#define BOOST_TEST_MODULE PCPVertexTest

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include "utils.hh"
#include <utopia/core/apply.hh>

#include "../PCPVertex.hh"
#include "../energy.hh"
#include "../algorithm.hh"
#include "../operations.hh"
#include "../PCPVertex_write_tasks.hh"

using namespace Utopia;
using namespace Utopia::Models::PCPVertex;

const double precision = 1e-12;

PCPVertex model_factory(bool periodic) {
    using Utopia::Models::PCPVertex::DataIO::time_energy_adaptor;

    if (periodic) {
        Utopia::PseudoParent pp("test_periodic.yml");
        return PCPVertex("PCPVertex", pp, {},
                         std::make_tuple(time_energy_adaptor));
    }
    else {
        Utopia::PseudoParent pp("test.yml");
        return PCPVertex("PCPVertex", pp, {},
                         std::make_tuple(time_energy_adaptor));
    }
}

class TEST_PCPVertex_energy_prediction : public PCPVertex
{
public:
    template<class ParentModel>
    TEST_PCPVertex_energy_prediction (
        const std::string name,
        ParentModel &parent_model,
        const Utopia::DataIO::Config custom_cfg = {})
    :
        PCPVertex(name, parent_model, custom_cfg,
             std::make_tuple(Utopia::Models::PCPVertex::DataIO::time_energy_adaptor))
    {
        this->prolog();

        // test the energy prediction throughout time
        for (int i = 0; i < 1000; i++) {
            test_energy_prediction();
            // NOTE involves iteration

            if ((i % 100) == 0) {
                this->jiggle_vertices(0.2);
            }
        }
        
        // compare dE with energy prediction
        this->jiggle_vertices(0.02);
        for (std::size_t i = 0; i < 500; i++) {
            test_energy_gradient();
        }
    }

    /// Test that the predicted energy and that obtained by iteration are equal
    /** Predicted energy: energy(beta = dt). Is compared to energy after
     *  iteration with step size dt.
     */
    void test_energy_prediction() {
        const double precision = 1e-10;

        std::string update_scheme(get_as<std::string>("update_scheme",
                                                        this->_cfg));
        BOOST_TEST(update_scheme == "steepest_gradient",
            "Test of energy prediction relies on fixed step size!");
        
        const auto dt = get_as<double>("dt", this->_cfg);

        // set the gradient
        this->init_minimization();

        // predict the energy terms
        auto energy = this->get_energy(dt);
        auto linetension = this->get_energy_linetension(dt);
        auto e_contr = this->get_energy_edge_contractility(dt);
        auto area_elast = this->get_energy_areaelasticity(dt);
        auto c_contr = this->get_energy_cell_contractility(dt);
        
        this->iterate();

        // test the tracing of energy
        BOOST_CHECK_CLOSE(this->_energy, this->get_energy(), precision);

        // energy prediction does not work in case of topological transitions
        if (this->get_num_T1s() == 0 and this->get_num_T2s() == 0) {
            // test the predictions
            BOOST_CHECK_CLOSE(energy, this->get_energy(),
                              precision);
            BOOST_CHECK_CLOSE(linetension, this->get_energy_linetension(),
                              precision);
            BOOST_CHECK_CLOSE(e_contr, this->get_energy_edge_contractility(),
                              precision);
            BOOST_CHECK_CLOSE(area_elast, this->get_energy_areaelasticity(),
                              precision);
            BOOST_CHECK_CLOSE(c_contr, this->get_energy_cell_contractility(),
                              precision);

            // test that energy is deterministic value
            std::vector<bool> energies_equal(100);
            std::transform(energies_equal.begin(), energies_equal.end(),
                        energies_equal.begin(),
                        [this, energy, precision](const auto&) {
                            return (  std::abs(this->get_energy() - energy)
                                    < precision);
                        });        

            BOOST_TEST(   energies_equal
                       == std::vector<bool>(energies_equal.size(), true));
        }
    }

    void test_energy_gradient() {
        const double precision = 1.e-2; // relative error in percent
        // NOTE only doing linear estimation here
        
        // chose a step size in steepest gradient update
        double dt = 0.00001;
        // NOTE error is expected to be smaller with small step size

        std::string update_scheme(get_as<std::string>("update_scheme",
                                                        this->_cfg));
        BOOST_TEST(update_scheme == "steepest_gradient",
            "Test of energy prediction relies on fixed step size!");
        
        // set the gradient
        this->init_minimization();

        // the current energy
        auto energy = this->get_energy();

        const auto& am = this->get_am();

        // A [0,N]-range uniform distribution used for evaluating probabilities
        std::uniform_int_distribution<> prob_distr(0, am.vertices().size()-1);

        // move a random vertex
        const auto& random_vertex = am.vertices()[prob_distr(*this->_rng)];
        const SpaceVec pos = am.position_of(random_vertex);
        this->get_am().move_by(random_vertex, random_vertex->state.f * dt);

        const SpaceVec d_pos = this->get_space()->displacement(
            pos, am.position_of(random_vertex));

        // integrate the gradient to obtain the energy change
        double dE = arma::dot(-1. * random_vertex->state.f, d_pos);

        // compare to the change in energy
        BOOST_CHECK_CLOSE(dE, this->get_energy() - energy, precision);
    }
};

BOOST_FIXTURE_TEST_SUITE (test_PCPVertex, ModelFixture)

    void test_model_initialization (PCPVertex& model) {
        test_custom_links(model);
    }

    // BOOST_AUTO_TEST_CASE(test_initialization_non_periodic) {
    //     auto model = model_factory(false);
    //     test_model_initialization(model);
    // }
    // FIXME requires activation

    BOOST_AUTO_TEST_CASE(test_initialization_periodic) {
        auto model = model_factory(true);
        test_model_initialization(model);
    }


    void test_model_minimization (PCPVertex& model)
    {
        model.prolog();

        Utopia::DataIO::Config minimization_cfg;
        minimization_cfg["tolerance"] = 1e-6;
        minimization_cfg["max_steps"] = 10000;
        minimization_cfg["num_repeat"] = 1;
        minimization_cfg["jiggle_tolerance"] = 5e-6;
        minimization_cfg["jiggle_intensity"] = 0.01;
        MinimizationParams minimization(minimization_cfg);
        model.minimize_energy(minimization);

        BOOST_TEST(model.get_time() > 0);
        BOOST_TEST(model.get_rel_energy_change() < minimization.tolerance);
    }

    // BOOST_AUTO_TEST_CASE(test_minimization_non_periodic)
    // {
    //     auto model = model_factory(false);
    //     test_model_minimization(model);
    // }
    // FIXME requires activation

    BOOST_AUTO_TEST_CASE(test_minimization_periodic)
    {
        auto model = model_factory(true);
        test_model_minimization(model);
    }

    /// Test the energy calculation
    BOOST_AUTO_TEST_CASE(test_energy_prediction)
    {
        Utopia::PseudoParent pp("test_periodic.yml");
        TEST_PCPVertex_energy_prediction test_model("PCPVertex", pp);
    }


BOOST_AUTO_TEST_SUITE_END()
