#define BOOST_TEST_MODULE PCPVertexTest

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include "utils.hh"

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
        test_energy_prediction();
    }

    void test_energy_prediction() {
        this->prolog();

        std::string update_scheme(get_as<std::string>("update_scheme",
                                                        this->_cfg));
        BOOST_TEST(update_scheme == "steepest_gradient",
            "Test of energy prediction relies on fixed step size!");
        
        const auto dt = get_as<double>("dt", this->_cfg);
        const auto& am = this->get_am();

        // predict the energy terms
        auto new_energy = this->get_energy(dt);
        auto linetension = this->get_energy_linetension(am.edges(), dt);
        auto e_contr = this->get_energy_edge_contractility(am.edges(), dt);
        auto area_elast = this->get_energy_areaelasticity(am.cells(), dt);
        auto c_contr = this->get_energy_cell_contractility(am.cells(), dt);
        
        this->iterate();

        // test the predictions
        BOOST_CHECK_CLOSE(linetension, this->get_energy_linetension(),
                        precision);
        BOOST_CHECK_CLOSE(e_contr, this->get_energy_edge_contractility(),
                        precision);
        BOOST_CHECK_CLOSE(area_elast, this->get_energy_areaelasticity(),
                        precision);
        BOOST_CHECK_CLOSE(c_contr, this->get_energy_cell_contractility(),
                        precision);
        BOOST_CHECK_CLOSE(new_energy, this->get_energy(),
                        precision);
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
