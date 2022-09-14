#define BOOST_TEST_MODULE PCPVertexTest

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "utils.hh"

#include <utopia/core/apply.hh>

#include "../PCPVertex.hh"
#include "../algorithm.hh"
#include "../operations.hh"
#include "../PCPVertex_write_tasks.hh"

using namespace Utopia;
using namespace Utopia::Models::PCPVertex;

using SpaceVec = typename PCPVertex::SpaceVec;

enum Cases {
    periodic,
    non_periodic,
    shape_elastic, // non-periodic with shape elasticity term
    columnar_periodic,
    columnar      // Columnar initialisation  
};

template<Cases C>
struct Fixture {    
    Models::PCPVertex::PCPVertex vertex_model;

    Fixture ()
    :
        vertex_model(model_factory())
    { }

    ~Fixture()
    {
        vertex_model.get_logger()->info("Tearing down ...");
        std::remove("test_data.h5");
        spdlog::drop_all();
    }
    
    PCPVertex model_factory() {
        using Utopia::Models::PCPVertex::DataIO::time_energy_adaptor;

        if constexpr (C == periodic) {
            Utopia::PseudoParent pp("test_periodic.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else if constexpr (C == non_periodic) {
            Utopia::PseudoParent pp("test.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else if constexpr (C == shape_elastic) {
            Utopia::PseudoParent pp("test_shape_elastic.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else if constexpr (C == columnar) {
            Utopia::PseudoParent pp("test_columnar.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else if constexpr (C == columnar_periodic) {
            Utopia::PseudoParent pp("test_columnar_periodic.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else {
            throw std::runtime_error(fmt::format("Testcase {} not implemented",
                                                 C));
        }
    }
};


template<bool periodic>
struct FixtureColumnar {    
    Models::PCPVertex::PCPVertex vertex_model;

    FixtureColumnar ()
    :
        vertex_model(model_factory())
    { }

    ~FixtureColumnar()
    {
        vertex_model.get_logger()->info("Tearing down ...");
        std::remove("test_data.h5");
        spdlog::drop_all();
    }
    
    PCPVertex model_factory() {
        using Utopia::Models::PCPVertex::DataIO::time_energy_adaptor;

        if constexpr (periodic) {
            Utopia::PseudoParent pp("test_column_periodic.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else {
            Utopia::PseudoParent pp("test_column.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
    }
};

class TEST_PCPVertex_energy : public PCPVertex
{
public:
    template<class ParentModel>
    TEST_PCPVertex_energy (
        const std::string name,
        ParentModel &parent_model,
        const Utopia::DataIO::Config custom_cfg = {})
    :
        PCPVertex(name, parent_model, custom_cfg,
             std::make_tuple(Utopia::Models::PCPVertex::DataIO::time_energy_adaptor))
    { }

    /** Runs the full test */
    void run_test() {
        this->prolog();

        this->perform_test();

        this->epilog();
    }

    /** Performs the full test */
    void perform_test() {        
        // compare dE with energy
        this->jiggle_vertices(0.02);
        for (std::size_t i = 0; i < 100; i++) {
            test_energy_gradient();
        }
    }

    void test_energy_gradient() {
        const double precision = 1.e-2; // relative error in percent
        // NOTE only doing linear estimation here
        
        // chose a step size in steepest gradient update
        double dt = 0.00001;
        // NOTE error is expected to be smaller with small step size
        
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
        this->get_am().move_by(
            random_vertex,
            random_vertex->state.get_force() * dt
        );

        const SpaceVec d_pos = this->get_space()->displacement(
            pos, am.position_of(random_vertex));

        // integrate the gradient to obtain the energy change
        double dE = arma::dot(-1. * random_vertex->state.get_force(), d_pos);

        // compare to the change in energy
        BOOST_CHECK_CLOSE(dE, this->get_energy() - energy, precision);
    }
};

typedef boost::mpl::vector< Fixture<Cases::periodic>,
                            Fixture<Cases::non_periodic>,
                            Fixture<Cases::shape_elastic>,
                            Fixture<Cases::columnar_periodic>,
                            Fixture<Cases::columnar>
                          > Fixtures;

BOOST_FIXTURE_TEST_SUITE (test_PCPVertex, ModelFixture)

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_model_initialization, F, Fixtures) {
        F fixture;
        test_custom_links(fixture.vertex_model);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_minimization, F, Fixtures) {
        F fixture;
        auto& model = fixture.vertex_model;
        model.prolog();

        Utopia::DataIO::Config minimization_cfg;

        std::cout << std::endl << "Beginning test of minimization in "
                                  "steepest gradient scheme .. \n\n";

        minimization_cfg["tolerance"] = 1e-6;
        minimization_cfg["dt"] = 1e-2;
        minimization_cfg["max_steps"] = 10000;
        minimization_cfg["num_repeat"] = 2;
        minimization_cfg["jiggle_tolerance"] = 5e-6;
        minimization_cfg["jiggle_intensity"] = 0.001;
        MinimizationParams minimization(minimization_cfg);

        model.minimize_energy(minimization);

        BOOST_TEST(model.get_time() >= 5);
        BOOST_TEST(model.get_energy_change() < minimization.tolerance);


        // std::cout << std::endl << "Beginning test of minimization in "
        //                           "steepest gradient noisy scheme .. \n\n";

        // minimization_cfg["num_steps"] = 100;
        // minimization_cfg["temperature"] = 0.0001;
        // minimization_cfg["linetension_fluctuation"] = 0.0001;
        // minimization_cfg["area_fluctuation"] = 0.0001;

        // auto time_start = model.get_time();
        // minimization = MinimizationParams(minimization_cfg);
        // model.minimize_energy(minimization);

        // BOOST_TEST(model.get_time() - time_start == 2 * (100 + 1) + 1);
    }

    BOOST_AUTO_TEST_CASE(test_energy_periodic)
    {
        Utopia::PseudoParent pp("test_periodic.yml");
        TEST_PCPVertex_energy test_model("PCPVertex", pp);

        test_model.run_test();
    }

    BOOST_AUTO_TEST_CASE(test_energy)
    {
        Utopia::PseudoParent pp("test.yml");
        TEST_PCPVertex_energy test_model("PCPVertex", pp);

        test_model.prolog();

        test_model.perform_test();


        // std::cout << std::endl << "Adding heterogeneities .. \n\n";

        // test_model.set_ppMLC_contractility(0.1, SpaceVec({1., 0.}), false);
        // test_model.set_pMLC_contractility(0.2, 1.);

        // test_model.perform_test();


        // std::cout << std::endl << "Adding curvature .. \n\n";

        // test_model.curve_boundary(0.005, true);
        // test_model.set_ppMLC_contractility(0.1, SpaceVec({1., 0.}), true);

        // test_model.perform_test();

        std::cout << std::endl << "Done. \n\n";

        test_model.epilog();
    }


BOOST_AUTO_TEST_SUITE_END()
