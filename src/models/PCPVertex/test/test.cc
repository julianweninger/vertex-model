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

using namespace Utopia::Models::PCPVertex;

PCPVertex model_factory(bool periodic) {
    using Utopia::Models::PCPVertex::DataIO::time_energy_adaptor;

    if (periodic) {
        Utopia::PseudoParent pp("test_periodic.yml");
        return PCPVertex("PCPVertex", pp, time_energy_adaptor);
    }
    else {
        Utopia::PseudoParent pp("test.yml");
        return PCPVertex("PCPVertex", pp, time_energy_adaptor);
    }
}

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

    BOOST_AUTO_TEST_CASE(test_minimization_periodic)
    {
        auto model = model_factory(true);
        test_model_minimization(model);
    }


BOOST_AUTO_TEST_SUITE_END()
