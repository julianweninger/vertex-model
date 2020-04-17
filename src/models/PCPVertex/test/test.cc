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

BOOST_FIXTURE_TEST_SUITE (test_PCPVertex_transitions, ModelFixture)

    void test_model_initialisation (PCPVertex& model) {
        test_custom_links(model);
    }

    // BOOST_AUTO_TEST_CASE(test_initialization_non_periodic) {
    //     auto model = model_factory(false);
    //     test_model_initialisation(model);
    // }
    // FIXME requires activation

    BOOST_AUTO_TEST_CASE(test_initialization_periodic) {
        auto model = model_factory(true);
        test_model_initialisation(model);
    }


BOOST_AUTO_TEST_SUITE_END()
