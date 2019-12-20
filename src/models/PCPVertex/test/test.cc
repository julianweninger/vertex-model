#define BOOST_TEST_MODULE PCPVertexTest

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include "../PCPVertex.hh"

using namespace Utopia;
using namespace Utopia::Models::PCPVertex;


// Factory for model
template<bool periodic_bc>
PCPVertex<periodic_bc> model_factory(std::string cfg) {
    // initialize pp inside to avoid problems with datamanager and memory 
    // access
    PseudoParent pp(cfg);

    auto model_cfg = pp.get_cfg()["PCPVertex"];

    BOOST_TEST(get_as<bool>("periodic_bc", model_cfg) == periodic_bc);

    return PCPVertex<periodic_bc>("PCPVertex", pp);
}

template<bool periodic_bc>
void test_T1_transition (std::string cfg)
{
    auto model = model_factory<periodic_bc>(cfg);
    auto model_cfg = model.get_cfg();

    auto edges = model.get_edges();

    // pick some inner edge
    auto e = edges[edges.size() / 3];

    // contract this inner edge
    double linetension = get_as<double>("linetension", model_cfg);
    e.lock()->linetension = 10*linetension;

    // run the model
    model.run();

    // check that the edge e has been removed
    BOOST_TEST(e.expired());

    // check that the number of egdes did not change
    BOOST_TEST(model.get_edges().size() == edges.size());

    model.get_logger()->info("Tear it down");
    spdlog::drop_all();
    std::remove("test_data.h5");
}

BOOST_AUTO_TEST_CASE(T1_transition_periodic)
{
    test_T1_transition<true>("test_transition_periodic.yml");
}

BOOST_AUTO_TEST_CASE(T1_transition_non_periodic)
{
    test_T1_transition<false>("test_transition_non_periodic.yml");
}



template<bool periodic_bc>
void test_T2_transition (std::string cfg)
{
    auto model = model_factory<periodic_bc>(cfg);
    auto model_cfg = model.get_cfg();

    auto cells = model.get_cells();

    // pick some inner cell
    auto c = cells[cells.size() / 2 + 1];

    // contract this inner cell
    double area_threshold = get_as<double>("area_threshold",
                                                    model_cfg);
    c.lock()->area_preferential = 0.5 * area_threshold;

    // run the model
    model.run();

    // check that the edge e has been removed
    BOOST_TEST(c.expired());

    // check that the number of egdes did not change
    BOOST_TEST(model.get_cells().size() == cells.size() - 1);

    model.get_logger()->info("Tear it down");
    spdlog::drop_all();
    std::remove("test_data.h5");
}

BOOST_AUTO_TEST_CASE(T2_transition_periodic)
{
    test_T2_transition<true>("test_transition_periodic.yml");
}

BOOST_AUTO_TEST_CASE(T2_transition_non_periodic)
{
    test_T2_transition<false>("test_transition_non_periodic.yml");
}


