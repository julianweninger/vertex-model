#define BOOST_TEST_MODULE PCPVertexTest

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include "../PCPVertex.hh"
#include "../PCPVertex_write_tasks.hh"

using namespace Utopia;
using namespace Utopia::Models::PCPVertex;


// Factory for model
template<bool periodic_bc>
PCPVertex<periodic_bc, false> model_factory(std::string cfg) {
    // initialize pp inside to avoid problems with datamanager and memory 
    // access
    PseudoParent pp(cfg);

    auto model_cfg = pp.get_cfg()["PCPVertex"];

    BOOST_TEST(get_as<bool>("periodic_bc", model_cfg) == periodic_bc);

    using Utopia::Models::PCPVertex::DataIO::time_energy_adaptor;

    return PCPVertex<periodic_bc, false>("PCPVertex", pp, time_energy_adaptor);
}

/// Destructor for the model
/** Besides destructing the model it is also necessary to destruct the logger
 *  and to free the datapath to construct a new model from the same config
 * 
 *  NOTE the models output_path needs to be "test_data.h5" 
 * 
 */
template<typename Model>
void destruct_model_factory(Model model) 
{
    model.get_logger()->info("Tear the model down");
    spdlog::drop_all();
    std::remove("test_data.h5");
}


template<bool periodic_bc>
void test_initialisation_hexagon (std::string cfg)
{
    auto model = model_factory<periodic_bc>(cfg);
    auto model_cfg = model.get_cfg();

    auto vertices = model.get_vertices();
    auto edges = model.get_edges();
    auto cells = model.get_cells();

    if constexpr (periodic_bc) {
        for (auto v : vertices) {
            assert(v.lock()->adj_edges.size() == 3);
            assert(v.lock()->adj_cells.size() == 3);
        }
        for (auto e : edges) {
            assert(not e.lock()->adj_cell_a.expired() and
                   not e.lock()->adj_cell_b.expired());
        }
        for (auto c : cells) {
            assert (c.lock()->vertices.size() == 6);
            assert (c.lock()->edges_ordered.size() == 6);
        }
    }
    // else {
    //     for (auto v : vertices) {
    //         assert(v.lock()->adj_edges.size() <= 3);
    //         assert(v.lock()->adj_cells.size() <= 3);
    //     }
    //     for (auto e : edges) {
    //         assert(not e.lock()->adj_cell_a.expired() or 
    //                not e.lock()->adj_cell_b.expired());
    //     }
    //     for (auto c : cells) {
    //         assert (c.lock()->vertices.size() == 6);
    //         assert (c.lock()->edges_ordered.size() == 6);
    //     }
    // }

    destruct_model_factory(model);
}

BOOST_AUTO_TEST_CASE(Initialisation_hexagon_periodic)
{
    test_initialisation_hexagon<true>("test_periodic.yml");
}

BOOST_AUTO_TEST_CASE(Initialisation_hexagon_non_periodic)
{
    test_initialisation_hexagon<false>("test_non_periodic.yml");
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

    destruct_model_factory(model);
}

BOOST_AUTO_TEST_CASE(T1_transition_periodic)
{
    test_T1_transition<true>("test_periodic.yml");
}

BOOST_AUTO_TEST_CASE(T1_transition_non_periodic)
{
    test_T1_transition<false>("test_non_periodic.yml");
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
    double area_threshold = get_as<double>("area_threshold", model_cfg);
    c.lock()->area_preferential = 0.5 * area_threshold;

    // run the model
    model.run();

    // check that the edge e has been removed
    BOOST_TEST(c.expired());

    // check that the number of egdes did not change
    BOOST_TEST(model.get_cells().size() == cells.size() - 1);

    destruct_model_factory(model);
}

BOOST_AUTO_TEST_CASE(T2_transition_periodic)
{
    test_T2_transition<true>("test_periodic.yml");
}

BOOST_AUTO_TEST_CASE(T2_transition_non_periodic)
{
    test_T2_transition<false>("test_non_periodic.yml");
}



template<bool periodic_bc>
void test_cell_division (std::string cfg)
{
    std::vector<int> test_cases;
    if constexpr (periodic_bc) {
        test_cases = {15, 1, 4, 5};
    }
    else { test_cases = {5}; }

    for (int cell_pos : test_cases) {
    auto model = model_factory<periodic_bc>(cfg);

    auto vertices = model.get_vertices();
    auto edges = model.get_edges();
    auto cells = model.get_cells();

    model.divide_cell(cells[cell_pos].lock(), 0.); 
    // this is a cell not at the boundary


    auto new_vertices = model.get_vertices();
    auto new_edges = model.get_edges();
    auto new_cells = model.get_cells();

    // assert that the right number of objects has been added during division
    assert(cells.size() == new_cells.size()-1);
    assert(edges.size() == new_edges.size()-3);
    assert(vertices.size() == new_vertices.size()-2);

    for (auto v : new_vertices) {
        assert(not v.expired());
        assert(not v.lock()->remove);            
    }
    for (auto e : new_edges) {
        assert(not e.expired());
        assert(not e.lock()->remove);
    }
    for (auto c : new_cells) {
        assert(not c.expired());
        assert(not c.lock()->remove);
    }

    int cnt_expired = 0;
    for (auto c : cells) {
        cnt_expired += c.expired();
    }
    assert(cnt_expired == 1);
    cnt_expired = 0;
    for (auto e : edges) {
        cnt_expired += e.expired();
    }
    assert(cnt_expired == 2);
    cnt_expired = 0;
    for (auto v : vertices) {
        cnt_expired += v.expired();
    }
    assert(cnt_expired == 0);

    if constexpr (periodic_bc)
    {
        for (auto v : new_vertices) {
            assert(v.lock()->adj_edges.size() == 3);
            assert(v.lock()->adj_cells.size() == 3);
        }
        for (auto e : new_edges) {
            assert(not e.lock()->adj_cell_a.expired() and
                   not e.lock()->adj_cell_b.expired());
        }
    }

    model.run();


    destruct_model_factory(model);
    }
}

BOOST_AUTO_TEST_CASE(cell_division_non_periodic)
{
    test_cell_division<false>("test_non_periodic.yml");
}

BOOST_AUTO_TEST_CASE(cell_division_periodic)
{
    test_cell_division<true>("test_periodic.yml");
}

