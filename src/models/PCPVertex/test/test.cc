#define BOOST_TEST_MODULE PCPVertexTest

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include "utils.hh"

template<bool periodic_bc>
void test_initialisation_hexagon (std::string cfg)
{
    auto model = model_factory<periodic_bc>(cfg);
    auto model_cfg = model.get_cfg();

    auto vertices = model.get_vertices();
    auto edges = model.get_edges();
    auto cells = model.get_cells();

    int num_rows = get_as<int>("lattice_rows", model_cfg);
    int num_columns = get_as<int>("lattice_columns", model_cfg);
    assert (cells.size() == num_rows * num_columns);

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
    else {
        for (auto v : vertices) {
            assert(v.lock()->adj_edges.size() <= 3);
            assert(v.lock()->adj_cells.size() <= 3);
        }
        for (auto e : edges) {
            assert(not e.lock()->adj_cell_a.expired() or 
                   not e.lock()->adj_cell_b.expired());
        }
        for (auto c : cells) {
            assert (c.lock()->vertices.size() == 6);
            assert (c.lock()->edges_ordered.size() == 6);
        }
    }

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

