#define BOOST_TEST_MODULE PCPVertexTestTransitions

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

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
        return PCPVertex("PCPVertex", pp, {},
                         std::make_tuple(time_energy_adaptor));
    }
    else {
        Utopia::PseudoParent pp("test.yml");
        return PCPVertex("PCPVertex", pp, {},
                         std::make_tuple(time_energy_adaptor));
    }
}

BOOST_FIXTURE_TEST_SUITE (test_PCPVertex_transitions, ModelFixture)

    void test_divide_cell(bool periodic) {
        auto model = model_factory(periodic);

        const auto& am = model.get_am();

        const auto cells = am.cells();
        const auto edges = am.edges();
        const auto vertices = am.vertices();

        auto cell = cells[cells.size() / 2];

        // set arbitrary values and check inheritance
        cell->state.area_preferential = 0.314;
        cell->state.area_preferential_var = 0.;
        cell->state.type = PCPVertex::CellType::support;
        cell->state.shape_index_preferential = 4.;
        cell->state.contractility = 0.114;
        model.divide_cell(cell, 0.);

        const auto cells_new = am.cells();
        const auto edges_new = am.edges();
        const auto vertices_new = am.vertices();

        BOOST_TEST(cells.size() + 1 == cells_new.size());
        BOOST_TEST(edges.size() + 3 == edges_new.size());
        BOOST_TEST(vertices.size() + 2 == vertices_new.size());

        for (auto new_cell : {cells_new[cells_new.size() - 2],
                              cells_new[cells_new.size() - 1]})
        {
            BOOST_TEST(new_cell->state.area_preferential
                       ==  cell->state.area_preferential);
            BOOST_TEST(new_cell->state.area_preferential_var
                       ==  cell->state.area_preferential_var);
            BOOST_TEST(new_cell->state.type
                       ==  cell->state.type);
            BOOST_TEST(new_cell->state.shape_index_preferential
                       ==  cell->state.shape_index_preferential);
            BOOST_TEST(new_cell->state.contractility
                       ==  cell->state.contractility);
        }

        test_custom_links(model);
        
        // check that it does not fail somewhere ...
        model.run();
    }

    // BOOST_AUTO_TEST_CASE(test_divide_cell_non_periodic) {
    //     test_divide_cell(false);
    // }
    // FIXME requires activation

    BOOST_AUTO_TEST_CASE(test_divide_cell_periodic) {
        test_divide_cell(true);
    }

    void test_T1_transition(bool periodic) {
        auto model = model_factory(periodic);
        model.prolog();

        const auto& am = model.get_am();

        const auto& cells = am.cells();
        const auto& edges = am.edges();
        const auto& vertices = am.vertices();

        auto num_cells = cells.size();
        auto num_edges = edges.size();
        auto num_vertices = vertices.size();

        auto edge = edges[edges.size() / 3];
        edge->state.linetension = 200.;

        auto a = edge->custom_links().a;
        auto b = edge->custom_links().b;

        PCPVertex::SpaceVec center = .5*(am.position_of(a) + am.position_of(b));
        am.move_to(a, 0.999*center);
        am.move_to(b, center);

        model.iterate();

        BOOST_TEST(edge.use_count() == 1);

        BOOST_TEST(cells.size() == num_cells);
        BOOST_TEST(edges.size() == num_edges);
        BOOST_TEST(vertices.size() == num_vertices);

        test_custom_links(model);
    }

    // BOOST_AUTO_TEST_CASE(test_T1_non_periodic) {
    //     test_T1_transition(false);
    // }
    // FIXME requires activation

    BOOST_AUTO_TEST_CASE(test_T1_periodic) {
        test_T1_transition(true);
    }

    void test_T2_transition(bool periodic) {
        auto model = model_factory(periodic);
        
        auto& am = model.get_am();

        auto& cells = am.cells();
        auto& edges = am.edges();
        auto& vertices = am.vertices();

        auto num_cells = cells.size();
        auto num_edges = edges.size();
        auto num_vertices = vertices.size();

        auto cell = cells[cells.size() / 2 + 4];
        cell->state.area_preferential = 0.1;
        
        model.prolog();
        for (int i = 0; i < 10; i++) {
            model.iterate();
        }

        bool found = std::find(cells.begin(), cells.end(), cell) != cells.end();
        BOOST_TEST(not found);

        BOOST_TEST(cell.use_count() == 1);

        BOOST_TEST(cells.size() == num_cells - 1);
        BOOST_TEST(edges.size() == (  num_edges
                                    - cell->custom_links().edges.size()));
        BOOST_TEST(vertices.size() == (  num_vertices
                                       - cell->custom_links().vertices.size()
                                       + 1));

        test_custom_links(model);
    }

    // BOOST_AUTO_TEST_CASE(test_T2_non_periodic) {
    //     test_T2_transition(false);
    // }
    // FIXME requires activation

    // BOOST_AUTO_TEST_CASE(test_T2_periodic) {
    //     test_T2_transition(true);
    // }
    // FIXME requires activation



BOOST_AUTO_TEST_SUITE_END()
