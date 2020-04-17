#define BOOST_TEST_MODULE PCPVertexTestNew

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "../PCPVertex.hh"
#include "../energy.hh"
#include "../algorithm.hh"
#include "../operations.hh"
#include "../PCPVertex_write_tasks.hh"

using namespace Utopia::Models::PCPVertex;

PCPVertex model_factory(bool periodic) {
    if (periodic) {
        Utopia::PseudoParent pp("test_periodic.yml");
        return PCPVertex("PCPVertex", pp);
    }
    else {
        Utopia::PseudoParent pp("test.yml");
        return PCPVertex("PCPVertex", pp);
    }
}

/// A fixture used in the test_PCPVertex test suite
/** Besides destructing the model it is also necessary to destruct the logger
 *  and to free the datapath to construct a new model from the same config
 *
 *  NOTE the models output_path needs to be "test_data.h5"
 *
 */
struct ModelFixture {
    std::shared_ptr<spdlog::logger> log;

    // Construct
    ModelFixture ()
    :
        log([]() {
            auto logger = spdlog::get("test");

            // Create it only if it does not already exist
            if (not logger) {
                logger = spdlog::stdout_color_mt("test");
            }

            // Set level and global logging pattern
            logger->set_level(spdlog::level::debug);
            spdlog::set_pattern("[%T.%e] [%^%l%$] [%n]  %v");
            // "[HH:MM:SS.mmm] [level(colored)] [logger]  <message>"

            return logger;
        }())
    { }

    // Teardown, invoked after each test
    ~ModelFixture () {
        log->info("Tearing down ...");
        std::remove("test_data.h5");
        spdlog::drop_all();
    }
};

BOOST_FIXTURE_TEST_SUITE (test_PCPVertex_transitions, ModelFixture)

    /// This tests the custom_links of all agents living in the model
    /** \details Tested are the neighborhood information stored in custom_links.
     *           These must be correct at all times.
     *  \note    Not tested is currently the order of vertices associated to a cell.
     *           This order is not maintained through some transitions.
     */
    void test_custom_links(PCPVertex model) {
        const auto& am = model.get_am();

        for (auto v : am.vertices()) {
            auto adj_edges = am.adjoint_edges_of(v);
            BOOST_TEST(adj_edges.size() == 3);
            for (auto e : adj_edges) {
                BOOST_TEST(e);
            }

            auto adj_cells = am.adjoint_cells_of(v);
            BOOST_TEST(adj_cells.size() == 3);
            for (auto c : adj_cells) {
                BOOST_TEST(c);
            }
        }
        for (auto e : am.edges()) {
            BOOST_TEST(e->custom_links().a);
            BOOST_TEST(e->custom_links().b);

            auto adj_edges = am.adjoint_edges_of(e->custom_links().a);
            BOOST_TEST((std::find(adj_edges.begin(), adj_edges.end(), e) != 
                       adj_edges.end()));
            adj_edges = am.adjoint_edges_of(e->custom_links().b);
            BOOST_TEST((std::find(adj_edges.begin(), adj_edges.end(), e) != 
                       adj_edges.end()));

            auto [adj_a, adj_b] = am.adjoints_of(e);
            BOOST_TEST(adj_a);
            BOOST_TEST(adj_b);
        }
        for (auto c : am.cells()) {
            BOOST_TEST(am.area_of(c) > 0);

            const auto& vertices = c->custom_links().vertices;
            const auto& edges = c->custom_links().edges;
            for (auto v : vertices) {
                BOOST_TEST(v);
            }

            BOOST_TEST(vertices.size() == edges.size());
            BOOST_TEST(vertices.size() >= 3);
            
            auto [e, flip] = edges[0];
            auto first = e->custom_links().a;
            auto iterator = e->custom_links().b;
            if (flip) {
                std::swap(first, iterator);
            }

            for (unsigned int i = 1; i < edges.size(); i++) {
                auto [e, flip] = edges[i];
                BOOST_TEST(e);

                if (not flip) {
                    BOOST_TEST(iterator == e->custom_links().a);
                    iterator = e->custom_links().b;
                }
                else {
                    BOOST_TEST(iterator == e->custom_links().b);
                    iterator = e->custom_links().a;
                }
            }
            BOOST_TEST(iterator == first);
        }
    }

    void test_divide_cell(bool periodic) {
        auto model = model_factory(periodic);

        const auto& am = model.get_am();

        const auto& cells = am.cells();
        const auto& edges = am.edges();
        const auto& vertices = am.vertices();

        int num_cells = cells.size();
        int num_edges = edges.size();
        int num_vertices = vertices.size();

        auto cell = cells[cells.size() / 2];
        model.divide_cell(cell, 0.);

        BOOST_TEST(cells.size() == num_cells + 1);
        BOOST_TEST(edges.size() == num_edges + 3);
        BOOST_TEST(vertices.size() == num_vertices + 2);

        test_custom_links(model);

        model.run();
    }

    // BOOST_AUTO_TEST_CASE(test_divide_cell_non_periodic) {
    //     test_divide_cell(false);
    // }

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

    BOOST_AUTO_TEST_CASE(test_T1_periodic) {
        test_T1_transition(true);
    }

    void test_T2_transition(bool periodic) {
        auto model = model_factory(periodic);
        model.prolog();
        
        const auto& am = model.get_am();

        const auto& cells = am.cells();
        const auto& edges = am.edges();
        const auto& vertices = am.vertices();

        auto num_cells = cells.size();
        auto num_edges = edges.size();
        auto num_vertices = vertices.size();

        auto cell = cells[cells.size() / 2];

        cell->state.area_preferential = 0.;
        for (int i = 0; i < 1000; i++) {
            model.iterate();
            if (cells.size() > num_cells) { 
                break;
            }
        }

        BOOST_TEST(cell.use_count() == 1);

        BOOST_TEST(cells.size() == num_cells - 1);
        BOOST_TEST(edges.size() == num_edges-cell->custom_links().edges.size());
        BOOST_TEST(vertices.size() == 
                   num_vertices - cell->custom_links().vertices.size() + 1);

        test_custom_links(model);
    }

    // BOOST_AUTO_TEST_CASE(test_T2_non_periodic) {
    //     test_T2_transition(false);
    // }

    BOOST_AUTO_TEST_CASE(test_T2_periodic) {
        test_T2_transition(true);
    }



BOOST_AUTO_TEST_SUITE_END()
