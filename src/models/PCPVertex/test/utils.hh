#ifndef UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH
#define UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH

#include <assert.h>
#include <math.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

/// A fixture used in the test_PCPVertex test suite
/** Besides destructing the model it is also necessary to destruct the logger
 *  and to free the datapath to construct a new model from the same config
 *
 *  \note the models output_path needs to be "test_data.h5"
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


/// This tests the custom_links of all agents living in the model
/** \details Tested are the neighborhood information stored in custom_links.
 *           These must be correct at all times.
 *  \note    Not tested is currently the order of vertices associated to a cell.
 *           This order is not maintained through some transitions.
 */
template <class PCPVertex>
void test_custom_links_periodic(PCPVertex& model) {
    const auto& am = model.get_am();

    for (auto v : am.vertices()) {
        BOOST_TEST(not v->state.remove);

        // test the adjacent edges
        auto adj_edges = am.adjoint_edges_of(v);
        BOOST_TEST(adj_edges.size() == 3);
        for (auto e : adj_edges) {
            BOOST_TEST(e);
            BOOST_TEST(not e->state.remove);
        }

        // test the adjacent cells
        auto adj_cells = am.adjoint_cells_of(v);
        BOOST_TEST(adj_cells.size() == 3);
        for (auto c : adj_cells) {
            BOOST_TEST(c);
            BOOST_TEST(not c->state.remove);
        }
    }

    for (auto e : am.edges()) {
        BOOST_TEST(not e->state.remove);

        // test the vertices
        BOOST_TEST(e->custom_links().a);
        BOOST_TEST(not e->custom_links().a->state.remove);
        BOOST_TEST(e->custom_links().b);
        BOOST_TEST(not e->custom_links().b->state.remove);

        auto adj_edges = am.adjoint_edges_of(e->custom_links().a);
        BOOST_TEST((std::find(adj_edges.begin(), adj_edges.end(), e) != 
                    adj_edges.end()));
        adj_edges = am.adjoint_edges_of(e->custom_links().b);
        BOOST_TEST((std::find(adj_edges.begin(), adj_edges.end(), e) != 
                    adj_edges.end()));

        // test the adjoint cells
        const auto [a, b] = am.adjoints_of(e);
        BOOST_TEST(a);
        BOOST_TEST(not a->state.remove);
        BOOST_TEST(b);
        BOOST_TEST(not b->state.remove);
    }

    for (auto c : am.cells()) {
        // test the cell
        BOOST_TEST(not c->state.remove);

        double area = am.area_of(c);
        BOOST_TEST(not isnan(area));
        BOOST_TEST(area > 0);

        BOOST_TEST(am.perimeter_of(c) > 0.);
        
        const auto& edges = c->custom_links().edges;
        const auto& vertices = c->custom_links().vertices;
        
        // test vertices
        BOOST_TEST(vertices.size() == edges.size());
        BOOST_TEST(vertices.size() >= 3);
        for (auto v : vertices) {
            BOOST_TEST(v);
            BOOST_TEST(not v->state.remove);
        }

        // test edges        
        auto [e, flip] = edges[0];
        auto first = e->custom_links().a;
        auto iterator = e->custom_links().b;
        if (flip) {
            std::swap(first, iterator);
        }

        for (unsigned int i = 1; i < edges.size(); i++) {
            auto [e, flip] = edges[i];
            BOOST_TEST(e);
            BOOST_TEST(not e->state.remove);

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

        // test cell neighbors
        BOOST_TEST(am.neighbors_of(c).size() == edges.size());
    }
}


/// This tests the custom_links of all agents living in the model
/** \details Tested are the neighborhood information stored in custom_links.
 *           These must be correct at all times.
 *  \note    Not tested is currently the order of vertices associated to a cell.
 *           This order is not maintained through some transitions.
 */
template <class PCPVertex>
void test_custom_links_non_periodic(PCPVertex& model) {
    const auto& am = model.get_am();

    for (const auto& v : am.vertices()) {
        BOOST_TEST(v);
        BOOST_TEST(not v->state.remove);

        // test the adjacent edges
        auto adj_edges = am.adjoint_edges_of(v);
        if (am.is_2_fold_boundary_vertex(v)) {
            BOOST_TEST(adj_edges.size() == 2);
        }
        else {
            BOOST_TEST(adj_edges.size() == 3);
        }
        for (auto e : adj_edges) {
            BOOST_TEST(e);
            BOOST_TEST(not e->state.remove);
        }

        // test the adjacent cells
        auto adj_cells = am.adjoint_cells_of(v);        
        if (am.is_2_fold_boundary_vertex(v)) {
            BOOST_TEST(adj_cells.size() == 1);
        }
        else if (am.is_3_fold_boundary_vertex(v)) {
            BOOST_TEST(adj_cells.size() == 2);
        }
        else {
            BOOST_TEST(adj_cells.size() == 3);
        }
        for (auto c : adj_cells) {
            BOOST_TEST(c);
            BOOST_TEST(not c->state.remove);
        }
    }

    for (auto e : am.edges()) {
        BOOST_TEST(e);
        BOOST_TEST(not e->state.remove);

        // test the vertices
        BOOST_TEST(e->custom_links().a);
        BOOST_TEST(not e->custom_links().a->state.remove);
        BOOST_TEST(e->custom_links().b);
        BOOST_TEST(not e->custom_links().b->state.remove);

        auto adj_edges = am.adjoint_edges_of(e->custom_links().a);
        BOOST_TEST((std::find(adj_edges.begin(), adj_edges.end(), e) != 
                    adj_edges.end()));
        adj_edges = am.adjoint_edges_of(e->custom_links().b);
        BOOST_TEST((std::find(adj_edges.begin(), adj_edges.end(), e) != 
                    adj_edges.end()));

        // test the adjoint cells
        const auto [a, b] = am.adjoints_of(e);
        if (am.is_1_cell_boundary_edge(e)) {
            BOOST_TEST((a or b));
            if (not a) {
                BOOST_TEST(b);
                BOOST_TEST(not b->state.remove);
            }
            else if (not b) {
                BOOST_TEST(a);
                BOOST_TEST(not a->state.remove);
            }
        }
        else {
            BOOST_TEST(a);
            BOOST_TEST(not a->state.remove);
            BOOST_TEST(b);
            BOOST_TEST(not b->state.remove);
        }
    }

    for (auto c : am.cells()) {
        BOOST_TEST(c);
        BOOST_TEST(not c->state.remove);

        const auto& edges = c->custom_links().edges;
        const auto& vertices = c->custom_links().vertices;
        
        // test the cell
        double area = am.area_of(c);
        BOOST_TEST(not isnan(area));
        BOOST_TEST(area > 0);

        BOOST_TEST(am.perimeter_of(c) > 0.);

        // test vertices
        BOOST_TEST(vertices.size() == edges.size());
        BOOST_TEST(vertices.size() >= 3);
        for (auto v : vertices) {
            BOOST_TEST(v);
            BOOST_TEST(not v->state.remove);
        }
        
        // test edges        
        auto [e, flip] = edges[0];
        BOOST_TEST(e);
        BOOST_TEST(not e->state.remove);
        auto first = e->custom_links().a;
        auto iterator = e->custom_links().b;
        BOOST_TEST(first);
        BOOST_TEST(not first->state.remove);
        BOOST_TEST(iterator);
        BOOST_TEST(not iterator->state.remove);
        if (flip) {
            std::swap(first, iterator);
        }

        for (unsigned int i = 1; i < edges.size(); i++) {
            auto [e, flip] = edges[i];
            BOOST_TEST(e);
            BOOST_TEST(not e->state.remove);
            BOOST_TEST(e->custom_links().a);
            BOOST_TEST(not e->custom_links().a->state.remove);
            BOOST_TEST(e->custom_links().b);
            BOOST_TEST(not e->custom_links().b->state.remove);

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

        // test cell neighbors
        if (not am.is_boundary(c)) {
            BOOST_TEST(am.neighbors_of(c).size() == edges.size());
        }
        else {
            BOOST_TEST(am.neighbors_of(c).size() < edges.size());
        }
    }
}


/// This tests the custom_links of all agents living in the model
/** \details Tested are the neighborhood information stored in custom_links.
 *           These must be correct at all times.
 *  \note    Not tested is currently the order of vertices associated to a cell.
 *           This order is not maintained through some transitions.
 */
template <class PCPVertex>
void test_custom_links(PCPVertex& model) {
    model.get_logger()->info("Testing model custom links ...");
    if (model.get_space()->periodic) {
        test_custom_links_periodic(model);
    }
    else {
        test_custom_links_non_periodic(model);
    }
    model.get_logger()->info("Ended test of model custom links.");

    model.get_logger()->info("Testing Vertex model energy terms ...");
    
    BOOST_TEST( std::isfinite(model.get_energy_linetension()) );
    BOOST_TEST( std::isfinite(model.get_energy_edge_contractility()) );
    BOOST_TEST( std::isfinite(model.get_energy_areaelasticity()) );
    BOOST_TEST( std::isfinite(model.get_energy_cell_contractility()) );
    
    BOOST_TEST( std::isfinite(model.get_boundary_shape_energy()) );
    BOOST_TEST( std::isfinite(model.get_boundary_area_energy()) );
    
    BOOST_TEST( std::isfinite(model.get_energy()) );

    model.get_logger()->info("Ended test of model energy terms.");
}

#endif // UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH
