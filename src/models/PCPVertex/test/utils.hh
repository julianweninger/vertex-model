#ifndef UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH
#define UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH

#include <assert.h>
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
void test_custom_links(PCPVertex& model) {
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

#endif // UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH
