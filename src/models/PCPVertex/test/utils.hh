#ifndef UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH
#define UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH

#include "../PCPVertex.hh"
#include "../initialisation.hh"
#include "../transitions.hh"
#include "../PCPVertex_write_tasks.hh"
#include <boost/test/unit_test.hpp>

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

template<typename Model>
void test_weak_links(Model model) {
    auto edges = model.get_edges();
    for (auto e_weak : edges) {
        auto e = e_weak.lock();
        for (auto v : {e->a, e->b}) {
            BOOST_TEST((std::find_if(v->adj_edges.begin(), v->adj_edges.end(),
                            [e](auto e_weak) {return e_weak.lock() == e; })
                        != v->adj_edges.end()));
        }
    }
    
    auto cells = model.get_cells();
    for (auto c_weak : cells) {
        auto c = c_weak.lock();
        for (auto v : c->vertices) {
            BOOST_TEST((std::find_if(v->adj_cells.begin(), v->adj_cells.end(),
                            [c](auto c_weak) {return c_weak.lock() == c; })
                       != v->adj_cells.end()));
        }
        for (auto [e, flip] : c->edges_ordered) {
            bool test = false;
            if (not e->adj_cell_a.expired() and e->adj_cell_a.lock() == c) {
                test = true;
            }
            if (not e->adj_cell_b.expired() and e->adj_cell_b.lock() == c) {
                test = true;
            }
            BOOST_TEST(test);
        }
    }
    
    auto vertices = model.get_vertices();
    // for (auto v : vertices) {
        // BOOST_TEST(v->adj_cells.size() <= 3); 
        // BOOST_TEST(v->adj_edges.size() <= 3);
        // TODO activate once T2 transition fixed
    // }
}

#endif