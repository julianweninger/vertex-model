#ifndef UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH
#define UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH

#include "../PCPVertex.hh"
#include "../initialisation.hh"
#include "../transitions.hh"
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

#endif