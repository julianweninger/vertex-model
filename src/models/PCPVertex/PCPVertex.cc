#include <iostream>

#include <utopia/data_io/data_manager/defaults.hh>

#include "PCPVertex.hh"
#include "algorithm.hh"
#include "operations.hh"
#include "PCPVertex_write_tasks.hh"

using namespace Utopia::Models::PCPVertex;
using namespace DataIO;
using Utopia::get_as;

/// Factory for model 
template<typename ParentType>
auto model_factory(ParentType &parent) {
    auto deciders = Utopia::DataIO::Default::default_deciders< PCPVertex >;
    deciders["equilibrium_decider"] =
    []() -> std::shared_ptr<Utopia::DataIO::Default::Decider<PCPVertex>> {
        return std::make_shared<DeciderEquilibriumCondition<PCPVertex>>();
    };
    return PCPVertex("PCPVertex", parent, {}, std::make_tuple(
        // the energy adaptors
        time_energy_adaptor, energy_adaptor,
        // transition adaptors
        transition_adaptor,
        // statistics
        interface_length_adaptor,
        // the position adaptors
        vertices_adaptor<typename PCPVertex::Space::SpaceVec>,
        cells_adaptor<typename PCPVertex::Space::SpaceVec>,
        edges_adaptor<typename PCPVertex::Space::SpaceVec>,
        cell_energies_adaptor, edge_energies_adaptor
        ),
        deciders);
}

int main (int, char** argv) {
    try {
        // Initialize the PseudoParent from config file path
        Utopia::PseudoParent pp(argv[1]);
        
        auto model = model_factory(pp);
        model.run(); 

        // Done
        return 0;
    }
    catch (Utopia::Exception& e) {
        return Utopia::handle_exception(e);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Exception occurred!" << std::endl;
        return 1;
    }
}
