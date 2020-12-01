#include <iostream>

#include "PCPVertex.hh"
#include "energy.hh"
#include "algorithm.hh"
#include "operations.hh"
#include "PCPVertex_write_tasks.hh"

using namespace Utopia::Models::PCPVertex;
using namespace DataIO;
using Utopia::get_as;

/// Factory for model 
template<typename ParentType>
auto model_factory(ParentType &parent) {
    return PCPVertex("PCPVertex", parent, {}, std::make_tuple(
        // the energy adaptors
        time_energy_adaptor, energy_adaptor,
        linetension_adaptor, areaelasticity_adaptor,
        cell_contractility_adaptor, edge_contractility_adaptor,
        cell_cell_polarity_adaptor, polarity_exclusion_adaptor,
        lagrange_net_polarisation_adaptor, lagrange_const_concentration_adaptor,
        // transition adaptors
        statistics_time_adaptor,
        T1_adaptor, T1_attempted_adaptor, T2_adaptor,
        // the position adaptors
        vertices_adaptor<typename PCPVertex::Space::SpaceVec>,
        cells_adaptor<typename PCPVertex::Space::SpaceVec,
                      typename PCPVertex::CellType>,
        edges_adaptor,
        cell_energies_adaptor, edge_energies_adaptor
        ));
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
