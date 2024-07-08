#include <iostream>

#include "PCPTopology.hh"
#include "PCPTopology_write_tasks.hh"
#include "../PlanarCellPolarity/PlanarCellPolarity_write_tasks.hh"

using namespace Utopia::Models::PCPTopology;

using TopologyModel = PCPTopology;

/// Factory for model 
template<typename ParentType>
auto model_factory(ParentType parent) {
    using VertexModel = Utopia::Models::PCPVertex::PCPVertex;

    using namespace Utopia::Models;
    using namespace Utopia::Models::PCPVertex::DataIO;
    
    return TopologyModel("PCPTopology", parent, {}, std::make_tuple(
        // energy adaptors
        continuous_time_adaptor, time_energy_adaptor, energy_adaptor,
        // transitions
        transition_adaptor,
        // statistics
        interface_length_adaptor,
        // the position adaptors
        vertices_adaptor<VertexModel::SpaceVec>,
        cells_adaptor<VertexModel::SpaceVec>,
        edges_adaptor<VertexModel::SpaceVec>,
        cell_energies_adaptor, edge_energies_adaptor,
        cell_cluster_adaptor,
        // PlanarCellPolarity spatial data
        PlanarCellPolarity::DataIO::pcp_cells_adaptor<
            VertexModel::SpaceVec,
            typename TopologyModel::ProteinVec >,
        PlanarCellPolarity::DataIO::pcp_edges_adaptor<
            typename TopologyModel::ProteinVec>
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
