#include <iostream>

#include "PCPTopology.hh"
#include "PCPTopology_write_tasks.hh"
#include "../PlanarCellPolarity/PlanarCellPolarity_write_tasks.hh"

using namespace Utopia::Models::PCPVertex;
using namespace DataIO;
using Utopia::get_as;

/// Factory for model 
template<typename ParentType>
auto model_factory(ParentType parent) {
    return PCPTopology("PCPTopology", parent, {}, std::make_tuple(
        // energy adaptors
        continuous_time_adaptor, time_energy_adaptor, energy_adaptor,
        linetension_adaptor, areaelasticity_adaptor,
        cell_contractility_adaptor, edge_contractility_adaptor,
        boundary_area_elasticity_adaptor, boundary_shape_elasticity_adaptor,
        // transitions
        transition_adaptor,
        // statistics
        statistics_time_adaptor, cell_neighbourhood_adaptor,
        cell_area_histogram_adaptor,
        cell_stats_adaptor<typename PCPTopology::Cell>,
        hair_cell_stats_adaptor<typename PCPTopology::Cell>,
        support_cell_stats_adaptor<typename PCPTopology::Cell>,
        bulk_cell_stats_adaptor<typename PCPTopology::Cell>,
        bulk_hair_cell_stats_adaptor<typename PCPTopology::Cell>,
        bulk_support_cell_stats_adaptor<typename PCPTopology::Cell>,
        interface_length_adaptor,
        // the position adaptors
        vertices_adaptor<typename PCPVertex::Space::SpaceVec>,
        cells_adaptor<typename PCPVertex::Space::SpaceVec,
                      typename PCPTopology::CellType>,
        edges_adaptor<typename PCPVertex::Space::SpaceVec>,
        cell_energies_adaptor, edge_energies_adaptor,
        cell_cluster_adaptor,
        // PlanarCellPolarity spatial data
        Utopia::Models::PlanarCellPolarity::DataIO::pcp_cells_adaptor
        <typename PCPVertex::Space::SpaceVec, typename PCPTopology::ProteinVec>,
        Utopia::Models::PlanarCellPolarity::DataIO::pcp_edges_adaptor
        <typename PCPTopology::ProteinVec>
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
