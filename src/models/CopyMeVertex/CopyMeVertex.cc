#include <iostream>

#include "../PCPTopology/PCPTopology.hh"
#include "../PCPTopology/PCPTopology_write_tasks.hh"
#include "../PlanarCellPolarity/PlanarCellPolarity_write_tasks.hh"

#include "operations.hh"
#include "work_function.hh"

using namespace Utopia::Models::PCPTopology;

using TopologyModel = PCPTopology;

using namespace Utopia::Models::CopyMeVertex::WorkFunction;
using namespace Utopia::Models::CopyMeVertex::OperationCollection;

/// Factory for model 
template<typename ParentType>
auto model_factory(ParentType parent) {
    using VertexModel = Utopia::Models::PCPVertex::PCPVertex;

    using namespace Utopia::Models;
    using namespace Utopia::Models::PCPVertex::DataIO;
    
    return TopologyModel("CopyMeVertex", parent, {}, std::make_tuple(
        // energy adaptors
        energy_adaptor,             // records the work-function
        time_energy_adaptor,        // the associated time
        continuous_time_adaptor,    // the associated time of PCPVertex
        // transitions
        transition_adaptor,         // Record the transitions
        // the position adaptors
        vertices_adaptor<VertexModel::SpaceVec>,    // Record vertex positions
        cells_adaptor<VertexModel::SpaceVec>,       // Record edge information
        edges_adaptor<VertexModel::SpaceVec>,       // Record cell information
        edge_energies_adaptor,                      // Record Edge energies
        cell_energies_adaptor                       // Record Cell energies
    ));
}


int main (int, char** argv) {
    try {
        // Initialize the PseudoParent from config file path
        Utopia::PseudoParent pp(argv[1]);

        auto model = model_factory(pp);
        
        // Register the builder of your work-function
        model.register_work_function_builder(
            "copy_me_const",
            copy_me_const_builder
        );

        // Register the builder of your operation
        model.register_operation_builder(
            "copy_me_operation",
            build_copy_me_operation
        );

        /// Run the model
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
