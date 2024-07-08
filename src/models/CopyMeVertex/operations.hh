#ifndef UTOPIA_MODELS_COPYMEVERTEX_OPERATIONS_HH
#define UTOPIA_MODELS_COPYMEVERTEX_OPERATIONS_HH

#include <typeinfo>
#include <algorithm>
#include <iterator>

#include <utopia/core/types.hh>

#include "../PCPVertex/utils.hh"
#include "../PCPVertex/PCPVertex.hh"
#include "../PCPVertex/minimization.hh"

#include "../PCPTopology/operations.hh"

namespace Utopia {
namespace Models {
namespace CopyMeVertex {
namespace OperationCollection {

using namespace Utopia::Models::PCPVertex;
using namespace Utopia::Models::PCPTopology::OperationCollection;

using VertexModel = Utopia::Models::PCPVertex::PCPVertex;

/// The type of builder that contains an OperationBundel
using OperationBundleBuilder = std::function<OperationBundle(
    std::string, const Config&, const MinimizationParams&
)>;

/// A sample operation writing to the logger and dividing a cell
OperationBundleBuilder build_copy_me_operation = [](
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params
)
{
    // The instruction when to apply the operation
    // and how to minimize the workfunction thereafter
    OperationParams params(name, cfg, default_minim_params);


    // TODO extract your parameters here from cfg
    double my_parameter = get_as<double>("my_parameter", cfg);
    auto my_vector = get_as< std::vector<int> >("my_vector", cfg);
    // auto -> compiler translated to std::vector<int>

    // The actual operation instruction
    // capture my_parameter and my_vector
    Operation operation = std::function<void(VertexModel&)>(
        [my_parameter, my_vector]
        (VertexModel& vertex_model)
        {
        // TODO write your own operation

        // write something to the logger of the vertex model
        vertex_model.get_logger()->info("This is operation foo.");
        vertex_model.get_logger()->info("This operation has parameter {}",
                                        my_parameter);
        
        // obtain the entities manager from vertex model
        const auto& am = vertex_model.get_am();

        // obtain the agents of the vertex model managed by agent-manager
        [[maybe_unused]] const auto& vertices = am.vertices();
        [[maybe_unused]] const auto& edges = am.edges();
        [[maybe_unused]] const auto& cells = am.cells();

        // count the cells
        vertex_model.get_logger()->info("The Vertex model has {} cells.",
                                        cells.size());

        // select the first cell
        const auto& cell = cells[0];

        // register a parameter with that cell
        vertex_model.get_logger()->debug(
            "Registering my parameter with the cell {}", cell->id()
        );
        if (not cell->state.has_parameter("my_parameter")) {
            cell->state.register_parameter("my_parameter", my_parameter);
        }

        // increase the vertex model's domain size in anticipiation of a
        // cell division
        vertex_model.get_logger()->debug("Increasing domain size by new cell's "
                                         "area={}",
                                         cell->state.area_preferential);
        vertex_model.increase_domain_size(cell->state.area_preferential);

        // Divide a cell at angle 0.
        vertex_model.get_logger()->debug("Dividing cell {}", cell->id());
        vertex_model.divide_cell(cell, 0.);

        // Count the number of cells
        vertex_model.get_logger()->info("The Vertex model has {} cells.",
                                        cells.size());
        }
    );

    return std::make_pair(operation, params);
};

} // namespace OperationCollection
} // namespace CopyMeVertex
} // namespace Models
} // namespace Utopia
#endif