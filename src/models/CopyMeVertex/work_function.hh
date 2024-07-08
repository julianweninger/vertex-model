#ifndef UTOPIA_MODELS_COPYMEVERTEX_WORKFUNCTION_HH
#define UTOPIA_MODELS_COPYMEVERTEX_WORKFUNCTION_HH


namespace Utopia {
namespace Models {
namespace CopyMeVertex {
namespace WorkFunction {

/// The type of the vertex model
using VertexModel = PCPVertex::PCPVertex;

/// The type of the work-function term of the vertex model
using WFTerm = typename PCPVertex::WorkFunction::WorkFunctionTerm<VertexModel>;

/// The type of a builder to a the work-function term
/** This type can be registered with the Vertex model (or Topology model) */
using WFTermBuilder = typename VertexModel::WFTermBuilder;


/// @brief An example term that adds a constant to the work-function
/** This term adds a constant value to the work-function and a 
 *  gaussian-distributed (constant) value for every cell.
 * 
 *  TODO copy this class and rename to your custom type of a WFTerm
 * 
 *  Parameters:
 *      - `constant`: the constant
 *      - `cell_parameter_mean`:    The mean of cell parameter value
 *      - `cell_parameter_std`:    The std. dev. of cell parameter value
 */
class CopyMeConst 
: 
    public WFTerm
{
    /// The type of a WFTerm the CopyMeConst is inherited from
    using Base = WFTerm;

    /// The type of Agent manager
    /** Here the Vertices, Edges, and Cells are stored, 
     *  as well as how to calculate their properties, e.g.
     *  - position of vertices
     *  - distance between vertices
     *  - length of edges
     *  - area or perimeter cells
     *  - adjoint cells to edges
     *  - neighbors to a cell
     */ 
    using AgentManager = typename Base::AgentManager;

    /// The type of {x, y} coordinates
    using SpaceVec = typename Base::SpaceVec;

    /// The type of a vertex
    using Vertex = typename Base::Vertex;

    /// The type of an edge
    using Edge = typename Base::Edge;

    /// The type of a cell
    using Cell = typename Base::Cell;

// --- The parameters needed to define this WFTerm ----------------------------
private:
    /// @brief  The constant parameter
    double _constant;

    /// @brief  The name under which to register the cell's parameter value
    std::string _cell_constant;

    // TODO replace with your custom parameters

    // --- Inherited from WorkFunctionTerm
    // const std::string _name;     // the name of the workfunction
    // const auto& _am;             // The entities manager containing edges, cells, and vertices


// --- The needed interface according to WFTerm -------------------------------
public:
    /// @brief  CopyMeConst Constructor
    /** How to initialise parameters, etc.
     */
    /// @param name         The name of the term
    /// @param params       The dictionary containing parameters 
    /// @param vertex_model The model from which to inherit the entities manager
    CopyMeConst (
        std::string name,
        const Config& params,
        const VertexModel& vertex_model
    )
    :
        Base(name, params, vertex_model),
        _constant(get_as<double>("constant", params)), 
        _cell_constant(this->_name + "__cell_constant")
        // TODO initialise your custom parameters from the params
    {
        // TODO register custom parameters that are unique to every cell, edge, 
        //      or vertex

        // register the cell constants drawn from a random distribution
        double mean = get_as<double>("cell_parameter_mean", params);
        double std = get_as<double>("cell_parameter_std", params);
        std::normal_distribution<double> distr(mean, std);
        
        // get the container of cells from the agent manager 
        // and iterate all cells
        for (const auto& cell : this->_am.cells()) {
            // register it with the state of the cell
            cell->state.register_parameter(
                _cell_constant,         // the name under which the parameter 
                                        // is registered
                distr(*vertex_model.get_rng()) // the value of the parameter
            );
        }
    }

    /// @brief  The instructions on how to compute the forces
    /** Forces need to be registered using  
     *  vertex->state.add_force(SpaceVec({Fx, Fy}))
     */
    void compute_and_set_forces () final {
        // TODO Define what forces arise from your WFTerm

        for (const auto& vertex : this->_am.vertices()) {
            vertex->state.add_force(SpaceVec({0., 0.})); // zero derivative
        }
    }

    // It is possible to define pressure and tensions here

    // It is possible to define instructions on how to compute energy for a
    // vertex or edge here

    /// How to compute the energy for a single cell
    /** Note, needs to be included in compute_energy(..., cells) to be 
     *  effective.
     */
    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        // TODO (optional) Define what energy each cell contributes
        //                 Can be done analogous for vertices and edges

        return cell->state.get_parameter(_cell_constant);
    }

    /// How compute the total energy for a set of vertices, edges, and cells
    double compute_energy (
        [[maybe_unused]] const Utopia::AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const Utopia::AgentContainer<Edge>& edges,
        const Utopia::AgentContainer<Cell>& cells
    ) const final
    {
        // TODO Define how to compute the total energy

        double energy = _constant; // the constant term
        for (const auto& cell : cells) {
            // accumulate from all cells
            energy += compute_energy(cell);
        }
        return energy;
    }

    /// @brief  How to update the parameter values of this WorkFunctionTerm
    /// @param params   The configuration
    void update_parameters (const Config& params) final {
        // TODO define how to change parameter values
        //      See PCPTopology::OperationCollection
        //      ::build_update_work_function_term in PCPTopology/operations.hh

        // change the value of constant if provided
        _constant  = get_as<double>("constant", params, _constant);

        // changing the cell's parameter is more involved and not implemented!
        if (params["cell_parameter_mean"] or params["cell_parameter_std"]) {
            throw std::runtime_error("No instructions implemented on how to "
                "change the mean and std of cell parameter values in "
                "CopyMeVertex::WorkFunction::CopyMeVertex"
            );
        }
    }

    // Data writers can be added for cell, edge, and vertex properties and
    // energies
};

/// The builder wrapping the CopyMeConst WFTerm
/** This builder can be registered with a PCPTopology instance:
 *  PCPTopology my_model (...);
 *  my_model.register_work_function_builder(
 *      "copy_me_const",        // a unique name to this WFTerm
 *      copy_me_const_builder   // the WFTermBuilder
 *  );
 */
WFTermBuilder copy_me_const_builder = [](
        std::string name, const Config& params,
        const VertexModel& model
) 
{
    return std::make_shared<CopyMeConst>(name, params, model);
    // TODO change CopyMeConst to the type of your custom WFTerm
};


} // namespace WorkFunction
} // namespace CopyMeVertex
} // namespace Models
} // namespace Utopia

#endif
