#ifndef UTOPIA_MODELS_PCPTOPOLOGY_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>

#include "../PCPVertex/geometry.hh"
#include "../PCPVertex/PCPVertex.hh"


namespace Utopia {
namespace Models {
namespace PCPVertex {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<>;

// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// The PCPTopology Model; the bare-basics a model needs
template <bool periodic_bc>
class PCPTopology:
    public Model<PCPTopology<periodic_bc>, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<PCPTopology<periodic_bc>, ModelTypes>;

    /// Data type of the group to write model data to, holding datasets
    using DataGroup = typename Base::DataGroup;

    /// Data type for a dataset
    using DataSet = typename Base::DataSet;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The Vertex model
    PCPVertex<periodic_bc> _vertex_model;
    
    /// VertexContainer
    /// NOTE insert and delete in VertexModel not updated automatically
    std::vector<std::weak_ptr<Vertex>> _vertices;
    
    /// EdgeContainer
    /// NOTE insert and delete in VertexModel not updated automatically    
    std::vector<std::weak_ptr<Edge>> _edges;

    /// CellContainer
    /// NOTE insert and delete in VertexModel not updated automatically
    std::vector<std::weak_ptr<Cell>> _cells;

    /// A tolerance value for equilibrium
    double _equilibration_tolerance;

    /// Number of steps performed in VertexModel per iteration
    int _num_equilibration_steps;

    /// Number of max iterations performed in VertexModel before aborting
    int _max_equilibration_iterations;

    double _cell_divisions_per_step;
    
    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::uniform_real_distribution<double> _prob_distr;

    // .. Temporary objects ...................................................


    // .. Datasets ............................................................
    // NOTE They should be named '_dset_<name>', where <name> is the
    //      dataset's actual name as set in its constructor. Ideally, do not
    //      hide them inside a struct ...
    // std::shared_ptr<DataSet> _dset_my_var;


public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the PCPTopology model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel>
    PCPTopology (const std::string name, ParentModel& parent)
    :
        // Initialize first via base model
        Base(name, parent),
        
        _vertex_model("PCPVertex", *this),

        _vertices(_vertex_model.get_vertices()),
        _edges(_vertex_model.get_edges()),
        _cells(_vertex_model.get_cells()),

        _equilibration_tolerance(get_as<double>("equilibration_tolerance",
                                             this->_cfg)),
        _num_equilibration_steps(get_as<int>("num_equilibration_steps",
                                             this->_cfg)),
        _max_equilibration_iterations(get_as<int>("max_equilibration_iterations",
                                                  this->_cfg)),
        _cell_divisions_per_step(get_as<double>("cell_divisions_per_step",
                                                this->_cfg)),
        _prob_distr(0.,1.)
    {
        if (_vertex_model.get_write_start() == 0) {
            _vertex_model.write_data();
            // NOTE Need this because env model is never run, but only iterated
        }

        this->equilibrate_vertex_model();

        // update the vertices, edges, and cell container
        this->update_object_containers();

        this->_log->debug("Initialised model.");
    }


private:
    // .. Setup functions .....................................................

    // .. Helper functions ....................................................
    /// Updates the Containers of vertices, edges, and cells
    /** NOTE insert and deletion opterations to these containers are not 
     *       updated automatically for changes within the vertex model
    */
    void update_object_containers() {
        _vertices = _vertex_model.get_vertices();
        _edges = _vertex_model.get_edges();
        _cells = _vertex_model.get_cells();
    }

    /// Equilibrates the vertex model
    /** Iterate the vertex model until it reaches an equilibrium state
     *  (see PCPVertex<periodic_bc>::equilibrium_state_reached(double tolerance) const)
     *  with tolerance = _equilibrium_tolerance
     * 
     *  More precisely, the _vertex_model is iterated for 
     *  _num_equilibration_steps, then the equilibrium condition (see above)
     *  is checked:
     * 
     *      If true, iteration is stopped and update_object_containers() called
     * 
     *      If false, _vertex_model is iterated again for _num_equilibration_steps 
     *          this is repeated for a maximum of _max_equilibration_iterations
     */
    void equilibrate_vertex_model() {
        bool equilibrated = false;
        int time_start = _vertex_model.get_time();
        while (not equilibrated) {
            this->_log->debug("Iterating vertex model for equilibration from "
                "time {} to {}", _vertex_model.get_time(),
                _vertex_model.get_time() + _num_equilibration_steps);
            for (int i = 0; i < _num_equilibration_steps; ++i) {
                _vertex_model.iterate();
            }
            equilibrated = _vertex_model.equilibrium_state_reached(
                                            _equilibration_tolerance);
            int max_steps = _num_equilibration_steps * _max_equilibration_iterations;
            if (_vertex_model.get_time() - time_start >= max_steps) {
                throw std::runtime_error("Equilibrium not reached within " +
                        std::to_string(_max_equilibration_iterations) + 
                        " iterations of " + 
                        std::to_string(_num_equilibration_steps) + " steps each "
                        "at a tolerance of " +
                        std::to_string(_equilibration_tolerance) + "!");
            }
        }

        // update the vertices, edges, and cell container
        this->update_object_containers();

        this->_log->debug("Vertex model equilibrated within {} steps", 
                         _vertex_model.get_time() - time_start);
    }

public:
    // -- Public Interface ----------------------------------------------------
    // .. Simulation Control ..................................................
    /// Perform N cell divisions
    /** \param num_cell_divisions number of cell divisions to be performed
     * 
     */
   void perform_cell_divisions(double num_cell_divisions) {
        int i;
        for (i = 1; i <= num_cell_divisions; ++i) {
            _vertex_model.divide_random_cell();
        }
        i -= 1;
        if (num_cell_divisions - i > 0 and 
                _prob_distr(*this->_rng) < num_cell_divisions - i)
        {
            _vertex_model.divide_random_cell();
        }

        equilibrate_vertex_model();
    }

    /// Iterate a single step
    void perform_step () {
        perform_cell_divisions(_cell_divisions_per_step);
    }


    /// Monitor model information
    void monitor () { }


    /// Write data
    void write_data () { }


    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

};

} // namespace PCPTopology
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPTOPOLOGY_HH
