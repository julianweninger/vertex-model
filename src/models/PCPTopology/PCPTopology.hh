#ifndef UTOPIA_MODELS_PCPTOPOLOGY_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_HH
// TODO Adjust above include guard (and at bottom of file)

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
/** TODO Add your class description here.
 *  ...
 */
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
    PCPVertex<periodic_bc> _vertex_model;
    
    std::shared_ptr<VertexContainer> _vertices;
    std::shared_ptr<EdgeContainer> _edges;
    std::shared_ptr<CellContainer> _cells;

    double _equilibration_tolerance;
    int _num_equilibration_steps;
    int _num_equilibration_iterations;

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

        _vertices(nullptr),
        _edges(nullptr),
        _cells(nullptr),
        // _vertices(_vertex_model.get_vertices_ptr()),
        // _edges(_vertex_model.get_edges_ptr()),
        // _cells(_vertex_model.get_cells_ptr()),

        _equilibration_tolerance(get_as<double>("equilibration_tolerance",
                                             this->_cfg)),
        _num_equilibration_steps(get_as<int>("num_equilibration_steps",
                                             this->_cfg)),
        _num_equilibration_iterations(get_as<int>("num_equilibration_iterations",
                                             this->_cfg))

        // Open the datasets
        // e.g. via _dset_state(this->create_dset("state", {})) <- 1d
        //      or  _dset_state(this->create_dset("state", {num_states})) <- 2d
    {
        if (_vertex_model.get_write_start() == 0) {
            _vertex_model.write_data();
            // NOTE Need this because env model is never run, but only iterated
        }

        this->equilibrate_vertex_model();

        this->_log->debug("Initialised model.");
    }


private:
    // .. Setup functions .....................................................

    // .. Helper functions ....................................................
    void equilibrate_vertex_model() {
        bool equilibrated = false;
        int time_start = _vertex_model.get_time();
        while (not equilibrated) {
            this->_log->debug("Iterating vertex model for equilibration from "
                             "time {} to {}",
                            _vertex_model.get_time(),
                            _vertex_model.get_time() + _num_equilibration_steps);
            for (int i = 0; i < _num_equilibration_steps; ++i) {
                _vertex_model.iterate();
            }
            equilibrated = _vertex_model.equilibrium_state_reached(
                                            _equilibration_tolerance);
            int max_steps = _num_equilibration_steps * _num_equilibration_iterations;
            if (_vertex_model.get_time() - time_start >= max_steps) {
                throw std::runtime_error("Equilibration not reached within {} "
                        "iterations of {} steps each at a tolerance of {}!",
                        _num_equilibration_iterations,
                        _num_equilibration_steps,
                        _equilibration_tolerance);
            }
        }
        this->_log->debug("Vertex model equilibrated within {} steps", 
                         _vertex_model.get_time() - time_start);
    }

public:
    // -- Public Interface ----------------------------------------------------
    // .. Simulation Control ..................................................

    /// Iterate a single step
    void perform_step () {
        // _vertex_model.iterate();
    }


    /// Monitor model information
    /** \details Here, functions and values can be supplied to the monitor that
     *          are then available to the frontend. The monitor() function is
     *          _only_ called if a certain emit interval has passed; thus, the
     *          performance hit is small.
     */
    void monitor () {
        // Can supply information to the monitor here in two ways:
        // this->_monitor.set_entry("key", value);
        // this->_monitor.set_entry("key", [this](){return 42.;});
    }


    /// Write data
    /** \details This function is called to write out data. It should be called
      *         at the end of the model constructor to write out the initial
      *         state. After that, the configuration determines at which times
      *         data is written.
      *         See \ref Utopia::DataIO::Dataset::write
      */
    void write_data () {
        // Example:
        // _dset_foo->write(it.begin(), it.end(),
        //     [](const auto& element) {
        //         return element.get_value();
        // });
    }


    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

};

} // namespace PCPTopology
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPTOPOLOGY_HH
