#ifndef UTOPIA_MODELS_PCPTOPOLOGY_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>

#include "../PCPVertex/PCPVertex.hh"
#include "../PCPVertex/space.hh"
#include "../PCPVertex/energy.hh"
#include "../PCPVertex/algorithm.hh"
#include "../PCPVertex/operations.hh"
#include "../PCPVertex/PCPVertex_write_tasks.hh"

#include "../NotchDelta/NotchDelta.hh"
#include "../NotchDelta/NotchDelta_write_tasks.hh"

#include "PCPTopology_write_tasks.hh"

#include <utopia/models/Environment/Environment.hh>

#include "operations.hh"

namespace Utopia {
namespace Models {
namespace PCPVertex {
// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
using namespace OperationCollection;

/// Type helper to define types used by the model
using PCPTopologyModelTypes = Utopia::ModelTypes<DefaultRNG, WriteMode::managed,
                                                 Space::CustomSpace<2>>;

// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// The PCPTopology Model; the bare-basics a model needs
class PCPTopology:
    public Model<PCPTopology, PCPTopologyModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<PCPTopology, PCPTopologyModelTypes>;

    /// Type of the config
    using typename Base::Config;

    /// The type of a cell
    using Cell = typename PCPVertex::Cell;

    /// The types a cell can take (support, hair)
    using CellType = typename PCPVertex::CellType;

    /// The type of a rule function acting on vertices of the agent manager
    using RuleFuncVertex = typename PCPVertex::RuleFuncVertex;
    
    /// The type of a rule function acting on edges of the agent manager
    using RuleFuncEdge = typename PCPVertex::RuleFuncEdge;
    
    /// The type of a rule function acting on cells of the agent manager
    using RuleFuncCell = typename PCPVertex::RuleFuncCell;

    /// The type of coordinates and vectors in space
    using SpaceVec = typename PCPVertex::SpaceVec;
                                    

private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The Vertex model
    PCPVertex _vertex_model;

    /// The parameter for energy minimization
    MinimizationParams _minimization_params;

    /// Container of operations
    /** The operations are applied in the order of registration at specified
     *  times (see `OperationBundle`) and energy minimization is performed as
     *  defined in `OperationBundle::MinimizationMode`. The parameter for
     *  minimization can be updated from default for every operation.
     *  Operations can duplicate with same or different parameter.
     */
    std::vector<OperationBundle> _operations;
    
    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::uniform_real_distribution<double> _prob_distr;

    // .. Temporary objects ...................................................

public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the PCPTopology model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel, typename... Taskargs>
    PCPTopology (const std::string name, ParentModel& parent, 
                 Taskargs&&... taskargs)
    :
        // Initialize first via base model
        Base(name, parent, std::forward<Taskargs>(taskargs)...),
        
        // construct the vertex model with an external maximum time stamp
        _vertex_model("PCPVertex", *this,
                    DataIO::time_energy_adaptor, DataIO::energy_adaptor,
                    DataIO::areaelasticity_adaptor,
                    DataIO::linetension_adaptor,
                    DataIO::contractility_adaptor,
                    DataIO::cell_cell_polarity_adaptor,
                    DataIO::polarity_exclusion_adaptor,
                    DataIO::lagrange_net_polarisation_adaptor,
                    DataIO::lagrange_const_concentration_adaptor,
                    DataIO::vertices_adaptor,  
                    DataIO::cells_adaptor<SpaceVec, CellType>,
                    DataIO::edges_adaptor),
        
        // the parameter
        _minimization_params(get_as<Config>("minimization", this->_cfg)),
        _operations{},
        _prob_distr(0.,1.)
    {
        this->_space = _vertex_model.get_space();

        setup_operations(get_as<Config>("operations", this->_cfg));

        this->_log->info("Model set up.");
    }


private:
    // .. Setup functions .....................................................
    void setup_operations(const Config& cfg) {
        this->_log->info("Setting up operations from {} configuration entr{} "
                         "...", cfg.size(), cfg.size() != 1 ? "ies" : "y");

        // For zombie or empty configurations, return empty container
        if (not cfg or not cfg.size()) {
            return;
        }
        // Otherwise, require a sequence
        if (not cfg.IsSequence()) {
            throw std::invalid_argument("The config for initializing the "
                "operations must be a sequence!");
        }

        // Iterate over the sequence of mappings
        for (const auto& ops : cfg) {
            // ops.IsMap() == true
            // The top `ops` keys are now the names of the desired environment
            // functions. Iterate over those ...
            for (const auto& op_pair : ops) {
                const auto name = op_pair.first.as<std::string>();
                const auto& op_cfg = op_pair.second;
                this->_log->trace("  Operation name:  {}", name);

                if (name == "differentiate_NotchDelta") {
                    auto notch_delta = std::make_shared<
                        NotchDelta::NotchDelta>(
                            "NotchDelta", *this,
                            NotchDelta::DataIO::density_time,
                            NotchDelta::DataIO::density_progenitor,
                            NotchDelta::DataIO::density_hair,
                            NotchDelta::DataIO::density_support,
                            NotchDelta::DataIO::density_ratio_hair_support,
                            NotchDelta::DataIO::number_hair_hair_contacts);
                    _operations.push_back(
                        build_differentiate_NotchDelta(name, op_cfg,
                            _minimization_params, notch_delta));
                }
                else if (name == "differentiate_random") {
                    _operations.push_back(
                        build_differentiate_random(name, op_cfg,
                                                   _minimization_params));
                }
                else if (name == "increment_area") {
                    _operations.push_back(
                        build_increment_area(name, op_cfg,
                                             _minimization_params));
                }
                else if (name == "increment_domain") {
                    _operations.push_back(
                        build_increment_domain(name, op_cfg,
                                               _minimization_params));
                }
                else if (name == "increment_edge_contractility") {
                    _operations.push_back(
                        build_increment_edge_contractility(name, op_cfg,
                            _minimization_params));
                }
                else if (name == "increment_linetension") {
                    _operations.push_back(
                        build_increment_linetension(name, op_cfg,
                                                    _minimization_params));
                }
                else if (name == "increment_shape_index") {
                    _operations.push_back(
                        build_increment_shape_index(name, op_cfg,
                                                    _minimization_params));
                }
                else if (name == "jiggle" or name == "void") {
                    _operations.push_back(
                        build_jiggle(name, op_cfg, _minimization_params));
                }
                else if (name == "proliferate") {
                    _operations.push_back(
                        build_proliferate(name, op_cfg,
                                            _minimization_params));
                }
                else {
                    throw std::invalid_argument(fmt::format(
                        "No operation '{}' available to construct! "
                        "Choose from: {}", name,
                            "differentiate_NotchDelta, "
                            "differentiate_random, "
                            "increment_area, "
                            "increment_domain, "
                            "increment_edge_contractility, "
                            "increment_linetension, "
                            "increment_shape_index, "
                            "jiggle, "
                            "proliferate, "
                            "void."));
                }

                this->_log->debug("Added '{}' operation.", name);
            }
        }
    }
    
    // .. Helper functions ....................................................
    void apply_operation(OperationBundle& operation_bundle,
                         bool prolog = false, bool epilog = false)
    {
        using MinimizationMode = OperationParams::MinimizationMode;

        auto& [operation, params] = operation_bundle;

        std::size_t iterates;
        if (prolog) {
            iterates = params.iterations_prolog;
        }
        else if (epilog) {
            iterates = params.iterations_epilog;
        }
        else if (params.times.size() and
                 *params.times.begin() == this->_time)
        {
            // Invoke at this time; pop element corresponding to this time
            params.times.erase(params.times.begin());

            if (params.probability < 1 and
                _prob_distr(*this->_rng) > params.probability)
            {
                iterates = 0;
            }
            else {
                iterates = params.iterations;
            }
        }
        else {
            iterates = 0;
        }

        auto emit_interval = params.emit_interval;
        if (iterates == 0) {
            this->_log->trace("Not invoking operation '{}' in iteration {}",
                              params.name, this->_time);
        }
        else if (emit_interval > 0) {
            this->_log->info("Applying operation '{}' with {} iterates ...",
                             params.name, iterates);
        }
        else {
            this->_log->debug("Applying operation '{}' with {} iterates ...",
                              params.name, iterates);
            emit_interval = iterates + 1;
        }

        for (std::size_t it = 0; it < iterates; it++) {
            operation(_vertex_model);

            if (params.minimization_mode == MinimizationMode::Every) {
                this->_log->debug("   Minimizing energy ...");
                _vertex_model.minimize_energy(params.minimization_params);
            }
            else {
                this->_log->debug("   NOT minimizing energy.");
            }

            if ((it+1) % emit_interval == 0) {
                this->_log->info("   Performed iterate {} / {} of operation "
                                 "'{}'.", it+1, iterates, params.name);
            }
            else {
                this->_log->debug("   Performed iterate {} / {} of operation "
                                  "'{}'.", it+1, iterates, params.name);
            }
        }
        if (iterates > 0 and
            params.minimization_mode == MinimizationMode::Once)
        {
            this->_log->debug("   Minimizing energy ...");
            _vertex_model.minimize_energy(params.minimization_params);
        }
        else if (iterates > 0 and
                 params.minimization_mode == MinimizationMode::Manual)
        {
            this->_log->debug("   Energy was NOT minimized!");
        }
    }

public:
    // -- Public Interface ----------------------------------------------------

    // .. Simulation Control ..................................................
    /// Iterate a single step
    void perform_step () {
        for (auto& operation_bundle : _operations) {
            apply_operation(operation_bundle);
        }
    }

    /// Monitor model information
    void monitor () {        
        this->_monitor.set_entry("num cells",
                                 _vertex_model.get_am().cells().size());
    }

    /// The prolog
    /** Performs the following tasks:
     *      1. initialize vertex model
     *      2. equilibrate vertex model
     *      3. prolog operations
     *      4. default prolog tasks
     */
    void prolog () {
        _vertex_model.prolog();
        _vertex_model.minimize_energy(_minimization_params);

        for (auto& operation : _operations) {
            apply_operation(operation, true, false);
        }
        
        return this->__prolog();
    }

    /// The epilog
    /** Performs the following tasks:
     *      1. epilog operations
     *      2. default epilog tasks
     */
    void epilog () {
        for (auto& operation : _operations) {
            apply_operation(operation, false, true);
        }

        _vertex_model.epilog();

        return this->__epilog();
    }

    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

    /// Getter for vertices
    const auto& get_am () const {
        return _vertex_model.get_am();
    }
    
    /// Add an operation
    void register_operation (OperationBundle& operation) {
        _operations.push_back(operation);
    }
};

} // namespace PCPTopology
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPTOPOLOGY_HH
