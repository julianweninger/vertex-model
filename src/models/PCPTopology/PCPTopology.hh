#ifndef UTOPIA_MODELS_PCPTOPOLOGY_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>
#include <utopia/data_io/data_manager/defaults.hh>

// model for energy minimization
#include "../PCPVertex/PCPVertex.hh"
#include "../PCPVertex/space.hh"
#include "../PCPVertex/energy.hh"
#include "../PCPVertex/algorithm.hh"
#include "../PCPVertex/operations.hh"
#include "../PCPVertex/PCPVertex_write_tasks.hh"

#include "PCPTopology_write_tasks.hh"
#include "operations.hh"

// coupled models
#include "../NotchDelta/NotchDelta.hh"
#include "../Collier/Collier.hh"
#include "../NotchDelta/Differentiation_write_tasks.hh"
#include "../Collier/Collier_write_tasks.hh"

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

    /// The type of a vertex
    using Vertex = typename PCPVertex::Vertex;

    /// The type of an edge
    using Edge = typename PCPVertex::Edge;

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

    /// The instance of the Collier model used by proliferation tasks
    /** \note Only initialized when needed
     */
    std::shared_ptr<Collier::Collier> _collier;

    /// Whether the _notch_delta model's prolog was performed
    std::shared_ptr<bool> _collier_prolog;

    /// The instance of the notch_delta model used by proliferation tasks
    /** \note Only initialized when needed
     */
    std::shared_ptr<NotchDelta::NotchDelta> _notch_delta;

    /// Whether the _notch_delta model's prolog was performed
    std::shared_ptr<bool> _notch_delta_prolog;
    
    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::uniform_real_distribution<double> _prob_distr;

    // .. Temporary objects ...................................................
    /// The number of T1 transitions
    std::size_t _num_T1s;

    /// The total number of T1 transitions
    std::size_t _num_T1s_total;

    /// The number of T1 transitions attempted
    std::size_t _num_T1s_attempted;

    /// The total number of T1 transitions attempted
    std::size_t _num_T1s_attempted_total;

    /// The number of T2 transitions
    std::size_t _num_T2s;

    /// The total number of T2 transitions
    std::size_t _num_T2s_total;

    std::size_t _estimate_minimizations;

public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the PCPTopology model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel, typename... WriterArgs>
    PCPTopology (const std::string name, ParentModel &parent_model,
                 const Config& custom_cfg = {},
                 std::tuple<WriterArgs...> &&writer_args = {})
    :
        // Initialize first via base model
        Base(name, parent_model, custom_cfg, writer_args),
        
        // construct the vertex model with an external maximum time stamp
        _vertex_model("PCPVertex", *this, {},
            std::make_tuple(
                // energy adaptors
                DataIO::time_energy_adaptor, DataIO::energy_adaptor,
                DataIO::areaelasticity_adaptor,
                DataIO::linetension_adaptor,
                DataIO::cell_contractility_adaptor,
                DataIO::edge_contractility_adaptor,
                DataIO::cell_cell_polarity_adaptor,
                DataIO::polarity_exclusion_adaptor,
                DataIO::lagrange_net_polarisation_adaptor,
                DataIO::lagrange_const_concentration_adaptor,
                // transition adaptors
                DataIO::statistics_time_adaptor,
                DataIO::T1_adaptor, DataIO::T1_attempted_adaptor,
                DataIO::T2_adaptor,
                // position adaptors
                DataIO::vertices_adaptor<SpaceVec>,
                DataIO::cells_adaptor<SpaceVec, CellType>,
                DataIO::edges_adaptor,
                DataIO::cell_energies_adaptor,
                DataIO::edge_energies_adaptor
                ),
            DataIO::build_custom_deciders<PCPVertex>()),
        
        // the parameter
        _minimization_params(get_as<Config>("minimization", this->_cfg)),
        _operations{},
        _prob_distr(0.,1.),
        _num_T1s(0),
        _num_T1s_total(0),
        _num_T1s_attempted(0),
        _num_T1s_attempted_total(0),
        _num_T2s(0),
        _num_T2s_total(0),
        _estimate_minimizations(0)
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

                if (name == "differentiate_Collier") {
                    this->setup_collier(
                            get_as<Config>("Collier", op_cfg, {}));
                    _operations.push_back(
                        build_differentiate_Collier(name, op_cfg,
                            _minimization_params, _collier,
                            _collier_prolog));
                }
                else if (name == "differentiate_NotchDelta") {
                    this->setup_notch_delta(
                            get_as<Config>("NotchDelta", op_cfg, {}));
                    _operations.push_back(
                        build_differentiate_NotchDelta(name, op_cfg,
                            _minimization_params, _notch_delta,
                            _notch_delta_prolog));
                }
                else if (name == "differentiate_random") {
                    _operations.push_back(
                        build_differentiate_random(name, op_cfg,
                                                   _minimization_params));
                }
                else if (name == "differentiate_hair_cluster") {
                    _operations.push_back(
                        build_differentiate_hair_cluster(name, op_cfg,
                                                         _minimization_params));
                }
                else if (name == "increment_area") {
                    _operations.push_back(
                        build_increment_area(name, op_cfg,
                                             _minimization_params));
                }
                else if (name == "increment_cell_contractility") {
                    _operations.push_back(
                        build_increment_cell_contractility(name, op_cfg,
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
                        build_proliferate(name, op_cfg, _minimization_params));
                }
                else if (name == "relax_area") {
                    _operations.push_back(
                        build_relax_area(name, op_cfg, _minimization_params));
                }
                else if (name == "set_area") {
                    _operations.push_back(
                        build_set_area(name, op_cfg,
                                       _minimization_params));
                }
                else {
                    throw std::invalid_argument(fmt::format(
                        "No operation '{}' available to construct! "
                        "Choose from: {}", name,
                            "differentiate_Collier, "
                            "differentiate_NotchDelta, "
                            "differentiate_random, "
                            "differentiate_hair_cluster, "
                            "increment_area, "
                            "increment_cell_contractility, "
                            "increment_domain, "
                            "increment_edge_contractility, "
                            "increment_linetension, "
                            "increment_shape_index, "
                            "jiggle, "
                            "proliferate, "
                            "relax_area, "
                            "set_area, "
                            "void."));
                }

                this->_log->debug("Added '{}' operation.", name);
                const auto& params = std::get<1>(_operations.back());
                _estimate_minimizations += params.get_num_minimizations(
                                                        this->_time_max);
            }
        }
    }
    
    /// Setup a Collier model
    void setup_collier (const Config& cfg = {})
    {
        if (_collier) {
            return;
        }

        this->_log->debug("Setting up Collier model from {}",
            cfg.size() ? 
                "custom configuration."
                : 
                fmt::format("configuration within {} model.", this->_name));

        _collier = std::shared_ptr<Collier::Collier>(
            new Collier::Collier("Collier", *this, cfg, 
                std::make_tuple(
                    Differentiation::DataIO::density_time,
                    Differentiation::DataIO::density_progenitor,
                    Differentiation::DataIO::density_hair,
                    Differentiation::DataIO::density_support,
                    Differentiation::DataIO::density_rosettes,
                    Collier::DataIO::density_notch,
                    Collier::DataIO::density_delta,
                    Collier::DataIO::density_nicd
                )
            )
        );

        _collier_prolog = std::make_shared<bool>(false);        
    }
    
    /// Setup a notch delta model
    void setup_notch_delta (const Config& cfg = {})
    {
        if (_notch_delta) {
            return;
        }

        this->_log->debug("Setting up NotchDelta model from {}",
            cfg.size() ? 
                "custom configuration."
                : 
                fmt::format("configuration within {} model.", this->_name));

        _notch_delta = std::shared_ptr<NotchDelta::NotchDelta>(
            new NotchDelta::NotchDelta("NotchDelta", *this, cfg, 
                std::make_tuple(
                    Differentiation::DataIO::density_time,
                    Differentiation::DataIO::density_progenitor,
                    Differentiation::DataIO::density_hair,
                    Differentiation::DataIO::density_support,
                    Differentiation::DataIO::density_rosettes
                )
            )
        );

        _notch_delta_prolog = std::make_shared<bool>(false);        
    }
    
    // .. Helper functions ....................................................
    /** Apply an operation from a bundle
     *  \param operation_bundle collection of operation and parameters
     *  \param prolog           Whether it is called during the prolog of the
     *                          model. `iterations_prolog` are performed.
     *  \param epilog   	    Whether it is called during the epilog of the
     *                          model. `iterations_epilog` are performed.
     *                          Whenever the energy is minimized the time of the
     *                          model is incremented and the datamanager is
     *                          called.
     */
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
                 *params.times.begin() == this->_time + 1)
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

        std::function<void()> monitor_mngr = [this](){
            this->_monitor.get_monitor_manager()->check_timer();
            this->__monitor();
            this->_monitor.get_monitor_manager()->emit_if_enabled();
        };

        for (std::size_t it = 0; it < iterates; it++) {
            if (not params.disable) {
                operation(_vertex_model);
            }
            // else: skip 

            if (params.minimization_mode == MinimizationMode::Every) {
                this->_log->debug("   Minimizing energy ...");
                _vertex_model.minimize_energy(params.minimization_params,
                                              monitor_mngr);
                
                // write data during epilog
                if (epilog) {
                    this->increment_time();
                    this->_datamanager(static_cast<PCPTopology&>(*this));
                }
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
            _vertex_model.minimize_energy(params.minimization_params,
                                          monitor_mngr);
                
            // write data during epilog
            if (epilog) {
                this->increment_time();
                this->_datamanager(static_cast<PCPTopology&>(*this));
            }
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

        _num_T1s = _vertex_model.get_num_T1s_total() - _num_T1s_total;
        _num_T1s_attempted = _vertex_model.get_num_T1s_attempted_total() - 
                             _num_T1s_attempted_total;
        _num_T2s = _vertex_model.get_num_T2s_total() - _num_T2s_total;

        _num_T1s_total = _vertex_model.get_num_T1s_total();
        _num_T1s_attempted_total = _vertex_model.get_num_T1s_attempted_total();
        _num_T2s_total = _vertex_model.get_num_T2s_total();
    }

    /// Monitor model information
    void monitor () {
        // overwrite the progress with progress estimate
        this->_monitor.get_monitor_manager()->set_time_entries(
            _vertex_model.get_num_minimizations(),
            std::max(_estimate_minimizations,
                     _vertex_model.get_num_minimizations() + 1));

        this->_monitor.set_entry("num_cells",
                                 _vertex_model.get_am().cells().size());
        this->_monitor.set_entry("num_T1_transitions",
                                 _vertex_model.get_num_T1s_total());
        this->_monitor.set_entry("num_T1_transitions_attempted",
                                 _vertex_model.get_num_T1s_attempted_total());
        this->_monitor.set_entry("num_T2_transitions",
                                 _vertex_model.get_num_T2s_total());

        this->_monitor.set_entry("time", this->get_time());
        this->_monitor.set_entry("progress",   float(this->get_time())
                                             / float(this->get_time_max()));

        _vertex_model.monitor();
    }

    /// The prolog
    /** Performs the following tasks:
     *      1. prolog of vertex model
     *      1. prolog operations
     *      1. default prolog tasks
     */
    void prolog () {
        _vertex_model.prolog();

        this->_log->info("Running prolog operations ...");

        for (auto& operation : _operations) {
            apply_operation(operation, true, false);
        }

        _num_T1s_total = _vertex_model.get_num_T1s_total();
        _num_T1s_attempted_total = _vertex_model.get_num_T1s_attempted_total();
        _num_T2s_total = _vertex_model.get_num_T2s_total();
        _num_T1s = _num_T1s_total;
        _num_T1s_attempted = _num_T1s_attempted_total;
        _num_T2s = _num_T2s_total;
        
        return this->__prolog();
    }

    /// The epilog
    /** Performs the following tasks:
     *      1. epilog operations
     *      1. default epilog tasks
     */
    void epilog () {
        for (auto& operation : _operations) {
            apply_operation(operation, false, true);
        }

        _num_T1s = _vertex_model.get_num_T1s_total() - _num_T1s_total;
        _num_T1s_attempted = _vertex_model.get_num_T1s_attempted_total() - 
                             _num_T1s_attempted_total;
        _num_T2s = _vertex_model.get_num_T2s_total() - _num_T2s_total;

        _num_T1s_total = _vertex_model.get_num_T1s_total();
        _num_T1s_attempted_total = _vertex_model.get_num_T1s_attempted_total();
        _num_T2s_total = _vertex_model.get_num_T2s_total();

        _vertex_model.epilog();

        if (_notch_delta) {
            _notch_delta->epilog();
        }

        return this->__epilog();
    }

    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model
    /// Getter for the continuous time
    /** \details This is the time of the vertex model used for energy
     *           minimization.
     */
    std::size_t get_continuous_time() const {
        return _vertex_model.get_time();
    }
    /// Getter for the total energy
    /** See PCPVertex::get_energy
     */
    double get_energy() const {
        return _vertex_model.get_energy();
    }
    /// Getter for the linetension energy
    /** See PCPVertex::get_energy_linetension
     */
    double get_energy_linetension() const {
        return _vertex_model.get_energy_linetension();
    }
    /// Getter for the edge contractility energy
    /** See PCPVertex::get_energy_edge_contractility
     */
    double get_energy_edge_contractility() const {
        return _vertex_model.get_energy_edge_contractility();
    }
    /// Getter for the areaelasticity energy
    /** See PCPVertex::get_energy_areaelasticity
     */
    double get_energy_areaelasticity () const {
        return _vertex_model.get_energy_areaelasticity();
    }
    /// Getter for the areaelasticity energy
    /** See PCPVertex::get_energy_cell_contractility
     */
    double get_energy_cell_contractility () const {
        return _vertex_model.get_energy_cell_contractility();
    }

    /// Getter for the linetension energy on a set of entities
    /** See PCPVertex::get_energy_linetension(const AgentContainer<Edge>& es) const
     */
    double get_energy_linetension(const AgentContainer<Edge>& es) const {
        return _vertex_model.get_energy_linetension(es);
    }

    /// Getter for the edge contractility energy on a set of entities
    /** See PCPVertex::get_energy_edge_contractility(const AgentContainer<Edge>& es) const
     */
    double get_energy_edge_contractility(const AgentContainer<Edge>& es) const {
        return _vertex_model.get_energy_edge_contractility(es);
    }

    /// Getter for the areaelasticity energy on a set of entities
    /** See PCPVertex::get_energy_areaelasticity(const AgentContainer<Cell>& cs) const
     */
    double get_energy_areaelasticity (const AgentContainer<Cell>& cs) const {
        return _vertex_model.get_energy_areaelasticity(cs);
    }

    /// Getter for the areaelasticity energy on a set of entities
    /** See PCPVertex::get_energy_cell_contractility(const AgentContainer<Cell>& cs) const
     */
    double get_energy_cell_contractility (const AgentContainer<Cell>& cs) const
    {
        return _vertex_model.get_energy_cell_contractility(cs);
    }
    
    // double get_energy_cell_cell_polarity() const {
    //     return _vertex_model.get_energy_cell_cell_polarity();
    // }
    // double get_energy_polarity_exclusion () const {
    //     return _vertex_model.get_energy_polarity_exclusion();
    // }
    // double get_energy_lagrange_net_polarisation() const {
    //     return _vertex_model.get_energy_lagrange_net_polarisation();
    // }
    // double get_energy_lagrange_const_concentration() const {
    //     return _vertex_model.get_energy_lagrange_const_concentration();
    // }

    std::size_t get_num_T1s() const {
        return _num_T1s;
    }

    std::size_t get_num_T1s_attempted() const {
        return _num_T1s_attempted;
    }

    std::size_t get_num_T2s() const {
        return _num_T2s;
    }

    /// Getter for vertices
    const auto& get_am () const {
        return _vertex_model.get_am();
    }
    
    /// Add an operation
    void register_operation (OperationBundle op_bundle) {
        _operations.push_back(op_bundle);
        auto [_, params] = op_bundle;
        this->_log->trace("Registered operation '{}'", params.name);

        _estimate_minimizations += params.get_num_minimizations(
                                                this->_time_max);
    }
};

} // namespace PCPTopology
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPTOPOLOGY_HH
