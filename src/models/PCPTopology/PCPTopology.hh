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
#include "../PCPVertex/work_function.hh"
#include "../PCPVertex/algorithm.hh"
#include "../PCPVertex/operations.hh"
#include "../PCPVertex/PCPVertex_write_tasks.hh"

#include "PCPTopology_write_tasks.hh"

#include "work_function.hh"

#include "operations.hh"
#include "operations_differentiation.hh"

// coupled models
#include "../NotchDelta/NotchDelta.hh"
#include "../NotchDelta/Differentiation_write_tasks.hh"

#include "../Collier/Collier.hh"
#include "../Collier/Collier_write_tasks.hh"

#include "../PlanarCellPolarity/PlanarCellPolarity.hh"
#include "../PlanarCellPolarity/PlanarCellPolarity_write_tasks.hh"



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

    using AgentManager = typename PCPVertex::AgentManager;

    /// The type of a vertex
    using Vertex = typename PCPVertex::Vertex;

    /// The type of an edge
    using Edge = typename PCPVertex::Edge;

    /// The type of a cell
    using Cell = typename PCPVertex::Cell;

    /// The type of a rule function acting on vertices of the agent manager
    using RuleFuncVertex = typename PCPVertex::RuleFuncVertex;
    
    /// The type of a rule function acting on edges of the agent manager
    using RuleFuncEdge = typename PCPVertex::RuleFuncEdge;
    
    /// The type of a rule function acting on cells of the agent manager
    using RuleFuncCell = typename PCPVertex::RuleFuncCell;

    /// The type of coordinates and vectors in space
    using SpaceVec = typename PCPVertex::SpaceVec;

    /// The type of Proteins
    using ProteinVec = typename PlanarCellPolarity::\
                                PlanarCellPolarity<PCPVertex::AgentManager>::\
                                ProteinVec;
                                    

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

    /// The instance of the notch_delta model used by proliferation tasks
    /** \note Only initialized when needed
     */
    std::shared_ptr<
        PlanarCellPolarity::PlanarCellPolarity<PCPVertex::AgentManager>> _pcp;

    /// Whether the _notch_delta model's prolog was performed
    std::shared_ptr<bool> _pcp_prolog;
    
    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::uniform_real_distribution<double> _prob_distr;

    // .. Temporary objects ...................................................
    /// The number of T1 transitions
    std::size_t _num_T1s;

    /// The total number of T1 transitions
    std::size_t _num_T1s_total;

    /// The frequency of T1 transitions per edge since initialisation
    double _T1_frequency_acc;

    /// The number of T1 transitions attempted
    std::size_t _num_T1s_attempted;

    /// The total number of T1 transitions attempted
    std::size_t _num_T1s_attempted_total;

    /// The frequency of attempted T1 transitions per edge since initialisation
    double _T1_attempt_frequency_acc;

    /// The number of T2 transitions
    std::size_t _num_T2s;

    /// The total number of T2 transitions
    std::size_t _num_T2s_total;

    /// The frequency of T2 transitions per cell since initialisation
    double _T2_frequency_acc;


    std::size_t _estimate_minimizations;

    /// A manager for monitor activity
    std::function<void()> _monitor_mngr = [this](){
        this->_monitor.get_monitor_manager()->check_timer();
        this->__monitor();
        this->_monitor.get_monitor_manager()->emit_if_enabled();
    };

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
                // transition adaptors
                DataIO::transition_adaptor,
                // statistics
                DataIO::interface_length_adaptor,
                // position adaptors
                DataIO::vertices_adaptor<SpaceVec>,
                DataIO::cells_adaptor<SpaceVec>,
                DataIO::edges_adaptor<SpaceVec>,
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
        _T1_frequency_acc(0.),
        _num_T1s_attempted(0),
        _num_T1s_attempted_total(0),
        _T1_attempt_frequency_acc(0.),
        _num_T2s(0),
        _num_T2s_total(0),
        _T2_frequency_acc(0.),
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
                this->_log->trace("  Adding operation '{}' ...", name);

                if (false) { }
                else if (name == "differentiate_Collier") {
                    this->setup_collier(
                            get_as<Config>("Collier", op_cfg, {}));
                    _operations.push_back(
                        build_differentiate_Collier(name, op_cfg,
                            _minimization_params, _collier,
                            _collier_prolog));
                }
                else if (name == "differentiate_domain") {
                    _operations.push_back(
                        build_differentiate_domain(name, op_cfg,
                                                   _minimization_params));
                }
                else if (name == "differentiate_cluster") {
                    _operations.push_back(
                        build_differentiate_cluster(name, op_cfg,
                                                         _minimization_params));
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
                else if (name == "minimize_cell_contacts") {
                    _operations.push_back(
                        build_minimize_cell_contacts(name, op_cfg, 
                                                _minimization_params));
                }
                else if (name == "proliferate") {
                    _operations.push_back(
                        build_proliferate(name, op_cfg, _minimization_params,
                                _log, _monitor_mngr));
                }
                else if (name == "proliferate_generations") {
                    _operations.push_back(
                        build_proliferate_generations(name, op_cfg,
                                _minimization_params,
                                _log, _monitor_mngr));
                }
                else if (name == "register_work_function_term") {
                    _operations.push_back(
                        build_register_work_function_term(
                            name, op_cfg, _minimization_params
                        )
                    );
                }
                else if (name == "simple_shear") {
                    _operations.push_back(
                        build_simple_shear(name, op_cfg, _minimization_params)
                    );
                }
                else if (name == "unregister_work_function_term") {
                    _operations.push_back(
                        build_unregister_work_function_term(
                            name, op_cfg, _minimization_params
                        )
                    );
                }
                else if (name == "update_work_function_term") {
                    _operations.push_back(
                        build_update_work_function_term(
                            name, op_cfg, _minimization_params
                        )
                    );
                }
                else if (name == "void") {
                    _operations.push_back(
                        build_void(name, op_cfg, _minimization_params)
                    );
                }
                else {
                    throw std::invalid_argument(fmt::format(
                        "No operation '{}' available to construct! "
                        "Choose from:\n"
                        "{}",
                        name,
                        "- differentiate_Collier\n"
                        "- differentiate_domain\n"
                        "- differentiate_hair_cluster\n"
                        "- differentiate_NotchDelta\n"
                        "- differentiate_random\n"
                        "- minimize_cell_contacts\n"
                        "- proliferate\n"
                        "- proliferate_generations\n"
                        "- register_work_function_term\n"
                        "- simple_shear\n"
                        "- unregister_work_function_term\n"
                        "- update_work_function_term\n"
                        "- void\n"
                    ));
                }

                const auto& params = std::get<1>(_operations.back());
                auto estimate = params.get_num_minimizations(this->_time_max);
                _estimate_minimizations += estimate;
                this->_log->debug("Added '{}' operation with an estimate of "
                                  "{} iterations.", name, estimate);
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
    
    /// Setup a notch delta model
    void setup_pcp (const Config& cfg = {})
    {
        if (_pcp) {
            return;
        }

        this->_log->debug("Setting up PlanarCellPolarity model from {}",
            cfg.size() ? 
                "custom configuration."
                : 
                fmt::format("configuration within {} model.", this->_name));

        using PCP = PlanarCellPolarity::PlanarCellPolarity<AgentManager>;

        _pcp = std::shared_ptr<PCP>(
            new PCP("PlanarCellPolarity", *this, this->get_am(), cfg,
                std::make_tuple(
                    PlanarCellPolarity::DataIO::time_energy_adaptor,
                    PlanarCellPolarity::DataIO::energy_adaptor
                )
            )
        );

        _pcp_prolog = std::make_shared<bool>(false);        
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

        for (std::size_t it = 0; it < iterates; it++) {
            if (not params.disable) {
                operation(_vertex_model);
            }
            // else: skip 

            if (params.minimization_mode == MinimizationMode::Every) {
                this->_log->debug("   Minimizing energy ...");
                _vertex_model.minimize_energy(params.minimization_params,
                                              _monitor_mngr);
                
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
                                          _monitor_mngr);
                
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
        _num_T1s_total = _vertex_model.get_num_T1s_total();
        _T1_frequency_acc = _vertex_model.\
                get_T1_frequency_accumulated();

        _num_T1s_attempted = _vertex_model.get_num_T1s_attempted_total() - 
                             _num_T1s_attempted_total;
        _num_T1s_attempted_total = _vertex_model.get_num_T1s_attempted_total();
        _T1_attempt_frequency_acc = _vertex_model.\
                get_T1_attempt_frequency_accumulated();
        
        _num_T2s = _vertex_model.get_num_T2s_total() - _num_T2s_total;
        _num_T2s_total = _vertex_model.get_num_T2s_total();
        _T2_frequency_acc = _vertex_model.\
                get_T2_frequency_accumulated();
    }

    /// Monitor model information
    void monitor () {
        // overwrite the progress with progress estimate
        if (_vertex_model.get_num_minimizations() < _estimate_minimizations) {
            this->_monitor.get_monitor_manager()->set_time_entries(
                _vertex_model.get_num_minimizations(),
                std::max(_estimate_minimizations,
                        _vertex_model.get_num_minimizations() + 1));
            this->_monitor.set_entry("time", this->get_time());
            this->_monitor.set_entry("progress",  float(this->get_time())
                                                / float(this->get_time_max()));
        }

        this->_monitor.set_entry("num_cells",
                                 _vertex_model.get_am().cells().size());
        this->_monitor.set_entry("num_T1_transitions",
                                 _vertex_model.get_num_T1s_total());
        this->_monitor.set_entry("num_T1_transitions_attempted",
                                 _vertex_model.get_num_T1s_attempted_total());
        this->_monitor.set_entry("num_T2_transitions",
                                 _vertex_model.get_num_T2s_total());


        _vertex_model.monitor();
        if (_pcp_prolog and *_pcp_prolog) {
            _pcp->monitor();
        }
        if (_collier_prolog and *_collier_prolog) {
            _collier->monitor();
        }
        if (_notch_delta_prolog and *_notch_delta_prolog) {
            _notch_delta->monitor();
        }
    }

    /// The prolog
    /** Performs the following tasks:
     *      1. prolog of vertex model
     *      1. Proliferation of hexagonal lattice
     *      1. prolog operations
     *      1. default prolog tasks
     */
    void prolog () {
        _vertex_model.prolog();

        const auto init_cfg = get_as<Config>("initialisation", this->_cfg);
        auto jiggle_cfg = get_as<Config>("0_add_noise", init_cfg);
        auto prolif_cfg = get_as<Config>("1_by_proliferation", init_cfg);

        std::unique_ptr<OperationBundle> jiggle_op(nullptr);
        std::unique_ptr<OperationBundle> prolif_op(nullptr);

        if (get_as<bool>("enabled", jiggle_cfg)) {
            jiggle_cfg["times"] = std::vector<std::size_t>({});
            jiggle_cfg["iterations_prolog"] = get_as<std::size_t>(
                "iterations_prolog", jiggle_cfg, 1
            );

            jiggle_op = std::make_unique<OperationBundle>(
                build_void("Initial_minimisation",
                           jiggle_cfg, _minimization_params));
            _estimate_minimizations += std::get<1>(*jiggle_op).get_num_minimizations(0);
        }

        if (get_as<bool>("enabled", prolif_cfg)) {
            prolif_cfg["times"] = std::vector<std::size_t>({});
            prolif_cfg["iterations_prolog"] = get_as<std::size_t>(
                "iterations_prolog", prolif_cfg, 1
            );
        
            prolif_op = std::make_unique<OperationBundle>(
                build_proliferate_generations(
                    "initial_proliferation", prolif_cfg,
                    _minimization_params, _log, _monitor_mngr));
            
            _estimate_minimizations += std::get<1>(*prolif_op).get_num_minimizations(0);
        }

        if (jiggle_op != nullptr) {
            this->_log->info("Running initial jiggle ...");
            apply_operation(*jiggle_op, true, false);
        }

        if (prolif_op != nullptr) {
            this->_log->info("Running initial proliferation ...");
            apply_operation(*prolif_op, true, false);
        }


        this->_log->info("Running prolog operations ...");

        for (auto& operation : _operations) {
            apply_operation(operation, true, false);
        }

        _num_T1s = _num_T1s_total;
        _num_T1s_total = _vertex_model.get_num_T1s_total();
        _T1_frequency_acc = _vertex_model.\
                get_T1_frequency_accumulated();

        _num_T1s_attempted = _num_T1s_attempted_total;
        _num_T1s_attempted_total = _vertex_model.get_num_T1s_attempted_total();
        _T1_attempt_frequency_acc = _vertex_model.\
                get_T1_attempt_frequency_accumulated();

        _num_T2s = _num_T2s_total;
        _num_T2s_total = _vertex_model.get_num_T2s_total();
        _T2_frequency_acc = _vertex_model.\
                get_T2_frequency_accumulated();
        
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
        _num_T1s_total = _vertex_model.get_num_T1s_total();
        _T1_frequency_acc = _vertex_model.\
                get_T1_frequency_accumulated();

        _num_T1s_attempted = _vertex_model.get_num_T1s_attempted_total() - 
                             _num_T1s_attempted_total;
        _num_T1s_attempted_total = _vertex_model.get_num_T1s_attempted_total();
        _T1_attempt_frequency_acc = _vertex_model.\
                get_T1_attempt_frequency_accumulated();

        _num_T2s = _vertex_model.get_num_T2s_total() - _num_T2s_total;
        _num_T2s_total = _vertex_model.get_num_T2s_total();
        _T2_frequency_acc = _vertex_model.\
                get_T2_frequency_accumulated();

        _vertex_model.epilog();

        if (_collier) {
            _collier->epilog();
        }

        if (_notch_delta) {
            _notch_delta->epilog();
        }

        if (_pcp) {
            _pcp->epilog();
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

    /// Get energy from a term implementation
    template <typename WFTerm>
    double get_energy (const std::shared_ptr<WFTerm>& functor) const {
        return _vertex_model.get_energy(functor);
    }
    /// Getter for energy of a container of edges and cells, resp.
    double get_energy (
            const AgentContainer<Vertex>& vs,
            const AgentContainer<Edge>& es,
            const AgentContainer<Cell>& cs
    ) const
    {
        return _vertex_model.get_energy(vs, es, cs);
    }

    /// Getter for energy all edges and cells
    double get_energy () const {
        return _vertex_model.get_energy();
    }

    const auto& get_work_function_terms () const {
        return _vertex_model.get_work_function_terms();
    }

    const auto& get_work_function_term_register () const {
        return _vertex_model.get_work_function_term_register();
    }


    /// Getter for the relative energy change from previous to last step
    double get_energy_change () const {
        return _vertex_model.get_energy_change();
    }

    /// Whether the relative change in energy fulfills the equilibrium condition
    bool equilibrium_condition() const {
        return _vertex_model.equilibrium_condition();
    }


    // .. Counter for transitions, etc.. ......................................

    /// Counter for the T1 neighborhood exchange transitions in last iteration
    std::size_t get_num_T1s() const {
        return _num_T1s;
    }

    /// Number of T1s per edge in last iteration
    double get_T1_frequency() const {
        return _num_T1s / static_cast<double>(this->get_am().edges().size());
    }

    /// Total counter for the T1 neighborhood exchange transitions
    std::size_t get_num_T1s_total() const {
        return _num_T1s_total;
    }

    /// Number of T1 per edge
    double get_T1_frequency_accumulated() const {
        return _T1_frequency_acc;
    }


    /// Counter for the attempted T1 neighborhood exchange transitions
    /// in last iteration
    std::size_t get_num_T1s_attempted() const {
        return _num_T1s_attempted;
    }

    /// Number of T1s attempted per edge in last iteration
    double get_T1_attempt_frequency() const {
        return (  _num_T1s_attempted
                / static_cast<double>(this->get_am().edges().size()));
    }

    /// Total counter for the attempted T1 neighborhood exchange transitions
    std::size_t get_num_T1s_attempted_total() const {
        return _num_T1s_attempted_total;
    }

    double get_T1_attempt_frequency_accumulated() const {
        return _T1_attempt_frequency_acc;
    }


    /// Counter for the T2 cell extrusion transitions in last iteration
    std::size_t get_num_T2s() const {
        return _num_T2s;
    }

    /// Number of T2s per cell in last iteration
    double get_T2_frequency() const {
        return _num_T2s / static_cast<double>(this->get_am().cells().size());
    }

    /// Total counter for the T2 cell extrusion transitions
    std::size_t get_num_T2s_total() const {
        return _num_T2s_total;
    }

    double get_T2_frequency_accumulated() const {
        return _T2_frequency_acc;
    }


    /// Getter for vertices
    const AgentManager& get_am () const {
        return _vertex_model.get_am();
    }

    /// PCPVertex::get_cluster_ids
    auto get_cluster_ids (const std::size_t& type) {
        return _vertex_model.get_cluster_ids(type);
    }
    
    /// Add an operation
    void register_operation (OperationBundle op_bundle) {
        _operations.push_back(op_bundle);
        auto [_, params] = op_bundle;
        this->_log->trace("Registered operation '{}'", params.name);

        _estimate_minimizations += params.get_num_minimizations(
                                                this->_time_max);
    }

    /// Return polarity proteins of this edge
    /** See PlanarCellPolarity::get_polarity_proteins */
    std::pair<ProteinVec, ProteinVec> get_polarity_proteins
    (const std::shared_ptr<Edge>& edge, bool flip = false) const
    {
        if (not _pcp) {
            ProteinVec null(1, arma::fill::zeros);
            return std::make_pair(null, null);
        }
        return _pcp->get_polarity_proteins(edge, flip);
    }

    std::size_t proteins_per_edge() const {
        if (not _pcp) {
            return 1;
        }

        return _pcp->proteins_per_edge();
    }

    /// Forward to PlanarCellPolarity::pcp_polarity
    SpaceVec pcp_polarity(const std::shared_ptr<Cell>& cell) const {
        if (not _pcp) {
            return SpaceVec({0., 0.});
        }

        return _pcp->pcp_polarity(cell);
    }

    double get_curvature_boundary() const {
        return _vertex_model.get_curvature_boundary();
    }

    auto fix_number_work_function_terms() {
        return _vertex_model.fix_number_work_function_terms();
    }
};

} // namespace PCPTopology
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPTOPOLOGY_HH
