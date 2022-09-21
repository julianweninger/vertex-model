#ifndef UTOPIA_MODELS_PCPVERTEX_HH
#define UTOPIA_MODELS_PCPVERTEX_HH

// standard library includes
#include <random>
#include <math.h>
#include <queue>

// third-party library includes
#include <boost/circular_buffer.hpp>

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/apply.hh>
#include <utopia/core/types.hh>
#include <utopia/core/agent_manager.hh>
#include <utopia/data_io/data_manager/defaults.hh>

#include "space.hh"

#include "entities.hh"
#include "entities_manager.hh"
#include "initialisation.hh"
#include "transitions.hh"
#include "utils.hh"

#include "work_function.hh"
#include "minimization.hh"
#include "boundary.hh"

namespace Utopia {
namespace Models {
namespace PCPVertex {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<DefaultRNG, WriteMode::managed,
                                      Space::CustomSpace<2>>;


// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// The PCPVertex Model
/** This model implements the relaxation of an energy function of a planar
 *  cell polarity system towards a local minimum.
 * 
 *  The energy function currently includes the following terms
 *      - area elasticity
 *      - line tension, alias surface tension
 * 
 *  The energy relaxation is performed in one of the following ways
 *      - along steepest descent with fixed step size
 * 
 *  The model accounts for the following topological changes
 *      - T1 transition: cell intercalation changes neighbourhood of cells,
 *          i.e. an edge shrinks to neglecting length and is replaced with 
 *          an orthogonal edge connecting the two next neighbour cells.
 *          Thus the cells adjacent to the removed edge are no longer neighbours
 *      - T2 transition: cell extrusion when cell area shrinks below threshold
 *          value
 */
class PCPVertex:
    public Model<PCPVertex, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<PCPVertex, ModelTypes>;

    /// Data type for the model time
    using Time = typename ModelTypes::Time;

    /// Type of the space the model is living in
    using Space = typename ModelTypes::Space;

    /// The type of a coordinate / vector in space
    using SpaceVec = typename Space::SpaceVec;

    /// The type of a config
    using Config = Utopia::DataIO::Config;

    /// The manager of the entities (agents) of this model
    using AgentManager = EntitiesManager<Model<PCPVertex, ModelTypes>>;

    /// The type of a Vertex
    using Vertex = typename AgentManager::Vertex;

    /// The type of an Edge
    using Edge = typename AgentManager::Edge;

    /// The type of a Cell
    using Cell = typename AgentManager::Cell;

    /// The types of a cell
    using CellState = typename Cell::State;

    using OrderedEdgeContainer = typename AgentManager::OrderedEdgeContainer;

    /// The type of a rule function acting on vertices of the agent manager
    using RuleFuncVertex = typename AgentManager::RuleFuncVertex;

    /// The type of a rule function acting on edges of the agent manager
    using RuleFuncEdge = typename AgentManager::RuleFuncEdge;

    /// The type of a rule function acting on cells of the agent manager
    using RuleFuncCell = typename AgentManager::RuleFuncCell;


    using WFVertexTerm = WorkFunction::WorkFunctionVertexTerm<PCPVertex>;

    using WFEdgeTerm = WorkFunction::WorkFunctionEdgeTerm<PCPVertex>;

    using WFCellTerm = WorkFunction::WorkFunctionCellTerm<PCPVertex>;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space

    // -- Members -------------------------------------------------------------
    /// The manager of the model's entities
    AgentManager _am;


    // -- Minimization parameters ---------------------------------------------

    // PARAMETERS
    MinimizationParams _default_minimization_params;

    /// timestep scaling
    double _dt;

    /// The tolerance during minimization
    double _minimization_tolerance;



    // -- Mechanical parameters -----------------------------------------------
    std::unordered_map<std::string,
                       std::shared_ptr<WFVertexTerm>> _vertex_terms;
    std::unordered_map<std::string, std::shared_ptr<WFEdgeTerm>> _tensions;
    std::unordered_map<std::string, std::shared_ptr<WFCellTerm>> _pressures;

    std::set<std::string> _work_function_terms_disabled;


    enum BoundaryType {
        Periodic,
        PeriodicSkewed,
        PeriodicCurved,
        Free,
        Fixed,
        FixedCurved,
        Potential
    } _boundary_type;

    bool boundary_is_fixed () const {
        if (_boundary_type == BoundaryType::Fixed) {
            return true;
        }
        if (_boundary_type == BoundaryType::FixedCurved) {
            return true;
        }
        return false;
    }


    // -- transition parameters -----------------------------------------------

    /// Whether topological transitions are enabled
    bool _enable_transitions;

    /// Whether T1 transitions are enabled
    bool _enable_T1_transitions;

    /// Edges shorter than this value are replaced in a T1 transition
    const double _T1_threshold;

    /// The factor by which new edges are longer than threshold
    /** \note Value from configuration is factor k_sep, but value stored is the
     *        length of separation: \f$ d_{sep} = k_{sep} * d_{min} \f$
     */
    const double _T1_separation;

    /// The probability, that a T1 transition occurs
    const double _T1_probability;

    /// The characteristic height of the energy barrier at T1 transitions
    /** The probability to perform a T1 transition is 
     *  \f$ p = \exp(-\Delta E / _T1_barrier) \f$, which is 1 for 
     *  \f$\Delta E < 0\f$.
     */
    const double _T1_barrier;

    /// A timeout after attempted unsuccessful T1 transition
    /** Default is 0
     */
    const std::size_t _T1_timeout;
    
    /// Whether T2 transitions are enabled
    bool _enable_T2_transitions;

    /// Cells with area smaller than this value are removed in T2 transition
    double _T2_threshold;

    
    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::uniform_real_distribution<double> _prob_distr;

    /// A [0,1]-range normal distribution
    std::normal_distribution<double> _normal_distr;

    // .. Temporary objects ...................................................
protected:
    boost::circular_buffer<double> _energy_buffer;

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

    /// A counter for the completed minimizations
    std::size_t _num_minimizations;

    /// The status of the model in minimize_energy for datamanager
    enum Status {
        /// In minimization
        Minimization,

        /// Has been externally perturbed
        Perturbed,

        /// Has been jiggled
        Jiggled
    } _status;

public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the PCPVertex model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel, typename... WriterArgs>
    PCPVertex (const std::string name,
        ParentModel &parent_model,
        const DataIO::Config& custom_cfg = {},
        std::tuple<WriterArgs...> &&writer_args = {},
        const DataIO::Default::DefaultDecidermap< PCPVertex >&
            w_deciders = Utopia::DataIO::Default::default_deciders< PCPVertex >,
        const DataIO::Default::DefaultTriggermap< PCPVertex >&
            w_triggers = Utopia::DataIO::Default::default_triggers< PCPVertex >)
    :
        // Initialize first via base model
        Base(name, parent_model, custom_cfg,
             writer_args, w_deciders, w_triggers),

        _am(*this),
        
        // Get member paramters from cfg
        _default_minimization_params(get_as<Config>("minimization",
                                                    this->_cfg)),
        _dt(_default_minimization_params.dt),
        _minimization_tolerance(_default_minimization_params.tolerance),
        _vertex_terms({}),
        _tensions({}),
        _pressures({}),
        _work_function_terms_disabled({}),
        
        _enable_transitions(
            get_as<bool>("enable_transitions", this->_cfg, true)),
        _enable_T1_transitions(
            get_as<bool>("enable_T1_transitions", this->_cfg, true)),
        _T1_threshold(get_as<double>("T1_threshold", this->_cfg)),
        _T1_separation(
            _T1_threshold * get_as<double>("T1_separation_factor",this->_cfg)),
        _T1_probability(get_as<double>("T1_probability", this->_cfg)),
        _T1_barrier(get_as<double>("T1_barrier", this->_cfg)),
        _T1_timeout(get_as<std::size_t>("T1_timeout", this->_cfg, 0)),
        _enable_T2_transitions(
            get_as<bool>("enable_T2_transitions", this->_cfg, true)),
        _T2_threshold(get_as<double>("T2_threshold", this->_cfg)),
        _prob_distr(0.,1.),
        _normal_distr(0.,1.),
        _energy_buffer(get_as<std::size_t>("Energy_buffer_size", this->_cfg)),
        _num_T1s(0),
        _num_T1s_total(0),
        _T1_frequency_acc(0.),
        _num_T1s_attempted(0),
        _num_T1s_attempted_total(0),
        _T1_attempt_frequency_acc(0.),
        _num_T2s(0),
        _num_T2s_total(0),
        _T2_frequency_acc(0.),
        _num_minimizations(0)
    {
        this->setup_work_function(get_as<Config>(
            "work_function_terms",
            this->_cfg
        ));

        double initial_jiggle(get_as<double>("initial_jiggle", this->_cfg, 0.));
        if (initial_jiggle > 1.e-12) {
            jiggle_vertices(initial_jiggle);
        }
        
        this->_log->info("Model initialized.");
    }


private:
    // .. Setup functions .....................................................
    void setup_work_function (const Config& cfg) {
        using namespace WorkFunction;

        this->_log->info("Setting up work-function from {} configuration entr{}"
                         " ...", cfg.size(), cfg.size() != 1 ? "ies" : "y");

        if (not cfg.size()) {
            throw std::runtime_error("No term registered to work-function! "
                "Note that also initially disabled terms must be registered!");
        }

        // Otherwise, require a sequence
        if (not cfg.IsSequence()) {
            throw std::invalid_argument("The config for initializing the "
                "work-function must be a sequence!");
        }

        // Iterate over the sequence of mappings
        for (const auto& terms : cfg) {
            // ops.IsMap() == true
            // The top `ops` keys are now the names of the desired environment
            // functions. Iterate over those ...
            for (const auto& term_pair : terms) {
                const std::string term = term_pair.first.as<std::string>();
                const Config& params = term_pair.second;

                if (not get_as<bool>("enabled", params, true)) {
                    this->_log->debug("Pre-registering work-function term "
                        "'{}', but disabling it.");

                    _work_function_terms_disabled.insert(term);                    
                }

                this->_log->debug("Registering work-function term '{}' ...",
                                  term);
                _work_function_terms_disabled.insert(term);                    

                if (term == "area_elasticity") {
                    register_pressure(
                        term,
                        std::make_shared<AreaElasticity<PCPVertex>>(
                            params, *this
                        )
                    );
                }
                else if (term == "boundary_stripe_potential") {
                    register_force(
                        term,
                        std::make_shared<BoundaryStripePotential<PCPVertex>>(
                            params, *this
                        )
                    );
                }
                else if (term == "cell_contractility") {
                    register_pressure(
                        term,
                        std::make_shared<CellContractility<PCPVertex>>(
                            params, *this
                        )
                    );
                }
                else if (term == "edge_contractility") {                  
                    register_tension(
                        term,
                        std::make_shared<EdgeContractility<PCPVertex>>(
                            params, *this
                        )
                    );
                }
                else if (term == "linetension") {
                    register_tension(
                        term,
                        std::make_shared<Linetension<PCPVertex>>(params, *this)
                    );
                }
                else if (term == "linetension_fluctuations") {
                    register_tension(
                        term,
                        std::make_shared<LinetensionFluctuations<PCPVertex>>(
                            params, *this
                        )
                    );
                }
                else if (term == "linetension_heterotypic") {
                    register_tension(
                        term,
                        std::make_shared<LinetensionHeterotypic<PCPVertex>>(
                            params, *this
                        )
                    );
                }
                else if (term == "shape_elasticity") {
                    register_pressure(
                        term,
                        std::make_shared<ShapeElasticity<PCPVertex>>(
                            params, *this
                        )
                    );
                }
                else {
                    throw std::runtime_error(fmt::format(
                        "No term `{}` known in PCPVertex namespace. "
                        "Use the `register_work_function_term` interface, or "
                        "choose one of the following available terms:\n"
                        " - area_elasticity\n"
                        " - boundary_stripe_potential\n"
                        " - cell_contractility\n"
                        " - edge_contractility\n"
                        " - linetension\n"
                        " - linetension_fluctuations\n"
                        " - linetension_heterotypic\n"
                        " - shape_elasticity\n"
                        "", term
                    ));
                }

            }
        }
    }

    // .. Force setter functions ..............................................

    /// Set the gradient of energy within the vertices
    /** \details This function applies the forces arising from the different
     *           energy terms.
     *  \note    When adding energy terms, remember to add their gradient here!
     */
    void compute_forces () {

        compute_tensions_and_pressures();

        for (const auto& vertex : _am.vertices()) {
            vertex->state.reset_force();
        }

        for (const auto& [name, functor] : _vertex_terms) {
            for (const auto& vertex : _am.vertices()) {
                vertex->state.add_force(functor->compute_force(vertex));
            }
        }

        // apply forces from tensions and pressures
        for (const auto& edge : _am.edges()) {
            SpaceVec displ = _am.displacement(edge);
            SpaceVec director = displ / arma::norm(displ);

            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;

            const double& T = edge->state.tension;
            a->state.add_force(+ T * director);
            b->state.add_force(- T * director);


            // left and right cell when facing along edge
            const auto& [cl, cr] = _am.adjoints_of<true>(edge);

            // get pressure difference between left and right cell
            double delta_P = 0.;
            if (cl != nullptr) {
                delta_P += cl->state.pressure;
            }
            if (cr != nullptr) {
                delta_P -= cr->state.pressure;
            }

            // force is normal to edge
            double Fa = 0.5 * delta_P;
            a->state.add_force(-1. * Fa * SpaceVec({displ[1], -displ[0]}));
            b->state.add_force(-1. * Fa * SpaceVec({displ[1], -displ[0]}));
        }

        // fix the boundary
        if (not _space->periodic and boundary_is_fixed()) {
            apply_rule<Update::sync>(
                [this](const auto& vertex) {
                    auto state = vertex->state;
                    if (this->_am.is_boundary(vertex)) {
                        state.fix_in_space = true;
                    }
                    else {
                        state.fix_in_space = false;
                    }
                    return state;
                },
                _am.vertices()
            );
        }

        // decider on whether to fix vertex
        std::function<bool(const std::shared_ptr<Vertex>& vertex)> fix_vertex;
        fix_vertex = [](const auto& vertex) {
            if (vertex->state.fix_in_space) {
                return true;
            }
            return false;
        };
        if (not _space->periodic and boundary_is_fixed()) {
            fix_vertex = [this, fix_vertex](const auto& vertex) {
                return (fix_vertex(vertex) or this->_am.is_boundary(vertex));
            };
        }
        
        // set forces on fixed vertices (e.g. boundary) to zero
        apply_rule<Update::sync>(
            [fix_vertex](const auto& vertex) {
                auto state = vertex->state;
                if (fix_vertex(vertex)) { 
                    state.reset_force();
                }
                return state;
            },
            _am.vertices()
        );
    }
    
    void compute_tensions_and_pressures () {
        for (const auto& edge : _am.edges()) {
            double tension = 0.;
            for (const auto& [name, functor] : _tensions) {
                tension += functor->compute_tension(edge);
            }
            edge->state.tension = tension;
        }

        for (const auto& cell : _am.cells()) {
            double tension = 0.;
            for (const auto& [name, functor] : _pressures) {
                tension += functor->compute_tension(cell);
            }
            for (const auto& [edge, f] : cell->custom_links().edges) {
                edge->state.tension += tension;
            }

            double pressure = 0.;
            for (const auto& [name, functor] : _pressures) {
                pressure += functor->compute_pressure(cell);
            }
            cell->state.pressure = pressure;
        }
    }

    /** The update of position
     * 
     *  Move vertex proportional to the gradient of energy (force)
     * 
     *  @param v    The pointer to the vertex to update
     */
    const RuleFuncVertex update_position = [this](const auto& vertex) {
        _am.move_by(vertex, vertex->state.get_force() * this->_dt);
        return vertex->state;
    };
    

    // -- The algorithm    ----------------------------------------------------
    // see algorithm.hh
    double steepest_gradient_step ();
    double perform_update_step();
    bool perform_transitions(bool enabled);


    // -- Helper functions ----------------------------------------------------

public:
    // -- Public Interface ----------------------------------------------------
    void jiggle_vertices(double intensity);
    void differentiate_hair_cells_random(double fraction);
    template <class NotchDelta>
    void differentiate_hair_cells_NotchDelta(
        std::shared_ptr<NotchDelta> notch_delta, int steps);
        
    /// Perform a cell division on specific cell
    /** Divides a specific cell into two identical cells with properties derived
     *  from the common parent cell. 
     *  The division is performed at a given angle through the parent cell's
     *  center. This defines the axis of division that will form a new edge 
     *  between the two new cells.
     * 
     *  \param cell         the cell that is to be divided
     *  \param division_angle   angle (in rad) at which the cell is divided
     */
    auto divide_cell(std::shared_ptr<Cell> cell, double division_angle) {
        return _am.divide_cell(cell, division_angle);
    }

    void increase_domain_size(double area, bool deform_plastic);
    double stretch_domain(SpaceVec stretch, bool compensate,
        bool fix_hc_area, bool fix_sc_area, bool deform_plastic);
    SpaceVec skew_domain(SpaceVec add_skew, bool absolute, bool deform_plastic);
    
    void curve_boundary(double add_curvature, bool deform_plastic);
    double get_curvature_boundary() const {
        if (_space->periodic) {
            return _space->get_curvature();
        }
        return 0.;
    }


    // .. Simulation Control ..................................................
    void init_minimization ();

    /// Iterate a single step
    /** \details Rules applied
     *      -# perform_transitions()
     *      -# perform_update_step()
     *      -# tracking of variables
     */
    void perform_step () {
        perform_transitions(_enable_transitions);

        for (const auto& [name, F] : _vertex_terms) {
            F->update(this->_dt);
        }
        for (const auto& [name, T] : _tensions) {
            T->update(this->_dt);
        }
        for (const auto& [name, P] : _pressures) {
            P->update(this->_dt);
        }

        if (_space->get_curvature() > 1.e-8) {
            throw std::runtime_error(fmt::format("Cannot perform step with "
                "curved periodic boundary conditions. Curvature {} > 0",
                _space->get_curvature()));
        }

        double E = perform_update_step();

        this->_log->trace("Energy changed by {}", E - _energy_buffer.back());
        _energy_buffer.push_back(E);
    }

    void prolog () {
        this->init_minimization();

        return this->__prolog();
    }
    
    /// Monitor model information
    void monitor () {
        this->_monitor.set_entry("time", this->_time);

        double mean_energy = 0.;
        if (_energy_buffer.size() >= 1) {
            mean_energy = std::accumulate(
                _energy_buffer.begin(),
                _energy_buffer.end(),
                0.
            );
            mean_energy /= _energy_buffer.size();

            this->_monitor.set_entry("energy", mean_energy);
        }
        if (_energy_buffer.size() >= 2) {
            this->_monitor.set_entry("energy_gradient",  get_energy_change());
        }
    }

    /// Minimize the energy
    /** Iterate this model until the change of energy is smaller than a given 
     *  tolerance.
     * 
     *  \param  params  The collection of parameter required.
     * 
     *  \details 1. Jiggle vertices using jiggle_vertices()
     *           2. Minimize energy to given tolerance by iterating this model.
     *              The tolerance may be reduced 
     *              See perform_step() for more details.
     *           3. Repeat 1. and 2. `num_repeat` times.
     * 
     *  \return num steps performed
     */
    std::size_t minimize_energy(const MinimizationParams& params,
                                std::function<void()> monitor_mngr = [](){ 
                                    return; })
    {
        _status = Status::Perturbed;
        this->increment_time();
        this->_datamanager(*this);            
        this->_log->debug("Incremented time (initial condition after external "
                          "perturbation): {:7d}",
                          this->_time);

        const auto time_0 = this->get_time();
        this->_log->debug("Minimizing energy from step {} with {} repeats ...",
                          time_0, params.num_repeat);
        _status = Status::Minimization;


        _dt = params.dt;

        if (std::get<1>(params.linetension_fluctuations) > 0.) {
            auto [registered, active] = is_registered_term(
                "linetension_fluctuations");
            if (not registered) {
                std::runtime_error("Cannot introduce `linetension_fluctuations`"
                    " because no such term is registered! Register it at setup "
                    "of vertex-model!");
            }

            Config cfg;
            cfg["timescale"] = std::get<0>(params.linetension_fluctuations);
            cfg["amplitude"] = std::get<1>(params.linetension_fluctuations);
            if (not active) {
                using T = WorkFunction::LinetensionFluctuations<PCPVertex>;

                register_tension(
                    "linetension_fluctuations",
                    std::make_shared<T>(cfg, *this)
                );
            }
            else {
                _tensions["linetension_fluctuations"]->update_parameters(cfg);
            }
        }
        else if (std::get<1>(is_registered_term("linetension_fluctuations"))) {
            Config cfg;
            cfg["amplitude"] = 0.;
            _tensions["linetension_fluctuations"]->update_parameters(cfg);
        }

        for (std::size_t i = 0; i < params.num_repeat; i++)
        {
            // jiggle vertices if required
            if (params.jiggle_intensity > 1.e-12) {
                _status = Status::Jiggled;
                this->jiggle_vertices(params.jiggle_intensity);
                this->increment_time();
                this->_datamanager(*this);            
                this->_log->debug("Incremented time after jiggling: {:7d}",
                                  this->_time);
                
                // reset status
                _status = Status::Minimization;
            }

            // initialize minimization
            const auto time_start = this->get_time();
            bool minimum_reached = false;

            this->init_minimization();
            _dt = params.dt;

            double tolerance;
            if (i+1 == params.num_repeat) {
                tolerance = params.tolerance;
            }
            else {
                tolerance = params.jiggle_tolerance;
            }
            _minimization_tolerance = tolerance;

            const bool tmp_enable_transitions = _enable_transitions;
            // NOTE save status and restore at the end

            // iterate a fixed number of steps
            if (params.num_steps > 0) {
                this->_log->debug("Iterating vertex model for {} steps",
                                  params.num_steps);
                for (std::size_t step = 0; step < params.num_steps; step++) {
                    // disable topological transitions in first iteration
                    if (step == 0) { _enable_transitions = false; }
                    else { _enable_transitions = tmp_enable_transitions; }

                    this->iterate();
                    monitor_mngr();

                    if (stop_now.load()) {
                        this->_log->warn("Was told to stop. Not iterating "
                            "further ...");
                        throw GotSignal(received_signum.load());
                    }
                }

                // end here after fixed number of steps
                minimum_reached = true;
            }
            else {
                this->_log->debug("Minimizing energy from step {} "
                                  "in {:d} / {:d} repeat ...",
                                  time_start, i+1, params.num_repeat);
            }

            // iterate until minimum reached
            while (not minimum_reached) {
                std::size_t step = this->get_time() - time_start;
                // disable topological transitions in first iteration
                if (step == 0) { _enable_transitions = false; }
                else { _enable_transitions = tmp_enable_transitions; }

                this->iterate();
                monitor_mngr();

                if (step >= params.min_steps) {
                    minimum_reached = equilibrium_condition();
                }
                double energy_change = get_energy_change();

                if (not minimum_reached) {
                    this->_log->trace("Energy changed by {} in last {} steps",
                                      energy_change, _energy_buffer.size());
                }
                
                if (not minimum_reached 
                    and this->get_time() - time_start >= params.max_steps)
                {
                    throw std::runtime_error(fmt::format(
                        "Equilibrium not reached within {} steps at a "
                        "tolerance of {}! Relative energy change in last {} "
                        "steps was {}.",
                        params.max_steps, tolerance, _energy_buffer.size(),
                        energy_change));
                }

                if (stop_now.load()) {
                    this->_log->warn("Was told to stop. Not iterating "
                        "further ...");
                    throw GotSignal(received_signum.load());
                }
            }
            this->_log->debug("  Energy minimized in {} steps.",
                              this->get_time() - time_start);

            _num_minimizations++;

            _enable_transitions = tmp_enable_transitions;
            // NOTE restore initial state of transitions allowed
        }

        return this->get_time() - time_0;
    }

    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

private:
    /// Checks that a term with this name is disabled and can thus be registered
    void prepare_register_term (const std::string& name) {
        if (not _work_function_terms_disabled.erase(name)) {
            this->_log->error("A list of pre-registered work-function terms "
                "available for registration:");
            for (const auto& n : _work_function_terms_disabled) {
                this->_log->error("  {}", n);
            }
            throw std::runtime_error(fmt::format(
                "Cannot register work-function term with name `{}`, because it "
                "was not pre-registered! See above for available "
                "pre-registered terms.",
                name
            ));
        }

        if (   _tensions.find(name) != _tensions.end()
            or _pressures.find(name) != _pressures.end()
            or _vertex_terms.find(name) != _vertex_terms.end())
        {
            throw std::runtime_error(fmt::format(
                "Cannot register work-function term with name `{}`, because a "
                "term with that name is already registered!",
                name
            ));
        }
    }

public:
    // .. Public energy terms .................................................
    void register_force(
        std::string name,
        const std::shared_ptr<WFVertexTerm>& force
    )
    {
        prepare_register_term(name);

        _vertex_terms.emplace(name, force);

        this->_log->info("Successfully registered force `{}`.", name);
    }
    
    void register_tension(
        std::string name,
        const std::shared_ptr<WFEdgeTerm>& tension
    )
    {
        prepare_register_term(name);

        _tensions.emplace(name, tension);

        this->_log->info("Successfully registered tension `{}`.", name);
    }

    void register_pressure(
        std::string name,
        const std::shared_ptr<WFCellTerm>& pressure
    )
    {
        prepare_register_term(name);

        _pressures[name] = pressure;

        this->_log->info("Successfully registered pressure `{}`.", name);
    }
    
    /// Check work function term register
    /**
     *  Returns:
     *      - bool: Whether term is registered
     *      - bool: Whether term is actively registered
     */
    std::pair<bool, bool> is_registered_term (const std::string& name) const {
        if (_vertex_terms.find(name) != _vertex_terms.end()) {
            return std::make_pair(true, true);
        }
        if (_tensions.find(name) != _tensions.end()) {
            return std::make_pair(true, true);
        }
        if (_pressures.find(name) != _pressures.end()) {
            return std::make_pair(true, true);
        }

        else if (   _work_function_terms_disabled.find(name)
                 != _work_function_terms_disabled.end())
        {
            return std::make_pair(true, false);
        }

        return std::make_pair(false, false);
    }


    /// Remove a term from the work function register
    /** Adds the term to the disabled terms (can be reused).
     */
    bool unregister_term(const std::string& name) {
        bool found = false;
        if (not found) {
            found = found or _vertex_terms.erase(name);
        }
        if (not found) {
            found = found or _tensions.erase(name);
        }
        if (not found) {
            found = found or _pressures.erase(name);
        }

        if (found) {
            _work_function_terms_disabled.insert(name);

            this->_log->info("Removed `{}` term from work function.", name);
        }
        else {
            this->_log->debug("Unable t remove `{}` term from work function. "
                              "It is not registered. ", name);
        }

        return found;
    }

    /// Get energy for an edge
    double get_energy (const std::shared_ptr<WFVertexTerm>& functor) const {
        double E = 0.;
        for (const auto& vertex : _am.vertices()) {
            E += functor->compute_energy(vertex);
        }
        return E;
    }

    /// Get energy for an edge
    double get_energy (const std::shared_ptr<WFEdgeTerm>& functor) const {
        double E = 0.;
        for (const auto& edge : _am.edges()) {
            E += functor->compute_energy(edge);
        }
        return E;
    }

    /// Get energy for a cell
    double get_energy (const std::shared_ptr<WFCellTerm>& functor) const {
        double E = 0.;
        for (const auto& cell : _am.cells()) {
            E += functor->compute_energy(cell);
        }
        return E;
    }

    /// Getter for energy of a container of edges and cells, resp.
    double get_energy (
            const AgentContainer<Vertex>& vs,
            const AgentContainer<Edge>& es,
            const AgentContainer<Cell>& cs
    ) const
    {
        double E = 0.;
        for (const auto& [name, functor] : _vertex_terms) {
            for (const auto& vertex : vs) {
                E += functor->compute_energy(vertex);
            }
        }
        for (const auto& [name, functor] : _tensions) {
            for (const auto& edge : es) {
                E += functor->compute_energy(edge);
            }
        }
        for (const auto& [name, functor] : _pressures) {
            for (const auto& cell : cs) {
                E += functor->compute_energy(cell);
            }
        }

        return E;
    }

    /// Getter for energy all edges and cells
    double get_energy () const {
        double E = 0.;
        for (const auto& [name, functor] : _vertex_terms) {
            E += get_energy(functor);
        }
        for (const auto& [name, functor] : _tensions) {
            E += get_energy(functor);
        }
        for (const auto& [name, functor] : _pressures) {
            E += get_energy(functor);
        }

        return E;
    }


    /// Getter for the relative energy change from previous to last step
    double get_energy_change () const {
        if (_energy_buffer.size() < 2) {
            return std::numeric_limits<double>::lowest();
        }
        double mean_energy = std::accumulate(
            _energy_buffer.begin(),
            _energy_buffer.end(),
            0.
        );
        mean_energy /= _energy_buffer.size();

        double mid = (_energy_buffer.size() - 1) / 2.;
        double trend = 0.;
        double norm = 0.;
        for (std::size_t i = 0; i < _energy_buffer.size(); i++) {
            double val = _energy_buffer[_energy_buffer.size() - 1 - i];
            trend += (val - mean_energy) * (mid - i);
            norm += std::pow(i - mid, 2.);
        }

        return trend / norm;
    }

    /// Whether the relative change in energy fulfills the equilibrium condition
    bool equilibrium_condition() const {
        if (_status != Status::Minimization) {
            return false;
        }

        double energy_change = get_energy_change();
        return (energy_change >= -_minimization_tolerance);
    }


    // .. Public energy terms for subset of entities ..........................
    auto get_work_function_terms () const {
        return std::make_tuple(
            _vertex_terms,
            _tensions,
            _pressures,
            _work_function_terms_disabled
        );
    }

    std::shared_ptr<WFCellTerm> get_work_function_cell_term
    (const std::string& name) const
    {
        const auto it = _pressures.find(name);

        if (it != _pressures.end()) {
            return std::get<1>(*it);
        }
        if (   _work_function_terms_disabled.find(name)
            != _work_function_terms_disabled.end())
        {
            return nullptr;
        }

        throw std::runtime_error(fmt::format(
            "Cannot find work-function term with name `{}`, because a "
            "term with that name is already registered!",
            name
        ));
    }
        
    // .. Counter for transitions, etc.. ......................................

    /// Counter for the T1 neighborhood exchange transitions in last iteration
    std::size_t get_num_T1s() const {
        return _num_T1s;
    }

    /// Number of T1s per edge in last iteration
    double get_T1_frequency() const {
        return _num_T1s / static_cast<double>(_am.edges().size());
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
        return _num_T1s_attempted / static_cast<double>(_am.edges().size());
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
        return _num_T2s / static_cast<double>(_am.cells().size());
    }

    /// Total counter for the T2 cell extrusion transitions
    std::size_t get_num_T2s_total() const {
        return _num_T2s_total;
    }

    double get_T2_frequency_accumulated() const {
        return _T2_frequency_acc;
    }


    /// Counter for the energy minimizations completed
    std::size_t get_num_minimizations() const {
        return _num_minimizations;
    }

    /// The entities manager - vertices, edges, cells
    const AgentManager& get_am () const {
        return _am;
    }

    /// Label all clusters of same type cells
    std::unordered_map<std::shared_ptr<Cell>, std::size_t> get_cluster_ids (
        const std::size_t type
    ) const
    {
        std::unordered_map<std::shared_ptr<Cell>, std::size_t> cluster_ids;
        cluster_ids.reserve(_am.cells().size());
        
        std::size_t cluster_id = 0;
        for (const auto& cell : _am.cells()) {
            cluster_ids[cell] = 0;
        }

        std::set<std::shared_ptr<Cell>> cluster_members;
        std::queue<std::shared_ptr<Cell>> cluster_members_candidates;

        for (const auto& cell : _am.cells()) {
            // skip assigned cells
            if (cluster_ids.at(cell) != 0)
            {
                continue;
            }

            if (cell->state.type != type) {
                continue;
            }

            cluster_id++;
            cluster_ids[cell] = cluster_id;
            cluster_members_candidates.push(cell);

            while (cluster_members_candidates.size()) {
                const auto& candidate = cluster_members_candidates.front();
                if (candidate->state.type == type) {
                    const auto [it, insert] = cluster_members.insert(candidate);
                    if (insert) {
                        for (const auto& n : _am.neighbors_of(candidate)) {
                            cluster_members_candidates.push(n);
                        }
                    }
                }
                cluster_members_candidates.pop();
            }

            for (const auto& m : cluster_members) {
                cluster_ids[m] = cluster_id;
            }

            cluster_members.clear();
        }

        return cluster_ids;
    }


    // .. Model properties ....................................................

    /// Enable or disable transitions
    void enable_transitions (bool enable_T1_transitions = true,
                             bool enable_T2_transitions = true)
    {
        _enable_T1_transitions = enable_T1_transitions;
        _enable_T2_transitions = enable_T2_transitions;
        _enable_transitions = (enable_T1_transitions or enable_T2_transitions);
    }
}; // class PCPVertex


} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPVERTEX_HH
