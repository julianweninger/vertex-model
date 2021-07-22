#ifndef UTOPIA_MODELS_PCPVERTEX_HH
#define UTOPIA_MODELS_PCPVERTEX_HH

// standard library includes
#include <random>
#include <math.h>

// third-party library includes

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

namespace Utopia {
namespace Models {
namespace PCPVertex {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// The collection of parameters used for energy minimization
struct MinimizationParams {
    /// The tolerance in energy change 
    double tolerance;

    enum UpdateScheme {
        SteepestGradient,
        SteepestGradientAdaptive,
        ConjugateGradient
    } update_scheme;

    /// The default timestep
    double dt;

    /// The maximum number of steps per minimization
    std::size_t max_steps;

    /// Iterate for a fixed number of steps
    /** If num_steps = 0, energy minimized for tolerance
     */
    std::size_t num_steps;

    /// The temperature for random brownian motion
    /** With mean 0 and variance \f$ \sigma^2 = 2 T \f$.
     */
    std::normal_distribution<double> temperature;

    /// Linetension fluctuation parameter
    /** Parameter fluctuations are implemented as Ornstein-Uhlenbeck process
     * 
     *  \f$  \frac{d\Lambda_{mn}}{dt} = - \frac{1}{\tau_\Lambda}
     *      (\Lambda_{mn}(t) - \Lambda_0)
     *      + \Delta \Lambda \sqrt{2 / \tau_\Lambda} \Theta_{mn}(t)
     *  \f$
     * 
     *  with the first value \f$ \tau \f$ and the second value
     *  \f$ \Delta \Lambda \f$.
     */
    std::pair<double, double> linetension_fluctuations;

    /// Area fluctuation parameter
    /** Parameter fluctuations are implemented as Ornstein-Uhlenbeck process
     * 
     *  \f$  \frac{dA_{i}}{dt} = - \frac{1}{\tau_A}
     *      (A_{i}(t) - A_0)
     *      + \Delta A * A_0 \sqrt{2 / \tau_A} \Theta_{i}(t)
     *  \f$
     * 
     *  with the first value \f$ \tau \f$,
     *  the second value the \f$ \Delta A \f$ relative fluctuations, 
     *  and the third value the minimum area (using a lognormal distribution).
     */
    std::tuple<double, double, double> area_fluctuations;

    /// How to evolve the activity of edge contractility
    /** Defines the activation and deactivation rate of contractility on every
     *  edge.
     * 
     * Global steady state density expected as
     * \f$
     *      \rho_c^\star = \frac{a}{a + b}
     * \f$.
     */
    std::pair<double, double> contractility_activity;

    /// The number of jiggling the vertices
    /** \details The first jiggle is applied before the first minimization,
     *           then the energy is minimized up to `jiggle_tolerance`.
     *           This is repeated up to the final minimization, which is 
     *           done up `tolerance`.
     * 
     *  \note Must be > 0
     */
    std::size_t num_repeat;

    /// The reduced tolerance for the preliminary minimizations
    /** \details Only the final minimization is done with `tolerance`.
     *  \note   Optional parameter, default is `tolerance`.
     *  \note   The value of `jiggle_tolerance` cannot be smaller than
     *          `tolerance`.
     */
    double jiggle_tolerance;

    /// The intensity of the jiggling
    /** \details The intensity scales relative to the typical lengthscale of a
     *           cell \f$ l = \sqrt{A_{domain} / num_cells} \f$
     */
    double jiggle_intensity;

    template <typename Config>
    MinimizationParams(const Config& cfg)
    :
        tolerance(get_as<double>("tolerance", cfg)),
        update_scheme(setup_update_scheme(cfg)),
        dt(get_as<double>("dt", cfg)),
        max_steps(get_as<std::size_t>("max_steps", cfg)),
        num_steps(get_as<std::size_t>("num_steps", cfg, 0)),
        temperature(0., sqrt(2 * get_as<double>("temperature", cfg, 0.))),
        linetension_fluctuations(
            std::make_pair(
                get_as<double>("linetension_fluctuation_tau", cfg, 1.),
                get_as<double>("linetension_fluctuation", cfg, 0.)
            )
        ),
        area_fluctuations(
            std::make_tuple(
                get_as<double>("area_fluctuation_tau", cfg, 1.),
                get_as<double>("area_fluctuation", cfg, 0.),
                get_as<double>("area_fluctuation_A_min", cfg, 0.)
            )
        ),
        contractility_activity(
            std::make_pair(
                get_as<double>("contractility_activation", cfg, 1.),
                get_as<double>("contractility_deactivation", cfg, 0.)
            )
        ),
        num_repeat(get_as<std::size_t>("num_repeat", cfg, 1)),
        jiggle_tolerance(get_as<double>("jiggle_tolerance", cfg, tolerance)),
        jiggle_intensity(get_as<double>("jiggle_intensity", cfg, 0.))
    {
        if (num_repeat == 0) {
            throw Utopia::KeyError("num_repeat", cfg, fmt::format(
                "Value must be larger than 0, but was {}", num_repeat));
        }
        if (num_repeat > 1 and jiggle_tolerance < tolerance) {
            throw Utopia::KeyError("jiggle_tolerance", cfg, fmt::format(
                "Value must be larger or equal to 'tolerance', but was {} < {}",
                jiggle_tolerance, tolerance));
        }

        if (    temperature.param().stddev() > 1.e-12
            and update_scheme != UpdateScheme::SteepestGradient)
        {
            throw Utopia::KeyError("temperature", cfg, fmt::format("Random "
                "brownian motion from temperature ({} > 0) is only allowed in "
                "fixed stepsize update scheme, such as `SteepestGradient`. "
                "Either set temperature to 0 or select a fixed stepsize update "
                "scheme (selected scheme: {}).", temperature.param().stddev(),
                get_update_scheme(update_scheme)));
        }
    }

    /// Initialize from config and inherit not defined values from default
    template <typename Config>
    MinimizationParams(const Config& cfg, const MinimizationParams& defaults)
    :
        tolerance(get_as<double>("tolerance", cfg, defaults.tolerance)),
        update_scheme(
            setup_update_scheme(cfg,
                                get_update_scheme(defaults.update_scheme))),
        dt(get_as<double>("dt", cfg, defaults.dt)),
        max_steps(get_as<std::size_t>("max_steps", cfg, defaults.max_steps)),
        num_steps(get_as<std::size_t>("num_steps", cfg, defaults.num_steps)),
        temperature(0., 0.),
        linetension_fluctuations(
            std::make_pair(
                get_as<double>("linetension_fluctuation_tau", cfg,
                    std::get<0>(defaults.linetension_fluctuations)),
                get_as<double>("linetension_fluctuation", cfg,
                    std::get<1>(defaults.linetension_fluctuations))
            )
        ),
        area_fluctuations(
            std::make_tuple(
                get_as<double>("area_fluctuation_tau", cfg,
                    std::get<0>(defaults.area_fluctuations)),
                get_as<double>("area_fluctuation", cfg,
                    std::get<1>(defaults.area_fluctuations)),
                get_as<double>("area_fluctuation_A_min", cfg,
                    std::get<2>(defaults.area_fluctuations))
            )
        ),
        contractility_activity(
            std::make_pair(
                get_as<double>("contractility_activation", cfg,
                    std::get<0>(defaults.contractility_activity)),
                get_as<double>("contractility_deactivation", cfg,
                    std::get<1>(defaults.contractility_activity))
            )
        ),
        num_repeat(get_as<std::size_t>("num_repeat", cfg, defaults.num_repeat)),
        jiggle_tolerance(get_as<double>("jiggle_tolerance", cfg,
                                        defaults.jiggle_tolerance)),
        jiggle_intensity(get_as<double>("jiggle_intensity", cfg,
                                        defaults.jiggle_intensity))
    {
        if (cfg["temperature"]) {
            temperature = std::normal_distribution<double>(0.,
                get_as<double>("temperature", cfg));
        }
        else {
            temperature = defaults.temperature;
        }
        
        if (num_repeat == 0) {
            throw Utopia::KeyError("num_repeat", cfg, fmt::format(
                "Value must be larger than 0, but was {}", num_repeat));
        }
        if (jiggle_tolerance < tolerance) {
            throw Utopia::KeyError("jiggle_tolerance", cfg, fmt::format(
                "Value must be larger or equal to 'tolerance', but was {} < {}",
                jiggle_tolerance, tolerance));
        }
    }

    /// Setup function for the update scheme
    /** Currently implemented update schemes:
     *      -# steepest_gradient : Steepest gradient update at fixed step size
     *      -# steepest_gradient_adaptive : Steepest gradient update at adaptive
     *              step size. Step size is to next minimum of energy in the 
     *              direction of steepest gradient
     *      -# conjugate_gradient : Conjugate gradient update method
     */
    template <typename Config>
    UpdateScheme setup_update_scheme(const Config& cfg,
                                     std::string default_scheme = "") {
        const auto update_scheme = get_as<std::string>("update_scheme", cfg,
                                                       default_scheme);

        if (update_scheme == "steepest_gradient") {
            return SteepestGradient;
        }
        if (update_scheme == "steepest_gradient_adaptive") {
            return SteepestGradientAdaptive;
        }
        if (update_scheme == "conjugate_gradient") {
            return ConjugateGradient;
        }

        throw KeyError("update_scheme", cfg, 
            "Update scheme must be one of the following: "
                "'steepest_gradient', "
                "'steepest_gradient_adaptive', "
                "'conjugate_gradient'.");
    }
    
    /// Transform the update scheme into a human readable string
    std::string get_update_scheme(UpdateScheme scheme) const {
        if (scheme == SteepestGradient) {
            return "steepest_gradient";
        }
        if (scheme == SteepestGradientAdaptive) {
            return "steepest_gradient_adaptive";
        }
        if (scheme == ConjugateGradient) {
            return "conjugate_gradient";
        }
        else {
            throw std::runtime_error(fmt::format("Unknown update scheme {}!",
                                                 scheme));
        }
    }
};

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

    /// The types of a cell
    using CellType = typename CellState::CellType;

    using OrderedEdgeContainer = typename AgentManager::OrderedEdgeContainer;

    /// The type of a rule function acting on vertices of the agent manager
    using RuleFuncVertex = typename AgentManager::RuleFuncVertex;

    /// The type of a rule function acting on edges of the agent manager
    using RuleFuncEdge = typename AgentManager::RuleFuncEdge;

    /// The type of a rule function acting on cells of the agent manager
    using RuleFuncCell = typename AgentManager::RuleFuncCell;

    /// A matrix for properties depending on the state of 2 cells 
    using CellCellPropertyMatrix = 
            arma::Mat<double>::fixed<CellType::num_cell_types + 1,
                                     CellType::num_cell_types + 1>;

    using UpdateScheme = MinimizationParams::UpdateScheme;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The manager of the model's entities
    AgentManager _am;


    // -- Minimization parameters ---------------------------------------------

    // PARAMETERS
    MinimizationParams _default_minimization_params;

    /// timestep scaling
    double _dt;

    /// The update scheme for energy minimization
    /** Currently implemented update schemes:
     *      -# steepest_gradient : Steepest gradient update at fixed step size
     *      -# steepest_gradient_adaptive : Steepest gradient update at adaptive
     *              step size. Step size is to next minimum of energy in the 
     *              direction of steepest gradient
     *      -# conjugate_gradient : Conjugate gradient update method
     */
    UpdateScheme _update_scheme;
    
    /// The tolerance during minimization
    double _minimization_tolerance;

    /// The temperature for random brownian motion
    /** With mean 0 and variance \f$ \sigma^2 = 2 T \f$.
     */
    std::normal_distribution<double> _distr_temperature;

    /// The timescale and amplitude of Ornstein-Uhlenbeck fluctuations on
    /// linetension
    std::pair<double, double> _linetension_fluctuations;

    /// The activation and deactivation rate of edge contractility
    std::pair<double, double> _contractility_activity;

    /// Whether edge contractility only aplies to apical boundary edges
    bool _apical_contractility;

    /// The timescale, amplitude, and lower area limit of Ornstein-Uhlenbeck
    /// fluctuations on preferential area using a lognormal distribution
    std::tuple<double, double, double> _area_fluctuations;


    // -- Mechanical parameters -----------------------------------------------

    /// Linetension constant Lambda
    /** The entries are linetensions at interfaces between two cells of types i 
     *  and j.
     *  \note The values of this matrix are used for new edges, e.g. in T1
     *        transitions. There is no other generic way to determine the
     *        surface tension between two cells from the previous edge.
     *  \note This is a symmetric matrix
     */
    CellCellPropertyMatrix _linetension;

    /// Linetension constant Lambda
    /** The entries are contractility at interfaces between two cells of types
     *  i and j
     *  \note The values of this matrix are used for new edges, e.g. in T1
     *        transitions. There is no other generic way to determine the
     *        surface tension between two cells from the previous edge.
     *  \note This is a symmetric matrix
     */
    CellCellPropertyMatrix _edge_contractility;
    
    /// Area elasticity constant K
    double _area_elasticity;

    /// Cell contractility implementation
    enum CellContractility {
        /// Normalize shape index to target cell area
        Contractile,

        /// Normalize shape index to cell area
        Shape_elastic
    } _cell_contractility_impl;

    /// Mechanical parameters for boundary
    /** The boundary of the domain is treated as one cell
     */
    struct BoundaryParam {
        double area_elasticity;
        double _area_preferential; /// Per cell averaged target area
        double area_preferential(std::size_t N) const {
            return _area_preferential * static_cast<double>(N);
        };

        double contractility;
        double shape_index_preferential;

        // -- A quadratic boundary potential in the shape of a stripe of
        //    a circle with radius R and width W ------------------------------

        /// The constant for the quadratic stripe potential
        double stripe_potential_constant;

        /// The width of the stripe
        double stripe_width;

        /// The curvature of the circle
        double stripe_curvature;

        /// The relative coordinate in horiz. axis where the place the origin
        double stripe_curvature_center;

        /// the center of the circle stripe
        /** originally the tissue center offsetted by 1. / curvature
         *  NOTE it is fixed and updated with changes in curvature to
         *       prevent macroscopic cell flows when always defining wrt 
         *       cell center
         */
        std::shared_ptr<SpaceVec> stripe_origin = nullptr;

        /// Whether to fix all vertices of the boundary in space
        bool fix_boundary;

        BoundaryParam (const Config& cfg)
        :
            area_elasticity(get_as<double>("area_elasticity", cfg)),
            _area_preferential(get_as<double>("area_preferential", cfg)),
            contractility(get_as<double>("contractility", cfg)),
            shape_index_preferential(get_as<double>("shape_index_preferential",
                                                    cfg)),
            
            // NOTE use operations to update
            stripe_potential_constant(0.),
            stripe_width(0.),
            stripe_curvature(0.),

            fix_boundary(get_as<bool>("fix_boundary", cfg, false))
        { }

    } _boundary_param;


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

    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::normal_distribution<double> _normal_distr;

    // .. Temporary objects ...................................................
protected:
    /// The total energy in the last step
    double _energy_previous_step;

    /// Current energy
    double _energy;

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
        _update_scheme(_default_minimization_params.update_scheme),
        _minimization_tolerance(_default_minimization_params.tolerance),
        _distr_temperature(_default_minimization_params.temperature),
        _linetension_fluctuations(
            _default_minimization_params.linetension_fluctuations),
        _contractility_activity(
            _default_minimization_params.contractility_activity),
        _apical_contractility(
            get_as<bool>("apical_edge_contractility", this->_cfg, false)),
        _area_fluctuations(
            _default_minimization_params.area_fluctuations),
        _linetension(this->setup_linetension(this->_cfg)),
        _edge_contractility(this->setup_edge_contractility(this->_cfg)),
        _area_elasticity(get_as<double>("area_elasticity", this->_cfg)),
        _cell_contractility_impl(setup_cell_contractility_impl(
            get_as<std::string>("cell_shape_implementation", this->_cfg)
        )),
        _boundary_param(get_as<Config>("boundary_parameter", this->_cfg)),
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
        _energy_previous_step(std::numeric_limits<double>::max()),
        _energy(0.),
        _num_T1s(0),
        _num_T1s_total(0),
        _num_T1s_attempted(0),
        _num_T1s_attempted_total(0),
        _num_T2s(0),
        _num_T2s_total(0),
        _num_minimizations(0)
    {
        double initial_jiggle(get_as<double>("initial_jiggle", this->_cfg, 0.));
        if (initial_jiggle > 1.e-12) {
            jiggle_vertices(initial_jiggle);
        }
        
        this->_log->info("Model initialized.");
    }


private:
    // .. Setup functions .....................................................
    /// Setup up the linetension from config
    CellCellPropertyMatrix setup_linetension(const Config& cfg)
    {
        const auto edge_cfg = cfg["agent_manager"]["edge_manager"];
        const double linetension = get_as<double>("linetension",
                                                  edge_cfg["agent_params"]);
        // NOTE not checking paths, because happened in constructor of _am
        
        CellCellPropertyMatrix matrix;
        return matrix.fill(linetension);
    }

    /// Setup up the edge contractility from config
    CellCellPropertyMatrix setup_edge_contractility(const Config& cfg)
    {
        const auto edge_cfg = cfg["agent_manager"]["edge_manager"];
        const double contractility = get_as<double>("contractility",
                                                    edge_cfg["agent_params"]);
        // NOTE not checking paths, because happened in constructor of _am
        
        CellCellPropertyMatrix matrix;
        return matrix.fill(contractility);
    }

    CellContractility setup_cell_contractility_impl (const std::string& impl)
    {
        if (impl == "contractile") {
            return CellContractility::Contractile;
        }
        else if (impl == "shape_elastic") {
            return CellContractility::Shape_elastic;
        }
        else {
            throw std::invalid_argument(fmt::format(
                "While setting up cell contractility implementation, received "
                "invalid configuration. {} is not a valid implementation. "
                "Choose one of the following: 'contractile', 'shape_elastic'.",
                impl
            ));
        }
    }


    // ..energy terms .........................................................
    // See energy.hh for implementation

    double line_tension_energy (const std::shared_ptr<Edge>& edge,
                                double beta = 0.) const;
    double edge_contractility_energy (const std::shared_ptr<Edge>& edge,
                                      double beta = 0.) const;
    double area_elasticity_energy (const std::shared_ptr<Cell>& cell,
                                   double beta = 0.) const;
    double cell_contractility_energy (const std::shared_ptr<Cell>& cell,
                                      double beta = 0.) const;

    // .. Force setter functions ..............................................
    /// Resets the forces of this vertex
    const RuleFuncVertex reset_forces = [] (const auto& vertex)
    {
        auto state = vertex->state;
        state.f.zeros();
        state.virtual_pos = std::make_pair(0, state.f);
        return state;
    };

    /** Calculates the forces from linetension
     * 
     *  Contractive force of the edge where energy is proportional to the edge's
     *  length.
     * 
     *  @param  e   Pointer to the edge for which to calculate the forces
     *              NOTE that forces only act on vertices
     * 
     *  \return energy associated with this edge
     */
    const RuleFuncEdge set_grad_linetension = [this](const auto& edge) {
        if (fabs(edge->state.linetension()) < 1.e-12) {
            return edge->state;
        }

        auto a = edge->custom_links().a;
        auto b = edge->custom_links().b;
        
        SpaceVec displ = this->_am.displacement(a, b);
        auto length = arma::norm(displ);

        SpaceVec force = edge->state.linetension() * displ / length;

        a->state.f += force;
        b->state.f -= force;

        return edge->state;
    };

    /** Calculates the forces from edge contractility
     * 
     *  Contractive force of the edge where energy is proportional to the edge's
     *  length.
     * 
     *  @param  e   Pointer to the edge for which to calculate the forces
     *              NOTE that forces only act on vertices
     * 
     *  \return energy associated with this edge
     */
    const RuleFuncEdge set_grad_edge_contractility = [this](const auto& edge) {
        if (fabs(edge->state.contractility()) < 1.e-12) {
            return edge->state;
        }

        auto a = edge->custom_links().a;
        auto b = edge->custom_links().b;

        SpaceVec force = edge->state.contractility() *
                         this->_am.displacement(a, b);

        a->state.f += force;
        b->state.f -= force;

        return edge->state;
    };

    /** Calculates the forces from area elasticity
     * 
     *  Response force to a deformation in area where the energy is 
     *  K/2 * (A - A0)**2
     * 
     *  @param c    Pointer to the cell for which to calculate the forces
     *              NOTE that forces only act on vertices
     * 
     *  \return energy associated with this edge
     */
    const RuleFuncCell set_grad_area_elasticity = [this](const auto& cell) {
        const auto state = cell->state;

        const auto rel_cell_area = (  this->_am.area_of(cell)
                                    / state.area_preferential());
        const SpaceVec cell_center = this->_am.barycenter_of(cell);
        
        const auto& edges = cell->custom_links().edges;
        for (unsigned int edges_it = 0; edges_it < edges.size(); edges_it++) {
            std::shared_ptr<Edge> e0; bool e0_flip;
            if (edges_it > 0) { 
                std::tie(e0, e0_flip) = edges[edges_it - 1];
            }
            else {
                std::tie(e0, e0_flip) = edges.back();
            }

            const auto [e1, e1_flip] = edges[edges_it];
            
            // vertices in ordering
            auto v_center = e0->custom_links().b;
            auto v_prior = e0->custom_links().a;
            if (e0_flip) {
                std::swap(v_center, v_prior);
            }

            std::shared_ptr<Vertex> v_post;
            if (not e1_flip) {
                v_post = e1->custom_links().b;
            }
            else {
                v_post = e1->custom_links().a;
            }

            // get positions relative to cell center
            SpaceVec prior = cell_center + 
                             this->_space->displacement(cell_center,
                                                    _am.position_of(v_prior));
            SpaceVec post = cell_center + 
                            this->_space->displacement(cell_center,
                                                    _am.position_of(v_post));

            SpaceVec displ = post - prior;

            // derivative of A to x_i, i.e. the position of v_center 
            SpaceVec dA_dx({0.5 * displ[1], -0.5 * displ[0]});

            SpaceVec force = (  -1. * this->_area_elasticity
                              * (rel_cell_area - 1) * dA_dx
                              / state.area_preferential());

            v_center->state.f += force;
        }
        
        return state;
    };

    const RuleFuncCell set_grad_shape_elasticity = [this](const auto& cell)
    {
        auto state = cell->state;

        if (fabs(state.contractility) < 1.e-12) {
            return state;
        }

        const double area = _am.area_of(cell);
        const double perimeter = this->_am.perimeter_of(cell);
        const double shape_index = perimeter / sqrt(area);

        const SpaceVec cell_center = this->_am.barycenter_of(cell);
        
        const auto& edges = cell->custom_links().edges;
        for (unsigned int edges_it = 0; edges_it < edges.size(); edges_it++) {
            std::shared_ptr<Edge> e0; bool e0_flip;
            if (edges_it > 0) { 
                std::tie(e0, e0_flip) = edges[edges_it - 1];
            }
            else {
                std::tie(e0, e0_flip) = edges.back();
            }

            const auto [e1, e1_flip] = edges[edges_it];
            
            // vertices in ordering
            auto v_center = e0->custom_links().b;
            auto v_prior = e0->custom_links().a;
            if (e0_flip) {
                std::swap(v_center, v_prior);
            }

            std::shared_ptr<Vertex> v_post;
            if (not e1_flip) {
                v_post = e1->custom_links().b;
            }
            else {
                v_post = e1->custom_links().a;
            }

            // get positions relative to cell center
            SpaceVec center = cell_center + 
                              this->_space->displacement(cell_center,
                                                    _am.position_of(v_center));
            SpaceVec prior = cell_center + 
                             this->_space->displacement(cell_center,
                                                    _am.position_of(v_prior));
            SpaceVec post = cell_center + 
                            this->_space->displacement(cell_center,
                                                    _am.position_of(v_post));

            // calculate dA / dx_i
            SpaceVec curv = post - prior;
            SpaceVec dA_dx({0.5 * curv[1], -0.5 * curv[0]});

            // calculate dP / dx_i
            SpaceVec displ_e_0 = this->_space->displacement(prior, center);
            SpaceVec displ_e_1 = this->_space->displacement(center, post);
            SpaceVec dP_dx = (  displ_e_0 / arma::norm(displ_e_0)
                              - displ_e_1 / arma::norm(displ_e_1));

            SpaceVec dE_dx = (  state.contractility
                              * (shape_index - state.shape_index_preferential)
                              * (  dP_dx / sqrt(area)
                                 - 0.5 * dA_dx * shape_index / area));

            v_center->state.f -= dE_dx;
        }

        return state;
    };

    /// Constriction on the cell's shape index
    const RuleFuncCell set_grad_cell_contractility = [this](const auto& cell)
    {
        auto state = cell->state;

        if (fabs(state.contractility) < 1.e-12) {
            return state;
        }

        double perimeter = this->_am.perimeter_of(cell);

        for (auto [e, flip] : cell->custom_links().edges) {
            auto a = e->custom_links().a;
            auto b = e->custom_links().b;
            if (flip) { std::swap(a, b); }

            SpaceVec displ = this->_am.displacement(a, b);
            double length = arma::norm(displ);

            SpaceVec force = (  state.contractility
                              * (  perimeter / sqrt(state.area_preferential())
                                 - state.shape_index_preferential)
                              * displ / length / sqrt(state.area_preferential())
                             );

            a->state.f += force;
            b->state.f -= force;
        }

        return state;
    };

    /// Apply a torque to the vertices of the cell
    const RuleFuncCell set_grad_torque = [this](const auto& cell) {
        if (not cell->custom_links().rotation_state) {
            return cell->state;
        }

        const SpaceVec cell_center = _am.barycenter_of(cell);
        double torque = cell->custom_links().rotation_state->torque;

        for (const auto& vertex : cell->custom_links().vertices) {
            const SpaceVec pos = _am.position_of(vertex);
            SpaceVec displ = _space->displacement(pos, cell_center);
            double distance = arma::norm(displ);

            SpaceVec force = torque / distance * SpaceVec({-displ[1],
                                                           displ[0]});

            vertex->state.f += force;
        }

        return cell->state;
    };



    void set_grad_boundary_area_elasticity
            (const OrderedEdgeContainer& boundary)
    {
        if (   _space->periodic
            or fabs(_boundary_param.area_elasticity) < 1.e-12)
        {
            return;
        }

        double area = _am.area_of(boundary);
        double N = _am.cells().size();

        const auto rel_area = area / _boundary_param.area_preferential(N);
        
        for (unsigned int edges_it = 0; edges_it < boundary.size(); edges_it++) {
            std::shared_ptr<Edge> e0; bool e0_flip;
            if (edges_it > 0) { 
                std::tie(e0, e0_flip) = boundary[edges_it - 1];
            }
            else {
                std::tie(e0, e0_flip) = boundary.back();
            }

            const auto [e1, e1_flip] = boundary[edges_it];
            
            // vertices in ordering
            auto v_center = e0->custom_links().b;
            auto v_prior  = e0->custom_links().a;
            if (e0_flip) {
                std::swap(v_center, v_prior);
            }

            std::shared_ptr<Vertex> v_post;
            if (not e1_flip) {
                v_post = e1->custom_links().b;
            }
            else {
                v_post = e1->custom_links().a;
            }

            // get positions relative to cell center
            SpaceVec prior = _am.position_of(v_prior);
            SpaceVec post  = _am.position_of(v_post);

            SpaceVec displ = post - prior;

            // derivative of A to x_i, i.e. the position of v_center 
            SpaceVec dA_dx({0.5 * displ[1], -0.5 * displ[0]});

            SpaceVec force = (  -1. * _boundary_param.area_elasticity
                              * (rel_area - 1) * dA_dx
                              / _boundary_param.area_preferential(N));

            v_center->state.f += force;
        }
        
        return;
    };

    void set_grad_boundary_shape_elasticity 
            (const OrderedEdgeContainer& boundary)
    {
        if (   _space->periodic
            or fabs(_boundary_param.contractility) < 1.e-12)
        {
            return;
        }

        const double area = _am.area_of(boundary, 0.);
        const double perimeter = this->_am.perimeter_of(boundary, 0.);
        const double shape_index = perimeter / sqrt(area);
        
        for (unsigned int edges_it = 0; edges_it < boundary.size(); edges_it++)
        {
            std::shared_ptr<Edge> e0; bool e0_flip;
            if (edges_it > 0) { 
                std::tie(e0, e0_flip) = boundary[edges_it - 1];
            }
            else {
                std::tie(e0, e0_flip) = boundary.back();
            }

            const auto [e1, e1_flip] = boundary[edges_it];
            
            // vertices in ordering
            auto v_center = e0->custom_links().b;
            auto v_prior  = e0->custom_links().a;
            if (e0_flip) {
                std::swap(v_center, v_prior);
            }

            std::shared_ptr<Vertex> v_post;
            if (not e1_flip) {
                v_post = e1->custom_links().b;
            }
            else {
                v_post = e1->custom_links().a;
            }

            // get positions relative to cell center
            SpaceVec center = _am.position_of(v_center);
            SpaceVec prior  = _am.position_of(v_prior);
            SpaceVec post   = _am.position_of(v_post);

            // calculate dA / dx_i
            SpaceVec displ = post - prior;
            SpaceVec dA_dx({0.5 * displ[1], -0.5 * displ[0]});

            // calculate dP / dx_i
            SpaceVec displ_2 = this->_space->displacement(prior, center);
            SpaceVec displ_3 = this->_space->displacement(center, post);
            SpaceVec dP_dx = (  displ_2 / arma::norm(displ_2)
                              - displ_3 / arma::norm(displ_3));

            SpaceVec dE_dx = (  _boundary_param.contractility
                              * (  shape_index
                                 - _boundary_param.shape_index_preferential)
                              * (  dP_dx / sqrt(area)
                                 - 0.5 * dA_dx * shape_index / area));

            v_center->state.f -= dE_dx;
        }

        return;
    };

    void set_grad_boundary_contractility
            (const OrderedEdgeContainer& boundary)
    {
        if (   _space->periodic
            or fabs(_boundary_param.contractility) < 1.e-12)
        {
            return;
        }

        const double perimeter = this->_am.perimeter_of(boundary, 0.);
        const auto N = _am.cells().size();
        
        for (auto [e, flip] : boundary) {
            auto a = e->custom_links().a;
            auto b = e->custom_links().b;
            if (flip) { std::swap(a, b); }

            SpaceVec displ = this->_am.displacement(a, b);
            double length = arma::norm(displ);

            SpaceVec force = (  _boundary_param.contractility
                              * (    perimeter
                                   / sqrt(_boundary_param.area_preferential(N))
                                 - _boundary_param.shape_index_preferential)
                              * displ / length
                              / sqrt(_boundary_param.area_preferential(N))
                             );

            a->state.f += force;
            b->state.f -= force;
        }

        return;
    };

    /// Derivative of a quadratic boundary potential
    void set_grad_boundary_stripe () {
        if (   _space->periodic
            or fabs(_boundary_param.stripe_potential_constant) < 1.e-12)
        {
            return;
        }

        double curvature = std::max(_boundary_param.stripe_curvature, 1.e-10);
        double R = 1. / curvature;

        // the (fixed) center of the circle stripe
        SpaceVec origin = *_boundary_param.stripe_origin - SpaceVec({0., R});

        double inner_radius = (R - _boundary_param.stripe_width / 2.);
        double outer_radius = (R + _boundary_param.stripe_width / 2.);
        
        // apply to all vertices outside the domain
        for (const auto &v : _am.vertices()) {
            // the position wrt origin
            SpaceVec pos = _am.position_of(v) - origin;

            double radius = arma::norm(pos);
            if (radius < inner_radius) {
                v->state.f += fabs(radius - inner_radius) * pos / radius;
            }
            else if (radius > outer_radius) {
                v->state.f -= fabs(radius - outer_radius) * pos / radius;
            }
        }

        return;
    };


    /// Set the gradient of energy within the vertices
    /** \details This function applies the forces arising from the different
     *           energy terms.
     *  \note    When adding energy terms, remember to add their gradient here!
     */
    void set_gradient () {
        // reset forces
        apply_rule<Update::sync>(reset_forces, _am.vertices());

        // apply new forces
        apply_rule<Update::async, Shuffle::off>(set_grad_linetension,
                                                _am.edges());
        apply_rule<Update::async, Shuffle::off>(set_grad_edge_contractility,
                                                _am.edges());
        apply_rule<Update::async, Shuffle::off>(set_grad_area_elasticity,
                                                _am.cells());
        if (this->_cell_contractility_impl == Contractile) {
            apply_rule<Update::async, Shuffle::off>(set_grad_cell_contractility,
                                                    _am.cells());
        }
        else if (this->_cell_contractility_impl == Shape_elastic) {
            apply_rule<Update::async, Shuffle::off>(set_grad_shape_elasticity,
                                                    _am.cells());
        }
        else {
            throw std::runtime_error(fmt::format(
                "Not Implemented Error: "
                "Cell Contractility Implementation {}!",
                this->_cell_contractility_impl));
        }

        // apply boundary forces
        if (not _space->periodic) {
            const auto boundary = _am.get_boundary_edges();
            set_grad_boundary_area_elasticity(boundary);
            if (this->_cell_contractility_impl == Contractile) {
                set_grad_boundary_contractility(boundary);
            }
            else if (this->_cell_contractility_impl == Shape_elastic) {
                set_grad_boundary_shape_elasticity(boundary);
            }
            else {
                throw std::runtime_error(fmt::format(
                    "Not Implemented Error: "
                    "Cell Contractility Implementation {}!",
                    this->_cell_contractility_impl));
            }
            set_grad_boundary_stripe();
        }

        if (_update_scheme == UpdateScheme::SteepestGradient) {
            apply_rule<Update::async, Shuffle::off>(set_grad_torque, 
                                                    _am.cells());
        }

        // fix the boundary
        if (not _space->periodic and _boundary_param.fix_boundary) {
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
        // set forces on fixed vertices (e.g. boundary) to zero
        // NOTE this is always done, as particular vertices can be manually
        //      fixed 
        apply_rule<Update::sync>(
            [](const auto& vertex) {
                auto state = vertex->state;
                if (state.fix_in_space) {  state.f = SpaceVec({0., 0.}); }
                return state;
            },
            _am.vertices()
        );

        // reset the virtual position
        apply_rule<Update::sync>(
            [this](const auto& vertex) {
                vertex->state.virtual_pos = 
                    std::make_pair(0., _am.position_of(vertex));
                return vertex->state;
            },
            _am.vertices()
        );
    }

    /** The update of position
     * 
     *  Move vertex proportional to the gradient of energy (force)
     * 
     *  @param v    The pointer to the vertex to update
     */
    const RuleFuncVertex update_position = [this](const auto& vertex) {
        _am.move_by(vertex, vertex->state.f * this->_dt);
        return vertex->state;
    };

    const RuleFuncVertex update_brownian_motion = [this](const auto& vertex) {
        SpaceVec temp = {_distr_temperature(*this->_rng),
                         _distr_temperature(*this->_rng)};
        _am.move_by(vertex, temp * this->_dt);
        
        return vertex->state;
    };

    /// Update linetension fluctuation in an Ornstein-Uhlenbeck process
    const RuleFuncEdge update_linetension_ornstein =
    [this](const auto& edge)
    {
        double linetension = edge->state._linetension_fluctuation;

        auto [tau, dL] = this->_linetension_fluctuations;
        
        double rand_l = dL * sqrt(2. * _dt / tau) * _normal_distr(*this->_rng);

        linetension += rand_l - _dt / tau * linetension;

        edge->state._linetension_fluctuation = linetension;
        
        return edge->state;
    };

    const RuleFuncEdge update_edge_contractility =
    [this](const auto& edge)
    {
        auto state = edge->state;

        const auto [act, deact] = this->_contractility_activity;

        if (not state.contractility_on) {
            state.contractility_on = (_prob_distr(*this->_rng) < act);
        }
        else {
            state.contractility_on = (_prob_distr(*this->_rng) > deact);
        }

        return state;
    };

    /// Update area preferential fluctuation in an Ornstein-Uhlenbeck process
    /** \note A^(0) > 0 required. Hence using a lognormal distribution enforcing
     *        A^(0) > A_min
     */
    const RuleFuncCell update_area_preferential_ornstein =
    [this](const auto& cell)
    {
        double A_fluc = cell->state._area_preferential_fluctuations;
        double A_targ = cell->state._area_preferential;

        // the timescale, relative fluctuation amplitude, and minimum area
        auto [tau, dA, A_min] = this->_area_fluctuations;
        
        // shift lognormal distribution to account for A_min
        double mean_distr = cell->state.area_preferential() - A_min;

        auto distr = get_lognormal_distribution(mean_distr, dA * A_targ);
        double rn = distr(*this->_rng) - mean_distr;
        // NOTE the rns have mean 0, min A_min, and stddev dA * A_targ

        // the Ornstein Uhlenbeck step for fluctuations
        A_fluc += sqrt(2. * _dt / tau) * rn - _dt / tau * A_fluc;

        cell->state._area_preferential_fluctuations = A_fluc;
        
        return cell->state;
    };
    

    // -- The algorithm    ----------------------------------------------------
    // see algorithm.hh
    std::tuple<bool, double, double> determine_timestep (
        double dt,
        const double energy_0,
        const double tolerance) const;
    double steepest_gradient_step (bool adaptive_step);
    double conjugate_gradient_step ();
    double perform_update_step(UpdateScheme update_scheme);
    bool perform_transitions(bool enabled);


    // -- Helper functions ----------------------------------------------------

    /// Update the rotation tracking of a cell
    /** Tracks the average angular velocity of the cell's vertices
     */
    const RuleFuncCell track_rotation = [this](const auto& cell)
    {
        if (not cell->custom_links().rotation_state) {
            return cell->state;
        }

        const auto& rot_state = cell->custom_links().rotation_state;

        double angular_vel = rot_state->angular_velocity(cell, _am);
        angular_vel *= _dt;

        double tracked_rotation = rot_state->tracked_rotation;
        double alpha = rot_state->tracking_persistence;

        rot_state->tracked_rotation += (  angular_vel
                                        - tracked_rotation * alpha);

        return cell->state;
    };

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
    void divide_cell(std::shared_ptr<Cell> cell, double division_angle) {
        return _am.divide_cell(cell, division_angle, _linetension,
                               _edge_contractility);
    }

    void increase_domain_size(double area);
    double stretch_domain(SpaceVec stretch, bool compensate,
        bool fix_hc_area, bool fix_sc_area);


    // .. Simulation Control ..................................................
    void init_minimization ();

    /// Iterate a single step
    /** \details Rules applied
     *      -# perform_transitions()
     *      -# perform_update_step()
     *      -# tracking of variables
     */
    void perform_step () {
        bool transition_occurred = perform_transitions(_enable_transitions);

        if (transition_occurred) {
            // restart the conjugate gradient update
            this->init_minimization();
        }

        _energy_previous_step = _energy;
        _energy = perform_update_step(_update_scheme);

        apply_rule<Update::sync>(track_rotation, _am.cells());

        this->_log->trace("Energy changed by {}",
                          get_rel_energy_change(_energy,
                                                _energy_previous_step));
    }

    void prolog () {
        this->init_minimization();
        return this->__prolog();
    }
    
    /// Monitor model information
    void monitor () {
        this->_monitor.set_entry("time", this->_time);
        this->_monitor.set_entry("energy", _energy);
        this->_monitor.set_entry("energy_change",
                                 get_rel_energy_change(_energy,
                                                       _energy_previous_step));
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
                          "perturbation): {:7d} / {:d}",
                          this->_time, this->_time_max);

        const auto time_0 = this->get_time();
        this->_log->debug("Minimizing energy from step {} with {} repeats ...",
                          time_0, params.num_repeat);
        _status = Status::Minimization;


        _update_scheme = params.update_scheme;
        _dt = params.dt;
        _distr_temperature = params.temperature;
        _linetension_fluctuations = params.linetension_fluctuations;
        _area_fluctuations = params.area_fluctuations;
        if (std::get<1>(_linetension_fluctuations) < 1.e-12) {
            for (const auto& e : _am.edges()) {
                e->state._linetension_fluctuation = 0.;
            }
        }
        if (std::get<1>(_area_fluctuations) < 1.e-12) {
            for (const auto& c : _am.cells()) {
                c->state._area_preferential_fluctuations = 0.;
            }
        }
        _contractility_activity = params.contractility_activity;
        if (std::get<1>(_contractility_activity) < 1.e-10) {
            for (const auto& e : _am.edges()) {
                e->state.contractility_on = true;
            }
        }


        const double energy_after_perturbation = this->get_energy();

        for (std::size_t i = 0; i < params.num_repeat; i++)
        {
            // jiggle vertices if required
            if (params.jiggle_intensity > 1.e-12) {
                _status = Status::Jiggled;
                this->jiggle_vertices(params.jiggle_intensity);
                this->increment_time();
                this->_datamanager(*this);            
                this->_log->debug("Incremented time after jiggling: "
                                  "{:7d} / {:d}",
                                  this->_time, this->_time_max);
                
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
                // disable topological transitions in first iteration
                if (this->get_time() - time_start == 0) {
                    _enable_transitions = false; }
                else { _enable_transitions = tmp_enable_transitions; }

                this->iterate();
                monitor_mngr();

                minimum_reached = equilibrium_condition();
                double energy_change = get_rel_energy_change(
                    _energy, _energy_previous_step);
                
                if (not minimum_reached 
                    and this->get_time() - time_start >= params.max_steps)
                {
                    throw std::runtime_error(fmt::format(
                        "Equilibrium not reached within {} steps at a "
                        "tolerance of {}! Relative energy change in last step "
                        "was {}.",
                        params.max_steps, tolerance, energy_change));
                }

                if (stop_now.load()) {
                    this->_log->warn("Was told to stop. Not iterating "
                        "further ...");
                    throw GotSignal(received_signum.load());
                }
            }
            this->_log->trace("  Energy minimized in {} steps.",
                              this->get_time() - time_start);

            _num_minimizations++;

            _enable_transitions = tmp_enable_transitions;
            // NOTE restore initial state of transitions allowed
        }

        auto num_steps = this->get_time() - time_0;
        if (get_rel_energy_change(_energy, energy_after_perturbation) > -1.e-3)
        {
            this->_log->info("WARN Energy minimized in {} steps and changed by "
                "{} (rel.).",
                num_steps,
                get_rel_energy_change(_energy, energy_after_perturbation));
        }
        else {
            this->_log->info("Energy minimized in {} steps and changed by "
                "{} (rel.).",
                num_steps,
                get_rel_energy_change(_energy, energy_after_perturbation));
        }

        return num_steps;
    }

    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

protected:

    // .. Prediction of energy terms ..........................................
    // see energy.hh for implementation
    // NOTE when adding energy terms remember to add them to get_energy(..)

    double get_energy_linetension(const AgentContainer<Edge>& es,
                                  double beta) const;
    double get_energy_edge_contractility(const AgentContainer<Edge>& es,
                                         double beta) const;
    double get_energy_areaelasticity (const AgentContainer<Cell>& cs,
                                      double beta) const;    
    double get_energy_cell_contractility (const AgentContainer<Cell>& cs,
                                          double beta) const;
                                          
    double get_boundary_area_energy (double beta) const;
    double get_boundary_shape_energy(double beta) const;
    double get_boundary_stripe_energy(double beta) const;

    double get_energy(const AgentContainer<Edge>& es,
                      const AgentContainer<Cell>& cs,
                      double beta = 0.) const;
    
    // .. Derived and global energy term predictions  .........................

    /// Predict the energy from linetension term
    /** by moving the vertices towards energy minimum at stepsize beta.
     *  This uses steepest gradient direction or conjugate gradient direction
     *  depending on update scheme.
     * 
     *  \warning The gradient is not updated! Must be done manually.
     */    
    double get_energy_linetension(double beta) const {
        return get_energy_linetension(_am.edges(), beta);
    }
    
    /// Predict the energy from edge contractility
    /** by moving the vertices towards energy minimum at stepsize beta.
     *  This uses steepest gradient direction or conjugate gradient direction
     *  depending on update scheme.
     * 
     *  \warning The gradient is not updated! Must be done manually.
     */        
    double get_energy_edge_contractility(double beta) const {
        return get_energy_edge_contractility(_am.edges(), beta);
    }
    
    /// Predict the energy from area elasticity
    /** by moving the vertices towards energy minimum at stepsize beta.
     *  This uses steepest gradient direction or conjugate gradient direction
     *  depending on update scheme.
     * 
     *  \warning The gradient is not updated! Must be done manually.
     */        
    double get_energy_areaelasticity (double beta) const {
        return get_energy_areaelasticity(_am.cells(), beta);
    }
    
    /// Predict the energy from cell contractility
    /** by moving the vertices towards energy minimum at stepsize beta.
     *  This uses steepest gradient direction or conjugate gradient direction
     *  depending on update scheme.
     * 
     *  \warning The gradient is not updated! Must be done manually.
     */
    double get_energy_cell_contractility (double beta) const {
        return get_energy_cell_contractility(_am.cells(), beta);
    }

    /// Predict the energy contribution of tbe boundary
    /** Includes shape and area elasticity
     */
    double get_boundary_energy(double beta) const {
        if (_space->periodic) {
            return 0.;
        }
        return (  get_boundary_area_energy(beta)
                + get_boundary_shape_energy(beta)
                + get_boundary_stripe_energy(beta));
    }
    
    /// Predict the energy
    /** by moving the vertices towards energy minimum at stepsize beta.
     *  This uses steepest gradient direction or conjugate gradient direction
     *  depending on update scheme.
     * 
     *  \warning The gradient is not updated! Must be done manually.
     */        
    double get_energy(double beta) const {
        return get_energy(_am.edges(), _am.cells(), beta);
    }

public:
    // .. Public energy terms .................................................

    /// Getter for energy associated with linetension
    /** Sums PCPVertex::line_tension_energy for all entities
     */
    double get_energy_linetension() const {
        return get_energy_linetension(_am.edges(), 0.);
    }
    
    /// Getter for energy associated with contractility of junctions
    /** Sums PCPVertex::edge_contractility_energy for all entities
     */
    double get_energy_edge_contractility() const {
        return get_energy_edge_contractility(_am.edges(), 0.);
    }
    
    /// Getter for energy associated with area elasticity
    /** Sums PCPVertex::area_elasticity_energy for all entities
     */
    double get_energy_areaelasticity () const {
        return get_energy_areaelasticity(_am.cells(), 0.);
    }

    /// Getter for energy associated with contractility of cells
    /** Sums PCPVertex::cell_contractility_energy for all entities
     */
    double get_energy_cell_contractility () const {
        return get_energy_cell_contractility(_am.cells(), 0.);
    }

    double get_boundary_shape_energy() const {
        return get_boundary_shape_energy(0.);
    }
    double get_boundary_area_energy () const {
        return get_boundary_area_energy(0.);
    }
    double get_boundary_stripe_energy () const {
        return get_boundary_stripe_energy(0.);
    }

    double get_boundary_energy() const {
        if (_space->periodic) {
            return 0.;
        }
        return (  get_boundary_area_energy(0.)
                + get_boundary_shape_energy(0.)
                + get_boundary_stripe_energy(0.));
    }

    /// Getter for the total energy
    double get_energy() const {
        return get_energy(_am.edges(), _am.cells(), 0.);
    }

    double get_rel_energy_change () const;
    double get_rel_energy_change (double energy, double energy_0) const;

    /// Whether the relative change in energy fulfills the equilibrium condition
    bool equilibrium_condition() const {
        if (_status != Status::Minimization) {
            return false;
        }

        double energy_change = get_rel_energy_change(_energy,
                                                     _energy_previous_step);
        return (energy_change > - _minimization_tolerance);
    }


    // .. Public energy terms for subset of entities ..........................

    /// Getter for energy associated with linetension
    /** Sums PCPVertex::line_tension_energy for provided entities
     */
    double get_energy_linetension(const AgentContainer<Edge>& es) const {
        return get_energy_linetension(es, 0.);
    }
    
    /// Getter for energy associated with contractility of junctions
    /** Sums PCPVertex::edge_contractility_energy for provided entities
     */
    double get_energy_edge_contractility(const AgentContainer<Edge>& es) const {
        return get_energy_edge_contractility(es, 0.);
    }
    
    /// Getter for energy associated with area elasticity
    /** Sums PCPVertex::area_elasticity_energy for provided entities
     */
    double get_energy_areaelasticity (const AgentContainer<Cell>& cs) const {
        return get_energy_areaelasticity(cs, 0.);
    }

    /// Getter for energy associated with contractility of cells
    /** Sums PCPVertex::cell_contractility_energy for provided entities
     */
    double get_energy_cell_contractility (const AgentContainer<Cell>& cs) const
    {
        return get_energy_cell_contractility(cs, 0.);
    }

    /// Getter for energy associated with contractility of cells
    /** Sums the following energies for provided entities
     *      -# PCPVertex::get_energy_linetension
     *      -# PCPVertex::get_energy_edge_contractility
     *      -# PCPVertex::get_energy_areaelasticity
     *      -# PCPVertex::get_energy_cell_contractility
     *      -# PCPVertex::get_boundary_area_energy
     *      -# PCPVertex::get_boundary_shape_energy
     *      -# PCPVertex::get_boundary_stripe_energy
     */
    double get_energy(const AgentContainer<Edge>& es,
                      const AgentContainer<Cell>& cs) const {
        return get_energy(es, cs, 0.);
    }

    
    // .. Counter for transitions, etc.. ......................................

    /// Counter for the T1 neighborhood exchange transitions in last iteration
    std::size_t get_num_T1s() const {
        return _num_T1s;
    }

    /// Counter for the attempted T1 neighborhood exchange transitions
    /// in last iteration
    std::size_t get_num_T1s_attempted() const {
        return _num_T1s_attempted;
    }

    /// Counter for the T2 cell extrusion transitions in last iteration
    std::size_t get_num_T2s() const {
        return _num_T2s;
    }

    /// Total counter for the T1 neighborhood exchange transitions
    std::size_t get_num_T1s_total() const {
        return _num_T1s_total;
    }

    /// Total counter for the attempted T1 neighborhood exchange transitions
    std::size_t get_num_T1s_attempted_total() const {
        return _num_T1s_attempted_total;
    }

    /// Total counter for the T2 cell extrusion transitions
    std::size_t get_num_T2s_total() const {
        return _num_T2s_total;
    }

    /// Counter for the energy minimizations completed
    std::size_t get_num_minimizations() const {
        return _num_minimizations;
    }

    /// The entities manager - vertices, edges, cells
    const AgentManager& get_am () const {
        return _am;
    }


    // .. Model properties ....................................................
    
    /// Get linetension matrix
    CellCellPropertyMatrix get_linetension () const {
        return _linetension;
    }
    
    /// Set linetension matrix
    void set_linetension (CellCellPropertyMatrix linetension,
                          bool update_edges)
    {
        for (int i = 0; i < CellType::num_cell_types + 1; i++) {
            for (int j = i + 1; j < CellType::num_cell_types + 1; j++) {
                if (linetension(i, j) != linetension(j, i)) {
                    throw std::invalid_argument(fmt::format(
                            "Cannot set edge linetension! "
                            "Linetension matrix needs to be symmetric, "
                            "but entry ({}, {})={} and ({}, {})={}",
                            i, j, linetension(i, j),
                            j, i, linetension(j, i)));
                }
            }
        }
        _linetension = linetension;

        if (update_edges) {
            RuleFuncEdge update = [this] (const auto& edge)
            {
                auto state = edge->state;
                auto [a, b] = this->_am.adjoints_of(edge);
                if (a and b) {
                    state._linetension = this->_linetension(a->state.type,
                                                           b->state.type);
                }
                else {
                    if (not a) {
                        std::swap(a, b);
                    }
                    if (not a) {
                        throw std::runtime_error(fmt::format(
                            "Adjoints of edge {} are two nullptr!",
                            edge->id()));
                    }
                    state._linetension = this->_linetension(
                        a->state.type, CellType::num_cell_types);
                }
                return state;
            };

            apply_rule<Update::sync>(update, this->_am.edges());
        }
    }

    /// Getter for edge contractility matrix
    CellCellPropertyMatrix get_edge_contractility () const {
        return _edge_contractility;
    }
    
    /// Setter for edge contractility matrix
    void set_edge_contractility (CellCellPropertyMatrix contractility,
                                 bool update_edges) 
    {
        for (int i = 0; i < CellType::num_cell_types + 1; i++) {
            for (int j = i + 1; j < CellType::num_cell_types + 1; j++) {
                if (contractility(i, j) != contractility(j, i)) {
                    throw std::invalid_argument(fmt::format(
                            "Cannot set edge contractility! "
                            "Contractility matrix needs to be symmetric, "
                            "but entry ({}, {})={} and ({}, {})={}",
                            i, j, contractility(i, j),
                            j, i, contractility(j, i)));
                }
            }
        }
        _edge_contractility = contractility;

        if (update_edges) {
            RuleFuncEdge update = [this] (const auto& edge)
            {
                auto state = edge->state;
                auto [a, b] = this->_am.adjoints_of(edge);
                if (a and b) {
                    state._contractility = this->_edge_contractility(
                        a->state.type, b->state.type);
                }
                else {
                    if (not a) {
                        std::swap(a, b);
                    }
                    if (not a) {
                        throw std::runtime_error(fmt::format(
                            "Adjoints of edge {} are two nullptr!",
                            edge->id()));
                    }
                    state._contractility = this->_edge_contractility(
                        a->state.type, CellType::num_cell_types);
                }
                return state;
            };

            apply_rule<Update::sync>(update, this->_am.edges());
        }
    }

    /// Setter for boundary parameter
    void set_boundary_parameter (const Config& cfg) {
        _boundary_param = BoundaryParam(cfg);
    }

    /// Setter for fixed boundary
    void fix_boundary (bool fix_boundary = true) {
        _boundary_param.fix_boundary = fix_boundary;

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

    /// Enable or disable transitions
    void enable_transitions (bool enable_T1_transitions = true,
                             bool enable_T2_transitions = true)
    {
        _enable_T1_transitions = enable_T1_transitions;
        _enable_T2_transitions = enable_T2_transitions;
        _enable_transitions = (enable_T1_transitions or enable_T2_transitions);
    }

    void init_stripe_boundary(bool force_update = false) {
        if (_boundary_param.stripe_origin) {
            if (not force_update) {
                return;
            }

            this->_log->warn("Moving the origin of the boundary stripe!");
        }
        
        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::lowest();
        double y_min = std::numeric_limits<double>::max();
        double y_max = std::numeric_limits<double>::lowest();
        for (const auto &v : _am.vertices()) {
            SpaceVec pos = _am.position_of(v);
            x_min = std::min(x_min, pos[0]);
            x_max = std::max(x_max, pos[0]);
            y_min = std::min(y_min, pos[1]);
            y_max = std::max(y_max, pos[1]);
        }
        double L = x_max - x_min;
        double H = y_max - y_min;

        double rel = _boundary_param.stripe_curvature_center;

        _boundary_param.stripe_origin = std::make_shared<SpaceVec>(
            SpaceVec({rel * L + x_min, 0.5 * H + y_min}));
    }

    /// Set new parameter for a stripe boundary potential
    /** Set parameters for a straight stripe.
     * 
     *  Adds a quadratic potential on all cells that are outside a stripe of
     *  fixed width.
     */
    void set_stripe_boundary_width(double potential_constant,
                                   double stripe_width)
    {
        if (not _boundary_param.stripe_origin) {
            throw std::runtime_error("The origin of the boundary stripe has "
                "not been initialized!");
        }

        _boundary_param.stripe_potential_constant = potential_constant;
        _boundary_param.stripe_width = stripe_width;
    }

    /// Set new parameter for a stripe boundary potential
    /** Set parameters for a curved stripe of fixed width.
     * 
     *  Adds a quadratic potential on all cells outside of a stripe of a circle
     *  with mean radius = 1. / curvature of a fixed width.
     * 
     *  \note the width has to be set using
     *        set_stripe_boundary_parameters(potential_const, width)
     */
    void set_stripe_boundary_curvature(double potential_constant,
                                       double rel_curvature)
    {
        if (_boundary_param.stripe_width > 1.e12) {
            throw std::runtime_error(fmt::format(
                "To set stripe boundary parameters with rel. curvature {}, the "
                "stripe width cannot be infinite ({} > 1.e12). Set the "
                "stripe width before using a convergence and extension like "
                "process.", rel_curvature, _boundary_param.stripe_width));
            // NOTE use before:
            // set_stripe_boundary_parameters(potential_constant, stripe_width)
        }

        _boundary_param.stripe_potential_constant = potential_constant;

        // determine the length of the tissue
        double min_x = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        for (const auto& vertex : _am.vertices()) {
            SpaceVec pos = _am.position_of(vertex);
            min_x = std::min(min_x, pos[0]);
            max_x = std::max(max_x, pos[0]);
        }
        double origin_x = (*_boundary_param.stripe_origin)[0];
        double kappa_max = 1. / std::max(fabs(max_x - origin_x),
                                         fabs(min_x - origin_x));

        _boundary_param.stripe_curvature = kappa_max * rel_curvature;
    }


}; // class PCPVertex


} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPVERTEX_HH
