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

#include "space.hh"
#include "entities.hh"
#include "entities_manager.hh"
#include "initialisation.hh"
#include "transitions.hh"

namespace Utopia {
namespace Models {
namespace PCPVertex {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// The collection of parameters used for energy minimization
struct MinimizationParams {
    /// The tolerance in energy change 
    double tolerance;

    /// The maximum number of steps per minimization
    std::size_t max_steps;

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
     *           cell \f$ l = \sqrt{A_{domain} / #cells} \f$
     */
    double jiggle_intensity;

    template <typename Config>
    MinimizationParams(const Config& cfg)
    :
        tolerance(get_as<double>("tolerance", cfg)),
        max_steps(get_as<std::size_t>("max_steps", cfg)),
        num_repeat(get_as<std::size_t>("num_repeat", cfg)),
        jiggle_tolerance(get_as<double>("jiggle_tolerance", cfg, tolerance)),
        jiggle_intensity(get_as<double>("jiggle_intensity", cfg))
    {
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

    /// Initialize from config and inherit not defined values from default
    template <typename Config>
    MinimizationParams(const Config& cfg, const MinimizationParams& defaults)
    :
        tolerance(get_as<double>("tolerance", cfg, defaults.tolerance)),
        max_steps(get_as<std::size_t>("max_steps", cfg, defaults.max_steps)),
        num_repeat(get_as<std::size_t>("num_repeat", cfg, defaults.num_repeat)),
        jiggle_tolerance(get_as<double>("jiggle_tolerance", cfg,
                                        defaults.jiggle_tolerance)),
        jiggle_intensity(get_as<double>("jiggle_intensity", cfg,
                                        defaults.jiggle_intensity))
    {
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

    /// The type of a rule function acting on vertices of the agent manager
    using RuleFuncVertex = typename AgentManager::RuleFuncVertex;

    /// The type of a rule function acting on edges of the agent manager
    using RuleFuncEdge = typename AgentManager::RuleFuncEdge;

    /// The type of a rule function acting on cells of the agent manager
    using RuleFuncCell = typename AgentManager::RuleFuncCell;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The manager of the model's entities
    AgentManager _am;

    // PARAMETERS
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
    const enum UpdateScheme {
        SteepestGradient,
        SteepestGradientAdaptive,
        ConjugateGradient
    } _update_scheme;
    
    /// The tolerance during minimization
    double _minimization_tolerance;

    /// timestep scaling for polarity
    double _gamma;

    /// Linetension constant Lambda
    /** The entries are linetensions at interfaces between two cells of types i 
     *  and j.
     *  \note The values of this matrix are used for new edges, e.g. in T1
     *        transitions. There is no other generic way to determine the
     *        surface tension between two cells from the previous edge.
     *  \note This is a symmetric matrix
     */
    arma::Mat<double>::fixed<CellType::num_cell_types,
                             CellType::num_cell_types> _linetension;

    /// Linetension constant Lambda
    /** The entries are contractility at interfaces between two cells of types
     *  i and j
     *  \note The values of this matrix are used for new edges, e.g. in T1
     *        transitions. There is no other generic way to determine the
     *        surface tension between two cells from the previous edge.
     *  \note This is a symmetric matrix
     */
    arma::Mat<double>::fixed<CellType::num_cell_types,
                             CellType::num_cell_types> _edge_contractility;

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
    
    /// Area elasticity constant K
    double _area_elasticity;

    /// Cells with area smaller than this value are removed in T2 transition
    double _T2_threshold;

    /// Interaction parameter of cell-cell polarity interaction
    double _cell_cell_polarity_interaction;

    /// Interaction parameter of cell-internal polarity interaction
    double _cell_polarity_exclusion;
    
    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::uniform_real_distribution<double> _prob_distr;

    // .. Temporary objects ...................................................
    /// The total energy in the last step
    double _energy_previous_step;

    /// Current energy
    double _energy;

public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the PCPVertex model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel, typename... Taskargs>
    PCPVertex (const std::string name, ParentModel& parent, 
               Taskargs&&... taskargs)
    :
        // Initialize first via base model
        Base(name, parent, std::forward<Taskargs>(taskargs)...),

        _am(*this),
        
        // Get member paramters from cfg
        _dt(get_as<double>("dt", this->_cfg)),
        _update_scheme(this->setup_update_scheme(this->_cfg)),
        _minimization_tolerance(get_as<double>(
            "minimization_tolerance", this->_cfg)),
        _gamma(get_as<double>("gamma", this->_cfg)),
        _linetension(this->setup_linetension(this->_cfg)),
        _edge_contractility(this->setup_edge_contractility(this->_cfg)),
        _T1_threshold(get_as<double>("T1_threshold", this->_cfg)),
        _T1_separation(_T1_threshold *
                       get_as<double>("T1_separation_factor",this->_cfg)),
        _T1_probability(get_as<double>("T1_probability", this->_cfg)),
        _T1_barrier(get_as<double>("T1_barrier", this->_cfg)),
        _area_elasticity(get_as<double>("area_elasticity", this->_cfg)),
        _T2_threshold(get_as<double>("T2_threshold", this->_cfg)),
        _cell_cell_polarity_interaction(get_as<double>(
            "cell_cell_polarity_interaction",this->_cfg)),
        _cell_polarity_exclusion(get_as<double>(
            "cell_polarity_exclusion", this->_cfg)),
        _prob_distr(0.,1.),
        _energy_previous_step(0.),
        _energy(0.)
    {
        // this->initialise_polarity_random(get_as<double>(
        //         "cell_initialisation_protein_level", this->_cfg));

        double initial_jiggle(get_as<double>("initial_jiggle", this->_cfg, 0.));
        if (initial_jiggle > 0.) {
            jiggle_vertices(initial_jiggle);
        }
        
        this->_log->info("Model initialized.");
    }


private:
    // .. Setup functions .....................................................
    /// Setup function for the update scheme
    /** Currently implemented update schemes:
     *      -# steepest_gradient : Steepest gradient update at fixed step size
     *      -# steepest_gradient_adaptive : Steepest gradient update at adaptive
     *              step size. Step size is to next minimum of energy in the 
     *              direction of steepest gradient
     *      -# conjugate_gradient : Conjugate gradient update method
     */
    UpdateScheme setup_update_scheme(const Config& cfg) {
        const auto update_scheme = get_as<std::string>("update_scheme", cfg);

        if (update_scheme == "steepest_gradient") {
            this->_log->info("Steepest gradient chosen as update scheme.");
            return SteepestGradient;
        }
        if (update_scheme == "steepest_gradient_adaptive") {
            this->_log->info("Steepest gradient with adaptive step size "
                             "chosen as update scheme.");
            return SteepestGradientAdaptive;
        }
        if (update_scheme == "conjugate_gradient") {
            this->_log->info("Conjugate gradient chosen as update scheme.");
            return ConjugateGradient;
        }

        throw KeyError("update_scheme", cfg, 
            "Update scheme must be one of the following: "
                "'steepest_gradient', "
                "'steepest_gradient_adaptive', "
                "'conjugate_gradient'.");
    }

    /// Setup up the linetension from config
    arma::Mat<double>::fixed<CellType::num_cell_types,
                             CellType::num_cell_types> setup_linetension(
            const Config& cfg)
    {
        const auto edge_cfg = cfg["agent_manager"]["edge_manager"];
        const double linetension = get_as<double>("linetension",
                                                  edge_cfg["agent_params"]);
        // NOTE not checking paths, because happened in constructor of _am
        
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> matrix;
        return matrix.fill(linetension);
    }

    /// Setup up the edge contractility from config
    arma::Mat<double>::fixed<CellType::num_cell_types,
                             CellType::num_cell_types> setup_edge_contractility(
            const Config& cfg)
    {
        const auto edge_cfg = cfg["agent_manager"]["edge_manager"];
        const double contractility = get_as<double>("contractility",
                                                    edge_cfg["agent_params"]);
        // NOTE not checking paths, because happened in constructor of _am
        
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> matrix;
        return matrix.fill(contractility);
    }

    // /// Initialise the polarity proteins with random levels
    // /** This initialisation fulfills the polarity constrains of zero net
    //  *  polarisation and constant level of proteins per cell
    //  */
    // void initialise_polarity_random (double initialisation_protein_level);


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
    // double cell_cell_polarity_energy (Edge_ptr &e) const;
    // double polarity_exclusion_energy (Edge_ptr &a, Edge_ptr &b,
    //                                   const Cell_ptr &cell) const;
    // double cell_polarity_exclusion_energy (Cell_ptr &c) const;
    // double lagrange_net_polarisation_energy (Cell_ptr &c) const;
    // double lagrange_const_concentration_energy (Cell_ptr &c) const;

    // .. Force setter functions ..............................................
    /// Resets the forces of this vertex
    const RuleFuncVertex reset_forces = [] (const auto& vertex)
    {
        auto state = vertex->state;
        state.f.zeros();
        return state;
    };

    /// Reset the forces associated with polarity-proteins
    const RuleFuncEdge reset_polarity_change = [](const auto& edge) {
        auto state = edge->state;
        
        state.d_sigma_a = 0.;
        state.d_sigma_b = 0.;
        
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
        auto a = edge->custom_links().a;
        auto b = edge->custom_links().b;
        
        SpaceVec displ = this->_am.displacement(a, b);
        auto length = arma::norm(displ);

        SpaceVec force = edge->state.linetension * displ / length;

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
        auto a = edge->custom_links().a;
        auto b = edge->custom_links().b;

        SpaceVec force = edge->state.contractility *
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

        const auto cell_area = this->_am.area_of(cell);
        const auto cell_center = this->_am.barycenter_of(cell);
        
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

            double dA_dx = 0.5 * displ[1];
            double dA_dy = -0.5 * displ[0];

            SpaceVec force = -1. * this->_area_elasticity * 
                             (cell_area - state.area_preferential) * 
                             SpaceVec({dA_dx, dA_dy});

            v_center->state.f += force;
        }
        
        return state;
    };

    /// Contractility of the cell perimeter
    /** Associated energy per cell is 
     *  \Gamma / 2 L_cell^2, with L_cell the cell perimeter
     */
    const RuleFuncCell set_grad_cell_contractility = [this](const auto& cell)
    {
        auto state = cell->state;
        double perimeter = this->_am.perimeter_of(cell);

        for (auto [e, flip] : cell->custom_links().edges) {
            auto a = e->custom_links().a;
            auto b = e->custom_links().b;
            if (flip) { std::swap(a, b); }

            SpaceVec displ = this->_am.displacement(a, b);
            double length = arma::norm(displ);

            SpaceVec force = state.contractility * displ / length *
                             (perimeter - state.perimeter_preferential());

            a->state.f += force;
            b->state.f -= force;
        }

        return state;
    };


    // /// Set forces from cell-cell polarity interaction
    // /** 
    //  *  \return PCPVertex::cell_cell_polarity_energy
    //  */
    // std::function<void(Edge_ptr&)> set_cell_cell_polarity = [this](Edge_ptr &e) {
    //     e->d_sigma_a -= _cell_cell_polarity_interaction * e->sigma_b;
    //     e->d_sigma_b -= _cell_cell_polarity_interaction * e->sigma_a;
    // };

    // /// Set forces from cell intrinsic exclusion of polarity proteins
    // /** \return PCPVertex::polarity_exclusion_energy
    //  */
    // std::function<void(Edge_ptr&, Edge_ptr&,
    //                      const Cell_ptr&)> set_polarity_exclusion = 
    //         [this](Edge_ptr &a, Edge_ptr &b, const Cell_ptr &cell)
    // {
    //     double sigma_a = a->get_sigma(cell);
    //     double sigma_b = b->get_sigma(cell);

    //     double d_sigma_a = _cell_polarity_exclusion * sigma_b;
    //     double d_sigma_b = _cell_polarity_exclusion * sigma_a;

    //     a->set_d_sigma(cell, a->get_d_sigma(cell) + d_sigma_a);
    //     b->set_d_sigma(cell, b->get_d_sigma(cell) + d_sigma_b);
    // };

    // /// Applies cell intrinsic exclusion of polarity proteins for a cell
    // /** Pairwise calculation of PCPVertex::polarity_exclusion for all edges of c
    //  */
    // std::function<void(Cell_ptr&)> apply_polarity_exclusion = 
    //         [this](Cell_ptr &c)
    // {
    //     EdgeContainer edges;
    //     for (auto [e, flip] : c->edges_ordered) {
    //         edges.push_back(e);
    //     }
    //     std::function<void(int, int)> factory = [c, edges, this](
    //             int pos_a, int pos_b)
    //     {
    //         auto a = edges[pos_a];
    //         auto b = edges[pos_b];

    //         this->set_polarity_exclusion(a, b, c);
    //     };
    //     for (int i = 1; i < edges.size(); i++) {
    //         factory(i-1, i);
    //     }
    //     factory(edges.size()-1, 0);
    // };

    // /// Set forces from constraint of zero net polarisation
    // /** Within a cell the net polarisation is zero.
    //  *  Included via lagrange multiplier I
    //  * 
    //  *  \f$ E = -\lambda_I^\alpha \sum_i \sigma_i^\alpha \f$
    //  * 
    //  *  \warning There is no argument available why 
    //  *           \f$\dot{\lambda} = \gamma dE/d\lambda\f$
    //  *           instead of \f$\dot{\lambda} = -\gamma dE/d\lambda\f$
    //  */
    // std::function<void(Cell_ptr&)> set_lagrange_net_polarisation =
    //         [this](Cell_ptr &c)
    // {
    //     double cum_sigma = 0.;
    //     double lagrange = c->lagrange_net_polarisation;
    //     for (auto [e, flip] : c->edges_ordered) {
    //         double sigma = e->get_sigma(c);
    //         e->set_d_sigma(c, e->get_d_sigma(c) + lagrange);      

    //         cum_sigma += sigma;
    //     }
    //     // update the lagrangian
    //     c->lagrange_net_polarisation -= cum_sigma * this->_gamma;
    // };
    
    // /// Set forces from constraint of constant protein level
    // /** Within a cell the concentration of proteins is constant
    //  *  Included via lagrange multiplier II     * 
    //  * 
    //  *  \f$ E = -\lambda_{II}^\alpha (\sum_i (\sigma_i^\alpha)^2 - c^\alpha) \f$
    //  */
    // std::function<void(Cell_ptr&)> set_lagrange_const_concentration =
    //         [this](Cell_ptr &c)
    // {
    //     double concentration = 0.;
    //     double lagrange = c->lagrange_const_concentration;
    //     for (auto [e, flip] : c->edges_ordered) {
    //         double sigma = e->get_sigma(c);
    //         double d_sigma = 2 * sigma * lagrange;
    //         e->set_d_sigma(c, e->get_d_sigma(c) + d_sigma);         

    //         concentration += std::pow(sigma, 2);
    //     }
    //     // update the lagrangian
    //     c->lagrange_const_concentration -= (concentration - 
    //                                      c->protein_concentration) * this->_gamma;
    // };

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
        apply_rule<Update::async, Shuffle::off>(set_grad_cell_contractility,
                                                _am.cells());
    }

    /** The update of position
     * 
     *  Move vertex proportional to the gradient of energy (force)
     * 
     *  @param v    The pointer to the vertex to update
     */
    const RuleFuncVertex update_position = [this](const auto& vertex) {
        const auto state = vertex->state;
        // if (arma::norm(state.f) * this->_dt > 0.1) {
        //     throw std::runtime_error(fmt::format("Timestep badly chosen. " 
        //         "Vertex ({}, {}) would move by ({}, {}), i.e. by more than {}",
        //         _am.position_of(vertex)[0], _am.position_of(vertex)[1], 
        //         state.f[0], state.f[1],
        //         arma::norm(state.f)));
        // }
        _am.move_by(vertex, state.f * this->_dt);
        return state;
    };

    /** The update of polarity protein levels
     * 
     *  Change polarity level proportional to the gradient of energy (force)
     * 
     *  @param e    The pointer to the edge to update
     */
    // std::function<void(Edge_ptr&)> update_polarity = [this](Edge_ptr &e) {
    //     e->sigma_a += e->d_sigma_a * _gamma;
    //     e->sigma_b += e->d_sigma_b * _gamma;
    // };
    
    // -- The algorithm    ----------------------------------------------------
    // see algorithm.hh
    std::pair<double, double> determine_timestep (double dt,
                                                  const double energy_0,
                                                  const double tolerance) const;
    double steepest_gradient_step (bool adaptive_step);
    double conjugate_gradient_step ();
    double perform_update_step(UpdateScheme update_scheme);

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
    void divide_cell(std::shared_ptr<Cell> cell, double division_angle) {
        this->_log->info("Dividing cell..");
        return _am.divide_cell(cell, division_angle, _linetension,
                               _edge_contractility);
    }

    void increase_domain_size(double area);
    double stretch_domain(SpaceVec stretch, bool compensate,
                          bool fix_hc_volume);


    // .. Simulation Control ..................................................
    void init_minimization ();

    /// Iterate a single step
    /** \details Rules applied
     *      -# perform T2 transitions on cells
     *      -# perform T1 transitions on edges 
     *      -# perform minimization step
     */
    void perform_step () {
        bool transition_occurred = false;

        // T2 transitions -- cell extrusion
        for (int i = _am.cells().size() - 1; i >= 0; i--) {
            double area = _am.area_of(_am.cells()[i]);
            if (area < _T2_threshold) {
                this->_log->info("Removing cell in T2 transition in step {}..",
                                 this->_time);
                bool T2 = _am.remove_cell_T2(_am.cells()[i]);
                transition_occurred = transition_occurred or T2;
            }
        }

        // T1 transition -- neighborhood change
        for (int i = _am.edges().size() - 1; i >= 0; i--) {
            double length = _am.length_of(_am.edges()[i]);
            if (length < _T1_threshold
                and _prob_distr(*this->_rng) < _T1_probability)
            {
                this->_log->info("Removing edge in T1 transition in step {}..",
                                 this->_time);
                bool T1 = _am.remove_edge_T1(_am.edges()[i],
                        _linetension, _edge_contractility,
                        [this](const AgentContainer<Edge>& es,
                               const AgentContainer<Cell>& cs) { 
                                    return this->get_energy(es, cs, 0.); },
                        _T1_separation, _T1_barrier, _prob_distr(*this->_rng));
                transition_occurred = transition_occurred or T1;
            }
        }

        if (transition_occurred) {
            // restart the conjugate gradient update
            this->init_minimization();
        }

        _energy = this->get_energy();
        // NOTE need this, because operations might have changed it
        //      since last update
        _energy_previous_step = _energy;
        _energy = perform_update_step(_update_scheme);

        // if (_gamma > 0) {
        //     throw std::logic_error("Polarity proteins update not implemented!");
        //     apply_rule<Update::sync>(reset_polarity_change, _am.edges());

        //     std::for_each(_edges.begin(), _edges.end(),
        //                   set_cell_cell_polarity);
        //     std::for_each(_cells.begin(), _cells.end(),
        //         apply_polarity_exclusion);
        //     std::for_each(_cells.begin(), _cells.end(),
        //         set_lagrange_net_polarisation);
        //     std::for_each(_cells.begin(), _cells.end(),
        //         set_lagrange_const_concentration);

        //     std::for_each(_edges.begin(), _edges.end(), update_polarity);
        // }
    }
    
    /// Monitor model information
    void monitor () {
        this->_monitor.set_entry("energy", _energy);
        this->_monitor.set_entry("energy_change",
                                 _energy - _energy_previous_step);
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
     */
    std::size_t minimize_energy(const MinimizationParams& params)
    {
        const auto time_0 = this->get_time();

        this->_log->debug("Minimizing energy from step {}", time_0);

        for (std::size_t i = 0; i < params.num_repeat; i++)
        {
            this->jiggle_vertices(params.jiggle_intensity);

            double tolerance;
            if (i+1 == params.num_repeat) {
                tolerance = params.tolerance;
            }
            else {
                tolerance = params.jiggle_tolerance;
            }
            const auto time_start = this->get_time();
            bool minimum_reached = false;
            while (not minimum_reached) {
                _minimization_tolerance = tolerance;
                this->iterate();

                double energy_change = this->get_rel_energy_change();
                minimum_reached = (fabs(energy_change) < tolerance);
                
                if (not minimum_reached 
                    and this->get_time() - time_start >= params.max_steps)
                {
                    throw std::runtime_error(fmt::format(
                        "Equilibrium not reached within {} steps at a "
                        "tolerance of {}! Energy change in last step was {}.", 
                        params.max_steps, tolerance, energy_change));
                }

                if (stop_now.load()) {
                    this->_log->warn("Was told to stop. Not iterating "
                        "further ...");
                    throw GotSignal(received_signum.load());
                }
            }
        }

        this->_log->debug("Energy minimized within {} steps", 
                          this->get_time() - time_0);
        return this->get_time() - time_0;
    }

    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

    // NOTE when adding energy terms remember to add them to get_energy(..)
    double get_energy_linetension(const AgentContainer<Edge>& es,
                                  double beta = 0.) const;
    double get_energy_linetension() const {
        return get_energy_linetension(_am.edges(), 0.); }
    
    double get_energy_edge_contractility(const AgentContainer<Edge>& es,
                                         double beta = 0.) const;
    double get_energy_edge_contractility() const {
        return get_energy_edge_contractility(_am.edges(), 0.);
    }
    
    double get_energy_areaelasticity (const AgentContainer<Cell>& cs,
                                      double beta = 0.) const;
    double get_energy_areaelasticity () const {
        return get_energy_areaelasticity(_am.cells(), 0.);
    }
    
    double get_energy_cell_contractility (const AgentContainer<Cell>& cs,
                                          double beta = 0.) const;
    double get_energy_cell_contractility () const {
        return get_energy_cell_contractility(_am.cells(), 0.);
    }
    
    // double get_energy_cell_cell_polarity(const EdgeContainer& es) const;
    // double get_energy_cell_cell_polarity() const {
    //     return get_energy_cell_cell_polarity(_edges);
    // }
    
    // double get_energy_polarity_exclusion (const CellContainer& cs) const;
    // double get_energy_polarity_exclusion () const {
    //     return get_energy_polarity_exclusion(_cells);
    // }
    
    // double get_energy_lagrange_net_polarisation(const CellContainer& cs) const;
    // double get_energy_lagrange_net_polarisation() const {
    //     return get_energy_lagrange_net_polarisation(_cells);
    // }
    
    // double get_energy_lagrange_const_concentration(
    //         const CellContainer& cs) const;
    // double get_energy_lagrange_const_concentration() const {
    //     return get_energy_lagrange_const_concentration(_cells);
    // }
    
    double get_energy(const AgentContainer<Edge>& es,
                      const AgentContainer<Cell>& cs,
                      double beta = 0.) const;
    double get_energy(double beta = 0.) const {
        return get_energy(_am.edges(), _am.cells(), beta);
    }
    double get_rel_energy_change () const;

    const AgentManager& get_am () const {
        return _am;
    }
    
    const auto get_linetension () const {
        return _linetension;
    }
    
    /// Set 
    void set_linetension (
            arma::Mat<double>::fixed<CellType::num_cell_types,
                                     CellType::num_cell_types> linetension,
            bool update_edges)
    {
        for (int i = 0; i < CellType::num_cell_types; i++) {
            for (int j = i + 1; j < CellType::num_cell_types; j++) {
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
                const auto& [a, b] = this->_am.adjoints_of(edge);
                state.linetension = this->_linetension(a->state.type,
                                                         b->state.type);
                return state;
            };

            apply_rule<Update::sync>(update, this->_am.edges());
        }
    }

    const auto get_edge_contractility () const {
        return _edge_contractility;
    }

    void set_edge_contractility (
            arma::Mat<double>::fixed<CellType::num_cell_types,
                                     CellType::num_cell_types> contractility,
            bool update_edges) 
    {
        for (int i = 0; i < CellType::num_cell_types; i++) {
            for (int j = i + 1; j < CellType::num_cell_types; j++) {
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
                const auto& [a, b] = this->_am.adjoints_of(edge);
                state.contractility = this->_edge_contractility(a->state.type,
                                                                b->state.type);
                return state;
            };

            apply_rule<Update::sync>(update, this->_am.edges());
        }
    }
}; // class PCPVertex



// void PCPVertex::initialise_polarity_random (
//         double initialisation_protein_level)
// {
//     for (auto c : _cells) {
//         std::vector<double> rn(6);
//         double sum = 0.;
//         double sum_squares = 0.;
//         for (int i = 0; i < 5; i++) {
//             rn[i] = 2*_prob_distr(*this->_rng) - 1.;
//             sum += rn[i];
//             sum_squares += std::pow(rn[i], 2);
//         }
//         rn[5] = -sum;
//         sum_squares += std::pow(rn[5], 2);

//         std::shuffle(rn.begin(), rn.end(), *this->_rng);

//         int it = 0;
//         for (auto [e, flip] : c->edges_ordered) {
//             double norm = initialisation_protein_level / sqrt(sum_squares); 
//             e->set_sigma(c, rn[it++] * norm) ;
//         }
//     }
// }

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPVERTEX_HH
