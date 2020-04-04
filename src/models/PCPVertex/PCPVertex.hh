#ifndef UTOPIA_MODELS_PCPVERTEX_HH
#define UTOPIA_MODELS_PCPVERTEX_HH

// standard library includes
#include <random>
#include <math.h>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>
#include <utopia/core/agent_manager.hh>

#include "geometry.hh"
#include "entities.hh"
#include "entities_manager.hh"
#include "initialisation.hh"
#include "space.hh"

#ifndef PI
#define PI 3.14159265
#endif


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
template<bool periodic_bc>
class PCPVertex:
    public Model<PCPVertex<periodic_bc>, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<PCPVertex<periodic_bc>, ModelTypes>;

    /// Data type for the model time
    using Time = typename ModelTypes::Time;

    /// The manager of the entities (agents) of this model
    using AgentManager = CustomAgentManager<Model<PCPVertex<periodic_bc>, 
                                                  ModelTypes>>;

    /// The type of a Vertex
    using VertexNew = typename AgentManager::Vertex;

    /// The types of a cell
    using CellType = typename Cell::CellType;

    /// The type of a config
    using Config = Utopia::DataIO::Config;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The manager of the model's entities
    AgentManager _am;
    
    /// Container of vertices
    VertexContainer _vertices;
    /// Container of edges
    EdgeContainer _edges;
    /// Container of cells
    CellContainer _cells;

    // PARAMETERS
    /// timestep scaling
    double _dt;

    /// The update scheme to choose
    /** See PCPVertex::perform_update_step .
     */
    enum UpdateScheme {
        SteepestGradient,
        SteepestGradientAdaptive,
        ConjugateGradient
    } _update_scheme;
    
    /// The precision during minimisation
    double _minimisation_precision;

    /// timestep scaling for polarity
    double _gamma;

    /// The domain size in x
    double _Lx;
    /// The domain size in y
    double _Ly;

    /// Linetension constant Lambda
    /**  The entries are linetensions at interfaces between two cells of types i 
     *  and j
     * 
     *  \note This is a symmetric matrix
     */
    arma::Mat<double>::fixed<CellType::num_cell_types,
                             CellType::num_cell_types> _linetension;

    /// Linetension constant Lambda
    /**  The entries are linetensions at interfaces between two cells of types i 
     *  and j
     * 
     *  \note This is a symmetric matrix
     */
    arma::Mat<double>::fixed<CellType::num_cell_types,
                             CellType::num_cell_types> _edge_contractility;

    /// Edges shorter than this value are replaced in a T1 transition
    /** using absolute length
     */
    double _T1_threshold;

    /// The probability, that a T1 transition occurs
    double _T1_probability;

    /// The characteristic height of the energy barrier at T1 transitions
    /** The probability to perform a T1 transition is 
     *  \f$ p = \exp(-\Delta E / _T1_barrier) \f$, which is 1 for 
     *  \f$\Delta E < 0\f$.
     */
    double _T1_barrier;
    
    /// Area elasticity constant K
    double _area_elasticity;

    /// Cells with area smaller than this value are removed in T2 transition
    double _T2_threshold;

    /// Contractility of cell perimeter
    /** This reflects mechanics and contractility of the actin-myosin ring
     */
    double _contractility;

    /// Noise applied in order zero of forces
    double _noise_constant;

    /// Noise applied in first order of forces
    double _noise_linear;

    /// Interaction parameter of cell-cell polarity interaction
    double _cell_cell_polarity_interaction;

    /// Interaction parameter of cell-internal polarity interaction
    double _cell_polarity_exclusion;
    
    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::uniform_real_distribution<double> _prob_distr;
    
    /// A normal distribution used for evaluating noise of order zero
    std::normal_distribution<double> _distr_noise_const;

    /// A normal distribution used for evaluating noise of order one
    std::normal_distribution<double> _distr_noise_linear;

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
        
        _vertices(),
        _edges(),
        _cells(),
        
        // Get member paramters from cfg
        _dt(get_as<double>("dt", this->_cfg)),
        _minimisation_precision(1e-8),
        _gamma(get_as<double>("gamma", this->_cfg)),
        _Lx(100.), _Ly(100.),
        _linetension(),
        _edge_contractility(),
        _T1_threshold(get_as<double>("T1_threshold", this->_cfg)),
        _T1_probability(get_as<double>("T1_probability", this->_cfg)),
        _T1_barrier(get_as<double>("T1_barrier", this->_cfg)),
        _area_elasticity(get_as<double>("area_elasticity", this->_cfg)),
        _T2_threshold(get_as<double>("T2_threshold", this->_cfg)),
        _contractility(get_as<double>("contractility", this->_cfg)),
        _cell_cell_polarity_interaction(get_as<double>(
            "cell_cell_polarity_interaction",this->_cfg)),
        _cell_polarity_exclusion(get_as<double>("cell_polarity_exclusion", 
                                                this->_cfg)),
        _prob_distr(0.,1.),
        _distr_noise_const(0., get_as<double>("noise_constant", this->_cfg)),
        _distr_noise_linear(0., get_as<double>("noise_linear", this->_cfg)),
        _energy_previous_step(0.),
        _energy(0.)
    {
        _linetension.fill(get_as<double>("linetension", this->_cfg));
        _edge_contractility.fill(get_as<double>("edge_contractility",
                                                this->_cfg));

        if (get_as<std::string>("update_scheme",
                                this->_cfg) == "steepest_gradient") {
            _update_scheme = SteepestGradient;
            this->_log->info("Steepest gradient chosen as update scheme.");
        }
        else if (get_as<std::string>("update_scheme",
                                this->_cfg) == "steepest_gradient_adaptive") {
            _update_scheme = SteepestGradientAdaptive;
            this->_log->info("Steepest gradient with adaptive step size "
                             "chosen as update scheme.");
        }
        else if (get_as<std::string>("update_scheme",
                                this->_cfg) == "conjugate_gradient") {
            _update_scheme = ConjugateGradient;
            this->_log->info("Conjugate gradient chosen as update scheme.");
        }
        else {
            throw KeyError("update_scheme", this->_cfg, "Update scheme must "
                "be one of the following: "
                "'steepest_gradient', "
                "'steepest_gradient_adaptive', "
                "'conjugate_gradient'.");
        }

        this->initialise_polarity_random(get_as<double>(
                "cell_initialisation_protein_level", this->_cfg));

        jiggle_vertices(get_as<double>("initial_jiggle", this->_cfg, 0.));
        
        this->_log->info("Model initialized.");
    }


private:
    // .. Setup functions .....................................................

    /// Initialise the polarity proteins with random levels
    /** This initialisation fulfills the polarity constrains of zero net
     *  polarisation and constant level of proteins per cell
     */
    void initialise_polarity_random (double initialisation_protein_level);


    // The energy terms
    // See energy.hh for implementation

    double line_tension_energy (Edge e, double beta = 0.) const;
    double edge_contractility_energy (Edge e, double beta = 0.) const;
    double area_elasticity_energy (Cell c, double beta = 0.) const;
    double cell_contractility_energy (Cell c, double beta = 0.) const;
    double cell_cell_polarity_energy (Edge_ptr &e) const;
    double polarity_exclusion_energy (Edge_ptr &a, Edge_ptr &b,
                                      const Cell_ptr &cell) const;
    double cell_polarity_exclusion_energy (Cell_ptr &c) const;
    double lagrange_net_polarisation_energy (Cell_ptr &c) const;
    double lagrange_const_concentration_energy (Cell_ptr &c) const;

    // .. Force setter functions ..............................................
    /// Resets the forces of this vertex
    std::function<void(Vertex_ptr&)> reset_forces = [](Vertex_ptr &v) {
        v->fx = 0.;
        v->fy = 0.;
    };

    /// Reset the forces associated with polarity-proteins
    std::function<void(Edge_ptr&)> reset_polarity_change = [](Edge_ptr &e) {
        e->d_sigma_a = 0.;
        e->d_sigma_b = 0.;
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
    std::function<void(Edge_ptr&)> set_linetension = [this](Edge_ptr &e) {
        auto [dx, dy] = displacement_absolute<periodic_bc>(*e->a, *e->b, 
                                                           _Lx, _Ly);
        const double length = e->template length<periodic_bc>(_Lx, _Ly);
        const double fx = e->linetension*dx/length;
        const double fy = e->linetension*dy/length;

        e->a->fx += fx;
        e->a->fy += fy;
        e->b->fx -= fx;
        e->b->fy -= fy;
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
    std::function<void(Edge_ptr&)> set_edge_contractility = [this](Edge_ptr &e) {
        auto [dx, dy] = displacement_absolute<periodic_bc>(*e->a, *e->b, 
                                                           _Lx, _Ly);
        const double length = e->template length<periodic_bc>(_Lx, _Ly);
        const double fx = e->contractility*dx;
        const double fy = e->contractility*dy;

        e->a->fx += fx;
        e->a->fy += fy;
        e->b->fx -= fx;
        e->b->fy -= fy;
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
    std::function<void(Cell_ptr&)> set_area_elasticity = [this](Cell_ptr &c) {   
        double area_abs = c->template area_abs<periodic_bc>(_Lx, _Ly);     
        
        Site centre_site = *c->template centre_site<periodic_bc>();
        
        for (int edges_it = 0; edges_it < c->edges_ordered.size(); edges_it++) {
            Vertex_ptr v_center, v_prior, v_post;

            auto [e0, e0_flip]  = c->edges_ordered[std::max(0, edges_it - 1)];
            if (edges_it == 0) { 
                std::tie(e0, e0_flip) = c->edges_ordered.back();
            }

            const auto [e1, e1_flip] = c->edges_ordered[edges_it];
            
            // vertices in ordering
            v_center = e0->b, v_prior = e0->a;
            if (e0_flip) {
                std::swap(v_center, v_prior);
            }

            if (e1_flip) {
                v_post = e1->a;
            }
            else {
                v_post = e1->b;
            }

            auto center = periodic_copy<periodic_bc>(*v_center, centre_site);
            auto prior = periodic_copy<periodic_bc>(*v_prior, centre_site);
            auto post = periodic_copy<periodic_bc>(*v_post, centre_site);

            double dA_dx = 0.5 * (post.y - prior.y);
            double dA_dy = 0.5 * (prior.x - post.x);
            
            // this is the force on this vertex
            v_center->fx -= _area_elasticity * (
                    area_abs - c->area_preferential)*dA_dx;
            v_center->fy -= _area_elasticity * (
                    area_abs - c->area_preferential)*dA_dy;
        }
    };

    /// Contractility of the cell perimeter
    /** Associated energy per cell is 
     *  \Gamma / 2 L_cell^2, with L_cell the cell perimeter
     */
    std::function<void(Cell_ptr&)> set_cell_contractility = [this](Cell_ptr &c)
    {
        double perimeter = 0.;

        for (auto [e, flip] : c->edges_ordered) {
            perimeter += e->template length<periodic_bc>(_Lx, _Ly);;
        }

        for (auto [e, flip] : c->edges_ordered) {
            auto a = e->a; auto b = e->b;
            if (flip) {
                std::swap(a, b);
            }
            auto [dx, dy] = displacement_absolute<periodic_bc>(*a, *b, _Lx, _Ly);
            const double length = e->template length<periodic_bc>(_Lx, _Ly);
            double fx = c->contractility * perimeter * dx / length;
            double fy = c->contractility * perimeter * dy / length;

            a->fx += fx;
            a->fy += fy;
            b->fx += -fx;
            b->fy += -fy;
        }
    };


    /// Set forces from cell-cell polarity interaction
    /** 
     *  \return PCPVertex::cell_cell_polarity_energy
     */
    std::function<void(Edge_ptr&)> set_cell_cell_polarity = [this](Edge_ptr &e) {
        e->d_sigma_a -= _cell_cell_polarity_interaction * e->sigma_b;
        e->d_sigma_b -= _cell_cell_polarity_interaction * e->sigma_a;
    };

    /// Set forces from cell intrinsic exclusion of polarity proteins
    /** \return PCPVertex::polarity_exclusion_energy
     */
    std::function<void(Edge_ptr&, Edge_ptr&,
                         const Cell_ptr&)> set_polarity_exclusion = 
            [this](Edge_ptr &a, Edge_ptr &b, const Cell_ptr &cell)
    {
        double sigma_a = a->get_sigma(cell);
        double sigma_b = b->get_sigma(cell);

        double d_sigma_a = _cell_polarity_exclusion * sigma_b;
        double d_sigma_b = _cell_polarity_exclusion * sigma_a;

        a->set_d_sigma(cell, a->get_d_sigma(cell) + d_sigma_a);
        b->set_d_sigma(cell, b->get_d_sigma(cell) + d_sigma_b);
    };

    /// Applies cell intrinsic exclusion of polarity proteins for a cell
    /** Pairwise calculation of PCPVertex::polarity_exclusion for all edges of c
     */
    std::function<void(Cell_ptr&)> apply_polarity_exclusion = 
            [this](Cell_ptr &c)
    {
        EdgeContainer edges;
        for (auto [e, flip] : c->edges_ordered) {
            edges.push_back(e);
        }
        std::function<void(int, int)> factory = [c, edges, this](
                int pos_a, int pos_b)
        {
            auto a = edges[pos_a];
            auto b = edges[pos_b];

            this->set_polarity_exclusion(a, b, c);
        };
        for (int i = 1; i < edges.size(); i++) {
            factory(i-1, i);
        }
        factory(edges.size()-1, 0);
    };

    /// Set forces from constraint of zero net polarisation
    /** Within a cell the net polarisation is zero.
     *  Included via lagrange multiplier I
     * 
     *  \f$ E = -\lambda_I^\alpha \sum_i \sigma_i^\alpha \f$
     * 
     *  \warning There is no argument available why 
     *           \f$\dot{\lambda} = \gamma dE/d\lambda\f$
     *           instead of \f$\dot{\lambda} = -\gamma dE/d\lambda\f$
     */
    std::function<void(Cell_ptr&)> set_lagrange_net_polarisation =
            [this](Cell_ptr &c)
    {
        double cum_sigma = 0.;
        double lagrange = c->lagrange_net_polarisation;
        for (auto [e, flip] : c->edges_ordered) {
            double sigma = e->get_sigma(c);
            e->set_d_sigma(c, e->get_d_sigma(c) + lagrange);      

            cum_sigma += sigma;
        }
        // update the lagrangian
        c->lagrange_net_polarisation -= cum_sigma * this->_gamma;
    };
    
    /// Set forces from constraint of constant protein level
    /** Within a cell the concentration of proteins is constant
     *  Included via lagrange multiplier II     * 
     * 
     *  \f$ E = -\lambda_{II}^\alpha (\sum_i (\sigma_i^\alpha)^2 - c^\alpha) \f$
     */
    std::function<void(Cell_ptr&)> set_lagrange_const_concentration =
            [this](Cell_ptr &c)
    {
        double concentration = 0.;
        double lagrange = c->lagrange_const_concentration;
        for (auto [e, flip] : c->edges_ordered) {
            double sigma = e->get_sigma(c);
            double d_sigma = 2 * sigma * lagrange;
            e->set_d_sigma(c, e->get_d_sigma(c) + d_sigma);         

            concentration += std::pow(sigma, 2);
        }
        // update the lagrangian
        c->lagrange_const_concentration -= (concentration - 
                                         c->protein_concentration) * this->_gamma;
    };

    void set_gradient () {
        // reset forces
        std::for_each(_vertices.begin(), _vertices.end(), reset_forces);

        // apply new forces
        std::for_each(this->_edges.begin(), this->_edges.end(),
                        this->set_linetension);
        std::for_each(this->_edges.begin(), this->_edges.end(),
                        this->set_edge_contractility);
        std::for_each(this->_cells.begin(), this->_cells.end(),
                        this->set_area_elasticity);
        std::for_each(this->_cells.begin(), this->_cells.end(),
                        this->set_cell_contractility);

        // NOTE remember to add additional terms also to this->get_energy()
    }

    /** The update of position
     * 
     *  Move vertex proportional to the gradient of energy (force)
     * 
     *  @param v    The pointer to the vertex to update
     */
    std::function<void(Vertex_ptr&)> update_position = [this](Vertex_ptr &v) {
        double Df_lin = _distr_noise_linear(*this->_rng);
        double Df_const = _distr_noise_const(*this->_rng);
        v->x += ((1 + Df_lin) * v->fx + Df_const) * _dt / _Lx;

        Df_lin = _distr_noise_linear(*this->_rng);
        Df_const = _distr_noise_const(*this->_rng);
        v->y += ((1 + Df_lin) * v->fy + Df_const) * _dt / _Ly;
        
        correct_periodic_bc<periodic_bc>(v);
    };

    /** The update of polarity protein levels
     * 
     *  Change polarity level proportional to the gradient of energy (force)
     * 
     *  @param e    The pointer to the edge to update
     */
    std::function<void(Edge_ptr&)> update_polarity = [this](Edge_ptr &e) {
        e->sigma_a += e->d_sigma_a * _gamma;
        e->sigma_b += e->d_sigma_b * _gamma;
    };
    
    // .. Transitions ....................................................
    // See transitions.hh 

    std::pair<EdgeContainer::iterator,
              bool> T1_transition (EdgeContainer::iterator edge_it);
    
    std::pair<CellContainer::iterator,
              bool> T2_transition (CellContainer::iterator &cell_it);
    
    CellContainer::iterator divide_cell(CellContainer::iterator cell_it,
                                        double division_angle);

    // see algorithm.hh
    std::pair<double, double> determine_timestep (double dt,
                                                  const double energy_0) const;
    double steepest_gradient_step (bool adaptive_step);
    double conjugate_gradient_step ();
    double perform_update_step(UpdateScheme update_scheme);

    // -- Helper functions ----------------------------------------------------
    void differentiate_hair_cells_hlpr(
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> linetension,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> edge_contractility,
        arma::Col<double>::fixed<CellType::num_cell_types> area_preferential);


public:
    // -- Public Interface ----------------------------------------------------
    void jiggle_vertices(double intensity);
    void differentiate_hair_cells_random(double fraction,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                    CellType::num_cell_types> linetension,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                    CellType::num_cell_types> edge_contractility,
        arma::Col<double>::fixed<CellType::num_cell_types> area_preferential);
    template <class NotchDelta>
    void differentiate_hair_cells_NotchDelta(
        std::shared_ptr<NotchDelta> notch_delta, int steps,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> linetension,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> edge_contractility,
        arma::Col<double>::fixed<CellType::num_cell_types> area_preferential);
    void divide_cell(Cell_ptr cell, double division_angle);
    void increase_domain_size(double area);
    double stretch_domain(double dx, double dy, bool compensate,
                          bool fix_hc_volume);


    // .. Simulation Control ..................................................
    void init_minimisation ();

    /// Iterate a single step
    /** \details Rules applied
     *      -# perform T2 transitions on cells
     *      -# perform T1 transitions on edges 
     *      -# perform minimisation step
     */
    void perform_step () {
        bool transition_occurred = false;

        // T2 transitions -- cell extrusion
        for (auto c_it = _cells.begin(); c_it != _cells.end(); /*void*/) {
            double area = (*c_it)->template area_abs<periodic_bc>(_Lx, _Ly);
            if (area < _T2_threshold) {
                std::tie(c_it, transition_occurred) = T2_transition(c_it);
                transition_occurred = true;
            }
            else {
                ++c_it;
            }
        }

        // T1 transition -- neighborhood change
        for (auto e_it = _edges.begin(); e_it != _edges.end(); /*void*/) {
            double length = (*e_it)->template length<periodic_bc>(_Lx, _Ly);
            if (length < _T1_threshold
                and _prob_distr(*this->_rng) < _T1_probability)
            {
                std::tie(e_it, transition_occurred) = T1_transition(e_it);
                transition_occurred = true;
            }
            else {
                ++e_it;
            }
        }

        if (transition_occurred) {
            // restart the conjugate gradient update
            this->init_minimisation();
        }

        _energy = this->get_energy();
        // NOTE need this, because operations might have changed it
        //      since last update
        _energy_previous_step = _energy;
        _energy = perform_update_step(_update_scheme);

        if (_gamma > 0) {
            throw std::logic_error("Polarity proteins update not implemented!");
            std::for_each(_edges.begin(), _edges.end(),
                          set_cell_cell_polarity);
            std::for_each(_cells.begin(), _cells.end(),
                apply_polarity_exclusion);
            std::for_each(_cells.begin(), _cells.end(),
                set_lagrange_net_polarisation);
            std::for_each(_cells.begin(), _cells.end(),
                set_lagrange_const_concentration);

            std::for_each(_edges.begin(), _edges.end(), update_polarity);
        }
    }
    
    /// Monitor model information
    void monitor () {
        this->_monitor.set_entry("energy", _energy);
        this->_monitor.set_entry("energy_change",
                                 _energy - _energy_previous_step);
    }

    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

    // NOTE when adding energy terms remember to add them to get_energy(..)
    double get_energy_linetension(EdgeContainer es = {},
                                  double beta = 0.) const;
    double get_energy_linetension_normalised () const;
    double get_energy_edge_contractility(EdgeContainer es = {},
                                         double beta = 0.) const;
    double get_energy_areaelasticity (CellContainer cs = {},
                                      double beta = 0.) const;
    double get_energy_areaelasticity_normalised () const;
    double get_energy_cell_contractility (CellContainer cs = {},
                                          double beta = 0.) const;
    double get_energy_contractility () const;
    double get_energy_contractility_normalised () const;
    double get_energy_cell_cell_polarity(EdgeContainer es = {}) const;
    double get_energy_cell_cell_polarity_normalised() const;
    double get_energy_polarity_exclusion (CellContainer cs = {}) const;
    double get_energy_polarity_exclusion_normalised() const;
    double get_energy_lagrange_net_polarisation(CellContainer cs = {}) const;
    double get_energy_lagrange_net_polarisation_normalised() const;
    double get_energy_lagrange_const_concentration(CellContainer cs = {}) const;
    double get_energy_lagrange_const_concentration_normalised() const;
    
    double get_energy(EdgeContainer es = {}, CellContainer cs = {},
                      double beta = 0.) const;
    double get_energy_normalised () const;
    double get_rel_energy_change () const;

    /// Getter for the domain size
    const std::pair<double, double> get_domain_size () const {
        return std::make_pair(_Lx, _Ly);
    }

    /// Getter for the periodic bc
    const bool get_periodic_bc () {
        return periodic_bc;
    }

    AgentManager& get_am () {
        return _am;
    }

    /// Change the const noise distribution
    void set_noise_const (double stddev) {
        _distr_noise_const = std::normal_distribution<double>(0., stddev);
    }

    /// Change the linear noise distribution
    void set_noise_linear (double stddev) {
        _distr_noise_linear = std::normal_distribution<double>(0., stddev);
    }

    /** Criterion for the equilibrium state
     * 
     *  Equilibrium if PCPVertex::get_mean_energy_change() change is smaller
     *  than threshold value
     *  
     *  \param threshold    The equilibrium threshold
     * 
     *  \return {whether equilibrium reached, relative energy change}
     *  \note relative energy change is not normalised for step size
     */
    std::pair<bool, double> equilibrium_state_reached() const {
        double energy_change = this->get_rel_energy_change ();
        return std::make_pair(fabs(energy_change) < _minimisation_precision,
                              energy_change);
    };

    /// Set a precision for minimisation
    void set_minimisation_precision (double precision) {
        _minimisation_precision = precision;
    }
}; // class PCPVertex



template <bool periodic_bc>
void PCPVertex<periodic_bc>::initialise_polarity_random (
        double initialisation_protein_level)
{
    for (auto c : _cells) {
        std::vector<double> rn(6);
        double sum = 0.;
        double sum_squares = 0.;
        for (int i = 0; i < 5; i++) {
            rn[i] = 2*_prob_distr(*this->_rng) - 1.;
            sum += rn[i];
            sum_squares += std::pow(rn[i], 2);
        }
        rn[5] = -sum;
        sum_squares += std::pow(rn[5], 2);

        std::shuffle(rn.begin(), rn.end(), *this->_rng);

        int it = 0;
        for (auto [e, flip] : c->edges_ordered) {
            double norm = initialisation_protein_level / sqrt(sum_squares); 
            e->set_sigma(c, rn[it++] * norm) ;
        }
    }
}

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPVERTEX_HH
