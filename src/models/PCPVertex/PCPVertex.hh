#ifndef UTOPIA_MODELS_PCPVERTEX_HH
#define UTOPIA_MODELS_PCPVERTEX_HH

// standard library includes
#include <random>
#include <math.h>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>

#include "geometry.hh"

#ifndef PI
#define PI 3.14159265
#endif


namespace Utopia {
namespace Models {
namespace PCPVertex {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<DefaultRNG, WriteMode::managed>;


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
template<bool periodic_bc, bool polarity_proteins>
class PCPVertex:
    public Model<PCPVertex<periodic_bc, polarity_proteins>, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<PCPVertex<periodic_bc, polarity_proteins>, ModelTypes>;

    /// Data type for the model time
    using Time = typename ModelTypes::Time;

    /// The types of a cell
    using CellType = typename Cell::CellType;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// Container of vertices
    VertexContainer _vertices;
    /// Container of edges
    EdgeContainer _edges;
    /// Container of cells
    CellContainer _cells;

    // PARAMETERS
    /// timestep scaling
    double _dt;
    
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

    /// The change of energy in last update
    double _energy_change;

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
        _energy_change(0.)
    {
        _linetension.fill(get_as<double>("linetension", this->_cfg));
        _edge_contractility.fill(get_as<double>("edge_contractility",
                                                this->_cfg));

        this->initialise_hexagonal(
            get_as<double>("hexagon_size", this->_cfg),
            get_as<int>("lattice_rows", this->_cfg),
            get_as<int>("lattice_columns", this->_cfg)
        );

        this->initialise_polarity_random(get_as<double>(
                "cell_initialisation_protein_level", this->_cfg));

        jiggle_vertices(get_as<double>("initial_jiggle", this->_cfg, 0.));
        
        this->_log->info("Model initialized.");
    }


private:
    // .. Setup functions .....................................................
    /** Initialiser for randomly distributed cells.
     * Use a Voronoi decomposition from randomly distributed cell centres
     * 
     */
    void initialise_voronoi(int num_cells);

    /// Initialiser for a hexagonal arrangement of cells.
    void initialise_hexagonal(double size, int num_rows, int num_columns);

    /// Initialise the polarity proteins with random levels
    /** This initialisation fulfills the polarity constrains of zero net
     *  polarisation and constant level of proteins per cell
     */
    void initialise_polarity_random (double initialisation_protein_level);

    // .. Helper functions ....................................................
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

    /// The energy associated with linetension per edge
    /** \f$ E = \sum_{ij} \lambda_{ij} l_ij \f$
     */
    double line_tension_energy (Edge e, double beta = 0.) const {
        if (beta > 0) {
            e = displace_edge_steepest_gradient<periodic_bc>(e, beta, _Lx, _Ly);
        }
        const double length = e.template length<periodic_bc>(_Lx, _Ly);
        return e.linetension * length;
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

    double edge_contractility_energy (Edge e, double beta = 0.) const {
        if (beta > 0) {
            e = displace_edge_steepest_gradient<periodic_bc>(e, beta, _Lx, _Ly);
        }
        const double length = e.template length<periodic_bc>(_Lx, _Ly);
        return 0.5 * e.contractility * pow(length, 2);
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

    /// The energy associated with area elasticity
    /** \f$ E = K/2 * (A - A0)**2 \f$
     */
    double area_elasticity_energy (Cell c, double beta = 0.) const {
        if (beta > 0.) {
            c = displace_cell_steepest_gradient<periodic_bc>(c, beta, _Lx, _Ly);
        }
        double area_abs = c.template area_abs<periodic_bc>(_Lx, _Ly);
        return 0.5 * _area_elasticity * pow(area_abs - c.area_preferential, 2);
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

    /// The energy associated with cell contractility
    double cell_contractility_energy (Cell c, double beta = 0.) const {
        if (beta > 0.) {
            c = displace_cell_steepest_gradient<periodic_bc>(c, beta, _Lx, _Ly);
        }

        double perimeter = 0.;
        for (auto [e, flip] : c.edges_ordered) {
            perimeter += e->template length<periodic_bc>(_Lx, _Ly);;
        }
        return 0.5 *c.contractility * std::pow(perimeter, 2);
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

    /// The energy associated with cell-cell polarity
    /** \f$ E = J_1 \sum_i \sigma_i^\alpha \sigma_i^\beta \f$, where
     *  \f$ i \f$ is naming an edge, \f$ \alpha \f$ and \f$ \beta \f$ naming the
     *  two neighbouring cells to edge \f$ i \f$. \f$ J_1 \f$ is an interaction
     *  parameter PCPVertex::_cell_cell_polarity_interaction.
     */
    double cell_cell_polarity_energy (Edge_ptr &e) const {
        return _cell_cell_polarity_interaction * e->sigma_a * e->sigma_b;
    };


    /// Set forces from cell-cell polarity interaction
    /** 
     *  \return PCPVertex::cell_cell_polarity_energy
     */
    std::function<void(Edge_ptr&)> set_cell_cell_polarity = [this](Edge_ptr &e) {
        e->d_sigma_a -= _cell_cell_polarity_interaction * e->sigma_b;
        e->d_sigma_b -= _cell_cell_polarity_interaction * e->sigma_a;
    };

    /// The energy associated with cell intrinsic exclusion of polarity proteins
    /** Interaction of proteins on neighbouring edges.
     *  For `cell_polarity_exclusion` > 0 accumulation of opposite sign polarity
     *  proteins is disfavoured on neighbouring edges
     *
     *  \f$ E = - J_2 \sum_{<i, j>} \sigma_i^\alpha \sigma_j^\alpha \f$, where
     *  edges \f$ i \f$ and \f$ j \f$ are adjacent bonds of cell \f$ \alpha \f$.
     *  \f$ J_2 \f$ is interaction parameter PCPVertex::_cell_polarity_exclusion.
     * 
     *  \param a    edge a, alias \f$ i \f$
     *  \param b    edge b, that is to be an interaction partner of a,
     *              for instance a neighbour, alias edge \f$ j \f$
     *  \param cell The cell \f$ \alpha \f$ within which exclusion is applied
     */
     double polarity_exclusion_energy (Edge_ptr &a, Edge_ptr &b,
                                       const Cell_ptr &cell) const
    {        
        double sigma_a = a->get_sigma(cell);
        double sigma_b = b->get_sigma(cell);
        
        return - _cell_polarity_exclusion * sigma_a * sigma_b;
    };

    /// The energy of polarity exclusion for the entire cell
    /** Cumulative for the pairwise interactions of edges within cell, 
     *  see PCPVertex::polarity_exclusion energy
     */
    double cell_polarity_exclusion_energy (Cell_ptr &c) const {
        EdgeContainer edges;
        double energy = 0.;
        for (auto [e, flip] : c->edges_ordered) {
            edges.push_back(e);
        }
        std::function<double(int, int)> factory = [c, edges, this](
                int pos_a, int pos_b)
        {
            auto a = edges[pos_a];
            auto b = edges[pos_b];

            return this->polarity_exclusion_energy(a, b, c);
        };
        for (int i = 1; i < edges.size(); i++) {
            energy += factory(i-1, i);
        }
        energy += factory(edges.size()-1, 0);
        
        return energy;
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
    
    /// The energy associated with constraint of zero net polarisation
    /** Within a cell the net polarisation is zero.
     *  Included via lagrange multiplier I
     * 
     *  \f$ E = -\lambda_I^\alpha \sum_i \sigma_i^\alpha \f$
     * 
     *  \warning There is no argument available why 
     *           \f$\dot{\lambda} = \gamma dE/d\lambda\f$
     *           instead of \f$\dot{\lambda} = -\gamma dE/d\lambda\f$
     */
    double lagrange_net_polarisation_energy (Cell_ptr &c) const {
        double cum_sigma = 0.;
        for (auto [e, flip] : c->edges_ordered) {
            cum_sigma += e->get_sigma(c);
        }

        return - c->lagrange_net_polarisation * cum_sigma;
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
    
    /// Energy associated with constraint of constant protein level
    /** Within a cell the concentration of proteins is constant
     *  Included via lagrange multiplier II     * 
     * 
     *  \f$ E = -\lambda_{II}^\alpha (\sum_i (\sigma_i^\alpha)^2 - c^\alpha) \f$
     */
    double lagrange_const_concentration_energy (Cell_ptr &c) const {
        double concentration = 0.;
        for (auto [e, flip] : c->edges_ordered) {
            concentration += std::pow(e->get_sigma(c), 2);
        }

        double lagrange = c->lagrange_const_concentration;
        return - lagrange * (concentration - c->protein_concentration);
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

    /** Perform a T1 transition on the edge edge_it in _edges
     * 
     *  edge_it will be removed and a new edge orthogonal to edge will be created
     */
    std::pair<EdgeContainer::iterator,
              bool> T1_transition (EdgeContainer::iterator edge_it);
    
    std::pair<CellContainer::iterator,
              bool> T2_transition (CellContainer::iterator &cell_it);
    

    /// Perform a cell division on specific cell
    /** Divides a specific cell into two identical cells with properties derived
     *  from the common parent cell. 
     *  The division is performed at a given angle through the parent cell's
     *  center. This defines the axis of division that will form a new edge 
     *  between the two new cells.
     * 
     *  The new edge has properties as given for initialisation.
     * 
     *  \param cell_it      iterator to the cell within _cells that is to be 
     *                      divided
     *  \param division_angle   angle (in rad) at which the cell
     * 
     *  \return iterator to the element following cell_it
     */
    CellContainer::iterator divide_cell(CellContainer::iterator cell_it,
                                        double division_angle);

    /// The step size to reach the line minimum along steepest gradient
    /** See www.acclab.helsinki.fi/~aakurone/atomistiset/lecturenotes/lecture12_2up.pdf
     */
    std::pair<double, double> determine_timestep (const double energy_0) const
    {
        // Bracket the minimum
        double dt = std::min(std::max(_dt, 1e-5), 1e-2);

        double pos_1 = 0.; // left boundary
        bool calculate_energy_min = true;
        double energy_1 = energy_0;
        bool calculate_energy_2 = true;
        double pos_2 = dt; // right boundary
        double energy_2 = this->get_energy(_edges, _cells, pos_2);
        double pos_min = dt/2.;
        double energy_min = this->get_energy(_edges, _cells, pos_min);
        while (true) {
            if (dt > 3.) {
                if (fabs(energy_2 - energy_1) < _minimisation_precision) {
                    // the energy function is flat
                    return std::make_pair(0., energy_0);
                }
                this->_log->warn("At step of {} along direction of "
                    "update the energy is still decreasing by {}", dt/2, 
                    energy_2 - energy_1);
                for (double dt = 1e-11; dt < 1000.; dt *= 2) {
                    this->_log->warn("At step of {} along direction of "
                        "update the energy is changing by {}", dt, 
                        this->get_energy(_edges, _cells, dt) - energy_1);
                }
                // NOTE only if energy_2 < energy_1: dt -> 2*dt
                throw std::runtime_error("Unable to find bracket to "
                    "local energy minimum!");
            }
            if (dt < 1e-10) {
                // the energy is increasing, hence already at local minimum
                // NOTE only if energy_min > energy_1: dt -> dt / 2
                return std::make_pair(0., energy_0);
            }

            if (energy_2 < energy_1) {
                dt *= 2.;
                energy_min = energy_2;
                energy_2 = this->get_energy(_edges, _cells, dt);
                continue; // dt was not a right boundary to minimum
            }

            if (energy_min > energy_1) {
                dt = dt/2;
                energy_2 = energy_min;
                energy_min = this->get_energy(_edges, _cells, dt/2);
                continue; // we know that close to pos_1 there is a minimum
            }

            // Both conditions fulfilled!
            // there is a minimum within interval [0, pos_2]
            pos_2 = dt;
            pos_min = dt / 2.;
            break;
        }

        if (fabs(energy_min - energy_0) < _minimisation_precision) {
            // the found minimum fulfills our condition 
            return std::make_pair(pos_min, energy_min);
        }

        // Find the minimum between brackets with parabolic interpolation
        while(true) {
            double term_3 = (pos_2 - pos_1)*(energy_2 - energy_min);
            double term_4 = (pos_2 - pos_min)*(energy_2 - energy_1);
            double term_1 = (pos_2 - pos_1)*term_3;
            double term_2 = (pos_2 - pos_min)*term_4;
            
            // the minimum of parabola through 1, 2, min
            double pos_4 = pos_2 - 0.5*(term_1 - term_2)/(term_3 - term_4);
            if (pos_4 < pos_1 or pos_4 > pos_2) {
                this->_log->warn("pos_1={}, pos_min={}, pos_2={}",
                                 pos_1, pos_min, pos_2);
                this->_log->warn("D_energy_1={}, D_energy_min={}, "
                                 "D_energy_2={}. wrt energy_0",
                                 energy_0-energy_1, energy_0-energy_min,
                                 energy_0-energy_2);
                this->_log->warn("Fitted minimum: x={}", pos_4);
                throw std::runtime_error("Energy minimisation failed! "
                                         "Parabola fit outside brackets");
            }

            double energy_4 = this->get_energy(_edges, _cells, pos_4);
            if (fabs(energy_4 - energy_min) < _minimisation_precision) {
                dt = pos_4;
                break; // found the minimum with required precision
            }
            
            if (energy_4 - energy_0 > 0.) {
                throw std::runtime_error("Energy minimisation failed! "
                                         "Found closer local maximum, hence "
                                         "overestimated size of brackets.");
            }

            // ** choose new brackets
            // pos_4 is to the left of minimum and is smaller
            // pos_4 becomes new minimum between 1 and 2->min
            if (pos_4 - pos_min < 0. and energy_4 - energy_min < 0.) {
                pos_2 = pos_min; energy_2 = energy_min;
                pos_min = pos_4; energy_min = energy_4;
            }
            // pos_4 is to the left of minimum but is larger
            // pos_4 becomes new pos 1
            else if (pos_4 - pos_min < 0.) {
                pos_1 = pos_4; energy_1 = energy_4;
            }
            // pos_4 is to the right of minimum and is smaller
            // pos_4 becomes new minimum between 1->min and 2
            else if (energy_4 - energy_min < 0.) {
                pos_1 = pos_min; energy_1 = energy_min;
                pos_min = pos_4; energy_min = energy_4;
            }
            // pos_4 is to the right of minimum but is larger
            // pos_4 becomes new pos 2
            else {
                pos_2 = pos_4; energy_2 = energy_4;
            }
        }

        return std::make_pair(dt, energy_min);
    }

public:
    // -- Public Interface ----------------------------------------------------
    
    /// Apply a perturbation to the position of vertices
    /** Move the x and y position by a random value in [-intensity, intensity]
     *  using a uniform distribution.
     * 
     * TODO write test
     */
    void jiggle_vertices(double intensity) {
        this->_log->debug("Jiggling the vertices on a length scale of "
                          "{} ..", intensity);
        for (auto v : _vertices) {
            v->x += 2*intensity * _prob_distr(*this->_rng) - intensity;
            v->y += 2*intensity * _prob_distr(*this->_rng) - intensity;
            correct_periodic_bc<periodic_bc>(v);
        }
    }

    /// Differentiates progenitor cells with random hair cell distribution
    /** \param fraction     fraction of hair cells. Others are support cells
     *  \param linetension  The symmetric matrix of linetension interactions
     *                      between two cells of same or different type
     *  \param area_preferential    The preferential cell area of the different
     *                              cell types
     */
    void differentiate_hair_cells(double fraction,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> linetension,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> edge_contractility,
        arma::Col<double>::fixed<CellType::num_cell_types> area_preferential)
    {
        this->_log->debug("Differentiating progenitor cells to {}% hair cells "
            "and {}% support cells ...", fraction, 1-fraction);

        for (int i = 0; i < CellType::num_cell_types; i++) {
            for (int j = i+1; j < CellType::num_cell_types; j++) {
                if (linetension(j, i) == linetension(i, j)) {
                    continue;
                }
                
                this->_log->warn("Invalid argument in differentiate_hair_cells."
                    "Got non symmetric 'linetension' matrix!");
                throw std::invalid_argument("Non symmetric 'linetension'"
                    "matrix");
            }
        }
        for (int i = 0; i < CellType::num_cell_types; i++) {
            for (int j = i+1; j < CellType::num_cell_types; j++) {
                if (edge_contractility(j, i) == edge_contractility(i, j)) {
                    continue;
                }
                
                this->_log->warn("Invalid argument in differentiate_hair_cells."
                    "Got non symmetric 'edge_contractility' matrix!");
                throw std::invalid_argument("Non symmetric 'edge_contractility'"
                    "matrix");
            }
        }
        
        // Set the cell type
        for (auto &c : _cells) {
            if (_prob_distr(*this->_rng) < fraction) { 
                c->type = CellType::hair; }
            else { c->type = CellType::support; }

            // set the preferential area
            c->area_preferential = area_preferential(c->type);
        }

        // Set the surface tension
        for (auto &e : _edges) {
            Cell::CellType cell_a_type, cell_b_type;
            if (e->adj_cell_a.expired()) {
                cell_a_type = Cell::CellType::support; }
            else { cell_a_type = e->adj_cell_a.lock()->type; }
            if (e->adj_cell_b.expired()) {
                cell_b_type = Cell::CellType::support; }
            else { cell_b_type = e->adj_cell_b.lock()->type; }

            e->linetension = linetension(cell_a_type, cell_b_type);
            e->contractility = edge_contractility(cell_a_type, cell_b_type);
        }

        _linetension = linetension;
        _edge_contractility = edge_contractility;
    }

    /// Perform a cell division on specific cell
    /** Divides a specific cell into two identical cells with properties derived
     *  from the common parent cell. 
     *  The division is performed at a given angle through the parent cell's
     *  center. This defines the axis of division that will form a new edge 
     *  between the two new cells.
     * 
     *  \param cell     pointer to the cell that is to be divided
     *  \param division_angle   angle (in rad) at which the cell
     */
    void divide_cell(Cell_ptr cell, double division_angle)
    {
        auto cell_it = std::find(_cells.begin(), _cells.end(), cell);

        if (cell_it == _cells.end()) {
            throw std::invalid_argument("Cannot divide cell at position "
                "({}, {}), because its not a member of cells in vertex model!");
        }

        this->divide_cell(cell_it, division_angle);
    }

    /// Increase the domain size by a certain area
    /** This remaps the domain of size A to size A + dA while keeping the 
     *  relation of Lx to Ly constant.
     * 
     *  Thereby proliferation of cells can be performed in a periodic setup
     *  without changing the parameters of the system. 
     */
    void increase_domain_size(double area) {
        if (-1. * area > _Lx * _Ly) {
            throw std::invalid_argument("Cannot decrease the domain size by "
                    "an area larger than the domain size. dA = " + 
                    std::to_string(area) + " and A = " +
                    std::to_string(_Lx * _Ly));
        }
        double ratio = _Lx / double(_Ly);
        _Ly = std::sqrt(_Ly*_Ly + area / ratio);
        _Lx = ratio * _Ly;
    }

    /// Stretch the domain size
    /** \param dx   The stretching distance in x
     *  \param dy   The stretching distance in y
     *  \param compensate   Whether to compensate the growth of the dissue by 
     *                      increase of preferential area
     *  \param fix_hc_volume   Whether to fix the volume of type CellType::hair
     * 
     *  \return The total change in area
     */
    double stretch_domain(double dx, double dy, bool compensate,
                          bool fix_hc_volume)
    {
        this->_log->debug("stretching domain by ({}, {}). Compensate {}, "
                          "fix hair cell volume {}", dx, dy, compensate, 
                          fix_hc_volume);
        _Lx += dx;
        _Ly += dy;

        if (not compensate) {
            return dx * _Ly + dy * _Lx;
        }

        int num_cells = _cells.size();
        if (fix_hc_volume) {
            for (auto c : _cells) {
                num_cells -= c->type == CellType::hair;
            }
            if (num_cells == 0) {
                throw std::runtime_error("All support cells eliminated!");
            }
        }

        double dA = (dx * _Ly + dy * _Lx) / num_cells;

        std::function<void(Cell_ptr&)> compensate_dA = [dA](Cell_ptr& cell) {
            cell->area_preferential += dA;
            return;
        };
        std::function<void(Cell_ptr&)> compensate_dA_non_hc = [dA](
                Cell_ptr& cell)
        {
            if (cell->type != CellType::hair) {
                cell->area_preferential += dA;
            }
            return;
        };

        if (not fix_hc_volume) {
            std::for_each(_cells.begin(), _cells.end(), compensate_dA);
        }
        else {
            std::for_each(_cells.begin(), _cells.end(), compensate_dA_non_hc);
        }

        return dx * _Ly + dy * _Lx;
    }

    // .. Simulation Control ..................................................

    void init_minimisation () {
        set_gradient();

        for (auto v : _vertices) {
            v->gx = v->x; v->gy = v->y;
            v->hx = v->x; v->hy = v->y;
        }
    }
    
    /// Iterate a single step
    /** \details Rules applied
     *      -# reset vertex forces, calculate cell area and edge length
     *      -# perform T2 transitions on cells
     *      -# perform T1 transitions on edges 
     *      -# Linetension on edges
     *      -# Area elasticity on cells
     *      -# Update vertex positions on vertices
     */
    void perform_step () {
        _energy_previous_step = this->get_energy();

        // line minimisation along direction of update h
        double new_energy;
        std::tie(_dt, new_energy) = determine_timestep(_energy_previous_step);
        _energy_change = (new_energy - _energy_previous_step) / new_energy;
        
        if (_energy_change < -1e-14) {
            this->_log->debug("Updating with timestep {} at energy change {}",
                              _dt, _energy_change);
            std::for_each(_vertices.begin(), _vertices.end(),
                          update_position);
        }
        else {
            this->_log->debug("NOT updating with step size {} along direction "
                              "of update at energy change {}", _dt,
                              _energy_change);
            return;
        }        

        // determine the conjugate gradient direction
        set_gradient();
        double gamma = 0.;
        double g_square = 0.;
        for (auto v : _vertices) {
            gamma += std::pow(v->fx, 2) + std::pow(v->fy, 2);
            g_square += std::pow(v->gx, 2) + std::pow(v->gy, 2);
            
            v->gx = v->fx; v->gy = v->fy;
        }
        gamma /= g_square;
        for (auto v : _vertices) {
            v->hx = v->gx + gamma * v->hx;
            v->hy = v->gy + gamma * v->hy;

            v->fx = v->hx; v->fy = v->hy; 
        }


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

        // if constexpr (polarity_proteins) {
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
        double energy = this->get_energy();
        this->_monitor.set_entry("energy", energy);
        this->_monitor.set_entry("energy_change",
                                 energy - _energy_previous_step);
    }


    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model
    /// Getter for energy associated with linetension
    double get_energy_linetension(EdgeContainer es = {}, double beta = 0.) const
    {
        if (es.empty()) { es = this->_edges; }
        double energy = 0.;
        for (auto &&e : es) {
            energy += line_tension_energy(*e, beta);
        }
        return energy;
    }

    /// Getter for normalised energy associated with linetension
    /** The energy is normalised to the number of edges
     */
    double get_energy_linetension_normalised () const {
        return get_energy_linetension() / double(_edges.size());
    }
    
    /// Getter for energy associated with contractility of junctions
    double get_energy_edge_contractility(EdgeContainer es = {},
                                         double beta = 0.) const
    {
        if (es.empty()) { es = this->_edges; }
        double energy = 0.;
        for (auto &&e : es) {
            energy += edge_contractility_energy(*e, beta);
        }
        return energy;
    }

    /// Getter for energy associated with area elasticity
    double get_energy_areaelasticity (CellContainer cs = {},
                                      double beta = 0.) const
    {
        if (cs.empty()) { cs = this->_cells; }
        double energy = 0.;
        for (auto &&c : cs) {
            energy += area_elasticity_energy(*c, beta);
        }
        return energy;
    }

    /// Getter for normalised energy associated with area elasticity
    /** The energy is normalised to the number of cells
     */
    double get_energy_areaelasticity_normalised () const {
        return get_energy_areaelasticity() / double(_cells.size());
    }

    /// Getter for energy associated with contractility of cells
    double get_energy_cell_contractility (CellContainer cs = {},
                                          double beta = 0.) const
    {
        if (cs.empty()) { cs = this->_cells; }
        double energy = 0.;
        for (auto &&c : cs) {
            energy += cell_contractility_energy(*c, beta);
        }
        return energy;
    }

    /// Getter for normalised energy associated with contractility of cells
    /** The energy is normalised to the number of cells
     */
    double get_energy_contractility () const {
        return get_energy_cell_contractility() + 
            get_energy_edge_contractility();
    }

    /// Getter for normalised energy associated with contractility of cells
    /** The energy is normalised to the number of cells
     */
    double get_energy_contractility_normalised () const {
        return get_energy_cell_contractility() / _cells.size() + 
            get_energy_edge_contractility() / _edges.size();
    }


    /// Getter for energy associated with cell-cell polarity
    double get_energy_cell_cell_polarity(
            EdgeContainer es = {}) const
    {
        if (es.empty()) { es = this->_edges; }
        double energy = 0.;
        for (auto &&e : es) {
            energy += cell_cell_polarity_energy(e);
        }
        return energy;
    }


    /// Getter for energy normalised associated with cell-cell polarity
    /** The energy is normalised to the number of cells
     */
    double get_energy_cell_cell_polarity_normalised() const {
        return get_energy_cell_cell_polarity() / double(_cells.size());
    }

    /// Getter for energy associated with polarity exclusion
    double get_energy_polarity_exclusion (
            CellContainer cs = {}) const
    {
        if (cs.empty()) { cs = this->_cells; }
        double energy = 0.;
        for (auto &&c : cs) {
            energy += cell_polarity_exclusion_energy(c);
        }
        return energy;
    }

    /// Getter for energy associated with polarity exclusion
    /** The energy is normalised to the number of cells
     */
    double get_energy_polarity_exclusion_normalised() const {
        return get_energy_polarity_exclusion() / double(_cells.size());
    }

    /// Getter for energy associated with lagrange multiplier I
    /** Langrange multiplier I is the constrain of zero net polarisation within
     *  a cell.
     */
    double get_energy_lagrange_net_polarisation(
            CellContainer cs = {}) const
    {
        if (cs.empty()) { cs = this->_cells; }
        double energy = 0.;
        for (auto &&c : cs) {
            energy += lagrange_net_polarisation_energy(c);
        }
        return energy;
    }

    /// Getter for normalised energy associated with lagrange multiplier I
    /** Langrange multiplier I is the constrain of zero net polarisation within
     *  a cell.
     * 
     *  The energy is normalised to the number of cells
     */
    double get_energy_lagrange_net_polarisation_normalised() const {
        return get_energy_lagrange_net_polarisation() / double(_cells.size());
    }

    /// Getter for energy associated with lagrange multiplier II
    /** Langrange multiplier II is the constrain of constant level of proteins
     */
    double get_energy_lagrange_const_concentration(
            CellContainer cs = {}) const
    {
        if (cs.empty()) { cs = this->_cells; }
        double energy = 0.;
        for (auto &&c : cs) {
            energy += lagrange_const_concentration_energy(c);
        }
        return energy;
    }

    /// Getter for normalised energy associated with lagrange multiplier II
    /** Langrange multiplier II is the constrain of constant level of proteins
     * 
     *  The energy is normalised to the number of cells
     */
    double get_energy_lagrange_const_concentration_normalised() const {
        return get_energy_lagrange_const_concentration() / _cells.size();
    }

    /// Getter for energy
    double get_energy (EdgeContainer es = {}, CellContainer cs = {}, 
                       double beta = 0.) const
    {
        if (beta > 10.) {
            throw std::runtime_error("Cannot calculate energy with "
                "epsilon multiplicator larger than 10.!");
        }

        if (es.empty()) { es = this->_edges; }
        if (cs.empty()) { cs = this->_cells; }
        return get_energy_linetension(es, beta) +
            get_energy_edge_contractility(es, beta) +
            get_energy_areaelasticity(cs, beta) +
            get_energy_cell_contractility(cs, beta) +
            get_energy_cell_cell_polarity(es) +
            get_energy_polarity_exclusion(cs) +
            get_energy_lagrange_net_polarisation(cs) +
            get_energy_lagrange_const_concentration(cs);
    }

    /// Getter for normalised energy
    /** The energy is normalised wrt number of vertices, edges, or cells, 
     *  respectively.
     */
    double get_energy_normalised () const {
        return get_energy_linetension_normalised() +
            get_energy_areaelasticity_normalised() +
            get_energy_contractility_normalised() +
            get_energy_cell_cell_polarity_normalised() +
            get_energy_polarity_exclusion_normalised() +
            get_energy_lagrange_net_polarisation_normalised() +
            get_energy_lagrange_const_concentration_normalised();
    }

    /// Getter for the relative energy change from previous to last step
    double get_rel_energy_change () const {
        const double energy = get_energy();
        double energy_change = energy - _energy_previous_step;
        return energy_change / energy;
    }

    /// Getter for the domain size
    const std::pair<double, double> get_domain_size () const {
        return std::make_pair(_Lx, _Ly);
    }

    /// Getter for the periodic bc
    const bool get_periodic_bc () {
        return periodic_bc;
    }

    /// Getter for vertices
    std::vector<std::weak_ptr<Vertex>> get_vertices () {
        std::vector<std::weak_ptr<Vertex>> vs;
        for (auto &v : _vertices) {
            vs.push_back(v);
        }
        return vs;
    }

    /// Getter for edges
    std::vector<std::weak_ptr<Edge>> get_edges () {
        std::vector<std::weak_ptr<Edge>> es;
        for (auto &e : _edges) {
            es.push_back(e);
        }
        return es;
    }

    /// Getter for cells
    std::vector<std::weak_ptr<Cell>> get_cells () {
        std::vector<std::weak_ptr<Cell>> cs;
        for (auto &c : _cells) {
            cs.push_back(c);
        }
        return cs;
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
     */
    bool equilibrium_state_reached() const {
        return fabs(_energy_change) < _minimisation_precision;
    };

    /// Set a precision for minimisation
    void set_minimisation_precision (double precision) {
        _minimisation_precision = precision;
    }
}; // class PCPVertex

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPVERTEX_HH
