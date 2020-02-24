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

    /// Edges shorter than this value are replaced in a T1 transition
    /** using relative length
     */
    double _length_threshold; 
    
    /// Area elasticity constant K
    double _area_elasticity;

    /// The preferential area of a cell
    arma::Col<double>::fixed<CellType::num_cell_types> _area_preferential;

    /// Cells with area smaller than this value are removed in T2 transition
    double _area_threshold;

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
    /// The current energy of area elasticity
    double _energy_linetension;

    /// The current energy of area elasticity
    double _energy_areaelasticity;

    /// The current energy of contractility
    double _energy_contractility;

    /// The current energy of cell-cell polarity
    double _energy_cell_cell_polarity;

    /// The current energy of polarity exclusion
    double _energy_polarity_exclusion;

    /// The current energy of associated with constraint of zero net polarisation
    double _energy_lagrange_net_polarisation;

    /// The current energy of associated with constraint of const protein level
    double _energy_lagrange_const_concentration;

    /// The total energy in the last step
    double _energy_previous_step;
    
    /// The time over which to average energe change
    /** in PCPVertex::get_mean_energy_change
     */
    int _energy_change_history_length;

    /// Relative energy change history
    /** of length PCPVertex::_energy_change_history_length
     */
    std::list<double> _energy_change_history;

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
        _gamma(get_as<double>("gamma", this->_cfg)),
        _Lx(100.), _Ly(100.),
        _linetension(),
        _length_threshold(get_as<double>("length_threshold", this->_cfg)),
        _area_elasticity(get_as<double>("area_elasticity", this->_cfg)),
        _area_preferential(),
        _area_threshold(get_as<double>("area_threshold", this->_cfg)),
        _contractility(get_as<double>("contractility", this->_cfg)),
        _cell_cell_polarity_interaction(get_as<double>(
            "cell_cell_polarity_interaction",this->_cfg)),
        _cell_polarity_exclusion(get_as<double>("cell_polarity_exclusion", 
                                                this->_cfg)),
        _prob_distr(0.,1.),
        _distr_noise_const(0., get_as<double>("noise_constant", this->_cfg)),
        _distr_noise_linear(0., get_as<double>("noise_linear", this->_cfg)),
        _energy_linetension(0.),
        _energy_areaelasticity(0.),
        _energy_contractility(0.),
        _energy_cell_cell_polarity(0.),
        _energy_polarity_exclusion(0.),
        _energy_lagrange_net_polarisation(0.),
        _energy_lagrange_const_concentration(0.),
        _energy_previous_step(0.),
        _energy_change_history_length(get_as<int>(
            "energy_change_history_length", this->_cfg)),
        _energy_change_history(_energy_change_history_length, 0.)
    {
        static_assert(CellType::num_cell_types == 3, "Initialisation of "
            "interaction matrices `linetension` and `area_preferential` only "
            "defined for 3 cell types.");
        
        this->_log->debug("Extracting area-preferential (expecting {} entries) ..",
                          CellType::num_cell_types);
        if (not this->_cfg["area_preferential"]) {
            throw std::invalid_argument("Expected cfg dict 'area_preferential' "
                "not available!");
        }
        _area_preferential(CellType::progenitor) = get_as<double>(
            "progenitor", this->_cfg["area_preferential"]);
        _area_preferential(CellType::hair) = get_as<double>(
            "hair", this->_cfg["area_preferential"]);
        _area_preferential(CellType::support) = get_as<double>(
            "support", this->_cfg["area_preferential"]);

        this->_log->debug("Extracting linetension (expecting {}! entries, "
                          "i.e. the upper diagonal matrix of a {}x{} matrix) ..",
                          CellType::num_cell_types, CellType::num_cell_types,
                          CellType::num_cell_types);
        if (not this->_cfg["linetension"]) {
            throw std::invalid_argument("Expected cfg dict 'linetension' "
                "not available!");
        }
        _linetension(CellType::progenitor, CellType::progenitor) = get_as<double>(
            "progenitor_progenitor", this->_cfg["linetension"]);
        _linetension(CellType::progenitor, CellType::hair) = get_as<double>(
            "progenitor_hair", this->_cfg["linetension"]);
        _linetension(CellType::progenitor, CellType::support) = get_as<double>(
            "progenitor_support", this->_cfg["linetension"]);
        _linetension(CellType::hair, CellType::hair) = get_as<double>(
            "hair_hair", this->_cfg["linetension"]);
        _linetension(CellType::hair, CellType::support) = get_as<double>(
            "hair_support", this->_cfg["linetension"]);
        _linetension(CellType::support, CellType::support) = get_as<double>(
            "support_support", this->_cfg["linetension"]);
        for (int i = 0; i < CellType::num_cell_types; i++) {
            for (int j = i+1; j < CellType::num_cell_types; j++) {
                _linetension(j, i) = _linetension(i, j);
            }
        }

        this->initialise_hexagonal(
            get_as<double>("hexagon_size", this->_cfg),
            get_as<int>("lattice_rows", this->_cfg),
            get_as<int>("lattice_columns", this->_cfg)
        );

        this->initialise_polarity_random(get_as<double>(
                "cell_initialisation_protein_level", this->_cfg));

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
    std::function<double(Edge_ptr&)> line_tension = [this](Edge_ptr &e) {
        auto [dx, dy] = displacement_absolute<periodic_bc>(*e->a, *e->b, 
                                                           _Lx, _Ly);
        const double length = e->template length<periodic_bc>(_Lx, _Ly);
        const double fx = e->linetension*dx/length;
        const double fy = e->linetension*dy/length;

        e->a->fx += fx;
        e->a->fy += fy;
        e->b->fx -= fx;
        e->b->fy -= fy;

        return e->linetension * length;
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
    std::function<double(Cell_ptr&)> area_elasticity = [this](Cell_ptr &c) {        
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

            auto center = periodic_copy<periodic_bc>(*v_center, *c->s);
            auto prior = periodic_copy<periodic_bc>(*v_prior, *c->s);
            auto post = periodic_copy<periodic_bc>(*v_post, *c->s);

            double dA_dx = 0.5 * (post.y - prior.y) * c->area_sgn;
            double dA_dy = 0.5 * (prior.x - post.x) * c->area_sgn;
            
            // this is the force on this vertex
            v_center->fx -= _area_elasticity * (
                    c->area_abs(_Lx, _Ly) - c->area_preferential)*dA_dx;
            v_center->fy -= _area_elasticity * (
                    c->area_abs(_Lx, _Ly) - c->area_preferential)*dA_dy;
        }

        return 0.5 * _area_elasticity * pow(c->area_abs(_Lx, _Ly) - c->area_preferential, 2);
    };

    /// Contractility of the cell perimeter
    /** Associated energy per cell is 
     *  \Gamma / 2 L_cell^2, with L_cell the cell perimeter
     */
    std::function<double(Cell_ptr&)> contractility = [this](Cell_ptr &c) {
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

        return 0.5 *c->contractility * std::pow(perimeter, 2);
    };

    /// Calculates forces of cell-cell polarity interaction
    /** Associated energy is
     *  
     *  \f$ E = J_1 \sum_i \sigma_i^\alpha \sigma_i^\beta \f$, where
     *  \f$ i \f$ is naming an edge, \f$ \alpha \f$ and \f$ \beta \f$ naming the
     *  two neighbouring cells to edge \f$ i \f$. \f$ J_1 \f$ is an interaction
     *  parameter PCPVertex::_cell_cell_polarity_interaction.
     *  
     *  \return energy associated with this edge
     */
    std::function<double(Edge_ptr&)> cell_cell_polarity = [this](Edge_ptr &e) {
        e->d_sigma_a -= _cell_cell_polarity_interaction * e->sigma_b;
        e->d_sigma_b -= _cell_cell_polarity_interaction * e->sigma_a;

        return _cell_cell_polarity_interaction * e->sigma_a * e->sigma_b;
    };

    /// Calculates forces on cell intrinsic exclusion of polarity proteins
    /** Interaction of proteins on neighbouring edges.
     *  For `cell_polarity_exclusion` > 0 accumulation of opposite sign polarity
     *  proteins is disfavoured on neighbouring edges
     * 
     *  Associated energy is
     *  \f$ E = - J_2 \sum_{<i, j>} \sigma_i^\alpha \sigma_j^\alpha \f$, where
     *  edges \f$ i \f$ and \f$ j \f$ are adjacent bonds of cell \f$ \alpha \f$.
     *  \f$ J_2 \f$ is interaction parameter PCPVertex::_cell_polarity_exclusion.
     * 
     *  \param a    edge a, alias \f$ i \f$
     *  \param b    edge b, that is to be an interaction partner of a,
     *              for instance a neighbour, alias edge \f$ j \f$
     *  \param cell The cell \f$ \alpha \f$ within which exclusion is applied
     */
    std::function<double(Edge_ptr&, Edge_ptr&,
                         const Cell_ptr&)> polarity_exclusion = 
            [this](Edge_ptr &a, Edge_ptr &b, const Cell_ptr &cell)
    {
        double sigma_a = a->get_sigma(cell);
        double sigma_b = b->get_sigma(cell);

        double d_sigma_a = _cell_polarity_exclusion * sigma_b;
        double d_sigma_b = _cell_polarity_exclusion * sigma_a;

        a->set_d_sigma(cell, a->get_d_sigma(cell) + d_sigma_a);
        b->set_d_sigma(cell, b->get_d_sigma(cell) + d_sigma_b);
        
        return - _cell_polarity_exclusion * sigma_a * sigma_b;
    };

    /// Applies cell intrinsic exclusion of polarity proteins for a cell
    /** Pairwise calculation of PCPVertex::polarity_exclusion for all edges of c
     */
    std::function<double(Cell_ptr&)> apply_polarity_exclusion = 
            [this](Cell_ptr &c)
    {
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
            bool flip_a = (a->adj_cell_b.lock() == c);
            bool flip_b = (b->adj_cell_b.lock() == c);

            return this->polarity_exclusion(a, b, c);
        };
        for (int i = 1; i < edges.size(); i++) {
            energy += factory(i-1, i);
        }
        energy += factory(edges.size()-1, 0);
        return energy;
    };

    /// Constraint of zero net polarisation
    /** Within a cell the net polarisation is zero.
     *  Included via lagrange multiplier I
     * 
     *  \f$ E = -\lambda_I^\alpha \sum_i \sigma_i^\alpha \f$
     * 
     *  \warning There is no argument available why 
     *           \f$\dot{\lambda} = \gamma dE/d\lambda\f$
     *           instead of \f$\dot{\lambda} = -\gamma dE/d\lambda\f$
     */
    std::function<double(Cell_ptr&)> lagrange_net_polarisation =
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

        return - lagrange * cum_sigma;
        
        // double cum_sigma = 0.;
        // double lagrange = c->lagrange_net_polarisation;
        // for (auto [e, flip] : c->edges_ordered) {
        //     double sigma = e->get_sigma(c);   

        //     cum_sigma += sigma;
        // }
        // for (auto [e, flip] : c->edges_ordered) {
        //     double d_sigma = -2. * lagrange * cum_sigma;
        //     e->set_d_sigma(c, e->get_d_sigma(c) + d_sigma);   
        // }

        // // update the lagrangian
        // c->lagrange_net_polarisation -= std::pow(cum_sigma, 2) * this->_dt;

        // return lagrange * std::pow(cum_sigma, 2);
    };
    
    /// Constraint of constant protein level
    /** Within a cell the concentration of proteins is constant
     *  Included via lagrange multiplier II     * 
     * 
     *  \f$ E = -\lambda_{II}^\alpha (\sum_i (\sigma_i^\alpha)^2 - c^\alpha) \f$
     */
    std::function<double(Cell_ptr&)> lagrange_const_concentration =
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

        return - lagrange * (concentration - c->protein_concentration);

        // double concentration = 0.;
        // double lagrange = c->lagrange_const_concentration;
        // for (auto [e, flip] : c->edges_ordered) {
        //     double sigma = e->get_sigma(c);       

        //     concentration += std::pow(sigma, 2);
        // }
        // for (auto [e, flip] : c->edges_ordered) {
        //     double sigma = e->get_sigma(c);
        //     double d_sigma = -4. * (concentration - c->protein_concentration) * sigma * lagrange;

        //     e->set_d_sigma(c, e->get_d_sigma(c) + d_sigma);
        // }
        // // update the lagrangian
        // c->lagrange_const_concentration -= std::pow(concentration - 
        //                                  c->protein_concentration, 2) * this->_dt;

        // return lagrange * std::pow(concentration - c->protein_concentration, 2);
    };

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

    /** Perform a T1 transition on the edge edge_it in _edges
     * 
     *  edge_it will be removed and a new edge orthogonal to edge will be created
     */
    EdgeContainer::iterator T1_transition (EdgeContainer::iterator edge_it);
    
    /** Erase cell from _cells while updating the topology (T2 transition)
     * 
     *  Removes the cell, its edges and its vertices and sets up a vertex at its
     *  centre.
     * 
     *  returns _cells.erase(cell_it)
     */
    CellContainer::iterator T2_transition (CellContainer::iterator &cell_it);
    

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

public:
    // -- Public Interface ----------------------------------------------------
    
    /// Differentiates progenitor cells with random hair cell distribution
    /** \param fraction     fraction of hair cells. Others are support cells
     */
    void differentiate_hair_cells(double fraction)
    {
        this->_log->debug("Differentiating progenitor cells to {}% hair cells "
            "and {}% support cells ...", fraction, 1-fraction);
        
        // Set the cell type
        for (auto &c : _cells) {
            if (_prob_distr(*this->_rng) < fraction) { 
                c->type = CellType::hair; }
            else { c->type = CellType::support; }

            c->area_preferential = _area_preferential(c->type);
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

            e->linetension = _linetension(cell_a_type,
                                          cell_b_type);
        }
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
     */
    void stretch_domain(double dx, double dy, bool compensate) {        
        _Lx += dx;
        _Ly += dy;

        if (not compensate) {
            return;
        }

        double dA = dx * _Ly + dy * _Lx;
        dA /= _cells.size();

        std::function<void(Cell_ptr&)> compensate_dA = [dA](Cell_ptr& cell) {
            cell->area_preferential += dA;
            return;
        };

        std::for_each(_cells.begin(), _cells.end(), compensate_dA);

        return;
    }

    // .. Simulation Control ..................................................

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
        // store previous energy
        _energy_change_history.pop_front();
        _energy_change_history.push_back(get_rel_energy_change ());
        _energy_previous_step = this->get_energy();

        // reset forces
        std::for_each(_vertices.begin(), _vertices.end(), reset_forces);
        if constexpr (polarity_proteins) {
            std::for_each(_edges.begin(), _edges.end(), reset_polarity_change);
        }

        // calculate cell area
        for (auto &c : _cells) {
            c->template cell_area<periodic_bc>();
        }

        // T2 transitions -- cell extrusion
        for (auto c_it = _cells.begin(); c_it != _cells.end(); /*void*/) {
            if ((*c_it)->area_abs(_Lx, _Ly) > _area_threshold) {
                ++c_it;
            }
            else {
                // erase c_it and update topology
                c_it = T2_transition(c_it);
            }
        }

        // T1 transition -- neighborhood change
        for (auto e_it = _edges.begin(); e_it != _edges.end(); /*void*/) {
            if ((*e_it)->template length<periodic_bc>(_Lx, _Ly) < _length_threshold)
            {
                e_it = T1_transition(e_it);
            }
            else {
                ++e_it;
            }
        }

        // line tension
        _energy_linetension = 0.;
        for (auto &e : _edges) {
            _energy_linetension += line_tension(e);
        }
        
        // area elasticity
        _energy_areaelasticity = 0.;
        for (auto &c : _cells) {
            _energy_areaelasticity += area_elasticity(c);
        }
        
        // contractility of perimeter
        _energy_contractility = 0.;
        for (auto &c : _cells) {
            _energy_contractility += contractility(c);
        }
        
        // update vertex positions from forces
        std::for_each(_vertices.begin(), _vertices.end(), update_position);

        if constexpr (polarity_proteins) {
            _energy_cell_cell_polarity = 0.;
            for (auto &e : _edges) {
                _energy_cell_cell_polarity += cell_cell_polarity(e);
            }

            _energy_polarity_exclusion = 0.;
            for (auto &c : _cells) {
                _energy_polarity_exclusion += apply_polarity_exclusion(c);
            }

            _energy_lagrange_net_polarisation = 0.;
            for (auto &c : _cells) {
                _energy_lagrange_net_polarisation += lagrange_net_polarisation(c);
            }

            _energy_lagrange_const_concentration = 0.;
            for (auto &c : _cells) {
                _energy_lagrange_const_concentration += lagrange_const_concentration(c);
            }

            for (auto& e : _edges) {
                update_polarity(e);
            }
        }
    }
    
    /// Monitor model information
    void monitor () { }

    /// The prolog
    /** Performs the following tasks:
     *      1. calculate the energies
     *      2. default prolog tasks
     */
    void prolog ();


    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

    /// Getter for energy associated with linetension
    /** The energy is normalised to the number of edges
     */
    double get_energy_linetension_normalised () const {
        return _energy_linetension / double(_edges.size());
    }

    /// Getter for energy associated with area elasticity
    /** The energy is normalised to the number of cells
     */
    double get_energy_areaelasticity_normalised () const {
        return _energy_areaelasticity / double(_cells.size());
    }

    /// Getter for energy associated with contractility
    /** The energy is normalised to the number of cells
     */
    double get_energy_contractility_normalised () const {
        return _energy_contractility / double(_cells.size());
    }


    /// Getter for energy associated with cell-cell polarity
    /** The energy is normalised to the number of cells
     */
    double get_energy_cell_cell_polarity_normalised() const {
        return _energy_cell_cell_polarity / double(_cells.size());
    }

    /// Getter for energy associated with polarity exclusion
    /** The energy is normalised to the number of cells
     */
    double get_energy_polarity_exclusion_normalised() const {
        return _energy_polarity_exclusion / double(_cells.size());
    }

    /// Getter for energy associated with lagrange multiplier I
    /** Langrange multiplier I is the constrain of zero net polarisation within
     *  a cell.
     * 
     *  The energy is normalised to the number of cells
     */
    double get_energy_lagrange_net_polarisation_normalised() const {
        return _energy_lagrange_net_polarisation / double(_cells.size());
    }

    /// Getter for energy associated with lagrange multiplier II
    /** Langrange multiplier II is the constrain of constant level of proteins
     * 
     *  The energy is normalised to the number of cells
     */
    double get_energy_lagrange_const_concentration_normalised() const {
        return _energy_lagrange_const_concentration / double(_cells.size());
    }

    /// Getter for energy
    double get_energy () const {
        return _energy_linetension + _energy_areaelasticity + 
            _energy_contractility + _energy_cell_cell_polarity +
            _energy_polarity_exclusion + _energy_lagrange_net_polarisation +
            _energy_lagrange_const_concentration;
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
        double energy_change = get_energy() - _energy_previous_step;
        return energy_change / get_energy();
    }

    /// Getter for mean energy change
    /** Mean is taken over last m steps, with 
     *  m = PCPVertex::_energy_change_history_length
     */
    double get_mean_energy_change() const {
        return std::accumulate(_energy_change_history.begin(),
                    _energy_change_history.end(), 0.0
                    ) / _energy_change_history.size();
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
    bool equilibrium_state_reached(double threshold) const {
        return fabs(get_mean_energy_change()) < threshold;
    };
}; // class PCPVertex

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPVERTEX_HH
