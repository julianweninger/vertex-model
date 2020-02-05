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
template<bool periodic_bc>
class PCPVertex:
    public Model<PCPVertex<periodic_bc>, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<PCPVertex<periodic_bc>, ModelTypes>;

    /// Data type for the model time
    using Time = typename ModelTypes::Time;


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

    /// The domain size in x
    double _Lx;
    /// The domain size in y
    double _Ly;

    /// Linetension constant Lambda
    double _linetension;
    /// Edges shorter than this value are replaced in a T1 transition
    double _length_threshold; 
    
    /// Area elasticity constant K
    double _area_elasticity;
    /// The prefrentrial area of a cell
    double _area_preferential;
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

    /// The total energy in the last step
    double _energy_previous_step;

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
        _Lx(100.), _Ly(100.),
        _linetension(get_as<double>("linetension", this->_cfg)),
        _length_threshold(get_as<double>("length_threshold", this->_cfg)),
        _area_elasticity(get_as<double>("area_elasticity", this->_cfg)),
        _area_preferential(get_as<double>("area_preferential", this->_cfg)),
        _area_threshold(get_as<double>("area_threshold", this->_cfg)),
        _contractility(get_as<double>("contractility", this->_cfg)),
        _prob_distr(0.,1.),
        _distr_noise_const(0., get_as<double>("noise_constant", this->_cfg)),
        _distr_noise_linear(0., get_as<double>("noise_linear", this->_cfg)),
        _energy_linetension(0.),
        _energy_areaelasticity(0.),
        _energy_contractility(0.),
        _energy_previous_step(0.)
    {
        this->initialise_hexagonal(
            get_as<double>("hexagon_size", this->_cfg),
            get_as<int>("lattice_rows", this->_cfg),
            get_as<int>("lattice_columns", this->_cfg)
        );

        // calculate edge lengths
        for (auto &e : _edges) {
            e->template update_length<periodic_bc>(_Lx, _Ly);
        }

        // calculate cell area
        for (auto &c : _cells) {
            c->template cell_area<periodic_bc>();
        }

        // line tension
        _energy_linetension = 0;
        for (auto &e : _edges) {
            _energy_linetension += line_tension(e);
        }
        
        // area elasticity
        _energy_areaelasticity = 0;
        for (auto &c : _cells) {
            _energy_areaelasticity += area_elasticity(c);
        }
        
        // area elasticity
        _energy_contractility = 0;
        for (auto &c : _cells) {
            _energy_contractility += contractility(c);
        }

        this->_log->info("Model initialized.");
    }


private:
    // .. Setup functions .....................................................
    
    // https://courses.cs.washington.edu/courses/cse326/00wi/projects/voronoi.html

    /** Initialiser for randomly distributed cells.
     * Use a Voronoi decomposition from randomly distributed cell centres
     * 
     */
    void initialise_voronoi(int num_cells) {
        throw std::logic_error("Initialise Voronoi function not yet "
            "implemented.");
    }

    /** Initialiser for a hexagonal arrangement of cells.
     *  Initializes a odd-r horizontal layout of hexagonal cells
     * 
     *  @param  size        The size (area) of a single cell
     *  @param  num_rows    The number of cells in a row.
     *                      Note, that this number has to be even in a periodic
     *                      setup
     *  @param  num_columns The number of cells in a column
     *  
     *  Cells
     *      - The cells of even row have the following locations
     *          - center = [(q + 1)*width, (1.5*r + 0.5) * height],
     *              with q and r the column and row id     * 
     *      - The cells of odd rows have the follwing locations
     *          - center = [(q+0.5)*width, (1.5*r + 1.25) * height]     * 
     *      - The id of any cell is given as
     *          - c_id = q + r*num_columns;
     *              notice, that in non-periodic boundary condition one additional
     *              column and row is simulated
     *  
     *  For every cell 2 vertices are created
     *      - for even rows these are
     *          - [q*width, (0.75 * r + 0.25) * height))]
     *          - [(q+0.5)*width, 0.75 * r * height))]
     * 
     *      - for odd rows
     *          - [q*width, r * 0.75 * height) )]
     *          - [(q+0.5)*width, (r * 0.75 + 0.25) * height) )]
     * 
     *      - Thereby the ids of the 6 vertices of every cell have the following ids
     *          - bottom left     2*c_id(q, r)
     *          - bottom          2*c_id(q, r) + 1
     *          - bottom right    2*c_id(q+1, r)
     *          - top left        2*c_id(q, r+1)
     *          - top             2*c_id(q, r+1) + 1
     *          - top right       2*c_id(q+1, r+1)
     * 
     *      - where e.g. c_id((q+1)%lim_columns, r) denotes the id of the cell
     *          to the right, which might be a periodic boundary
     *          lim_columns = num_columns in periodic bc, but 
     *          lim_columns = num_columns+1 otherwise, with the additional
     *          last cell of the row having periodic bc
     * 
     *  For every cell 3 edges are created
     *      - pair rows
     *          - The edges are created as the lower left, lower right and
     *          left edge of every cell.
     *          - Their enumeration is in this order, thus
     *              e_bottom_left = 3 * (q + r*num_cols)
     *              e_bottom_right = 3 * (q + r*num_cols) + 1
     *              e_left = 3 * (q + r*num_cols) + 2
     *          - Thereby the remaining edges have the following ids
     *              e_top_left = 3 * (q + (r+1)*num_cols)
     *              e_top_right = 3 * (q + (r+1)*num_cols) + 1
     *              e_right = 3 * ((q+1) + r*num_cols) + 2
     *                  i.e. the left edge of the cell to the right
     *          - Note that the last to edges originate from impair rows
     *      - impair rows
     *          - The edges are created forming a T around the lower left vertex
     *          of every cell in the order of lower left, lower right, and left
     *          - Thus
     *              - e_bottom_left = 3 * (q + r*num_cols) + 1
     *              - e_left = 3 * (q + r*num_cols) + 2
     *              - e_bottom_right = 3 * ((q+1) + r*num_cols)
     *          - And the remaining edges have the following ids
     *              - e_top_left = 3 * (q + (r+1)*num_cols) + 1
     *              - e_right = 3 * ((q+1) + r*num_cols) + 2
     *              - e_top_right = 3 * ((q+1) + (r+1)*num_cols)
     * 
     *  for additional details see https://www.redblobgames.com/grids/hexagons/
     */
    void initialise_hexagonal(double size, int num_rows, int num_columns) {
        double width = sqrt(3) * size;
        double height = 2 * size;

        if constexpr (periodic_bc) {
            _Lx = num_columns * width;
            _Ly = 0.75 * num_rows * height;

            if (num_rows % 2 != 0) {
                throw std::invalid_argument( "\nERROR with periodic boundary "
                    "conditions the hexagonal cell lattice needs pair number "
                    "of rows, but received impair number. Requested "
                    "number of rows was " + std::to_string(num_columns) + "!");
            }

            // set up with relative coordinates
            width = width / _Lx;
            height = height / _Ly;
        }
        else {
            _Lx = (num_columns + 0.5) * width;
            _Ly = (0.75 * num_rows + 0.5) * height;
        }
            
        width = 1. / double(num_columns);
        height = 4. / (3. * num_rows);

        // Add vertices
        // NOTE one more row and column of vertices initialized in
        //      non-periodic bc this will be undone before initializing edges  
        //      and cells
        // NOTE some of these vertices have to be removed eventually
        int lim_rows = num_rows;
        int lim_columns = num_columns;
        if constexpr (not periodic_bc) {
            lim_columns += 1;
            lim_rows += 1;
        }        
        for (int r = 0; r < lim_rows; r += 1) {
            // pair rows
            if (r % 2 == 0) {
            for (int q = 0; q < lim_columns; ++q) {
                _vertices.push_back(std::make_shared<Vertex>(
                    q*width, (0.75 * r + 0.25) * height));
                _vertices.push_back(std::make_shared<Vertex>(
                    (q+0.5)*width, 0.75 * r * height));
            }
            }
            // impair rows
            else {
            for (int q = 0; q < lim_columns; ++q) {
                _vertices.push_back(std::make_shared<Vertex>(
                    q*width, r * 0.75 * height) );
                _vertices.push_back(std::make_shared<Vertex>(
                    (q+0.5)*width, (r * 0.75 + 0.25) * height) );
            }
            }
        }
        // these vertices are not needed 
        if constexpr (not periodic_bc) {
            _vertices[2*(lim_columns - 1) + 1] = nullptr;
            if (num_rows % 2 == 1) {
                _vertices.back() = nullptr;
            }
            else {
                _vertices[2*(lim_rows-1)*lim_columns] = nullptr;
            }
        }

        // Add edges
        for (int r = 0; r < lim_rows; r++) {
            /** pair rows
             *  The edges are created as the lower left, lower right and
             *  left edge of every cell.
             *  Their enumeration is in this order, thus
             *      e_bottom_left = 3 * (q + r*num_cols)
             *      e_bottom_right = 3 * (q + r*num_cols) + 1
             *      e_left = 3 * (q + r*num_cols) + 2
             *  Thereby the remaining edges have the following ids
             *      e_top_left = 3 * (q + (r+1)*num_cols)
             *      e_top_right = 3 * (q + (r+1)*num_cols) + 1
             *      e_right = 3 * ((q+1) + r*num_cols) + 2
             *          i.e. the left edge of the cell to the right
             *  Note that the last to edges originate from impair rows
             */
            if (r % 2 == 0) {
                int c_id; // id of the resp cell
                for (int q = 0; q < lim_columns; ++q) {
                    // id of the cell
                    c_id = q + r * lim_columns;
                    // lower left edge
                    _edges.push_back(std::make_shared<Edge>(
                        _vertices[2 * c_id], _vertices[2*c_id + 1],
                        _linetension));
                    // lower right edge
                    _edges.push_back(std::make_shared<Edge>(
                        _vertices[2 * c_id + 1],
                        _vertices[2*((q+1)%lim_columns + r*lim_columns)],
                        _linetension));
                    // left edge
                    _edges.push_back(std::make_shared<Edge>(
                        _vertices[2 * c_id],
                        _vertices[2*(q + ((r+1)%lim_rows)*lim_columns)],
                        _linetension));
                }
                // delete not needed edges
                if constexpr (not periodic_bc) {
                    if (r == 0) {
                        _edges[3*(lim_columns - 1)] = nullptr;
                    }
                    
                    _edges[3*((r+1)*lim_columns - 1) + 1] = nullptr;

                    if (r == lim_rows - 1) {
                        _edges[3*(r*lim_columns)] = nullptr;
                        for (int q = 0; q < lim_columns; ++q) {
                            _edges[3*(q + r*lim_columns) + 2] = nullptr;
                        }
                    }
                }
            }            
            /** impair rows
             *  The edges are created forming a T around the lower left vertex
             *  of every cell in the order of lower left, lower right, and left
             *  Thus
             *      e_bottom_left = 3 * (q + r*num_cols) + 1
             *      e_left = 3 * (q + r*num_cols) + 2
             *      e_bottom_right = 3 * ((q+1) + r*num_cols)
             *  And the remaining edges have the following ids
             *      e_top_left = 3 * (q + (r+1)*num_cols) + 1
             *      e_right = 3 * ((q+1) + r*num_cols) + 2
             *      e_top_right = 3 * ((q+1) + (r+1)*num_cols)
             */
            else {
                // handle last column separately
                int c_id;
                for (int q = 0; q < lim_columns; ++q) {
                    // id of the cell
                    c_id = q + r * lim_columns;
                    // lower left edge
                    _edges.push_back(std::make_shared<Edge>(
                        _vertices[2*c_id],
                        _vertices[2*c_id + 1],
                        _linetension));
                    // lower right edge
                    // connects the lower vertex with
                    // the lower left vertex of the cell to the right
                    _edges.push_back(std::make_shared<Edge>(
                        _vertices[2 * c_id + 1],
                        _vertices[2*((q+1)%lim_columns + r*lim_columns)],
                        _linetension));
                    // left edge connects the lower left vertex with
                    // the upper vertex of the cell in the row above
                    _edges.push_back(std::make_shared<Edge>(
                        _vertices[2 * c_id + 1],
                        _vertices[2*(q + ((r+1)%lim_rows)*lim_columns) + 1],
                        _linetension));
                }
                // delete not needed edges
                if constexpr (not periodic_bc) {
                    _edges[3*((r+1)*lim_columns - 1) + 1] = nullptr;

                    if (r == lim_rows - 1) {
                        _edges[3*(lim_rows*lim_columns - 1)] = nullptr;

                        for (int q = 0; q < lim_columns; ++q) {
                            _edges[3*(q + r*lim_columns) + 2] = nullptr;
                        }
                    }
                }
            }
        }

        // ** add cells
        // handle last row separately
        for (int r = 0; r < num_rows; r++) {
            // pair rows
            if (r % 2 == 0) {
                for (int q = 0; q < num_columns; q++) {
                    auto center = std::make_shared<Site>(
                                        (q + 0.5)*width,
                                        (0.75*r + 0.5) * height);
                    EdgeContainer edges;
                    int c_id = q + r * lim_columns;
                    edges.push_back(_edges[3*c_id]); // lower left
                    edges.push_back(_edges[3*c_id + 1]); // lower right
                    edges.push_back(_edges[3*c_id + 2]); // left
                    // right
                    edges.push_back(_edges[3*((q+1)%lim_columns + r*lim_columns) + 2]);

                    edges.push_back(_edges[3*(c_id + lim_columns)]); // upper left
                    edges.push_back(_edges[3*(c_id + lim_columns) + 1]); // upper right
                    // NOTE a pair row cannot be at periodic boundary
                    //      hence no need to handle periodicity

                    _cells.push_back(std::make_shared<Cell>(*center, edges,
                                                            _area_preferential,
                                                            _contractility));
                    _cells.back()->link_members();                    
                }
            }
            else { // impare rows
                for (int q = 0; q < num_columns; q++) {
                    auto center = std::make_shared<Site>(
                                        (q+1)*width,
                                        (0.75*r + 0.5) * height);
                    EdgeContainer edges;
                    int c_id = q + r * lim_columns;
                    edges.push_back(_edges[3*c_id + 1]); // lower left
                    edges.push_back(_edges[3*c_id + 2]); // left
                    // lower right
                    edges.push_back(_edges[3*((q+1)%lim_columns + r*lim_columns)]); 
                    
                    // right
                    edges.push_back(_edges[3*((q+1)%lim_columns + r*lim_columns) + 2]);

                    // upper left
                    edges.push_back( _edges[3*(q + (r+1)%lim_rows * lim_columns) + 1]);

                    // upper right
                    edges.push_back(_edges[3*((q+1)%lim_columns + (r+1)%lim_rows * lim_columns)]);

                    _cells.push_back(std::make_shared<Cell>(*center, edges,
                                                            _area_preferential,
                                                            _contractility));
                    _cells.back()->link_members();
                }
            }
        }

        // remove expired objects
        _vertices.erase(
            std::remove_if(_vertices.begin(), _vertices.end(),
                           [](auto v) { return v == nullptr; }),
            _vertices.end());
        _edges.erase(
            std::remove_if(_edges.begin(), _edges.end(),
                           [](auto e) { return e == nullptr; }),
            _edges.end());

        for (auto e : _edges) {
            e->link_members();
        }

        this->_log->info("Initialised hexagonal cells.");
    }

    // .. Helper functions ....................................................
    /// Resets the forces of this vertex
    std::function<void(Vertex_ptr&)> reset_forces = [](Vertex_ptr &v) {
        v->fx = 0.;
        v->fy = 0.;
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

        const double fx = e->linetension*dx/e->length;
        const double fy = e->linetension*dy/e->length;

        e->a->fx += fx;
        e->a->fy += fy;
        e->b->fx -= fx;
        e->b->fy -= fy;

        return e->linetension * e->length;
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
            // the edge prior the vertex
            auto [e0, e0_flip]  = c->edges_ordered[std::max(0, edges_it - 1)];
            if (edges_it == 0) { std::tie(e0, e0_flip) = c->edges_ordered.back(); }

            // the edge post the vertex
            const auto [e1, e1_flip] = c->edges_ordered[edges_it];
            
            // vertices in ordering
            Vertex_ptr v_center = e0->b, v_prior = e0->a;
            if (e0_flip) {
                std::swap(v_center, v_prior);
            }

            Vertex_ptr v_post;
            if (e1_flip) {
                v_post = e1->a;
            }
            else {
                v_post = e1->b;
            }

            auto [dx1, dy1] = displacement_absolute<periodic_bc>(*v_center,
                                    *v_prior, _Lx, _Ly);
            auto [dx2, dy2] = displacement_absolute<periodic_bc>(*v_post,
                                    *v_center, _Lx, _Ly);

            double dA_dx = -0.5 * (dy1 + dy2) * c->area_sgn;
            double dA_dy = 0.5 * (dx1 + dx2) * c->area_sgn;
            
            // this is the force on this vertex
            v_center->fx -= _area_elasticity * (c->area_abs(_Lx, _Ly) - c->area_preferential)*dA_dx;
            v_center->fy -= _area_elasticity * (c->area_abs(_Lx, _Ly) - c->area_preferential)*dA_dy;
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
            perimeter += e->length;
        }

        for (auto [e, flip] : c->edges_ordered) {
            auto a = e->a; auto b = e->b;
            if (flip) {
                std::swap(a, b);
            }
            auto [dx, dy] = displacement_absolute<periodic_bc>(*a, *b, _Lx, _Ly);
            
            double fx = c->contractility * perimeter * dx / e->length;
            double fy = c->contractility * perimeter * dy / e->length;

            a->fx += fx;
            a->fy += fy;
            b->fx += -fx;
            b->fy += -fy;
        }

        return 0.5 *c->contractility * std::pow(perimeter, 2);
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

    /** Erase cell from _cells while updating the topology (T2 transition)
     * 
     *  Removes the cell, its edges and its vertices and sets up a vertex at its
     *  centre.
     * 
     *  returns _cells.erase(cell_it)
     */
    CellContainer::iterator T2_transition (CellContainer::iterator &cell_it) 
    {
        this->_log->debug("Removing cell in T2 transition..");

        auto cell = *cell_it;
        
        // create a new vertex at the center of c
        cell->template cell_area<periodic_bc>(); // updates the center of c        
        auto new_v = std::make_shared<Vertex>(cell->s);
        _vertices.push_back(new_v);

        // Tag the objects that are to be removed
        cell->remove = true;
        for (auto &v : cell->vertices) {
            v->remove = true;
        }
        for (auto &e_pair : cell->edges_ordered) {
            std::get<Edge_ptr>(e_pair)->remove = true;
        }

        // update cells
        for (auto c_it = _cells.begin(); c_it != _cells.end(); /*void*/) {
            auto c = *c_it;
            
            if (c->remove) {
                ++c_it;
                continue;
            }

            int num_vertices = c->vertices.size();

            // remove links to removed vertexes
            c->vertices.erase(
                std::remove_if(
                    c->vertices.begin(), c->vertices.end(),
                    [](auto &v) { return v->remove; }),
                c->vertices.end()
            );

            // remove links to removed edges
            c->edges_ordered.erase(
                std::remove_if(
                    c->edges_ordered.begin(), c->edges_ordered.end(),
                    [](auto &e_pair) { 
                        return std::get<Edge_ptr>(e_pair)->remove; }),
                c->edges_ordered.end()
            );
            // No new edge to add
            
            // NOTE edges will be still ordered once the edges are updates

            // add new_v and update area
            if (num_vertices > c->vertices.size()) {
                c->vertices.push_back(new_v);
                c->template cell_area<periodic_bc>();
            }

            // continue iteration
            ++c_it;
        }

        // update edges
        for (auto e_it = _edges.begin(); e_it != _edges.end(); /*void*/) {
            auto e = *e_it;
            
            if (e->remove) {
                e_it = _edges.erase(e_it);
                continue;
            }

            // replace vertices that have been removed
            if (e->a->remove) {
                e->a = new_v;
                e->template update_length<periodic_bc>(_Lx, _Ly);
            }
            else if (e->b->remove) {
                e->b = new_v;
                e->template update_length<periodic_bc>(_Lx, _Ly);
            }

            ++e_it;
        }

        // remove vertices
        _vertices.erase(
            std::remove_if(
                _vertices.begin(), _vertices.end(),
                [](auto v) { return v->remove; }),
            _vertices.end()
        );

        // remove cell
        return _cells.erase(cell_it);
    }

    /** Perform a T1 transition on the edge edge_it in _edges
     * 
     *  edge_it will be removed and a new edge orthogonal to edge will be created
     */
    EdgeContainer::iterator T1_transition (EdgeContainer::iterator edge_it)
    {
        this->_log->debug("Removing edge in T1 transition..");

        auto edge = *edge_it;

        // tag objects to be removed
        edge->remove = true;
        edge->a->remove = true;
        edge->b->remove = true;

        // the adj cells that are currently neighbours sharing edge
        // NOTE a, b arbitrary
        Cell_ptr adj_cell_a = edge->adj_cell_a.lock();
        Cell_ptr adj_cell_b = edge->adj_cell_b.lock();

        // the two cells that become neighbours in this T1 transition
        // NOTE c is the cell adjoint to vertex edge->a (arbitrary)
        //      d is the cell adjoint to vertex edge->b
        Cell_ptr adj_cell_c, adj_cell_d;
        if (edge->a->adj_cells.size() != 3) {
            throw std::runtime_error("Not implemented error: In T1 transition "
                "(edge length < threshold) expected 3 adj_cells to a vertex, "
                "but encountered " + 
                std::to_string(edge->a->adj_cells.size()) + "! "
                "This would be e.g. an edge orthogonal to the outer boundary of "
                "a collection of cells.");
        }
        if (edge->b->adj_cells.size() != 3) {
            throw std::runtime_error("Not implemented error: In T1 transition "
                "(edge length < threshold) expected 3 adj_cells to a vertex, "
                "but encountered " + 
                std::to_string(edge->b->adj_cells.size()) + "! "
                "This would be e.g. an edge orthogonal to the outer boundary of "
                "a collection of cells.");
        }
        for (auto &c_weak : edge->a->adj_cells) {
            auto c = c_weak.lock();
            if (c != adj_cell_a and c != adj_cell_b) {
                adj_cell_c = c;
            }
        }
        for (auto &c_weak : edge->b->adj_cells) {
            auto c = c_weak.lock();
            if (c != adj_cell_a and c != adj_cell_b) {
                adj_cell_d = c;
            }
        }

        // The involved edges
        // (a, b) and (c, d) currently share a vertex
        // (a, c) and (b, c) currently share an adjoint cell a, resp. b
        // they will share a vertex a (resp. b) after transition
        Edge_ptr adj_edge_a, adj_edge_b, adj_edge_c, adj_edge_d;
        if (edge->a->adj_edges.size() != 3 or edge->b->adj_edges.size() != 3) {
            throw std::runtime_error("Expected 3 adj_edges to a vertex, but "
                    "encountered " + std::to_string(edge->a->adj_edges.size())
                    + "!");
        }
        if (edge->b->adj_edges.size() != 3) {
            throw std::runtime_error("Expected 3 adj_edges to a vertex, but "
                    "encountered " + std::to_string(edge->b->adj_edges.size())
                    + "!");
        }
        for (auto &e_weak : edge->a->adj_edges) {
            auto e = e_weak.lock();
            if (e != edge) {
                if (e->adj_cell_a.lock() == adj_cell_a or
                    e->adj_cell_b.lock() == adj_cell_a)
                {
                    adj_edge_a = e;
                }
                else {
                    adj_edge_b = e;
                }
            }
        }
        for (auto &e_weak : edge->b->adj_edges) {
            auto e = e_weak.lock();
            if (e != edge) {
                if (e->adj_cell_a.lock() == adj_cell_a or
                    e->adj_cell_b.lock() == adj_cell_a)
                {
                    adj_edge_c = e;
                }
                else {
                    adj_edge_d = e;
                }
            }
        }

        // create two new vertices that create an edge of threshold length 
        // pointing from cell a to b
        auto [dx, dy] = displacement_absolute<periodic_bc>(*adj_cell_a->s,
                                                           *adj_cell_b->s, 
                                                           _Lx, _Ly);
        auto length = distance<periodic_bc>(*adj_cell_a->s, 
                                            *adj_cell_b->s, _Lx, _Ly);
        dx = dx / length * _length_threshold / _Lx;
        dy = dy / length * _length_threshold / _Ly;
        auto new_v_a = std::make_shared<Vertex>(
                                0.5 * (edge->a->x + edge->b->x) - dx / 2.,
                                0.5 * (edge->a->y + edge->b->y) - dy / 2.,
                                EdgeContainer({adj_edge_a, adj_edge_c}),
                                CellContainer({adj_cell_a, adj_cell_c,
                                               adj_cell_d}));
        auto new_v_b = std::make_shared<Vertex>(
                                0.5 * (edge->a->x + edge->b->x) + dx / 2.,
                                0.5 * (edge->a->y + edge->b->y) + dy / 2.,
                                EdgeContainer({adj_edge_b, adj_edge_d}),
                                CellContainer({adj_cell_b, adj_cell_c,
                                               adj_cell_d}));
        correct_periodic_bc<periodic_bc>(new_v_a);
        correct_periodic_bc<periodic_bc>(new_v_b);
        _vertices.push_back(new_v_a);
        _vertices.push_back(new_v_b);

        // create a new edge
        Edge_ptr new_edge = std::make_shared<Edge>(new_v_a, new_v_b,
                                                   _linetension, adj_cell_c,
                                                   adj_cell_d);
        new_edge->link_members();

        // remove objects
        _vertices.erase(
            std::remove_if(
                _vertices.begin(), _vertices.end(),
                [](auto v) { return v->remove; }),
            _vertices.end()
        );
        for (auto &c : {adj_cell_a, adj_cell_b, adj_cell_c, adj_cell_d}) {
            c->vertices.erase(
                std::remove_if(
                    c->vertices.begin(), c->vertices.end(), 
                    [](auto v) { return v->remove; }),
                c->vertices.end()
            );
        }
        // NOTE edge_it will be removed at the very end
        // NOTE a and b will be replaced within edges

        // remove edge from adj_cells a and b
        // NOTE the neighbouring edges now have a common vertex, hence order of
        //      edges is maintained
        for (auto &c : {adj_cell_a, adj_cell_b}) {
            c->edges_ordered.erase(
                std::remove_if(c->edges_ordered.begin(),
                               c->edges_ordered.end(), 
                               [](auto e_pair) {
                                    return std::get<Edge_ptr>(e_pair)->remove; }),
                c->edges_ordered.end()
            );
        }
        
        // replace vertices in edges
        // NOTE a, c share cell a; b, d share cell b
        //      hence, new_v_a associated with cell a
        //      and new_v_b associated with cell b
        for (auto &e : {adj_edge_a, adj_edge_c}) {
            if (e->a->remove) { e->a = new_v_a; }
            else { e->b = new_v_a; }
        }
        for (auto &e : {adj_edge_b, adj_edge_d}) {
            if (e->a->remove) { e->a = new_v_b; }
            else { e->b = new_v_b; }
        }

        // add new vertices to cells
        // NOTE vertex a is associated with cell a; 
        //      vertex b is associated with cell b
        //      both associated with cells c and d
        for (auto c : {adj_cell_a, adj_cell_c, adj_cell_d}) {
            c->vertices.push_back(new_v_a);
        }
        for (auto c : {adj_cell_b, adj_cell_c, adj_cell_d}) {
            c->vertices.push_back(new_v_b);
        }

        // make a new edge in adj_cells c and d
        // NOTE this is symmetric, directionality is given by void order_edges()
        for (auto &c : {adj_cell_c, adj_cell_d}) {
            c->edges_ordered.push_back(
                std::make_pair(new_edge, false) );
            c->order_edges();
        }

        // update the objects
        adj_cell_a->template cell_area<periodic_bc>();
        adj_cell_b->template cell_area<periodic_bc>();
        adj_cell_c->template cell_area<periodic_bc>();
        adj_cell_d->template cell_area<periodic_bc>();
        
        new_edge->template update_length<periodic_bc>(_Lx, _Ly);
        adj_edge_a->template update_length<periodic_bc>(_Lx, _Ly);
        adj_edge_b->template update_length<periodic_bc>(_Lx, _Ly);
        adj_edge_c->template update_length<periodic_bc>(_Lx, _Ly);
        adj_edge_d->template update_length<periodic_bc>(_Lx, _Ly);

        // replace the edge at adge_it
        edge_it = _edges.erase(edge_it);
        edge_it = _edges.insert(edge_it, new_edge);

        return ++edge_it;
    }

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
                                        double division_angle)
    {
        this->_log->debug("Dividing cell..");

        // The cell to be divided
        auto cell = *cell_it;
        cell->remove = true;
        cell_it = _cells.erase(cell_it);
        auto cell_center = std::make_shared<Vertex>(cell->s);

        // generate the axis of division
        double dx = cos(division_angle) / _Lx;
        double dy = sin(division_angle) / _Ly;
        auto tmp_vertex = std::make_shared<Vertex>(cell_center->x + dx, 
                                                   cell_center->y + dy);
        auto division_axis = Edge(cell_center, tmp_vertex, 0.);

        // determine the new vertices from this axis
        // These are the intersections of the division axis with edges of cell
        // NOTE there must be exactly 2 intersections for nicely shaped cells
        std::vector<Vertex_ptr> new_vertices;
        for (auto e_pair : cell->edges_ordered) {
            auto e = std::get<Edge_ptr>(e_pair);

            // calculate the intersection
            auto inter = intersection<periodic_bc>(*e, division_axis,
                                                   true, false);
            
            // if no intersection end here
            if (not inter) { continue; }

            // create a new vertex at the intersection site
            auto new_v = std::make_shared<Vertex>(inter, EdgeContainer({e}),
                                CellContainer({e->adj_cell_a.lock(),
                                               e->adj_cell_b.lock()}));
            // NOTE the edge link is used later and removed thereafter
            
            // keep track of this vertex
            _vertices.push_back(new_v);
            new_vertices.push_back(new_v);

            // the edge itself will be divided and removed
            e->remove = true;
        }
        if (new_vertices.size() != 2) {
            this->_log->warn("Cannot perform cell division on cell at ({}, {}), "
                "with a division angle of {}.", cell->s->x, cell->s->y,
                division_angle / 2 / PI * 360);
            this->_log->warn("division axis is {}, {} -> {}, {}",
                division_axis.a->x, division_axis.a->y,
                division_axis.b->x, division_axis.b->y);
            this->_log->warn("Edges of cell are:");
            for (auto [e, flip] : cell->edges_ordered) {
                if (flip) {
                    this->_log->warn("   {}, {} -> {}, {}", e->b->x, e->b->y,
                        e->a->x, e->a->y);
                }
                else {
                    this->_log->warn("   {}, {} -> {}, {}", e->a->x, e->a->y,
                        e->b->x, e->b->y);
                }
            }
            this->_log->warn("Intersections are:");
            for (auto v : new_vertices) {
                this->_log->warn("   {}, {}", v->x, v->y);
            }
            throw std::runtime_error("During cell division, expected 2 new "
                "vertices, but got " + std::to_string(new_vertices.size()) + "!");
        }

        // the 2 edges that are divided by the new edge
        auto edge_a = new_vertices[0]->adj_edges.back().lock();
        auto edge_b = new_vertices[1]->adj_edges.back().lock();
        // remove them, they will be replaced by 2 new edges each
        _edges.erase(std::remove_if(_edges.begin(), _edges.end(), 
                                    [](auto e){ return e->remove;}),
                     _edges.end());
        new_vertices[0]->adj_edges.clear();
        new_vertices[1]->adj_edges.clear();

        // create a new edge connecting the 2 new vertices
        // NOTE it has properties as at model initialisation
        auto new_edge = std::make_shared<Edge>(new_vertices[0], 
                                               new_vertices[1], _linetension);
        new_edge->link_members();
        _edges.push_back(new_edge);

        // this is how to divide an edge at a pivot vertex
        /* \param e     The Edge to divide
         * \param flip  The direction of e
         * \param pivot The Vertex where to divide e
         * 
         * \return  the first half of e from a to pivot, the second half of e 
         *          from pivot to b
         */
        auto divide_edge = [cell](Edge_ptr e, bool flip, Vertex_ptr pivot)
        {
            // The start and end of e, considering the direction within the 
            // iteration
            Vertex_ptr a, b;
            if (flip) {
                a = e->b;
                b = e->a;
            }
            else {
                a = e->a;
                b = e->b;
            }

            for (auto v : {a, b}) {
                v->adj_edges.erase(
                    std::remove_if(
                        v->adj_edges.begin(),
                        v->adj_edges.end(),
                        [](auto e){ return e.lock()->remove; }),
                    v->adj_edges.end());
            }

            // The first half of the edge from a to pivot
            auto new_edge_0 = std::make_shared<Edge>(a, pivot, e->linetension,
                                        e->adj_cell_a.lock(),
                                        e->adj_cell_b.lock());
            new_edge_0->link_members();

            // the second half of the edge from pivot to b
            auto new_edge_1 = std::make_shared<Edge>(b, pivot, e->linetension,
                                        e->adj_cell_a.lock(),
                                        e->adj_cell_b.lock());
            new_edge_1->link_members();

            // update the crosslinks in the adj cells
            for (auto c_weak : {e->adj_cell_a, e->adj_cell_b}) {
                if (c_weak.lock() == cell) { continue; }
                else {
                    auto c = c_weak.lock();
                    c->edges_ordered.erase(
                        std::remove_if(
                            c->edges_ordered.begin(), c->edges_ordered.end(),
                            [](auto e_pair){ 
                                return std::get<Edge_ptr>(e_pair)->remove; }
                        ),
                        c->edges_ordered.end()
                    );
                    c->edges_ordered.push_back(std::make_pair(new_edge_0, 
                                                              false));
                    c->edges_ordered.push_back(std::make_pair(new_edge_1, 
                                                              false));
                    
                    c->order_edges();
                }
            }

            return std::make_pair(new_edge_0, new_edge_1);
        };
                
        // separate the edges to form 2 cells
        EdgeContainer new_edges_cell_0, new_edges_cell_1;
        
        // the new edge (division axis) is in both cells
        new_edges_cell_0.push_back(new_edge);
        new_edges_cell_1.push_back(new_edge);

        /* 1. start iteration of ordered edges of cell at random point, add
         *    edges to 1st cell.
         * 2. stop the iteration when edge_a or _b is reached, divide this edge
         *    into 2 new edges at its intersection with the division axis.
         * 3. continue iteration from second half of this divided edge, but add
         *    these edges to the 2nd cell.
         * 4. stop iteration at other edge_b or _a, resp. divide this edge as
         *    before.
         * 5. finish iteration from second half of this divided edge, but again 
         *    add to 1st cell.
         */
        int it_edges;
        // 1. start iteration
        for (it_edges = 0; it_edges < cell->edges_ordered.size(); it_edges++) {
            auto e_pair = cell->edges_ordered[it_edges];
            auto e = std::get<Edge_ptr>(e_pair);
            // check whether edge_a or _b reached
            if (e == edge_a or e == edge_b) {
                // the intersection with the division axis
                Vertex_ptr edge_pivot;
                if (e == edge_a) { edge_pivot = new_vertices[0]; }
                else { edge_pivot = new_vertices[1]; }

                Edge_ptr new_edge_0, new_edge_1;
                tie(new_edge_0, new_edge_1) = divide_edge(e,
                                                    std::get<bool>(e_pair), 
                                                    edge_pivot);

                _edges.push_back(new_edge_0);
                _edges.push_back(new_edge_1);

                new_edges_cell_0.push_back(new_edge_0);
                new_edges_cell_1.push_back(new_edge_1);

                // 2. stop iteration here
                break;
            }
            else {
                // (1.) this edge is in 1st cell.
                new_edges_cell_0.push_back(e);
            }
        }
        // 3. continue iteration for edges in 2nd cell
        for (it_edges += 1; it_edges < cell->edges_ordered.size(); it_edges++) {
            auto e_pair = cell->edges_ordered[it_edges];
            auto e = std::get<Edge_ptr>(e_pair);
            // check whether edge_a or _b reached
            if (e == edge_a or e == edge_b) {
                // the intersection with the division axis
                Vertex_ptr edge_pivot;
                if (e == edge_a) { edge_pivot = new_vertices[0]; }
                else { edge_pivot = new_vertices[1]; }

                Edge_ptr new_edge_0, new_edge_1;
                tie(new_edge_0, new_edge_1) = divide_edge(e,
                                                    std::get<bool>(e_pair), 
                                                    edge_pivot);

                _edges.push_back(new_edge_0);
                _edges.push_back(new_edge_1);

                new_edges_cell_1.push_back(new_edge_0);
                new_edges_cell_0.push_back(new_edge_1);

                // 4. stop iteration here
                break;
            }
            else {
                // (3.) this edge is in 2nd cell.
                new_edges_cell_1.push_back(e);
            }
        }
        // 5. finish iteration for cells in 1st cell
        for (it_edges += 1; it_edges < cell->edges_ordered.size(); it_edges++) {
            auto e_pair = cell->edges_ordered[it_edges];
            auto e = std::get<Edge_ptr>(e_pair);
            new_edges_cell_0.push_back(e);
        }

        // create 2 new cells
        auto new_cell_0 = std::make_shared<Cell>(*cell_center, new_edges_cell_0, 
                                                 cell->area_preferential,
                                                 cell->contractility);
        new_cell_0->link_members();
        auto new_cell_1 = std::make_shared<Cell>(*cell_center, new_edges_cell_1, 
                                                 cell->area_preferential,
                                                 cell->contractility);
        new_cell_1->link_members();
        _cells.push_back(new_cell_0);
        _cells.push_back(new_cell_1);
        
        // remove expired crosslinks
        for (auto new_c : {new_cell_0, new_cell_1}) {
            for (auto v : new_c->vertices) {
                v->adj_cells.erase(std::remove_if(v->adj_cells.begin(),
                                                v->adj_cells.end(),
                                                [cell](auto c){
                                                    return c.lock() == cell;}),
                                v->adj_cells.end());
            }
            for (auto e_pair : new_c->edges_ordered) {
                auto e = std::get<Edge_ptr>(e_pair);
                if (e->adj_cell_a.lock() == cell) {
                    e->adj_cell_a = new_c;
                }
                else if (e->adj_cell_b.lock() == cell) {
                    e->adj_cell_b = new_c;
                }
            }
        }

        return cell_it;
    }

public:
    // -- Public Interface ----------------------------------------------------
    // .. Simulation Control ..................................................
    
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
        _energy_previous_step = this->get_energy();

        // reset forces
        std::for_each(_vertices.begin(), _vertices.end(), reset_forces);

        // calculate edge lengths
        for (auto &e : _edges) {
            e->template update_length<periodic_bc>(_Lx, _Ly);

            // breakpoint in debug mode
            #ifndef NDEBUG
                // periodic bc only defined if edge defines shortest distance
                // between two points 
                if (periodic_bc and e->length > 0.9 * 0.25 * sqrt(_Lx*_Ly)) {
                    auto [dx, dy] = displacement_absolute<periodic_bc>(*e->a,
                                        *e->b, _Lx, _Ly);
                    if (dx > 0.9 * _Lx) {
                        throw std::runtime_error("Edge is supposed to connect "
                            "two vertices, that are separated by " +
                            std::to_string(dx) + " on a periodic domain, "
                            "that measures " + std::to_string(_Lx) + " in x. "
                            "This is only defined as long as the edge defines "
                            "the shortest path between the two vertices. "
                            "Hence it is no longer defined for dx -> Lx / 2");
                    }
                    if (dy > 0.9 * _Ly) {
                        throw std::runtime_error("Edge is supposed to connect "
                            "two vertices, that are separated by " +
                            std::to_string(dy) + " on a periodic domain, "
                            "that measures " + std::to_string(_Ly) + " in y. "
                            "This is only defined as long as the edge defines "
                            "the shortest path between the two vertices. "
                            "Hence it is no longer defined for dy -> Ly / 2");
                    }
                }
            #endif
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
            if ((*e_it)->length < _length_threshold)
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
    }


    /// Monitor model information
    void monitor () { }


    /// Writes the initial data
    /** This function is taken care of in the run mode, but not in an iterated
     *  mode
     */ 
    void write_data_initial () {
        this->_datamanager(*this);
    }


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

    /// Getter for energy
    double get_energy () const {
        return _energy_linetension + _energy_areaelasticity + 
               _energy_contractility;
    }

    /// Getter for normalised energy
    /** The energy is normalised wrt number of vertices, edges, or cells, 
     *  respectively.
     */
    double get_energy_normalised () const {
        return get_energy_linetension_normalised() + 
               get_energy_areaelasticity_normalised() +
               get_energy_contractility_normalised();
    }

    /// Getter for the relative energy change from previous to last step
    double get_rel_energy_change () const {
        double energy_change = fabs(get_energy() - _energy_previous_step);
        return energy_change / get_energy();
    }

    /// Getter for the domain size
    const std::pair<double, double> get_domain_size () const {
        return std::make_pair(_Lx, _Ly);
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

    /// Increase the domain size by a certain area
    /** This remaps the domain of size A to size A + dA while keeping the 
     *  relation of Lx to Ly constant.
     * 
     *  Thereby proliferation of cells can be performed in a periodic setup
     *  without changing the parameters of the system. 
     */
    void increase_domain_size(double area) {
        double ratio = _Lx / double(_Ly);
        _Ly = std::sqrt(_Ly*_Ly + area / ratio);
        _Lx = ratio * _Ly;
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
     *  Equilibrium if relative energy change is smaller than threshold
     *  
     *  \param threshold    The equilibrium threshold
     */
    bool equilibrium_state_reached(double threshold) const {
        return get_rel_energy_change() < threshold;
    };
}; // class PCPVertex

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPVERTEX_HH
