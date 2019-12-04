#ifndef UTOPIA_MODELS_PCPVERTEX_HH
#define UTOPIA_MODELS_PCPVERTEX_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>

#include "geometry.hh"


namespace Utopia {
namespace Models {
namespace PCPVertex {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<>;


// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// The PCPVertex Model; the bare-basics a model needs
/** TODO Add your class description here.
 *  ...
 */
template<bool periodic_bc>
class PCPVertex:
    public Model<PCPVertex<periodic_bc>, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<PCPVertex<periodic_bc>, ModelTypes>;

    /// Data type of the group to write model data to, holding datasets
    using DataGroup = typename Base::DataGroup;

    /// Data type for a dataset
    using DataSet = typename Base::DataSet;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    VertexContainer _vertices;
    EdgeContainer _edges;
    CellContainer _cells;

    // PARAMETERS
    double _dt = 0.1;    /// timestep scaling

    double _linetension = 10.;
    double _area_elasticity = 0.01;
    double _area_preferential = 1000.;

    // .. Temporary objects ...................................................


    // .. Datasets ............................................................
    // NOTE They should be named '_dset_<name>', where <name> is the
    //      dataset's actual name as set in its constructor. Ideally, do not
    //      hide them inside a struct ...
    const std::shared_ptr<DataGroup> _grp_vertices;
    const std::shared_ptr<DataGroup> _grp_edges;
    const std::shared_ptr<DataGroup> _grp_cells;    
    const std::shared_ptr<DataSet> _dset_forces;

public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the PCPVertex model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel>
    PCPVertex (const std::string name, ParentModel& parent)
    :
        // Initialize first via base model
        Base(name, parent),
        // Get member paramters from cfg
        _vertices(),
        _edges(),
        _cells(),
        
        _dt(get_as<double>("dt", this->_cfg)),

        _linetension(get_as<double>("linetension", this->_cfg)),
        _area_elasticity(get_as<double>("area_elasticity", this->_cfg)),
        _area_preferential(get_as<double>("area_preferential", this->_cfg)),

        // Open the datasets
        // e.g. via _dset_state(this->create_dset("state", {})) <- 1d
        //      or  _dset_state(this->create_dset("state", {num_states})) <- 2d
        _grp_vertices(this->_hdfgrp->open_group("vertices")),
        _grp_edges(this->_hdfgrp->open_group("edges")),
        _grp_cells(this->_hdfgrp->open_group("cells")),
        _dset_forces(this->create_dset("forces", {}))
    {
        this->initialise_hexagonal(
            get_as<double>("hexagon_size", this->_cfg),
            get_as<int>("lattice_rows", this->_cfg),
            get_as<int>("lattice_columns", this->_cfg)
        );

        this->_log->info("Model initialized.");
    }


private:
    // .. Setup functions .....................................................
/// https://courses.cs.washington.edu/courses/cse326/00wi/projects/voronoi.html
    void initialise()
    {
        this->_log->debug("Initialising ..");
    }

    /** Initializes a odd-r horizontal layout
     *  https://www.redblobgames.com/grids/hexagons/
     *  
     *  odd-r horizontal layout
     * 
     *  The cells of even row have the following locations
     *      center = [(q + 1)*width, (1.5*r + 0.5) * height]
     *      with q and r the column and row id
     *  The cells of odd rows have the follwing locations
     *      center = [(q+0.5)*width, (1.5*r + 1.25) * height]
     *  The id of any cell is given as
     *      c_id = q + r*num_columns
     *      notice, that in non-periodic boundary condition one additional
     *      column and row is simulated
     *  
     *  For every cell 2 vertices are created, for even rows these are
     *      [q*width, (0.75 * r + 0.25) * height))]
     *      [(q+0.5)*width, 0.75 * r * height))]
     *  and for odd rows
     *      [q*width, r * 0.75 * height) )]
     *      [(q+0.5)*width, (r * 0.75 + 0.25) * height) )]
     *  Thereby the ids of the 6 vertices of every cell have the following ids
     *      bottom left     2*c_id(q, r)
     *      bottom          2*c_id(q, r) + 1
     *      bottom right    2*c_id(q+1, r)
     *      top left        2*c_id(q, r+1)
     *      top             2*c_id(q, r+1) + 1
     *      top right       2*c_id(q+1, r+1)
     * 
     *          where e.g. c_id((q+1)%lim_columns, r) denotes the id of the cell
     *          to the right, which might be a periodic boundary
     *          lim_columns = num_columns in periodic bc, but 
     *          lim_columns = num_columns+1 otherwise, with the additional
     *          last cell of the row having periodic bc
     * 
     *  For every cell 3 edges are created
     *      pair rows
     *          The edges are created as the lower left, lower right and
     *          left edge of every cell.
     *          Their enumeration is in this order, thus
     *              e_bottom_left = 3 * (q + r*num_cols)
     *              e_bottom_right = 3 * (q + r*num_cols) + 1
     *              e_left = 3 * (q + r*num_cols) + 2
     *          Thereby the remaining edges have the following ids
     *              e_top_left = 3 * (q + (r+1)*num_cols)
     *              e_top_right = 3 * (q + (r+1)*num_cols) + 1
     *              e_right = 3 * ((q+1) + r*num_cols) + 2
     *                  i.e. the left edge of the cell to the right
     *          Note that the last to edges originate from impair rows
     *      impair rows
     *          The edges are created forming a T around the lower left vertex
     *          of every cell in the order of lower left, lower right, and left
     *          Thus
     *              e_bottom_left = 3 * (q + r*num_cols) + 1
     *              e_left = 3 * (q + r*num_cols) + 2
     *              e_bottom_right = 3 * ((q+1) + r*num_cols)
     *          And the remaining edges have the following ids
     *              e_top_left = 3 * (q + (r+1)*num_cols) + 1
     *              e_right = 3 * ((q+1) + r*num_cols) + 2
     *              e_top_right = 3 * ((q+1) + (r+1)*num_cols)
     *      
     */
    void initialise_hexagonal(double size, int num_rows, int num_columns) {
        this->_log->debug("Initialising hexagonal cells ..");

        double height = 2 * size;
        double width = sqrt(3) * size;

        if constexpr (periodic_bc) {
            this->_log->warn("Resizing the domain!!");
            
            Lx = num_columns * width;
            Ly = 0.75 * num_rows * height;

            if (num_rows % 2 != 0) {
                throw std::invalid_argument( "\nERROR with periodic boundary "
                    "conditions the hexagonal cell lattice needs pair number "
                    "of rows, but received impair number. Requested "
                    "number of rows was " + std::to_string(num_columns) + "!");
            }
        }

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
                    // y = ((r+1)*0.75 + 0.25) * h
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
                        _vertices[2 * c_id], _vertices[2*c_id + 1] ));
                    // lower right edge
                    _edges.push_back(std::make_shared<Edge>(
                        _vertices[2 * c_id + 1],
                        _vertices[2*((q+1)%lim_columns + r*lim_columns)]));
                    // left edge
                    _edges.push_back(std::make_shared<Edge>(
                        _vertices[2 * c_id],
                        _vertices[2*(q + ((r+1)%lim_rows)*lim_columns)]));
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
                        _vertices[2*c_id + 1]));
                    // lower right edge
                    // connects the lower vertex with
                    // the lower left vertex of the cell to the right
                    _edges.push_back(std::make_shared<Edge>(
                        _vertices[2 * c_id + 1],
                        _vertices[2*((q+1)%lim_columns + r*lim_columns)]));
                    // left edge connects the lower left vertex with
                    // the upper vertex of the cell in the row above
                    _edges.push_back(std::make_shared<Edge>(
                        _vertices[2 * c_id + 1],
                        _vertices[2*(q + ((r+1)%lim_rows)*lim_columns) + 1]));
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

                    _cells.push_back(std::make_shared<Cell>(center, edges));
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

                    _cells.push_back(std::make_shared<Cell>(center, edges));
                }
            }
        }

        for (auto it = _edges.begin(); it != _edges.end(); /*void*/) {
            if (!it->get()) { // nullptr here
                it = _edges.erase(it);
            }
            else {
                ++it; 
            }
        }
        for (auto it = _vertices.begin(); it != _vertices.end(); /*void*/) {
            if (!it->get()) { // nullptr here
                it = _vertices.erase(it);
            }
            else {
                ++it; 
            }
        }

        this->_log->info("Initialised hexagonal cells.");
    }

    // .. Helper functions ....................................................

public:
    // -- Public Interface ----------------------------------------------------
    // .. Simulation Control ..................................................

    /// Iterate a single step
    /** \details Here you can add a detailed description what exactly happens 
      *         in a single iteration step
      */
    void perform_step () {
        // reset forces
        for (auto v : _vertices) {
            v->fx = 0.;
            v->fy = 0.;
        }

        // line tension
        for (auto e : _edges) {
            double dx = e->b->x - e->a->x;
            double dy = e->b->y - e->a->y;

            if constexpr (periodic_bc) {
                if (dx <= -Lx / 2.) { dx = dx + Lx; }
                else if (dx > Lx / 2.) { dx = dx - Lx; }

                if (dy <= -Ly / 2.) { dy = dy + Ly; }
                else if (dy > Ly / 2.) { dy = dy - Ly; }
            }

            e->length = sqrt(pow(dx, 2.) + pow(dy, 2.));

            const double fx = _linetension*dx/e->length;
            const double fy = _linetension*dy/e->length;

            e->a->fx += fx;
            e->a->fy += fy;
            e->b->fx -= fx;
            e->b->fy -= fy;
        }
        
        // area elasticity
        int cnt = 0;
        for (auto c : _cells) {
            double area = c->cell_area();

            for (int edges_it = 0; edges_it < c->edges_ordered.size(); edges_it++) {
                auto e0_pair  = c->edges_ordered[std::max(0, edges_it - 1)];
                if (edges_it == 0) { e0_pair = c->edges_ordered.back(); }
                const auto e1_pair = c->edges_ordered[edges_it];
                
                // vertices in the anti-clockwise ordering
                Vertex_ptr v_center, v_prior, v_post;
                if (std::get<bool>(e0_pair)) {
                    v_center = std::get<Edge_ptr>(e0_pair)->a;
                    v_prior = std::get<Edge_ptr>(e0_pair)->b; 
                }
                else {
                    v_center = std::get<Edge_ptr>(e0_pair)->b;
                    v_prior = std::get<Edge_ptr>(e0_pair)->a;
                }
                if (std::get<bool>(e1_pair)) {
                    v_post = std::get<Edge_ptr>(e1_pair)->a;
                }
                else {
                    v_post = std::get<Edge_ptr>(e1_pair)->b;
                }

                double dA_dx = 0.5 * (v_post->y - v_prior->y) * c->area_sgn;
                double dA_dy = 0.5 * (v_prior->x - v_post->x) * c->area_sgn;
                // TODO periodicity and fix
                
                v_center->fx -= _area_elasticity * (c->area - _area_preferential)*dA_dx;
                v_center->fy -= _area_elasticity * (c->area - _area_preferential)*dA_dy;
            }
        }
        
        // update vertex positions from forces
        for (auto v : _vertices) {
            v->x += v->fx * _dt;
            v->y += v->fy * _dt;
        }
    }


    /// Monitor model information
    /** \details Here, functions and values can be supplied to the monitor that
     *          are then available to the frontend. The monitor() function is
     *          _only_ called if a certain emit interval has passed; thus, the
     *          performance hit is small.
     */
    void monitor () {
        // Can supply information to the monitor here in two ways:
        // this->_monitor.set_entry("key", value);
        // this->_monitor.set_entry("key", [this](){return 42.;});
    }


    /// Write data
    /** \details This function is called to write out data. It should be called
      *         at the end of the model constructor to write out the initial
      *         state. After that, the configuration determines at which times
      *         data is written.
      *         See \ref Utopia::DataIO::Dataset::write
      */
    void write_data () {
        const auto num_vertices = _vertices.size();
        const auto num_edges = _edges.size();
        const auto num_cells = _cells.size();

        // Get a logger to use here (Note: needs to have been setup beforehand)
        auto log = spdlog::get("data_io");
        log->info("Saving graph with {} vertices and {} edges ...",
                num_vertices,
                num_edges);

        // Store additional metadata in the group attributes
        auto dset_vs = _grp_vertices->open_dataset(std::to_string(this->_time),
                                                   {2, num_vertices});
        dset_vs->add_attribute("num_vertices", num_vertices);

        int running_id = 0;
        dset_vs->write(_vertices.begin(), _vertices.end(), [&](auto v) {
            v->current_id = running_id++;
            return v->x;
        });
        dset_vs->write(_vertices.begin(), _vertices.end(), [&](auto v) {
            return v->y;
        });

        auto dset_es = _grp_edges->open_dataset(std::to_string(this->_time), {2, num_edges});
        dset_es->add_attribute("num_edges", num_edges);
        
        dset_es->write(_edges.begin(), _edges.end(), [&](auto e) {
            return e->a->current_id;
        });
        dset_es->write(_edges.begin(), _edges.end(), [&](auto e) {
            return e->b->current_id;
        });
        
        auto dset_cs = _grp_cells->open_dataset(std::to_string(this->_time), {2, num_cells});
        dset_cs->add_attribute("num_cells", num_cells);

        dset_cs->write(_cells.begin(), _cells.end(), [&](auto c) {
            return c->cell_center()->x;
        });
        dset_cs->write(_cells.begin(), _cells.end(), [&](auto c) {
            return c->s->y;
        });

        double tot_forces_2 = 0;
        for (const auto &v : _vertices) {
            tot_forces_2 += pow(v->fx, 2) + pow(v->fy, 2);
        }
        _dset_forces->write(sqrt(tot_forces_2)/_vertices.size());

    }


    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

};

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPVERTEX_HH
