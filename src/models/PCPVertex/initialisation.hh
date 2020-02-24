#ifndef UTOPIA_MODELS_SAVANNAHETEROGENEOUS_INITIALISATION_HH
#define UTOPIA_MODELS_SAVANNAHETEROGENEOUS_INITIALISATION_HH

namespace Utopia {
namespace Models {
namespace PCPVertex {
    
// https://courses.cs.washington.edu/courses/cse326/00wi/projects/voronoi.html
template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc,
               polarity_proteins>::initialise_voronoi (int num_cells)
{
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
template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc,
               polarity_proteins>::initialise_hexagonal (double size,
                                                        int num_rows,
                                                        int num_columns)
{
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
    double linetension = _linetension(CellType::progenitor,
                                      CellType::progenitor);
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
                    linetension, nullptr, nullptr, 0., 0.));
                // lower right edge
                _edges.push_back(std::make_shared<Edge>(
                    _vertices[2 * c_id + 1],
                    _vertices[2*((q+1)%lim_columns + r*lim_columns)],
                    linetension, nullptr, nullptr, 0., 0.));
                // left edge
                _edges.push_back(std::make_shared<Edge>(
                    _vertices[2 * c_id],
                    _vertices[2*(q + ((r+1)%lim_rows)*lim_columns)],
                    linetension, nullptr, nullptr, 0., 0.));
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
                    linetension, nullptr, nullptr, 0., 0.));
                // lower right edge
                // connects the lower vertex with
                // the lower left vertex of the cell to the right
                _edges.push_back(std::make_shared<Edge>(
                    _vertices[2 * c_id + 1],
                    _vertices[2*((q+1)%lim_columns + r*lim_columns)],
                    linetension, nullptr, nullptr, 0., 0.));
                // left edge connects the lower left vertex with
                // the upper vertex of the cell in the row above
                _edges.push_back(std::make_shared<Edge>(
                    _vertices[2 * c_id + 1],
                    _vertices[2*(q + ((r+1)%lim_rows)*lim_columns) + 1],
                    linetension, nullptr, nullptr, 0., 0.));
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
    double area_preferential = _area_preferential(CellType::progenitor);
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
                                        area_preferential, _contractility));
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
                                        area_preferential, _contractility));
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

    // calculate cell area
    for (auto &c : _cells) {
        c->template cell_area<periodic_bc>();
    }

    this->_log->info("Initialised hexagonal cells.");
}

template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc,
               polarity_proteins>::initialise_polarity_random (
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

template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc,
               polarity_proteins>::prolog ()
{
    // calculate cell area
    for (auto &c : _cells) {
        c->template cell_area<periodic_bc>();
    }

    // initialise the energy terms
    // NOTE since no update is performed, the position of the vertices and
    //      the level of polarity is not changed

    // line tension
    _energy_linetension = 0;
    for (auto &e : _edges) {
        _energy_linetension += this->line_tension(e);
    }
    
    // area elasticity
    _energy_areaelasticity = 0;
    for (auto &c : _cells) {
        _energy_areaelasticity += this->area_elasticity(c);
    }
    
    // area elasticity
    _energy_contractility = 0;
    for (auto &c : _cells) {
        _energy_contractility += this->contractility(c);
    }

    if constexpr (polarity_proteins) {
        _energy_cell_cell_polarity = 0.;
        for (auto &e : _edges) {
            _energy_cell_cell_polarity += this->cell_cell_polarity(e);
        }

        _energy_polarity_exclusion = 0.;
        for (auto &c : _cells) {
            _energy_polarity_exclusion += this->apply_polarity_exclusion(c);
        }

        _energy_lagrange_net_polarisation = 0.;
        for (auto &c : _cells) {
            _energy_lagrange_net_polarisation += 
                    this->lagrange_net_polarisation(c);
            c->lagrange_net_polarisation = 0.; // NOTE undo changes
        }

        _energy_lagrange_const_concentration = 0.;
        for (auto &c : _cells) {
            _energy_lagrange_const_concentration += 
                    this->lagrange_const_concentration(c);
            c->lagrange_const_concentration = 0.; // NOTE undo changes
        }
    }

    return this->__prolog();
}

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif