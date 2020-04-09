#ifndef UTOPIA_MODELS_PCPVERTEX_INITIALISATION_HH
#define UTOPIA_MODELS_PCPVERTEX_INITIALISATION_HH

namespace Utopia::Models::PCPVertex {

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
template<class Model>
void EntitiesManager<Model>::setup_agents_hexagonal_structure (
        const Config& cfg)
{
    if (_space->dim != 2) {
        throw std::invalid_argument("Initialisation of hexagonal arrangement "
            "only defined in 2 dimensional space!");
    }

    double size = get_as<double>("hexagon_size", cfg);
    int num_rows = get_as<int>("lattice_rows", cfg);
    int num_columns = get_as<int>("lattice_columns", cfg);

    SpaceVec cell_shape = SpaceVec({sqrt(3), 2}) * size;

    if (_space->periodic) {
        if (num_rows % 2 != 0) {
            throw std::invalid_argument( "\nERROR with periodic boundary "
                "conditions the hexagonal cell lattice needs pair number "
                "of rows, but received impair number. Requested "
                "number of rows was " + std::to_string(num_columns) + "!");
        }
        
        _space->set_domain_size(SpaceVec({double(num_columns),
                                           0.75 * num_rows}) % cell_shape);
    }
    else {
        _space->set_domain_size(SpaceVec({num_columns + 1.,
                                           0.75 * (num_rows + 1.)}) % 
                                           cell_shape);
    }

    // Add vertices
    // NOTE one more row and column of vertices initialized in
    //      non-periodic bc this will be undone before initializing edges  
    //      and cells
    // NOTE some of these vertices have to be removed eventually
    int lim_rows = num_rows;
    int lim_columns = num_columns;
    if (not _space->periodic) {
        lim_columns += 1;
        lim_rows += 1;
    }
    for (int r = 0; r < lim_rows; r += 1) {
        // pair rows
        if (r % 2 == 0) {
        for (int q = 0; q < lim_columns; ++q) {
            this->add_vertex(SpaceVec({double(q), (.75*r + .25)}) % cell_shape);
            this->add_vertex(SpaceVec({q + 0.5, 0.75 * r}) % cell_shape);
        }
        }
        // impair rows
        else {
        for (int q = 0; q < lim_columns; ++q) {
            this->add_vertex(SpaceVec({double(q), r * 0.75}) % cell_shape);
            this->add_vertex(SpaceVec({q + 0.5, r * 0.75 + 0.25}) % cell_shape);
        }
        }
    }

    // Add edges
    const auto vertices = this->vertices();
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
                this->add_edge(vertices[2 * c_id], vertices[2*c_id + 1]);
                // lower right edge
                this->add_edge(vertices[2 * c_id + 1],
                        vertices[2*((q+1)%lim_columns + r*lim_columns)]);
                // left edge
                this->add_edge(vertices[2 * c_id],
                        vertices[2*(q + ((r+1)%lim_rows)*lim_columns)]);
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
                this->add_edge(vertices[2*c_id], vertices[2*c_id + 1]);
                // lower right edge
                // connects the lower vertex with
                // the lower left vertex of the cell to the right
                this->add_edge(vertices[2 * c_id + 1],
                        vertices[2*((q+1)%lim_columns + r*lim_columns)]);
                // left edge connects the lower left vertex with
                // the upper vertex of the cell in the row above
                this->add_edge(vertices[2 * c_id + 1],
                        vertices[2*(q + ((r+1)%lim_rows)*lim_columns) + 1]);
            }
        }
    }

    const auto edges = this->edges();
    // ** add cells
    // handle last row separately
    for (int r = 0; r < num_rows; r++) {
        // pair rows
        if (r % 2 == 0) {
            for (int q = 0; q < num_columns; q++) {
                int c_id = q + r * lim_columns;
                this->add_cell(
                    SpaceVec({q + 0.5, 0.75 * r + 0.5}) % cell_shape,
                    { edges[3*c_id], // lower left
                      edges[3*c_id + 1], // lower right
                      edges[3*c_id + 2], // left
                      edges[3*((q+1)%lim_columns + r*lim_columns) + 2], // right
                      edges[3*(c_id + lim_columns)], // upper left
                      edges[3*(c_id + lim_columns) + 1] // upper right
                    });
            }
        }
        else { // impare rows
            for (int q = 0; q < num_columns; q++) {
                int c_id = q + r * lim_columns;
                this->add_cell(
                    SpaceVec({q + 1., 0.75 * r + 0.5}) % cell_shape,
                    { edges[3*c_id + 1], // lower left
                      edges[3*c_id + 2], // lower right
                      edges[3*((q+1)%lim_columns + r*lim_columns)], // left
                      edges[3*((q+1)%lim_columns + r*lim_columns) + 2], // right
                      // upper left
                      edges[3*(q + (r+1)%lim_rows * lim_columns) + 1],
                      // upper right
                      edges[3*((q+1)%lim_columns + (r+1)%lim_rows * lim_columns)]
                    });
            }
        }
    }

    // remove not needed objects
    // in non periodic bc more objects are initialised that needed
    if (not _space->periodic) {
        
        // vertices
        AgentContainer<Vertex> vertices_remove;
        vertices_remove.push_back(vertices[2*(lim_columns - 1) + 1]);
        if (num_rows % 2 == 1) {
            vertices_remove.push_back(vertices.back());
        }
        else {
            vertices_remove.push_back(vertices[2*(lim_rows-1)*lim_columns]);
        }

        // edges
        AgentContainer<Edge> edges_remove;
        for (int r = 0; r < lim_rows; r++) {
            if (r % 2 == 0) {
                if (r == 0) {
                    edges_remove.push_back(edges[3*(lim_columns - 1)]);
                }
                
                edges_remove.push_back(edges[3*((r+1)*lim_columns - 1) + 1]);

                if (r == lim_rows - 1) {
                    edges_remove.push_back(edges[3*(r*lim_columns)]);
                    for (int q = 0; q < lim_columns; ++q) {
                        edges_remove.push_back(
                            edges[3*(q + r*lim_columns) + 2]);
                    }
                }
            }
            else {
                edges_remove.push_back(edges[3*((r+1)*lim_columns - 1) + 1]);

                if (r == lim_rows - 1) {
                    edges_remove.push_back(edges[3*(lim_rows*lim_columns - 1)]);

                    for (int q = 0; q < lim_columns; ++q) {
                        edges_remove.push_back(
                            edges[3*(q + r*lim_columns) + 2]);
                    }
                }
            }
        }

        for (const auto& v : vertices_remove) {
            _vertex_manager.remove_agent(v);
        }
        for (const auto& e : edges_remove) {
            _edge_manager.remove_agent(e);
        }
    }

    for (const auto& e : this->edges()) {
        position_of(e->custom_links().a);
        position_of(e->custom_links().b);
    }

    this->_log->info("Initialised hexagonal cells.");
}

} // namespace Utopia::Models::PCPVertex
#endif