#ifndef UTOPIA_MODELS_PCPVERTEX_INITIALISATION_HH
#define UTOPIA_MODELS_PCPVERTEX_INITIALISATION_HH

#include "utopia/core/zip.hh"

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
    std::size_t num_rows = get_as<int>("lattice_rows", cfg);
    std::size_t num_columns = get_as<int>("lattice_columns", cfg);

    enum HexShape {
        PointyTop,
        FlatTop
    };
    std::function<HexShape(const Config&)> setup_hexshape =
        [this](const Config& cfg) {
            auto shape = get_as<std::string>("hexagon_shape", cfg);

            if (shape == "pointy" or shape == "pointy_top") {
                return HexShape::PointyTop;
            }
            else if (shape == "flat" or shape == "flat_top") {
                return HexShape::FlatTop;
            }
            else {
                this->_log->error("Failed to set up hexagonal cell structure "
                    "with shape {}", shape);
                throw std::invalid_argument(fmt::format(
                    "Invalid `hexagon_shape` {} when setting up cells in "
                    "Vertex model, chose one of the following: "
                    "pointy, "
                    "pointy_top, "
                    "flat, "
                    "flat_top.",
                    shape
                ));
            }
        };
    HexShape shape = setup_hexshape(cfg);
    
    // The following is for a PointyTop layout
    // The FlatTop layout is generated by a 90 deg rotation of the domain
    if (shape == HexShape::FlatTop) {
        std::swap(num_rows, num_columns);
    }

    SpaceVec cell_shape = SpaceVec({sqrt(3), 2}) * size;

    if (_space->periodic) {
        if (num_rows % 2 != 0) {
            throw std::invalid_argument(fmt::format(
                "ERROR with periodic boundary conditions the {} hexagonal "
                "cell lattice needs pair number of {}, but received impair "
                "number ({})!",
                get_as<std::string>("hexagon_shape", cfg),
                shape == HexShape::PointyTop ? "rows" : "columns",
                num_rows
            ));
        }
        
        _space->set_domain_size(SpaceVec({double(num_columns),
                                           0.75 * num_rows}) % cell_shape);
    }

    // Add vertices
    // NOTE one more row and column of vertices initialized in
    //      non-periodic bc this will be undone before initializing edges  
    //      and cells
    // NOTE some of these vertices have to be removed eventually
    std::size_t lim_rows = num_rows;
    std::size_t lim_columns = num_columns;
    if (not _space->periodic) {
        lim_columns += 1;
        lim_rows += 1;
    }
    for (std::size_t r = 0; r < lim_rows; r += 1) {
        // pair rows
        if (r % 2 == 0) {
        for (std::size_t q = 0; q < lim_columns; ++q) {
            SpaceVec position = SpaceVec({double(q),
                                         (.75*r + .25)}) % cell_shape;
            this->add_vertex(position);
            
            position = SpaceVec({q + 0.5, 0.75 * r}) % cell_shape;
            this->add_vertex(position);
        }
        }
        // impair rows
        else {
        for (std::size_t q = 0; q < lim_columns; ++q) {
            SpaceVec position = SpaceVec({double(q), r * 0.75}) % cell_shape;
            this->add_vertex(position);
            
            position = SpaceVec({q + 0.5, r * 0.75 + 0.25}) % cell_shape;
            this->add_vertex(position);
        }
        }
    }

    // Add edges
    const auto vertices = this->vertices();
    for (std::size_t r = 0; r < lim_rows; r++) {
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
            std::size_t c_id; // id of the resp cell
            for (std::size_t q = 0; q < lim_columns; ++q) {
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
            std::size_t c_id;
            for (std::size_t q = 0; q < lim_columns; ++q) {
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
    for (std::size_t r = 0; r < num_rows; r++) {
        // pair rows
        if (r % 2 == 0) {
            for (std::size_t q = 0; q < num_columns; q++) {
                std::size_t c_id = q + r * lim_columns;
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
            for (std::size_t q = 0; q < num_columns; q++) {
                std::size_t c_id = q + r * lim_columns;
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
        for (std::size_t r = 0; r < lim_rows; r++) {
            if (r % 2 == 0) {
                if (r == 0) {
                    edges_remove.push_back(edges[3*(lim_columns - 1)]);
                }
                
                edges_remove.push_back(edges[3*((r+1)*lim_columns - 1) + 1]);

                if (r == lim_rows - 1) {
                    edges_remove.push_back(edges[3*(r*lim_columns)]);
                    for (std::size_t q = 0; q < lim_columns; ++q) {
                        edges_remove.push_back(
                            edges[3*(q + r*lim_columns) + 2]);
                    }
                }
            }
            else {
                edges_remove.push_back(edges[3*((r+1)*lim_columns - 1) + 1]);

                if (r == lim_rows - 1) {
                    edges_remove.push_back(edges[3*(lim_rows*lim_columns - 1)]);

                    for (std::size_t q = 0; q < lim_columns; ++q) {
                        edges_remove.push_back(
                            edges[3*(q + r*lim_columns) + 2]);
                    }
                }
            }
        }

        for (const auto& v : vertices_remove) {
            v->state.remove = true;
            this->remove_vertex(v);
        }
        for (const auto& e : edges_remove) {
            e->state.remove = true;

            // remove the edge from weak links
            for (const auto& v : {e->custom_links().a, e->custom_links().b}) {
                if (not v->state.remove) {
                    auto adjoints = _vertices_adjoint_edges.at(v->id());
                    adjoints.erase(std::remove(adjoints.begin(), adjoints.end(),
                                               e),
                                   adjoints.end());
                    _vertices_adjoint_edges.erase(v->id());
                    _vertices_adjoint_edges[v->id()] = adjoints;
                }
            }
            
            this->remove_edge(e);
        }

        SpaceVec origin = this->barycenter_of(this->get_boundary_edges());

        for (const auto& v : this->vertices()) {
            this->move_by(v, -1. * origin);
        }
    }

    this->_log->info("Initialised hexagonal cells.");

    auto structure_str = get_as<std::string>("HC_structure", cfg, "");
    if (structure_str.length() > 0) {
        this->_log->info("Initialising HC in '{}' structure ...",
                         structure_str);

        if (structure_str == "ratio_1_to_2") {
            if (_space->periodic and num_columns % 3 != 0) {
                throw std::invalid_argument(fmt::format(
                    "Failed to set up HC structure '{}' on a hexagonal lattice "
                    "with {} columns. Columns needs to be a multiple of 3!",
                    structure_str, num_columns
                ));
            }
            for (std::size_t row = 0; row < num_rows; row++) {
                for (std::size_t col = 0; col < num_columns; col++) {
                    auto id = row * num_columns + col;
                    if (   (row % 2 == 0 and col % 3 == 0)
                        or (row % 2 == 1 and col % 3 == 1))
                    {
                        this->cells()[id]->state.type = CellType::hair;
                    }
                    else {
                        this->cells()[id]->state.type = CellType::support;
                    }
                }
            }
        }
        else if (structure_str == "ratio_1_to_3") {
            for (std::size_t row = 0; row < num_rows; row++) {
                for (std::size_t col = 0; col < num_columns; col++) {
                    auto id = row * num_columns + col;
                    if (row % 2 == 0 and col % 2 == 0) {
                        this->cells()[id]->state.type = CellType::hair;
                    }
                    else {
                        this->cells()[id]->state.type = CellType::support;
                    }
                }
            }
        }
        else if (structure_str == "ratio_1_to_4") {
            if (_space->periodic and num_columns % 5 != 0) {
                throw std::invalid_argument(fmt::format(
                    "Failed to set up HC structure '{}' on a hexagonal lattice "
                    "with {} columns. Columns needs to be a multiple of 5!",
                    structure_str, num_columns
                ));
            }
            for (std::size_t row = 0; row < num_rows; row++) {
                for (std::size_t col = 0; col < num_columns; col++) {
                    auto id = row * num_columns + col;
                    if (   (row % 2 == 0 and col % 5 == 0)
                        or (row % 2 == 1 and col % 5 == 2)) {
                        this->cells()[id]->state.type = CellType::hair;
                    }
                    else {
                        this->cells()[id]->state.type = CellType::support;
                    }
                }
            }
        }
        else if (structure_str == "ratio_1_to_5") {
            if (    _space->periodic
                and (num_columns % 3 != 0 or num_rows % 4 != 0)) {
                throw std::invalid_argument(fmt::format(
                    "Failed to set up HC structure '{}' on a hexagonal lattice "
                    "with {} columns and {} rows. Columns and rows need to be "
                    "a multiple of 3 and 4, respectively!",
                    structure_str, num_columns, num_rows
                ));
            }
            for (std::size_t row = 0; row < num_rows; row++) {
                for (std::size_t col = 0; col < num_columns; col++) {
                    auto id = row * num_columns + col;
                    if (row % 2 == 0 and col % 3 == (row % 6) / 2)
                    {
                        this->cells()[id]->state.type = CellType::hair;
                    }
                    else {
                        this->cells()[id]->state.type = CellType::support;
                    }
                }
            }
        }
        else if (structure_str == "ratio_1_to_6") {
            if (    _space->periodic
                and (num_columns % 14 != 0 or num_rows % 14 != 0)) {
                throw std::invalid_argument(fmt::format(
                    "Failed to set up HC structure '{}' on a hexagonal lattice "
                    "with {} columns and {} rows. Columns and rows need to be "
                    "both a multiple of 14!",
                    structure_str, num_columns, num_rows
                ));
            }
            for (int m = -num_rows; m < int(num_rows); m++) {
                for (int n = -num_columns; n < int(num_columns); n++) {
                    int row = 2 * m + n;
                    int col = 4 * n + m + row / 2;

                    if (    row >= 0 and row < int(num_rows)
                        and col >= 0 and col < int(num_columns))
                    {
                        auto id = row * num_columns + col;
                        this->cells()[id]->state.type = CellType::hair;
                    }
                }
            }
            for (const auto& cell : cells()) {
                if (cell->state.type == CellType::progenitor) {
                    cell->state.type = CellType::support;
                }
            }
        }
        else {
            throw std::invalid_argument(fmt::format(
                "Invalid HC structure '{}' in initialisation of hexagonal "
                "lattice! Available structures are: "
                    "ratio_1_to_2, "
                    "ratio_1_to_3, "
                    "ratio_1_to_4, "
                    "ratio_1_to_5, "
                    "ratio_1_to_6.",
                structure_str));
        }
    }

    // Perform the rotation
    if (shape == HexShape::FlatTop) {
        SpaceVec domain = _space->get_domain_size();
        SpaceVec origin = domain / 2.;

        std::vector<SpaceVec> displacements;
        for (const auto& vertex : this->vertices()) {
            SpaceVec displ = _space->displacement(origin, position_of(vertex));
            displacements.push_back(displ);
        }

        _space->set_domain_size(SpaceVec({domain[1], domain[0]}));
        for (const auto& [vertex, displ] :
                    Itertools::zip(this->vertices(), displacements))
        {
            SpaceVec new_pos = (  _space->get_domain_size() / 2.
                                + SpaceVec({displ[1], -displ[0]}));
            move_to(vertex, new_pos);
        }
    }
}


/// Create vertices, edges, and cells in a side view columnar arrangement
template<class Model>
void EntitiesManager<Model>::setup_agents_column (const Config& cfg)
{
    if (_space->dim != 2) {
        throw std::invalid_argument("Initialisation of columnar cell "
            "arrangement only defined in 2 dimensional space!");
    }

    this->_log->debug("Setting up cells in columnar structure ...");

    auto num_cells = get_as<std::size_t>("num_cells", cfg);

    // determine the cell shape
    auto size = get_as<double>("size", cfg);
    auto height = get_as<double>("height", cfg);
    double width = size / height;
    double _offset = get_as<double>("offset", cfg);

    SpaceVec cell_shape({width, height});

    if (_space->periodic) {        
        _space->set_domain_size(  SpaceVec({double(num_cells), 2 * _offset + 1})
                                % cell_shape);
    }

    // Add vertices
    // NOTE one more column of vertices initialized in
    //      non-periodic bc this will be undone before initializing edges  
    //      and cells
    // NOTE some of these vertices have to be removed eventually
    std::size_t num_columns = num_cells;
    std::size_t lim_columns = num_cells;
    if (not _space->periodic) {
        lim_columns += 1;
    }
    
    SpaceVec offset({0., 0.});
    if (_space->periodic) {
        offset = SpaceVec({0., _offset}) % cell_shape;
    }
    for (std::size_t q = 0; q < lim_columns; ++q) {
        SpaceVec position = SpaceVec({double(q), 0.}) % cell_shape;
        this->add_vertex(position + offset);
        
        position = SpaceVec({double(q), 1.}) % cell_shape;
        this->add_vertex(position + offset);
    }
    const auto& vertices = this->vertices();

    // Add edges
    for (std::size_t q = 0; q < num_columns; ++q) {
        // left edge
        this->add_edge(vertices[2 * q], vertices[2 * q + 1]);

        // lower edge
        this->add_edge(vertices[2 * q], vertices[2 * ((q + 1) % lim_columns)]);
        
        // upper edge
        this->add_edge(vertices[2 * q + 1],
                       vertices[2 * ((q + 1) % lim_columns) + 1]);
    }
    if (not _space->periodic) {
        // right edge
        this->add_edge(vertices[2 * num_columns],
                       vertices[2 * num_columns + 1]);
    }
    const auto& edges = this->edges();

    // Add cells
    for (std::size_t q = 0; q < num_columns; q++) {
        this->add_cell(
            SpaceVec({q + 0.5, 0.5}) % cell_shape,
            {
              edges[3 * q],       // left
              edges[3 * q + 1],   // lower
              edges[3 * ((q + 1) % lim_columns)], // right
              edges[3 * q + 2]    // upper
            });
    }

    this->_log->info("Initialised {} cells in columnar structure.",
                     this->cells().size());
}

} // namespace Utopia::Models::PCPVertex
#endif