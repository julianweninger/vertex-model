#ifndef UTOPIA_MODELS_PCPVERTEX_TRANSITIONS_HH
#define UTOPIA_MODELS_PCPVERTEX_TRANSITIONS_HH

#ifndef PI
#define PI 3.14159265
#endif

namespace Utopia::Models::PCPVertex {

/// Perform a cell division on specific cell
/** Divides a specific cell into two identical cells with properties derived
 *  from the common parent cell. 
 *  The division is performed at a given angle through the parent cell's
 *  center. This defines the axis of division that will form a new edge 
 *  between the two new cells.
 * 
 *  The new edge has properties as given for initialisation.
 * 
 *  \param cell         the cell that is to be divided
 *  \param division_angle   angle (in rad) at which the cell is divided
 *  \param linetension  The linetension of the dividing edge
 *  \param edge_contractility   The contractility of the dividing egde
 * 
 *  \note The remaining parameters of objects are inherited from existing 
 *        objects.
 */
template<class Model>
void CustomAgentManager<Model>::divide_cell(const std::shared_ptr<Cell> cell,
        double division_angle, double linetension, double edge_contractility)
{
    // The cell to be divided
    cell->state.remove = true;
    _cell_manager.remove_agent(cell);

    const auto cell_center = this->barycenter_of(cell);

    // generate the axis of division
    auto domain = this->_space->get_domain_size();
    const SpaceVec axis = SpaceVec({cos(division_angle), sin(division_angle)}) / 
                          domain;

    // determine the new vertices from this axis
    // These are the intersections of the division axis with edges of cell
    // NOTE there must be exactly 2 intersections for nicely shaped cells
    AgentContainer<Vertex> new_vertices;
    for (const auto& [e, flip] : cell->custom_links().edges) {
        // calculate the intersection
        auto a = this->position_of(e->custom_links().a);
        auto b = this->position_of(e->custom_links().b);
        const auto [inter, success] = this->_space->intersection(
                a, this->_space->displacement(a, b),
                cell_center, axis,
                true, false);
        
        // if no intersection end here
        if (not success) { continue; }

        // create a new vertex at the intersection site
        auto new_v = this->add_vertex(inter);
        const auto [adj_cell_a, adj_cell_b] = this->adjoints_of(e);
        _vertices_adjoint_cells[new_v->id()] = {adj_cell_a, adj_cell_b};
        _vertices_adjoint_edges[new_v->id()] = {e};
        // NOTE the edge link is used later and removed thereafter
        
        // keep track of this vertex
        new_vertices.push_back(new_v);

        // the edge itself will be divided and removed
        e->state.remove = true;
    }
    // check for valid setting
    if (new_vertices.size() != 2) {
        this->_log->warn("Cannot perform cell division on cell at ({}, {}), "
            "with a division angle of {}.", cell_center[0], cell_center[1],
            division_angle / 2 / PI * 360);
        this->_log->warn("division axis is {}, {} -> {}, {}",
            cell_center[0], cell_center[1],
            (axis + cell_center)[0], (axis + cell_center)[1]);
        this->_log->warn("Edges of cell are:");
        for (const auto& [e, flip] : cell->custom_links().edges) {
            auto a = this->position_of(e->custom_links().a);
            auto b = this->position_of(e->custom_links().b);
            if (flip) { std::swap(a, b); }
            this->_log->warn("   {}, {} -> {}, {}",a[0], a[1], b[0], b[1]);
        }
        this->_log->warn("Intersections are:");
        for (const auto& v : new_vertices) {
            auto pos = this->position_of(v);
            this->_log->warn("   {}, {}", pos[0], pos[1]);
        }
        throw std::runtime_error("During cell division, expected 2 new "
            "vertices, but got " + std::to_string(new_vertices.size()) + "!");
    }

    // the 2 edges that are divided by the new edge
    const auto edge_a = adjoint_edges_of(new_vertices[0]).back();
    const auto edge_b = adjoint_edges_of(new_vertices[1]).back();

    // remove them, they will be replaced by 2 new edges each
    _edge_manager.erase_agent_if([](const auto& edge){
                                    return edge->state.remove; });
    _vertices_adjoint_edges[new_vertices[0]->id()].clear();
    _vertices_adjoint_edges[new_vertices[1]->id()].clear();

    // create the new edge dividing the cell
    DataIO::Config edge_cfg;
    edge_cfg["linetension"] = linetension;
    edge_cfg["contractility"] = edge_contractility;
    const auto new_edge = add_edge(new_vertices[0], new_vertices[1], edge_cfg);

    // this is how to divide an edge at a pivot vertex
    /* \param e     The Edge to divide
        * \param flip  The direction of e
        * \param pivot The Vertex where to divide e
        * 
        * \return  the first half of e from a to pivot, the second half of e 
        *          from pivot to b
        */
    auto divide_edge = [this, cell](std::shared_ptr<Edge> edge, bool flip,
                              std::shared_ptr<Vertex> pivot)
    {
        // The start and end of e, considering the direction within the 
        // iteration
        auto a = edge->custom_links().a;
        auto b = edge->custom_links().b;
        if (flip) { std::swap(a, b); }

        for (auto v : {a, b}) {
            auto adj_edges = _vertices_adjoint_edges.at(v->id());
            adj_edges.erase(std::remove_if(
                    adj_edges.begin(), adj_edges.end(),
                    [](const auto& e){ return e->state.remove; }),
                adj_edges.end());
        }

        DataIO::Config edge_cfg;
        edge_cfg["linetension"] = edge->state.linetension;
        edge_cfg["contractility"] = edge->state.contractility;

        // Divide edge at pivot intow two new edges
        auto new_edge_0 = this->add_edge(a, pivot, edge_cfg);
        auto new_edge_1 = this->add_edge(pivot, b, edge_cfg);

        _edges_adjoint_cells[new_edge_0->id()] = adjoints_of(edge);
        _edges_adjoint_cells[new_edge_1->id()] = adjoints_of(edge);

        // update the crosslinks in the adjoint cell
        const auto [adj_cell_a, adj_cell_b] = adjoints_of(edge);
        for (const auto& c : {adj_cell_a, adj_cell_b}) {
            if (not c) { continue; }
            if (c == cell) { continue; }
            else {
                auto& adj_edges = c->custom_links().edges;
                auto it = adj_edges.erase(std::remove_if(
                        adj_edges.begin(), adj_edges.end(),
                        [](const auto& e){ 
                            return std::get<0>(e)->state.remove; }),
                    adj_edges.end());
                // edge is flipped in neighboring cell
                // add b->pivot (edge_1, true) then pivot->a (edge_0, true)
                adj_edges.insert(it, std::make_pair(new_edge_1, true));
                adj_edges.insert(it, std::make_pair(new_edge_0, true));
            }
            auto& vertices = c->custom_links().vertices;
            vertices.insert(std::find(vertices.begin(), vertices.end(), a),
                            pivot);
        }

        return std::make_pair(new_edge_0, new_edge_1);
    };
            
    // separate the edges to form 2 cells
    AgentContainer<Edge> new_edges_cell_0, new_edges_cell_1;
    
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
    const auto& edges = cell->custom_links().edges;
    // 1. start iteration
    for (it_edges = 0; it_edges < edges.size(); it_edges++) {
        const auto [e, flip] = edges[it_edges];
        // check whether edge_a or _b reached
        if (e == edge_a or e == edge_b) {
            // the intersection with the division axis
            std::shared_ptr<Vertex> pivot;
            if (e == edge_a) { pivot = new_vertices[0]; }
            else { pivot = new_vertices[1]; }

            std::shared_ptr<Edge> new_edge_0, new_edge_1;
            std::tie(new_edge_0, new_edge_1) = divide_edge(e, flip, pivot);

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
    for (it_edges += 1; it_edges < edges.size(); it_edges++) {
        const auto [e, flip] = edges[it_edges];
        // check whether edge_a or _b reached
        if (e == edge_a or e == edge_b) {
            // the intersection with the division axis
            std::shared_ptr<Vertex> pivot;
            if (e == edge_a) { pivot = new_vertices[0]; }
            else { pivot = new_vertices[1]; }

            std::shared_ptr<Edge> new_edge_0, new_edge_1;
            std::tie(new_edge_0, new_edge_1) = divide_edge(e, flip, pivot);

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
    for (it_edges += 1; it_edges < edges.size(); it_edges++) {
        auto [e, flip] = edges[it_edges];
        new_edges_cell_0.push_back(e);
    }

    // create 2 new cells
    DataIO::Config cell_cfg;
    cell_cfg["area_preferential"] = cell->state.area_preferential;
    cell_cfg["contractility"] = cell->state.contractility;
    cell_cfg["protein_concentration"] = cell->state.protein_concentration;
    
    using CellType = typename Cell::State::CellType;
    if (cell->state.type == CellType::progenitor) {
        cell_cfg["cell_type"] = "progenitor";
    }
    else if (cell->state.type == CellType::hair) {
        cell_cfg["cell_type"] = "hair";
    }
    else if (cell->state.type == CellType::support) {
        cell_cfg["cell_type"] = "support";
    }
    else { 
        cell_cfg["cell_type"] = "not implemented within division";
    }
    const auto new_cell_0 = this->add_cell(cell_center, new_edges_cell_0,
                                           cell_cfg);
    const auto new_cell_1 = this->add_cell(cell_center, new_edges_cell_1,
                                           cell_cfg);
    
    // remove expired crosslinks
    for (const auto& new_c : {new_cell_0, new_cell_1}) {
        for (const auto& v : new_c->custom_links().vertices) {
            auto adj_cells = _vertices_adjoint_cells.at(v->id());

            adj_cells.erase(std::remove_if(adj_cells.begin(), adj_cells.end(),
                                           [cell](const auto& c){
                                                return c == cell; }),
                            adj_cells.end());
        }
        for (const auto& [e, flip] : new_c->custom_links().edges) {
            auto [adj_cell_a, adj_cell_b] = adjoints_of(e);
            if (adj_cell_a == cell) {
                adj_cell_a = new_c;
            }
            else if (adj_cell_b == cell) {
                adj_cell_b = new_c;
            }
            _edges_adjoint_cells[e->id()] = std::make_pair(adj_cell_a,
                                                           adj_cell_b);
        }
    }
}

} // namespace Utopia::Models::PCPVertex
#endif