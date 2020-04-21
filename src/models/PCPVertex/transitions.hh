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
 *  
 *  \warning Does not conserve the ordering of vertices
 */
template<class Model>
void EntitiesManager<Model>::divide_cell(const std::shared_ptr<Cell> cell,
        double division_angle, double linetension, double edge_contractility)
{
    if (not _space->periodic) {
        throw std::runtime_error("Cell division not implemented in "
            "non-periodic boundary condition!");
    }

    // The cell to be divided
    cell->state.remove = true;
    remove_cell(cell);

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
        this->_log->error("Cannot perform cell division on cell at ({}, {}), "
            "with a division angle of {}.", cell_center[0], cell_center[1],
            division_angle / 2 / PI * 360);
        this->_log->error("division axis is {}, {} -> {}, {}",
            cell_center[0], cell_center[1],
            (axis + cell_center)[0], (axis + cell_center)[1]);
        this->_log->error("Edges of cell are:");
        for (const auto& [e, flip] : cell->custom_links().edges) {
            auto a = this->position_of(e->custom_links().a);
            auto b = this->position_of(e->custom_links().b);
            if (flip) { std::swap(a, b); }
            this->_log->error("   {}, {} -> {}, {}",a[0], a[1], b[0], b[1]);
        }
        this->_log->error("Intersections are:");
        for (const auto& v : new_vertices) {
            auto pos = this->position_of(v);
            this->_log->error("   {}, {}", pos[0], pos[1]);
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
            auto& adj_edges = _vertices_adjoint_edges[v->id()];
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
            if (c == cell) { continue; }

            auto& adj_edges = c->custom_links().edges;
            auto it = adj_edges.erase(std::find_if(adj_edges.begin(),
                            adj_edges.end(), [edge](auto& e_pair) {
                                return std::get<0>(e_pair) == edge; }));
            // edge is flipped in neighboring cell
            // add b->pivot (edge_1, true) then pivot->a (edge_0, true)
            it = adj_edges.insert(it, std::make_pair(new_edge_1, true));
            adj_edges.insert(it+1, std::make_pair(new_edge_0, true));

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
    unsigned int it_edges;
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
            auto& adj_cells = _vertices_adjoint_cells[v->id()];

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

    return;
} // divide cell


/// Removes an edge in a T1 neighborhood exchange
/** \details Adds a new edge separating the adjoint cells
 * 
 *  \param linetension  The linetension of the new edge
 *  \param contractility  The contractility of the new edge
 *  \param get_energy   Calculate the energy for container of edges and cells
 *  \param T1_barrier   The height of the energy barrier
 * 
 *  \returns whether edge was removed.
 *           The edge is not removed if one of the adjoint cells is triangular.
 *           Remove it rather in a T2 transition.
 *           The edge is also not removed if the transition increases energy.
 */
template<class Model>
bool EntitiesManager<Model>::remove_edge_T1 (const std::shared_ptr<Edge> edge,
        double linetension, double contractility,
        std::function<double(const AgentContainer<Edge>&,
                             const AgentContainer<Cell>&)> get_energy,
        double T1_threshold, double T1_barrier, double random_number)
{
    if (not _space->periodic) {
        throw std::runtime_error("T1 transition not implemented in "
            "non-periodic boundary condition!");
    }

    auto vertex_a = edge->custom_links().a;
    auto vertex_b = edge->custom_links().b;

    // create copies of the objects - T1 might be aborted
    const Edge edge_copy = *edge;
    const auto adjoint_edges_vertex_a = adjoint_edges_of(vertex_a);
    const auto adjoint_edges_vertex_b = adjoint_edges_of(vertex_b);

    // tag objects to be removed
    edge->state.remove = true;
    vertex_a->state.remove = true;
    vertex_b->state.remove = true;

    auto [adj_cell_a, adj_cell_b] = adjoints_of(edge);
    // NOTE a, b arbitrary
    
    for (auto c : {adj_cell_a, adj_cell_b}) {
        if (c->custom_links().edges.size() <= 3) {
            this->_log->warn("The triangular cell has area {}. Adapt the "
                "threshold area for T2 transitions, such that this cell "
                "is removed, rather than performing a T1 transition!",
                area_of(c));
            return false;
            // TODO check algorithm at this point
        }
    }

    // the two cells that become neighbours in this T1 transition
    // NOTE c is the cell adjoint to vertex edge->a (arbitrary)
    //      d is the cell adjoint to vertex edge->b
    std::shared_ptr<Cell> adj_cell_c, adj_cell_d;
    for (const auto& c : adjoint_cells_of(vertex_a)) {
        if (c != adj_cell_a and c != adj_cell_b) {
            adj_cell_c = c;
        }
    }
    for (const auto& c : adjoint_cells_of(vertex_b)) {
        if (c != adj_cell_a and c != adj_cell_b) {
            adj_cell_d = c;
        }
    }

    // The involved edges
    // (a, b) and (c, d) currently share a vertex
    // (a, c) and (b, c) currently share an adjoint cell a, resp. b
    // they will share a vertex a (resp. b) after transition
    std::shared_ptr<Edge> adj_edge_a, adj_edge_b, adj_edge_c, adj_edge_d;
    for (const auto& e : adjoint_edges_of(vertex_a)) {
        if (e != edge) {
            const auto& [e_adj_cell_a, e_adj_cell_b] =  adjoints_of(e);
            if (e_adj_cell_a == adj_cell_a or e_adj_cell_b == adj_cell_a)
            {
                adj_edge_a = e;
            }
            else {
                adj_edge_b = e;
            }
        }
    }
    for (const auto& e : adjoint_edges_of(vertex_b)) {
        if (e != edge) {
            const auto& [e_adj_cell_a, e_adj_cell_b] =  adjoints_of(e);
            if (e_adj_cell_a == adj_cell_a or e_adj_cell_b == adj_cell_a)
            {
                adj_edge_c = e;
            }
            else {
                adj_edge_d = e;
            }
        }
    }

    // make copies of the status quo
    const auto adj_edge_a_copy = *adj_edge_a;
    const auto adj_edge_b_copy = *adj_edge_b;
    const auto adj_edge_c_copy = *adj_edge_c;
    const auto adj_edge_d_copy = *adj_edge_d;

    const auto adj_cell_a_copy = *adj_cell_a;
    const auto adj_cell_b_copy = *adj_cell_b;
    const auto adj_cell_c_copy = *adj_cell_c;
    const auto adj_cell_d_copy = *adj_cell_d;

    double current_energy = get_energy({edge, adj_edge_a, adj_edge_b,
                                              adj_edge_c, adj_edge_d},
                                             {adj_cell_a, adj_cell_b,
                                              adj_cell_c, adj_cell_d});

    // create two new vertices that create an edge of threshold length 
    // pointing from cell a to b
    auto displ = _space->displacement(barycenter_of(adj_cell_a),
                                      barycenter_of(adj_cell_b));
    displ = displ / arma::norm(displ) * 2 * T1_threshold /
            _space->get_domain_size();
    
    auto center = position_of(vertex_a) + displacement(vertex_a, vertex_b) / 2;
    auto new_v_a = add_vertex(SpaceVec(center - displ / 2.));
    auto new_v_b = add_vertex(SpaceVec(center + displ / 2.));

    _vertices_adjoint_edges[new_v_a->id()] = {adj_edge_a, adj_edge_c};
    _vertices_adjoint_edges[new_v_b->id()] = {adj_edge_b, adj_edge_d};

    _vertices_adjoint_cells[new_v_a->id()] = {adj_cell_a, adj_cell_c,
                                              adj_cell_d};
    _vertices_adjoint_cells[new_v_b->id()] = {adj_cell_b, adj_cell_c,
                                              adj_cell_d};

    // create a new edge
    DataIO::Config edge_cfg;
    edge_cfg["linetension"] = linetension;
    edge_cfg["contractility"] = contractility;
    auto new_edge = add_edge(new_v_a, new_v_b, edge_cfg);
    _edges_adjoint_cells[new_edge->id()] = std::make_pair(adj_cell_c,
                                                          adj_cell_d);

    // remove objects
    for (auto &c : {adj_cell_a, adj_cell_b, adj_cell_c, adj_cell_d}) {
        auto& vertices = c->custom_links().vertices;
        vertices.erase(std::remove_if(
                vertices.begin(), vertices.end(), 
                [](const auto& v) { return v->state.remove; }),
            vertices.end());
    }
    // NOTE edge_it will be removed at the very end
    // NOTE a and b will be replaced within edges

    // remove edge from adj_cells a and b
    // NOTE the neighbouring edges now have a common vertex, hence order of
    //      edges is maintained
    for (auto &c : {adj_cell_a, adj_cell_b}) {
        auto& edges = c->custom_links().edges;
        edges.erase(std::remove_if(
                edges.begin(), edges.end(), 
                [](auto e_pair) { return std::get<0>(e_pair)->state.remove; }),
            edges.end());
    }
    
    // replace vertices in edges
    // NOTE a, c share cell a; b, d share cell b
    //      hence, new_v_a associated with cell a
    //      and new_v_b associated with cell b
    for (auto &e : {adj_edge_a, adj_edge_c}) {
        if (e->custom_links().a->state.remove) { 
            e->custom_links().a = new_v_a;
        }
        else {
            e->custom_links().b = new_v_a;
        }
    }
    for (auto &e : {adj_edge_b, adj_edge_d}) {
        if (e->custom_links().a->state.remove) {
            e->custom_links().a = new_v_b;
        }
        else {
            e->custom_links().b = new_v_b;
        }
    }

    // add new vertices to cells
    // NOTE vertex a is associated with cell a; 
    //      vertex b is associated with cell b
    //      both associated with cells c and d
    for (auto c : {adj_cell_a, adj_cell_c, adj_cell_d}) {
        c->custom_links().vertices.push_back(new_v_a);
    }
    for (auto c : {adj_cell_b, adj_cell_c, adj_cell_d}) {
        c->custom_links().vertices.push_back(new_v_b);
    }

    // make a new edge in adj_cells c and d
    // NOTE this is symmetric, directionality is given by void order_edges()
    for (auto &c : {adj_cell_c, adj_cell_d}) {
        auto& edges = c->custom_links().edges;
        for (auto it = edges.begin(); it != edges.end(); it++) {
            auto [e, flip] = *it;
            std::shared_ptr<Vertex> a; // the start of e
            if (not flip) {
                a = e->custom_links().a;
            }
            else {
                a = e->custom_links().b;
            }

            if (a == new_edge->custom_links().b) {
                // the start is end of new_edge -> insert before
                edges.insert(it, std::make_pair(new_edge, false));
                break;
            }
            else if (a == new_edge->custom_links().a) {
                // the start is start of new_edge -> insert before, but flip
                edges.insert(it, std::make_pair(new_edge, true));
                break;
            }
        }
    }

    double new_energy = get_energy({new_edge, adj_edge_a, adj_edge_b,
                                    adj_edge_c, adj_edge_d},
                                   {adj_cell_a, adj_cell_b,
                                    adj_cell_c, adj_cell_d});

    double probability = exp(-(new_energy - current_energy)/T1_barrier);
    if (random_number > probability)
    {
        this->_log->info("Aborting T1 transition, because energy increased by "
                "{} .. The probability to do this T1 transition is {}.",
                new_energy - current_energy, probability);
        
        // undo the changes
        edge->state.remove = false;
        vertex_a->state.remove = false;
        vertex_b->state.remove = false;

        _vertices_adjoint_edges[vertex_a->id()] = adjoint_edges_vertex_a;
        _vertices_adjoint_edges[vertex_b->id()] = adjoint_edges_vertex_b;
        // NOTE these are changed with add_edge(..)

        adj_edge_a->custom_links().a = adj_edge_a_copy.custom_links().a;
        adj_edge_a->custom_links().b = adj_edge_a_copy.custom_links().b;
        adj_edge_b->custom_links().a = adj_edge_b_copy.custom_links().a;
        adj_edge_b->custom_links().b = adj_edge_b_copy.custom_links().b;
        adj_edge_c->custom_links().a = adj_edge_c_copy.custom_links().a;
        adj_edge_c->custom_links().b = adj_edge_c_copy.custom_links().b;
        adj_edge_d->custom_links().a = adj_edge_d_copy.custom_links().a;
        adj_edge_d->custom_links().b = adj_edge_d_copy.custom_links().b;

        adj_cell_a->custom_links().vertices = adj_cell_a_copy.custom_links().vertices;
        adj_cell_a->custom_links().edges = adj_cell_a_copy.custom_links().edges;
        adj_cell_b->custom_links().vertices = adj_cell_b_copy.custom_links().vertices;
        adj_cell_b->custom_links().edges = adj_cell_b_copy.custom_links().edges;
        adj_cell_c->custom_links().vertices = adj_cell_c_copy.custom_links().vertices;
        adj_cell_c->custom_links().edges = adj_cell_c_copy.custom_links().edges;
        adj_cell_d->custom_links().vertices = adj_cell_d_copy.custom_links().vertices;
        adj_cell_d->custom_links().edges = adj_cell_d_copy.custom_links().edges;

        remove_edge(new_edge);
        remove_vertex(new_v_a);
        remove_vertex(new_v_b);

        return false;
    }
    else if (new_energy >= current_energy) {
        this->_log->debug("This T1 transition increases energy by {}, but "
            "given the parameter 'T1_barrier'={} the probability is {} do "
            "perform this transition non the less.",
            new_energy - current_energy, T1_barrier, probability);
    }

    // remove the entities
    remove_edge(edge);
    remove_vertex(vertex_a);
    remove_vertex(vertex_b);

    return true;
} // remove edge T1

/** Remove cell updating the topology (T2 transition)
 * 
 *  \details Removes the cell, its edges and its vertices and sets up a vertex at its
 *           centre.
 * 
 *  \warning Does not conserve the ordering of vertices
 */
template<class Model>
bool EntitiesManager<Model>::remove_cell_T2 (const std::shared_ptr<Cell> cell)
{
    if (cell->custom_links().edges.size() != 3) {
        this->_log->info("Delaying T2 transition, because the cell has more "
            "3 vertices (has {} vertices).", cell->custom_links().edges.size());
        this->_log->debug("Since correct implementation is "
            "missing, a T2 transition on this cell would violate the condition, "
            "that a vertex has 3 (or less at boundary) adjoint edges and cells. "
            "Waiting, that T1 transitions occur so that T2 becomes possible..");
        return false;
    }

    cell->state.remove = true;
    for (auto &v : cell->custom_links().vertices) {
        v->state.remove = true;
    }
    for (auto [e, flip] : cell->custom_links().edges) {
        e->state.remove = true;
    }

    auto new_vertex = add_vertex(barycenter_of(cell));

    for (auto v : cell->custom_links().vertices) {
        remove_vertex(v); // remove it from the vertex manager
    }

    for (auto [e, flip] : cell->custom_links().edges) {
        remove_edge(e); // remove it from the edge manager
    }

    // update the custom_links to objects that will be removed
    for (auto& e : edges()) {
        if (e->custom_links().a->state.remove) {
            e->custom_links().a = new_vertex;
            _vertices_adjoint_edges[new_vertex->id()].push_back(e);
        }
        else if (e->custom_links().b->state.remove) {
            e->custom_links().b = new_vertex;
            _vertices_adjoint_edges[new_vertex->id()].push_back(e);
        }
    }

    remove_cell(cell); // remove it from the cell manager

    // update the custom_links to objects that will be removed
    for (auto& c : cells()) {
        auto& vertices = c->custom_links().vertices;
        unsigned int num_vertices = vertices.size();
        vertices.erase(std::remove_if(vertices.begin(), vertices.end(),
                                      [](auto& v) { return v->state.remove; }),
                       vertices.end());
        if (num_vertices > vertices.size()) {
            vertices.push_back(new_vertex);
            _vertices_adjoint_cells[new_vertex->id()].push_back(c);
        }

        auto& edges = c->custom_links().edges;
        edges.erase(std::remove_if(edges.begin(), edges.end(),
                            [](auto& e_pair) { 
                                return std::get<0>(e_pair)->state.remove; }),
                    edges.end());
        // NOTE the edges are still ordered
    }

    return true;
} // remove edge T1

} // namespace Utopia::Models::PCPVertex
#endif