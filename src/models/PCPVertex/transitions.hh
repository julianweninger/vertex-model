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
template <class Model>
template <typename EdgeParamMatrix>
void EntitiesManager<Model>::divide_cell(const std::shared_ptr<Cell> cell,
        double division_angle,
        EdgeParamMatrix linetension, EdgeParamMatrix edge_contractility)
{
    this->_log->debug("Dividing cell ...");

    // The cell to be divided
    cell->state.remove = true;
    remove_cell(cell);

    const SpaceVec cell_center = this->barycenter_of(cell);

    // generate the axis of division
    const SpaceVec axis = SpaceVec({cos(division_angle), sin(division_angle)});

    // determine the new vertices from this axis
    // These are the intersections of the division axis with edges of cell
    // NOTE there must be exactly 2 intersections for nicely shaped cells
    AgentContainer<Vertex> new_vertices;
    for (const auto& [e, flip] : cell->custom_links().edges) {
        // calculate the intersection
        SpaceVec a = this->position_of(e->custom_links().a);
        SpaceVec b = this->position_of(e->custom_links().b);
        SpaceVec inter;
        bool success;
        std::tie(inter, success) = this->_space->intersection(
                a, this->_space->displacement(a, b),
                cell_center, axis,
                true, false);
        
        // if no intersection end here
        if (not success) { continue; }

        // create a new vertex at the intersection site
        auto new_v = this->add_vertex(inter);
        const auto [adj_cell_a, adj_cell_b] = this->adjoints_of(e);
        
        AgentContainer<Cell> adjoint_cells;
        if (adj_cell_a) { adjoint_cells.push_back(adj_cell_a); }
        if (adj_cell_b) { adjoint_cells.push_back(adj_cell_b); }
        
        _vertices_adjoint_cells[new_v->id()] = adjoint_cells;
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
            SpaceVec a = this->position_of(e->custom_links().a);
            SpaceVec b = this->position_of(e->custom_links().b);
            if (flip) { std::swap(a, b); }
            this->_log->error("   {}, {} -> {}, {}",a[0], a[1], b[0], b[1]);
        }
        this->_log->error("Intersections are:");
        for (const auto& v : new_vertices) {
            SpaceVec pos = this->position_of(v);
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
    edge_cfg["linetension"] = linetension.at(cell->state.type,
                                             cell->state.type);
    edge_cfg["contractility"] = edge_contractility.at(cell->state.type,
                                                      cell->state.type);
    const auto new_edge = add_edge(new_vertices[0], new_vertices[1]);

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

        DataIO::Config edge_cfg = edge->state.get_cfg();

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
    const auto new_cell_0 = this->add_cell(cell_center, new_edges_cell_0, cell);
    const auto new_cell_1 = this->add_cell(cell_center, new_edges_cell_1, cell);
    
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

            if (adj_cell_a == new_c and adj_cell_b == new_c) {
                adj_cell_b = nullptr;
            }
            
            _edges_adjoint_cells[e->id()] = std::make_pair(adj_cell_a,
                                                           adj_cell_b);
        }
    }

    this->_log->trace("Successfully divided cell.");

    return;
} // divide cell


/// Removes an edge in a T1 neighborhood exchange
/** \details Adds a new edge separating the adjoint cells
 * 
 *  \param edge         The edge to remodel
 *  \param linetension  The linetension of the new edge
 *  \param contractility  The contractility of the new edge
 *  \param get_energy   Calculate the energy for container of edges and cells
 *  \param separation   The separation of the new vertices 
 *  \param T1_barrier   The height of the energy barrier
 *  \param random_number    A random number in [0, 1]
 *  \param time         Timepoint this function is called. Sets last_T1_attempt
 *                      in new junction
 * 
 *  \returns whether edge was removed.
 *           The edge is not removed if one of the adjoint cells is triangular.
 *           Remove it rather in a T2 transition.
 *           The edge is also not removed if the transition increases energy.
 */
template <class Model>
template <typename EdgeParamMatrix>
bool EntitiesManager<Model>::remove_edge_T1 (const std::shared_ptr<Edge> edge,
        EdgeParamMatrix linetension, EdgeParamMatrix contractility,
        std::function<double(const AgentContainer<Edge>&,
                             const AgentContainer<Cell>&)> get_energy,
        double separation, double T1_barrier, double random_number,
        std::size_t time)
{
    auto vertex_a = edge->custom_links().a;
    auto vertex_b = edge->custom_links().b;

    if (   is_2_fold_boundary_vertex(vertex_a)
        or is_2_fold_boundary_vertex(vertex_b))
    {
        return remove_boundary_edge(edge, get_energy, T1_barrier,
                                    random_number);
    }

    if (is_1_cell_boundary_edge(edge)) {
        this->_log->debug("Remodelling 1 cell boundary edge {} ({} -> {}) in "
                          "T1 transition ...", edge->id(),
                          vertex_a->id(), vertex_b->id());

    }
    else if (is_2_cell_boundary_edge<false>(edge)) {
        if (    is_3_fold_boundary_vertex(vertex_a)
            and is_3_fold_boundary_vertex(vertex_b))
        {
            this->_log->error("Cannot remove edge {} with two 3-fold boundary "
                              "vertices {} and {}. This would separate the "
                              "tissue into two non-connected parts. Not "
                              "implemented!", edge->id(), vertex_a->id(),
                              vertex_b->id());
            return false;
        }

        this->_log->debug("Remodelling 2 cell boundary edge {} in "
                          "T1 transition ...", edge->id());
    }
    else {
        this->_log->debug("Remodelling edge {} in T1 transition ...",
                          edge->id());
    }

    
    SpaceVec displ = displacement(vertex_a, vertex_b);

    // the vector separating the two new vertices orthogonal to displ
    // new edge will be of length `separation`
    SpaceVec sep = SpaceVec({displ[1], -displ[0]})/arma::norm(displ)*separation;

    // choose a the cell right of edge and b left
    auto [adj_cell_a, adj_cell_b] = adjoints_of(edge);
    // is a boundary edge ?
    if (not adj_cell_a or not adj_cell_b) {
        // make sure the void (nullptr) cell is adj_cell_b
        if (not adj_cell_a) {
            std::swap(adj_cell_a, adj_cell_b);
        }
        // NOTE adj_cell_b is nullptr now
        //      adj_cell_a is always valid pointer

        // make sure that adj_cell_a is to the right of edge
        SpaceVec center = position_of(vertex_a) + 0.5 * displ;
        SpaceVec pos_a = barycenter_of(adj_cell_a);
        if (arma::dot(_space->displacement(pos_a, center), sep) < 0) {
            // the vector A->B is anti-parallel to edge rotated by 90deg
            // anti-clockwise, hence flip edge
            std::swap(vertex_a, vertex_b);  
        
            displ = displacement(vertex_a, vertex_b);
            sep = SpaceVec({displ[1], -displ[0]})/arma::norm(displ)*separation;          
        }
    }
    else if (arma::dot(displacement(adj_cell_a, adj_cell_b), sep) < 0) {
        // the vector A->B is anti-parallel to edge rotated by 90deg
        // anti-clockwise, hence flip edge
        std::swap(vertex_a, vertex_b);
        
        displ = displacement(vertex_a, vertex_b);
        sep = SpaceVec({displ[1], -displ[0]})/arma::norm(displ)*separation;
    }
    
    
    for (auto c : {adj_cell_a, adj_cell_b}) {
        if (not c) { continue; }
        if (c->custom_links().edges.size() <= 3) {
            this->_log->warn("The triangular cell has area {}. Adapt the "
                "threshold area for T2 transitions, such that this cell "
                "is removed, rather than performing a T1 transition!",
                area_of(c));
            return false;
        }
    }

    // create copies of the objects - T1 might be aborted
    const Edge edge_copy = *edge;
    const auto adjoint_edges_vertex_a = adjoint_edges_of(vertex_a);
    const auto adjoint_edges_vertex_b = adjoint_edges_of(vertex_b);

    // tag objects to be removed
    edge->state.remove = true;
    vertex_a->state.remove = true;
    vertex_b->state.remove = true;

    // the two cells that become neighbours in this T1 transition
    // NOTE c is the cell adjoint to vertex edge->a
    //      d is the cell adjoint to vertex edge->b
    //      is nullptr if vertex is boundary vertex
    std::shared_ptr<Cell> adj_cell_c = nullptr;
    std::shared_ptr<Cell> adj_cell_d = nullptr;
    for (const auto& c : adjoint_cells_of(vertex_a)) {
        if (c != adj_cell_a and c != adj_cell_b) {
            adj_cell_c = c;
            break;
        }
    }
    for (const auto& c : adjoint_cells_of(vertex_b)) {
        if (c != adj_cell_a and c != adj_cell_b) {
            adj_cell_d = c;
            break;
        }
    }

    AgentContainer<Cell> cells = {adj_cell_a};
    if (adj_cell_b) { cells.push_back(adj_cell_b); }
    if (adj_cell_c) { cells.push_back(adj_cell_c); }
    if (adj_cell_d) { cells.push_back(adj_cell_d); }

    // The involved edges
    // (a, b) and (c, d) currently share vertex a, resp. b
    // (a, c) and (b, d) currently share an adjoint cell a, resp. b
    // they will share a vertex a (resp. b) after transition
    std::shared_ptr<Edge> adj_edge_a, adj_edge_b, adj_edge_c, adj_edge_d;
    for (const auto& e : adjoint_edges_of(vertex_a)) {
        if (e == edge) {
            continue;
        }
        // e is either adj_edge_a or _b 

        const auto& [e_adj_cell_a, e_adj_cell_b] =  adjoints_of(e);
        if (e_adj_cell_a == adj_cell_a or e_adj_cell_b == adj_cell_a)
        {
            adj_edge_a = e;
        }
        else {
            adj_edge_b = e;
        }
    }
    for (const auto& e : adjoint_edges_of(vertex_b)) {
        if (e == edge) {
            continue;
        }
        // e is either adj_edge_c or _d

        const auto& [e_adj_cell_a, e_adj_cell_b] =  adjoints_of(e);
        if (e_adj_cell_a == adj_cell_a or e_adj_cell_b == adj_cell_a)
        {
            adj_edge_c = e;
        }
        else {
            adj_edge_d = e;
        }
    }
    AgentContainer<Edge> edges = {adj_edge_a, adj_edge_c};
    if (adj_edge_b) { edges.push_back(adj_edge_b); }
    if (adj_edge_d) { edges.push_back(adj_edge_d); }


    // make copies of the status quo
    std::unordered_map<std::size_t, const Edge> adj_edge_copies;
    for (const auto& e : edges) {
        adj_edge_copies.insert({e->id(), *e});
    }

    std::unordered_map<std::size_t,
                       const AgentContainer<Vertex>> adj_cells_vs_copies;
    std::unordered_map<std::size_t,
                       const OrderedEdgeContainer> adj_cells_es_copies;
    for (const auto& c : cells) {
        adj_cells_vs_copies.insert({c->id(), c->custom_links().vertices});
        adj_cells_es_copies.insert({c->id(), c->custom_links().edges});
    }


    // the current energy
    AgentContainer<Edge> edges_tmp = edges;
    edges_tmp.push_back(edge);
    double current_energy = get_energy(edges_tmp, cells);
    // WARN agents are not yet removed at this point!
    //      Any component if get_energy that works globally on the agent
    //      container must be dealt with care!

    // create two new vertices using the separation
    SpaceVec center = (  position_of(vertex_a)
                       + 0.5 * displacement(vertex_a, vertex_b));
    auto new_v_a = add_vertex(center - sep / 2.);
    auto new_v_b = add_vertex(center + sep / 2.);

    _vertices_adjoint_edges[new_v_a->id()] = {adj_edge_a, adj_edge_c};    
    AgentContainer<Edge> adj_edges_tmp = {};
    if (adj_edge_b) { adj_edges_tmp.push_back(adj_edge_b); }
    if (adj_edge_d) { adj_edges_tmp.push_back(adj_edge_d); }
    _vertices_adjoint_edges[new_v_b->id()] = adj_edges_tmp;

    AgentContainer<Cell> adjoint_cells_a = {adj_cell_a};
    AgentContainer<Cell> adjoint_cells_b = {};
    if (adj_cell_b) {
        adjoint_cells_b.push_back(adj_cell_b);
    }
    for (const auto& c: {adj_cell_c, adj_cell_d}) {
        if (not c) { continue; }
        
        adjoint_cells_a.push_back(c);
        adjoint_cells_b.push_back(c);
    }
    _vertices_adjoint_cells[new_v_a->id()] = adjoint_cells_a;
    _vertices_adjoint_cells[new_v_b->id()] = adjoint_cells_b;

    // create a new edge
    DataIO::Config edge_cfg;
    if (adj_cell_c and adj_cell_d) {
        edge_cfg["linetension"] = linetension.at(adj_cell_c->state.type,
                                                 adj_cell_d->state.type);
        edge_cfg["contractility"] = contractility.at(adj_cell_c->state.type,
                                                     adj_cell_d->state.type);
    }
    else if (adj_cell_c) {
        edge_cfg["linetension"] = linetension.at(adj_cell_c->state.type,
                                                 CellType::num_cell_types);
        edge_cfg["contractility"] = contractility.at(adj_cell_c->state.type,
                                                     CellType::num_cell_types);
    }
    else if (adj_cell_d) {
        edge_cfg["linetension"] = linetension.at(adj_cell_d->state.type,
                                                 CellType::num_cell_types);
        edge_cfg["contractility"] = contractility.at(adj_cell_d->state.type,
                                                     CellType::num_cell_types);
    }

    auto new_edge = add_edge(new_v_a, new_v_b, edge_cfg);
    _edges_adjoint_cells[new_edge->id()] = std::make_pair(adj_cell_c,
                                                          adj_cell_d);
    // NOTE OK also if adj_cell_c or _d nullptr


    // remove objects
    for (auto &c : cells) {
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
        if (not c) { continue; } 

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
        if (not e) { continue; }

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
        if (not c) { continue; }
        c->custom_links().vertices.push_back(new_v_a);
    }
    for (auto c : {adj_cell_b, adj_cell_c, adj_cell_d}) {
        if (not c) { continue; }
        c->custom_links().vertices.push_back(new_v_b);
    }

    // make a new edge in adj_cells c and d
    // NOTE this is symmetric, directionality is given by void order_edges()
    for (auto &c : {adj_cell_c, adj_cell_d}) {
        if (not c) { continue; }

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

    edges_tmp = edges;
    edges_tmp.push_back(new_edge);
    double new_energy = get_energy(edges_tmp, cells);
    // WARN agents are not yet removed at this point!
    //      Any component if get_energy that works globally on the agent
    //      container must be dealt with care!

    double probability = exp(-(new_energy - current_energy)/T1_barrier);
    if (random_number > probability)
    {
        this->_log->debug("Aborting T1 transition, because energy increased by "
                "{} .. The probability to do this T1 transition is {}.",
                new_energy - current_energy, probability);
        
        // undo the changes
        edge->state.remove = false;
        vertex_a->state.remove = false;
        vertex_b->state.remove = false;

        _vertices_adjoint_edges[vertex_a->id()] = adjoint_edges_vertex_a;
        _vertices_adjoint_edges[vertex_b->id()] = adjoint_edges_vertex_b;
        // NOTE these are changed with add_edge(..)

        // reset the custom links
        for (const auto& e : edges) {
            e->custom_links().a = adj_edge_copies.at(e->id()).custom_links().a;
            e->custom_links().b = adj_edge_copies.at(e->id()).custom_links().b;
        }
        for (const auto& c : cells) {
            c->custom_links().vertices = adj_cells_vs_copies.at(c->id());
            c->custom_links().edges = adj_cells_es_copies.at(c->id());
        }

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

    new_edge->state.last_T1_attempt = time;

    return true;
} // remove edge T1

/// Removes a 1 cell boundary edge
/** Removes the edge and merges the two vertices in a single new vertex.
 * 
 *  \param edge         The edge to remodel
 *  \param get_energy   Calculate the energy for container of edges and cells
 *  \param T1_barrier   The height of the energy barrier
 *  \param random_number    A random number in [0, 1]
 * 
 *  \returns whether edge was removed.
 *           The edge is not removed if one of the adjoint cells is triangular.
 *           Remove it rather in a T2 transition.
 *           The edge is also not removed if the transition increases energy.
 */
template<class Model>
bool EntitiesManager<Model>::remove_boundary_edge(
        const std::shared_ptr<Edge> edge,
        std::function<double(const AgentContainer<Edge>&,
                             const AgentContainer<Cell>&)> get_energy,
        double T1_barrier, double random_number)
{
    auto vertex_a = edge->custom_links().a;
    auto vertex_b = edge->custom_links().b;

    // if a 3 fold vertex involved, call this one a
    if (is_3_fold_boundary_vertex(vertex_b))
    {
        std::swap(vertex_a, vertex_b);
    }

    this->_log->debug("Removing a boundary edge {} in merging "
                      "vertices {} and {} ...",
                      edge->id(), vertex_a->id(), vertex_b->id());
    
    // get the adjacent cells to edge. `adj_cell_b` will be the void cell
    auto [adj_cell_a, adj_cell_b] = adjoints_of(edge);
    if (not adj_cell_a) {
        std::swap(adj_cell_a, adj_cell_b);
    }
    if (not adj_cell_a) {
        throw std::runtime_error(fmt::format(
                "The 1 cell boundary edge {} to be removed in T1 transition "
                "has no adjoint cell! It is expected to have 1 adjacent cell "
                "and 1 void neighbor (nullptr).", edge->id()));
    }
    
    if (adj_cell_a->custom_links().edges.size() <= 3) {
        this->_log->warn("The triangular cell has area {}. Adapt the "
            "threshold area for T2 transitions, such that this cell "
            "is removed, rather than performing a T1 transition!",
            area_of(adj_cell_a));
        return false;
    }

    // tag objects to be removed
    edge->state.remove = true;
    vertex_a->state.remove = true;
    vertex_b->state.remove = true;

    // the two cells that become neighbours in this T1 transition
    // NOTE c is the cell adjoint to vertex edge->a
    //      is nullptr if vertex is boundary vertex
    std::shared_ptr<Cell> adj_cell_c = nullptr;
    if (is_3_fold_boundary_vertex(vertex_a)) {
        for (const auto& c : adjoint_cells_of(vertex_a)) {
            if (c != adj_cell_a and c != adj_cell_b) {
                adj_cell_c = c;
                break;
            }
        }
    }

    AgentContainer<Cell> adj_cells = {adj_cell_a};
    if (adj_cell_c) { adj_cells.push_back(adj_cell_c); }

    // The involved edges
    // (a, b) and (c) currently share vertex a, resp. b
    // b is nullptr for a 2_fold_vertex
    // (a, c) currently share an adjoint cell a, b is boundary edge
    std::shared_ptr<Edge> adj_edge_a, adj_edge_b, adj_edge_c;
    adj_edge_b = nullptr; // NOTE adj_edge_b can be nullptr
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
        // there are 2 adjoint edges, the one that is not `edge` is adj_edge_c
        if (e != edge) {
            adj_edge_c = e;
        }
    }

    AgentContainer<Edge> adj_edges = {edge, adj_edge_a, adj_edge_c};
    if (adj_edge_b) { adj_edges.push_back(adj_edge_b); }


    // create copies of the objects - T1 might be aborted
    const Edge edge_copy = *edge;
    const auto adjoint_edges_vertex_a = adjoint_edges_of(vertex_a);
    const auto adjoint_edges_vertex_b = adjoint_edges_of(vertex_b);

    std::vector<Edge> adj_edges_copy{};
    for (const auto& edge : adj_edges) {
        adj_edges_copy.push_back(*edge);
    }
    std::vector<AgentContainer<Vertex>> adj_cells_vs_copy{};
    std::vector<OrderedEdgeContainer> adj_cells_es_copy{};
    for (const auto& cell : adj_cells) {
        adj_cells_vs_copy.push_back(cell->custom_links().vertices);
        adj_cells_es_copy.push_back(cell->custom_links().edges);
    }

    double current_energy = get_energy(adj_edges, adj_cells);
    // WARN agents are not yet removed at this point!
    //      Any component if get_energy that works globally on the agent
    //      container must be dealt with care!

    // create a new vertex in the edge's center
    SpaceVec center = (  position_of(vertex_a)
                       + 0.5 * displacement(vertex_a, vertex_b));
    auto new_v = add_vertex(center);

    // assign the new adjoint edges and cells
    AgentContainer<Edge> new_adjoint_edges = {adj_edge_a, adj_edge_c};
    if (adj_edge_b) { new_adjoint_edges.push_back(adj_edge_b); }    
    _vertices_adjoint_edges[new_v->id()] = new_adjoint_edges;

    _vertices_adjoint_cells[new_v->id()] = adj_cells;


    // remove objects
    for (auto &c : adj_cells) {
        auto& vertices = c->custom_links().vertices;
        vertices.erase(std::remove_if(
                vertices.begin(), vertices.end(), 
                [](const auto& v) { return v->state.remove; }),
            vertices.end());
    }
    // NOTE edge_it will be removed at the very end
    // NOTE a and b will be replaced within edges

    // remove edge from adj_cell_a
    // NOTE the neighbouring edges now have a common vertex, hence order of
    //      edges is maintained
    adj_cell_a->custom_links().edges.erase(
        std::remove_if(
            adj_cell_a->custom_links().edges.begin(),
            adj_cell_a->custom_links().edges.end(), 
            [](auto e_pair) { return std::get<0>(e_pair)->state.remove; }),
        adj_cell_a->custom_links().edges.end());
    
    // replace vertices in edges
    for (auto &e : adj_edges) {
        if (e->custom_links().a->state.remove) { 
            e->custom_links().a = new_v;
        }
        else {
            e->custom_links().b = new_v;
        }
    }

    // add new vertices to cells
    // NOTE vertex a is associated with cell a; 
    //      vertex b is associated with cell b
    //      both associated with cells c and d
    for (auto c : adj_cells) {
        c->custom_links().vertices.push_back(new_v);
    }

    double new_energy = get_energy(new_adjoint_edges, adj_cells);
    // WARN agents are not yet removed at this point!
    //      Any component if get_energy that works globally on the agent
    //      container must be dealt with care!

    double probability = exp(-(new_energy - current_energy)/T1_barrier);
    if (random_number > probability)
    {
        this->_log->debug("Aborting T1 transition, because energy increased by "
                "{} .. The probability to do this T1 transition is {}.",
                new_energy - current_energy, probability);
        
        // undo the changes        
        edge->state.remove = false;
        vertex_a->state.remove = false;
        vertex_b->state.remove = false;

        _vertices_adjoint_edges[vertex_a->id()] = adjoint_edges_vertex_a;
        _vertices_adjoint_edges[vertex_b->id()] = adjoint_edges_vertex_b;
        // NOTE these are changed with add_edge(..)

        for (std::size_t i = 0; i < adj_edges.size(); i++) {
            auto edge = adj_edges[i];
            auto edge_copy = adj_edges_copy[i];
            edge->custom_links().a = edge_copy.custom_links().a;
            edge->custom_links().b = edge_copy.custom_links().b;

        }
        for (std::size_t i = 0; i < adj_cells.size(); i++) {
            auto cell = adj_cells[i];
            cell->custom_links().vertices = adj_cells_vs_copy[i];
            cell->custom_links().edges = adj_cells_es_copy[i];
        }

        remove_vertex(new_v);

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
} // remove boundary edge T1


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
    if (neighbors_of(cell).size() > 3) {
        this->_log->debug("Delaying T2 transition, because the cell {} has "
            "more than 3 neighbors (has {} neighbors and {} edges).",
            cell->id(), neighbors_of(cell).size(),
            cell->custom_links().edges.size());
        return false;
    }
    if (cell->custom_links().edges.size() > 3) {
        this->_log->debug("Removing boundary cell {} with {} neighbor(s) and "
            "{} edges in T2 transition ...", cell->id(), 
            neighbors_of(cell).size(), cell->custom_links().edges.size());
    }
    else if (not this->is_boundary(cell)) {
        this->_log->debug("Removing cell {} in T2 transition ...", cell->id());
    }
    else {
        this->_log->debug("Removing boundary cell {} in T2 transition ...",
                          cell->id());
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