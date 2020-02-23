#ifndef UTOPIA_MODELS_SAVANNAHETEROGENEOUS_TRANSITIONS_HH
#define UTOPIA_MODELS_SAVANNAHETEROGENEOUS_TRANSITIONS_HH

namespace Utopia {
namespace Models {
namespace PCPVertex {

template <bool periodic_bc, bool polarity_proteins>
EdgeContainer::iterator PCPVertex<periodic_bc,
    polarity_proteins>::T1_transition (EdgeContainer::iterator edge_it)
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
    auto [dx, dy] = displacement<periodic_bc>(*adj_cell_a->s, *adj_cell_b->s);
    auto length = sqrt(std::pow(dx, 2) + std::pow(dy, 2));
    dx = dx / length * _length_threshold;
    dy = dy / length * _length_threshold;
    auto tmp_periodic_copy = periodic_copy<periodic_bc>(*edge->b, *edge->a);
    auto centre_site = Site(0.5 * (edge->a->x + tmp_periodic_copy.x),
                            0.5 * (edge->a->y + tmp_periodic_copy.y));
    auto new_v_a = std::make_shared<Vertex>(
                            centre_site.x - dx / 2., centre_site.y - dy / 2.,
                            EdgeContainer({adj_edge_a, adj_edge_c}),
                            CellContainer({adj_cell_a, adj_cell_c,
                                            adj_cell_d}));
    auto new_v_b = std::make_shared<Vertex>(
                            centre_site.x + dx / 2., centre_site.y + dy / 2.,
                            EdgeContainer({adj_edge_b, adj_edge_d}),
                            CellContainer({adj_cell_b, adj_cell_c,
                                            adj_cell_d}));
    correct_periodic_bc<periodic_bc>(new_v_a);
    correct_periodic_bc<periodic_bc>(new_v_b);
    _vertices.push_back(new_v_a);
    _vertices.push_back(new_v_b);

    // create a new edge
    double linetension = _linetension(adj_cell_c->type, adj_cell_d->type);
    Edge_ptr new_edge = std::make_shared<Edge>(new_v_a, new_v_b,
                                               linetension, adj_cell_c,
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

template <bool periodic_bc, bool polarity_proteins>
CellContainer::iterator PCPVertex<periodic_bc,
    polarity_proteins>::T2_transition (CellContainer::iterator &cell_it) 
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
        // NOTE No new edge to add
        // NOTE Edges will be still ordered        

        // add new_v and update area
        if (num_vertices > c->vertices.size()) {
            c->vertices.push_back(new_v);
            // NOTE Vertices not ordered
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

template <bool periodic_bc, bool polarity_proteins>
CellContainer::iterator PCPVertex<periodic_bc,
    polarity_proteins>::divide_cell(CellContainer::iterator cell_it,
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
    double linetension = _linetension(cell->type, cell->type);
    auto new_edge = std::make_shared<Edge>(new_vertices[0], 
                                           new_vertices[1], linetension);
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
            if (c_weak.expired()) { continue; }
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
                                                cell->contractility,
                                                cell->type,
                                                cell->protein_concentration);
    new_cell_0->link_members();
    auto new_cell_1 = std::make_shared<Cell>(*cell_center, new_edges_cell_1, 
                                                cell->area_preferential,
                                                cell->contractility,
                                                cell->type,
                                                cell->protein_concentration);
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

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif