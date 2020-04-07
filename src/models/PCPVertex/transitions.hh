#ifndef UTOPIA_MODELS_SAVANNAHETEROGENEOUS_TRANSITIONS_HH
#define UTOPIA_MODELS_SAVANNAHETEROGENEOUS_TRANSITIONS_HH

namespace Utopia {
namespace Models {
namespace PCPVertex {

/** Perform a T1 transition on the edge edge_it in _edges
 * 
 *  edge_it will be removed and a new edge orthogonal to edge will be created
 */
template <bool periodic_bc>
std::pair<EdgeContainer::iterator, bool> PCPVertex<periodic_bc>::T1_transition (
        EdgeContainer::iterator edge_it)
{
    this->_log->info("Removing edge in T1 transition in step {}..",
                      this->_time);
    if constexpr (not periodic_bc) {
        if ((*edge_it)->adj_cell_a.expired() or (*edge_it)->adj_cell_b.expired()) {
            throw std::runtime_error("In T1 transition: Removing edge "
                "with less than 2 adj cells not implemented!");
        }
    }

    auto edge = *edge_it;

    Edge edge_copy = *edge;
    Vertex vertex_a_copy = *edge_copy.a;
    Vertex vertex_b_copy = *edge_copy.b;

    // tag objects to be removed
    edge->remove = true;
    edge->a->remove = true;
    edge->b->remove = true;

    // the adj cells that are currently neighbours sharing edge
    // NOTE a, b arbitrary
    Cell_ptr adj_cell_a = edge->adj_cell_a.lock();
    Cell_ptr adj_cell_b = edge->adj_cell_b.lock();
    
    for (auto c : {adj_cell_a, adj_cell_b}) {
        if (c->edges_ordered.size() <= 3) {
            this->_log->warn("The triangular cell has area {}. Adapt the "
                "threshold area for T2 transitions, such that this cell "
                "is removed, rather than performing a T1 transition!",
                c->template area<periodic_bc>());
            return std::make_pair(++edge_it, false);
            // TODO check algorithm at this point
            // throw std::runtime_error("Cannot perform T1 transition of an "
            //     "edge adjacent to a cell with only 3 edges (triangle). "
            //     "The cell would become a (n-1) dimensional object.");
        }
    }

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

    // make copies of the status quo
    Edge adj_edge_a_copy = Edge(*adj_edge_a), adj_edge_b_copy = Edge(*adj_edge_b);
    Edge adj_edge_c_copy = Edge(*adj_edge_c), adj_edge_d_copy = Edge(*adj_edge_d);
    Cell adj_cell_a_copy = Cell(*adj_cell_a), adj_cell_b_copy = Cell(*adj_cell_b);
    Cell adj_cell_c_copy = Cell(*adj_cell_c), adj_cell_d_copy = Cell(*adj_cell_d);

    double current_energy = 0.;
    // double current_energy = this->get_energy({edge, adj_edge_a, adj_edge_b,
    //                                           adj_edge_c, adj_edge_d},
    //                                          {adj_cell_a, adj_cell_b,
    //                                           adj_cell_c, adj_cell_d});

    // create two new vertices that create an edge of threshold length 
    // pointing from cell a to b
    auto [dx, dy] = displacement<periodic_bc>(
        *adj_cell_a->template centre_site<periodic_bc>(),
        *adj_cell_b->template centre_site<periodic_bc>());
    auto length = sqrt(std::pow(dx, 2) + std::pow(dy, 2));
    dx = dx / length * _T1_threshold / _Lx;
    dy = dy / length * _T1_threshold / _Ly;
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

    // create a new edge
    double linetension = _linetension(adj_cell_c->type, adj_cell_d->type);
    double contractility = _edge_contractility(adj_cell_c->type,
                                               adj_cell_d->type);
    Edge_ptr new_edge = std::make_shared<Edge>(new_v_a, new_v_b,
                                               linetension, contractility,
                                               adj_cell_c, adj_cell_d);
    new_edge->link_members();

    // remove objects
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
        c->template area<periodic_bc>();
    }

    double new_energy = 0.;
    // double new_energy = this->get_energy({new_edge, adj_edge_a, adj_edge_b,
    //                                       adj_edge_c, adj_edge_d},
    //                                      {adj_cell_a, adj_cell_b,
    //                                       adj_cell_c, adj_cell_d});

    double probability = exp(-(new_energy - current_energy)/_T1_barrier);
    if (new_energy > current_energy
        and _prob_distr(*this->_rng) > probability)
    {
        this->_log->info("Aborting T1 transition, because energy increased by "
                "{} .. The probability to do this T1 transition is {}.",
                new_energy - current_energy, probability);

        *edge = edge_copy;
        *(edge->a) = vertex_a_copy;
        *(edge->b) = vertex_b_copy;

        *adj_edge_a = adj_edge_a_copy;
        *adj_edge_b = adj_edge_b_copy;
        *adj_edge_c = adj_edge_c_copy;
        *adj_edge_d = adj_edge_d_copy;

        *adj_cell_a = adj_cell_a_copy;
        *adj_cell_b = adj_cell_b_copy;
        *adj_cell_c = adj_cell_c_copy;
        *adj_cell_d = adj_cell_d_copy;

        return std::make_pair(++edge_it, false);
    }
    else if (new_energy >= current_energy) {
        this->_log->debug("This T1 transition increases energy by {}, but "
            "given the parameter 'T1_barrier'={} the probability is {} do "
            "perform this transition non the less.",
            new_energy - current_energy, _T1_barrier, probability);
    }

    // replace the edge at adge_it
    edge_it = _edges.erase(edge_it);
    edge_it = _edges.insert(edge_it, new_edge);

    _vertices.erase(
        std::remove_if(
            _vertices.begin(), _vertices.end(),
            [](auto v) { return v->remove; }),
        _vertices.end()
    );

    _vertices.push_back(new_v_a);
    _vertices.push_back(new_v_b);

    return std::make_pair(++edge_it, true);
}
/** Erase cell from _cells while updating the topology (T2 transition)
 * 
 *  Removes the cell, its edges and its vertices and sets up a vertex at its
 *  centre.
 * 
 *  returns _cells.erase(cell_it)
 */
template <bool periodic_bc>
std::pair<CellContainer::iterator, bool> PCPVertex<
        periodic_bc>::T2_transition (CellContainer::iterator &cell_it) 
{
    this->_log->info("Removing cell in T2 transition..");

    auto cell = *cell_it;

    if (cell->edges_ordered.size() != 3) {
        this->_log->info("Delaying T2 transition, because the cell has more "
            "3 vertices. Since correct implementation is missing, a T2 "
            "transition on this cell would violate the condition, that a "
            "vertex has 3 (or less at boundary) adj_edges and _cells. "
            "Hope, that T1 transitions occur so that T2 becomes possible.");
        return std::make_pair(++cell_it, false);
    }
    
    // create a new vertex at the center of c
    cell->template area<periodic_bc>(); // updates the center of c        
    auto new_v = std::make_shared<Vertex>(cell->template centre_site<periodic_bc>());
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
            new_v->adj_cells.push_back(c);
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
            new_v->adj_edges.push_back(e);
        }
        else if (e->b->remove) {
            e->b = new_v;
            new_v->adj_edges.push_back(e);
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
    return std::make_pair(_cells.erase(cell_it), true);
}

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif