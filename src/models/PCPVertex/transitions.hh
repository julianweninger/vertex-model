#ifndef UTOPIA_MODELS_SAVANNAHETEROGENEOUS_TRANSITIONS_HH
#define UTOPIA_MODELS_SAVANNAHETEROGENEOUS_TRANSITIONS_HH

namespace Utopia {
namespace Models {
namespace PCPVertex {

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