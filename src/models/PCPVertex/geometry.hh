#ifndef UTOPIA_MODELS_SAVANNAHETEROGENEOUS_GEOMETRY_HH
#define UTOPIA_MODELS_SAVANNAHETEROGENEOUS_GEOMETRY_HH

#include <memory>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>


namespace Utopia {
namespace Models {
namespace PCPVertex {

double Lx = 100., Ly = 100.;

struct Site;
struct Vertex;
struct Edge;
struct Cell;

using Site_ptr = std::shared_ptr<Site>;
using Vertex_ptr = std::shared_ptr<Vertex>;
using Edge_ptr = std::shared_ptr<Edge>;
using Cell_ptr = std::shared_ptr<Cell>;

using SiteContainer = std::vector<Site_ptr>;
using VertexContainer = std::vector<Vertex_ptr>;
using EdgeContainer = std::vector<Edge_ptr>;
using CellContainer = std::vector<Cell_ptr>;

/// The Site described by its coordinates (x,y)
struct Site {
    /// The position
    double x, y;

    /// An id unique within the current container of vertices
    int current_id;

    /// Whether this object is to be removed 
    bool remove;

    /// Constructor
    Site(double x, double y)
    :
        x(x),
        y(y),
        current_id(0),
        remove(false)
    { }
};

/// The Vertex as a positional node within the cell boundary network 
struct Vertex : Site {
    // double x, y
    double fx, fy;

    /// Container of the adjacent edges
    std::vector<std::weak_ptr<Edge>> adj_edges;

    /// Container of the the adjacent cells
    std::vector<std::weak_ptr<Cell>> adj_cells;

    /// Constructor
    Vertex(double x, double y, EdgeContainer adj_es = {},
           CellContainer adj_cs = {})
    :
        Site(x, y),
        fx(0.),
        fy(0.),
        adj_edges(),
        adj_cells()
    {
        for (auto e : adj_es) {
            adj_edges.push_back(e);
        }
        for (auto c : adj_cs) {
            adj_cells.push_back(c);
        }
    }

    /// Derive a vertex from a Site_ptr
    Vertex(Site_ptr s, EdgeContainer adj_es = {},
           CellContainer adj_cs = {})
    : Vertex(s->x, s->y, adj_es, adj_cs)
    { }

    /// Derive a vertex from a Site
    Vertex(Site s, EdgeContainer adj_es = {},
           CellContainer adj_cs = {})
    : Vertex(s.x, s.y, adj_es, adj_cs) 
    { }
};

/// Corrects the position of s wrt to boundaries
template <bool periodic_bc>
void correct_periodic_bc(Site_ptr s)
{
    if constexpr (periodic_bc) {
        if (s->x >= Lx) { s->x -= Lx; }
        else if (s->x < 0) { s->x += Lx; }
        if (s->y >= Ly) { s->y -= Ly; }
        else if (s->y < 0) { s->y += Ly; }
    }

    return;
}

/// Creates an object-copy, that is at true distance to another object
template <bool periodic_bc>
Site periodic_copy(const Site s, const Site s_fixed)
{
    if constexpr (not periodic_bc) {
        return Site(s.x, s.y);
    }

    double dx = s.x - s_fixed.x;
    double dy = s.y - s_fixed.y;

    if (dx <= -Lx / 2.) { dx += Lx; }
    else if (dx > Lx / 2.) { dx -= Lx; }

    if (dy <= -Ly / 2.) { dy += Ly; }
    else if (dy > Ly / 2.) { dy -= Ly; }

    return Site(s_fixed.x + dx, s_fixed.y + dy);
}

/// The displacement vector from a to b
template <bool periodic_bc>
std::pair<double, double> displacement(const Site &a, const Site &b) 
{
    double dx = b.x - a.x;
    double dy = b.y - a.y;

    if constexpr (periodic_bc) {
        if (dx <= -Lx / 2.) { dx += Lx; }
        else if (dx > Lx / 2.) { dx -= Lx; }

        if (dy <= -Ly / 2.) { dy += Ly; }
        else if (dy > Ly / 2.) { dy -= Ly; }
    }

    return std::make_pair(dx, dy);
}

/// The distance of two sites
template <bool periodic_bc>
double distance(const Site &a, const Site &b) {
    auto v_ab = displacement<periodic_bc>(a, b);

    return std::max(sqrt(pow(std::get<0>(v_ab), 2) + pow(std::get<1>(v_ab), 2)),
                    1e-10);
}

/// The Edge object
/** It is defined as the shortest connection between two vertices a and b
 */
struct Edge {
    /// Start and end vertices of this edge
    Vertex_ptr a, b;

    /// Length of the edge
    /** \warning manual update */
    double length;

    /// Update the length of the edge
    template <bool periodic_bc>
    double update_length() {
        length = distance<periodic_bc>(*a, *b);
        return length;
    }

    /// The linetension parameter property to this edge
    double linetension;

    /// Container of the the adjacent cells
    std::vector<std::weak_ptr<Cell>> adj_cells;

    /// Whether this object is to be removed 
    bool remove;

    /// Constructor
    /** \param a    Start Vertex
     *  \param b    End Vertex
     *  \param linetension  The linetension property
     *  \param adj_cs   The adjacent cells to this edge
     *  \param l    The distance from a to b, length of this edge
     */
    Edge(Vertex_ptr a, Vertex_ptr b, double linetension,
         CellContainer adj_cs = {}, double l = 0.)
    :
        a(a),
        b(b),
        length(l),
        linetension(linetension),
        adj_cells(),
        remove(false)
    {
        for (auto c : adj_cs) {
            adj_cells.push_back(c);
        }
    }
};

/// Intersection site of two Edges
/** \param e0           the first edge
 *  \param e1           the second edge
 *  \param finite_e0    if false, edge a extends beyond the vertices that define
 *                      its direction
 *  \param finite_e0    if false, edge b extends beyond the vertices that define
 *                      its direction
 *  
 *  \tparam periodic_bc the periodicity of space, periodic if true, non-periodic
 *                      otherwise.
 *                      
 *  \warning    The approach in periodic boundary conditions creates a periodic
 *              copy of e1, such that e1->a is considered a periodic copy of 
 *              itself wrt e0->a, the end vertices are then periodic copies wrt
 *              the start vertex of the resp edge. This approach should be
 *              possible as long as an edge is considered the shortest
 *              connection between the start and end vertex. To be double checked.
 *              See periodic_copy(const Site_ptr s, const Site_ptr s_fixed)
 * 
 *  \return shared pointer to the intersection of edges e0 and e1, resp. 
 *          infinite lines through their start and endpoint. Nullptr if no
 *          intersection.
 */
template <bool periodic_bc>
Site_ptr intersection(const Edge e0, const Edge e1, 
                     const bool finite_e0 = true, const bool finite_e1 = true) 
{
    // the start vertices of the resp. edges
    Site va = Site(e0.a->x, e0.a->y);
    Site vc = periodic_copy<periodic_bc>(Site(e1.a->x, e1.a->y), va);
    
    // the end vertices of the resp. edges
    Site vb = periodic_copy<periodic_bc>(Site(e0.b->x, e0.b->y), va);
    Site vd = periodic_copy<periodic_bc>(Site(e1.b->x, e1.b->y), vc);

    auto vector_e0 = displacement<periodic_bc>(va, *(e0.b));
    auto vector_e1 = displacement<periodic_bc>(*(e1.a), *(e1.b));

    // line defined by e0: y0 = a + bx
    double a, b;
    if (std::get<0>(vector_e0) == 0.) {
        vb.x += 1e-8;
        vector_e0 = displacement<periodic_bc>(va, vb);
    }
    b = std::get<1>(vector_e0) / std::get<0>(vector_e0);
    a = va.y - b * va.x;

    // line defined by e1: y1 = c + dx
    double c, d;
    if (std::get<0>(vector_e1) == 0.) {
        vd.x += 1e-8;
        vector_e1 = displacement<periodic_bc>(vc, vd);
    }
    d = std::get<1>(vector_e1) / std::get<0>(vector_e1);
    c = vc.y - d * vc.x;

    // parallel lines
    if (b == d) {
        return nullptr;
    }

    // coordinates of line intersection
    double x = (c-a) / (b-d);
    double y = a + b*x;

    // check validity with edges
    if (finite_e0) {
        double ax, bx;
        ax = va.x;
        bx = ax + std::get<0>(vector_e0);

        if (ax <= bx and (x < ax - 1e-10 or x > bx + 1e-10)) {
            return nullptr;
        }
        if (ax > bx and (x > ax + 1e-10 or x < bx - 1e-10)) {
            return nullptr;
        }
    }
    if (finite_e1) {
        double ax, bx;
        ax = vc.x;
        bx = ax + std::get<0>(vector_e1);

        if (ax <= bx and (x < ax - 1e-10 or x > bx + 1e-10)) {
            return nullptr;
        }
        if (ax > bx and (x > ax + 1e-10 or x < bx - 1e-10)) {
            return nullptr;
        }
    }

    auto site = std::make_shared<Site>(x, y);
    correct_periodic_bc<periodic_bc>(site);

    return site;
}

/// The Cell defined by its id, its vertices and its area
struct Cell {
    /// The vertices that bound the cell
    VertexContainer vertices;

    /// The cell edges in order (clockwise or anti-clockwise)
    /** The edges start with a random edge e0.
     *  From the endpoint e0->b, the following edge connects e0->b = e1->a with
     *  b' = e1->b, etc.
     *  Edges are thus ordered clockwise or anti-clockwise from e0->a to 
     *  eN->b = e0->a.
     * 
     *  The edges have a priori random assignment of start and endpoint a and b,
     *  which are thus exchangable.
     *  Therefore the boolian refers to the flip condition: if flip, then
     *  exchange a and b in this edge to obtain a closed loop of vertices.
     * 
     *  The order of the edges is eventually clockwise, or anti-clockwise.
     *  When calculating the area of a polygon clockwise-ordered vertices 
     *  result in a negative value for the area, which absolute value
     *  corresponds to the area of the anti-clockwise ordered equivalent.
     * 
     *  NOTE manual update with cell_area
     */
    std::vector<std::pair<Edge_ptr, bool>> edges_ordered;

    /// The area of the cell
    /** NOTE manual update */
    double area;
    
    /// The order of the edges
    /** If 1, then the edges are ordered anti-clockwise. Otherwise ordered 
     *  clockwise.
     * 
     *  NOTE updated together with area
     */
    char area_sgn;

    /// center of the cell
    /** NOTE updated together with area */
    Site_ptr s;

    double area_preferential;

    /// Whether this object is to be removed 
    bool remove;

    Cell(Site_ptr s, EdgeContainer &es, double area_preferential)
    :
        vertices(),
        edges_ordered(),
        area(0.),
        area_sgn(0),
        s(s),
        area_preferential(area_preferential),
        remove(false)
    {
        order_edges(es);
    }

    void order_edges(EdgeContainer es) {
        edges_ordered.clear();
        vertices.clear();

        if (es.empty()) { return; }

        Edge_ptr e = es.back();
        Vertex_ptr a, b;

        // put b to the right of a
        double dx = e->b->x - e->a->x;
        if ((dx < 0 and dx > -Lx / 2.) or dx > Lx / 2.) {
            edges_ordered.push_back(std::make_pair(e, true));
            a = e->b;
            b = e->a;
        }
        else {
            edges_ordered.push_back(std::make_pair(e, false));
            a = e->a;
            b = e->b;
        }
        es.pop_back();
        const auto a0 = a;

        int cnt_iteration = 0;
        int cnt_max = 2*(es.size()+1);
        while (a0 != b)
        {
            for (auto e_it = es.begin(); e_it != es.end(); /*void*/) {
                e = *e_it;
                if (e->a == b) {
                    edges_ordered.push_back(std::make_pair(e, false));
                    a = e->a;
                    b = e->b;

                    e_it = es.erase(e_it);
                    break;
                }
                else if (e->b == b) {
                    edges_ordered.push_back(std::make_pair(e, true));
                    a = e->b;
                    b = e->a;

                    e_it = es.erase(e_it);
                    break;
                }
                else {
                    ++e_it;
                }
            }

            #ifndef NDEBUG
            if (cnt_iteration++ > cnt_max) {
                std::cout << "\nFailed to order these edges: \n";
                for (auto e_pair : edges_ordered) {
                    auto e = std::get<Edge_ptr>(e_pair);
                    if (not std::get<bool>(e_pair)) {
                        std::cout << " (" << e->a->x << ", " << e->a->y << ") "
                            "to (" << e->b->x << ", " << e->b->y << ")\n";
                    }
                    else {
                        std::cout << " (" << e->b->x << ", " << e->b->y << ") "
                            "to (" << e->a->x << ", " << e->a->y << ")\n";
                    }
                }
                std::cout << " stopped here.\n";
                for (auto e : es) {
                    std::cout << " (" << e->a->x << ", " << e->a->y << ") "
                        "to (" << e->b->x << ", " << e->b->y << ")\n";
                }
                std::cout << std::flush;

                throw std::runtime_error("Could not order edges of cell at (" +
                    std::to_string(s->x) + ", " + std::to_string(s->y) + "). "
                    "Edges to order are listed above.");
            }
            #endif
        }

        // add vertices from edges
        for (const auto e : edges_ordered) {
            if (std::get<bool>(e)) { // flip
                vertices.push_back(std::get<Edge_ptr>(e)->b);
            }
            else {
                vertices.push_back(std::get<Edge_ptr>(e)->a);
            }
        }
    }

    void order_edges() {
        EdgeContainer es;
        for (const auto e_pair : edges_ordered) {
            es.push_back(std::get<Edge_ptr>(e_pair));
        }
        order_edges(es);
    }

    /// The area of a polygon cell
    template <bool periodic_bc>
    double cell_area () {
        area = 0;
        double center_x = 0;
        double center_y = 0;

        double ref_x, ref_y; // the reference vertex
        if (std::get<bool>(edges_ordered.front())) {
            auto v_ref = std::get<Edge_ptr>(edges_ordered.front())->b;
            ref_x = v_ref->x;
            ref_y = v_ref->y;
        }
        else {
            auto v_ref = std::get<Edge_ptr>(edges_ordered.front())->a;
            ref_x = v_ref->x;
            ref_y = v_ref->y;
        }

        for (const auto e_pair : edges_ordered) {
            const auto e = std::get<Edge_ptr>(e_pair);
            Vertex_ptr a, b;
            if (std::get<bool>(e_pair)) {
                a = e->b;
                b = e->a;
            }
            else {
                a = e->a;
                b = e->b;
            }

            double bx = b->x;
            double by = b->y;
            double dx = bx - ref_x;
            double dy = by - ref_y;

            if constexpr (periodic_bc) {
                if (dx <= -Lx / 2.) { dx = dx + Lx; }
                else if (dx > Lx / 2.) { dx = dx - Lx; }

                if (dy <= -Ly / 2.) { dy = dy + Ly; }
                else if (dy > Ly / 2.) { dy = dy - Ly; }
            }

            // ref is a, ref + (dx, dy) is b
            double da = ref_x * (ref_y + dy) - (ref_x + dx) * ref_y;
            area += da;
            center_x -= (ref_x + ref_x + dx) * da;
            center_y -= (ref_y + ref_y + dy) * da;

            // set ref duplicate of b -- b will be a in next step
            ref_x = ref_x + dx;
            ref_y = ref_y + dy;
        }

        if (area < 0) { area_sgn = -1; }
        // else if (area == 0) { area_sgn = 0;} // NOTE don't care
        else { area_sgn = 1;}

        area = 0.5 * area * area_sgn;
        if constexpr (periodic_bc) {
            double x = -area_sgn * center_x / 6 / area;
            double y = -area_sgn * center_y / 6 / area;
            
            if (x >= Lx) {
                x -= Lx;
            }
            else if (x < 0) {
                x += Lx;
            }
            if (y >= Ly) {
                x -= Ly;
            }
            else if (y < 0) {
                y += Ly;
            }

            s = std::make_shared<Site>(x, y);
        }
        else {
            s = std::make_shared<Site>(-area_sgn * center_x / 6 / area, 
                                       -area_sgn * center_y / 6 / area);
        }
        correct_periodic_bc<periodic_bc>(s);

        return area;
    }
};

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif