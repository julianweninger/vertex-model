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
double DX = 0.01, DY = 0.01;

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


struct Site {
    double x, y;
    int current_id;

    /// Whether this object is to be removed 
    bool remove;

    Site(double x, double y)
    :
        x(x),
        y(y),
        current_id(0),
        remove(false)
    { }
};

/// The Vertex described by its coordinates (x,y)
struct Vertex : Site {
    // double x, y
    double fx, fy;

    /// Container of the adjacent edges
    std::vector<std::weak_ptr<Edge>> adj_edges;
    // TODO use weak ptr here

    /// Container of the the adjacent cells
    std::vector<std::weak_ptr<Cell>> adj_cells;

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

    Vertex(Site_ptr s) : Vertex(s->x, s->y)
    { }
};

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

template <bool periodic_bc>
std::pair<double, double> displacement(const Site_ptr &a, const Site_ptr &b) 
{
    double dx = b->x - a->x;
    double dy = b->y - a->y;

    if constexpr (periodic_bc) {
        if (dx <= -Lx / 2.) { dx += Lx; }
        else if (dx > Lx / 2.) { dx -= Lx; }

        if (dy <= -Ly / 2.) { dy += Ly; }
        else if (dy > Ly / 2.) { dy -= Ly; }
    }

    return std::make_pair(dx, dy);
}

template <bool periodic_bc>
double distance(const Site_ptr &a, const Site_ptr &b) {
    auto v_ab = displacement<periodic_bc>(a, b);

    return std::max(sqrt(pow(std::get<0>(v_ab), 2) + pow(std::get<1>(v_ab), 2)),
                    1e-10);
}

/// The Edge described by the ids of the Vertixes
struct Edge {    
    Vertex_ptr a, b;

    double length;
    /** NOTE manual update */

    double linetension;

    /// Container of the the adjacent cells
    std::vector<std::weak_ptr<Cell>> adj_cells;

    /// Whether this object is to be removed 
    bool remove;

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

/// The length of an edge
template <bool periodic_bc>
double edge_length (Edge_ptr e) 
{
    return distance<periodic_bc>(e->a, e->b);
}

/// The Cell defined by its id, its vertices and its area
struct Cell {
    VertexContainer vertices;
    /// The cell edges in order with the boolian flip
    /** The edges start at a random point, such that e0->a->x , e0->b->b
     *  From the endpoint b, the following edge fill continue to b', etc.
     *  Edges are thus ordered anti-clockwise from e0->a to eN->b = e0->a
     */
    std::vector<std::pair<Edge_ptr, bool>> edges_ordered;
    /** NOTE manual update with cell_area*/

    double area;
    /** NOTE manual update */
    char area_sgn;
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