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

const bool periodic_bc = false;
double Lx = 100., Ly = 100.;
double DX = 0.01, DY = 0.01;

struct Site {
    double x, y;
    int current_id;

    Site(double x, double y)
    :
        x(x),
        y(y),
        current_id(0)
    { }
};

using Site_ptr = std::shared_ptr<Site>;

/// The Vertex described by its coordinates (x,y)
struct Vertex : Site {
    // double x, y
    double fx, fy;

    Vertex(double x, double y)
    :
        Site(x, y),
        fx(0.),
        fy(0.)
    { }

    Vertex(const Vertex &b)
    :
        Vertex(b.x, b.y)
    {}

    ~Vertex() { }
};

using Vertex_ptr = std::shared_ptr<Vertex>;

/// The map of Vertexes to their Id
using VertexContainer = std::vector<Vertex_ptr>;

double distance(const Site_ptr &a, const Site_ptr &b) {
    sqrt(pow(b->x - b->x, 2) + pow(b->y - a->y, 2));
}

/// The Edge described by the ids of the Vertixes
struct Edge {
    const std::shared_ptr<Vertex> a, b;
    double length;
    /** NOTE manual update */

    std::shared_ptr<Edge> replace_edge;

    Edge(std::shared_ptr<Vertex> a, std::shared_ptr<Vertex> b, double l = 0.)
    :
        a(a),
        b(b),
        length(l),
        replace_edge(nullptr)
    { }

    ~Edge() { }
};


using Edge_ptr = std::shared_ptr<Edge>;

/// The Edges between two vertexes
/** described by the Ids of the two Vertexes and its length */
using EdgeContainer = std::vector<Edge_ptr>;

/// The length of an edge
double edge_length (Edge_ptr e) 
{
    double dx = e->b->x - e->a->x;
    double dy = e->b->y - e->a->y;

    if constexpr (periodic_bc) {
        if (dx <= -Lx / 2.) { dx = dx + Lx; }
        else if (dx > Lx / 2.) { dx = dx - Lx; }

        if (dy <= -Ly / 2.) { dy = dy + Ly; }
        else if (dy > Ly / 2.) { dy = dy - Ly; }
    }

    return (sqrt(pow(dx, 2.) + pow(dy, 2.)));
}

/// The Cell defined by its id, its vertexes and its area
struct Cell {
    VertexContainer vertexes;
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

    std::shared_ptr<Site> s; // for Voronoi construction

    Cell(std::shared_ptr<Site> s, EdgeContainer &es, double a = 0.)
    :
        edges_ordered(),
        area(a),
        s(s)
    {
        order_edges(es);
        return;
    }

    Cell(EdgeContainer &es, double a = 0.)
    :
        Cell(nullptr, es, a)
    { }

    ~Cell() { }

    void order_edges(EdgeContainer es) {
        edges_ordered.clear();
        vertexes.clear();

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

        // add vertexes from edges
        for (const auto e : edges_ordered) {
            if (std::get<bool>(e)) { // flip
                vertexes.push_back(std::get<Edge_ptr>(e)->b);
            }
            else {
                vertexes.push_back(std::get<Edge_ptr>(e)->a);
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
    double cell_area () {
        area = 0;

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
            area += ref_x * (ref_y + dy) - (ref_x + dx) * ref_y;

            // set ref duplicate of b -- b will be a in next step
            ref_x = ref_x + dx;
            ref_y = ref_y + dy;
        }

        if (area < 0) { area_sgn = -1; }
        // else if (area == 0) { area_sgn = 0;} // NOTE don't care
        else { area_sgn = 1;}

        area = 0.5 * abs(area);

        return area;
    }
};

using Cell_ptr = std::shared_ptr<Cell>;

/// The collection of cells
using CellContainer = std::vector<Cell_ptr>;

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif