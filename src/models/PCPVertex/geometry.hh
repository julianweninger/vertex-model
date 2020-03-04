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

/// The Site described by its relative coordinates (x,y)
struct Site {
    /// The relative x position within the domain
    double x;

    /// The relative y position within the domain
    double y;

    /// An id unique within the current container of vertices
    int current_id;

    /// Whether this object is to be removed 
    bool remove;

    /// Constructor
    /** \param x, y     relative position within the domain
     */
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
    // double x, y  relative positions
    /// the steepest gradient slope
    double fx, fy;

    /// Container of the adjacent edges
    std::vector<std::weak_ptr<Edge>> adj_edges;

    /// Container of the the adjacent cells
    std::vector<std::weak_ptr<Cell>> adj_cells;

    /// Constructor
    /** \param x, y     relative position within the domain
     */
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
        if (s->x >= 1.) { s->x -= 1.; }
        else if (s->x < 0.) { s->x += 1.; }
        if (s->y >= 1.) { s->y -= 1.; }
        else if (s->y < 0.) { s->y += 1.; }
    }

    return;
}
/// Absolute position
std::pair<double, double> absolute_position(const Site s,
                       const std::pair<double, double> domain_size)
{
    return std::make_pair(std::get<0>(domain_size)*s.x,
                          std::get<1>(domain_size)*s.y);
}

/// The relative displacement vector from a to b
template <bool periodic_bc>
std::pair<double, double> displacement(const Site &a, const Site &b) 
{
    double dx = b.x - a.x;
    double dy = b.y - a.y;

    if (dx >= 0.5) { dx -= 1.; }
    else if (dx < -0.5) { dx += 1.; }

    if (dy >= 0.5) { dy -= 1.; }
    else if (dy < -0.5) { dy += 1.; }

    #ifndef NDEBUG
    if constexpr (periodic_bc) {
        double tolerance = 0.01;
        if (dx > 0.5 - tolerance or dy > 0.5 - tolerance) {
            throw std::runtime_error("ERROR The displacement between Site a "
                    "and b is not unique, because of periodic boundary "
                    "conditions! a is (" + std::to_string(a.x) + ", " +
                    std::to_string(a.y) + ") and b is (" + std::to_string(b.x) +
                    ", " + std::to_string(b.y) + "). The vector b - a is (" + 
                    std::to_string(dx) + ", " + std::to_string(dy) + "). ");
        }
    }
    #endif

    return std::make_pair(dx, dy);
}

/// The absolute displacement vector from a to b
template <bool periodic_bc>
std::pair<double, double> displacement_absolute(const Site &a, const Site &b,
                                                double Lx, double Ly)
{
    auto [dx, dy] = displacement<periodic_bc>(a, b);

    return std::make_pair(Lx*dx, Ly*dy);
}

/// The absolute distance of two sites
template <bool periodic_bc>
double distance(const Site &a, const Site &b, double Lx, double Ly) {
    auto [dx, dy] = displacement_absolute<periodic_bc>(a, b, Lx, Ly);

    return std::max(sqrt(pow(dx, 2) + pow(dy, 2)), 1e-10);
}

/// Creates an object-copy, that is at true distance to another object
template <bool periodic_bc>
Site periodic_copy(const Site s, const Site s_fixed)
{
    if constexpr (not periodic_bc) {
        return Site(s.x, s.y);
    }

    auto [dx, dy] = displacement<periodic_bc>(s_fixed, s);

    return Site(s_fixed.x + dx, s_fixed.y + dy);
}

/// The Edge object
/** It is defined as the shortest connection between two vertices a and b
 * 
 *  Every edge has two adjoint cells, one to either side, except boundary edges
 *  have only one.
 */
struct Edge : public std::enable_shared_from_this<Edge> {
    /// Start and end vertices of this edge
    Vertex_ptr a, b;

    /// Absolute length of the edge
    /** \param Lx, Ly   Domain size 
     */
    template <bool periodic_bc>
    double length(double Lx, double Ly) {
        return distance<periodic_bc>(*a, *b, Lx, Ly);
    }

    /// The linetension parameter property to this edge
    double linetension;

    /// The contractility parameter
    double contractility;

    /// The first adjoint cell
    std::weak_ptr<Cell> adj_cell_a;

    /// The second adjoint cell
    std::weak_ptr<Cell> adj_cell_b;

    /// Polartity protein level on side to cell a
    double sigma_a;

    /// Polartity protein level on side to cell b
    double sigma_b;

    double get_sigma(Cell_ptr c) const {
        if (c == adj_cell_a.lock()) {
            return sigma_a;
        }
        else if (c == adj_cell_b.lock()) {
            return sigma_b;
        }
        else {
            throw std::runtime_error("The provided cell is not an adjoint cell"
                " of this edge!");
        }
    }

    void set_sigma(Cell_ptr c, double value) {
        if (c == adj_cell_a.lock()) {
            sigma_a = value;
        }
        else if (c == adj_cell_b.lock()) {
            sigma_b = value;
        }
        else {
            throw std::runtime_error("The provided cell is not an adjoint cell"
                " of this edge!");
        }
    }

    /// Steepest descent in protein level (side a)
    double d_sigma_a;

    /// Steepest descent in protein level (side b)
    double d_sigma_b;

    double get_d_sigma(Cell_ptr c) const {
        if (c == adj_cell_a.lock()) {
            return d_sigma_a;
        }
        else if (c == adj_cell_b.lock()) {
            return d_sigma_b;
        }
        else {
            throw std::runtime_error("The provided cell is not an adjoint cell"
                " of this edge!");
        }
    }

    void set_d_sigma(Cell_ptr c, double value) {
        if (c == adj_cell_a.lock()) {
            d_sigma_a = value;
        }
        else if (c == adj_cell_b.lock()) {
            d_sigma_b = value;
        }
        else {
            throw std::runtime_error("The provided cell is not an adjoint cell"
                " of this edge!");
        }
    }

    /// Whether this object is to be removed 
    bool remove;

    /// Constructor
    /** \param a    Start Vertex
     *  \param b    End Vertex
     *  \param linetension  The linetension property
     *  \param contractility    The contractility parameter
     *  \param adj_cs   The adjacent cells to this edge
     */
    Edge(Vertex_ptr a, Vertex_ptr b, double linetension, double contractility,
         Cell_ptr adj_cell_a = nullptr, Cell_ptr adj_cell_b = nullptr,
         double sigma_a = 0., double sigma_b = 0.)
    :
        a(a),
        b(b),
        linetension(linetension),
        contractility(contractility),
        adj_cell_a(adj_cell_a), adj_cell_b(adj_cell_b),
        sigma_a(sigma_a), sigma_b(sigma_b),
        d_sigma_a(0.), d_sigma_b(0.),
        remove(false)
    { }

    /// Set the weak pointer to this for members
    void link_members () {
        a->adj_edges.push_back(weak_from_this());
        b->adj_edges.push_back(weak_from_this());
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
    Site va = *e0.a;
    Site vc = periodic_copy<periodic_bc>(*e1.a, va);
    
    // the end vertices of the resp. edges
    Site vb = periodic_copy<periodic_bc>(*e0.b, va);
    Site vd = periodic_copy<periodic_bc>(*e1.b, vc);

    auto vector_e0 = displacement<periodic_bc>(va, vb);
    auto vector_e1 = displacement<periodic_bc>(vc, vd);

    // line defined by e0: y0 = a + bx
    double a, b;
    if (fabs(std::get<0>(vector_e0)) < 1e-8) {
        vb.x += 1e-8;
        vector_e0 = displacement<periodic_bc>(va, vb);
    }
    b = std::get<1>(vector_e0) / std::get<0>(vector_e0);
    a = va.y - b * va.x;

    // line defined by e1: y1 = c + dx
    double c, d;
    if (fabs(std::get<0>(vector_e1)) < 1e-8) {
        vd.x += 1e-8;
        vector_e1 = displacement<periodic_bc>(vc, vd);
    }
    d = std::get<1>(vector_e1) / std::get<0>(vector_e1);
    c = vc.y - d * vc.x;

    // parallel lines
    if (fabs(d - b) < 1e-8) {
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
struct Cell : public std::enable_shared_from_this<Cell> {
    /// The type of a cell
    enum CellType {
        progenitor,
        hair,
        support,
        num_cell_types,
    } type;

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
    
    /// The order of the edges
    /** If 1, then the edges are ordered anti-clockwise. Otherwise ordered 
     *  clockwise.
     * 
     *  NOTE updated together with area
     */
    char area_sgn;

    /// site of the cell
    /** NOTE access only centre_site! */
    Site_ptr site;

    /// Site at centre of cell
    template<bool periodic_bc>
    Site_ptr centre_site () {
        this->template area<periodic_bc>();
        return site;
    }

    /// The preferential area of the cell
    /** In absolute coordinates
     */
    double area_preferential;

    /// The contractility of the cell
    double contractility;

    /// Lagrange multiplier I
    /** Constraint of zero net polarisation
     */
    double lagrange_net_polarisation;

    /// Lagrange multiplier II
    /** Constraint of constant protein level
     */
    double lagrange_const_concentration;

    /// The initial protein concentration
    double protein_concentration;

    /// Whether this object is to be removed 
    bool remove;

    /// Constructor of a cell
    /** \param s    The center of the cell
     *  \param es   The edges defining the boundary of the cell in arbitrary
     *              order and direction
     *  \param area_preferential    The preferential size of this cell
     *  \param cell_type    The type of cell
     *  \param contractility        The contractility of the cell associated
     *                              with contractility of the actin-myosin ring
     */
    Cell(Site s, EdgeContainer es, double area_preferential, 
         double contractility, CellType cell_type = progenitor,
         double protein_concentration = 0.)
    :
        type(cell_type),
        vertices(),
        edges_ordered(),
        area_sgn(0),
        site(std::make_shared<Site>(s)),
        area_preferential(area_preferential),
        contractility(contractility),
        lagrange_net_polarisation(0.),
        lagrange_const_concentration(0.),
        protein_concentration(protein_concentration),
        remove(false)
    {
        order_edges(es);
    }

    /// Set the weak pointer to this for members
    void link_members () {
        for (auto v : vertices) {
            v->adj_cells.push_back(weak_from_this());
        }
        for (auto [e, flip] : edges_ordered) {
            if (e->adj_cell_a.expired()) {
                e->adj_cell_a = weak_from_this();
            }
            else if (e->adj_cell_b.expired()) {
                e->adj_cell_b = weak_from_this();
            }
        }
    }
    
    /// Reset the cell's boundary to new edges es
    /** Removes the current edges and vertices and sets them anew from es.
     *  The new edges are ordered.
     */
    void order_edges(EdgeContainer es) {
        edges_ordered.clear();
        vertices.clear();

        if (es.size() < 3) { 
            throw std::runtime_error("Cannot order egdes, because less "
                "than 3 received. No cell this is!"); }

        Edge_ptr e = es.back();
        Vertex_ptr a, b;

        // put b to the right of a
        double dx = e->b->x - e->a->x;
        if ((dx < 0 and dx > -0.5) or dx > 0.5) {
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
            
            // Breakpoint in debug mode
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
                    std::to_string(site->x) + ", " + std::to_string(site->y) + "). "
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

    /// Update the ordering of edges
    /** Edges are ordered to be clockwise or anti-clockwise
     */
    void order_edges() {
        EdgeContainer es;
        for (const auto e_pair : edges_ordered) {
            es.push_back(std::get<Edge_ptr>(e_pair));
        }
        order_edges(es);
    }

    /// The relative area of a polygon cell
    template <bool periodic_bc>
    double area () {
        // reset the area
        double area = 0.;

        // The new center of the cell
        Site center(0., 0.);

        // The starting vertex for iteration
        // periodic copies will be made wrt to this point
        Site start_vertex = *std::get<Edge_ptr>(edges_ordered.front())->a;
        if (std::get<bool>(edges_ordered.front())) {
            start_vertex = *std::get<Edge_ptr>(edges_ordered.front())->b;
        }

        for (const auto [e, flip] : edges_ordered) {            
            auto a = periodic_copy<periodic_bc>(*e->a, start_vertex);
            auto b = periodic_copy<periodic_bc>(*e->b, start_vertex);


            #ifndef NDEBUG
            displacement<periodic_bc>(a, b);
            #endif

            if (flip) { 
                std::swap(a, b);
            }

            double da = a.x * b.y - b.x * a.y;
            area += da;
            center.x += (a.x + b.x) * da;
            center.y += (a.y + b.y) * da;
        }

        if (area < 0) {
            // invert edge ordering
            auto tmp_es = edges_ordered;
            edges_ordered.clear();
            for (auto [e, flip] : tmp_es) {
                edges_ordered.insert(edges_ordered.begin(),
                                          std::make_pair(e, not flip));
            }

            area_sgn = 1;
        }
        else { area_sgn = 1; }

        area = 0.5 * fabs(area);

        center.x = area_sgn * center.x / 6 / area;
        center.y = area_sgn * center.y / 6 / area;
        
        // update the center site s
        site->x = center.x; site->y = center.y;
        correct_periodic_bc<periodic_bc>(site);

        return area;
    }

    /// The absolute area of the cell
    /** \param Lx, Ly   Domain size
     */
    template <bool periodic_bc>
    double area_abs(double Lx, double Ly) {
        return this->template area<periodic_bc>() * Lx * Ly;
    }
};

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif