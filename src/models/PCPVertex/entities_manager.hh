#ifndef UTOPIA_MODELS_PCPVERTEX_ENTITIESMANAGER_HH
#define UTOPIA_MODELS_PCPVERTEX_ENTITIESMANAGER_HH

#include <vector>
#include <numeric>

#include "../Collier/Collier.hh"
#include "../NotchDelta/NotchDelta.hh"

namespace Utopia::Models::PCPVertex {

template<class Model>
class EntitiesManager {

public:
    /// The type of the space
    using Space = typename Model::Space;

    /// The type of a vector in space
    using SpaceVec = typename Space::SpaceVec;

    /// The type of the configuration
    using Config = typename Model::Config;

    /// The traits of a Vertex
    /** Using no config for initialisation
     */
    using VertexTraits = Utopia::AgentTraits<VertexState, Update::manual, true>;

    /// The type of the managed vertices
    using Vertex = Utopia::Agent<VertexTraits, Space>;

    /// The links of an Edge
    template <class EdgeContainer>
    struct EdgeLinks {
        /// Start and end vertices of an Edge
        std::shared_ptr<Vertex> a, b;
    };

    /// The traits of an Edge
    using EdgeTraits = Utopia::AgentTraits<EdgeState, Update::manual, true,
                                           EmptyTag, EdgeLinks>;

    /// The type of the managed edges
    using Edge = Utopia::Agent<EdgeTraits, Space>;

    /// A container of edges associated with a flip boolean
    /** Every edge within the container has associated vertices a and b.
     *  The ordering is such that one can iterate the edges vertices a_i to b_i,
     *  where b_i the same vertex as a_{i+1} and b_N = a_0.
     *  If the edge's associated flip boolean is True, the ordering is possible
     *  with a_i and 
     */
    using OrderedEdgeContainer = std::vector<std::pair<std::shared_ptr<Edge>,
                                                       bool>>;

    /// The links of a Cell
    template <class CellContainer>
    struct CellLinks {
        /// The vertices forming this cell
        AgentContainer<Vertex> vertices;

        /// The edges forming the boundary of this cell
        /** \details The edges are maintained such that one can iterate from
         *           the first vertex back to itself in anti-clockwise manner
         *           using the edges two vertices.
         *           The boolian indicates whether the edge starts at
         *           vertex a or b. 
         */
        OrderedEdgeContainer edges;
        
        /// A cell in the NotchDelta model
        std::shared_ptr<Utopia::Models::NotchDelta::NotchDelta::Cell> nd_cell;
        
        /// A cell in the Collier model
        std::shared_ptr<Utopia::Models::Collier::Collier::Cell> c_cell;
    };
    using CellTraits = Utopia::AgentTraits<CellState, Update::manual, false,
                                           EmptyTag, CellLinks>;

    /// The type of the managed cells
    using Cell = Utopia::Agent<CellTraits, Space>;

    /// The type of a rule function acting on vertices of this agent manager
    /** This is a convenience type def that models can use to easily have this
      * type available.
      */
    using RuleFuncVertex = std::function<typename Vertex::State(
                                               const std::shared_ptr<Vertex>&)>;

    /// The type of a rule function acting on edges of this agent manager
    /** This is a convenience type def that models can use to easily have this
      * type available.
      */
    using RuleFuncEdge = std::function<typename Edge::State(
                                               const std::shared_ptr<Edge>&)>;

    /// The type of a rule function acting on cells of this agent manager
    /** This is a convenience type def that models can use to easily have this
      * type available.
      */
    using RuleFuncCell = std::function<typename Cell::State(
                                               const std::shared_ptr<Cell>&)>;

    /// The layout of the cellular agents
    enum Layout {
        /// A 2D apical view of a single layer of cells
        Apical,

        /// A 2D side view of a single layer of cells
        Columnar
    };


private:
    /// The logger (same as the model this manager resides in)
    const std::shared_ptr<spdlog::logger> _log;
    
    /// Agent manager configuration node
    const Config _cfg;

    /// The physical space the agents are to reside in
    const std::shared_ptr<Space> _space;

    const Layout _layout;

    /// The manager of the vertices
    Utopia::AgentManager<VertexTraits, Model> _vertex_manager;

    /// The manager of the edges
    Utopia::AgentManager<EdgeTraits, Model> _edge_manager;

    /// The manager of the cells
    Utopia::AgentManager<CellTraits, Model> _cell_manager;

    /// The container of adjoint cells to every cell
    std::unordered_map<int, std::pair<
            std::shared_ptr<Cell>, std::shared_ptr<Cell>>> _edges_adjoint_cells;

    /// The container of adjoint edges to every vertex
    std::unordered_map<int, AgentContainer<Edge>> _vertices_adjoint_edges;

    /// The container of adjoint edges to every vertex
    std::unordered_map<int, AgentContainer<Cell>> _vertices_adjoint_cells;

public:
    EntitiesManager (Model& model, const Config& custom_cfg = {})
    :
        _log(spdlog::stdout_color_mt(  model.get_logger()->name()
                                     + ".entities_manager")),
        _cfg(setup_cfg(model, custom_cfg)),
        _space(model.get_space()),
        _layout(setup_layout(_cfg)),
        _vertex_manager(model, setup_vertex_cfg()),
        _edge_manager(model, setup_edge_cfg()),
        _cell_manager(model, setup_cell_cfg())
    {
        // Set this model instance's log level
        if (_cfg["log_level"]) {
            // Via value given in configuration
            const auto lvl = get_as<std::string>("log_level", _cfg);
            _log->debug("Setting log level to '{}' ...", lvl);
            _log->set_level(spdlog::level::from_str(lvl));
        }
        else {
            // No config value given; use the level of the parent's logger
            _log->set_level(model.get_logger()->level());
        }

        setup_agents();
        _log->info("EntitiesManager is all set up.");
    }

    // -- Public interface ----------------------------------------------------
    /// Getter for space object pointer
    const auto& get_space() const {
        return _space;
    }
    
    /// Getter for the layout structure
    const auto& get_layout() const {
        return _layout;
    }

    /// Getter for the extent of the domain
    /** Returns the rectangle enclosing all vertices in non-periodic BC as the
     *  lower left and upper right corner positions.
     */
    const std::pair<SpaceVec, SpaceVec> get_extent() const {
        if (_space->periodic) {
            return std::make_pair(SpaceVec({0, 0}), _space->get_domain_size());
        }

        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::min();
        double y_min = std::numeric_limits<double>::max();
        double y_max = std::numeric_limits<double>::min();

        for (const auto& v : this->vertices()) {
            SpaceVec pos = this->position_of(v);
            x_min = std::min(x_min, pos[0]);
            x_max = std::max(x_max, pos[0]);
            y_min = std::min(x_min, pos[1]);
            y_max = std::max(x_max, pos[1]);
        }

        return std::make_pair(SpaceVec({x_min, y_min}),
                              SpaceVec({x_max, y_max}));
    }
    
    /// Return const reference to the managed vertices
    const auto& vertices () const {
        return _vertex_manager.agents();
    }

    auto vertices_id_counter () const {
        return _vertex_manager.id_counter();
    }

    /// Return const reference to the managed edges
    const auto& edges () const {
        return _edge_manager.agents();
    }

    auto edges_id_counter () const {
        return _edge_manager.id_counter();
    }

    /// Return const reference to the managed cells
    const auto& cells () const {
        return _cell_manager.agents();
    }

    auto cells_id_counter () const {
        return _cell_manager.id_counter();
    }

    SpaceVec position_of (const Vertex& vertex) const {
        return vertex.position();
    }

    SpaceVec position_of (const std::shared_ptr<Vertex>& vertex) const {
        return position_of(*vertex);
    }

    /// Move an vertex to a new position in the space
    void move_to(Vertex& vertex,
                 const SpaceVec& pos) const
    {
        return _vertex_manager.move_to(vertex, pos);
    }
    
    /// Move an vertex to a new position in the space
    void move_to(const std::shared_ptr<Vertex>& vertex,
                 const SpaceVec& pos) const
    {
        return move_to(*vertex, pos);
    }

    /// Move an vertex relative to its current position
    void move_by(Vertex& vertex,
                 const SpaceVec& move_vec) const
    {
        return move_to(vertex, position_of(vertex) + move_vec);
    }

    /// Move an vertex relative to its current position
    void move_by(const std::shared_ptr<Vertex>& vertex,
                 const SpaceVec& move_vec) const
    {
        return move_by(*vertex, move_vec);
    }

    /// The displacement between two vertices
    /** \details see Utopia::Space::displacement
     */
    SpaceVec displacement (const Vertex& a, const Vertex& b) const {
        return _space->displacement(position_of(a), position_of(b));
    }

    /// The displacement between two vertices
    /** \details see Utopia::Space::displacement
     */
    SpaceVec displacement (const std::shared_ptr<Vertex>& a,
                           const std::shared_ptr<Vertex>& b) const {
        return _space->displacement(position_of(a), position_of(b));
    }

    /// The displacement between the barycenter of two cells
    /** \details see Utopia::Space::displacement
     */
    SpaceVec displacement (const Cell& a, const Cell& b) const {
        return _space->displacement(barycenter_of(a), barycenter_of(b));
    }

    /// The displacement between the barycenter of two cells
    /** \details see Utopia::Space::displacement
     */
    SpaceVec displacement (const std::shared_ptr<Cell>& a,
                           const std::shared_ptr<Cell>& b) const {
        return _space->displacement(barycenter_of(a), barycenter_of(b));
    }

    /// The displacement between the barycenter defined by an edge
    /** \details see Utopia::Space::displacement
     */
    SpaceVec displacement (const std::shared_ptr<Edge>& e,
                           bool flip = false) const
    {
        if (not flip) {
            return displacement(e->custom_links().a, e->custom_links().b);
        }
        else {
            return displacement(e->custom_links().b, e->custom_links().a);
        }
    }

    /// The displacement between the barycenter defined by an edge
    /** \details see Utopia::Space::displacement
     */
    SpaceVec displacement (const std::pair<std::shared_ptr<Edge>,
                           bool> e_pair) const
    {
        auto [e, flip] = e_pair;
        return displacement(e, flip);
    }

    /// The distance between two vertices
    /** \details see Utopia::Space::distance
     */
    auto distance (const Vertex& a, const Vertex& b) const {
        return _space->distance(position_of(a), position_of(b));
    }

    /// The distance between two vertices
    /** \details see Utopia::Space::distance
     */
    auto distance (const std::shared_ptr<Vertex>& a,
                   const std::shared_ptr<Vertex>& b) const {
        return _space->distance(position_of(a), position_of(b));
    }

    /// The distance between two cells
    /** \details see Utopia::Space::distance
     */
    auto distance (const Cell& a, const Cell& b) const {
        return _space->distance(barycenter_of(a), barycenter_of(b));
    }

    /// The distance between two cells
    /** \details see Utopia::Space::distance
     */
    auto distance (const std::shared_ptr<Cell>& a,
                   const std::shared_ptr<Cell>& b) const {
        return _space->distance(barycenter_of(a), barycenter_of(b));
    }

    auto length_of (const std::shared_ptr<Edge>& edge) const {
        return distance(edge->custom_links().a, edge->custom_links().b);
    }

    SpaceVec center_of (const std::shared_ptr<Edge>& edge) const {
        SpaceVec a = position_of(edge->custom_links().a);
        SpaceVec displ = displacement(edge);
        return _space->map_into_space(a + .5 * displ);
    }

    /// Calculate the perimeter of a cell
    double perimeter_of (const std::shared_ptr<Cell>& cell) const
    {
        return perimeter_of(cell->custom_links().edges);
    }

    /// Calculate the perimeter of a boundary.
    /** This can be a boundary of a cell or any other closed loop of edges.
     */
    double perimeter_of (const OrderedEdgeContainer& boundary) const
    {
        double perimeter = 0.;
        for (const auto& [e, flip] : boundary) {
            perimeter += length_of(e);
        }
        return perimeter;
    }

    /// Calculate the area of a cell
    /** \tparam get_sign    If false the area can have a sign. It will be
     *                      negative when the cell's edges are ordered clockwise
     *                      instead of anti-clockwise.
     *                      If true, throws on counter-clockwise ordering
     * 
     *  \note Does not check whether the edges are ordered and in which
     *        orientation they are ordered.
     */
    template <bool get_sign = false>
    double area_of (const std::shared_ptr<Cell>& cell) const
    {
        return area_of<get_sign>(cell->custom_links().edges);
    }

    /// Calculate the area of a boundary
    /** \param boundary The edges defining a boundary. This can be of a cell
     *                  or any other closed loop.
     *  \tparam get_sign    If false the area can have a sign. It will be
     *                      negative when the cell's edges are ordered clockwise
     *                      instead of anti-clockwise.
     *                      If true, throws on counter-clockwise ordering
     * 
     *  \note Does not check whether the edges are ordered and in which
     *        orientation they are ordered.
     */
    template <bool get_sign = false>
    double area_of(const OrderedEdgeContainer& boundary) const 
    {
        static_assert(Space::dim == 2, "Area of a cell is only implemented for "
                      "2 dimensional space!");

        if (boundary.size() < 3) {
            throw std::runtime_error(fmt::format("Cannot calculate area for a "
                "boundary with {} < 3 edge(s)", boundary.size()));
        }

        // define a reference in space
        auto [e, flip] = boundary.front();
        SpaceVec ref;
        if (not flip) { ref = position_of(e->custom_links().a); }
        else { ref = position_of(e->custom_links().b); }

        double area = 0.;
        for (const auto& [e, flip] : boundary) {
            SpaceVec a = position_of(e->custom_links().a);
            SpaceVec b = position_of(e->custom_links().b);

            if (flip) { std::swap(a, b); }

            // define the vertices positions relative to the reference
            /* this is important in periodic space to calculate with 
                * "real" coordinates */
            a = ref + _space->displacement(ref, a);
            b = ref + _space->displacement(ref, b);
            ref = a;

            area += a[0] * b[1] - b[0] * a[1];
        }

        area /= 2.;

        if constexpr (not get_sign) {
            if (area < 1.e-12) {
                throw std::runtime_error(fmt::format(
                    "Negative area ({}) of a boundary defining a cell with "
                    "{} edges!", area, boundary.size()));
            }
        }

        return area;
    }

    /// Calculate the shape index of a cell
    /** shape index \f$ p = P / \sqrt(A) \f$ with the cell's perimeter \f$ P \f$
     *  and area \f$ A \f$.
     */ 
    double shape_index_of (const std::shared_ptr<Cell>& cell) const
    {
        return shape_index_of(cell->custom_links().edges);
    }

    /// Calculate the shape index of a boundary
    /** shape index \f$ p = P / \sqrt(A) \f$ with the cell's perimeter \f$ P \f$
     *  and area \f$ A \f$.
     */ 
    double shape_index_of (const OrderedEdgeContainer& boundary) const
    {
        return perimeter_of(boundary) / sqrt( area_of(boundary) );
    }
    
    /// Returns the barycenter of the given cell
    /** \note   It is assumed that the edges are ordered anti-clockwise or 
     *          clock-wise
     * 
     *  \param cell     The considered cell
     */
    SpaceVec barycenter_of (const std::shared_ptr<Cell>& cell) const {
        return barycenter_of(cell->custom_links().edges);
    }
    
    /// Returns the barycenter of the given boundary
    /** \note   It is assumed that the edges are ordered anti-clockwise or 
     *          clock-wise
     * 
     *  \param cell     The considered cell
     */
    SpaceVec barycenter_of (const OrderedEdgeContainer& boundary) const {
        static_assert(Space::dim == 2, "Center of a cell is only implemented "
                      "for 2 dimensional space!");

        if (boundary.size() < 3) {
            throw std::runtime_error("Cannot calculate barycenter of an"
                "edge collection with less than 3 edges!");
        }

        // define a reference in space
        auto [e, flip] = boundary.front();
        SpaceVec ref;
        if (not flip) { ref = position_of(e->custom_links().a); }
        else { ref = position_of(e->custom_links().b); }

        double area = 0.;
        SpaceVec center(arma::fill::zeros);
        for (const auto& [e, flip] : boundary) {
            SpaceVec a = position_of(e->custom_links().a);
            SpaceVec b = position_of(e->custom_links().b);

            if (flip) { std::swap(a, b); }

            // define the vertices positions relative to the reference
            /* this is important in periodic space to calculate with 
                * "real" coordinates */
            a = ref + _space->displacement(ref, a);
            b = ref + _space->displacement(ref, b);
            ref = a; 

            double da = a[0] * b[1] - b[0] * a[1];
            area += da;
            center += (a + b) * da;
        }

        area /= 2;
        if (area < 1.e-12) { area += 1e-12; }
        center /= (6 * area);

        return _space->map_into_space(center);
    }

    /// The adjoint edges of a vertex
    /** The container of edges starting or ending in this vertex
     * 
     *  \param vertex   The considered vertex
     */
    AgentContainer<Edge> adjoint_edges_of(
            const std::shared_ptr<Vertex>& vertex) const
    {
        return _vertices_adjoint_edges.at(vertex->id());
    }

    /// The adjoint cells of a vertex
    /** The container of cells that include this vertex in their boundary
     * 
     *  \param vertex   The considered vertex
     */
    AgentContainer<Cell> adjoint_cells_of(
            const std::shared_ptr<Vertex>& vertex) const
    {
        return _vertices_adjoint_cells.at(vertex->id());
    }

    /// The adjoint cells of an edge
    /** The cells on either side of this edge, i.e. those that include the edge
     *  in their boundary. If aligned, the first cell is the cell on the left
     *  of the edge (vertices a->b). In this cell this edge is oriented
     *  anti-clockwise.
     * 
     *  \param edge   The considered edge
     */
    template<bool aligned = false>
    std::pair<std::shared_ptr<Cell>, std::shared_ptr<Cell>> adjoints_of(
            const std::shared_ptr<Edge>& edge) const
    {
        if constexpr (aligned) {
            auto [a, b] = _edges_adjoint_cells.at(edge->id());
            SpaceVec _a, _b;
            if (a) { _a = barycenter_of(a); }
            else   { _a = center_of(edge); }
            if (b) { _b = barycenter_of(b); }
            else   { _b = center_of(edge); }
            
            // is a left of edge: a->b clockwise wrt edge?
            SpaceVec displ = this->displacement(edge);
            SpaceVec connector = _space->displacement(_a, _b);
            if ((displ[0] * connector[1] - displ[1] * connector[0]) > 0) {
                std::swap(a, b);
            }
            return std::make_pair(a, b);
        }
        else {
            return _edges_adjoint_cells.at(edge->id());
        }
    }

    /// The neighboring cells to a cell
    /** The cells that are separated from this cell by a common edge
     * 
     *  \param cell   The considered cell
     */
    AgentContainer<Cell> neighbors_of(const std::shared_ptr<Cell>& cell) const
    {
        AgentContainer<Cell> neighbors;

        for (const auto& [e, flip] : cell->custom_links().edges) {
            const auto [adj_cell_a, adj_cell_b] = adjoints_of(e);
            if (adj_cell_a and adj_cell_a != cell) {
                neighbors.push_back(adj_cell_a);
            }
            else if (adj_cell_b and adj_cell_b != cell) {
                neighbors.push_back(adj_cell_b);
            }
        }

        return neighbors;
    }

    /// The neighbors of a cell given a Mannhatten-distance
    /** The larger neighborhood of a cell.
     * 
     *  For distance = 1, returns neighbors_of(cell). Otherwise, iterates the 
     *  neighbours and includes their neighbors too. Repeat N = distance times.
     * 
     *  \param cell     The considered cell
     *  \param distance The manhatten-distance to cells to include in
     *                  neighborhood. For distance = 1, equiv to 
     *                  neighbors_of(cell). Otherwise, iterates the 
     *                  neighbours and includes their neighbors too. 
     *                  Repeat N = distance times.
     *  \tparam check_confluence    Whether to check for confluence of tissue.
     *                              If false, for distance larger the number 
     *                              of cells in the tissue, returns all cells.
     *                              In a non-confluent tissue, this could be a 
     *                              subset of all cells.
     *  
     *  \returns A container with all neighbors of cell. This does not include 
     *      cell itself. For distance=0, this is an empty container.
     */
    template<bool check_confluence=false>
    AgentContainer<Cell> neighbors_of(
            const std::shared_ptr<Cell>& cell,
            const std::size_t distance
    ) const
    {
        if (distance == 1) {
            return neighbors_of(cell);
        }
        if (distance == 0) {
            return AgentContainer<Cell>({});
        }
        if constexpr (not check_confluence) {
            if (distance >= cells().size()) {
                return cells();
            }
        }

        std::unordered_set<std::shared_ptr<Cell>> neighbors({cell});
        for (std::size_t d = 0; d < distance; d++) {
            AgentContainer<Cell> tmp(neighbors.begin(), neighbors.end());
            for (const auto& nb : tmp) {
                auto nnbs = neighbors_of(nb);
                for (const auto& nnb : nnbs) {
                    neighbors.insert(nnb);
                }
            }
        }

        neighbors.erase(cell);
        neighbors.erase(nullptr);
        
        return AgentContainer<Cell>(neighbors.begin(), neighbors.end());
    }

    /// The neighboring cells to a cell of type `hair`
    /** \param cell   The considered cell
     */
    AgentContainer<Cell> neighbors_of_type(
            const std::shared_ptr<Cell>& cell,
            const std::size_t& type,
            const std::size_t& distance=1
    ) const
    {
        std::set<std::shared_ptr<Cell>> subset_neighbors;

        auto neighbors = this->neighbors_of(cell, distance);

        for (const auto & nb : neighbors) {
            if (nb->state.type == type)
            {
                subset_neighbors.insert(nb);
            }
        }
        subset_neighbors.erase(nullptr);

        return AgentContainer<Cell>(subset_neighbors.begin(),
                                    subset_neighbors.end());
    }

    /// Hexatic order of a cell
    /** Returns the hexatic and the corrected hexatic order. 
     *  Correction rescales the position of next neighbour HC from an elongated
     *  arrangement to a circular configuration.
     */
    std::pair<double, double> hexatic_order_of (
        const std::shared_ptr<Cell>& cell,
            const std::size_t& type,
            const std::size_t& distance=1
    ) const 
    {
        const auto neighbors = neighbors_of_type(cell, type, distance);
        
        if (neighbors.size() < 3) {
            return std::make_pair(0., 0.);
        }
        
        SpaceVec center = this->barycenter_of(cell);

        // calculate elongation and rotation of that object
        arma::mat22 q = this->elongation_of(cell);
        double abs_q = sqrt(arma::trace(q * q.t())/2.);

        // the rotation angle
        double phi = atan2(q[2], q[0]) / 2.;
        // the ratio of long and short axis of the ellipse
        double ratio_WH = exp(2 * abs_q);


        // calculate corrected hexatic order
        // NOTE a sheared hexagon maps to 1
        using namespace std::complex_literals;
        std::complex<double> hex_order = std::accumulate(
            neighbors.begin(), neighbors.end(),
            std::complex<double>(0., 0.),
            [this, center]
            (std::complex<double> val, const auto& nb)
            {
                using namespace std::complex_literals;

                SpaceVec pos = this->barycenter_of(nb);
                SpaceVec displ = this->_space->displacement(center, pos);

                return val + std::exp(1i * 6. * atan2(displ[1], displ[0]));
            }
        );
        std::complex<double> hex_order_corr = std::accumulate(
            neighbors.begin(), neighbors.end(),
            std::complex<double>(0., 0.),
            [this, center, phi, ratio_WH]
            (std::complex<double> val, const auto& nb)
            {
                using namespace std::complex_literals;

                SpaceVec pos = this->barycenter_of(nb);
                SpaceVec displ = this->_space->displacement(center, pos);
                double dx_ = (  (  displ[0] * cos(-phi) - displ[1] * sin(-phi))
                              / ratio_WH);
                double dy_ =  displ[0] * sin(-phi) + displ[1] * cos(-phi);
                double dx = dx_ * cos(phi) - dy_ * sin(phi);
                double dy = dx_ * sin(phi) + dy_ * cos(phi);

                return val + std::exp(1i * 6. * atan2(dy, dx));
            }
        );

        std::complex<double> N(neighbors.size());
        return std::make_pair(std::norm(hex_order / N),
                              std::norm(hex_order_corr / N));
    }

    // see transitions.hh
    std::pair<std::shared_ptr<Cell>, std::shared_ptr<Cell> >
    divide_cell(const std::shared_ptr<Cell> cell, double division_angle);

    // see transitions.hh
    bool remove_edge_T1 (const std::shared_ptr<Edge> edge,
        std::function<double(const AgentContainer<Vertex>&,
                             const AgentContainer<Edge>&,
                             const AgentContainer<Cell>&)> get_energy,
        double separation, double T1_barrier, double random_number,
        std::size_t time);

    // see transitions.hh
    bool remove_cell_T2 (const std::shared_ptr<Cell> cell);


    // -- Identifiers for objects at a non-periodic boundary ------------------
    
    /// Whether the vertex is a boundary vertex with 3 adjoint edges
    /** Such cells have only 2 adjoint cells, but 3 adjoint edges
     * 
     *  \param vertex   The considered vertex
     */
    bool is_3_fold_boundary_vertex(const std::shared_ptr<Vertex>& vertex) const
    {
        if (adjoint_cells_of(vertex).size() == 2) {
            return true;
        }
        else {
            return false;
        }
    }

    /// Whether the vertex is a boundary vertex with 2 adjoint edges
    /** Such cells have only 1 adjoint cells and 2 adjoint edges
     * 
     *  \param vertex   The considered vertex
     */
    bool is_2_fold_boundary_vertex(const std::shared_ptr<Vertex>& vertex) const
    {
        if (adjoint_cells_of(vertex).size() == 1) {
            return true;
        }
        else {
            return false;
        }
    }
    
    /// Whether a vertex is part of the tissue boundary
    /** This can either be a is_3_fold_boundary_vertex() or
     *  a is_2_fold_boundary_vertex() vertex.
     */
    bool is_boundary(const std::shared_ptr<Vertex>& vertex) const {
        if (   is_3_fold_boundary_vertex(vertex)
            or is_2_fold_boundary_vertex(vertex))
        {
            return true;
        }
        else {
            return false;
        }
    }

    /// Whether an edge ends at a vertex part of the tissue boundary
    /** Such an edge starts (ends) at a vertex that is part of the boundary
     * 
     *  \tparam safe    Whether to distinguish from is_1_cell_boundary_edge().
     *                  A 1 cell boundary edge fullfills the condition, but is
     *                  considered different. If false, the difference is not
     *                  made.
     *  \param edge     The considered edge
     */
    template<bool safe = true>
    bool is_2_cell_boundary_edge(const std::shared_ptr<Edge>& edge) const {
        // exclude the is_1_cell_boundary_edge() edges.
        if constexpr (safe) {
            if (is_1_cell_boundary_edge(edge)) {
                return false;
            }
        }

        // is one vertex part of boundary?
        for (const auto& v : {edge->custom_links().a, edge->custom_links().b}) {
            if (is_boundary(v)) {
                return true;
            }
        }

        // bulk edge
        return false;
    }

    /// Whether an edge is the interface of a cell with the void
    /** Such an edge has only 1 adjoint cell. Both adjoint vertices are boundary
     *  vertices. The adjoint vertices can be 2 or 3 fold each.
     */
    bool is_1_cell_boundary_edge(const std::shared_ptr<Edge>& edge) const {
        const auto [a, b] = adjoints_of(edge);
        if (not a or not b) {
            return true;
        }
        else {
            return false;
        }
    }

    /// Whether the edge is part of the tissue boundary
    /** This can be a is_2_cell_boundary_edge() or is_1_cell_boundary_edge()
     *  edge.
     */
    bool is_boundary(const std::shared_ptr<Edge>& edge) const {
        if (   is_1_cell_boundary_edge(edge)
            or is_2_cell_boundary_edge<false>(edge))
        {
            return true;
        }
        else {
            return false;
        }
    }

    /// Whether a cell is part of the tissue boundary
    /** Such a cell has at least one vertex that is part of the boundary.
     */
    bool is_boundary(const std::shared_ptr<Cell>& cell) const {
        for (const auto& v : cell->custom_links().vertices) {
            if (is_boundary(v)) {
                return true;
            }
        }
        return false;
    }

    OrderedEdgeContainer get_boundary_edges () const {
        if (_space->periodic) {
            return OrderedEdgeContainer{};
        }

        const auto& edges = this->edges();

        OrderedEdgeContainer boundary{};
        boundary.reserve(edges.size());

        auto edge = *std::find_if(
            edges.begin(), edges.end(),
            [this](const auto& edge) {
                return (   this->is_1_cell_boundary_edge(edge)
                        and not edge->state.remove);
                /* WARN this is only a work for a unprecise removal in
                 *      remove_boundary_edge!
                 *      get_boundary_edges will be called before the
                 *      remove_boundary_edge terminates. At that point the 
                 *      edge has not been removed from the edges container,
                 *      but all weak links are set correctly.
                 *      Hence, once a true edge has been picked, iteration 
                 *      works fine!
                 */
            });

        const auto start = edge->custom_links().a;
        auto iter = edge->custom_links().b;
        boundary.push_back(std::make_pair(edge, false));

        while (iter != start) {
            const auto tmp_edges = adjoint_edges_of(iter);
            auto edge_it = std::find_if(
                tmp_edges.begin(), tmp_edges.end(),
                [this, edge](const auto& e_it) {
                    if (e_it == edge) {
                        return false;
                    }
                    return this->is_1_cell_boundary_edge(e_it);
                });
            if (edge_it == tmp_edges.end()) {
                throw std::runtime_error("No 1 cell boundary edge found to "
                    "continue iteration of boundary!");
            }
            edge = *edge_it;

            if (iter == edge->custom_links().a) {
                boundary.push_back(std::make_pair(edge, false));
                iter = edge->custom_links().b;
            }
            else if (iter == edge->custom_links().b) {
                boundary.push_back(std::make_pair(edge, true));
                iter = edge->custom_links().a;
            }
            else {
                throw std::runtime_error(fmt::format("Failed to iterate "
                    "boundary with edge {} ({} -> {}) not starting at vertex {}",
                    edge->id(),
                    edge->custom_links().a->id(), edge->custom_links().b->id(),
                    iter->id()));
            }

            if (boundary.size() > edges.size()) {
                throw std::runtime_error("Failed to determine tissue boundary!");
            }
        }

        double area = area_of<true>(boundary);

        if (area < 0.) {
            boundary = this->reverse_edge_ordering(boundary);
            area = fabs(area);
        }

        boundary.shrink_to_fit();

        return boundary;
    }

    std::pair<OrderedEdgeContainer,
              OrderedEdgeContainer> get_apical_basal_boundary_edges () const
    {
        if (_layout != Columnar) {
            return std::make_pair(OrderedEdgeContainer{},
                                  OrderedEdgeContainer{});
        }
        if (not _space->periodic) {
            throw std::runtime_error("Apical boundary edges not defined in "
                                     "non-periodic space!");
        }

        // A function to gather connected boundary edges starting at edge
        std::function<OrderedEdgeContainer(std::shared_ptr<Edge>)>
        gather_boundary = [this](auto edge) {
            OrderedEdgeContainer boundary{};

            auto start = edge->custom_links().a;
            auto iter = edge->custom_links().b;

            // determine orientation of edge in adjoint cell
            auto [adj_cell, _adj_cell] = this->adjoints_of(edge);
            if (not adj_cell) {
                std::swap(adj_cell, _adj_cell);
            }
            auto [_e, flip] = *std::find_if(
                adj_cell->custom_links().edges.begin(),
                adj_cell->custom_links().edges.end(),
                [edge](const auto& e_pair) {
                    return std::get<0>(e_pair) == edge;
                });
            if (flip) {
                std::swap(start, iter);
            }
            boundary.push_back(std::make_pair(edge, flip));

            while (iter != start) {
                const auto tmp_edges = adjoint_edges_of(iter);
                auto edge_it = std::find_if(
                    tmp_edges.begin(), tmp_edges.end(),
                    [this, edge](const auto& e_it) {
                        if (e_it == edge) {
                            return false;
                        }
                        return this->is_1_cell_boundary_edge(e_it);
                    });
                if (edge_it == tmp_edges.end()) {
                    throw std::runtime_error("No 1 cell boundary edge found to "
                        "continue iteration of boundary!");
                }
                edge = *edge_it;

                if (iter == edge->custom_links().a) {
                    boundary.push_back(std::make_pair(edge, false));
                    iter = edge->custom_links().b;
                }
                else if (iter == edge->custom_links().b) {
                    boundary.push_back(std::make_pair(edge, true));
                    iter = edge->custom_links().a;
                }
                else {
                    throw std::runtime_error(fmt::format("Failed to iterate "
                        "boundary with edge {} ({} -> {}) not starting at vertex {}",
                        edge->id(),
                        edge->custom_links().a->id(),
                        edge->custom_links().b->id(),
                        iter->id()));
                }

                if (boundary.size() > edges().size()) {
                    throw std::runtime_error("Failed to determine tissue boundary!");
                }
            }

            return boundary;
        };

        auto edges = this->edges();

        // get boundary for one boundary edge
        auto edge = *std::find_if(
            edges.begin(), edges.end(),
            [this](const auto& edge) {
                return this->is_1_cell_boundary_edge(edge);
            });
        auto boundary_0 = gather_boundary(edge);

        // get the boundary for the other boundary
        auto edge_it = std::find_if(
            edges.begin(), edges.end(),
            [this, boundary_0](const auto& edge) {
                return (    this->is_1_cell_boundary_edge(edge)
                        and std::find_if(boundary_0.begin(), boundary_0.end(),
                                         [edge](const auto& e_pair) {
                                             return std::get<0>(e_pair) == edge;
                                         }) == boundary_0.end()
                       );
            });
        if (edge_it == edges.end()) {
            throw std::runtime_error("No 1 cell boundary edge found to "
                                     "determine second boundary!");
        }
        edge = *edge_it;

        auto boundary_1 = gather_boundary(edge);

        // determine which boundary points left, which right
        double sign = std::accumulate(
            boundary_0.begin(), boundary_0.end(), 0.,
            [this](double val, const auto& e_pair) -> double {
                return val + this->displacement(e_pair)[0];
            });
        
        // apical boundary points left (sign < 0)
        if (sign < 0) {
            return std::make_pair(boundary_0, boundary_1);
        }
        else {
            return std::make_pair(boundary_1, boundary_0);
        }
    }


    // -- Dual lattice --------------------------------------------------------

public:
    using Triangle = std::tuple<SpaceVec, SpaceVec, SpaceVec>;

    std::tuple<double, arma::mat22, arma::mat22> decompose
    (const arma::mat22& mat) const {
        double tr = arma::trace(mat);
        arma::mat22 sym =  (0.5 * (mat + arma::trans(mat))
                            - tr / 2. * arma::eye(2, 2));
        arma::mat22 asym = 0.5 * (mat - arma::trans(mat));

        return std::make_tuple(tr, sym, asym);
    }

    arma::mat22 equilateral_triangle () const {
        // the basis vectors of a equilateral triangle
        SpaceVec e0({1., 0.});              // a -> b
        SpaceVec e1({-0.5, sqrt(3) / 2.});  // b -> c

        arma::mat22 mat({{e0[0], e1[0]},
                         {e0[1], e1[1]}});

        return mat;
        // NOTE if change, use upper triangular form
    }

    /// Express triangle in basis of equilateral triangle
    /** If safe, c is set to be upper most vertex and a and b are ordered 
     *  anti-clockwise. Otherwise, this is assumed
     */
    template<bool ordered = false, bool c_upper_vertex = false>
    arma::mat22 shape (const Triangle& T) const {
        SpaceVec a, b, c;
        std::tie(a, b, c) = T;

        if constexpr (not ordered) {
            double area = 0.5 * (  a[0] * b[1] - b[0] * a[1]
                                 + b[0] * c[1] - c[0] * b[1]
                                 + c[0] * a[1] - a[0] * c[1]);

            // is the ordering clockwise ?
            if (area < 1.e-12) {
                // reorder to anti-clockwise
                std::swap(b, c);
            }
        }

        if constexpr (not c_upper_vertex) {
            // set c to be upper vertex
            if (a[1] > c[1] and a[1] > b[1]) {
                // a is uppermost -> rotate anti-clokwise
                std::swap(a, c);
                std::swap(a, b);
            }
            else if (b[1] > c[1] and b[1] > a[1]) {
                // b is uppermost -> rotate clockwise
                std::swap(a, c);
                std::swap(b, c);
            }
            // now c is uppermost
        }

        SpaceVec v0 = b - a;
        SpaceVec v1 = c - b;
        arma::mat22 base_1({{v0[0], v1[0]},
                            {v0[1], v1[1]}});
        
        arma::mat22 base_0 = equilateral_triangle();
        return base_1 * base_0.i();
    }

    std::tuple<double, arma::mat22, double, double> interpretation
    (const arma::mat22& shape) const
    {
        auto [tr, sym, asym] = decompose(shape);

        double area = arma::det(shape);
        [[maybe_unused]] double area0 = sqrt(3) / 4;
        // NOTE area := area / area0

        arma::mat22 h = asym + arma::eye(2, 2) * tr / 2.;

        arma::mat22 rot = h / sqrt(arma::trace(h * h.t()) / 2.);
        double theta = acos(rot[0]);
        if (rot[1] < 0.) {
            theta = -theta;
        }

        double s = sqrt(arma::trace(sym * sym.t()) / 2.);

        arma::mat22 q;
        if (s < 1.e-12) {
            q = sym * rot.t();
        }
        else {
            q = 1 / s * asinh(s / sqrt(area)) * (sym * rot.t());
        }

        double r = exp(2 *  sqrt(arma::trace(q * q.t()) / 2.));

        return std::make_tuple(area, q, r, theta);
    }

    /** If safe, c is set to be upper most vertex and a and b are ordered 
     *  anti-clockwise. Otherwise, this is assumed.
     *  If map_into_space, the coordinates lie within the space, otherwise
     *  might lie outside (periodic boundaries).
     */
    template<bool map_into_space = true, bool order = true>
    Triangle dual_triangle (const std::shared_ptr<Vertex>& vertex) const
    {
        auto adjoints = adjoint_cells_of(vertex);
        
        if (adjoints.size() != 3) {
            throw std::runtime_error(
                "Cannot define dual triangle on boundary vertex!");
        }

        // the dual vertices
        SpaceVec a, b, c;
        a = this->barycenter_of(adjoints[0]);
        b = a + _space->displacement(a, this->barycenter_of(adjoints[1]));
        c = a + _space->displacement(a, this->barycenter_of(adjoints[2]));

        if constexpr (order) {
            double area = 0.5 * (  a[0] * b[1] - b[0] * a[1]
                                 + b[0] * c[1] - c[0] * b[1]
                                 + c[0] * a[1] - a[0] * c[1]);

            // is the ordering clockwise ?
            if (area < 1.e-12) {
                // reorder to anti-clockwise
                std::swap(b, c);
            }
        }
        
        if constexpr (map_into_space) {
            b = _space->map_into_space(b);
            c = _space->map_into_space(c);
        }

        return std::make_tuple(a, b, c);
    }

    /// Calculate the elongation tensor of the object
    /** \param corners      The corners of a polygon
     *  \param center       The center of the polygon
     *  \tparam sorted      Whether the corners are sorted anti-clockwise
     */
    template<bool sorted = false>
    arma::mat22 elongation_of (
        std::vector<SpaceVec> corners,
        SpaceVec center = SpaceVec({0., 0.})
    ) const 
    {
        if (arma::norm(center) > 1.e-8) {
            std::for_each(corners.begin(), corners.end(),
                          [this, center](const SpaceVec& pos) {
                                return this->_space->displacement(center, pos);
                          });
        }
        
        if constexpr (not sorted) {
            std::sort(corners.begin(), corners.end(), 
                     [](const SpaceVec& a, const SpaceVec& b) {
                            return atan2(a[1], a[0]) < atan2(b[1], b[0]);
                     }
            );
        }


        arma::mat22 shape = arma::zeros(2, 2);
        double dual_area = 0.;
        for (std::size_t i = 0; i < corners.size(); i++) {
            std::size_t j = (i + 1) % corners.size();
            arma::mat22 s = this->shape<true, true>(std::make_tuple(
                SpaceVec({0., 0.}), corners[i], corners[j]));

            auto [area, q, r, theta] = interpretation(s);

            shape += area * q;
            dual_area += area;
        }
        
        return shape / dual_area;
    }

    template <bool dual_lattice=true>
    arma::mat22 elongation_of (const std::shared_ptr<Cell>& cell) const {
        if (this->is_boundary(cell)) {
            arma::mat22 q;
            q.fill(arma::datum::nan);
            return q;
        }

        const auto& edges = cell->custom_links().edges;
        SpaceVec center = barycenter_of(cell);

        std::vector<SpaceVec> corners{};
        if constexpr (dual_lattice) {
            for (const auto& [edge, flip] : edges) {
                const auto& [adj_a, adj_b] = adjoints_of(edge);
                if (adj_a != cell) {
                    corners.push_back(
                        this->_space->displacement(center,
                                                   this->barycenter_of(adj_a)));
                }
                else {
                    corners.push_back(
                        this->_space->displacement(center,
                                                   this->barycenter_of(adj_b)));
                }
            }
        }
        else {
            for (const auto& [edge, flip] : edges) {
                auto a = edge->custom_links().a;
                auto b = edge->custom_links().b;
                if (flip) {
                    std::swap(a, b);
                }

                SpaceVec pos = position_of(a);
                corners.push_back(_space->displacement(pos, center));
            }
        }
        
        return elongation_of<true>(corners);
    }

private:
    // -- Setup functions -----------------------------------------------------
    
    /// Set up the custom agent manager configuration member
    /** \details This function determines whether to use a custom configuration
      *         or the one provided by the model this AgentManager belongs to
      */
    Config setup_cfg(Model& model, const Config& custom_cfg) {
        Config cfg;

        if (custom_cfg.size() > 0) {
            _log->debug("Using custom config for agent manager setup ...");
            cfg = custom_cfg;
        }
        else {
            _log->debug("Using '{}' model's configuration for agent manager "
                        "setup ... ", model.get_name());

            if (not model.get_cfg()["agent_manager"]) {
                throw std::invalid_argument("Missing config entry "
                    "'agent_manager' in model configuration! Either specify "
                    "that key or pass a custom configuration node to the "
                    "AgentManager constructor.");
            }
            cfg = model.get_cfg()["agent_manager"];
        }

        return cfg;
    }

    /// Setup the layout from cfg
    Layout setup_layout(const Config& cfg) {        
        auto method = get_as<std::string>("setup_method", cfg);

        // call the requested method
        if (method == "hexagonal") {
            return Apical;
        }
        else if (method == "single") {
            return Apical;
        }
        else if (method == "column") {
            return Columnar;
        }
        // method not found!
        else {
            throw KeyError(method, cfg, "Cannot extract layout of agents! "
                           "Please choose one of the following setup "
                           "methods: 'hexagonal', 'single', 'column'.");
        }
    }

    /// Setup vertex manager
    /** \details Connot use default agent constructor if initializing hexagon.
     *           Hence, set up agent manager with 0 agents and add manually
     */
    Config setup_vertex_cfg () {
        this->_log->debug("Setting up vertex manager ..");
        if (not _cfg["vertex_manager"]) {
            throw KeyError("vertex_manager", _cfg, "In 'agent_manager'");
        }
        Config cfg = _cfg["vertex_manager"];
        if (get_as<std::string>("setup_method", _cfg) != "") {
            cfg["initial_num_agents"] = 0;
        }
        return cfg;
    }

    /// Setup edge manager
    /** \details Connot use default agent constructor if initializing hexagon.
     *           Hence, set up agent manager with 0 agents and add manually
     */
    Config setup_edge_cfg () {
        this->_log->debug("Setting up edge manager ..");
        if (not _cfg["edge_manager"]) {
            throw KeyError("edge_manager", _cfg, "In 'agent_manager'");
        }
        Config cfg = _cfg["edge_manager"];
        if (get_as<std::string>("setup_method", _cfg) != "") {
            cfg["initial_num_agents"] = 0;
        }

        return cfg;
    }

    /// Setup vertex manager
    /** \details Connot use default agent constructor if initializing hexagon.
     *           Hence, set up agent manager with 0 agents and add manually
     */
    Config setup_cell_cfg () {
        this->_log->debug("Setting up cell manager ..");
        if (not _cfg["cell_manager"]) {
            throw KeyError("cell_manager", _cfg, "In 'agent_manager'");
        }

        Config cfg = _cfg["cell_manager"];
        if (get_as<std::string>("setup_method", _cfg) != "") {
            cfg["initial_num_agents"] = 0;
        }

        return cfg;
    }

    /// Setup the agents following a specified method
    /** Extracts the `setup_method` from _cfg
     * 
     *  The `setup_method` can be:
     *      - `hexagonal`: Initialize the agents in a hexagonal lattice using
     *                     setup_agents_hexagonal_structure().
     *      - `single`: Initialize a single hexagon using
     *                  setup_agents_single_hexagon().
     * 
     *  The parameters for the setup method have to be passed in
     *  `setup_params/<method>`.
     */
    void setup_agents() {
        // extract the method of how to arrange the initial agents
        auto method = get_as<std::string>("setup_method", _cfg);

        // check that the appropriate params have been passed
        if (not _cfg["setup_params"]) {
            throw KeyError ("setup_params", _cfg, "No parameters provided "
                            "for the setup of the agent manager!");
        }
        if (not _cfg["setup_params"][method]) {
            throw KeyError (method, _cfg["setup_params"], "No parameters "
                            "provided for the setup of agent manager with the "
                            "specified method!");
        }
        
        // call the requested method
        if (method == "hexagonal") {
            this->setup_agents_hexagonal_structure(
                    _cfg["setup_params"]["hexagonal"]);
        }
        else if (method == "single") {
            this->setup_agents_single_hexagon(_cfg["setup_params"]["single"]);
        }
        else if (method == "column") {
            this->setup_agents_column(_cfg["setup_params"]["column"]);
        }
        // method not found!
        else {
            throw KeyError(method, _cfg, "Method not implemented to setup agent "
                           "manager! Please choose one of the following setup "
                           "methods: 'hexagonal', 'single', 'column'.");
        }
    }

    // see initialisation.hh
    void setup_agents_hexagonal_structure(const Config& cfg);
    void setup_agents_column(const Config& cfg);

    /// Setup a single cell of hexagonal shape in center of space
    /** Requires the config entry `size`, the length of an edge of the hexagon.
     */
    void setup_agents_single_hexagon(const Config& cfg) {
        double size = get_as<double>("size", cfg);   
        double width = sqrt(3) * size;
        double height = 2 * size;

        SpaceVec pos = _space->extent / 2.;

        add_vertex(pos + SpaceVec({0., height/2.}));
        add_vertex(pos + SpaceVec({width/2., height/4.}));
        add_vertex(pos + SpaceVec({width/2., -height/4.}));
        add_vertex(pos + SpaceVec({0., -height/2.}));
        add_vertex(pos + SpaceVec({-width/2., -height/4.}));
        add_vertex(pos + SpaceVec({-width/2., height/4.}));

        auto vertices = this->vertices();
        add_edge(vertices[0], vertices[1]);
        add_edge(vertices[1], vertices[2]);
        add_edge(vertices[2], vertices[3]);
        add_edge(vertices[3], vertices[4]);
        add_edge(vertices[4], vertices[5]);
        add_edge(vertices[5], vertices[0]);

        auto cell = add_cell(this->edges());
    }


    // -- Agent creation -----------------------------------------------------

    /// Create a Vertex and associate it with the VertexManager
    const auto& add_vertex (const SpaceVec& pos)
    {
        return _vertex_manager.add_agent(pos);
    }

    /// Create a Edge and associate it with the EdgeManager
    /** Creates an edge between vertices a and b. The newly created edge is
     *  added to the adjoint edges of the involved vertices.
     * 
     *  \param a    The start vertex of the new edge
     *  \param b    The end vertex of the new edge
     */
    const auto& add_edge (const std::shared_ptr<Vertex>& a,
                          const std::shared_ptr<Vertex>& b)
    {
        const auto& edge = _edge_manager.add_agent({0., 0.});
        edge->custom_links().a = a;
        edge->custom_links().b = b;

        // Create the required weak links
        for (const auto& v : {a, b}) {
            _vertices_adjoint_edges[v->id()].push_back(edge);
        }

        return edge;
    }

    const auto& add_edge (
        const std::shared_ptr<Vertex>& a,
        const std::shared_ptr<Vertex>& b,
        const std::shared_ptr<Edge>& parent
    )
    {
        const auto& edge = add_edge(a, b);
        edge->state.inherit_from_state(parent->state);
        return edge;
    }


    /// Create a Cell and associate it with the CellManager
    /** Creates a cell using edges as the cell's boundary.
     *  The adjoint objects of the involved edges and vertices are updated.
     * 
     *  \param pos  Where to place the cell
     *  \param edges    The (unordered) container of edges defining the boundary
     *                  of the cell
     *  \param custom_cfg   (optional) A custom cfg used to initialize the
     *                      cell's state. If not provided the agent manager's
     *                      default is used.
     */
    const auto& add_cell (
        AgentContainer<Edge> edges,
        const Config& custom_cfg = {}
    )
    {
        const auto& cell = _cell_manager.add_agent(
            SpaceVec({0., 0.}),
            custom_cfg
        );
        cell->custom_links().edges = this->order_edges(edges);

        // check that edges are anti-clockwise
        if (area_of<true>(cell) < 0.) {
            // they are clockwise -> reverse ordering
            cell->custom_links().edges = reverse_edge_ordering(
                    cell->custom_links().edges);
        }
        
        AgentContainer<Vertex> vertices;
        for (auto [e, flip] : cell->custom_links().edges) {
            if (not flip) {
                vertices.push_back(e->custom_links().a);
            }
            else {
                vertices.push_back(e->custom_links().b);
            }
        }
        cell->custom_links().vertices = vertices;

        // Create the required weak links
        for (const auto& v : vertices) {
            _vertices_adjoint_cells[v->id()].push_back(cell);
        }
        for (const auto& [e, flip] : cell->custom_links().edges) {
            auto [adj_cell_a, adj_cell_b] = _edges_adjoint_cells[e->id()];

            if (not adj_cell_a) {
                adj_cell_a = cell;
            }
            else if (not adj_cell_b) {
                adj_cell_b = cell;
            }

            _edges_adjoint_cells[e->id()] = std::make_pair(adj_cell_a,
                                                           adj_cell_b);
        }

        return cell;
    }

    /// Create a cell and associate it with the CellManager
    /** See add_cell()
     */
    const auto& add_cell(
        AgentContainer<Edge> edges,
        const std::shared_ptr<Cell>& parent_cell
    )
    {
        const auto& cell = add_cell(edges);
        cell->state.inherit_from_state(parent_cell->state);
        return cell;
    }

    /// Remove a vertex
    /** \warning This does not remove the vertex from other entities adjoint
     *           objects
     */
    void remove_vertex (const std::shared_ptr<Vertex>& vertex) {
        vertex->state.remove = true;
        _vertices_adjoint_edges.erase(vertex->id());
        _vertices_adjoint_cells.erase(vertex->id());
        _vertex_manager.remove_agent(vertex);
    }

    /// Remove an edge
    /** \warning This does not remove the edge from other entities adjoint
     *           objects
     */
    void remove_edge (const std::shared_ptr<Edge>& edge) {
        edge->state.remove = true;
        _edges_adjoint_cells.erase(edge->id());
        _edge_manager.remove_agent(edge);
    }

    // see transitions.hh
    bool remove_boundary_edge(const std::shared_ptr<Edge> edge,
        std::function<double(const AgentContainer<Vertex>&,
                             const AgentContainer<Edge>&,
                             const AgentContainer<Cell>&)> get_energy,
        double T1_barrier, double random_number);

    /// Remove an edge
    /** \warning This does not remove the cell from other entities adjoint
     *           objects
     */
    void remove_cell (const std::shared_ptr<Cell>& cell) {
        cell->state.remove = true;
        if (cell->custom_links().nd_cell) {
            cell->custom_links().nd_cell->state.cell_type = 
                Utopia::Models::NotchDelta::CellState::StateType::inactive;
        }
        if (cell->custom_links().c_cell) {
            cell->custom_links().c_cell->state.cell_type = 
                Utopia::Models::Collier::CellState::StateType::inactive;
        }
        _cell_manager.remove_agent(cell);
    }

    /// Order a container of Edges
    /** The established order is either anti-clockwise or clockwise.
     * 
     *  \warning    It fails if the edges do not form a closed boundary
     */
    OrderedEdgeContainer order_edges (AgentContainer<Edge> edges) const
    {
        OrderedEdgeContainer ordered_edges;
       
        auto e_it = edges.begin();
        auto first_vertex = (*e_it)->custom_links().a;
        auto vertex_it = (*e_it)->custom_links().b;        
        ordered_edges.push_back(std::make_pair(*e_it, false));
        e_it = edges.erase(e_it);
        
        while (not edges.empty()) {
            bool found_next = false;
            for (e_it = edges.begin(); e_it != edges.end(); e_it++) {
                auto e = *e_it;
                if (e->custom_links().a == vertex_it) {
                    ordered_edges.push_back(
                        std::make_pair(e, false));
                    vertex_it = e->custom_links().b; // iterate
                    e_it = edges.erase(e_it);
                    found_next = true;
                    break;
                }
                if (e->custom_links().b == vertex_it) {
                    ordered_edges.push_back(
                        std::make_pair(e, true));
                    vertex_it = e->custom_links().a; // iterate
                    e_it = edges.erase(e_it);
                    found_next = true;
                    break;
                }
            }
            if (not found_next) {
                std::cout << "\nFailed to order these edges: \n";
                for (const auto& [e, flip] : ordered_edges) {
                    auto a = e->custom_links().a;
                    auto b = e->custom_links().b;
                    if (not flip) {
                        std::cout << position_of(a) << " to "
                                  << position_of(b) << "\n";
                    }
                    else {
                        std::cout << position_of(b) << " to "
                                  << position_of(a) << "\n";
                    }
                }
                std::cout << " stopped here.\n";
                for (const auto& e : edges) {
                    auto a = e->custom_links().a;
                    auto b = e->custom_links().b;
                    std::cout << position_of(a) << " to " 
                              << position_of(b) << "\n";
                }
                std::cout << std::flush;

                throw std::runtime_error("Could not order edges. Edges to "
                    "order are listed above.");
            }
        }

        return ordered_edges;        
    }

    /// Reverse the ordering of an ordered edge container
    OrderedEdgeContainer reverse_edge_ordering (
            OrderedEdgeContainer ordered_edges)
            const
    {
        OrderedEdgeContainer new_edges;
        // they are clockwise -> flip all edges and inverse order
        for (const auto& [e, flip] : ordered_edges) {
            new_edges.insert(new_edges.begin(), std::make_pair(e, not flip));
        }

        return new_edges;
    }

public:
    /// Return a pointer to the logger of this model
    std::shared_ptr<spdlog::logger> get_logger() const {
        return _log;
    }
}; // EntitiesManager

} // namespace Utopia::Models::PCPVertex

#endif