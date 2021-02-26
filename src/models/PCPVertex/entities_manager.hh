#ifndef UTOPIA_MODELS_PCPVERTEX_ENTITIESMANAGER_HH
#define UTOPIA_MODELS_PCPVERTEX_ENTITIESMANAGER_HH

#include "../Collier/Collier.hh"
#include "../NotchDelta/NotchDelta.hh"

namespace Utopia::Models::PCPVertex {

template<class Model>
struct EntitiesManager {

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
    using EdgeTraits = Utopia::AgentTraits<EdgeState, Update::manual, false,
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

        std::shared_ptr<RotationCellState> rotation_state;
    };
    using CellTraits = Utopia::AgentTraits<CellState, Update::manual, false,
                                           EmptyTag, CellLinks>;

    /// The type of the managed cells
    using Cell = Utopia::Agent<CellTraits, Space>;

    /// The cell-type (hair, support) of the cells
    using CellType = typename Cell::State::CellType;

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


private:
    /// The logger (same as the model this manager resides in)
    const std::shared_ptr<spdlog::logger> _log;
    
    /// Agent manager configuration node
    const Config _cfg;

    /// The physical space the agents are to reside in
    const std::shared_ptr<Space> _space;

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

    /// Function that will be used to prepare positions for adding an agent
    std::function<SpaceVec(const SpaceVec&)> _prepare_pos;

public:
    EntitiesManager (Model& model, const Config& custom_cfg = {})
    :
        _log(spdlog::stdout_color_mt(  model.get_logger()->name()
                                     + ".entities_manager")),
        _cfg(setup_cfg(model, custom_cfg)),
        _space(model.get_space()),
        _vertex_manager(model, setup_vertex_cfg()),
        _edge_manager(model, setup_edge_cfg()),
        _cell_manager(model, setup_cell_cfg()),
        _prepare_pos(setup_prepare_pos_func())
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
    const auto& get_space() const {
        return _space;
    }    
    
    /// Return const reference to the managed vertices
    const auto& vertices () const {
        return _vertex_manager.agents();
    }

    /// Return const reference to the managed edges
    const auto& edges () const {
        return _edge_manager.agents();
    }

    /// Return const reference to the managed cells
    const auto& cells () const {
        return _cell_manager.agents();
    }

    SpaceVec position_of (const Vertex& vertex) const {
        return _space->map_to_absolute_space(vertex.position());
    }

    SpaceVec position_of (const std::shared_ptr<Vertex>& vertex) const {
        return position_of(*vertex);
    }

    /// Move an vertex to a new position in the space
    void move_to(Vertex& vertex,
                 const SpaceVec& pos) const
    {
        return _vertex_manager.move_to(vertex, _prepare_pos(pos));
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

    /// Calculate to where a vertex would move
    /** Vertices move along the self-managed value Vertex::State::f.
     * 
     *  The result is also stored in vertex->state.virtual_position.
     * 
     *  \param beta     The step size along the direction of update
     */
    const SpaceVec& displace_virtual (const std::shared_ptr<Vertex>& vertex,
                                      double beta) const {
        if (std::get<double>(vertex->state.virtual_pos) == beta) {
            return std::get<SpaceVec>(vertex->state.virtual_pos);
        }

        SpaceVec pos = position_of(vertex) + beta * vertex->state.f;
        
        if (this->_space->periodic) {
            pos = this->_space->map_into_space(pos);
        }
        else {
            if (not this->_space->contains(pos)) {
                throw OutOfSpace(pos, this->_space, "Could not move agent!");
            }
        }

        vertex->state.virtual_pos = std::make_pair(beta, pos);
        return std::get<SpaceVec>(vertex->state.virtual_pos);
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

    auto length_of (const std::shared_ptr<Edge>& edge, double beta) const
    {
        if (beta == 0) {
            return distance(edge->custom_links().a, edge->custom_links().b);
        }
        else {
            const SpaceVec &a = displace_virtual(edge->custom_links().a, beta);
            const SpaceVec &b = displace_virtual(edge->custom_links().b, beta);
            return _space->distance(a, b);
        }
    }

    /// Calculate the perimeter of a cell
    double perimeter_of (const std::shared_ptr<Cell>& cell,
                         double beta = 0.) const
    {
        return perimeter_of(cell->custom_links().edges, beta);
    }

    /// Calculate the perimeter of a boundary.
    /** This can be a boundary of a cell or any other closed loop of edges.
     */
    double perimeter_of (const OrderedEdgeContainer& boundary,
                         double beta) const
    {
        double perimeter = 0.;
        for (const auto [e, flip] : boundary) {
            perimeter += length_of(e, beta);
        }
        return perimeter;
    }

    /// Calculate the area of a cell
    /** \param beta     The update step length at which to predict the area
     *  \tparam get_sign    If false the area can have a sign. It will be
     *                      negative when the cell's edges are ordered clockwise
     *                      instead of anti-clockwise.
     *                      If true, throws on counter-clockwise ordering
     * 
     *  \note Does not check whether the edges are ordered and in which
     *        orientation they are ordered.
     */
    template <bool get_sign = false>
    double area_of (const std::shared_ptr<Cell>& cell, double beta = 0.) const
    {
        return area_of<get_sign>(cell->custom_links().edges, beta);
    }

    /// Calculate the area of a boundary
    /** \param boundary The edges defining a boundary. This can be of a cell
     *                  or any other closed loop.
     *  \param beta     The update step length at which to predict the area
     *  \tparam get_sign    If false the area can have a sign. It will be
     *                      negative when the cell's edges are ordered clockwise
     *                      instead of anti-clockwise.
     *                      If true, throws on counter-clockwise ordering
     * 
     *  \note Does not check whether the edges are ordered and in which
     *        orientation they are ordered.
     */
    template <bool get_sign = false>
    double area_of(const OrderedEdgeContainer& boundary,
                   double beta = 0.) const 
    {
        static_assert(Space::dim == 2, "Area of a cell is only implemented for "
                      "2 dimensional space!");

        // define a reference in space
        auto [e, flip] = boundary.front();
        std::shared_ptr<Vertex> reference;
        if (not flip) { reference = e->custom_links().a; }
        else { reference = e->custom_links().b; }

        SpaceVec ref;
        if (beta == 0.) {
            ref = position_of(reference);
        }
        else {
            ref = displace_virtual(reference, beta);
        }

        double area = 0.;
        if (beta == 0.) {
            for (const auto [e, flip] : boundary) {
                SpaceVec a = position_of(e->custom_links().a);
                SpaceVec b = position_of(e->custom_links().b);

                // define the vertices positions relative to the reference
                /* this is important in periodic space to calculate with 
                    * "real" coordinates */
                a = ref + _space->displacement(ref, a);
                b = ref + _space->displacement(ref, b);

                if (flip) { std::swap(a, b); }

                area += a[0] * b[1] - b[0] * a[1];
            }
        }
        else {
            for (const auto [e, flip] : boundary) {
                SpaceVec a = this->displace_virtual(e->custom_links().a,
                                                    beta);
                SpaceVec b = this->displace_virtual(e->custom_links().b,
                                                    beta);

                a = ref + _space->displacement(ref, a);
                b = ref + _space->displacement(ref, b);
                
                if (flip) { std::swap(a, b); }

                area += a[0] * b[1] - b[0] * a[1];
            }
        }

        area /= 2.;

        if constexpr (not get_sign) {
            // With beta = 0, negative area not allowed
            if (area < 0. and beta == 0.) {
                throw std::runtime_error(fmt::format(
                    "Negative area ({}) of a boundary defining a cell!", area));
            }
            else if (area < 0.) {
                return std::nan("1");
            }
        }

        return area;
    }

    /// Calculate the shape index of a cell
    /** shape index \f$ p = P / \sqrt(A) \f$ with the cell's perimeter \f$ P \f$
     *  and area \f$ A \f$.
     * 
     *  \param beta     The update step length at which to predict the area
     */ 
    double shape_index_of (const std::shared_ptr<Cell>& cell,
                           double beta = 0.) const
    {
        return shape_index_of(cell->custom_links().edges, beta);
    }

    /// Calculate the shape index of a boundary
    /** shape index \f$ p = P / \sqrt(A) \f$ with the cell's perimeter \f$ P \f$
     *  and area \f$ A \f$.
     * 
     *  \param beta     The update step length at which to predict the area
     */ 
    double shape_index_of (const OrderedEdgeContainer& boundary,
                           double beta = 0.) const
    {
        return perimeter_of(boundary, beta) / sqrt( area_of(boundary, beta) );
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

        // define a reference in space
        auto [e, flip] = boundary.front();
        std::shared_ptr<Vertex> reference;
        if (not flip) { reference = e->custom_links().a; }
        else { reference = e->custom_links().b; }

        double area = 0.;
        SpaceVec center(arma::fill::zeros);
        for (const auto [e, flip] : boundary) {
            // define the vertices positions relative to the reference
            /* this is important in periodic space to calculate with "real"
             * coordinates */
            SpaceVec a = position_of(reference) +
                         displacement(reference, e->custom_links().a);
            SpaceVec b = position_of(reference) + 
                         displacement(reference, e->custom_links().b);

            if (flip) { std::swap(a, b); }

            double da = a[0] * b[1] - b[0] * a[1];
            area += da;
            center += (a + b) * da;
        }
        area /= 2;
        if (area == 0) { area += 1e-12; }
        return _space->map_into_space(center / (6 * area));
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
     *  in their boundary
     * 
     *  \param edge   The considered edge
     */
    std::pair<std::shared_ptr<Cell>, std::shared_ptr<Cell>> adjoints_of(
            const std::shared_ptr<Edge>& edge) const
    {
        return _edges_adjoint_cells.at(edge->id());
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

    /// The neighboring cells to a cell of type `hair`
    /** The neighbors and next neighbors of cell that are of CellType::hair.
     * 
     *  \param cell   The considered cell
     */
    AgentContainer<Cell> hair_neighbors_of(
            const std::shared_ptr<Cell>& cell) const
    {
        std::set<std::shared_ptr<Cell>> hair_neighbors;

        auto neighbors = this->neighbors_of(cell);

        for (const auto & nb : neighbors) {
            if (nb->state.type == CellType::hair)
            {
                hair_neighbors.insert(nb);
            }
            else {
                auto next_neighbors = this->neighbors_of(nb);
                for (const auto& nn : next_neighbors) {
                    if (nn->state.type == CellType::hair)
                    {
                        hair_neighbors.insert(nn);
                    }
                }
            }
        }

        hair_neighbors.erase(cell);
        hair_neighbors.erase(nullptr);

        return AgentContainer<Cell>(hair_neighbors.begin(),
                                    hair_neighbors.end());
    }

    // see transitions.hh
    template <typename EdgeParamMatrix>
    void divide_cell(const std::shared_ptr<Cell> cell, double division_angle,
            EdgeParamMatrix linetension, EdgeParamMatrix edge_contractility);

    // see transitions.hh
    template <typename EdgeParamMatrix>
    bool remove_edge_T1 (const std::shared_ptr<Edge> edge,
        EdgeParamMatrix linetension, EdgeParamMatrix contractility,
        std::function<double(const AgentContainer<Edge>&,
                             const AgentContainer<Cell>&)> get_energy,
        double separation, double T1_barrier, double random_number);

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

        const auto& cells = this->cells();
        const auto& edges = this->edges();

        OrderedEdgeContainer boundary{};
        boundary.reserve(cells.size());

        auto edge = *std::find_if(
            edges.begin(), edges.end(),
            [this](const auto& edge) {
                return this->is_1_cell_boundary_edge(edge);
            });

        const auto start = edge->custom_links().a;
        auto iter = edge->custom_links().b;
        boundary.push_back(std::make_pair(edge, false));

        while (iter != start) {
            const auto tmp_edges = adjoint_edges_of(iter);
            edge = *std::find_if(
                tmp_edges.begin(), tmp_edges.end(),
                [this, edge](const auto& e_it) {
                    if (e_it == edge) {
                        return false;
                    }
                    return this->is_1_cell_boundary_edge(e_it);
                });

            if (iter == edge->custom_links().a) {
                boundary.push_back(std::make_pair(edge, false));
                iter = edge->custom_links().b;
            }
            else {
                boundary.push_back(std::make_pair(edge, true));
                iter = edge->custom_links().a;
            }
        }

        double area = std::accumulate(
            boundary.begin(), boundary.end(), 0.,
            [this](double val, const auto e_pair) {
                const auto [e, flip] = e_pair;
                
                SpaceVec a = this->position_of(e->custom_links().a);
                SpaceVec b = this->position_of(e->custom_links().b);

                if (flip) { std::swap(a, b); }

                return val + (a[0] * b[1] - b[0] * a[1]);
            });
        area /= 2.;

        if (area < 0.) {
            boundary = this->reverse_edge_ordering(boundary);
            area = fabs(area);
        }

        boundary.shrink_to_fit();

        return boundary;
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
        // method not found!
        else {
            throw KeyError(method, _cfg, "Method not implemented to setup agent "
                           "manager! Please choose one of the following setup "
                           "methods: 'hexagonal', 'single'.");
        }
    }

    // see initialisation.hh
    void setup_agents_hexagonal_structure(const Config& cfg);

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

        auto cell = add_cell(pos, this->edges());
    }


    // -- Agent creation -----------------------------------------------------

    /// Create a Vertex and associate it with the VertexManager
    auto add_vertex (const SpaceVec& pos, const Config& custom_cfg = {})
    {
        return _vertex_manager.add_agent(_prepare_pos(pos), custom_cfg);
    }

    /// Create a Edge and associate it with the EdgeManager
    /** Creates an edge between vertices a and b. The newly created edge is
     *  added to the adjoint edges of the involved vertices.
     * 
     *  \param a    The start vertex of the new edge
     *  \param b    The end vertex of the new edge
     *  \param custom_cfg   (optional) A custom cfg used to initialize the
     *                      edge's state. If not provided the agent manager's
     *                      default is used.
     */
    auto add_edge (std::shared_ptr<Vertex> a, std::shared_ptr<Vertex> b,
                   const Config& custom_cfg = {})
    {
        auto e = _edge_manager.add_agent({0., 0.}, custom_cfg);
        e->custom_links().a = a;
        e->custom_links().b = b;

        // Create the required weak links
        for (const auto& v : {a, b}) {
            _vertices_adjoint_edges[v->id()].push_back(e);
        }

        return e;
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
    auto add_cell (const SpaceVec& pos,
                   AgentContainer<Edge> edges,
                   const Config& custom_cfg = {})
    {
        auto cell = _cell_manager.add_agent(_prepare_pos(pos), custom_cfg);
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
    auto add_cell(AgentContainer<Edge> edges,
                  const std::shared_ptr<Cell>& parent_cell)
    {
        return add_cell(SpaceVec({0., 0.}), edges,
                        parent_cell->state.create_cfg_from_props());
    }

    /// Create a cell with properties inherited from a parent cell
    /** Arguments passed on to add_cell() without precising the pos. The cell's
     *  barycenter will be calculated from the positions of the vertices
     *  defining the boundary.
     */
    auto add_cell(const SpaceVec& pos,
                  AgentContainer<Edge> edges,
                  const std::shared_ptr<Cell>& parent_cell)
    {
        return add_cell(pos, edges, parent_cell->state.create_cfg_from_props());
    }

    /// Remove a vertex
    /** \warning This does not remove the vertex from other entities adjoint
     *           objects
     */
    void remove_vertex (const std::shared_ptr<Vertex>& vertex) {
        _vertices_adjoint_edges.erase(vertex->id());
        _vertices_adjoint_cells.erase(vertex->id());
        _vertex_manager.remove_agent(vertex);
    }

    /// Remove an edge
    /** \warning This does not remove the edge from other entities adjoint
     *           objects
     */
    void remove_edge (const std::shared_ptr<Edge>& edge) {
        _edges_adjoint_cells.erase(edge->id());
        _edge_manager.remove_agent(edge);
    }

    // see transitions.hh
    bool remove_boundary_edge(const std::shared_ptr<Edge> edge,
        std::function<double(const AgentContainer<Edge>&,
                             const AgentContainer<Cell>&)> get_energy,
        double T1_barrier, double random_number);

    /// Remove an edge
    /** \warning This does not remove the cell from other entities adjoint
     *           objects
     */
    void remove_cell (const std::shared_ptr<Cell>& cell) {
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

    /// Depending on periodicity, return the function to prepare positions of
    /// agents before they are added
    /** The function that is used to prepare a position before an agent is
     *  added (if passed explicitly) depends on whether the space is periodic or
     *  not. In the case of a periodic space the position is automatically
     *  mapped into space again across the borders. For a nonperiodic space a
     *  position outside of the borders will throw an error.
     */
    std::function<SpaceVec(const SpaceVec&)> setup_prepare_pos_func() const {
        for (size_t i = 0; i < this->_space->dim; i++) {
            if (this->_space->get_domain_size()[i] < 1) {
                std::cout << "Received domain size: "
                          << this->_space->get_domain_size()
                          << std::endl << std::flush;
                throw std::invalid_argument("Domain size has to be larger than "
                    "1 in all entries!");
                // NOTE this can be remove if consistently integrated in
                //      agent_manager.
            }
        }

        // If periodic, map the position back into space
        if (_space->periodic) {
            return
                [this](const SpaceVec& pos){
                    return this->_space->map_to_relative_space(
                                this->_space->map_into_space(pos));
                };
        }
        // If non-periodic, check wether the position is valid
        else {
            return
                [this](const SpaceVec& pos) {
                    if (not _space->contains(pos)) {
                        throw OutOfSpace(pos, _space,
                                         "Given position is out of space!");
                    }
                    return this->_space->map_to_relative_space(pos);
                };
        }
    }
}; // EntitiesManager

} // namespace Utopia::Models::PCPVertex

#endif