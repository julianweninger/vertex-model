#ifndef UTOPIA_MODELS_PCPVERTEX_ENTITIESMANAGER_HH
#define UTOPIA_MODELS_PCPVERTEX_ENTITIESMANAGER_HH

namespace Utopia::Models::PCPVertex {

template<class Model>
struct CustomAgentManager {
public:
    /// The type of the space
    using Space = typename Model::Space;

    /// The type of a vector in space
    using SpaceVec = typename Space::SpaceVec;

    /// The type of the configuration
    using Config = typename Model::Config;

    /// The traits of a Vertex
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
        std::vector<std::pair<std::shared_ptr<Edge>, bool>> edges;
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

public:
    CustomAgentManager (Model& model, const Config& custom_cfg = {})
    :
        _log(model.get_logger()),
        _cfg(setup_cfg(model, custom_cfg)),
        _space(model.get_space()),
        _vertex_manager(model, setup_vertex_cfg(model)),
        _edge_manager(model, setup_edge_cfg(model)),
        _cell_manager(model, setup_cell_cfg(model))
    {
        setup_agents();
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

    /// The displacement between two vertices
    /** \details see Utopia::Space::displacement
     */
    auto displacement (const Vertex& a, const Vertex& b) const {
        return _space->displacement(a.position(), b.position());
    }

    /// The displacement between two vertices
    /** \details see Utopia::Space::displacement
     */
    auto displacement (const std::shared_ptr<Vertex>& a,
                       const std::shared_ptr<Vertex>& b) const {
        return _space->displacement(a->position(), b->position());
    }

    /// The distance between two vertices
    /** \details see Utopia::Space::distance
     */
    auto distance (const Vertex& a, const Vertex& b) const {
        return _space->distance(a.position(), b.position());
    }

    /// The distance between two vertices
    /** \details see Utopia::Space::distance
     */
    auto distance (const std::shared_ptr<Vertex>& a,
                   const std::shared_ptr<Vertex>& b) const {
        return _space->distance(a->position(), b->position());
    }

    /// Calculate the perimeter of a cell
    double perimeter_of (const Cell& cell) const {
        double perimeter = 0.;
        for (auto [e, flip] : cell.custom_links().edges) {
            perimeter += distance(e->custom_links().a, e->custom_links().b);
        }
        return perimeter;
    }

    /// Calculate the perimeter of a cell
    double perimeter_of (const std::shared_ptr<Cell>& cell) const {
        return perimeter_of(*cell);
    }

    /// Calculate the area of a cell
    /** \note   It is assumed that the edges are ordered anti-clockwise.
     *          If the edges are ordered clockwise, the area is correct but of 
     *          negative sign.
     */
    double area_of (const Cell& cell) const {
        static_assert(Space::dim == 2, "Area of a cell is only implemented for "
                      "2 dimensional space!");

        // the ordered edges using flip boolian: [edge, flip]
        const auto& edges = cell.custom_links().edges;

        // define a reference in space
        auto [e, flip] = edges.front();
        std::shared_ptr<Vertex> reference;
        if (not flip) { reference = e->custom_links().a; }
        else { reference = e->custom_links().b; }

        double area = 0.;
        for (const auto [e, flip] : edges) {
            // define the vertices positions relative to the reference
            /* this is important in periodic space to calculate with "real"
             * coordinates */
            SpaceVec a = reference->position() +
                         displacement(reference, e->custom_links().a);
            SpaceVec b = reference->position() + 
                         displacement(reference, e->custom_links().b);

            if (flip) { std::swap(a, b); }

            area += a(0) * b(1) - b(0) * a(1);
        }

        // since the edges have to be ordered anti-clockwise 
        return 0.5 * area;
    }

    /// Calculate the area of a cell
    /** \note   It is assumed that the edges are ordered anti-clockwise.
     *          If the edges are ordered clockwise, the area is correct but of 
     *          negative sign.
     */
    double area_of (const std::shared_ptr<Cell>& cell) const {
        return area_of(*cell);
    }

    /// Returns the barycenter of the given cell
    /** \note   It is assumed that the edges are ordered anti-clockwise or 
     *          clock-wise
     */
    SpaceVec barycenter_of (const Cell& cell) const {
        static_assert(Space::dim == 2, "Center of a cell is only implemented "
                      "for 2 dimensional space!");

        // the ordered edges using flip boolian: [edge, flip]
        const auto& edges = cell.custom_links().edges;

        // define a reference in space
        auto [e, flip] = edges.front();
        std::shared_ptr<Vertex> reference;
        if (not flip) { reference = e->custom_links().a; }
        else { reference = e->custom_links().b; }

        double area = 0.;
        SpaceVec center(arma::fill::zeros);
        for (const auto [e, flip] : edges) {
            // define the vertices positions relative to the reference
            /* this is important in periodic space to calculate with "real"
             * coordinates */
            SpaceVec a = reference->position() +
                         displacement(reference, e->custom_links().a);
            SpaceVec b = reference->position() + 
                         displacement(reference, e->custom_links().b);

            if (flip) { std::swap(a, b); }

            double da = a(0) * b(1) - b(0) * a(1);
            area += da;
            center += (a + b) * da;
        }
        area /= 2;
        return center / (6 * area);
    }
    
    /// Returns the barycenter of the given cell
    /** \note   It is assumed that the edges are ordered anti-clockwise or 
     *          clock-wise
     */
    SpaceVec barycenter_of (const std::shared_ptr<Cell>& cell) const {
        return barycenter_of(*cell);
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
    Config setup_vertex_cfg (Model& model) {
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
    Config setup_edge_cfg (Model& model) {
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
    Config setup_cell_cfg (Model& model) {
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

    void setup_agents() {
        auto method = get_as<std::string>("setup_method", _cfg);
        if (not _cfg["setup_params"]) {
            throw KeyError ("setup_params", _cfg, "No parameters provided "
                            "for the setup of agent manager!");
        }
        if (not _cfg["setup_params"][method]) {
            throw KeyError (method, _cfg["setup_params"], "No parameters "
                            "provided for the setup of agent manager with the "
                            "specified method!");
        }
        if (method == "hexagonal") {
            this->setup_agents_hexagonal_structure(
                    _cfg["setup_params"]["hexagonal"]);
        }
        else if (method == "single") {
            this->setup_agents_single_hexagon(_cfg["setup_params"]["single"]);
        }
        else {
            throw KeyError(method, _cfg, "Method not implemented to setup agent "
                           "manager! Please choose one of the following setup "
                           "methods: 'hexagonal'");
        }
    }

    void setup_agents_hexagonal_structure(const Config& cfg);

    /// Setup a single cell of hexagonal shape in center of space
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

    /// Create a Vertex and associate it with the VertexManager
    auto add_vertex (const SpaceVec& pos, const Config& custom_cfg = {})
    {
        return _vertex_manager.add_agent(pos, custom_cfg);
    }

    /// Create a Edge and associate it with the EdgeManager
    auto add_edge (std::shared_ptr<Vertex> a, std::shared_ptr<Vertex> b,
                   const Config& custom_cfg = {})
    {
        auto e = _edge_manager.add_agent({arma::datum::nan, arma::datum::nan},
                                         custom_cfg);
        e->custom_links().a = a;
        e->custom_links().b = b;
        return e;
    }

    /// Create a Cell and associate it with the CellManager
    auto add_cell (const SpaceVec& pos,
                   AgentContainer<Edge> edges,
                   const Config& custom_cfg = {})
    {
        auto c = _cell_manager.add_agent(pos);
        c->custom_links().edges = this->order_edges(edges);

        // check that edges are anti-clockwise
        if (area_of(*c) < 0.) {
            // they are clockwise -> flip all edges
            auto edges = c->custom_links().edges;
            for (int i = 0; i < edges.size(); i++) {
                edges[i] = std::make_pair(std::get<0>(edges[i]),
                                          not std::get<1>(edges[i]));
            }
            c->custom_links().edges = edges;
        }
        
        AgentContainer<Vertex> vertices;
        for (auto [e, flip] : c->custom_links().edges) {
            if (not flip) {
                vertices.push_back(e->custom_links().a);
            }
            else {
                vertices.push_back(e->custom_links().b);
            }
        }
        c->custom_links().vertices = vertices;

        return c;
    }

    /// Order a container of Edges
    /** \details    The established order is either anti-clockwise or clockwise.
     *  \warning    It fails if the edges do not form a closed boundary
     */
    std::vector<std::pair<std::shared_ptr<Edge>, bool>> order_edges (
            AgentContainer<Edge> edges)
    {
        std::vector<std::pair<std::shared_ptr<Edge>, bool>> ordered_edges;
       
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
                        std::cout << a->position() << " to "
                                  << b->position() << "\n";
                    }
                    else {
                        std::cout << b->position() << " to "
                                  << a->position() << "\n";
                    }
                }
                std::cout << " stopped here.\n";
                for (const auto& e : edges) {
                    auto a = e->custom_links().a;
                    auto b = e->custom_links().b;
                    std::cout << a->position() << " to " << b->position() << "\n";
                }
                std::cout << std::flush;

                throw std::runtime_error("Could not order edges. Edges to "
                    "order are listed above.");
            }
        }

        return ordered_edges;        
    }
}; // CustomAgentManager

} // namespace Utopia::Models::PCPVertex

#endif