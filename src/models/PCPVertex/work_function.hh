#ifndef UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION_HH
#define UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION_HH


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction {

/// @brief  The class of a tension term in the work function
/// @tparam Model   The model within which the AgentManager lives
/** A tension is the derivative of the energy to the length of an junction */
template <typename Model>
class WorkFunctionTerm {
public:
    /// The type of {x, y} coordinates
    using SpaceVec = typename Model::SpaceVec;

    /// The type of Agent manager
    /** Here the Vertices, Edges, and Cells are stored, 
     *  as well as how to calculate their properties, e.g.
     *  - position of vertices
     *  - distance between vertices
     *  - length of edges
     *  - area or perimeter cells
     *  - adjoint cells to edges
     *  - neighbors to a cell
     */ 
    using AgentManager = typename Model::AgentManager;

    /// The type of a vertex
    using Vertex = typename Model::Vertex;

    /// The type of an edge
    using Edge = typename Model::Edge;

    /// The type of a cell
    using Cell = typename Model::Cell;

protected:
    /// The name of this work-function term
    const std::string _name;

    /// The manager of vertices, edges, cells, and their relations
    /** See entities_manager.hh */
    const AgentManager& _am;

public:
    // --- Constructors and destructors ---------------------------------------
    WorkFunctionTerm (std::string name,
                      [[maybe_unused]] const DataIO::Config& cfg,
                      const Model& model)
    :
        _name(name),
        _am(model.get_am())
    { }

    virtual ~WorkFunctionTerm() { }

    // --- Force calculations  -------------------------------------------------
    /// @brief The instruction how to calculate forces acting on vertices
    ///        from this term, i.e the derivative of the energy.
    virtual void compute_and_set_forces () = 0;


    /// @brief  The force acting on a vertex
    /** i.e. the derivate of the work-function w.r.t the position of vertex i.
     * 
     *  Note, needs to be included in compute_and_set_forces () to be effective.
     *  
     *  @param vertex     The vertex i
     *  @return force (double, double)
    **/ 
    virtual SpaceVec compute_force
    ([[maybe_unused]] const std::shared_ptr<Vertex>& vertex) const
    {
        return SpaceVec({0., 0.});
    }

    /// @brief  The tension of an edge
    /** i.e. the derivate of the work-function w.r.t the length of edge <i,j>.
     * 
     *  Note, needs to be included in compute_and_set_forces () to be effective.
     *  
     *  @param edge     The edge <i,j>
     *  @return tension (double)
    **/ 
    virtual double compute_tension
    ([[maybe_unused]] const std::shared_ptr<Edge>& edge) const
    {
        return 0.;
    }

    /// @brief  The pressure of a cell
    /** i.e. the derivate of the work-function w.r.t the area of cell \alpha.
     * 
     *  Note, needs to be included in compute_and_set_forces () to be effective.
     *  
     *  @param cell     The cell \alpha
     *  @return pressure (double)
    **/ 
    virtual double compute_pressure
    ([[maybe_unused]] const std::shared_ptr<Cell>& cell) const
    {
        return 0.;
    }
    

    // --- Energy calculations  -----------------------------------------------
    /// @brief  The energy associated with vertex i.
    /** Note, needs to be included in compute_energy(vertices, ...) to be 
     *  effective.
     * 
     *  @param vertex     The vertex i
     *  @return energy (double)
    **/ 
    virtual double compute_energy
    ([[maybe_unused]] const std::shared_ptr<Vertex>& vertex) const 
    {
        return 0.;
    }
    
    /// @brief  The energy associated with edge <i,j>.
    /** Note, needs to be included in compute_energy(..., edges, ...) to be 
     *  effective.
     * 
     *  @param edge     The edge <i,j>
     *  @return energy (double)
    **/ 
    virtual double compute_energy
    ([[maybe_unused]] const std::shared_ptr<Edge>& edge) const 
    {
        return 0.;
    }
    
    /// @brief  The energy associated with cell \alpha.
    /** Note, needs to be included in compute_energy(..., cells) to be 
     *  effective.
     * 
     *  @param cell     The cell \alpha
     *  
     *  @return energy (double)
    **/ 
    virtual double compute_energy
    ([[maybe_unused]] const std::shared_ptr<Cell>& cell) const 
    {
        return 0.;
    }


    /// @brief  The total energy for a subset of entities.
    /** @param vertices  The container of vertices
     *  @param edges     The container of edges
     *  @param cells     The container of cells
     *  
     *  @return energy (double)
    **/ 
    virtual double compute_energy (
        const AgentContainer<Vertex>& vertices,
        const AgentContainer<Edge>& edges,
        const AgentContainer<Cell>& cells
    ) const = 0;

    /// @brief  The total energy for all entities.
    /** 
     *  @return energy (double)
    **/ 
    double compute_energy () const {
        return compute_energy(
            this->_am.vertices(),
            this->_am.edges(),
            this->_am.cells()
        );
    }


    // --- Continuous updates -------------------------------------------------
    /// Continuous update to WF-term
    /** This function is continuously called alongside compute_and_set_forces.
     */
    virtual void update ([[maybe_unused]] double dt) { return; }


    // --- Setters and Getters ------------------------------------------------
    /// Update the parameter values
    virtual void update_parameters (const DataIO::Config& cfg) = 0;

    /// The name of the WF-term
    const std::string& get_name () const {
        return _name;
    }

    // --- Data writers -------------------------------------------------------
    /// The names to provide the data writer task for each vertex
    virtual std::vector<std::string> write_task_vertex_properties_names () const {
        return std::vector<std::string>({});
    }

    /// The data for each vertex to be written by the data writer
    virtual std::vector<std::vector<double>> write_vertex_properties () const {
        return std::vector<std::vector<double>>({});
    }

    /// The names to provide the data writer task for each edge
    virtual std::vector<std::string> write_task_edge_properties_names () const {
        return std::vector<std::string>({});
    }

    /// The data for each edge to be written by the data writer
    virtual std::vector<std::vector<double>> write_edge_properties () const {
        return std::vector<std::vector<double>>({});
    }
    
    /// The names to provide the data writer task for each cell
    virtual std::vector<std::string> write_task_cell_properties_names () const {
        return std::vector<std::string>({});
    }
    
    /// The data for each cell to be written by the data writer
    virtual std::vector<std::vector<double>> write_cell_properties () const {
        return std::vector<std::vector<double>>({});
    }


    /// The names to provide the data writer task for each vertex-energy
    virtual std::vector<std::string> write_task_vertex_energies_names () const {
        return std::vector<std::string>({});
    }

    /// The energy-data for each vertex to be written by the data writer
    virtual std::vector<std::vector<double>> write_vertex_energies () const {
        return std::vector<std::vector<double>>({});
    }

    /// The names to provide the data writer task for each edge-energy
    virtual std::vector<std::string> write_task_edge_energies_names () const {
        return std::vector<std::string>({});
    }

    /// The energy-data for each edge to be written by the data writer
    virtual std::vector<std::vector<double>> write_edge_energies () const {
        return std::vector<std::vector<double>>({});
    }

    /// The names to provide the data writer task for each cell-energy
    virtual std::vector<std::string> write_task_cell_energies_names () const {
        return std::vector<std::string>({});
    }

    /// The energy-data for each cell to be written by the data writer
    virtual std::vector<std::vector<double>> write_cell_energies () const {
        return std::vector<std::vector<double>>({});
    }

    // --- Tests --------------------------------------------------------------
    virtual bool test_constraints 
    ([[maybe_unused]] const std::shared_ptr<spdlog::logger>& logger) const 
    {
        return true;
    }
};


/// @brief The linetension term
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `linetension`: the tension \f$ k \f$
 */
template <typename Model>
class Linetension : public WorkFunctionTerm<Model>
{
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

private:
    double _linetension;

public:
    Linetension (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _linetension(get_as<double>("linetension", cfg))
    { }

    void compute_and_set_forces () final {
        for (const auto& edge : this->_am.edges()) {
            SpaceVec displ = this->_am.displacement(edge);
            SpaceVec director = displ / arma::norm(displ);

            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;

            a->state.add_force(+ _linetension * director);
            b->state.add_force(- _linetension * director);
        }
    }

    SpaceVec compute_force(const std::shared_ptr<Vertex>& vertex) const final {
        SpaceVec force({0., 0.});
        for (const auto& edge : this->_am.adjoint_edges_of(vertex)) {
            SpaceVec displ = this->_am.displacement(edge);
            SpaceVec director = displ / arma::norm(displ);

            if (edge->custom_links().a == vertex) {
                force += _linetension * director;
            }
            else {
                force -= _linetension * director;
            }
        }
        return force;
    }

    double compute_tension([[maybe_unused]] const std::shared_ptr<Edge>& edge)
    const final
    {
        return _linetension;
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return _linetension * this->_am.length_of(edge);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& edge : edges) {
            energy += compute_energy(edge);
        }
        return energy;
    }

    void update_parameters (const DataIO::Config& cfg) final {
        _linetension  = get_as<double>("linetension", cfg, _linetension);
    }



    std::vector<std::string> write_task_edge_energies_names () const final {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const final {
        std::vector<double> energies({});
        std::vector<double> tensions({});
        
        for (const auto& edge : this->_am.edges()) {
            energies.push_back(compute_energy(edge));
            tensions.push_back(compute_tension(edge));
        }

        return std::vector<std::vector<double>>({
            energies,
            tensions
        });
    }
};


arma::mat setup_symmetric_matrix (
    const std::vector<std::vector<double>>& values
)
{
    // If the mask is empty, just return the default cost
    if (values.size() == 0){
        throw std::runtime_error("Matrix cannot be empty!");
    }

    std::set<std::size_t> dim_values;
    for (const auto& v : values){
        dim_values.insert(v.size());
    }
    dim_values.insert(values.size());
    if (dim_values.size() != 1){
        throw std::runtime_error("The provided value-matrix does not have "
                                 "quadratic shape!");
    }

    arma::mat mat(values.size(), values.size(), arma::fill::zeros);

    for (std::size_t i = 0; i < values.size(); i++) {
        for (std::size_t j = 0; j < values.size(); j++) {
            mat(i, j) = values[i][j];
            if (i > j and fabs(mat(i,j) - mat(j,i)) > 1.e-8) {
                std::cout << mat;
                throw std::runtime_error("Matrix is not symmetric!");
            }
        }
    }

    return mat;
}

/// @brief The linetension term with heterotypic values
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `linetension`: the tension \f$ k \f$. A symmetric matrix. 
 *              The i,j coordinates map to the type of cell on left and right
 *              side. 
 *      - `boundary_type`: To which value a boundary cell is mapped.
 */
template <typename Model>
class LinetensionHeterotypic : public WorkFunctionTerm<Model>
{
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

    typedef std::vector< std::vector<double> > stdmat;

private:
    arma::mat _linetension;

    std::size_t _boundary_type;

    const double& get_linetension (const std::shared_ptr<Edge>& edge) const {
        const auto& [ca, cb] = this->_am.adjoints_of(edge);
        std::size_t type_a, type_b;
        if (ca) { type_a = ca->state.type; }
        else { type_a = _boundary_type; }
        if (cb) { type_b = cb->state.type; }
        else { type_b = _boundary_type; }

        if (std::max(type_a, type_b) > _linetension.n_rows)
        {
            std::cout << _linetension << std::endl;

            throw std::runtime_error(fmt::format(
                "In WF-term LinetensionHeterotypic, no parameter registered "
                "for cells of type {}. Parameters for {} types registered.",
                std::max(type_a, type_b),
                _linetension.n_rows
            ));
        }

        return _linetension.at(type_a, type_b);
    }

public:
    LinetensionHeterotypic (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _linetension(setup_symmetric_matrix(
            get_as<stdmat>("linetension", cfg))),
        _boundary_type(get_as<std::size_t>("boundary_type", cfg))
    { }

    void compute_and_set_forces () final {
        for (const auto& edge : this->_am.edges()) {
            SpaceVec displ = this->_am.displacement(edge);
            SpaceVec director = displ / arma::norm(displ);

            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;

            a->state.add_force(+ get_linetension(edge) * director);
            b->state.add_force(- get_linetension(edge) * director);
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return get_linetension(edge);
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return get_linetension(edge) * this->_am.length_of(edge);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& edge : edges) {
            energy += compute_energy(edge);
        }
        return energy;
    }

    void update_parameters (const DataIO::Config& cfg) final {
        if (cfg["linetension"]) {
            _linetension = setup_symmetric_matrix(get_as<stdmat>(
                "linetension", cfg
            ));
        }
        _boundary_type = get_as<std::size_t>(
            "boundary_type", cfg, _boundary_type);
    }



    std::vector<std::string> write_task_edge_properties_names () const final {
        return std::vector<std::string>({
            "linetension"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const final {
        std::vector<double> tensions({});
        for (const auto& edge : this->_am.edges()) {
            tensions.push_back(get_linetension(edge));
        }
        return std::vector<std::vector<double>>({tensions});
    }

    std::vector<std::string> write_task_edge_energies_names () const final {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const final {
        std::vector<double> energies({});
        std::vector<double> tensions({});
        
        for (const auto& edge : this->_am.edges()) {
            energies.push_back(compute_energy(edge));
            tensions.push_back(compute_tension(edge));
        }

        return std::vector<std::vector<double>>({
            energies,
            tensions
        });
    }
};

/// @brief The linetension term with Ornstein-Uhlenbeck fluctuations
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `linetension`: the tension \f$ k \f$. A symmetric matrix. 
 *              The i,j coordinates map to the type of cell on left and right
 *              side. 
 *      - `boundary_type`: To which value a boundary cell is mapped.
 */
template <typename Model>
class LinetensionFluctuations : public WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

    using RNG = typename Model::Base::RNG;

private:
    const std::string _name;

    const std::shared_ptr<RNG> _rng;

    /// A [0,1]-range normal distribution
    std::normal_distribution<double> _normal_distr;

    double _tau;

    double _amplitude;

    double get_linetension (const std::shared_ptr<Edge>& edge) const {
        if (not edge->state.has_parameter(_name)) {
            edge->state.register_parameter(_name, 0.);
        }
        return edge->state.get_parameter(_name);
    }

public:
    LinetensionFluctuations (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _name(get_as<std::string>("name", cfg,
                                  "PCPVertex_Linetension_fluctuation")),
        _rng(model.get_rng()),
        _normal_distr(0., 1.),
        _tau(get_as<double>("timescale", cfg)),
        _amplitude(get_as<double>("amplitude", cfg))
    {
        for (const auto& edge : this->_am.edges()) {
            edge->state.register_parameter(_name, 0.);
        }
    }

    ~LinetensionFluctuations () {
        for (const auto& edge : this->_am.edges()) {
            edge->state.unregister_parameter(_name);
        }
    }

    void compute_and_set_forces () final {
        for (const auto& edge : this->_am.edges()) {
            SpaceVec displ = this->_am.displacement(edge);
            SpaceVec director = displ / arma::norm(displ);

            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;

            a->state.add_force(+ get_linetension(edge) * director);
            b->state.add_force(- get_linetension(edge) * director);
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return get_linetension(edge);
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return get_linetension(edge) * this->_am.length_of(edge);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& edge : edges) {
            energy += compute_energy(edge);
        }
        return energy;
    }

    void update (double dt) final {
        for (const auto& edge : this->_am.edges()) {
            double tension = get_linetension(edge);

            double rn = _amplitude *sqrt(2. * dt / _tau) * _normal_distr(*_rng);
            tension += rn - dt / _tau * tension;

            edge->state.update_parameter(_name, {tension});
        }
    }

    void update_parameters (const DataIO::Config& cfg) final {
        _tau = get_as<double>("timescale", cfg, _tau);
        _amplitude = get_as<double>("amplitude", cfg, _amplitude);
    }


    std::vector<std::string> write_task_edge_properties_names () const final {
        return std::vector<std::string>({
            "linetension"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const final {
        std::vector<double> tensions({});
        for (const auto& edge : this->_am.edges()) {
            tensions.push_back(get_linetension(edge));
        }
        return std::vector<std::vector<double>>({tensions});
    }

    std::vector<std::string> write_task_edge_energies_names () const final {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const final {
        std::vector<double> energies({});
        std::vector<double> tensions({});
        
        for (const auto& edge : this->_am.edges()) {
            energies.push_back(compute_energy(edge));
            tensions.push_back(compute_tension(edge));
        }

        return std::vector<std::vector<double>>({
            energies,
            tensions
        });
    }
};


/// @brief The edge contractility term's base class
/** \f$ E_{i,j} = \Gamma l_{i,j}^2\f$, a term quadratic in edge length with 
 *  contractility parameter \f$ \Gamma \f$.
 */
template <typename Model>
class EdgeContractilityBase : public WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// The contractility of a particular edge
    virtual double get_contractility(const std::shared_ptr<Edge>& edge) const=0;

public:
    EdgeContractilityBase (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model)
    { }

    void compute_and_set_forces () final {
        for (const auto& edge : this->_am.edges()) {
            SpaceVec displ = this->_am.displacement(edge);

            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;

            double contract = get_contractility(edge);

            a->state.add_force(+ contract * displ);
            b->state.add_force(- contract * displ);
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return get_contractility(edge) * this->_am.length_of(edge);
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return (
              0.5 * get_contractility(edge) 
            * std::pow(this->_am.length_of(edge), 2)
        );
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& edge : edges) {
            energy += compute_energy(edge);
        }
        return energy;
    }

    std::vector<std::string> write_task_edge_properties_names () const override
    {
        return std::vector<std::string>({
            "contractility"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const override {
        std::vector<double> contractilities({});
        for (const auto& edge : this->_am.edges()) {
            contractilities.push_back(get_contractility(edge));
        }
        return std::vector<std::vector<double>>({contractilities});
    }

    std::vector<std::string> write_task_edge_energies_names () const override {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const override {
        std::vector<double> energies({});
        std::vector<double> tensions({});
        
        for (const auto& edge : this->_am.edges()) {
            energies.push_back(compute_energy(edge));
            tensions.push_back(compute_tension(edge));
        }

        return std::vector<std::vector<double>>({
            energies,
            tensions
        });
    }
};



/// @brief A uniform-homotypic edge contractility term
/** 
 *  Parameters:
 *      - `contractility`: the contractility \f$ \Gamma \f$
 */
template <typename Model>
class EdgeContractility : public EdgeContractilityBase<Model>
{
public:
    using Base = EdgeContractilityBase<Model>;

    using Edge = typename Base::Edge;

private:
    /// @brief  The contractility parameter
    double _contractility;

protected:
    inline double 
    get_contractility([[maybe_unused]] const std::shared_ptr<Edge>& edge)
    const override 
    {
        return _contractility;
    }

public:
    EdgeContractility (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _contractility(get_as<double>("contractility", cfg))
    { }

    void update_parameters (const DataIO::Config& cfg) override {
        _contractility = get_as<double>("contractility", cfg, _contractility);
    }
};

/// @brief The edge contractility term with heterotypic values
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `contractility`: the contractility \f$ \Gamma \f$. A symmetric matrix. 
 *              The i,j coordinates map to the type of cell on left and right
 *              side.
 *      - `increment_contractility` (only update): Add incremental value to 
 *              `contractility`. Cannot be update together with `contractility`.
 *      - `boundary_type`: To which value a boundary cell is mapped.
 */
template <typename Model>
class EdgeContractilityHeterotypic : public EdgeContractilityBase<Model>
{
public:
    using Base = EdgeContractilityBase<Model>;

    using Edge = typename Base::Edge;

    typedef std::vector< std::vector<double> > stdmat;

protected:
    /// The cell-type dependent contractility
    arma::mat _contractility;

    /// A cell-type index associated with the boundary
    std::size_t _boundary_type;

    double get_contractility (const std::shared_ptr<Edge>& edge) const override
    {
        const auto& [ca, cb] = this->_am.adjoints_of(edge);
        std::size_t type_a, type_b;
        if (ca) { type_a = ca->state.type; }
        else { type_a = _boundary_type; }
        if (cb) { type_b = cb->state.type; }
        else { type_b = _boundary_type; }

        if (std::max(type_a, type_b) > _contractility.n_rows)
        {
            std::cout << _contractility << std::endl;

            throw std::runtime_error(fmt::format(
                "In WF-term EdgeContractilityHeterotypic, no parameter registered "
                "for cells of type {}. Parameters for {} types registered.",
                std::max(type_a, type_b),
                _contractility.n_rows
            ));
        }

        return _contractility.at(type_a, type_b);
    }

public:
    EdgeContractilityHeterotypic (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _contractility(setup_symmetric_matrix(
            get_as<stdmat>("contractility", cfg))),
        _boundary_type(get_as<std::size_t>("boundary_type", cfg))
    { }

    void update_parameters (const DataIO::Config& cfg) override {
        if (cfg["increment_contractility"]) {
            _contractility += setup_symmetric_matrix(get_as<stdmat>(
                "increment_contractility", cfg
            ));
        }
        else if (cfg["contractility"]) {
            _contractility = setup_symmetric_matrix(get_as<stdmat>(
                "contractility", cfg
            ));
        }
        _boundary_type = get_as<std::size_t>(
            "boundary_type", cfg, _boundary_type);
    }
};

/// Heterotypic edge contractility with additional spatial gradient
/** \f$ \Gamma = \Gamma_0 + x * \Gamma_1 \f$, where x is the x-coordinate of 
 *  the edge's center relative to the half-point of the domain in x. 
 *  \Gamma_0 gives the mean contractility and \Gamma_1 the difference between
 *  hight (left) and low (right) contractility.
 *  
 * Parameters:
 *      - `contractility`: The mean contractility \f$ \Gamma_0^{i,j} \f$ 
 *              dependent on the cell types adjacent to a junction. A symmetric 
 *              matrix. The i,j coordinates map to the type of cell on left and 
 *              right side.
 *      - `gradient_contractility`: The gradient in contractility 
 *              \f$ \Gamma_1^{i,j} \f$ between left and right side of the domain
 *              dependent on the cell types adjacent to a junction. A symmetric 
 *              matrix. The i,j coordinates map to the type of cell on left and 
 *              right side. In relative units, i.e. for x = Lx, 
 *              \f$ \Gamma = \Gamma_0 - \Gamma_1 / 2 \f$ and for x = 0,
 *              \f$ \Gamma = \Gamma_0 + \Gamma_1 / 2 \f$.
 *      - `boundary_type`: To which value a boundary cell is mapped.
*/
template <typename Model>
class EdgeContractilityHeterotypicSpatialBase
: 
    public EdgeContractilityBase<Model> 
{
public:
    using Base = EdgeContractilityBase<Model>;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

    using SpaceVec = typename Base::SpaceVec;

protected:
    /// The name of the edge's contractility property
    std::string _contractility_param;

    /// A cell-type index associated with the boundary
    std::size_t _boundary_type;

    /// The origin and length (horizontal) of the domain
    std::pair<SpaceVec, SpaceVec> _domain_info;


    double get_contractility (const std::shared_ptr<Edge>& edge) const final 
    {
        if (not edge->state.has_parameter(_contractility_param)) {
            edge->state.register_parameter(_contractility_param, 0.);
            this->update_contractility(edge);
        }
        return edge->state.get_parameter(_contractility_param);
    }

    virtual double contractility_at (
        const SpaceVec& rpos,
        std::size_t type_a,
        std::size_t type_b
    ) 
    const = 0;

    /// Caluclate the origin and length of the domain
    /** Based on the vertex positions and periodicity of space.
     *  It deals with open, semi-periodic and periodic BC.
     */
    std::pair<SpaceVec, SpaceVec> get_domain_info () const {
        // Establish domain size
        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::min();
        double y_min = std::numeric_limits<double>::max();
        double y_max = std::numeric_limits<double>::min();

        const auto& space = this->_am.get_space();
        SpaceVec domain = space->get_domain_size();
        for (const auto& v : this->_am.vertices()) {
            SpaceVec pos = this->_am.position_of(v);
            if (space->periodic) {
                // for semi-periodic space
                pos[0] -= std::round(pos[0]/domain[0]) * domain[0];
            }
            x_min = std::min(x_min, pos[0]);
            x_max = std::max(x_max, pos[0]);
            y_min = std::min(y_min, pos[1]);
            y_max = std::max(y_max, pos[1]);
        }
        double Lx = (x_max - x_min);
        double Ly = (y_max - y_min);

        if (space->periodic) {
            if (Lx > 0.9 * domain[0]) {
                // is not semi-periodic space
                Lx = domain[0];
            }
            Ly = domain[1];
        }

        return std::make_pair(
            SpaceVec({x_min + Lx/2., y_min + Ly/2.}), // origin at center
            SpaceVec({Lx, Ly})  // extent
        );
    }

    /// @brief  Update the contractility parameter of an edge
    /// @param edge 
    void update_contractility (const std::shared_ptr<Edge>& edge) const { 
        SpaceVec origin = std::get<0>(_domain_info);
        SpaceVec extent = std::get<1>(_domain_info);

        SpaceVec pos = (
            this->_am.position_of(edge->custom_links().a) + 
            0.5 * this->_am.displacement(edge)
        );
        if (this->_am.get_space()->periodic) {
            SpaceVec domain = this->_am.get_space()->get_domain_size();
            if (extent[0] > 0.9 * domain[0]) {
                // periodic space
                pos -= origin;
                pos -= arma::round(pos/domain) % domain;
                pos = -(2 * arma::abs(pos) - domain / 2.);
            }
            else {
                // semi-periodic space
                pos[0] -= std::round(pos[0]/domain[0]) * domain[0];
                pos[0] -= origin[0];

                // y is periodic
                pos[1] -= origin[1];
                pos[1] -= std::round(pos[1]/domain[1]) * domain[1];
                pos[1] = -(2 * std::abs(pos[1]) - domain[1] / 2.);
            }
        }
        else {
            // non-periodic space
            pos -= origin;
        }
        SpaceVec rpos = pos / extent;

        const auto& [ca, cb] = this->_am.adjoints_of(edge);
        std::size_t type_a, type_b;
        if (ca) { type_a = ca->state.type; }
        else { type_a = this->_boundary_type; }
        if (cb) { type_b = cb->state.type; }
        else { type_b = this->_boundary_type; }

        // NOTE rpos in [-0.5, 0.5] for x and y
        edge->state.update_parameter(
            _contractility_param,
            contractility_at(rpos, type_a, type_b)
        );
    }

    /// Update all edge contractility parameters
    /** Also updates domain info (origin and length)
     */
    void update_contractility () {
        _domain_info = get_domain_info();

        // update contractility
        for (const auto& edge : this->_am.edges()) {
            update_contractility(edge);
        }
    }

public:
    EdgeContractilityHeterotypicSpatialBase (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _contractility_param(name + "__contractility"),
        _boundary_type(get_as<std::size_t>("boundary_type", cfg)),
        _domain_info(get_domain_info())
    {
        for (const auto& edge : this->_am.edges()) {
            edge->state.register_parameter(_contractility_param, 0.);
        }
    }

    ~EdgeContractilityHeterotypicSpatialBase() {
        for (const auto& edge : this->_am.edges()) {
            edge->state.unregister_parameter(_contractility_param);
        }
    }

    void update ([[maybe_unused]] double dt) override {
        update_contractility();
    }

    void update_parameters (const DataIO::Config& cfg) override {
        _boundary_type = get_as<std::size_t>(
            "boundary_type", cfg, _boundary_type
        );
        update_contractility();
    }
};

/// Heterotypic edge contractility with additional spatial gradient
/** \f$ \Gamma = \Gamma_0 + x * \Gamma_1 \f$, where x is the x-coordinate of 
 *  the edge's center relative to the half-point of the domain in x. 
 *  \Gamma_0 gives the mean contractility and \Gamma_1 the difference between
 *  hight (left) and low (right) contractility.
 *  
 * Parameters:
 *      - `contractility`: The mean contractility \f$ \Gamma_0^{i,j} \f$ 
 *              dependent on the cell types adjacent to a junction. A symmetric 
 *              matrix. The i,j coordinates map to the type of cell on left and 
 *              right side.
 *      - `increment_contractility` (only update): Add incremental value to 
 *              `contractility`. Cannot be update together with `contractility`.
 *      - `gradient_contractility`: The gradient in contractility 
 *              \f$ \Gamma_1^{i,j} \f$ between left and right side of the domain
 *              dependent on the cell types adjacent to a junction. A symmetric 
 *              matrix. The i,j coordinates map to the type of cell on left and 
 *              right side. In relative units, i.e. for x = Lx, 
 *              \f$ \Gamma = \Gamma_0 - \Gamma_1 / 2 \f$ and for x = 0,
 *              \f$ \Gamma = \Gamma_0 + \Gamma_1 / 2 \f$.
 *      - `increment_gradient_contractility` (only update): Add incremental
 *              value to `gradient_contractility`. Cannot be update together
 *              with `gradient_contractility`.
 *      - `boundary_type`: To which value a boundary cell is mapped.
*/
template <typename Model>
class EdgeContractilityHeterotypicGraded 
: 
    public EdgeContractilityHeterotypicSpatialBase<Model>
{
public:
    using Base = EdgeContractilityHeterotypicSpatialBase<Model>;

    using SpaceVec = typename Base::SpaceVec;

    typedef std::vector< std::vector<double> > stdmat;

protected:
    /// The cell-type dependent contractility
    arma::mat _contractility;

    /// The cell-type dependent spatial gradient in contractility
    arma::mat _gradient_contractility;
    
    double contractility_at(
        const SpaceVec& rpos,
        std::size_t type_a,
        std::size_t type_b
    ) 
    const override 
    {
        #if UTOPIA_DEBUG
            std::size_t N = std::min(
                this->_contractility.n_rows,
                this->_gradient_contractility.n_rows
            );
            if (std::max(type_a, type_b) > N) {
                std::cout << this->_contractility << std::endl;
                std::cout << this->_gradient_contractility << std::endl;

                throw std::runtime_error(fmt::format(
                    "In WF-term EdgeContractilityHeterotypicGraded, no "
                    "parameter registered for cells of type {}. Parameters "
                    "for {} / {} types registered.",
                    std::max(type_a, type_b),
                    this->_contractility.n_rows,
                    this->_gradient_contractility.n_rows
                ));
            }
        #endif

        return (_contractility.at(type_a, type_b) 
                + rpos[0] * _gradient_contractility.at(type_a, type_b));
    }

public:
    EdgeContractilityHeterotypicGraded (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _contractility(setup_symmetric_matrix(
            get_as<stdmat>("contractility", cfg))),
        _gradient_contractility(setup_symmetric_matrix(
            get_as<stdmat>("gradient_contractility", cfg)))
    {
        this->update_contractility();
    }

    void update_parameters (const DataIO::Config& cfg) override {
        if (cfg["increment_contractility"]) {
            _contractility += setup_symmetric_matrix(get_as<stdmat>(
                "increment_contractility", cfg
            ));
        }
        else if (cfg["contractility"]) {
            _contractility = setup_symmetric_matrix(get_as<stdmat>(
                "contractility", cfg
            ));
        }


        if (cfg["increment_gradient_contractility"]) {
            _gradient_contractility += setup_symmetric_matrix(get_as<stdmat>(
                "increment_gradient_contractility", cfg
            ));
        }
        else if (cfg["gradient_contractility"]) {
            _gradient_contractility = setup_symmetric_matrix(get_as<stdmat>(
                "gradient_contractility", cfg
            ));
        }
        
        Base::update_parameters(cfg);
    }
};

/// Heterotypic edge contractility with additional spatial gradient
/** \f$ \Gamma = \Gamma_0 + x * \Gamma_1 \f$, where x is the x-coordinate of 
 *  the edge's center relative to the half-point of the domain in x. 
 *  \Gamma_0 gives the mean contractility and \Gamma_1 the difference between
 *  hight (left) and low (right) contractility.
 *  
 * Parameters:
 *      - `contractility`: The mean contractility \f$ \Gamma_0^{i,j} \f$ 
 *              dependent on the cell types adjacent to a junction. A symmetric 
 *              matrix. The i,j coordinates map to the type of cell on left and 
 *              right side.
 *      - `gradient_contractility`: The gradient in contractility 
 *              \f$ \Gamma_1^{i,j} \f$ between left and right side of the domain
 *              dependent on the cell types adjacent to a junction. A symmetric 
 *              matrix. The i,j coordinates map to the type of cell on left and 
 *              right side. In relative units, i.e. for x = Lx, 
 *              \f$ \Gamma = \Gamma_0 - \Gamma_1 / 2 \f$ and for x = 0,
 *              \f$ \Gamma = \Gamma_0 + \Gamma_1 / 2 \f$.
 *      - `boundary_type`: To which value a boundary cell is mapped.
*/
template <typename Model>
class EdgeContractilityHeterotypicGradedMulti
: 
    public EdgeContractilityHeterotypicSpatialBase<Model> 
{
public:
    using Base = EdgeContractilityHeterotypicSpatialBase<Model>;

    using SpaceVec = typename Base::SpaceVec;

    typedef std::vector< std::vector<double> > stdmat;

protected:
    /// The cell-type dependent contractility
    std::vector<arma::mat> _contractility;

    /// The cell-type dependent spatial gradient in contractility
    std::vector<arma::mat> _gradient_contractility;
    
    double contractility_at(
        const SpaceVec& rpos,
        std::size_t type_a,
        std::size_t type_b
    ) 
    const override 
    {
        double N = _contractility.size();
        double x = rpos[0] + 0.5;
        std::size_t block = std::floor(x * N);
        block = std::min(block, static_cast<std::size_t>(N - 1));

        arma::mat contractility = _contractility[block];
        arma::mat gradient_contractility = _gradient_contractility[block];

        #if UTOPIA_DEBUG
            std::size_t n = std::min(
                contractility.n_rows,
                gradient_contractility.n_rows
            );
            if (std::max(type_a, type_b) > n) {
                std::cout << contractility << std::endl;
                std::cout << gradient_contractility << std::endl;

                throw std::runtime_error(fmt::format(
                    "In WF-term EdgeContractilityHeterotypicGraded, no "
                    "parameter registered for cells of type {}. Parameters "
                    "for {} / {} types registered.",
                    std::max(type_a, type_b),
                    contractility.n_rows,
                    gradient_contractility.n_rows
                ));
            }
        #endif
        
        // find the relative position within block
        double min = double(block) / N;
        double max = (block + 1.) / N;
        x -= min;
        x /= (max - min);
        x -= 0.5;

        return (
            contractility.at(type_a, type_b) 
            + x * gradient_contractility.at(type_a, type_b)
        );
    }

    std::vector<arma::mat> setup_matrices(const std::vector<stdmat>& vecc) {
        std::vector<arma::mat> vecm{};
        std::size_t N = 0;
        for (const stdmat& __stdmat : vecc) {
            vecm.push_back(setup_symmetric_matrix(__stdmat));
            if (vecm.size() == 1) {
                N = vecm[0].n_rows;
            }
            else if (vecm.back().n_rows != N ) {
                std::cout << "Contractility matrices:\n";
                for (const arma::mat& matrix : vecm) {
                    std::cout << matrix << std::endl;
                }
                throw std::runtime_error(fmt::format(
                    "In WF-term EdgeContractilityHeterotypicGradedMulti, "
                    "list-entry {} of contractility- or "
                    "contractility-gradient-matrices (see above) did not have "
                    "the same dimensionality as previous matrices ({} vs {})!",
                    vecm.size(),
                    N,
                    vecm.back().n_rows
                ));
            }
        }
        return vecm;
    }

    void validate_matrices() const {
        if (_contractility.size() != _gradient_contractility.size()) {
            std::cout << "Contractility matrices:\n";
            for (const arma::mat& matrix : _contractility) {
                std::cout << matrix << std::endl;
            }

            std::cout << std::endl << "Contractility-gradient matrices:\n";
            for (const arma::mat& matrix : _gradient_contractility) {
                std::cout << matrix << std::endl;
            }

            throw std::runtime_error(fmt::format(
                "In WF-term EdgeContractilityHeterotypicGradedMulti, list "
                "of contractility- and contractility-gradient-matrices "
                "(see above) did not have the same number of entries "
                "({} and {})!",
                _contractility.size(),
                _gradient_contractility.size()
            ));
        }
        if (_contractility.size() == 0) {
            throw std::runtime_error(fmt::format(
                "In WF-term EdgeContractilityHeterotypicGradedMulti, list "
                "of contractility- and contractility-gradient-matrices "
                "were empty!"
            ));
        }
    }

public:
    EdgeContractilityHeterotypicGradedMulti (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _contractility(setup_matrices(
            get_as<std::vector<stdmat>>("contractility", cfg))),
        _gradient_contractility(setup_matrices(
            get_as<std::vector<stdmat>>("gradient_contractility", cfg)))
    {
        validate_matrices();
        this->update_contractility();
    }

    void update_parameters (const DataIO::Config& cfg) override {
        if (cfg["contractility"]) {
            _contractility = setup_matrices(
                get_as<std::vector<stdmat>>("contractility", cfg)
            );
        }

        if (cfg["gradient_contractility"]) {
            _gradient_contractility = setup_matrices(
                get_as<std::vector<stdmat>>("gradient_contractility", cfg)
            );
        }

        validate_matrices();
        
        Base::update_parameters(cfg);
    }
};


/// @brief The edge contractility term to apical and basal boundary
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 *  
 * Parameters:
 *      - `apical_contractility`: the contractility \f$ \Gamma \f$ applied to 
 *              apical junctions. A vector where the entries correspond to cell 
 *              types.
 *      - `basal_contractility`: the contractility \f$ \Gamma \f$ applied to 
 *              apical junctions. A vector where the entries correspond to cell 
 *              types.
 */
template <typename Model>
class EdgeContractilityHeterotypicBoundary : public EdgeContractilityBase<Model>
{
    using Base = EdgeContractilityBase<Model>;

    using SpaceVec = typename Base::SpaceVec;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

private:
    /// Apical contractility per cell-type
    std::vector<double> _apical_contractility;

    /// Basal contractility per cell-type
    std::vector<double> _basal_contractility;

    /// The entry whether none (0) / basal (1) / apical (2)
    std::string _classification;

    double get_contractility (const std::shared_ptr<Edge>& edge) const override
    {
        const auto& [ca, cb] = this->_am.adjoints_of(edge);

        std::size_t type;
        std::shared_ptr<Cell> cell;

        if (ca) { 
            type = ca->state.type;
            cell = ca;
        }
        else {
            type = cb->state.type;
            cell = cb;
        }

        if (not edge->state.has_parameter(_classification)) {
            if (ca and cb) {
                edge->state.register_parameter(_classification, 0);
            }
            else {
                throw std::runtime_error("T1 transitions not implemented");
            }
        }

        auto classification = edge->state.get_parameter(_classification);
        if (fabs(classification - 0) < 1.e-3) {
            return 0;
        }
        else if (fabs(classification - 3) < 1.e-3) {
            return 0;
        }
        else if (fabs(classification - 1) < 1.e-3) {
            if (type >= _basal_contractility.size()) {
                throw std::runtime_error(fmt::format(
                    "In WF-term EdgeContractilityHeterotypicBoundary, no "
                    "parameter registered for cells of type {}. "
                    "Parameters for {} types registered.",
                    type,
                    _basal_contractility.size()
                ));

            }

            return _basal_contractility[cell->state.type];
        }
        else if (fabs(classification - 2) < 1.e-3) {
            if (type >= _apical_contractility.size()) {
                throw std::runtime_error(fmt::format(
                    "In WF-term EdgeContractilityHeterotypicBoundary, no "
                    "parameter registered for cells of type {}. "
                    "Parameters for {} types registered.",
                    type,
                    _apical_contractility.size()
                ));
            }

            return _apical_contractility[cell->state.type];
        }
        else {
            throw std::runtime_error(fmt::format(
                "Unknown none (0) / basal (1) / apical (2) specification of "
                "junction {}!",
                classification
            ));
        }

    }
public:
    EdgeContractilityHeterotypicBoundary (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _apical_contractility(get_as<std::vector<double>>(
            "apical_contractility", cfg)),
        _basal_contractility(get_as<std::vector<double>>(
            "basal_contractility", cfg)),
        _classification(name + "__classification")
    {
        for (const auto& edge : this->_am.edges()) {
            auto [ca, cb] = this->_am.adjoints_of(edge);

            if (ca and cb) {
                // bulk
                edge->state.register_parameter(_classification, 0);
                continue;
            }

            if (not ca) { std::swap(ca, cb); }

            SpaceVec displ({0., 1.});
            for (const auto& [e, flip] : ca->custom_links().edges) {
                if (e == edge) {
                    displ = this->_am.displacement(edge, flip);
                    break;
                }
            }

            if (fabs(atan(displ[1] / displ[0])) > 75 / 180. * M_PI) {
                // vertical
                edge->state.register_parameter(_classification, 3);
            }
            else if (displ[0] >= 0) {
                // basal junction
                edge->state.register_parameter(_classification, 1);
            }
            else {
                // apical junction
                edge->state.register_parameter(_classification, 2);
            }
        }
    }

    ~EdgeContractilityHeterotypicBoundary () {
        for (const auto& edge : this->_am.edges()) {
            edge->state.unregister_parameter(_classification);
        }
    }

    void update_parameters (const DataIO::Config& cfg) final {
        _apical_contractility = get_as<std::vector<double>>(
            "apical_contractility", cfg);
        _basal_contractility = get_as<std::vector<double>>(
            "basal_contractility", cfg);
    }

    std::vector<std::string> write_task_edge_properties_names () const override 
    {
        return std::vector<std::string>({
            "contractility",
            "classification"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const override {
        std::vector<double> contractilities({});
        std::vector<double> classifications({});
        for (const auto& edge : this->_am.edges()) {
            contractilities.push_back(get_contractility(edge));
            classifications.push_back(
                edge->state.get_parameter(_classification)
            );
        }
        return std::vector<std::vector<double>>({
            contractilities,
            classifications
        });
    }
};


/// @brief The edge contractility term with axial-asymmetric values
/** \f$ E_{i,j} = k l_{i,j}^2 (t \cdot \tau)^2\f$, a term linear in 
 *  edge length \f$ l \f$, modulated with the angle between the junctions' 
 *  tangent \f$ t \f$ and the tissue axis. The tissue axis runs parallel 
 *  to the arc if curvature non-zero.
 * 
 *  Parameters:
 *      - `contractility`: the contractility \f$ \Gamma \f$. A symmetric matrix. 
 *              The i,j coordinates map to the type of cell on left and right
 *              side. 
 *      - `boundary_type`: To which value a boundary cell is mapped.
 *      - `angle` (double): Angle of the axis to the tissue axis
 *      - `curvature` (double): The curvature of the tissue axis
 *      - `origin` (SpaceVec, optional): The center of the tissue. 
 *          The tangent of tissue axis is pointing along 
 *          x axis at origin. If not provided, the center of the tissue is used
 *          as origin.
 */
template <typename Model>
class EdgeContractilityAxial : public WorkFunctionTerm<Model>
{
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

    typedef std::vector< std::vector<double> > stdmat;

private:
    /// The matrix of contractilities for different junction types
    arma::mat _contractility;

    /// The type of a boundary cell
    /** The contractility matrix needs to include this type if non-periodic BC*/
    std::size_t _boundary_type;

    /// The axis relative to which a cell's polarity is defined
    SpaceVec _axis;

    /// The curvature of the axis.
    /** The orientation of a junction is measured with respect to an axis. This
     *  axis is running at an angle to the tangent of the arc defined by the 
     *  curvature.
     */
    double _curvature;

    /// The origin of the tissue
    /** Specifically, the position at which the tangent of the arc is horizontal
     */
    SpaceVec _origin;

    double get_contractility (const std::shared_ptr<Edge>& edge) const {
        const auto& [ca, cb] = this->_am.adjoints_of(edge);
        std::size_t type_a, type_b;
        if (ca) { type_a = ca->state.type; }
        else { type_a = _boundary_type; }
        if (cb) { type_b = cb->state.type; }
        else { type_b = _boundary_type; }

        if (std::max(type_a, type_b) > _contractility.n_rows)
        {
            std::cout << _contractility << std::endl;

            throw std::runtime_error(fmt::format(
                "In WF-term EdgeContractilityModulated, no parameter registered "
                "for cells of type {}. Parameters for {} types registered.",
                std::max(type_a, type_b),
                _contractility.n_rows
            ));
        }

        return _contractility.at(type_a, type_b);
    }

    /// @brief The local axis to an edge
    /** @param edge     The edge considered
     *  @return     Determines the local axis at the center of an edge.
     *              If no curvature, axis is same everywhere.
     *              Else, the axis is measured as an angle to the tangent
     *              of the arc.
     */
    SpaceVec local_axis (std::shared_ptr<Edge> edge) const {
        if (fabs(_curvature) < 1.e-10) {
            return _axis;
        }
        else if (this->_am.get_space()->periodic) {
            throw std::runtime_error("No curvature in periodic BC implemented");
        }

        SpaceVec pos = 0.5 * (  this->_am.position_of(edge->custom_links().a)
                              + this->_am.position_of(edge->custom_links().b));

        SpaceVec origin = _origin - SpaceVec({0., 1. / _curvature});

        SpaceVec rpos = pos - origin;
        double theta = std::atan2(rpos[0], rpos[1]);
        double alpha = std::atan2(_axis[1], _axis[0] + 1.e-8);

        return SpaceVec({cos(alpha - theta), sin(alpha - theta)});
    }


public:
    EdgeContractilityAxial (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _contractility(setup_symmetric_matrix(
            get_as<stdmat>("contractility", cfg))),
        _boundary_type(get_as<std::size_t>("boundary_type", cfg)),
        _axis({
            cos(get_as<double>("angle", cfg)),
            sin(get_as<double>("angle", cfg))
        }),
        _curvature(get_as<double>("curvature", cfg, 0.)),
        _origin(SpaceVec({0., 0.}))
    {
        if (cfg["origin"]) {
            _origin = get_as_SpaceVec<2>("origin", cfg);
        }
        else if (not model.get_space()->periodic) {
            _origin = this->_am.barycenter_of(this->_am.get_boundary_edges());
        }
        else {
            _origin = model.get_space()->get_domain_size() / 2.;
        }
    }

    SpaceVec compute_force
    ([[maybe_unused]] const std::shared_ptr<Vertex>& vertex) const final
    {
        throw std::runtime_error("Compute force on single vertex not "
            "implemented in EgdeContractilityAxial!");
    }

    void compute_and_set_forces () final {
        for (const auto& edge : this->_am.edges()) {
            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;
            
            SpaceVec axis = local_axis(edge);
            SpaceVec displ = this->_am.displacement(edge);

            double k = get_contractility(edge);
            SpaceVec dE_dx;
            if (fabs(_curvature) < 1.e-10) {
                dE_dx = - 1. * k * arma::dot(displ, axis) * axis;

                a->state.add_force(- dE_dx);
                b->state.add_force(+ dE_dx);
            }
            else if (this->_am.get_space()->periodic) {
                throw std::runtime_error("No curvature in periodic BC implemented");
            }
            else {
                // derivative of axis to displacement of vertex
                SpaceVec _a = this->_am.position_of(edge->custom_links().a);
                SpaceVec _b = this->_am.position_of(edge->custom_links().b);

                SpaceVec origin = _origin - SpaceVec({0., 1. / _curvature});
                _a -= origin;
                _b -= origin;

                double sum_x = _a[0] + _b[0];
                double sum_y = _a[1] + _b[1];
                double arg = sum_x / sum_y;

                double tmp_0 = sum_y * (pow(arg, 2) + 1);
                double tmp_1 = pow(sum_y, 2) * (pow(arg, 2) + 1);

                double dpx_dx =  axis[1] / tmp_0;
                double dpy_dx = -axis[0] / tmp_0;

                double dpx_dy = -(sum_x * axis[1]) / tmp_1;
                double dpy_dy =  (sum_x * axis[0]) / tmp_1;

                SpaceVec tmp({
                    displ[0] * dpx_dx + displ[1] * dpy_dx,
                    displ[0] * dpx_dy + displ[1] * dpy_dy,
                });

                double _dE_dx = k * arma::dot(displ, axis);

                a->state.add_force(- _dE_dx * (tmp - axis));
                b->state.add_force(- _dE_dx * (tmp + axis));
            }
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        const auto& a = edge->custom_links().a;
        const auto& b = edge->custom_links().b;

        SpaceVec axis = local_axis(edge);
        SpaceVec displ = this->_am.displacement(a, b);

        double modulation = pow(arma::dot(displ/arma::norm(displ), axis), 2);
        double length = arma::norm(displ);

        return get_contractility(edge) * modulation * length;
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        const auto& a = edge->custom_links().a;
        const auto& b = edge->custom_links().b;

        SpaceVec axis = local_axis(edge);
        SpaceVec displ = this->_am.displacement(a, b);

        // Gamma -> Gamma * cos^2 (theta), where theta angle with p-d axis
        return 0.5 * get_contractility(edge) * pow(arma::dot(displ, axis), 2);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& edge : edges) {
            energy += compute_energy(edge);
        }
        return energy;
    }

    void update_parameters (const DataIO::Config& cfg) final {
        if (cfg["contractility"]) {
            _contractility = setup_symmetric_matrix(get_as<stdmat>(
                "contractility", cfg
            ));
        }
        _boundary_type = get_as<std::size_t>(
            "boundary_type", cfg, _boundary_type);
        double angle = atan2(_axis[1], _axis[0]);
        _axis = SpaceVec({
            cos(get_as<double>("angle", cfg, angle)),
            sin(get_as<double>("angle", cfg, angle))
        });
        _curvature = get_as<double>("curvature", cfg, _curvature);
    }


    std::vector<std::string> write_task_edge_properties_names () const final {
        return std::vector<std::string>({
            "contractility"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const final {
        std::vector<double> contractilities({});
        for (const auto& edge : this->_am.edges()) {
            double k = get_contractility(edge);
            SpaceVec axis = local_axis(edge);
            SpaceVec displ = this->_am.displacement(edge);

            double modulation = pow(
                arma::dot(displ, axis) / arma::norm(displ),2);
            contractilities.push_back(k * modulation);
        }
        return std::vector<std::vector<double>>({contractilities});
    }

    std::vector<std::string> write_task_edge_energies_names () const final {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const final {
        std::vector<double> energies({});
        std::vector<double> tensions({});
        
        for (const auto& edge : this->_am.edges()) {
            energies.push_back(compute_energy(edge));
            tensions.push_back(compute_tension(edge));
        }

        return std::vector<std::vector<double>>({
            energies,
            tensions
        });
    }
};

/// @brief The edge contractility term with polar-asymmetric values
/** \f$ E_{i,j} = k l_{i,j}^2 ((x_{i+1} - x_i)/2. \cdot p)^2\f$, a term linear in 
 *  edge length \f$ l \f$, modulated with the angle between the cell's polarity
 *  vector and the vector connecting the cell center to the midpoint of an edge
 * 
 *  Parameters:
 *      - `contractility`: the contractility \f$ \Gamma \f$. A symmetric matrix. 
 *              The i,j coordinates map to the type of cell on left and right
 *              side. 
 *      - `boundary_type`: To which value a boundary cell is mapped.
 *      - `axis` (double, optional, default: 0): Angle of the tissue axis to 
 *              the tangential of an arc.
 *      - `curvature` (double): The curvature of the tissue axis
 *      - `origin` (SpaceVec, optional): The center of the tissue. 
 *          The tangent of tissue axis is pointing along 
 *          x axis at origin. If not provided, the center of the tissue is used
 *          as origin.
 *      - `setup` (str): How to setup cellular polarity.
 *          - 'uniform': uniform real distribution [0, 2 * M_PI).
 */
template <typename Model>
class EdgeContractilityPolar : public WorkFunctionTerm<Model>
{
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

    using RNG = typename Model::Base::RNG;

    typedef std::vector< std::vector<double> > stdmat;

private:
    /// The matrix of contractilities for different junction types
    arma::mat _contractility;

    /// The type of a boundary cell
    /** The contractility matrix needs to include this type if non-periodic BC*/
    std::size_t _boundary_type;

    /// The axis relative to which a cell's polarity is defined
    /** The axis is a global axis, also with curvature */
    double _axis;

    /// The curvature used to calculate the polarity_to_axis
    /** Only used in data writing */
    double _curvature;

    /// The origin of the tissue
    /** Specifically, the position at which the tangent of the arc is horizontal
     */
    SpaceVec _origin;

    /// Cell types to exclude from calculation
    /** Modulation is applied to each side of a junction addively. */
    std::set<std::size_t> _exclude_types;

    /// The damping factor for energy minimization
    double _gamma;

    /// Polarity fluctuation parameter
    /** Parameter fluctuations are implemented as Ornstein-Uhlenbeck process
     * 
     *  \f$  \frac{dp_{mn}}{dt} = - \frac{1}{\tau_p}
     *      (p_{mn}(t) - \Lambda_0)
     *      + \Delta p \sqrt{2 / \tau_p} \Theta_{mn}(t)
     *  \f$
     * 
     *  with the first value \f$ \tau \f$ and the second value
     *  \f$ \Delta p \f$.
     */
    std::pair<double, double> _fluctuations;

    const std::shared_ptr<RNG> _rng;

    /// A [0,1]-range normal distribution
    std::normal_distribution<double> _normal_distr;

    /// How many iterations to perform per update step
    std::size_t _num_steps;

    /// @brief The local axis to a position
    /** @return     Determines the local axis at the center of an edge.
     *              If no curvature, axis is same everywhere.
     *              Else, the axis is measured as an angle to the tangent
     *              of the arc.
     */
    double local_axis (SpaceVec pos) const {
        if (fabs(_curvature) < 1.e-10) {
            return _axis;
        }
        else if (this->_am.get_space()->periodic) {
            throw std::runtime_error("No curvature in periodic BC implemented");
        }

        SpaceVec origin = _origin - SpaceVec({0., 1. / _curvature});

        SpaceVec rpos = pos - origin;
        double theta = std::atan2(rpos[0], rpos[1]);

        return _axis - theta;
    }

    /// Getter for cell polarity
    double get_polarity (const std::shared_ptr<Cell>& cell) const {
        return cell->state.get_parameter(this->_name+"__polarity");
    }

    /// Setter for cell polarity
    void update_polarity (const std::shared_ptr<Cell>& cell, double value) {
        return cell->state.update_parameter(
            this->_name + "__polarity",
            constrain_angle(value)
        );
    }

    /// Initializer for cell polarity
    void register_polarity (const std::shared_ptr<Cell>& cell, double value) {
        return cell->state.register_parameter(
            this->_name + "__polarity",
            constrain_angle(value)
        );
    }

    /// Getter for heterotypic junctional contractility
    double get_contractility (const std::shared_ptr<Edge>& edge) const {
        const auto& [ca, cb] = this->_am.adjoints_of(edge);
        std::size_t type_a, type_b;
        if (ca) { type_a = ca->state.type; }
        else { type_a = _boundary_type; }
        if (cb) { type_b = cb->state.type; }
        else { type_b = _boundary_type; }

        if (std::max(type_a, type_b) > _contractility.n_rows)
        {
            std::cout << _contractility << std::endl;

            throw std::runtime_error(fmt::format(
                "In WF-term EdgeContractilityModulated, no parameter registered "
                "for cells of type {}. Parameters for {} types registered.",
                std::max(type_a, type_b),
                _contractility.n_rows
            ));
        }

        return _contractility.at(type_a, type_b);
    }

    /// Getter for modulated junctional contractility
    double get_modulated_contractility (const std::shared_ptr<Edge>& edge) const
    {
        double k = get_contractility(edge);

        if (fabs(k) < 1.e-10) {
            return 0.;
        }

        double modulation = 0.;
        const auto [adj_a, adj_b] = this->_am.adjoints_of(edge);
        for (const auto& cell : {adj_a, adj_b}) {
            if (cell == nullptr) {
                continue;
            }

            if (_exclude_types.find(cell->state.type) != _exclude_types.end()) {
                continue;
            }

            // find orientation of edge within cell
            auto e_pair = *std::find_if(
                cell->custom_links().edges.begin(),
                cell->custom_links().edges.end(),
                [edge](const auto& ep) {
                    return std::get<0>(ep) == edge;
                });

            // displ in anti-clockwise within cell
            SpaceVec displ = this->_am.displacement(edge, std::get<bool>(e_pair));
            double length = arma::norm(displ);

            double polarity = get_polarity(cell);

            // The vector to align to is 90 deg rotated to polarity of cell
            SpaceVec align_to({sin(polarity), -cos(polarity)});

            // the modulation factor
            modulation += (arma::dot(displ/length, align_to) + 1.) / 2.;
        }

        return k * modulation;
    }


public:
    EdgeContractilityPolar (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _contractility(setup_symmetric_matrix(
            get_as<stdmat>("contractility", cfg))),
        _boundary_type(get_as<std::size_t>("boundary_type", cfg)),
        _axis(get_as<double>("axis", cfg)),
        _curvature(get_as<double>("curvature", cfg, 0.)),
        _origin(SpaceVec({0., 0.})),
        _exclude_types(setup_cell_exclusion(
            get_as<std::vector<std::size_t>>("exclude_cell_types", cfg)
        )),
        _gamma(get_as<double>("gamma", cfg)),
        _fluctuations(std::make_pair(
            get_as<double>("polarity_fluctuations_tau", cfg),
            get_as<double>("polarity_fluctuations", cfg)
        )),
        _rng(model.get_rng()),
        _normal_distr(0., 1.),
        _num_steps(get_as<std::size_t>("num_steps", cfg))
    {
        if (cfg["origin"]) {
            _origin = get_as_SpaceVec<2>("origin", cfg);
        }
        else if (not model.get_space()->periodic) {
            _origin = this->_am.barycenter_of(this->_am.get_boundary_edges());
        }
        else {
            _origin = model.get_space()->get_domain_size() / 2.;
        }

        auto setup = get_as<std::string>("setup", cfg);

        if (false) {}
        else if (setup == "gaussian") {
            std::normal_distribution<double> distr(
                get_as<double>("orientation_mean", cfg),
                get_as<double>("orientation_std", cfg)
            );
            for (const auto& cell : this->_am.cells()) {
                register_polarity(cell, distr(*model.get_rng()));
            }
        }
        else if (setup == "uniform") {
            std::uniform_real_distribution distr(0., 2 * M_PI);
            for (const auto& cell : this->_am.cells()) {
                register_polarity(cell, distr(*model.get_rng()));
            }
        }
        else if (setup == "value") {
            double value = get_as<double>("orientation", cfg);
            for (const auto& cell : this->_am.cells()) {
                register_polarity(cell, value);
            }
        }
        else {
            throw std::runtime_error(fmt::format(
                "Initialization of EdgeContractilityPolar WF-term `{}` with "
                "chosen `setup`: {} not implemented. Choose one of the "
                "following: {}.",
                name,
                setup,
                "- gaussian\n"
                "- uniform\n"
                "- value\n"
            ));
        }
    }

    ~EdgeContractilityPolar () {
        for (const auto& cell : this->_am.cells()) {
            cell->state.unregister_parameter(this->_name + "__polarity");
        }
    }

private:
    auto setup_cell_exclusion (std::vector<std::size_t> exclude) const {
        std::set<std::size_t> _exclude{};
        _exclude.insert(exclude.begin(), exclude.end());
        return _exclude;
    }

public:

    SpaceVec compute_force
    ([[maybe_unused]] const std::shared_ptr<Vertex>& vertex) const final
    {
        throw std::runtime_error("Compute force on single vertex not "
            "implemented in EgdeContractilityPolar!");
    }

    void compute_and_set_forces () final {
        for (const auto& edge : this->_am.edges()) {

            const auto [adj_a, adj_b] = this->_am.adjoints_of(edge);
            for (const auto& cell : {adj_a, adj_b}) {
                if (cell == nullptr) {
                    continue;
                }

                if (_exclude_types.find(cell->state.type)!=_exclude_types.end())
                {
                    continue;
                }

                // find orientation of edge within cell
                auto e_pair = *std::find_if(
                    cell->custom_links().edges.begin(),
                    cell->custom_links().edges.end(),
                    [edge](const auto& ep) {
                        return std::get<0>(ep) == edge;
                    });
                bool flip = std::get<bool>(e_pair);
                
                auto a = edge->custom_links().a;
                auto b = edge->custom_links().b;
                if (flip) {
                    std::swap(a, b);
                }
                SpaceVec displ = this->_am.displacement(a, b);
                double length = arma::norm(displ);

                double polarity = get_polarity(cell);

                SpaceVec align_to({sin(polarity), -cos(polarity)});

                SpaceVec T1 = -(arma::dot(displ/length, align_to) + 1)/2.*displ;
                SpaceVec T2 = -0.25 * length * align_to;
                SpaceVec T3 = +0.25 * arma::dot(displ/length, align_to) * displ;

                SpaceVec dE_dx = get_contractility(edge) * (T1 + T2 + T3);

                a->state.add_force(- dE_dx);
                b->state.add_force(+ dE_dx);
            }
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        double length = this->_am.length_of(edge);
        return get_modulated_contractility(edge) * length;
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        double length = this->_am.length_of(edge);
        return 0.5 * get_modulated_contractility(edge) * pow(length, 2);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& edge : edges) {
            energy += compute_energy(edge);
        }
        return energy;
    }

    /// The update step, i.e. one minimization step
    void perform_step (double dt) {
        for (const auto& cell : this->_am.cells()) {
            if (_exclude_types.find(cell->state.type) != _exclude_types.end()) {
                continue;
            }

            double polarity = get_polarity(cell);

            double dE_dp = 0.;
            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto a = edge->custom_links().a;
                auto b = edge->custom_links().b;
                if (flip) { std::swap(a, b); }

                SpaceVec displ = this->_am.displacement(a, b);
                double length = arma::norm(displ);

                double k = get_contractility(edge);
                dE_dp += 0.25 * k * std::pow(length, 2) * 
                            arma::dot(displ,
                                      SpaceVec({cos(polarity), sin(polarity)}));
            }

            const auto& [tau, Dp] = _fluctuations;
            double rn = Dp *sqrt(2. * dt / tau) * _normal_distr(*_rng);
            
            update_polarity(cell, polarity - dE_dp * _gamma + rn);
        }
    }

    void update([[maybe_unused]] double dt) {
        for (std::size_t i = 0; i < _num_steps; i++) {
            perform_step(dt);
        }
    }

    void update_parameters (const DataIO::Config& cfg) final {
        if (cfg["contractility"]) {
            _contractility = setup_symmetric_matrix(
                get_as<stdmat>("contractility", cfg)
            );
        }
        _boundary_type = get_as<std::size_t>(
            "boundary_type", cfg, _boundary_type);
        
        _axis = get_as<double>("axis", cfg, _axis);
        _curvature = get_as<double>("curvature", cfg, _curvature);

        _gamma = get_as<double>("gamma", cfg, _gamma);

        _fluctuations = std::make_pair(
            get_as<double>("polarity_fluctuations_tau", cfg, 
                           std::get<0>(_fluctuations)),
            get_as<double>("polarity_fluctuations", cfg,
                           std::get<1>(_fluctuations))
        ),
        _num_steps = get_as<std::size_t>("num_steps", cfg, _num_steps);
    }


    std::vector<std::string> write_task_edge_properties_names () const final {
        return std::vector<std::string>({
            "contractility"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const final {
        std::vector<double> contractilities({});
        for (const auto& edge : this->_am.edges()) {
            contractilities.push_back(get_modulated_contractility(edge));
        }
        return std::vector<std::vector<double>>({contractilities});
    }

    std::vector<std::string> write_task_edge_energies_names () const final {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const final {
        std::vector<double> energies({});
        std::vector<double> tensions({});
        
        for (const auto& edge : this->_am.edges()) {
            energies.push_back(compute_energy(edge));
            tensions.push_back(compute_tension(edge));
        }

        return std::vector<std::vector<double>>({
            energies,
            tensions
        });
    }

    std::vector<std::string> write_task_cell_properties_names () const final {
        return std::vector<std::string>({"polarity", "polarity_to_axis"});
    }

    std::vector<std::vector<double>> write_cell_properties () const final {
        std::vector<double> polarities({});
        std::vector<double> polarities_to_axis({});

        for (const auto& cell : this->_am.cells()) {
            double polarity = get_polarity(cell);

            SpaceVec center = this->_am.barycenter_of(cell);
            double axis = local_axis(center);

            polarities.push_back(polarity);
            polarities_to_axis.push_back(constrain_angle(polarity - axis));
        }
        return std::vector<std::vector<double>>({
            polarities, polarities_to_axis
        });
    }
};

/// @brief The cell contractility
/** \f$ E_\alpha = 1/ 2k (P_\alpha / (\sqrt{A_\alpha^{(0)}} P^{(0)}) - 1)^2\f$, 
 *  an elastic penalty on perimeter length \f$ P_\alpha \f$.
 * 
 *  Parameters:
 *      - `contractility`: the elastic modulus \f$ k \f$
 *      - `preferential_shape`: The target perimeter \f$ P^{(0)} \f$ 
 *              at which cell is tension free.
 */
template <typename Model>
class CellContractility : public WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunctionTerm<Model>;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

    using SpaceVec = typename Base::SpaceVec;

private:
    double _contractility;

    double _preferential_shape;

public:
    CellContractility (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _contractility(get_as<double>("contractility", cfg)),
        _preferential_shape(get_as<double>("preferential_shape", cfg))
    { }

    void compute_and_set_forces () final {
        for (const auto& cell : this->_am.cells()) {
            double S0 = _preferential_shape*sqrt(cell->state.area_preferential);
            double tension = (
                  _contractility
                * (this->_am.perimeter_of(cell) / S0 - 1)
                / S0
            );

            for (const auto& [edge, flip] : cell->custom_links().edges) {
                SpaceVec director = this->_am.displacement(edge);
                director /= arma::norm(director);

                edge->custom_links().a->state.add_force(+ tension * director);
                edge->custom_links().b->state.add_force(- tension * director);
            }
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        const auto& [c1, c2] = this->_am.adjoints_of(edge);

        double T = 0.;
        if (c1) {
            double S0 = _preferential_shape*sqrt(c1->state.area_preferential);
            T += _contractility * (this->_am.perimeter_of(c1) / S0 - 1) / S0;
        }
        if (c2) {
            double S0 = _preferential_shape*sqrt(c2->state.area_preferential);
            T += _contractility * (this->_am.perimeter_of(c2) / S0 - 1) / S0;
        }

        return T;
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        double S0 = _preferential_shape * sqrt(cell->state.area_preferential);
        const double P = this->_am.perimeter_of(cell);
        return 0.5 * _contractility * std::pow(P / S0 - 1, 2);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += compute_energy(cell);
        }
        return energy;
    }


    void update_parameters (const DataIO::Config& cfg) final {
        _contractility = get_as<double>("contractility", cfg, _contractility);
        _preferential_shape = get_as<double>("preferential_shape", cfg,
                                             _preferential_shape);
    }


    std::vector<std::string> write_task_edge_energies_names () const final {
        return std::vector<std::string>({
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const final {
        std::vector<double> tensions({});
        
        for (const auto& edge : this->_am.edges()) {
            tensions.push_back(compute_tension(edge));
        }

        return std::vector<std::vector<double>>({
            tensions
        });
    }

    std::vector<std::string> write_task_cell_energies_names () const final {
        return std::vector<std::string>({
            "energy"
        });
    }

    std::vector<std::vector<double>> write_cell_energies () const final {
        std::vector<double> energies({});
        
        for (const auto& cell : this->_am.cells()) {
            energies.push_back(compute_energy(cell));
        }

        return std::vector<std::vector<double>>({
            energies
        });
    }
};

/// @brief The cell's perimeter contractility
/** \f$ E_\alpha = 1/2 k P_\alpha^2\f$, a contractility of the cell's perimeter.
 * 
 *  Parameters:
 *      - `contractility`: the elastic modulus \f$ k \f$
 */
template <typename Model>
class PerimeterContractility : public WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunctionTerm<Model>;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

    using SpaceVec = typename Base::SpaceVec;

private:
    double _contractility;

public:
    PerimeterContractility (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _contractility(get_as<double>("contractility", cfg))
    { }

    void compute_and_set_forces () final {
        for (const auto& cell : this->_am.cells()) {
            double tension = _contractility * this->_am.perimeter_of(cell);

            for (const auto& [edge, flip] : cell->custom_links().edges) {
                SpaceVec director = this->_am.displacement(edge);
                director /= arma::norm(director);

                edge->custom_links().a->state.add_force(+ tension * director);
                edge->custom_links().b->state.add_force(- tension * director);
            }
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        const auto& [c1, c2] = this->_am.adjoints_of(edge);

        double T = 0.;
        if (c1) {
            T += _contractility * this->_am.perimeter_of(c1);
        }
        if (c2) {
            T += _contractility * this->_am.perimeter_of(c2);
        }

        return T;
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        return 0.5 * _contractility * std::pow(this->_am.perimeter_of(cell), 2);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += compute_energy(cell);
        }
        return energy;
    }


    void update_parameters (const DataIO::Config& cfg) final {
        _contractility = get_as<double>("contractility", cfg, _contractility);
    }


    std::vector<std::string> write_task_edge_energies_names () const final {
        return std::vector<std::string>({
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const final {
        std::vector<double> tensions({});
        
        for (const auto& edge : this->_am.edges()) {
            tensions.push_back(compute_tension(edge));
        }

        return std::vector<std::vector<double>>({
            tensions
        });
    }

    std::vector<std::string> write_task_cell_energies_names () const final {
        return std::vector<std::string>({
            "energy"
        });
    }

    std::vector<std::vector<double>> write_cell_energies () const final {
        std::vector<double> energies({});
        
        for (const auto& cell : this->_am.cells()) {
            energies.push_back(compute_energy(cell));
        }

        return std::vector<std::vector<double>>({
            energies
        });
    }
};

/// @brief The shape elasticity of a cell
/** \f$ E_\alpha = k (P_\alpha / \sqrt(A^{(0)}) - s^{(0)})\f$, 
 *  an elastic penalty on cell shape index 
 *  \f$ s_\alpha = P_\alpha / \sqrt(A^{(0)}) \f$,
 *  where \f$ P_\alpha \f$ is the cell's perimeter.
 * 
 *  Parameters:
 *      - `contractility`: the elastic modulus \f$ k \f$
 *      - `preferential_shape_index`: The target shape index \f$ s^{(0)} \f$ 
 *              at which cell is tension free.
 */
template <typename Model>
class ShapeElasticity : public WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

private:
    double _elastic_modulus;

    double _preferential_shape_index;

public:
    ShapeElasticity (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _elastic_modulus(get_as<double>("elastic_modulus", cfg)),
        _preferential_shape_index(
            get_as<double>("preferential_shape_index", cfg))
    { }

    void compute_and_set_forces () final {
        for (const auto& cell : this->_am.cells()) {
            const double& s0 = _preferential_shape_index;
            const double area = this->_am.area_of(cell);
            const double perimeter = this->_am.perimeter_of(cell);
            const double s = perimeter / sqrt(area);
            double tension = _elastic_modulus * (s - s0) / sqrt(area);
            double pressure = (
                  _elastic_modulus * (s - s0) * perimeter
                * (-0.5 * std::pow(area, -1.5))
            );
            pressure /= 2.;

            for (const auto& [edge, flip] : cell->custom_links().edges) {
                SpaceVec displ = this->_am.displacement(edge);
                SpaceVec director = displ / arma::norm(displ);
                
                SpaceVec normal({displ[1], -displ[0]});
                if (flip) {
                    normal *= -1;
                }

                edge->custom_links().a->state.add_force(+ tension * director);
                edge->custom_links().b->state.add_force(- tension * director);
                
                edge->custom_links().a->state.add_force(- pressure * normal);
                edge->custom_links().b->state.add_force(- pressure * normal);
            }
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        const auto& [c1, c2] = this->_am.adjoints_of(edge);

        double T = 0.;
        if (c1) {
            const double& s0 = _preferential_shape_index;
            double area = this->_am.area_of(c1);
            double perimeter = this->_am.perimeter_of(c1);
            double s = perimeter / sqrt(area);
            T += _elastic_modulus * (s - s0) / sqrt(area);
        }
        if (c2) {
            const double& s0 = _preferential_shape_index;
            double area = this->_am.area_of(c2);
            double perimeter = this->_am.perimeter_of(c2);
            double s = perimeter / sqrt(area);
            T += _elastic_modulus * (s - s0) / sqrt(area);
        }

        return T;
    }

    double compute_pressure(const std::shared_ptr<Cell>& cell) const final {
        const double& s0 = _preferential_shape_index;
        const double area = this->_am.area_of(cell);
        const double perimeter = this->_am.perimeter_of(cell);
        const double s = perimeter / sqrt(area);
        return (  _elastic_modulus * (s - s0) * perimeter
                * (-0.5 * std::pow(area, -1.5)));
    }


    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        const double& s0 = _preferential_shape_index;
        const double s = this->_am.shape_index_of(cell);
        return 0.5 * _elastic_modulus * std::pow(s - s0, 2);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += compute_energy(cell);
        }
        return energy;
    }

    void update_parameters (const DataIO::Config& cfg) final {
        _elastic_modulus = get_as<double>("elastic_modulus", cfg,
                                          _elastic_modulus);
        _preferential_shape_index = get_as<double>(
            "preferential_shape_index", cfg, _preferential_shape_index);
    }


    std::vector<std::string> write_task_edge_energies_names () const final {
        return std::vector<std::string>({
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const final {
        std::vector<double> tensions({});
        
        for (const auto& edge : this->_am.edges()) {
            tensions.push_back(compute_tension(edge));
        }

        return std::vector<std::vector<double>>({
            tensions
        });
    }

    std::vector<std::string> write_task_cell_energies_names () const final {
        return std::vector<std::string>({
            "energy",
            "pressure"
        });
    }

    std::vector<std::vector<double>> write_cell_energies () const final {
        std::vector<double> energies({});
        std::vector<double> pressures({});
        
        for (const auto& cell : this->_am.cells()) {
            energies.push_back(compute_energy(cell));
            pressures.push_back(compute_pressure(cell));
        }

        return std::vector<std::vector<double>>({
            energies,
            pressures
        });
    }
};

/// @brief The base class of area elasticity of a cell
/** \f$ E_\alpha = k/2 (A_\alpha / A^{(0)} - 1)^2\f$, an elastic penalty on 
 *  cell area  \f$ A_\alpha \f$ deviating from a target value f$ A^{(0)} \f$.
 *  Requires a derived class defining the target value A0.
 * 
 *  Parameters:
 *      - `elastic_modulus`: the elastic modulus \f$ k \f$
 */
template <typename Model>
class AreaElasticityBase : public WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    double _elastic_modulus;

    virtual double preferential_area(const std::shared_ptr<Cell>& cell) const=0;

public:
    AreaElasticityBase (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _elastic_modulus(get_as<double>("elastic_modulus", cfg))
    { }

    void compute_and_set_forces () final {
        for (const auto& cell : this->_am.cells()) {
            double pressure = compute_pressure(cell) / 2.;

            for (const auto& [edge, flip] : cell->custom_links().edges) {
                SpaceVec displ = this->_am.displacement(edge);
                
                SpaceVec normal({displ[1], -displ[0]});
                if (flip) {
                    normal *= -1;
                }
                
                edge->custom_links().a->state.add_force(- pressure * normal);
                edge->custom_links().b->state.add_force(- pressure * normal);
            }
        }
    }

    double compute_pressure(const std::shared_ptr<Cell>& cell) const final {
        const double& A0 = preferential_area(cell);
        return _elastic_modulus * (this->_am.area_of(cell) / A0 - 1.) / A0;
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        const double& A0 = preferential_area(cell);
        const double A = this->_am.area_of(cell);
        return 0.5 * _elastic_modulus * std::pow(A / A0 - 1, 2);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += compute_energy(cell);
        }
        return energy;
    }

    void update_parameters (const DataIO::Config& cfg) override {
        _elastic_modulus = get_as<double>("elastic_modulus", cfg, 
                                          _elastic_modulus);
    }
    

    std::vector<std::string> write_task_cell_energies_names () const final {
        return std::vector<std::string>({
            "energy",
            "pressure"
        });
    }

    std::vector<std::vector<double>> write_cell_energies () const final {
        std::vector<double> energies({});
        std::vector<double> pressures({});
        
        for (const auto& cell : this->_am.cells()) {
            energies.push_back(compute_energy(cell));
            pressures.push_back(compute_pressure(cell));
        }

        return std::vector<std::vector<double>>({
            energies,
            pressures
        });
    }
};

/** An area elasticity with uniform target area
 *  Derived from AreaElasticityBase with a uniform target area.
 *
 *  NOTE Relies on cell state variable `area_preferential` and hence cannot be
 *      combined with other AreaElasticity terms also relying on that variable.
 *
 *  Parameters:
 *      - `preferential_area`: the target value A0.
 */
template <typename Model>
class AreaElasticity : public AreaElasticityBase<Model> {
public:
    using Base = AreaElasticityBase<Model>;

    using Cell = typename Base::Cell;

protected:
    double preferential_area(const std::shared_ptr<Cell>& cell) const final {
        return cell->state.area_preferential;
    }

public: 
    AreaElasticity (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model)
    {
        double A0 = get_as<double>("preferential_area", cfg);

        for (const auto& cell : this->_am.cells()) {
            cell->state.area_preferential = A0;
        }
    }

    void update_parameters (const DataIO::Config& cfg) override {
        if (cfg["preferential_area"]) {                
            double A0 = get_as<double>("preferential_area", cfg);

            for (const auto& cell : this->_am.cells()) {
                cell->state.area_preferential = A0;
            }
        }

        Base::update_parameters(cfg);
    }
};


/// @brief The area elasticity of a cell with heterotypic target values
/** Derived from AreaElasticityBase with a target area dependent on cell type. 
 *
 *  NOTE Relies on cell state variable `area_preferential` and hence cannot be
 *      combined with other AreaElasticity terms also relying on that variable.
 * 
 *  Parameters:
 *      - `method` (str): The method how to set the heterotypic values for
 *              preferential area \f$ A^{(0)} \f$. Available:
 *          - `set_heterotypic`: Set \f$ A^{(0)} \f$ by cell type.
 *              Additional parameters:
 *              - `preferential_areas` (vector<double>): The values for 
 *                  \f$ A^{(0)} \f$ by type.
 * 
 *          - `set_uniform`: Sets all values of \f$ A^{(0)} \f$ indepentent of 
 *              cell type.
 *              Additional parameters:
 *                  - `preferential_area` (double): The new \f$ A^{(0)} \f$ for
 *                      all cells.
 * 
 *          - `increment_heterotypic`: Increments \f$ A^{(0)} \f$ by cell type.
 *              Additional parameters:
 *              - `increments` (vector<double>): Incremental values per 
 *                      cell type.
 *              - `restrain_types` (vector<uint>): Which types are compensating
 *                      area change. If not empty, the cells of these types 
 *                      change in \f$ A^{(0)} \f$ such that total 
 *                      \f$ \sum A^{(0)}_i = const \f$.
 */
template <typename Model>
class AreaElasticityHeterotypic : public AreaElasticityBase<Model> 
{
public:
    using Base = AreaElasticityBase<Model>;

    using Cell = typename Base::Cell;

protected:
    double preferential_area(const std::shared_ptr<Cell>& cell) const final {
        return cell->state.area_preferential;
    }

public:
    AreaElasticityHeterotypic (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model)
    {
        this->setup_preferential_area(cfg);
    }


private:
    void setup_preferential_area(const Config& cfg) {
        const auto method = get_as<std::string>("method", cfg);
        if (false) { }

        else if (method == "increment_heterotypic") {
            auto _increments = get_as<std::vector<double>>("increments", cfg);
            auto restrain = get_as<std::vector<std::size_t>>("restrain_types",
                                                             cfg);

            auto increments(_increments);
            if (restrain.size() > 0) {
                std::size_t types = 0;
                std::size_t cnt_restrain = 0;
                double area_change = 0.;
                for (const auto& cell : this->_am.cells()) {
                    types = std::max(types, cell->state.type);
                    area_change += increments[cell->state.type];
                    auto find_it = std::find(restrain.begin(),
                                            restrain.end(),
                                            cell->state.type);
                    if (find_it != restrain.end())
                    {
                        cnt_restrain++;
                    }
                }

                if (cnt_restrain == 0 and fabs(area_change) > 1.e-8) {
                    throw std::runtime_error(fmt::format("Heterotypic "
                        "increment of preferential area with restrain on types "
                        "failed, because no cells with that type!"));
                }

                area_change /= cnt_restrain;
                for (const auto& type : restrain) {
                    increments[type] -= area_change;
                }
            }

            for (const auto& cell : this->_am.cells()) {
                if (cell->state.type >= increments.size()) {
                    throw std::runtime_error(fmt::format("In update of WF "
                        "AreaElasticityHeterogeneous, no incremental area "
                        "value defined for cell of type {}.",
                        cell->state.type));
                }
                cell->state.area_preferential += increments[cell->state.type];
            }
        }

        else if (method == "set_heterotypic") {
            auto A0 = get_as<std::vector<double>>("preferential_areas", cfg);
            std::size_t types = 0;
            for (const auto& cell : this->_am.cells()) {
                types = std::max(types, cell->state.type);
            }

            if (A0.size() <= types) {
                throw std::runtime_error(fmt::format("For heterotypic "
                    "preferential area, one value for every type must be "
                    "provided, but got {} values for {} types!",
                    A0.size(), types+1));
            }

            for (const auto& cell : this->_am.cells()) {
                cell->state.area_preferential = A0[cell->state.type];
            }
        }

        else if (method == "set_uniform") {
            double A0 = get_as<double>("preferential_area", cfg);

            for (const auto& cell : this->_am.cells()) {
                cell->state.area_preferential = A0;
            }
        }

        else {
            throw std::runtime_error(fmt::format("Unknown method '{}' to "
                "set preferential area in AreaElasticity! Choose one of the "
                "following:"
                "- increment_heterotypic\n"
                "- set_heterotypic\n"
                "- set_uniform\n",
                method
            ));
        }
    }


public:
    /// @brief  Updates the parameters
    /// @param cfg The configuration forwarded to 
    ///            AreaElasticityHeterotypic::setup_preferential_area
    void update_parameters (const DataIO::Config& cfg) override {
        setup_preferential_area(cfg);
        Base::update_parameters(cfg);
    }
};


/// @brief The area elasticity with target values that follow a gradient
/** Derived from AreaElasticityBase with a target area dependent on cell type.
 *  The values of \f$ A^{(0)}\f$ follow a sinusoidal gradient. Over the length
 *  of the domain, \f$ A^{(0)}\f$ changes from `left` to `right` and back to 
 *  `left` with a \$ \cos(x)^2 \f$ profile.
 *
 *  NOTE Relies on cell state variable `area_preferential` and hence cannot be
 *      combined with other AreaElasticity terms also relying on that variable.
 * 
 *  Parameters:
 *      - `relaxation_time` (double): Over which timescale the area of non-
 *              regulated types relaxes to keep total area constant.
 *      - `relax_types` (list[int]): The cell types that are relaxed
 *      - `num_types` (int): The number of total cell types
 *      - `min_area` (double): The minium preferential area of a cell
 *      - `manhatten_distance` (int): The distance over which cells adapt.
 *          Distance is given as manhatten-distance, i.e. direct neighbours
 *          have distance 1. See also EntitiesManager::neighbors_of.
 */
template <typename Model>
class AreaElasticityGradient : public AreaElasticityBase<Model> 
{
public:
    using SpaceVec = typename Model::SpaceVec;

    using Base = AreaElasticityBase<Model>;

    using Cell = typename Base::Cell;

protected:
    double preferential_area(const std::shared_ptr<Cell>& cell) const final {
        return cell->state.area_preferential;
    }

private:
    /// The maximum and minimum value of the cos^2 profile
    std::vector<std::pair<double, double>> _regulate_gradient;

    /// Relative length of the plateau inserted at the minimum of cos^2
    double _plateau_length;

    /// List of types that are not regulated and are relaxed instead
    std::unordered_set<std::size_t> _relax_types;

    /// Manhatten distance in relaxation
    std::size_t _distance;

    /// Minimum area of a cell
    double _min_area;

    /// Time scale of relaxation
    double _tau;

    /// The targeted average area
    double _reference_area;

    const enum Shape {
        grad,
        cos,
        cos2
    } _shape;

    /// Update the preferential area of regulated cells
    void regulate_gradient () {
        std::function<double(double, double, double)> A0;
        if (_shape == Shape::grad) {
            A0 = [this](double x, double left, double right) {
                x /= 1. - this->_plateau_length;
                x = std::min(x, 1.);
                return left + (right - left) * x;
            };
        }
        else if (_shape == Shape::cos) {
            A0 = [this](double x, double left, double right) {
                x /= 1. - _plateau_length;
                x = std::min(x, 1.);
                return (left - right) * (1 + std::cos(x * M_PI))/ 2 + right;
            };
        }
        else if (_shape == Shape::cos2) {
            A0 = [this](double x, double left, double right) {
                x = std::min(x, 1.-x);
                x /= 1 - _plateau_length / 2.;
                x = std::min(x, 0.5);
                return (left - right) * std::pow(std::cos(x * M_PI), 2) + right;
            };
        }
        else {
            throw std::runtime_error(fmt::format(
                "Not Implemented shape ({}) in AreaElasticityGradient!",
                _shape
            ));
        }
        

        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::min();

        const auto& space = this->_am.get_space();
        double __Lx = space->get_domain_size()[0];
        for (const auto& v : this->_am.vertices()) {
            double x = this->_am.position_of(v)[0];
            if (space->periodic) {
                x -= std::round(x/__Lx) * __Lx;
            }
            x_min = std::min(x_min, x);
            x_max = std::max(x_max, x);
        }
        double L = (x_max - x_min);


        for (const auto& cell : this->_am.cells()) {
            if (cell->state.type >= _regulate_gradient.size()
                or std::find(
                        _relax_types.begin(), _relax_types.end(), 
                        cell->state.type
                   ) != _relax_types.end()
            )
            {
                continue;
            }

            double x = this->_am.barycenter_of(cell)[0];
            if (space->periodic) {
                x -= std::round(x/__Lx) * __Lx;
            }
            x -= x_min;
            x /= L;

            const auto& [left, right] = _regulate_gradient[cell->state.type];
            cell->state.area_preferential = A0(x, left, right);
        }
    }

    /// Update the preferential area of relaxing cells
    void relax_cells () {
        if (_relax_types.size() == 0) {
            return;
        }

        std::map<std::shared_ptr<Cell>, double> targets({});
        for (const auto& cell : this->_am.cells()) {
            auto find = std::find(_relax_types.begin(), _relax_types.end(),
                                  cell->state.type);
            if (find == _relax_types.end()) {
                continue;
            }

            const auto& neighbors = this->_am.neighbors_of(cell, _distance);
            
            std::size_t N_tot = neighbors.size() + 1;
            std::size_t N = 1;
            double average_area_tot = cell->state.area_preferential;
            double average_area = cell->state.area_preferential;
            for (const auto& n : neighbors) {
                average_area_tot += n->state.area_preferential;
                
                auto find = std::find(_relax_types.begin(), _relax_types.end(),
                                      n->state.type);
                if (find != _relax_types.end()) {
                    average_area += n->state.area_preferential;
                    N++;
                }
            }
            average_area_tot /= static_cast<double>(N_tot);
            average_area /= static_cast<double>(N);

            auto A = cell->state.area_preferential;
            double dA = (average_area_tot/_reference_area - 1.) * A;
            double A_target = A - dA / _tau;
            A_target += (average_area - A) / _tau;

            targets[cell] = std::max(A_target, _min_area);
        }

        for (const auto& [cell, A_target] : targets) {
            double dA = A_target - cell->state.area_preferential;
            cell->state.area_preferential += dA;
        }
    }

public:
    AreaElasticityGradient (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _regulate_gradient({}),
        _plateau_length(get_as<double>("plateau_length", cfg)),
        _relax_types({}),
        _distance(get_as<std::size_t>("manhatten_distance", cfg)),
        _min_area(get_as<double>("min_area", cfg)),
        _tau(get_as<double>("relaxation_time", cfg)),
        _reference_area(get_as<double>("preferential_area", cfg)),
        _shape(setup_shape(cfg))
    {
        for (std::size_t i = 0; i < get_as<double>("num_types", cfg); i++) {
            _regulate_gradient.push_back(std::make_pair(_reference_area,
                                                        _reference_area));
        }

        if (_plateau_length > 1 or _plateau_length < 0) {
            throw std::runtime_error(fmt::format("In AreaElasticityGradient, "
                "parameter `plateau_length` must be in [0, 1], but was {}", 
                _plateau_length));
        }

        auto relax_types = get_as<std::vector<std::size_t>>("relax_types", cfg);
        _relax_types = std::unordered_set<std::size_t>(relax_types.begin(),
                                                       relax_types.end());
    }

private:
    Shape setup_shape(const Config& cfg){
        auto shape = get_as<std::string>("shape", cfg);
        if (shape == "grad" or shape == "gradient") {
            return Shape::grad;
        }
        if (shape == "cos" or shape == "cosinus") {
            return Shape::cos;
        }
        if (shape == "cos2" or shape == "cosinus2" or shape == "cosinus_square")
        {
            return Shape::cos2;
        }
        throw std::runtime_error(fmt::format(
            "While setting up the AreaElasticityGradient work-function ('{}'), "
            "got unknown shape '{}'. Implemented shapes are: \n"
            " - grad / gradient (a linear gradient)\n",
            " - cos / cosinus (a cosinus shape, discontinuous in periodic space)\n",
            " - cos2 (a cos^2(x) shape, continuous in periodic space)\n",
            this->_name,
            shape
        ));
    }


public:

    /// @brief  Updates the parameters
    /// @param cfg The configuration forwarded to 
    ///            AreaElasticityHeterotypic::update_preferential_area
    /** Additional parameters in cfg:
     *      - `increment_areas` (list[(double, double)]): The incremental values
     *          for left and right values of gradient per cell type.
     * 
    */
    void update_parameters (const DataIO::Config& cfg) override {
        if (cfg["area_elasticity"]) {
            DataIO::Config _cfg{};
            _cfg["area_elasticity"] = cfg["area_elasticity"];
            Base::update_parameters(_cfg);
        }

        _plateau_length = get_as<double>("plateau_length", cfg,
                                              _plateau_length);
        if (_plateau_length > 1 or _plateau_length < 0) {
            throw std::runtime_error(fmt::format("In AreaElasticityGradient, "
                "parameter `plateau_length` must be in [0, 1], but was {}", 
                _plateau_length));
        }

        std::size_t num = get_as<double>("num_types", cfg,
                                         _regulate_gradient.size());
        for (std::size_t i = _regulate_gradient.size(); i < num; i++) {
            _regulate_gradient.push_back(std::make_pair(_reference_area,
                                                        _reference_area));
        }
        
        if (cfg["increment_areas"]) {
            auto increments = get_as<std::vector<std::pair<double, double>>>(
                "increment_areas", cfg);

            if (increments.size() > num) {
                throw std::runtime_error(fmt::format("Too many values provided "
                    "in update_parameters of WF AreaElasticityGradient. "
                    "Provided {} values, but `num_types` is {}!",
                    increments.size(), num));
            }

            for (std::size_t i = 0; i < increments.size(); i++) {
                const auto& [left, right] = increments[i];
                const auto& [_left, _right] = _regulate_gradient[i];
                _regulate_gradient[i] = std::make_pair(
                    _left + left,
                    _right + right
                );
            }

            regulate_gradient();
        }
        if (cfg["relax_types"]) {
            auto relax_types = get_as<std::vector<std::size_t>>("relax_types",
                                                                cfg);

            _relax_types = std::unordered_set<std::size_t>(relax_types.begin(),
                                                           relax_types.end());
        }

        Base::update_parameters(cfg);
    }

    void update ([[maybe_unused]] double dt) override {
        regulate_gradient();
        relax_cells();
    }
};


/// @brief The area elasticity of a cell to prevent extrusion
/** An area elasticity (see AreaElasticityBase), where A0 = A_min if A < A_min
 *  and A otherwise. It therefore only is non-zero for A < A_min, preventing
 *  further reduction of A.
 *
 *  NOTE This derived class does NOT rely on cell state variable 
 *      `area_preferential` and can be combined with other AreaElasticity terms.
 * 
 *  Parameters:
 *      - `area_minimum`: The minimal area, below which area elasticity is 
 *              applied.
 */
template <typename Model>
class AreaElasticityMinimum : public AreaElasticityBase<Model> 
{
public:
    using Base = AreaElasticityBase<Model>;

    using Cell = typename Base::Cell;

protected:
    double _minimum_area;

    double preferential_area(const std::shared_ptr<Cell>& cell) const final {
        double A = this->_am.area_of(cell);
        if (A < _minimum_area) {
            return _minimum_area;
        }
        return A;
    }

public:
    AreaElasticityMinimum (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _minimum_area(get_as<double>("minimum_area", cfg))
    { }

    
    /// @brief  Updates the parameters
    /// @param cfg The configuration forwarded to 
    ///            AreaElasticityHeterotypic::setup_preferential_area
    void update_parameters (const DataIO::Config& cfg) override {
        _minimum_area = get_as<double>("minimum_area", cfg, _minimum_area);

        Base::update_parameters(cfg);
    }
};

} // namespace WorkFunction
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
