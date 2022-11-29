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
    using SpaceVec = typename Model::SpaceVec;

    using AgentManager = typename Model::AgentManager;

    using Vertex = typename Model::Vertex;

    using Edge = typename Model::Edge;

    using Cell = typename Model::Cell;

protected:
    const std::string _name;

    const AgentManager& _am;

public:
    WorkFunctionTerm (std::string name,
                      [[maybe_unused]] const DataIO::Config& cfg,
                      const Model& model)
    :
        _name(name),
        _am(model.get_am())
    { }

    virtual ~WorkFunctionTerm() { }


    virtual void compute_and_set_forces () = 0;


    virtual SpaceVec compute_force
    ([[maybe_unused]] const std::shared_ptr<Vertex>& vertex) const
    {
        return SpaceVec({0., 0.});
    }

    virtual double compute_tension
    ([[maybe_unused]] const std::shared_ptr<Edge>& edge) const
    {
        return 0.;
    }

    virtual double compute_pressure
    ([[maybe_unused]] const std::shared_ptr<Cell>& cell) const
    {
        return 0.;
    }
    

    virtual double compute_energy
    ([[maybe_unused]] const std::shared_ptr<Vertex>& vertex) const 
    {
        return 0.;
    }
    
    virtual double compute_energy
    ([[maybe_unused]] const std::shared_ptr<Edge>& edge) const 
    {
        return 0.;
    }
    
    virtual double compute_energy
    ([[maybe_unused]] const std::shared_ptr<Cell>& cell) const 
    {
        return 0.;
    }

    virtual double compute_energy (
        const AgentContainer<Vertex>& vertices,
        const AgentContainer<Edge>& edges,
        const AgentContainer<Cell>& cells
    ) const = 0;

    double compute_energy () const {
        return compute_energy(
            this->_am.vertices(),
            this->_am.edges(),
            this->_am.cells()
        );
    }


    virtual void update ([[maybe_unused]] double dt) { return; }

    virtual void update_parameters (const DataIO::Config& cfg) = 0;


    const std::string& get_name () const {
        return _name;
    }


    virtual std::vector<std::string> write_task_vertex_properties_names () const {
        return std::vector<std::string>({});
    }

    virtual std::vector<std::vector<double>> write_vertex_properties () const {
        return std::vector<std::vector<double>>({});
    }

    virtual std::vector<std::string> write_task_edge_properties_names () const {
        return std::vector<std::string>({});
    }

    virtual std::vector<std::vector<double>> write_edge_properties () const {
        return std::vector<std::vector<double>>({});
    }

    virtual std::vector<std::string> write_task_cell_properties_names () const {
        return std::vector<std::string>({});
    }

    virtual std::vector<std::vector<double>> write_cell_properties () const {
        return std::vector<std::vector<double>>({});
    }


    virtual std::vector<std::string> write_task_vertex_energies_names () const {
        return std::vector<std::string>({});
    }

    virtual std::vector<std::vector<double>> write_vertex_energies () const {
        return std::vector<std::vector<double>>({});
    }

    virtual std::vector<std::string> write_task_edge_energies_names () const {
        return std::vector<std::string>({});
    }

    virtual std::vector<std::vector<double>> write_edge_energies () const {
        return std::vector<std::vector<double>>({});
    }

    virtual std::vector<std::string> write_task_cell_energies_names () const {
        return std::vector<std::string>({});
    }

    virtual std::vector<std::vector<double>> write_cell_energies () const {
        return std::vector<std::vector<double>>({});
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



    std::vector<std::string> write_task_edge_energies_names () const {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const {
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
    const std::vector<std::vector<double>>& values) 
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
        stdmat tension(_linetension.n_rows);
        for (size_t i = 0; i < _linetension.n_rows; ++i) {
            tension[i] = arma::conv_to< std::vector<double> >::from(
                _linetension.row(i)
            );
        };
        _linetension = setup_symmetric_matrix(get_as<stdmat>(
            "linetension", cfg, tension
        ));
        _boundary_type = get_as<std::size_t>(
            "boundary_type", cfg, _boundary_type);
    }



    std::vector<std::string> write_task_edge_properties_names () const {
        return std::vector<std::string>({
            "linetension"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const {
        std::vector<double> tensions({});
        for (const auto& edge : this->_am.edges()) {
            tensions.push_back(get_linetension(edge));
        }
        return std::vector<std::vector<double>>({tensions});
    }

    std::vector<std::string> write_task_edge_energies_names () const {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const {
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

    void update_parameters (const DataIO::Config& cfg) {
        _tau = get_as<double>("timescale", cfg, _tau);
        _amplitude = get_as<double>("amplitude", cfg, _amplitude);
    }


    std::vector<std::string> write_task_edge_properties_names () const {
        return std::vector<std::string>({
            "linetension"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const {
        std::vector<double> tensions({});
        for (const auto& edge : this->_am.edges()) {
            tensions.push_back(get_linetension(edge));
        }
        return std::vector<std::vector<double>>({tensions});
    }

    std::vector<std::string> write_task_edge_energies_names () const {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const {
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


/// @brief The edge contractility term
/** \f$ E_{i,j} = k l_{i,j}^2\f$, a term quadratic in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `contractility`: the contractility \f$ k \f$
 */
template <typename Model>
class EdgeContractility : public WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;
private:
    double _contractility;

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


    void compute_and_set_forces () final {
        for (const auto& edge : this->_am.edges()) {
            SpaceVec displ = this->_am.displacement(edge);

            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;

            a->state.add_force(+ _contractility * displ);
            b->state.add_force(- _contractility * displ);
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return _contractility * this->_am.length_of(edge);
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return 0.5 * _contractility * std::pow(this->_am.length_of(edge), 2);
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

    void update_parameters (const DataIO::Config& cfg) {
        _contractility = get_as<double>("contractility", cfg, _contractility);

    }


    std::vector<std::string> write_task_edge_energies_names () const {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const {
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

/// @brief The edge contractility term with heterotypic values
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `contractility`: the contractility \f$ \Gamma \f$. A symmetric matrix. 
 *              The i,j coordinates map to the type of cell on left and right
 *              side. 
 *      - `boundary_type`: To which value a boundary cell is mapped.
 */
template <typename Model>
class EdgeContractilityHeterotypic : public WorkFunctionTerm<Model>
{
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

    typedef std::vector< std::vector<double> > stdmat;

private:
    arma::mat _contractility;

    std::size_t _boundary_type;

    const double& get_contractility (const std::shared_ptr<Edge>& edge) const {
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

    void compute_and_set_forces () final {
        for (const auto& edge : this->_am.edges()) {
            SpaceVec displ = this->_am.displacement(edge);

            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;

            a->state.add_force(+ get_contractility(edge) * displ);
            b->state.add_force(- get_contractility(edge) * displ);
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return get_contractility(edge) * this->_am.length_of(edge);
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return (
            0.5
            * get_contractility(edge)
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

    void update_parameters (const DataIO::Config& cfg) final {
        stdmat tension(_contractility.n_rows);
        for (size_t i = 0; i < _contractility.n_rows; ++i) {
            tension[i] = arma::conv_to< std::vector<double> >::from(
                _contractility.row(i)
            );
        };
        _contractility = setup_symmetric_matrix(get_as<stdmat>(
            "contractility", cfg, tension
        ));
        _boundary_type = get_as<std::size_t>(
            "boundary_type", cfg, _boundary_type);
    }




    std::vector<std::string> write_task_edge_properties_names () const {
        return std::vector<std::string>({
            "contractility"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const {
        std::vector<double> contractilities({});
        for (const auto& edge : this->_am.edges()) {
            contractilities.push_back(get_contractility(edge));
        }
        return std::vector<std::vector<double>>({contractilities});
    }

    std::vector<std::string> write_task_edge_energies_names () const {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const {
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

/// @brief The edge contractility term with heterotypic values
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `contractility`: the contractility \f$ \Gamma \f$. A symmetric matrix. 
 *              The i,j coordinates map to the type of cell on left and right
 *              side. 
 *      - `boundary_type`: To which value a boundary cell is mapped.
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
    arma::mat _contractility;

    std::size_t _boundary_type;

    double _angle;

    SpaceVec _axis;

    double _curvature;

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
        else {
            _origin = this->_am.barycenter_of(this->_am.get_boundary_edges());
        }
    }

    void compute_and_set_forces () final {
        for (const auto& edge : this->_am.edges()) {
            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;
            
            SpaceVec axis = _axis;
            SpaceVec displ = this->_am.displacement(edge);
            
            // rotate axis angle wrt tangential of circle
            if (fabs(_curvature) > 1.e-10) {
                SpaceVec pos = (  this->_am.position_of(edge->custom_links().a)
                                + 0.5 * displ);

                SpaceVec origin = _origin - SpaceVec({0., -1. / _curvature});

                SpaceVec displ = pos - origin;
                double theta = std::atan2(displ[0], displ[1]);

                axis = SpaceVec({
                    axis[0]*cos(-theta) - axis[1]*sin(-theta),
                    axis[0]*sin(-theta) + axis[1]*cos(-theta)
                });
            }

            double k = get_contractility(edge);
            SpaceVec dE_dx = k * arma::dot(displ, axis) * axis;

            a->state.add_force(- dE_dx);
            b->state.add_force(+ dE_dx);
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return get_contractility(edge) * this->_am.length_of(edge);
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        const auto& a = edge->custom_links().a;
        const auto& b = edge->custom_links().b;

        SpaceVec axis = _axis;
        SpaceVec displ = this->_am.displacement(a, b);
        
        // rotate ppMLC_axis so that points along curved tissue axis
        if (fabs(_curvature) > 1.e-10) {
            SpaceVec pos = (  this->_am.position_of(edge->custom_links().a)
                            + 0.5 * displ);

            SpaceVec origin = _origin - SpaceVec({0., -1. / _curvature});

            SpaceVec displ = pos - origin;
            double theta = std::atan2(displ[0], displ[1]);

            axis = SpaceVec({
                axis[0]*cos(-theta) - axis[1]*sin(-theta),
                axis[0]*sin(-theta) + axis[1]*cos(-theta)
            });
        }

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
        stdmat tension = arma::conv_to<stdmat>::from(_contractility);
        _contractility = setup_symmetric_matrix(get_as<stdmat>(
            "contractility", cfg, tension
        ));
        _boundary_type = get_as<std::size_t>(
            "boundary_type", cfg, _boundary_type);
        double angle = atan2(_axis[1], _axis[0]);
        _axis = SpaceVec({
            cos(get_as<double>("angle", cfg, angle)),
            sin(get_as<double>("angle", cfg, angle))
        });
        _curvature = get_as<double>("curvature", cfg, _curvature);
    }


    std::vector<std::string> write_task_edge_properties_names () const {
        return std::vector<std::string>({
            "contractility"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const {
        std::vector<double> contractilities({});
        for (const auto& edge : this->_am.edges()) {
            contractilities.push_back(get_contractility(edge));
        }
        return std::vector<std::vector<double>>({contractilities});
    }

    std::vector<std::string> write_task_edge_energies_names () const {
        return std::vector<std::string>({
            "energy",
            "tension"
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const {
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

/// @brief The cell contractility
/** \f$ E_\alpha = k (P_\alpha / (\sqrt{A_\alpha^{(0)}} P^{(0)}) - 1)^2\f$, 
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

    void update_parameters (const DataIO::Config& cfg) {
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

    void update_parameters (const DataIO::Config& cfg) {
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

/// @brief The area elasticity of a cell
/** \f$ E_\alpha = k (A_\alpha / A^{(0)} - 1)\f$, an elastic penalty on 
 *  cell area  \f$ A_\alpha \f$..
 * 
 *  Parameters:
 *      - `contractility`: the elastic modulus \f$ k \f$
 *      - `preferential_area`: The target area \f$ A^{(0)} \f$ 
 *              at which cell is pressure free.
 */
template <typename Model>
class AreaElasticity : public WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    double _elastic_modulus;

public:
    AreaElasticity (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _elastic_modulus(get_as<double>("elastic_modulus", cfg))
    {
        double A0 = get_as<double>("preferential_area", cfg);

        for (const auto& cell : this->_am.cells()) {
            cell->state.area_preferential = A0;
        }
    }

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
        const double& A0 = cell->state.area_preferential;
        return _elastic_modulus * (this->_am.area_of(cell) / A0 - 1.) / A0;
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        const double& A0 = cell->state.area_preferential;
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

        if (cfg["preferential_area"]) {                
            double A0 = get_as<double>("preferential_area", cfg);

            for (const auto& cell : this->_am.cells()) {
                cell->state.area_preferential = A0;
            }
        }
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


template <typename Model>
class AreaElasticityHeterotypic : public AreaElasticity<Model> 
{
public:
    using Base = AreaElasticity<Model>;

    using Cell = typename Base::Cell;

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
                    throw std::runtime_error(fmt::format("Heterotypic increment of "
                        "preferential area with restrain on types failed, "
                        "because no cells with that type!"));
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
                "following", method,
                "- increment_heterotypic\n"
                "- set_heterotypic\n"
                "- set_uniform\n"
            ));
        }
    }


public:
    void update_parameters (const DataIO::Config& cfg) override {
        if (cfg["area_elasticity"]) {
            DataIO::Config _cfg{};
            _cfg["area_elasticit"] = cfg["area_elasticity"];
            Base::update_parameters(_cfg);
        }
        
        setup_preferential_area(cfg);
    }
};

} // namespace WorkFunction
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
