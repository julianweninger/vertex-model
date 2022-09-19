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
class WorkFunctionEdgeTerm {
public:
    using AgentManager = typename Model::AgentManager;
    using Edge = typename Model::Edge;

protected:
    const AgentManager& _am;

public:
    WorkFunctionEdgeTerm ([[maybe_unused]] const DataIO::Config& cfg,
                          const Model& model)
    :
        _am(model.get_am())
    { }

    virtual double compute_tension(const std::shared_ptr<Edge>& edge) const = 0;

    virtual double compute_energy(const std::shared_ptr<Edge>& edge) const = 0;

    virtual void update ([[maybe_unused]] double dt) { return; }

    virtual void update_parameters (const DataIO::Config& cfg) = 0;
};

/// @brief The class of a pressure term in the work function
/// @tparam Model   The model within which the AgentManager lives
/** A pressure is the derivative of the energy to the area of a cell. 
 */
template <typename Model>
class WorkFunctionCellTerm {
public:
    using AgentManager = typename Model::AgentManager;
    using Cell = typename Model::Cell;
    using Edge = typename Model::Edge;


protected:
    const AgentManager& _am;
    
public:
    WorkFunctionCellTerm ([[maybe_unused]] const DataIO::Config& cfg,
                          const Model& model)
    :
        _am(model.get_am())
    { }

    virtual double compute_tension(const std::shared_ptr<Cell>& cell) const=0;

    double compute_tension(const std::shared_ptr<Edge>& edge) const {
        const auto& [cl, cr] = this->_am.template adjoints_of<true>(edge);

        double T = 0.;
        if (cl) {
            T += compute_tension(cl);
        }
        if (cr) {
            T += compute_tension(cr);
        }
        return T;
    };

    virtual double compute_pressure(const std::shared_ptr<Cell>& cell) const=0;

    virtual double compute_energy(const std::shared_ptr<Cell>& cell) const = 0;

    virtual void update ([[maybe_unused]] double dt) { return; }

    virtual void update_parameters (const DataIO::Config& cfg) = 0;
};


/// @brief The linetension term
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `linetension`: the tension \f$ k \f$
 */
template <typename Model>
class Linetension : public WorkFunctionEdgeTerm<Model>
{
    using Base = WorkFunctionEdgeTerm<Model>;

    using Edge = typename Base::Edge;

private:
    double _linetension;

public:
    Linetension (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        _linetension(get_as<double>("linetension", cfg))
    { }

    double compute_tension([[maybe_unused]] const std::shared_ptr<Edge>& edge)
    const final
    {
        return _linetension;
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return _linetension * this->_am.length_of(edge);
    }

    void update_parameters (const DataIO::Config& cfg) final {
        _linetension  = get_as<double>("linetension", cfg, _linetension);
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

/// @brief The linetension term
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `linetension`: the tension \f$ k \f$. A symmetric matrix. 
 *              The i,j coordinates map to the type of cell on left and right
 *              side. 
 *      - `boundary_type`: To which value a boundary cell is mapped.
 */
template <typename Model>
class LinetensionHeterotypic : public WorkFunctionEdgeTerm<Model>
{
    using Base = WorkFunctionEdgeTerm<Model>;

    using Edge = typename Base::Edge;

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

        if (std::max(type_a, type_b) > _linetension.index_max())
        {
            std::cout << _linetension << std::endl;

            throw std::runtime_error(fmt::format(
                "In WF-term LinetensionHeterotypic, no parameter registered "
                "for cells of type {}. Parameters for {} types registered.",
                std::max(type_a, type_b),
                _linetension.index_max()
            ));
        }

        return _linetension.at(type_a, type_b);
    }

public:
    LinetensionHeterotypic (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        _linetension(setup_symmetric_matrix(
            get_as<stdmat>("linetension", cfg))),
        _boundary_type(get_as<std::size_t>("boundary_type", cfg))
    { }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return get_linetension(edge);
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return get_linetension(edge) * this->_am.length_of(edge);
    }

    void update_parameters (const DataIO::Config& cfg) final {
        stdmat tension = arma::conv_to<stdmat>::from(_linetension);
        _linetension = setup_symmetric_matrix(get_as<stdmat>(
            "linetension", cfg, tension
        ));
        _boundary_type = get_as<std::size_t>(
            "boundary_type", cfg, _boundary_type);
    }
};



/// @brief The linetension term
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `linetension`: the tension \f$ k \f$. A symmetric matrix. 
 *              The i,j coordinates map to the type of cell on left and right
 *              side. 
 *      - `boundary_type`: To which value a boundary cell is mapped.
 */
template <typename Model>
class LinetensionFluctuations : public WorkFunctionEdgeTerm<Model>
{
    using Base = WorkFunctionEdgeTerm<Model>;

    using Edge = typename Base::Edge;

    using RNG = typename Model::Base::RNG;

private:
    const std::string _name;

    const std::shared_ptr<RNG> _rng;

    /// A [0,1]-range normal distribution
    std::normal_distribution<double> _normal_distr;

    double _tau;

    double _amplitude;

    const double& get_linetension (const std::shared_ptr<Edge>& edge) const {
        if (not edge->state.has_parameter(_name)) {
            edge->state.register_parameter(_name, {0.});
        }
        return edge->state.get_parameter(_name)[0];
    }

public:
    LinetensionFluctuations (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        _name(get_as<std::string>("name", cfg,
                                  "PCPVertex_Linetension_fluctuation")),
        _rng(model.get_rng()),
        _normal_distr(0., 1.),
        _tau(get_as<double>("timescale", cfg)),
        _amplitude(get_as<double>("amplitude", cfg))
    {
        for (const auto& edge : this->_am.edges()) {
            edge->state.register_parameter(_name, {0.});
        }
    }

    ~LinetensionFluctuations () {
        for (const auto& edge : this->_am.edges()) {
            edge->state.unregister_parameter(_name);
        }
    }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return get_linetension(edge);
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return get_linetension(edge) * this->_am.length_of(edge);
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
};

/// @brief The linetension term
/** \f$ E_{i,j} = k l_{i,j}^2\f$, a term quadratic in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `contractility`: the contractility \f$ k \f$
 */
template <typename Model>
class EdgeContractility : public WorkFunctionEdgeTerm<Model>
{
    using Base = WorkFunctionEdgeTerm<Model>;

    using Edge = typename Base::Edge;

private:
    double _contractility;

public:
    EdgeContractility (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        _contractility(get_as<double>("contractility", cfg))
    { }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return _contractility * this->_am.length_of(edge);
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return 0.5 * _contractility * std::pow(this->_am.length_of(edge), 2);
    }

    void update_parameters (const DataIO::Config& cfg) {
        _contractility = get_as<double>("contractility", cfg, _contractility);

    }
};

/// @brief The cell contractility
/** \f$ E_\alpha = k (P_\alpha / P^{(0)} - 1)\f$, an elastic penalty on 
 *  perimeter length \f$ P_\alpha \f$.
 * 
 *  Parameters:
 *      - `contractility`: the elastic modulus \f$ k \f$
 *      - `preferential_perimeter`: The target perimeter \f$ P^{(0)} \f$ 
 *              at which cell is tension free.
 */
template <typename Model>
class CellContractility : public WorkFunctionCellTerm<Model>
{
    using Base = WorkFunctionCellTerm<Model>;

    using Cell = typename Base::Cell;
    using Edge = typename Base::Edge;

private:
    double _contractility;

    double _preferential_perimeter;

public:
    CellContractility (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        _contractility(get_as<double>("contractility", cfg)),
        _preferential_perimeter(get_as<double>("preferential_perimeter", cfg))
    { }

    double compute_tension(const std::shared_ptr<Cell>& cell) const final {
        const double& P0 = _preferential_perimeter;
        return _contractility * (this->_am.perimeter_of(cell) / P0 - 1) / P0;
    }

    double compute_pressure([[maybe_unused]] const std::shared_ptr<Cell>& cell)
    const final
    {
        return 0.;
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        const double& P0 = _preferential_perimeter;
        const double P = this->_am.perimeter_of(cell);
        return 0.5 * _contractility * std::pow(P / P0 - 1, 2);
    }

    void update_parameters (const DataIO::Config& cfg) {
        _contractility = get_as<double>("contractility", cfg, _contractility);
        _preferential_perimeter = get_as<double>("preferential_perimeter", cfg,
                                                 _preferential_perimeter);
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
 *      - `preferential_perimeter`: The target perimeter \f$ P^{(0)} \f$ 
 *              at which cell is tension free.
 */
template <typename Model>
class ShapeElasticity : public WorkFunctionCellTerm<Model>
{
    using Base = WorkFunctionCellTerm<Model>;

    using Cell = typename Base::Cell;

private:
    double _elastic_modulus;

    double _preferential_shape_index;

public:
    ShapeElasticity (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        _elastic_modulus(get_as<double>("elastic_modulus", cfg)),
        _preferential_shape_index(
            get_as<double>("preferential_shape_index", cfg))
    { }

    double compute_tension(const std::shared_ptr<Cell>& cell) const final {
        const double& s0 = _preferential_shape_index;
        const double area = this->_am.area_of(cell);
        const double perimeter = this->_am.perimeter_of(cell);
        const double s = perimeter / sqrt(area);
        return _elastic_modulus * (s - s0) / sqrt(area);
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

    void update_parameters (const DataIO::Config& cfg) {
        _elastic_modulus = get_as<double>("elastic_modulus", cfg,
                                          _elastic_modulus);
        _preferential_shape_index = get_as<double>(
            "preferential_shape_index", cfg, _preferential_shape_index);
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
class AreaElasticity : public WorkFunctionCellTerm<Model>
{
    using Base = WorkFunctionCellTerm<Model>;

    using Cell = typename Base::Cell;

private:
    double _elastic_modulus;

    double _preferential_area;

public:
    AreaElasticity (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        _elastic_modulus(get_as<double>("elastic_modulus", cfg)),
        _preferential_area(get_as<double>("preferential_area", cfg))
    { }

    double compute_tension([[maybe_unused]] const std::shared_ptr<Cell>& cell)
    const final 
    {
        return 0.;
    }

    double compute_pressure(const std::shared_ptr<Cell>& cell) const final {
        const double& A0 = _preferential_area;
        return _elastic_modulus * (this->_am.area_of(cell) / A0 - 1.) / A0;
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        const double& A0 = _preferential_area;
        const double A = this->_am.area_of(cell);
        return 0.5 * _elastic_modulus * std::pow(A / A0 - 1, 2);
    }

    void update_parameters (const DataIO::Config& cfg) {
        _elastic_modulus = get_as<double>("elastic_modulus", cfg, 
                                          _elastic_modulus);
        _preferential_area = get_as<double>("preferential_area", cfg, 
                                            _preferential_area);
    }
};

/// @brief The area elasticity of a cell
/** \f$ E_\alpha = k (A_\alpha / A^{(0)}_\alpha - 1)\f$, an elastic penalty on 
 *  cell area \f$ A_\alpha \f$ wrt target value \f$ A^{(0)}_\alpha \f$.
 * 
 *  Parameters:
 *      - `contractility`: the elastic modulus \f$ k \f$
 *      - `preferential_area`: The target area \f$ A^{(0)} \f$ 
 *              at which cell is pressure free. Positions map to cell type.
 */
template <typename Model>
class AreaElasticityHeterotypic : public WorkFunctionCellTerm<Model>
{
    using Base = WorkFunctionCellTerm<Model>;

    using Cell = typename Base::Cell;

private:
    double _elastic_modulus;

    std::vector<double> _preferential_area;

    const auto& get_preferential_area (const std::shared_ptr<Cell>& cell) const
    {
        if (cell->state.type >= _preferential_area.size()) {
            std::cout << _preferential_area.size();
            throw std::runtime_error(fmt::format(
                "In WF-term AreaElasticityHeterotypic, no parameter registered "
                "for cells of type {}. Parameters for {} types registered.",
                cell->state.type, _preferential_area.size()
            ));
        }
        return _preferential_area.at(cell->state.type);
    }

public:
    AreaElasticityHeterotypic (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        _elastic_modulus(get_as<double>("elastic_modulus", cfg)),
        _preferential_area(
            get_as<std::vector<double>>("preferential_area", cfg))
    { }

    double compute_tension([[maybe_unused]] const std::shared_ptr<Cell>& cell)
    const final 
    {
        return 0.;
    }

    double compute_pressure(const std::shared_ptr<Cell>& cell) const final {
        const double& A0 = get_preferential_area(cell);
        return _elastic_modulus * (this->_am.area_of(cell) / A0 - 1.) / A0;
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        const double& k = _elastic_modulus;
        const double& A0 = get_preferential_area(cell);
        return 0.5 * k * std::pow(this->_am.area_of(cell) / A0 - 1, 2);
    }

    void update_parameters (const DataIO::Config& cfg) {
        _elastic_modulus = get_as<double>("elastic_modulus", cfg,
                                          _elastic_modulus);
        _preferential_area = get_as<std::vector<double>>("preferential_area",
                                                         cfg,
                                                         _preferential_area);
    }
};

/// @brief The area elasticity of a cell
/** \f$ E_\alpha = k (A_\alpha / A^{(0)}_\alpha - 1)\f$, an elastic penalty on 
 *  cell area \f$ A_\alpha \f$ wrt target value \f$ A^{(0)}_\alpha \f$.
 * 
 *  Parameters:
 *      - `contractility`: the elastic modulus \f$ k \f$
 *      - `preferential_area`: The target area \f$ A^{(0)} \f$ 
 *              at which cell is pressure free. Positions map to cell type.
 */
template <typename Model>
class AreaElasticityIndividual : public WorkFunctionCellTerm<Model>
{
    using Base = WorkFunctionCellTerm<Model>;

    using Cell = typename Base::Cell;

private:
    const std::string _name;

public:
    AreaElasticityIndividual (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        _name(get_as<std::string>("name", cfg, "AreaElasticityIndividual"))
    {
        double A0 = get_as<double>("area_preferential", cfg);
        double k = get_as<double>("elastic_modulus", cfg);
        for (const auto& cell : this->_am.cells()) {
            cell->state.register_parameter(_name, {k, A0});
        }
    }

    ~AreaElasticityIndividual() {
        for (const auto& cell : this->_am.cells()) {
            cell->state.unregister_parameter(_name);
        }
    }

    double compute_tension([[maybe_unused]] const std::shared_ptr<Cell>& cell)
    const final 
    {
        return 0.;
    }

    double compute_pressure(const std::shared_ptr<Cell>& cell) const final {
        const auto& params = get_parameters(cell);
        return params[0] * (this->_am.area_of(cell)/params[1] - 1.) / params[1];
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        const auto& params = get_parameters(cell);
        return 0.5*params[0]*std::pow(this->_am.area_of(cell)/params[1] - 1, 2);
    }
    
    const auto& get_parameters (const std::shared_ptr<Cell>& cell) const
    {
        if (not cell->state.has_parameter(_name)) {
            throw std::runtime_error(fmt::format("Cell {} has no parameter "
                "`{}` registered. Required from work function term `{}` "
                "(AreaElasticityIndividual)!",
                cell->id(), _name, _name));
        }
        return cell->state.get_parameter(_name);
    }

    void update_parameters (const DataIO::Config& cfg) override {
        return;
    }

    void set_parameters (const std::shared_ptr<Cell>& cell, 
                         double k, double A0)
    {
        const auto& params = get_parameters(cell);
        cell->state.update_parameter(_name, {k, A0});
    }
};


} // namespace WorkFunction
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
