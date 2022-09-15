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
    WorkFunctionEdgeTerm (const DataIO::Config& cfg, const Model& model)
    :
        _am(model.get_am())
    { }

    virtual double compute_tension(const std::shared_ptr<Edge>& edge) const = 0;

    virtual double compute_energy(const std::shared_ptr<Edge>& edge) const = 0;
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
    WorkFunctionCellTerm (const DataIO::Config& cfg, const Model& model)
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
    const double linetension;

public:
    Linetension (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        linetension(get_as<double>("linetension", cfg))
    { }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return linetension;
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return linetension * this->_am.length_of(edge);
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
    const double contractility;

public:
    EdgeContractility (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        contractility(get_as<double>("contractility", cfg))
    { }

    double compute_tension(const std::shared_ptr<Edge>& edge) const final {
        return contractility * this->_am.length_of(edge);
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return 0.5 * contractility * std::pow(this->_am.length_of(edge), 2);
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
    const double _contractility;

    const double _preferential_perimeter;

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
    const double _elastic_modulus;

    const double _preferential_shape_index;

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
};

/// @brief The area elasticity of a cell
/** \f$ E_\alpha = k (A_\alpha / A^{(0)} - 1)\f$, an elastic penalty on 
 *  cell area  \f$ A_\alpha \f$..
 * 
 *  Parameters:
 *      - `contractility`: the elastic modulus \f$ k \f$
 *      - `preferential_perimeter`: The target area \f$ A^{(0)} \f$ 
 *              at which cell is pressure free.
 */
template <typename Model>
class AreaElasticity : public WorkFunctionCellTerm<Model>
{
    using Base = WorkFunctionCellTerm<Model>;

    using Cell = typename Base::Cell;

private:
    const double _elastic_modulus;

    const double _preferential_area;

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
        const double& A0 = _elastic_modulus;
        return _elastic_modulus * (this->_am.area_of(cell) / A0 - 1.) / A0;
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        const double& A0 = _elastic_modulus;
        const double A = this->_am.area_of(cell);
        return 0.5 * _elastic_modulus * std::pow(A / A0 - 1, 2);
    }
};


} // namespace WorkFunction
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
