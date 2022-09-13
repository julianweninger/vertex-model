#ifndef UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION_HH
#define UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION_HH


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction {

template <class Model>
std::pair<typename Model::Tension,
          typename Model::TensionEnergy> setup_linetension \
(const DataIO::Config& cfg, const Model& model)
{
    using Edge = typename Model::Edge; 
    using Tension = typename Model::Tension; 
    using TensionEnergy = typename Model::TensionEnergy;

    double k = get_as<double>("linetension", cfg);
    const auto& am = model.get_am();

    Tension linetension = \
    [k]
    ([[maybe_unused]] const std::shared_ptr<Edge>& edge)
    {
        return k;
    };
    TensionEnergy linetension_E = \
    [k, am]
    (const std::shared_ptr<Edge>& edge)
    {
        return k * am.length_of(edge);
    };

    return std::make_pair(linetension, linetension_E);
}

template <class Model>
std::pair<typename Model::Tension,
          typename Model::TensionEnergy>  setup_edge_contractility
(const DataIO::Config& cfg, const Model& model)
{
    using Edge = typename Model::Edge; 
    using Tension = typename Model::Tension; 
    using TensionEnergy = typename Model::TensionEnergy;

    double k = get_as<double>("contractility", cfg);
    const auto& am = model.get_am();

    Tension contractility = [k, am](const std::shared_ptr<Edge>& edge) {
        return k * am.length_of(edge);
    };
    TensionEnergy contractility_E = \
    [k, am]
    (const std::shared_ptr<Edge>& edge)
    {
        return 0.5 * k * std::pow(am.length_of(edge), 2);
    };

    return std::make_pair(contractility, contractility_E);
}

template <class Model>
std::pair<typename Model::TensionCell,
          typename Model::TensionCellEnergy>  setup_cell_contractility
(const Config& cfg, const Model& model)
{
    using Cell = typename Model::Cell; 
    using TensionCell = typename Model::TensionCell; 
    using TensionCellEnergy = typename Model::TensionCellEnergy;

    double k = get_as<double>("contractility", cfg);
    double P0 = get_as<double>("preferential_perimeter", cfg);
    const auto& am = model.get_am();

    TensionCell contractility = \
    [k, P0, am]
    (const std::shared_ptr<Cell>& cell)
    {
        return k * (am.perimeter_of(cell) - P0) / P0;
    };
    TensionCellEnergy contractility_E = \
    [k, P0, am]
    (const std::shared_ptr<Cell>& cell)
    {
        return 0.5 * k * std::pow(am.perimeter_of(cell) / P0 - 1, 2);
    };

    return std::make_pair(contractility, contractility_E);
}

template <class Model>
std::pair<typename Model::Pressure,
          typename Model::PressureEnergy>  setup_area_elasticity
(const Config& cfg, const Model& model)
{
    using Cell = typename Model::Cell; 
    using Pressure = typename Model::Pressure; 
    using PressureEnergy = typename Model::PressureEnergy;
    
    double k = get_as<double>("elastic_modulus", cfg);
    double A0 = get_as<double>("preferential_area", cfg);
    const auto& am = model.get_am();

    Pressure area_elasticity = \
    [k, A0, am]
    (const std::shared_ptr<Cell>& cell)
    {
        return k * (am.area_of(cell) / A0 - 1.) / A0;
    };

    PressureEnergy area_elasticity_E = \
    [k, A0, am]
    (const std::shared_ptr<Cell>& cell)
    {
        return 0.5 * k * std::pow(am.area_of(cell) / A0 - 1, 2);
    };

    return std::make_pair(area_elasticity, area_elasticity_E);
}


} // namespace WorkFunction
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
