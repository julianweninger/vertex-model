#ifndef UTOPIA_MODELS_PCPTOPOLOGY_OPERATIONS_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_OPERATIONS_HH

#include <typeinfo>

#include <utopia/core/types.hh>

#include "../PCPVertex/PCPVertex.hh"

namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace OperationCollection {

/// Parameters on how and when to apply operations
struct OperationParams {
    using Time = typename PCPVertex::Time;

    /// Name of the operation
    const std::string name;

    /// The times when to invoke
    /** \note This are the times of the model,
     *        hence at iteration 1 the time is 0
     */
    std::set<Time> times;

    /// The number of iterations of the operations
    /** \details default: 1
     */
    std::size_t iterations;

    /// The interval of logging
    /** \details default: 0
     */
    std::size_t emit_interval;

    /// The propability to invoke
    /** \details Default: 1; If evaluated false, all iterations are skipped.
     */
    double probability;

    /// The number of iterations during prolog
    /** \details Default: 0
     *  \note no probability applied here
     */
    std::size_t iterations_prolog;

    /// The number of iterations during epilog
    /** \details Default: 0
     *  \note no probability applied here
     */
    std::size_t iterations_epilog;

    /// When to minimize the energy
    /** \details default: every
     */
    enum MinimizationMode {
        Every,  /// After every iteration
        Once,   /// After all iterations at specific time
        Manual  /// Do not automatically minimize
    } minimization_mode;

    /// The parameters for energy minimization
    /** \details Can be updated by passing 'minimization' dict
     */
    MinimizationParams minimization_params;

    /// Initialize from cfg
    /** \details minimization_params are updated from cfg
     */
    template <typename Config>
    OperationParams (std::string name, const Config& cfg,
                     const MinimizationParams& default_minim_params)
    :
        name(name),
        times(),
        iterations(get_as<std::size_t>("iterations", cfg, 1)),
        emit_interval( get_as<std::size_t>("emit_interval", cfg, 0)),
        probability(get_as<double>("probability", cfg, 1)),
        iterations_prolog(get_as<std::size_t>("iterations_prolog", cfg, 0)),
        iterations_epilog(get_as<std::size_t>("iterations_epilog", cfg, 0)),
        minimization_mode(setup_minimization_mode(get_as<std::string>(
            "mode", get_as<Config>("minimization", cfg, Config()), "every"))),
        minimization_params((get_as<Config>("minimization", cfg, Config()),
                            default_minim_params))
    {
        auto times_list = get_as<std::vector<Time>>("times", cfg);
        // TODO Consider wrapping negative values around

        // Make sure negative times are not included
        times_list.erase(
            std::remove_if(times_list.begin(), times_list.end(),
                            [](auto& t){ return (t < 0); }),
            times_list.end()
        );

        // Populate the set; this will impose ordering
        times.insert(times_list.begin(), times_list.end());
    }

    MinimizationMode setup_minimization_mode(std::string mode) {
        if (mode == "every") {
            return MinimizationMode::Every;
        }
        else if (mode == "once") {
            return MinimizationMode::Once;
        }
        else if (mode == "manual" or mode == "off") {
            return MinimizationMode::Manual;
        }
        else {
            throw std::invalid_argument(fmt::format(
                "Unknown mode for minimization. Was {}. Known modes are: {}",
                mode, "'every', 'once', 'manual', 'off'"));
        }
    }
};

using Operation = std::function<void(PCPVertex& vertex_model)>;

using OperationBundle = typename std::pair<Operation, OperationParams>;

/// The operation to differentiate the types of cells randomly
/** Progenitor cells turn to hair cell with given probability and to support
 *  cell otherwise.
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `probability` (double): The probability for a cell to differentiate
 *              to a hair cell; otherwise support cell. Gives a fraction of 
 *              hair cells for large enough number of cells.
 */
OperationBundle build_differentiate_random (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);
    
    double probability(get_as<double>("probability", cfg));
    
    if (params.iterations_prolog + params.iterations_epilog +
        params.iterations * (params.times.size()) > 1)
    {
        throw std::invalid_argument(fmt::format(
            "Differentiate can be applied only once, because terminal "
            "process. Iterations in prolog: {}, in run {}, in epilog",
            params.iterations_prolog,
            params.iterations * (params.times.size()),
            params.iterations_epilog));
    }

    std::uniform_real_distribution<double> prob_distr;

    Operation operation = [probability, prob_distr{std::move(prob_distr)}]
            (PCPVertex& vertex_model) mutable
    {
        using CellType = PCPVertex::CellType;

        PCPVertex::RuleFuncCell update = [vertex_model, probability,
                                          prob_distr{std::move(prob_distr)}]
                (const auto& cell) mutable
        {
            auto state = cell->state;
            if (prob_distr(*vertex_model.get_rng()) < probability) {
                state.type = CellType::hair;
            }
            else {
                state.type = CellType::support;
            }
            return state;
        };

        apply_rule<Update::sync>(update, vertex_model.get_am().cells());
    };

    return std::make_pair(operation, params);
}

/// The operation to increment
/** The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `steps` (uint): The numer of steps performed in the NotchDelta model.
 * 
 *  \param notch_delta  Pointer to a notch delta model
 */
template <class NotchDelta>
OperationBundle build_differentiate_NotchDelta (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params,
        std::shared_ptr<NotchDelta> notch_delta)
{
    OperationParams params(name, cfg, default_minim_params);

    std::size_t steps(get_as<std::size_t>("steps", cfg));
    
    if (params.iterations_prolog + params.iterations_epilog +
        params.iterations * (params.times.size()) > 1)
    {
        throw std::invalid_argument(fmt::format(
            "Differentiate can be applied only once, because terminal "
            "process. Iterations in prolog: {}, in run {}, in epilog",
            params.iterations_prolog,
            params.iterations * (params.times.size()),
            params.iterations_epilog));
    }

    Operation operation = [notch_delta, steps] (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;
        using NDCellType = typename NotchDelta::CellType;

        const auto& nd_cells = notch_delta->get_cm().cells();
        const auto& cells = vertex_model.get_am().cells();

        std::unordered_map<std::shared_ptr<PCPVertex::Cell>,
                           std::shared_ptr<typename NotchDelta::Cell>> cell_map;
        if (cells.size() > nd_cells.size()) {
            throw std::runtime_error(fmt::format("Cannot link cells of "
                "NotchDelta and Vertex models because ore cells in Vertex "
                "({} cells) than in NotchDelta model ({} cells)!",
                nd_cells.size(), cells.size()));
        }
        unsigned int iterator;
        cell_map.reserve(cells.size());
        for (iterator = 0; iterator < cells.size(); iterator++) {
            cell_map.insert({cells[iterator], nd_cells[iterator]});
        }
        for (void(); iterator < nd_cells.size(); iterator++) {
            nd_cells[iterator]->state.cell_type = NDCellType::inactive;
            nd_cells[iterator]->custom_links().neighbors.clear();
        }

        for (const auto& c : cells) {
            auto mapped_cell = cell_map.at(c);
            mapped_cell->custom_links().neighbors.clear();
            for (auto n : vertex_model.get_am().neighbors_of(c)) {
                mapped_cell->custom_links().neighbors.push_back(cell_map.at(n));
            }
        }

        notch_delta->prolog();
        for (std::size_t i = 0; i < steps; i++) {
            notch_delta->iterate();
        }
        notch_delta->epilog();

        for (const auto [cell, nd_cell] : cell_map) {
            auto type = nd_cell->state.cell_type;
            if (type == NDCellType::hair) {
                cell->state.type = CellType::hair;
            }
            else if (type == NDCellType::support) {
                cell->state.type = CellType::support;
            }
            else {
                cell->state.type = CellType::progenitor;
            }
        }
    };

    return std::make_pair(operation, params);
}

/// The operation to increment preferential area
/** \details The following parameter are extracted from cfg 
 *           (besides those passed to `OperationParams`):
 *               - `progenitor` (double, default: 0.): incremental value for
 *                      cells of type progenitor
 *               - `var_progenitor` (double, default: 0.): variance for cells
 *                      of type progenitor; using normal distribution
 *               - `hair` (double, default: 0.): incremental value for cells
 *                      of type hair
 *               - `var_hair` (double, default: 0.): variance for cells
 *                      of type hair; using normal distribution
 *               - `support` (double, default: 0.): incremental value for cells
 *                      of type support
 *               - `var_support` (double, default: 0.): variance for cells
 *                      of type support; using normal distribution
 *               - `adapt_support` (bool, default: false): If true, the total
 *                      area of hair and support cells remains constant, hence
 *                      support cells compensate area changes from hair cells.
 */
OperationBundle build_increment_area (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    double incr_prog(get_as<double>("progenitor", cfg, 0.));
    double var_prog(get_as<double>("var_progenitor", cfg, 0.));
    double incr_hair(get_as<double>("hair", cfg, 0.));
    double var_hair(get_as<double>("var_hair", cfg, 0.));
    double incr_support(get_as<double>("support", cfg, 0.));
    double var_support(get_as<double>("var_support", cfg, 0.));
    bool adapt_support(get_as<bool>("adapt_support", cfg, false));
    
    if (adapt_support and incr_support != 0) {
        throw std::invalid_argument(fmt::format(
            "In cfg {}: If `adapt_support = True`, "
            "then the increment value for the support cells "
            "`support` must be zero, but was {}!", name, incr_support));
    }

    Operation operation = [incr_prog, var_prog,
                           incr_hair, var_hair,
                           incr_support, var_support, adapt_support]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        std::normal_distribution<> dist_prog{incr_prog, var_prog};
        std::normal_distribution<> dist_hair{incr_hair, var_hair};
        std::normal_distribution<> dist_support;
        if (adapt_support) {
            double incr_support;
            std::size_t num_hcs = std::count_if(
                cells.begin(), cells.end(), 
                [](const auto& c) {
                    return c->state.type == CellType::hair;
                });

            if (num_hcs < cells.size()) {
                incr_support = -1. * incr_hair * num_hcs /
                               (cells.size() - num_hcs);
            }
            else {
                incr_support = 0;
            }
            
            dist_support = std::normal_distribution<>(incr_support,
                                                      var_support);
        }
        else {
            dist_support = std::normal_distribution<>(incr_support,
                                                      var_support);
        }

        PCPVertex::RuleFuncCell update = [vertex_model,
                                          dist_prog{std::move(dist_prog)},
                                          dist_hair{std::move(dist_hair)},
                                          dist_support{std::move(dist_support)}]
                (const auto& cell) mutable
        {
            auto state = cell->state;
            if (state.type == CellType::progenitor) {
                state.area_preferential += dist_prog(*vertex_model.get_rng());
            }
            else if (state.type == CellType::support) {
                state.area_preferential += dist_support(*vertex_model.get_rng());
            }
            else if (state.type == CellType::hair) {
                state.area_preferential += dist_hair(*vertex_model.get_rng());
            }
            return state;
        };
        
        apply_rule<Update::sync>(update, cells);
    };

    return std::make_pair(operation, params);
}

/// The operation to increment the domain size
/** The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `value` (SpaceVec<2>): The incremental increase in x and y
 *      - `compensate` (bool): Whether to compensate the increase in area in the
 *              preferential area of the cells. If true, keeps the mechanical
 *              properties constant
 *      - `fix_hair_cell_volume` (bool, default: false): Whether to fix the area
 *              of cells of type hair.
 */
OperationBundle build_increment_domain (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    using SpaceVec = PCPVertex::SpaceVec;

    OperationParams params(name, cfg, default_minim_params);
    
    SpaceVec increment(get_as_SpaceVec<2>("value", cfg));
    bool compensate(get_as<bool>("compensate", cfg));
    bool fix_hair_cell_volume(get_as<bool>("fix_hair_cell_volume", cfg, false));

    Operation operation = [increment, compensate,
                           fix_hair_cell_volume] (PCPVertex& vertex_model)
    {
        vertex_model.stretch_domain(increment, compensate,
                                    fix_hair_cell_volume);
    };

    return std::make_pair(operation, params);
}

/// The operation to increment edge contractility
/** The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `progenitor_progenitor` (double, default: 0.): Incremental value for 
 *              edges between cells of types progenitor and progenitor.
 *      - `progenitor_hair` (double, default: 0.): Incremental value for 
 *              edges between cells of types progenitor and hair.
 *      - `progenitor_support` (double, default: 0.): Incremental value for 
 *              edges between cells of types progenitor and support.
 *      - `hair_hair` (double, default: 0.): Incremental value for 
 *              edges between cells of types hair and hair.
 *      - `hair_support` (double, default: 0.): Incremental value for 
 *              edges between cells of types hair and support.
 *      - `support_support` (double, default: 0.): Incremental value for 
 *              edges between cells of types support and support.
 */
OperationBundle build_increment_edge_contractility (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);
    
    double incr_prog_prog(get_as<double>("progenitor_progenitor", cfg, 0.));
    double incr_prog_hair(get_as<double>("progenitor_hair", cfg, 0.));
    double incr_prog_supp(get_as<double>("progenitor_support", cfg, 0.));
    double incr_hair_hair(get_as<double>("hair_hair", cfg, 0.));
    double incr_hair_supp(get_as<double>("hair_support", cfg, 0.));
    double incr_supp_supp(get_as<double>("support_support", cfg, 0.));

    if (cfg["hair_progenitor"]) {
        double incr_hair_prog = get_as<double>("hair_progenitor", cfg);
        if (cfg["progenitor_hair"] and incr_hair_prog != incr_prog_hair) {
            throw KeyError("", cfg, fmt::format("Increment has to be "
                "symmetric matrix."));
        }
        incr_prog_hair = incr_hair_prog;
    }
    if (cfg["support_progenitor"]) {
        double incr_supp_prog = get_as<double>("hair_progenitor", cfg);
        if (cfg["progenitor_support"] and incr_supp_prog != incr_prog_supp) {
            throw KeyError("", cfg, fmt::format("Increment has to be "
                "symmetric matrix."));
        }
        incr_prog_supp = incr_supp_prog;
    }
    if (cfg["support_hair"]) {
        double incr_supp_hair = get_as<double>("hair_progenitor", cfg);
        if (cfg["hair_support"] and incr_supp_hair != incr_hair_supp) {
            throw KeyError("", cfg, fmt::format("Increment has to be "
                "symmetric matrix."));
        }
        incr_hair_supp = incr_supp_hair;
    }

    Operation operation = [incr_prog_prog, incr_prog_hair, incr_prog_supp,
                           incr_hair_hair, incr_hair_supp, incr_supp_supp]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        auto contractility = vertex_model.get_edge_contractility();

        contractility(CellType::progenitor,
                    CellType::progenitor) += incr_prog_prog;
        contractility(CellType::progenitor, CellType::hair) += incr_prog_hair;
        contractility(CellType::progenitor,
                      CellType::support) += incr_prog_supp;
        
        contractility(CellType::hair, CellType::progenitor) += incr_prog_hair;
        contractility(CellType::hair, CellType::hair) += incr_hair_hair;
        contractility(CellType::hair, CellType::support) += incr_hair_supp;
        
        contractility(CellType::support,
                      CellType::progenitor) += incr_prog_supp;
        contractility(CellType::support, CellType::hair) = incr_hair_supp;
        contractility(CellType::support, CellType::support) += incr_supp_supp;

        // set the contractility and update the edge properties
        vertex_model.set_edge_contractility(contractility, true);
    };

    return std::make_pair(operation, params);
}

/// The operation to increment linetension
/** The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `progenitor_progenitor` (double, default: 0.): Incremental value for 
 *              edges between cells of types progenitor and progenitor.
 *      - `progenitor_hair` (double, default: 0.): Incremental value for 
 *              edges between cells of types progenitor and hair.
 *      - `progenitor_support` (double, default: 0.): Incremental value for 
 *              edges between cells of types progenitor and support.
 *      - `hair_hair` (double, default: 0.): Incremental value for 
 *              edges between cells of types hair and hair.
 *      - `hair_support` (double, default: 0.): Incremental value for 
 *              edges between cells of types hair and support.
 *      - `support_support` (double, default: 0.): Incremental value for 
 *              edges between cells of types support and support.
 */
OperationBundle build_increment_linetension (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);
    
    double incr_prog_prog(get_as<double>("progenitor_progenitor", cfg, 0.));
    double incr_prog_hair(get_as<double>("progenitor_hair", cfg, 0.));
    double incr_prog_supp(get_as<double>("progenitor_support", cfg, 0.));
    double incr_hair_hair(get_as<double>("hair_hair", cfg, 0.));
    double incr_hair_supp(get_as<double>("hair_support", cfg, 0.));
    double incr_supp_supp(get_as<double>("support_support", cfg, 0.));

    if (cfg["hair_progenitor"]) {
        double incr_hair_prog = get_as<double>("hair_progenitor", cfg);
        if (cfg["progenitor_hair"] and incr_hair_prog != incr_prog_hair) {
            throw KeyError("", cfg, fmt::format("Increment has to be "
                "symmetric matrix."));
        }
        incr_prog_hair = incr_hair_prog;
    }
    if (cfg["support_progenitor"]) {
        double incr_supp_prog = get_as<double>("hair_progenitor", cfg);
        if (cfg["progenitor_support"] and incr_supp_prog != incr_prog_supp) {
            throw KeyError("", cfg, fmt::format("Increment has to be "
                "symmetric matrix."));
        }
        incr_prog_supp = incr_supp_prog;
    }
    if (cfg["support_hair"]) {
        double incr_supp_hair = get_as<double>("hair_progenitor", cfg);
        if (cfg["hair_support"] and incr_supp_hair != incr_hair_supp) {
            throw KeyError("", cfg, fmt::format("Increment has to be "
                "symmetric matrix."));
        }
        incr_hair_supp = incr_supp_hair;
    }

    Operation operation = [incr_prog_prog, incr_prog_hair, incr_prog_supp,
                           incr_hair_hair, incr_hair_supp, incr_supp_supp]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        auto linetension = vertex_model.get_linetension();

        linetension(CellType::progenitor,
                    CellType::progenitor) += incr_prog_prog;
        linetension(CellType::progenitor, CellType::hair) += incr_prog_hair;
        linetension(CellType::progenitor, CellType::support) += incr_prog_supp;
        
        linetension(CellType::hair, CellType::progenitor) += incr_prog_hair;
        linetension(CellType::hair, CellType::hair) += incr_hair_hair;
        linetension(CellType::hair, CellType::support) += incr_hair_supp;
        
        linetension(CellType::support, CellType::progenitor) += incr_prog_supp;
        linetension(CellType::support, CellType::hair) = incr_hair_supp;
        linetension(CellType::support, CellType::support) += incr_supp_supp;

        // set the contractility and update the edge properties
        vertex_model.set_linetension(linetension, true);
    };

    return std::make_pair(operation, params);
}

/// The operation to increment the shape index
/** The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `progenitor` (double, default: 0.): Incremental value for cells of
 *              type progenitor.
 *      - `hair` (double, default: 0.): Incremental value for cells of
 *              type hair.
 *      - `support` (double, default: 0.): Incremental value for cells of
 *              type support.
 */
OperationBundle build_increment_shape_index (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);
    
    /// Incremental value for progenitor cells.
    /** \details default: 0 . pass as 'progenitor'. */
    double incr_prog(get_as<double>("progenitor", cfg, 0.));

    /// Incremental value for hair cells.
    /** \details default: 0 . pass as 'hair'. */
    double incr_hair(get_as<double>("hair", cfg, 0.));

    /// Incremental value for support cells.
    /** \details default: 0 . pass as 'support'. */
    double incr_supp(get_as<double>("support", cfg, 0.));

    Operation operation = [incr_prog, incr_hair, incr_supp]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        PCPVertex::RuleFuncCell update = [incr_prog, incr_hair, incr_supp]
                (const auto& cell)
        {
            auto state = cell->state;
            if (state.type == CellType::progenitor) {
                state.shape_index_preferential += incr_prog;
            }
            else if (state.type == CellType::support) {
                state.shape_index_preferential += incr_supp;
            }
            else if (state.type == CellType::hair) {
                state.shape_index_preferential += incr_hair;
            }
            return state;
        };

        apply_rule<Update::sync>(update, vertex_model.get_am().cells());
    };

    return std::make_pair(operation, params);
}

/// The operation to proliferate cells
/** Cell division happens as follows:
 *      #. Choose random cell of oldest generation.
 *      #. Increment area of cell to 2x preferential area. In 'num_increases'
 *         iterations. Minimize energy after every iteration.
 *      #. If area < threshold * 2x area_preferential: throw
 *      #. Divide cell in 2 daughter cells using the configuration of mother
 *         cell.
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `num_increases` (uint, > 0): In how many steps the area is to be increased
 *              to 2x A0 before division.
 *      - `area_threshold` (double, default: 0.): Threshold factor of how much 
 *              cell must increase: If area < threshold * 2x area_preferential
 *              throws.
 */
OperationBundle build_proliferate (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);
    
    auto num_increases(get_as<std::size_t>("num_increases", cfg));
    if (num_increases == 0) {
        throw std::invalid_argument("Parameter 'num_increases' in "
            "proliferation must be larger than 0. It specifies in how "
            "many steps the area of the cell should be increased.");
    }

    auto threshold(get_as<double>("area_threshold", cfg, 0.5));
    
    auto generation_max_id = std::make_shared<std::size_t>(0);
    
    std::uniform_real_distribution<double> prob_distr(0.,1.);

    Operation operation = [num_increases, threshold, generation_max_id,
                           prob_distr{std::move(prob_distr)}, params]
            (PCPVertex& vertex_model) mutable
    {
        using Cell = typename PCPVertex::Cell;

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        AgentContainer<Cell> current_generation(cells.size());
        auto it = std::copy_if (cells.begin(), cells.end(),
                                current_generation.begin(),
                                [generation_max_id](const auto& cell){
                                    return cell->id() <
                                        *generation_max_id; } );
        current_generation.resize(
            std::distance(current_generation.begin(), it));
        if (current_generation.size() == 0) {
            for (const auto& c : cells) {
                *generation_max_id = std::max(*generation_max_id, c->id());
            }
            *generation_max_id = *generation_max_id + 1;
            current_generation = cells;
        }

        std::uniform_int_distribution<> int_dist(
            0, current_generation.size() - 1);
        
        auto cell = current_generation[int_dist(*vertex_model.get_rng())];
        if (not vertex_model.get_space()->periodic) {
            for (auto [e, flip] : cell->custom_links().edges) {
                const auto [adj_cell_a, adj_cell_b] = am.adjoints_of(e);
                if (not adj_cell_a or not adj_cell_b) {
                    throw std::runtime_error(
                        "Cannot divide randomly chosen cell, because it is a "
                        "boundary cell. Division of boundary cells in "
                        "non-periodic boundary conditions is not implemented.");
                }
            }
        }

        double area_preferential = cell->state.area_preferential;
        double dA = cell->state.area_preferential / num_increases;
        for (unsigned int i = 0; i < num_increases; i++) {            
            vertex_model.increase_domain_size(dA);
            cell->state.area_preferential += dA;
            vertex_model.minimize_energy(params.minimization_params);
        }

        if (am.area_of(cell) < threshold * cell->state.area_preferential)
        {
            throw std::runtime_error(fmt::format("Cell division failed! "
                "Cell did not grow to area larger than threshold. "
                "For division requested minimal area: {}. \n"
                "For division preferred area: {}. \n"
                "Area reached: {}.", threshold * cell->state.area_preferential,
                cell->state.area_preferential, am.area_of(cell)));
        }

        cell->state.area_preferential = area_preferential;

        double angle = prob_distr(*vertex_model.get_rng()) * PI;
        vertex_model.divide_cell(cell, angle);        
    };

    return std::make_pair(operation, params);
}

/// The operation to do nothing but jiggle and minimize
OperationBundle build_jiggle (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    Operation operation = [] ([[maybe_unused]] PCPVertex& vertex_model) { };

    return std::make_pair(operation, params);
}

} // namespace OperationCollection
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif