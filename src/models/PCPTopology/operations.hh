#ifndef UTOPIA_MODELS_PCPTOPOLOGY_OPERATIONS_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_OPERATIONS_HH

#include <typeinfo>
#include <algorithm>
#include <iterator>

#include <utopia/core/types.hh>

#include "../PCPVertex/utils.hh"
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
    /** \note This are the times of the model as per result,
     *        hence at iteration 1 the time is resulting in 1.
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

    /// Additional minimisation estimate per iteration
    std::size_t add_num_minimisations;

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

    /// Whether to disable this operation
    /** Iterate the operation, but as void function.
     */
    bool disable;

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
        add_num_minimisations(0),
        probability(get_as<double>("probability", cfg, 1)),
        iterations_prolog(get_as<std::size_t>("iterations_prolog", cfg, 0)),
        iterations_epilog(get_as<std::size_t>("iterations_epilog", cfg, 0)),
        disable(get_as<bool>("disable", cfg, false)),
        minimization_mode(setup_minimization_mode(get_as<std::string>(
            "mode", get_as<Config>("minimization", cfg, Config()), "every"))),
        minimization_params(get_as<Config>("minimization", cfg, Config()),
                            default_minim_params)
    {
        if (not cfg["times"]) {
            throw std::runtime_error(fmt::format("Could not extract key "
                "'times' from operation '{}'!", name));
        }
        else if (cfg["times"].IsSequence()) {
            auto times_list = get_as<std::vector<Time>>("times", cfg);

            // Make sure negative times are not included
            times_list.erase(
                std::remove_if(times_list.begin(), times_list.end(),
                                [](auto& t){ return (t <= 0); }),
                times_list.end()
            );

            // Populate the set; this will impose ordering
            times.insert(times_list.begin(), times_list.end());
        }
        else if (cfg["times"].IsMap()) {
            auto start = get_as<Time>("start", cfg["times"]);
            auto stop = get_as<Time>("stop", cfg["times"]);
            auto step = get_as<Time>("step", cfg["times"]);

            for (Time t = start; t < stop; t += step) {
                if (t <= 0) {
                    continue;
                }
                times.insert(t);
            }
        }
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

    /// Get the number of minimizations triggered for this operations
    std::size_t get_num_minimizations(std::size_t time_max) const {
        std::size_t iters;
        if (minimization_mode == MinimizationMode::Manual) {
            iters = 0;
        }
        else if (minimization_mode == MinimizationMode::Once) {
            iters = std::accumulate(times.begin(), times.end(), 0,
                                    [time_max](auto val, auto t) {
                                        return val + (t <= time_max);
                                    });
            if (iterations_prolog > 0) {
                iters++;
            }
            if (iterations_epilog > 0) {
                iters++;
            }
        }
        else if (minimization_mode == MinimizationMode::Every) {
            iters = std::accumulate(times.begin(), times.end(), 0,
                                    [time_max](auto val, auto t) {
                                        return val + (t <= time_max);
                                    });
            iters *= iterations;
            iters += iterations_prolog + iterations_epilog;
        }
        else {
            throw std::invalid_argument(fmt::format("Cannot estimate "
                "iterations on unknown MinimizationMode {}",
                minimization_mode));
        }

        return iters * (  minimization_params.num_repeat
                        + add_num_minimisations);
    }
};

using Operation = std::function<void(PCPVertex& vertex_model)>;

using OperationBundle = typename std::pair<Operation, OperationParams>;


/// The operation to do nothing
OperationBundle build_void (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    Operation operation = [] ([[maybe_unused]] PCPVertex& vertex_model) { };

    return std::make_pair(operation, params);
}

/// Update the paramters of a registered work function
/** Takes configs:
 *      - `name' (str): The name of the WF-term to update
 *      - `parameters' (Config): Parameters forwarded to
 *              update_work_function_term
 */
OperationBundle build_update_work_function_term (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    auto wf_name = get_as<std::string>("name", cfg);
    auto update_params = get_as<Config>("parameters", cfg);

    Operation operation = [wf_name, update_params] (PCPVertex& vertex_model)
    {
        auto registered = vertex_model.is_registered_term(wf_name);
        if (registered) {
            vertex_model.update_work_function_term(wf_name, update_params);
        }
        else {
            vertex_model.get_logger()->warn("Cannot update unregistered "
                "work-function term {}!", wf_name);
        }
    };

    return std::make_pair(operation, params);
}

/// Update the paramters of a registered work function
/** Takes configs:
 *      - `term' (str): The reference to a work-function term
 *      - `name' (str, optional): The name given to this WF-term
 *      - `parameters' (Config): Parameters forwarded to the constructor
 *              of the selected term
 */
OperationBundle build_register_work_function_term (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    auto term = get_as<std::string>("term", cfg);
    auto wf_name = get_as<std::string>("name", cfg, term);
    auto wf_params = get_as<Config>("parameters", cfg);

    Operation operation = [term, wf_name, wf_params] (PCPVertex& vertex_model)
    {
        vertex_model.register_work_function_term(term, wf_name, wf_params);
    };

    return std::make_pair(operation, params);
}

/// Remove a registered work function
/** Takes configs:
 *      - `name' (str, optional): The name given to this WF-term
 */
OperationBundle build_unregister_work_function_term (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    auto wf_name = get_as<std::string>("name", cfg);

    Operation operation = [wf_name] (PCPVertex& vertex_model)
    {
        vertex_model.erase_work_function_term(wf_name);
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
 *      - `angle_distribution` (pair(double, double), optional): The mean and 
 *              stddev of a normal distribution for the angle of cell division 
 *              axis. In rad units [0, PI]. If not provided (or stddev < 1.e-11), a
 *              uniform distribution [0, Pi] is used, i.e. no preferential
 *              axis.
 */
template<typename Logger>
OperationBundle build_proliferate (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params,
        std::shared_ptr<Logger> logger,
        std::function<void()> monitor = [](){ })
{
    OperationParams params(name, cfg, default_minim_params);

    // Cells are picked from oldest generation at random
    // As the id of any new cell is incremented wrt youngest cell
    // this keeps track of the youngest cell of the oldest generation
    auto generation_max_id = std::make_shared<std::size_t>(0);

    // In how many steps to double area
    auto num_increases(get_as<std::size_t>("num_increases", cfg));
    MinimizationParams minimization_after_increase = params.minimization_params;
    if (num_increases > 0) {
        minimization_after_increase = MinimizationParams(
            get_as<Config>("minimization_after_increase", cfg),
            params.minimization_params);
    }
    
    // The minimum area for division (relative to 2 * A_0)
    auto threshold(get_as<double>("area_threshold", cfg));

    // A distribution of division angles in rad
    auto [mean, stddev] = get_as<std::pair<double, double>>(
        "angle_distribution", cfg, std::make_pair(0., 1.e-12));
    stddev = std::max(stddev, 1.e-12);
    std::normal_distribution<double> normal_distr(mean, stddev);
    std::uniform_real_distribution<double> uniform_distr(0., M_PI);

    Operation divide_cell = [num_increases, minimization_after_increase,
                             threshold, generation_max_id,
                             uniform_distr{std::move(uniform_distr)},
                             normal_distr{std::move(normal_distr)}, params]
            (PCPVertex& vertex_model) mutable
    {
        using Cell = typename PCPVertex::Cell;

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        // The cells of the oldest generation
        // i.e. all cells that have an id smaller than the max_id
        AgentContainer<Cell> current_generation(cells.size());
        auto it = std::copy_if (cells.begin(), cells.end(),
                                current_generation.begin(),
                                [generation_max_id, am](const auto& cell){
                                    return (cell->id() < *generation_max_id);
                                } );
        current_generation.resize(
            std::distance(current_generation.begin(), it));
        if (current_generation.size() == 0) {
            for (const auto& c : cells) {
                *generation_max_id = std::max(*generation_max_id, c->id());
            }
            *generation_max_id = *generation_max_id + 1;
            
            current_generation.clear();
            current_generation.resize(cells.size());
            std::copy(cells.begin(), cells.end(), current_generation.begin());
        }

        std::uniform_int_distribution<> int_dist(
            0, current_generation.size() - 1);
        
        auto cell = current_generation[int_dist(*vertex_model.get_rng())];

        if (num_increases > 0) {
            const double A0 = cell->state.area_preferential;
            double dA = A0 / num_increases;
            for (std::size_t i = 1; i <= num_increases; i++) {            
                vertex_model.increase_domain_size(dA);
                cell->state.area_preferential += dA;
                vertex_model.minimize_energy(minimization_after_increase);
            }

            if (am.area_of(cell) < threshold * A0)
            {
                throw std::runtime_error(fmt::format("Cell division failed! "
                    "Cell did not grow to area larger than threshold. "
                    "For division requested minimal area: {}. \n"
                    "For division preferred area: {}. \n"
                    "Area reached: {}.",
                    threshold * A0, A0, am.area_of(cell)
                ));
            }

            // reset parameters, but don't update
            cell->state.area_preferential = A0;
        }
    else {
        vertex_model.increase_domain_size(cell->state.area_preferential);
    }

        double angle;
        if (normal_distr.stddev() > 1.e-11) {
            angle = normal_distr(*vertex_model.get_rng());
        }
        else {
            angle = uniform_distr(*vertex_model.get_rng());
        }
        auto [ca, cb] = vertex_model.divide_cell(cell, angle);
    };
    std::size_t add_minimizations = (  minimization_after_increase.num_repeat
                                     * num_increases);


    Operation operation;
    auto num_divs = get_as<std::size_t>("num_cells", cfg, 1);
    if (num_divs > 1) {
        bool minimize = get_as<bool>("minimize_after_division", cfg);
        MinimizationParams minimization_after_division = params.minimization_params;
        if (minimize) {
            minimization_after_division = MinimizationParams(
                get_as<Config>("minimization_after_division", cfg),
                params.minimization_params);
            add_minimizations += minimization_after_division.num_repeat;
        }
        operation = [num_divs, minimization_after_division, minimize,
                     divide_cell, logger, monitor]
                    (PCPVertex& vertex_model)
        {
            logger->debug(" Proliferating {} cells ...", num_divs);
            for (std::size_t i = 0; i < num_divs; i++) {
                divide_cell(vertex_model);
                if (minimize) {
                    vertex_model.minimize_energy(minimization_after_division,
                                                 monitor);
                }
            }
        };
        add_minimizations *= num_divs;
    }
    else {
        operation = divide_cell;
    }

    params.add_num_minimisations = add_minimizations;
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
template<typename Logger>
OperationBundle build_proliferate_generations (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params,
        std::shared_ptr<Logger> logger = nullptr,
        std::function<void()> monitor = [](){ })
{
    OperationParams params(name, cfg, default_minim_params);

    auto num_generations = get_as<std::size_t>("num_generations", cfg);

    MinimizationParams minimization_between_generations(
        get_as<Config>("minimization_between_generations", cfg),
        params.minimization_params);

    // How to do the single cell division
    Config cfg_proliferation;
    cfg_proliferation["num_increases"] = 0;
    cfg_proliferation["area_threshold"] = 0.;

    Config minimization_tmp;
    minimization_tmp["mode"] = "manual";
    cfg_proliferation["minimization_after_increase"] = minimization_tmp;
    cfg_proliferation["times"] = std::vector<std::size_t>({});

    auto op_pair = build_proliferate(
        "proliferate", cfg_proliferation, minimization_between_generations,
        logger, monitor);
    auto proliferate = std::get<0>(op_pair);


    Operation operation = [proliferate, num_generations,
                           minimization_between_generations, logger, monitor]
            (PCPVertex& vertex_model)
    {
        if (not vertex_model.get_space()->periodic) {
            throw std::runtime_error("Proliferation of cells in generations "
                                     "is not suitable for non-periodic space!");
        }

        const auto& cells = vertex_model.get_am().cells();

        logger->info("Proliferating {} cells in {} generations. Expecting {} "
                     "cells after proliferation.", cells.size(),
                     num_generations,
                     cells.size() * std::pow(2, num_generations));
        if (cells.size() * std::pow(2, num_generations) > 512) {
            logger->warn("Expecting {} > 512 cells after proliferation "
                         "in {} generations. This might take a while ..",
                         cells.size() * std::pow(2, num_generations),
                         num_generations);

        }
        for (std::size_t gen = 0; gen < num_generations; gen++) {
            auto num_cells = cells.size();
            logger->debug(" Proliferating generation {} of {} cells, "
                          "then minimizing energy ...", gen, num_cells);
            for (std::size_t i = 0; i < num_cells; i++)
            {
                proliferate(vertex_model);
            }

            if (gen + 1 < num_generations) {
                logger->trace("  There are {} cells after proliferation. "
                            "Minimizing energy now", cells.size());
                vertex_model.minimize_energy(minimization_between_generations,
                                             monitor);
            }
        }

        logger->info("There are {} cells after proliferation.", cells.size());
    };

    params.add_num_minimisations = (
          minimization_between_generations.num_repeat
        * std::max(static_cast<int>(num_generations) - 1, 0));

    return std::make_pair(operation, params);
}


} // namespace OperationCollection
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif