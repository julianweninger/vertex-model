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
        auto times_list = get_as<std::vector<Time>>("times", cfg);
        // TODO Consider wrapping negative values around

        // Make sure negative times are not included
        times_list.erase(
            std::remove_if(times_list.begin(), times_list.end(),
                            [](auto& t){ return (t <= 0); }),
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

/// Minimize using alternating brownian motion and minimization
/** Iterates the vertex model using a brownian motion update scheme,
 *  then minimizes energy using the usual configuration.
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `num_repeat` (std::size_t): How often to repeat the cycle
 *      - `brownian_iteration` (minimization Config): The configuration for
 *              brownian motion minimization update.
 *      - `minimization` (minimization Config): The update to the default
 *              minimization config. Used to minimize after brownian motion.
 */
template<typename Logger>
OperationBundle build_brownian_noise (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params,
        std::shared_ptr<Logger> logger = nullptr,
        std::function<void()> monitor = [](){ })
{
    using MinimizationMode = OperationParams::MinimizationMode;

    OperationParams params(name, cfg, default_minim_params);

    MinimizationParams brownian_iteration(
        get_as<Config>("brownian_iteration", cfg),
        params.minimization_params);
    auto minimization = params.minimization_params;
    bool minimize = (    params.minimization_mode == MinimizationMode::Every
                     and not params.disable);

    auto num_repeat = get_as<std::size_t>("num_repeat", cfg);
    params.add_num_minimisations = 2 * num_repeat;

    Operation operation = [num_repeat, brownian_iteration, minimization,
                           minimize,
                           monitor, logger]
                          (PCPVertex& vertex_model)
    {
        logger->debug("Repeating cycle of minimization and brownian motion {} "
                      "times ..", num_repeat);
        for (std::size_t i = 0; i < num_repeat; i++) {
            logger->trace("   Minimizing energy {} / {} ..", i, num_repeat);
            if (minimize) {
                vertex_model.minimize_energy(minimization, monitor);
            }
            
            logger->trace("   Iterating brownian motion {} / {} ..",
                          i, num_repeat);
            vertex_model.minimize_energy(brownian_iteration, monitor);
        }
    };

    return std::make_pair(operation, params);
}


/// The operation to differentiate the types of cells randomly
/** Progenitor cells turn to hair cell with given probability and to support
 *  cell otherwise.
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `probability` (double): The probability for a cell to differentiate
 *              to a hair cell; otherwise support cell. Gives a fraction of 
 *              hair cells for large enough number of cells.
 *      - `fraction` (double): Similar to `probability`, but a fraction of the 
 *              cells is selected for HC differentiation. Resulting HC fraction
 *              is fixed. Number of cells is rounded to next smaller integer
 * 
 *  \note The entities properties do not change during differentiation.
 */
OperationBundle build_differentiate_domain (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);
    
    double x(get_as<double>("relative_width", cfg));
    double y(get_as<double>("relative_height", cfg, 1.));

    Operation operation = [x, y](PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;
        using SpaceVec = PCPVertex::SpaceVec;
        
        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        double abs_x, abs_y;
        if (vertex_model.get_space()->periodic) {
            SpaceVec domain = vertex_model.get_space()->get_domain_size();
            abs_x = x * domain[0];
            abs_y = y * domain[1];
        }
        else {
            double min_x = std::numeric_limits<double>::max();
            double max_x = std::numeric_limits<double>::min();
            double min_y = std::numeric_limits<double>::max();
            double max_y = std::numeric_limits<double>::min();
            for (const auto& vertex : am.vertices()) {
                SpaceVec pos = am.position_of(vertex);
                min_x = std::min(min_x, pos[0]);
                max_x = std::max(max_x, pos[0]);

                min_y = std::min(min_y, pos[1]);
                max_y = std::max(max_y, pos[1]);
            }
            abs_x = x * (max_x - min_x) - min_x;
            abs_y = y * (max_y - min_y) - max_x;
        }

        // Pick cells in lower left rectangle         
        apply_rule<Update::sync>(
            [am, abs_x, abs_y](const auto& cell)
            {
                auto state = cell->state;
                SpaceVec pos = am.barycenter_of(cell);

                if (pos[0] < abs_x and pos[1] < abs_y) {
                    state.type = CellType::hair;
                }
                else {
                    state.type = CellType::support;
                }

                return state;
            },
            cells
        );
    };

    return std::make_pair(operation, params);
}

/// The operation differentiate cells using the Collier model
/** \details The differentiation occurs as in the Collier model.
 *  The geometric properties and cell states are synchronized after n (`steps`)
 *  iterations of the Collier model.
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `steps` (uint, default: 1): The numer of steps performed in the
 *              Collier model without synchronisation to the Vertex model.
 * 
 *  \note The entities properties do not change during differentiation.
 * 
 *  \note The association of Vertex cells to Collier cells is arbitrary and
 *      might change with every call of operation. The cell's neighborhood is
 *      identical anyways.
 * 
 *  \param collier  Pointer to a Collier model
 *  \param prolog       Whether the prolog of Collier model already performed
 */
template <class Collier>
OperationBundle build_differentiate_Collier (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params,
        std::shared_ptr<Collier> collier,
        std::shared_ptr<bool> prolog)
{
    if (not collier) {
        throw std::runtime_error("Received nullptr in "
            "build_operation_Collier!");
    }
    
    OperationParams params(name, cfg, default_minim_params);

    std::size_t steps(get_as<std::size_t>("steps", cfg, 1));

    Operation operation = [collier, prolog, steps] (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;
        // Collier model cells
        using CCellType = typename Collier::CellType;

        collier->get_logger()->debug(
            "Setting up the custom neighbourhood of the cells as per the "
            "neighbourhood in the Vertex model ...");
            

        const auto& c_cells = collier->get_cm().cells();
        const auto& cells = vertex_model.get_am().cells();
        if (cells.size() > c_cells.size()) {
            throw std::runtime_error(fmt::format("Cannot link cells of "
                "Collier and Vertex models because more cells in Vertex "
                "({} cells) than in Collier model ({} cells)!",
                cells.size(), c_cells.size()));
        }

        unsigned int iterator;
        for (iterator = 0; iterator < cells.size(); iterator++) {
            cells[iterator]->custom_links().c_cell = c_cells[iterator];
            
            auto type = cells[iterator]->state.type;
            if (not *prolog) {
                void(); // use initialization of Collier Model
            }
            else if (type == CellType::hair) {
                c_cells[iterator]->state.cell_type = CCellType::hair;
            }
            else if (type == CellType::support) {
                c_cells[iterator]->state.cell_type = CCellType::support;
            }
            else if (type == CellType::progenitor) {
                c_cells[iterator]->state.cell_type = CCellType::progenitor;
            }
            else {
                throw std::runtime_error(fmt::format("Cell type {} from vertex "
                    "model not available in Collier model!", type));
            }
        }
        for (void(); iterator < c_cells.size(); iterator++) {
            c_cells[iterator]->state.cell_type = CCellType::inactive;
            c_cells[iterator]->custom_links().neighbors.clear();
        }

        for (const auto& cell : cells) {
            auto& c_cell = cell->custom_links().c_cell;
            c_cell->custom_links().neighbors.clear();
            for (auto n : vertex_model.get_am().neighbors_of(cell)) {
                c_cell->custom_links().neighbors.push_back(
                    n->custom_links().c_cell);
            }
        }
    
        if (not *prolog) {
            collier->prolog();
            *prolog = true;
        }


        collier->get_logger()->debug(
            "Itterating the Collier model for {} steps ...", steps);
        for (std::size_t i = 0; i < steps; i++) {
            collier->iterate();
        }
        collier->identify_clusters();

        collier->get_logger()->debug(
            "Synchronizing Collier changes to Vertex model ...");
        for (const auto& cell : vertex_model.get_am().cells()) {
            const auto& c_cell = cell->custom_links().c_cell;

            auto type = c_cell->state.cell_type;
            if (type == CCellType::hair) {
                cell->state.type = CellType::hair;
            }
            else if (type == CCellType::support) {
                cell->state.type = CellType::support;
            }
            else if (type == CCellType::progenitor){
                cell->state.type = CellType::progenitor;
            }
            else {
                throw std::runtime_error(fmt::format("Cell type {} from "
                    "Collier not available in vertex model!", type));
            }
        }
    };

    return std::make_pair(operation, params);
}

/// The operation to differentiate two hair cells in contact
/** The remaining cells become type progenitor
 * 
 *  \note mainly for testing purposes, i.e. which mechanism can separate two
 *        hair cells?
 *  
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `num_seeds` (std::size_t): The number of seeds for new clusters
 *      - `cluster_size` (std::size_t): How many next neighbours to add to the 
 *              cluster
 * 
 *  \note The entities properties do not change during differentiation.
 * 
 *  \warning not continuously tested
 */
OperationBundle build_differentiate_hair_cluster (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);
    
    auto num_seeds(get_as<std::size_t>("num_seeds", cfg));
    auto cluster_size(get_as<std::size_t>("cluster_size", cfg));
    
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


    Operation operation = [num_seeds, cluster_size](PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;
        
        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();
        AgentContainer<PCPVertex::Cell> cluster_seeds{};
        cluster_seeds.reserve(num_seeds);
        std::sample(cells.begin(), cells.end(),
                    std::back_inserter(cluster_seeds),
                    num_seeds, *vertex_model.get_rng());

        for (const auto& cell : cluster_seeds) {
            AgentContainer<PCPVertex::Cell> cluster({cell});
            std::size_t iterations = 0;
            while (    (cluster.size() < cluster_size)
                   and (iterations++ < 10 * cluster_size))
            {
                std::uniform_int_distribution<> int_distr(0, cluster.size() - 1);
                const auto& n = cluster[int_distr(*vertex_model.get_rng())];
                auto nbs = am.neighbors_of(n);
                int_distr = std::uniform_int_distribution<>(0, nbs.size() - 1);
                const auto& new_cell = nbs[int_distr(*vertex_model.get_rng())];

                if ((   std::find(cluster.begin(), cluster.end(), new_cell)
                     == cluster.end()))
                {
                    cluster.push_back(new_cell);
                }
            }

            PCPVertex::RuleFuncCell differentiate = [](const auto& cell)
            {
                cell->state.type = CellType::hair;
                return cell->state;
            };

            apply_rule<Update::sync>(differentiate, cluster);
        }

        PCPVertex::RuleFuncCell differentiate_others = [](const auto& cell)
        {
            auto state = cell->state;
            if (state.type == CellType::progenitor) {
                state.type = CellType::support;
            }
            return state;
        };

        apply_rule<Update::sync>(differentiate_others, cells);
    };

    return std::make_pair(operation, params);
}


/// The operation differentiate cells using the NotchDelta model
/** \details The differentiation occurs as in the NotchDelta model.
 *  The geometric properties and cell states are synchronized after n (`steps`)
 *  iterations of the NotchDelta model.
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `steps` (uint, default: 1): The numer of steps performed in the
 *              NotchDelta model without synchronisation to the Vertex model.
 *      - `discard_initialisation` (bool, default: false): Discard the initial
 *              state of the NotchDelta model. Overwrite it with the state of 
 *              vertex-model cells.
 * 
 *  \note The entities properties do not change during differentiation.
 * 
 *  \note The association of Vertex cells to NotchDelta cells is arbitrary and
 *      might change with every call of operation. The cell's neighborhood is
 *      identical anyways.
 * 
 *  \param notch_delta  Pointer to a notch delta model
 *  \param prolog       Whether the prolog of ND is already finished
 */
template <class NotchDelta>
OperationBundle build_differentiate_NotchDelta (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params,
        std::shared_ptr<NotchDelta> notch_delta,
        std::shared_ptr<bool> prolog)
{
    if (not notch_delta) {
        throw std::runtime_error("Received nullptr in "
            "build_operation_NotchDelta!");
    }
    
    OperationParams params(name, cfg, default_minim_params);

    auto steps(get_as<std::size_t>("steps", cfg, 1));
    auto discard_initialisation(get_as<bool>("discard_initialisation", cfg,
                                             false));

    Operation operation = [notch_delta, prolog, steps,
                           discard_initialisation] (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;
        using NDCellType = typename NotchDelta::CellType;

        notch_delta->get_logger()->debug(
            "Setting up the custom neighbourhood of the cells as per the "
            "neighbourhood in the Vertex model ...");            

        const auto& nd_cells = notch_delta->get_cm().cells();
        const auto& cells = vertex_model.get_am().cells();
        if (cells.size() > nd_cells.size()) {
            throw std::runtime_error(fmt::format("Cannot link cells of "
                "NotchDelta and Vertex models because more cells in Vertex "
                "({} cells) than in NotchDelta model ({} cells)!",
                cells.size(), nd_cells.size()));
        }

        unsigned int iterator;
        for (iterator = 0; iterator < cells.size(); iterator++) {
            cells[iterator]->custom_links().nd_cell = nd_cells[iterator];
            
            auto type = cells[iterator]->state.type;
            if (not *prolog and not discard_initialisation) {
                void(); // use initialization of ND Model
            }
            else if (type == CellType::hair) {
                nd_cells[iterator]->state.cell_type = NDCellType::hair;
            }
            else if (type == CellType::support) {
                nd_cells[iterator]->state.cell_type = NDCellType::support;
            }
            else if (type == CellType::progenitor) {
                nd_cells[iterator]->state.cell_type = NDCellType::progenitor;
            }
            else {
                throw std::runtime_error(fmt::format("Cell type {} from vertex "
                    "model not available in NotchDelta model!", type));
            }
        }
        for (void(); iterator < nd_cells.size(); iterator++) {
            nd_cells[iterator]->state.cell_type = NDCellType::inactive;
            nd_cells[iterator]->custom_links().neighbors.clear();
        }

        for (const auto& cell : cells) {
            auto& nd_cell = cell->custom_links().nd_cell;
            nd_cell->custom_links().neighbors.clear();
            for (auto n : vertex_model.get_am().neighbors_of(cell)) {
                nd_cell->custom_links().neighbors.push_back(
                    n->custom_links().nd_cell);
            }
        }
    
        if (not *prolog) {
            notch_delta->prolog();
            *prolog = true;
        }

        notch_delta->get_logger()->debug(
            "Itterating the NotchDelta model for {} steps ...", steps);   
        for (std::size_t i = 0; i < steps; i++) {
            notch_delta->iterate();
        }
        notch_delta->identify_clusters();

        notch_delta->get_logger()->debug(
            "Synchronizing NotchDelta changes to Vertex model ...");   
        for (const auto& cell : vertex_model.get_am().cells()) {
            const auto& nd_cell = cell->custom_links().nd_cell;

            auto type = nd_cell->state.cell_type;
            if (type == NDCellType::hair) {
                cell->state.type = CellType::hair;
            }
            else if (type == NDCellType::support) {
                cell->state.type = CellType::support;
            }
            else if (type == NDCellType::progenitor) {
                cell->state.type = CellType::progenitor;
            }
            else {
                throw std::runtime_error(fmt::format("Cell type {} from "
                    "NotchDelta model not available in Vertex model!", type));
            }
        }
    };

    return std::make_pair(operation, params);
}


/// The operation to differentiate the types of cells randomly
/** Progenitor cells turn to hair cell with given probability and to support
 *  cell otherwise.
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `probability` (double): The probability for a cell to differentiate
 *              to a hair cell; otherwise support cell. Gives a fraction of 
 *              hair cells for large enough number of cells.
 *      - `fraction` (double): Similar to `probability`, but a fraction of the 
 *              cells is selected for HC differentiation. Resulting HC fraction
 *              is fixed. Number of cells is rounded to next smaller integer.
 * 
 *  \note The entities properties do not change during differentiation.
 */
OperationBundle build_differentiate_random (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);
    
    double probability(get_as<double>("probability", cfg, 0.));
    double fraction(get_as<double>("fraction", cfg, 0.));

    if (probability > 1.e-12 and fraction > 1.e-12) {
        throw std::invalid_argument(fmt::format(
            "In operation `differentiate random: "
            "Received probability and fraction > 0, choose either! "
            "Probability was {} and fraction was {}",
            probability, fraction));
    }
    if (fraction > 1 or fraction < 0) {
        throw std::invalid_argument(fmt::format(
            "In operation `differentiate random: `fraction` needs to be in "
            "[0, 1], but was {}!",
            fraction));
    }
    
    if (params.iterations_prolog + params.iterations_epilog +
        params.iterations * (params.times.size()) > 1)
    {
        throw std::invalid_argument(fmt::format(
            "In operation `differentiate random: "
            "Differentiate can be applied only once, because terminal "
            "process. Iterations in prolog: {}, in run {}, in epilog",
            params.iterations_prolog,
            params.iterations * (params.times.size()),
            params.iterations_epilog));
    }

    std::uniform_real_distribution<double> prob_distr(0., 1.);

    if (fraction > 1.e-12) {
        Operation operation = [fraction, prob_distr{std::move(prob_distr)}]
                (PCPVertex& vertex_model) mutable
        {
            using CellType = PCPVertex::CellType;

            // Differentiate cell as HCs
            PCPVertex::RuleFuncCell update_HCs = [
                            vertex_model,
                            prob_distr{std::move(prob_distr)}]
                    (const auto& cell) mutable
            {
                cell->state.type = CellType::hair;

                return cell->state;
            };

            // Pick n = f*N random cells 
            auto cells = vertex_model.get_am().cells();
            std::shuffle(cells.begin(), cells.end(), *vertex_model.get_rng());
            cells.erase(cells.begin() + std::size_t(fraction*cells.size()),
                        cells.end());

            // Differentiate n cells as HCs
            apply_rule<Update::sync>(update_HCs, cells);

            // Differentiate others as SCs
            PCPVertex::RuleFuncCell update_SCs = [
                            vertex_model,
                            prob_distr{std::move(prob_distr)}]
                    (const auto& cell) mutable
            {
                auto state = cell->state;
                if (state.type == CellType::progenitor) {
                    state.type = CellType::support;
                }

                return state;
            };
            apply_rule<Update::sync>(update_SCs, vertex_model.get_am().cells());

        };

        return std::make_pair(operation, params);
    }
    else {
        Operation operation = [probability, prob_distr{std::move(prob_distr)}]
                (PCPVertex& vertex_model) mutable
        {
            using CellType = PCPVertex::CellType;

            // differentiate cell HC with probability p
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
}

/// Enable or disable topological transitions
/** The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `enable_T1s` (bool): Whether to enable T1 transitions
 *      - `enable_T2s` (bool): Whether to enable T2 transitions
 */
OperationBundle build_enable_transitions (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    auto enable_T1s = get_as<bool>("enable_T1s", cfg);
    auto enable_T2s = get_as<bool>("enable_T2s", cfg);

    Operation operation = [enable_T1s, enable_T2s](PCPVertex& vertex_model)
    {
        vertex_model.enable_transitions(enable_T1s, enable_T2s);
    };

    return std::make_pair(operation, params);
}

/// An operation to fix boundary vertices in space
/** Sets the fix_boundary in the vertex model and is updated continuously
 *  in update
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `fix_boundary` (bool, default: true): Whether to fix the boundary
 *                                              in space
 */
OperationBundle build_fix_boundary (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    bool fix_boundary = get_as<bool>("fix_boundary", cfg, true);

    Operation operation = [fix_boundary] (PCPVertex& vertex_model)
    {
        vertex_model.fix_boundary(fix_boundary);
    };

    return std::make_pair(operation, params);
}

/// The operation to increment preferential area
/** \details The following parameter are extracted from cfg 
 *           (besides those passed to `OperationParams`):
 *               - `progenitor` (double, default: 0.): incremental value for
 *                      cells of type progenitor
 *               - `hair` (double, default: 0.): incremental value for cells
 *                      of type hair
 *               - `support` (double, default: 0.): incremental value for cells
 *                      of type support
 *               - `adapt_support` (bool, default: false): If true, the total
 *                      area of hair and support cells remains constant, hence
 *                      support cells compensate area changes from hair cells.
 *               - `relax_domain` (bool, default: false): If true, the domain
 *                      size is adapted to fit the cells as per preferential 
 *                      area.
 *               - `gradient_hair_max` and `gradient_hair_min`
 *                  (double, default: 0.): The difference between
 *                      the cells that grow the most and the least. Using a
 *                      negative slope for non-periodic and a cos(x)**2 in 
 *                      periodic BC along the horizontal (P-D) axis.
 *                      If `adapt_support`, the total A0 of all cells remains
 *                      constant (uniform substraction of the change)
 *               - `hair_gradient_center` (double, default: 0.5): The relative 
 *                      position of dA = 0 in non-periodic BC and the relative 
 *                      position of dA = max in periodic BC.
 * 
 *  \note This affects the current entities properties, it does not overwrite
 *        changes in the past.
 */
OperationBundle build_increment_area (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    double prog(get_as<double>("progenitor", cfg, 0.));
    double hair(get_as<double>("hair", cfg, 0.));
    double support(get_as<double>("support", cfg, 0.));
    bool adapt_support(get_as<bool>("adapt_support", cfg, false));
    bool relax_domain(get_as<bool>("relax_domain", cfg, false));
    bool relax_domain_PD_axis(get_as<bool>("relax_domain_PD_axis",
                                           cfg, false));
    double gradient_hair_max(get_as<double>("gradient_hair_max", cfg, 0.));
    double gradient_hair_min(get_as<double>("gradient_hair_min", cfg, 0.));
    double hair_gradient_center(get_as<double>("hair_gradient_center",
                                               cfg, 0.5));
    
    if (adapt_support and fabs(support) > 1.e-12) {
        throw std::invalid_argument(fmt::format(
            "In cfg {}: If `adapt_support = True`, "
            "then the increment value for the support cells "
            "`support` must be zero, but was {}!", name, support));
    }

    if (adapt_support + relax_domain + relax_domain_PD_axis > 1) {
        throw std::invalid_argument(fmt::format(
            "In cfg {}: `adapt_support` and `relax_domain` cannot both be "
            "true!", name));
    }

    Operation operation = [prog, hair, support, adapt_support,
                           gradient_hair_max, gradient_hair_min,
                           hair_gradient_center, 
                           relax_domain, relax_domain_PD_axis]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;
        using SpaceVec = PCPVertex::SpaceVec;

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        double _support = support;

        if (adapt_support) {
            std::size_t num_hcs = std::count_if(
                cells.begin(), cells.end(), 
                [](const auto& c) {
                    return c->state.type == CellType::hair;
                });

            if (num_hcs < cells.size()) {
                _support = -1. * hair * num_hcs / (cells.size() - num_hcs);
            }
        }

        PCPVertex::RuleFuncCell update = [vertex_model, prog, hair, _support]
                (const auto& cell)
        {
            auto state = cell->state;
            
            double dA0 = 0.;
            if (state.type == CellType::progenitor) {
                dA0 = prog;
            }
            else if (state.type == CellType::hair) {
                dA0 = hair;
            }
            else if (state.type == CellType::support) {
                dA0 = _support;
            }

            state._area_preferential += dA0;
            return state;
        };
        
        apply_rule<Update::sync>(update, cells);

        if (fabs(gradient_hair_max) > 1.e-12) {
            auto [min, max] = am.get_extent();

            double L = (max - min)[0];
            double x_min = min[0];

            std::function<double(const SpaceVec&)> map;
            if (not vertex_model.get_space()->periodic) {
                map = 
                    [gradient_hair_max, gradient_hair_min,
                     hair_gradient_center, L, x_min]
                    (const SpaceVec& pos)
                    {
                        double hair_gradient = (  gradient_hair_max
                                                - gradient_hair_min);
                        double rel_x = (pos[0] - x_min) / L;
                        double fac = (hair_gradient_center - rel_x);
                        return hair_gradient * fac + gradient_hair_min;
                    };
            }
            else {
                map = 
                    [gradient_hair_max, gradient_hair_min,
                     hair_gradient_center, L, x_min]
                    (const SpaceVec& pos)
                    {
                        double hair_gradient = (  gradient_hair_max
                                                - gradient_hair_min);
                        double rel_x = (  (pos[0] - x_min) / L
                                        - hair_gradient_center);
                        double fac = std::pow(cos(rel_x * M_PI), 2);
                        return hair_gradient * fac + gradient_hair_min;
                    };
            }

            double area = std::accumulate(
                cells.begin(), cells.end(), 0.,
                [](const double& val, const auto& cell){
                    return val + cell->state._area_preferential;
                });

            for (const auto& cell : cells) {
                if (cell->state.type == CellType::hair) {
                    SpaceVec pos = am.barycenter_of(cell);
                    double dA = map(pos);
                    cell->state._area_preferential += dA;
                }
            }
            if (adapt_support) {
                double new_area = std::accumulate(
                    cells.begin(), cells.end(), 0.,
                    [](const double& val, const auto& cell){
                        return val + cell->state._area_preferential;
                    });

                for (const auto& cell : cells) {
                    cell->state._area_preferential -= (  (new_area - area)
                                                       / cells.size());
                }
            }
        }

        if (relax_domain) {
            double area = std::accumulate(cells.begin(), cells.end(), 0.,
                                [](const double& val, const auto& cell) {
                                    return val + cell->state._area_preferential;
                                });

            PCPVertex::SpaceVec domain_size = 
                    vertex_model.get_space()->get_domain_size();
            vertex_model.increase_domain_size(
                area - domain_size[0]*domain_size[1]);
        }
        else if (relax_domain_PD_axis) {
            double area = std::accumulate(cells.begin(), cells.end(), 0.,
                            [](const double& val, const auto& cell) {
                                return val + cell->state._area_preferential; });

            PCPVertex::SpaceVec domain_size = 
                    vertex_model.get_space()->get_domain_size();
            double dA = area - domain_size[0]*domain_size[1];

            vertex_model.stretch_domain({dA/domain_size[1], 0.},
                                        false, false, false);
        }

        return;
    };

    return std::make_pair(operation, params);
}

/// The operation to increment cell contractility
/** \details The following parameter are extracted from cfg 
 *           (besides those passed to `OperationParams`):
 *               - `progenitor` (double, default: 0.): incremental value for
 *                      cells of type progenitor
 *               - `hair` (double, default: 0.): incremental value for cells
 *                      of type hair
 *               - `support` (double, default: 0.): incremental value for cells
 *                      of type support
 * 
 *  \note This affects the current entities properties, it does not overwrite
 *        changes in the past.
 */
OperationBundle build_increment_cell_contractility (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    double incr_prog(get_as<double>("progenitor", cfg, 0.));
    double incr_hair(get_as<double>("hair", cfg, 0.));
    double incr_support(get_as<double>("support", cfg, 0.));

    Operation operation = [incr_prog, incr_hair, incr_support]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();
        
        PCPVertex::RuleFuncCell update = [incr_prog, incr_hair, incr_support]
                (const auto& cell)
        {
            auto state = cell->state;

            if (state.type == CellType::progenitor) {
                state.contractility += incr_prog;
            }
            else if (state.type == CellType::hair) {
                state.contractility += incr_hair;
            }
            else if (state.type == CellType::support) {
                state.contractility += incr_support;
            }

            return state;
        };
        
        apply_rule<Update::sync>(update, cells);

        return;
    };

    return std::make_pair(operation, params);
}

/// The operation to increment the domain size
/** The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `value` (SpaceVec<2> or float): The incremental increase in x and y,
 *              or the incremental area (ratio Lx to Ly constant)
 *      - `compensate` (bool): Whether to compensate the increase in area in the
 *              preferential area of the cells. If true, keeps the mechanical
 *              properties constant
 *      - `fix_hair_cell_area` (bool, default: false): Whether to fix the area
 *              of cells of type hair.
 */
OperationBundle build_increment_domain (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    using SpaceVec = PCPVertex::SpaceVec;

    OperationParams params(name, cfg, default_minim_params);

    if (cfg["value"] and cfg["value"].IsSequence()) {
        SpaceVec increment(get_as_SpaceVec<2>("value", cfg));
        bool compensate(get_as<bool>("compensate", cfg));
        bool fix_hc_area(get_as<bool>("fix_hair_cell_area", cfg, false));
        bool fix_sc_area(get_as<bool>("fix_support_cell_area", cfg, false));

        if (compensate and fix_hc_area and fix_sc_area) {
            throw std::invalid_argument("Cannot compensate area while fixing "
                "hair and support cell area in operation 'increment_domain'!");
        }
        
        Operation operation = [increment, compensate, fix_hc_area,
                               fix_sc_area] (PCPVertex& vertex_model)
        {
            vertex_model.stretch_domain(increment, compensate,
                                        fix_hc_area, fix_sc_area);
        };

        return std::make_pair(operation, params);
    }
    else {
        double increment(get_as<double>("value", cfg));

        Operation operation = [increment] (PCPVertex& vertex_model)
        {
            vertex_model.increase_domain_size(increment);
        };

        return std::make_pair(operation, params);
    }
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
 * 
 *  The following paramter are excluded and passed on to the Vertex model
 *      - `set_ppMLC` (double, default: 0.): The contractility of SC-SC
 *              junctions with a specified orientation.
 *              See PCPVertex::set_ppMLC_contractility.
 *      - `set_ppMLC_angle` (double, default: 0.): The angle of the ppMLC
 *              contractility orientation, where 0 is the (1, 0) axis and
 *              Pi the (-1, 0) axis.
 *              See PCPVertex::set_ppMLC_contractility.
 *      - `set_pMLC` (double, default: 0.): The contractility of HC-SC junctions
 *              where the orientation is specified by the HC's polarity vector.
 *              See PCPVertex::set_pMLC_contractility.
 *              
 * 
 *  \note This increments the matrix of contractility, hence affects current and
 *        future edges between cells of corresponding type.
 */
OperationBundle build_increment_edge_contractility (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    using SpaceVec = PCPVertex::SpaceVec;

    OperationParams params(name, cfg, default_minim_params);
    
    double incr_prog_prog(get_as<double>("progenitor_progenitor", cfg, 0.));
    double incr_prog_hair(get_as<double>("progenitor_hair", cfg, 0.));
    double incr_prog_supp(get_as<double>("progenitor_support", cfg, 0.));
    double incr_prog_bnd(get_as<double>("progenitor_boundary", cfg, 0.));
    double incr_hair_hair(get_as<double>("hair_hair", cfg, 0.));
    double incr_hair_supp(get_as<double>("hair_support", cfg, 0.));
    double incr_hair_bnd(get_as<double>("hair_boundary", cfg, 0.));
    double incr_supp_supp(get_as<double>("support_support", cfg, 0.));
    double incr_supp_bnd(get_as<double>("support_boundary", cfg, 0.));

    double ppMLC(get_as<double>("set_ppMLC", cfg, 0.));
    double ppMLC_angle(get_as<double>("set_ppMLC_angle", cfg, 0.));
    auto ppMLC_axis = SpaceVec({cos(ppMLC_angle), sin(ppMLC_angle)});

    double pMLC(get_as<double>("set_pMLC", cfg, 0.));
    auto reset_HC_polarity(get_as<std::pair<bool, double>>(
        "reset_HC_polarity", cfg, std::make_pair(false, -M_PI_2)));
    std::size_t reset_HC_polarity_antispin_nshift = get_as<std::size_t>(
        "reset_HC_polarity_antispin_nshift", cfg, 0);
    std::size_t reset_HC_polarity_antispin_n = get_as<std::size_t>(
        "reset_HC_polarity_antispin_n", cfg, 0);
    auto reset_HC_polarity_shuffle(get_as<double>(
        "reset_HC_polarity_shuffle", cfg, 0.));

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
                           incr_hair_hair, incr_hair_supp, incr_supp_supp,
                           incr_prog_bnd, incr_hair_bnd, incr_supp_bnd,
                           ppMLC, ppMLC_axis, pMLC,
                           reset_HC_polarity, 
                           reset_HC_polarity_antispin_nshift,
                           reset_HC_polarity_antispin_n,
                           reset_HC_polarity_shuffle]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        auto contractility = vertex_model.get_edge_contractility();

        contractility(CellType::progenitor,
                      CellType::progenitor) += incr_prog_prog;
        contractility(CellType::progenitor, CellType::hair) += incr_prog_hair;
        contractility(CellType::progenitor,
                      CellType::support) += incr_prog_supp;
        contractility(CellType::progenitor,
                      CellType::num_cell_types) += incr_prog_bnd;
        
        contractility(CellType::hair, CellType::progenitor) += incr_prog_hair;
        contractility(CellType::hair, CellType::hair) += incr_hair_hair;
        contractility(CellType::hair, CellType::support) += incr_hair_supp;
        contractility(CellType::hair,
                      CellType::num_cell_types) += incr_hair_bnd;
        
        contractility(CellType::support,
                      CellType::progenitor) += incr_prog_supp;
        contractility(CellType::support, CellType::hair) += incr_hair_supp;
        contractility(CellType::support, CellType::support) += incr_supp_supp;
        contractility(CellType::support,
                      CellType::num_cell_types) += incr_supp_bnd;
                      
        contractility(CellType::num_cell_types,
                      CellType::progenitor) += incr_prog_bnd;
        contractility(CellType::num_cell_types,
                      CellType::hair) += incr_hair_bnd;
        contractility(CellType::num_cell_types,
                      CellType::support) += incr_supp_bnd;

        // set the contractility and update the edge properties
        vertex_model.set_edge_contractility(contractility, true);

        vertex_model.set_ppMLC_contractility(ppMLC, ppMLC_axis);
        vertex_model.set_pMLC_contractility(pMLC);

        if (std::get<bool>(reset_HC_polarity)) {
            const auto& cells = vertex_model.get_am().cells();

            std::size_t nshift = reset_HC_polarity_antispin_nshift;
            std::size_t n = reset_HC_polarity_antispin_n;

            std::size_t cell_cnt = 0;
            for (const auto& cell : cells) {
                if (cell->state.type == CellType::support) {
                    continue;
                }
                auto state = cell->state;
                state._polarity = std::get<double>(reset_HC_polarity);
                state._polarity_fluctuations = 0.;

                std::size_t flip = 0;
                // in a flip row
                if (nshift > 0 and (cell_cnt / nshift) % 2 == 1) {
                    flip++;
                }
                // in a row flip each nth
                if (n > 0 and cell_cnt % n == n - 1) {
                    flip++;
                }
                if (flip % 2 == 1) {
                    state._polarity = std::fmod(state._polarity + M_PI, 
                                                2 * M_PI);
                }
                cell->state = state;

                cell_cnt += 1;
            }
        }

        if (fabs(reset_HC_polarity_shuffle) > 1.e-12) {
            auto distr = std::uniform_real_distribution<double>(
                -reset_HC_polarity_shuffle, reset_HC_polarity_shuffle);
            const auto& cells = vertex_model.get_am().cells();
            const auto& rng = vertex_model.get_rng();
            for (const auto& cell : cells) {
                cell->state._polarity = std::fmod(
                    cell->state._polarity + distr(*rng) + 2 * M_PI, 
                    2 * M_PI);
                cell->state._polarity_fluctuations = 0.;
            }
        }


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
 * 
 *  \note This increments the matrix of linetension, hence affects current and
 *        future edges between cells of corresponding type.
 */
OperationBundle build_increment_linetension (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);
    
    double incr_prog_prog(get_as<double>("progenitor_progenitor", cfg, 0.));
    double incr_prog_hair(get_as<double>("progenitor_hair", cfg, 0.));
    double incr_prog_supp(get_as<double>("progenitor_support", cfg, 0.));
    double incr_prog_bnd(get_as<double>("progenitor_boundary", cfg, 0.));
    double incr_hair_hair(get_as<double>("hair_hair", cfg, 0.));
    double incr_hair_supp(get_as<double>("hair_support", cfg, 0.));
    double incr_hair_bnd(get_as<double>("hair_boundary", cfg, 0.));
    double incr_supp_supp(get_as<double>("support_support", cfg, 0.));
    double incr_supp_bnd(get_as<double>("support_boundary", cfg, 0.));

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
                           incr_hair_hair, incr_hair_supp, incr_supp_supp,
                           incr_prog_bnd, incr_hair_bnd, incr_supp_bnd]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        auto linetension = vertex_model.get_linetension();

        linetension(CellType::progenitor,
                    CellType::progenitor) += incr_prog_prog;
        linetension(CellType::progenitor, CellType::hair) += incr_prog_hair;
        linetension(CellType::progenitor, CellType::support) += incr_prog_supp;
        linetension(CellType::progenitor,
                    CellType::num_cell_types) += incr_prog_bnd;
        
        linetension(CellType::hair, CellType::progenitor) += incr_prog_hair;
        linetension(CellType::hair, CellType::hair) += incr_hair_hair;
        linetension(CellType::hair, CellType::support) += incr_hair_supp;
        linetension(CellType::hair,
                    CellType::num_cell_types) += incr_hair_bnd;
        
        linetension(CellType::support, CellType::progenitor) += incr_prog_supp;
        linetension(CellType::support, CellType::hair) += incr_hair_supp;
        linetension(CellType::support, CellType::support) += incr_supp_supp;
        linetension(CellType::support,
                    CellType::num_cell_types) += incr_supp_bnd;
                    
        linetension(CellType::num_cell_types,
                    CellType::progenitor) += incr_prog_bnd;
        linetension(CellType::num_cell_types,
                    CellType::hair) += incr_hair_bnd;
        linetension(CellType::num_cell_types,
                    CellType::support) += incr_supp_bnd;

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
 * 
 *  \note This affects the current entities properties, it does not overwrite
 *        changes in the past.
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
            else if (state.type == CellType::hair) {
                state.shape_index_preferential += incr_hair;
            }
            else if (state.type == CellType::support) {
                state.shape_index_preferential += incr_supp;
            }
            return state;
        };

        apply_rule<Update::sync>(update, vertex_model.get_am().cells());
    };

    return std::make_pair(operation, params);
}

/// Initialises a stripe that fits all vertices
OperationBundle build_initialise_stripe_boundary(
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    Operation operation = [cfg](PCPVertex& vertex_model) {
        vertex_model.initialise_stripe_boundary(cfg);
    };

    return std::make_pair(operation, params);
}

/// Increment the width of the stripe
/** Stripe needs to be initialised!
 */
template<typename Logger>
OperationBundle build_increment_stripe_width (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params,
        std::shared_ptr<Logger> logger = nullptr,
        std::function<void()> monitor = [](){ })
{
    OperationParams params(name, cfg, default_minim_params);

    // the incremental change in height
    double dH = get_as<double>("dH", cfg);

    // the target height (if > 0)
    auto H_target = get_as<double>("H_target", cfg, 0.);
    bool relative_target = get_as<bool>("relative_target", cfg, false);
    std::shared_ptr<double> H_final = std::make_shared<double>(0.);

    std::shared_ptr<double> update_potential(nullptr);
    if (cfg["update_potential_constant"]) {
        update_potential = std::make_shared<double>(
            get_as<double>("update_potential_constant", cfg));
    }

    MinimizationParams minimization(
        get_as<Config>("minimization", cfg, Config()),
                       params.minimization_params);

    Operation operation =
    [dH, H_target, relative_target, H_final, update_potential,
     minimization, logger, monitor]
    (PCPVertex& vertex_model)
    {
        if (update_potential) {
            vertex_model.update_stripe_boundary_potential(*update_potential);
        }

        double H = vertex_model.get_stripe_boundary_width();
        if (H > 1.e12) {
            throw std::runtime_error("Stripe boundary needs to be initialised "
                "before it can be changed in `increment_stipe_width`!");
        }
        if (*H_final < 1.e-12) {
            *H_final = H_target;
            if (relative_target) {
                *H_final *= H;
            }
        }

        // Reduce H by dH until H_final
        if (*H_final > 1.e-12) {
            while (H > *H_final + 1.e-12) {
                H = vertex_model.increment_stripe_boundary_width(
                    std::max(dH, *H_final - H));
                if (H > *H_final + 1.e-12) {
                    vertex_model.minimize_energy(minimization, monitor);
                }
            }
        }

        // reduce by dH
        else {
            vertex_model.increment_stripe_boundary_width(dH);
        }
    };

    return std::make_pair(operation, params);
}

/// Increment the width of the stripe
/** Stripe needs to be initialised!
 */
OperationBundle build_increment_stripe_curvature (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    double dk = get_as<double>("increment_curvature", cfg);

    std::shared_ptr<double> update_potential(nullptr);
    if (cfg["update_potential_constant"]) {
        update_potential = std::make_shared<double>(
            get_as<double>("update_potential_constant", cfg));
    }

    Operation operation = [dk, update_potential](PCPVertex& vertex_model)
    {
        double H = vertex_model.get_stripe_boundary_width();
        if (H > 1.e12) {
            throw std::runtime_error("Stripe boundary needs to be initialised "
                "before it can be changed in `increment_stipe_curvature`!");
        }

        if (update_potential) {
            vertex_model.update_stripe_boundary_potential(*update_potential);
        }

        vertex_model.increment_stripe_boundary_curvature(dk);
    };

    return std::make_pair(operation, params);
}

template <class PlanarCellPolarity, typename Logger>
OperationBundle build_iterate_pcp (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params,
        std::shared_ptr<PlanarCellPolarity> pcp,
        std::shared_ptr<bool> prolog,
        std::shared_ptr<Logger> logger = nullptr,
        std::function<void()> monitor = [](){ })
{
    if (not pcp) {
        throw std::runtime_error("Received nullptr in build_iterate_pcp!");
    }
    
    OperationParams params(name, cfg, default_minim_params);

    auto num_steps(get_as<std::size_t>("num_steps", cfg, 1));
    auto pcp_cfg(get_as<Config>("PlanarCellPolarity", cfg, Config()));

    Operation operation = [pcp, prolog, pcp_cfg, num_steps, monitor, logger]
            ([[maybe_unused]] PCPVertex& vertex_model)
    {
        pcp->update_parameters(pcp_cfg);

        if (not *prolog) {
            pcp->prolog();
            *prolog = true;
        }

        for (std::size_t step = 0; step < num_steps; step++) {
            pcp->iterate();
            monitor();

            if (stop_now.load()) {
                pcp->get_logger()->warn("Was told to stop. Not iterating "
                    "further ...");
                throw GotSignal(received_signum.load());
            }
        }
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

        double area_preferential = cell->state.area_preferential();
        if (num_increases > 0) {
            double dA = cell->state.area_preferential() / num_increases;
            for (unsigned int i = 0; i < num_increases; i++) {            
                vertex_model.increase_domain_size(dA);
                cell->state._area_preferential += dA;
                vertex_model.minimize_energy(minimization_after_increase);
            }

            if (am.area_of(cell) < threshold * cell->state.area_preferential())
            {
                throw std::runtime_error(fmt::format("Cell division failed! "
                    "Cell did not grow to area larger than threshold. "
                    "For division requested minimal area: {}. \n"
                    "For division preferred area: {}. \n"
                    "Area reached: {}.",
                    threshold * cell->state.area_preferential(),
                    cell->state.area_preferential(), am.area_of(cell)));
            }

            cell->state._area_preferential = area_preferential;
        }
        else {
            vertex_model.increase_domain_size(area_preferential);
        }

        double angle;
        if (normal_distr.stddev() > 1.e-11) {
            angle = normal_distr(*vertex_model.get_rng());
        }
        else {
            angle = uniform_distr(*vertex_model.get_rng());
        }
        vertex_model.divide_cell(cell, angle);
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


/// A pure shear experiment -- convergence and extension
OperationBundle build_pure_shear (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    double v0 = get_as<double>("shear_velocity", cfg);
    bool deform_plastic = get_as<bool>("deform_plastic", cfg);

    auto initialised = std::make_shared<bool>(false);
    auto rho = std::make_shared<double>(0.);
    auto Lx0 = std::make_shared<double>(0.);
    auto Ly0 = std::make_shared<double>(0.);

    Operation operation =
        [rho, v0, Lx0, Ly0,initialised, deform_plastic]
        (PCPVertex& vertex_model)
    {
        if (not vertex_model.get_space()->periodic) {
            throw std::runtime_error("Space needs to be periodic to apply "
                "`pure shear` operation");
        }

        using SpaceVec = typename PCPVertex::SpaceVec;
        auto space = vertex_model.get_space();

        if (not *initialised) {
            auto domain = space->get_domain_size();
            *Lx0 = domain[0];
            *Ly0 = domain[1];
            vertex_model.get_logger()->info("Initializing the domain shape "
                "for pure shear operation. Calling the pure ;shear operation "
                "again will overwrite changes.");
            
            *initialised = true;
        }

        auto domain = space->get_domain_size();
        double Lx = *Lx0 * exp(*rho);
        double Ly = *Ly0 * exp(-*rho);
        if (fabs(domain[0] - Lx) > 1.e-10 or fabs(domain[1] - Ly) > 1.e-10) {
            vertex_model.get_logger()->error("Vertex model domain: {}, {}",
                                             domain[0], domain[1]);
            vertex_model.get_logger()->error("Expected domain: {}, {}",
                                             Lx, Ly);
            throw std::runtime_error("Failed to apply `pure shear` operation, "
                "because the domain size was changed between two operations. "
                "Calling this operation would overwrite changes!");
        }

        *rho += v0;
        Lx = *Lx0 * exp(*rho);
        Ly = *Ly0 * exp(-*rho);

        const auto& cells = vertex_model.get_am().cells();
        double max_A0 = 0.;
        for (const auto& cell : cells) {
            max_A0 = std::max(max_A0, cell->state.area_preferential());
        }
        double L = 2 * sqrt(max_A0 / M_PI); // diameter of a cell

        if (Lx < 3 * L or Ly < 3 * L) {
            throw std::runtime_error("Cannot apply `pure shear` operation as "
                "the shorter tissue axis does not fit 3 diameters of a cell!");
        }

        vertex_model.stretch_domain(SpaceVec({Lx, Ly}) - domain,
                                    false, false, false,
                                    deform_plastic);
        // NOTE area preserving
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

/// The operation to relax preferential area
/** \details The following parameter are extracted from cfg 
 *           (besides those passed to `OperationParams`):
 *               - `factor` (double): incremental factor c. Should be in [0, 1].
 *               - `factor_type` (double, default: 0.): incremental factor c'.
 *                      Should be in [0, 1]. A factor that relaxes all cells
 *                      of one type to the average cell area
 *               - `factor_tissue` (double, default: 0.): incremental factor c''.
 *                      Should be in [0, 1]. A factor that relaxes all cells
 *                      to aim for a target average cell area.
 *               - `average_cell_area` (double, default: 1): The target
 *                      average cell area.
 *               - `minimum` (double, optional): Don't set preferential area 
 *                                               below this value.
 *               - `minimum` (double, optional): Don't set preferential area 
 *                                               above this value.
 *               - `minimum_hair` (double, optional): Limit to the preferential 
 *                          area of hair cells.
 *               - `maximum_hair` (double, optional): Limit to the preferential 
 *                          area of hair cells.
 *               - `minimum_support` (double, optional): Limit to the 
 *                          preferential area of cells.
 *               - `maximum_support` (double, optional): Limit to the 
 *                          preferential area of cells.
 * 
 *  Preferential area is relaxed towards the actual cell area:
 *  \f$ A_0(t+1) - A_0(t) =   c (A - A_0) + c\prime (<A>_{type} - A0)
 *                          + c\prime\prime (<A> - A_{target}) \f$
 * 
 *  Operation is applied to all cells with respective \f$A\f$ and \f$A_0\f$.
 */
OperationBundle build_relax_area (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    double factor(get_as<double>("factor", cfg));
    double factor_type(  get_as<double>("factor_type", cfg, 0.));
    double factor_tissue(get_as<double>("factor_tissue", cfg, 0.));
    double target_area(  get_as<double>("average_cell_area", cfg, 1.));

    double minimum(get_as<double>("minimum", cfg, 0.));
    double maximum(get_as<double>("maximum", cfg, 0.));

    double minimum_HC(get_as<double>("minimum_hair", cfg, 0.));
    double maximum_HC(get_as<double>("maximum_hair", cfg, 0.));

    double minimum_SC(get_as<double>("minimum_support", cfg, 0.));
    double maximum_SC(get_as<double>("maximum_support", cfg, 0.));

    Operation operation = [factor, factor_type, factor_tissue, target_area,
                           minimum, maximum,
                           minimum_HC, maximum_HC,
                           minimum_SC, maximum_SC]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        double average_HC = 0.;
        std::size_t num_HC = 0;
        double average_SC = 0.;
        std::size_t num_SC = 0;
        double average = 0.;
        for (const auto& cell : cells) {
            double area = am.area_of(cell);
            if (cell->state.type == CellType::hair) {
                average_HC += area;
                num_HC++;
            }
            else if (cell->state.type == CellType::support) {
                average_SC += area;
                num_SC++;
            }
            average += area;
        }
        average_HC /= static_cast<double>(num_HC); 
        average_SC /= static_cast<double>(num_SC); 
        average /= static_cast<double>(cells.size());

        PCPVertex::RuleFuncCell update = 
        [factor, factor_type, factor_tissue, target_area,
         average_HC, average_SC, average,
         minimum, maximum, minimum_HC, maximum_HC, minimum_SC, maximum_SC,
         am]
        (const auto& cell)
        {
            auto state = cell->state;
            double A0 = state.area_preferential();
            state._area_preferential += factor * (am.area_of(cell) - A0);
            if (state.type == CellType::hair) {
                state._area_preferential += factor_type * (average_HC - A0);
            }
            else if (state.type == CellType::support) {
                state._area_preferential += factor_type * (average_SC - A0);
            }
            state._area_preferential += factor_tissue * (target_area - average);

            // general maximum and minimum
            state._area_preferential = std::max(state._area_preferential,
                                                minimum);
            if (maximum > 1.e-12) {
                state._area_preferential = std::min(state._area_preferential,
                                                    maximum);
            }

            // hair cell maximum and minimum 
            if (state.type == CellType::hair) {
                state._area_preferential = std::max(state._area_preferential,
                                                    minimum_HC);
                if (maximum_HC > 1.e-12) {
                    state._area_preferential = 
                        std::min(state._area_preferential, maximum_HC);
                }
            }
            // hair cell maximum and minimum 
            else if (state.type == CellType::support) {
                state._area_preferential = std::max(state._area_preferential,
                                                    minimum_SC);
                if (maximum_SC > 1.e-12) {
                    state._area_preferential = 
                        std::min(state._area_preferential, maximum_SC);
                }
            }

            return state;
        };
        
        apply_rule<Update::sync>(update, cells);

        return;
    };

    return std::make_pair(operation, params);
}

/// The operation to set preferential area
/** \details The following parameter are extracted from cfg 
 *           (besides those passed to `OperationParams`):
 *               - `progenitor` (double, default: 0.): new value for
 *                      cells of type progenitor. Value 0. is ignored and
 *                      previous value of preferential area kept.
 *               - `stddev_progenitor` (double, default: 0.): stddev for cells
 *                      of type progenitor; using lognormal distribution
 *               - `hair` (double, default: 0.): new value for cells
 *                      of type hair. Value 0. is ignored and
 *                      previous value of preferential area kept.
 *               - `stddev_hair` (double, default: 0.): stddev for cells
 *                      of type hair; using lognormal distribution
 *               - `support` (double, default: 0.): new value for cells
 *                      of type support. Value 0. is ignored and
 *                      previous value of preferential area kept.
 *               - `stddev_support` (double, default: 0.): stddev for cells
 *                      of type support; using lognormal distribution
 *               - `relax_domain` (bool, default: false): If true, the domain
 *                      size is adapted to fit the cells as per preferential 
 *                      area.
 *               - `relax_domain_PD_axis` (bool, default: false): If true,
 *                      the domain size is adapted to fit the cells as per
 *                      preferential area. Domain size changes only in PD (x)
 *                      axis.
 */
OperationBundle build_set_area (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    double prog(get_as<double>("progenitor", cfg, 0.));
    double stddev_prog(get_as<double>("stddev_progenitor", cfg, 0.));
    double hair(get_as<double>("hair", cfg, 0.));
    double stddev_hair(get_as<double>("stddev_hair", cfg, 0.));
    double support(get_as<double>("support", cfg, 0.));
    double stddev_support(get_as<double>("stddev_support", cfg, 0.));
    bool relax_domain(get_as<bool>("relax_domain", cfg, false));
    bool relax_domain_PD_axis(get_as<bool>("relax_domain_PD_axis",
                                           cfg, false));

    if (relax_domain and relax_domain_PD_axis) {
        throw std::invalid_argument("In operation `set area` only one of cfgs "
            "`relax_domain` and `relax_domain_PD_axis` can be true!");
    }

    Operation operation = [prog, stddev_prog,
                           hair, stddev_hair,
                           support, stddev_support,
                           relax_domain, relax_domain_PD_axis]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        std::lognormal_distribution<> dist_prog, dist_hair, dist_support;
        if (prog > 1.e-12 and stddev_prog > 1.e-12) {
            dist_prog = get_lognormal_distribution(prog, stddev_prog);
        }
        if (hair > 1.e-12 and stddev_hair > 1.e-12) {
            dist_hair = get_lognormal_distribution(hair, stddev_hair);
        }
        if (support > 1.e-12 and stddev_support > 1.e-12) {
            dist_support = get_lognormal_distribution(support, stddev_support);
        }

        PCPVertex::RuleFuncCell update = [vertex_model,
                                          dist_prog{std::move(dist_prog)},
                                          dist_hair{std::move(dist_hair)},
                                          dist_support{std::move(dist_support)},
                                          prog, stddev_prog,
                                          hair, stddev_hair,
                                          support, stddev_support]
                (const auto& cell) mutable
        {
            std::function<double(double, double,
                                 std::lognormal_distribution<>&)> new_area = 
                [vertex_model](double mean, double stddev,
                               std::lognormal_distribution<>& distr)
            {
                if (stddev > 1.e-12) {
                    return distr(*vertex_model.get_rng());
                }
                else {
                    return mean;
                }
            };

            auto state = cell->state;
            double new_A = 0.;
            double mean = 0.;
            if (state.type == CellType::progenitor and prog > 1.e-12) {
                mean = prog;
                new_A = new_area(prog, stddev_prog, dist_prog);
            }
            else if (state.type == CellType::hair and hair > 1.e-12) {
                mean = hair;
                new_A = new_area(hair, stddev_hair, dist_hair);
            }
            else if (state.type == CellType::support and support > 1.e-12) {
                mean = support;
                new_A = new_area(support, stddev_support, dist_support);
            }

            if (mean > 1.e-12) {
                state._area_preferential = mean;
                state._area_preferential_fluctuations = new_A - mean;
            }

            return state;
        };
        
        apply_rule<Update::sync>(update, cells);

        if (relax_domain) {
            double area = std::accumulate(cells.begin(), cells.end(), 0.,
                            [](const double& val, const auto& cell) {
                                return val + cell->state.area_preferential();
                            });

            PCPVertex::SpaceVec domain_size = 
                    vertex_model.get_space()->get_domain_size();
            vertex_model.increase_domain_size(
                    area - domain_size[0]*domain_size[1]
            );
        }
        else if (relax_domain_PD_axis) {
            double area = std::accumulate(cells.begin(), cells.end(), 0.,
                            [](const double& val, const auto& cell) {
                                return val + cell->state.area_preferential();
                            });

            PCPVertex::SpaceVec domain_size = 
                    vertex_model.get_space()->get_domain_size();
            double dA = area - domain_size[0]*domain_size[1];

            vertex_model.stretch_domain({dA/domain_size[1], 0.},
                                        false, false, false);
        }

        return;
    };

    return std::make_pair(operation, params);
}

OperationBundle build_update_boundary_parameter (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    auto boundary_cfg = get_as<Config>("boundary_parameter", cfg);

    Operation operation = [boundary_cfg] (PCPVertex& vertex_model)
    {
        vertex_model.update_boundary_parameter(boundary_cfg);
    };

    return std::make_pair(operation, params);
}

/// Set a torque to cells
/** The configuration is passed on to the constructor of RotationCellState
 *  The config following the keys to the respective cell type:
 *      - `hair`
 *      - `support`
 *      - `progenitor`
 * 
 *  If a configuration is not given, torque zero is set.
 */
OperationBundle build_set_torque (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    Config cfg_default;
    cfg_default["torque"] = 0.;

    Config cfg_hair = get_as<Config>("hair", cfg, cfg_default);
    Config cfg_support = get_as<Config>("support", cfg, cfg_default);
    Config cfg_progenitor = get_as<Config>("progenitor", cfg, cfg_default);

    Operation operation = [cfg_hair, cfg_support, cfg_progenitor]
                          (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        const auto& cells = vertex_model.get_am().cells();
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::hair) {
                auto rot_state = std::make_shared<RotationCellState>(cfg_hair);
                cell->custom_links().rotation_state = rot_state;
            }
            else if (cell->state.type == CellType::support) {
                auto rot_state = std::make_shared<RotationCellState>(
                                                cfg_support);
                cell->custom_links().rotation_state = rot_state;
            }
            else if (cell->state.type == CellType::progenitor) {
                auto rot_state = std::make_shared<RotationCellState>(
                                                cfg_progenitor);
                cell->custom_links().rotation_state = rot_state;
            }
        }
    };

    return std::make_pair(operation, params);
}

/// A simple shear model
/** Move all boundary vertices in a simple shear way:
 *  \f$ dx = const * y \f$
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `max_shear` (double): the maximum displacement at \f$ y_{max} \f$
 *                              with const = max_shear / W, where H the width
 *                              of the tissue
 */
OperationBundle build_simple_shear (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    using SpaceVec = PCPVertex::SpaceVec;

    OperationParams params(name, cfg, default_minim_params);

    SpaceVec v0 = get_as_SpaceVec<2>("shear_velocity", cfg);
    bool deform_plastic = get_as<bool>("deform_plastic", cfg);

    Operation operation = [v0, deform_plastic](PCPVertex& vertex_model)
    {
        const auto& space = vertex_model.get_space();

        SpaceVec skew = space->get_skew();
        SpaceVec domain = space->get_domain_size();

        skew += v0;

        space->set_skew(skew);

        const auto& am = vertex_model.get_am();
        if (deform_plastic) {
            if (fabs(space->get_curvature()) > 1.e-12) {
                double curvature = space->get_curvature();
                for (const auto& vertex : am.vertices()) {
                    SpaceVec pos = am.position_of(vertex);
                    auto [rho, theta] = space->transform_radial(pos);
                    double skew_theta = v0[0] * curvature;
                    double max_theta = (domain[0] / 2.) * curvature;

                    double r = (rho - 1. / curvature) / domain[1] + 0.5;
                    double r_theta = theta / (2 * max_theta);

                    theta += skew_theta * r;
                    rho += v0[1] * (r_theta + max_theta);

                    pos = space->transform_cartesian(rho, theta);
                    am.move_to(vertex, pos);
                }
            }
            else {
                for (const auto& vertex : am.vertices()) {
                    SpaceVec rpos = am.position_of(vertex) / domain;
                    am.move_by(vertex, v0 % SpaceVec({rpos[1], rpos[0]}));
                }
            }
        }
    };

    return std::make_pair(operation, params);
}


/// A bending model in periodic boundary conditions
OperationBundle build_bending_box_bc (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    using SpaceVec = PCPVertex::SpaceVec;

    OperationParams params(name, cfg, default_minim_params);

    double v0 = get_as<double>("increment_curvature", cfg);
    bool deform_plastic = get_as<bool>("deform_plastic", cfg);

    Operation operation = [v0, deform_plastic](PCPVertex& vertex_model)
    {
        const auto& space = vertex_model.get_space();

        double curvature = space->get_curvature();
        curvature += v0;

        space->set_curvature(curvature);

        if (deform_plastic) {
            const SpaceVec domain = space->get_domain_size();

            double k0 = std::max(curvature - v0, 1.e-8);
            double R0 = 1. / k0;

            double k1 = std::max(curvature, 1.e-8);
            double R1 = 1. / k1;
            
            double max_theta_0 = (domain[0] / 2.) * k0;
            double max_theta_1 = (domain[0] / 2.) * k1;

            SpaceVec origin_0({domain[0] / 2., domain[1] / 2. - R0});
            SpaceVec origin_1({domain[0] / 2., domain[1] / 2. - R1});

            const auto& am = vertex_model.get_am();
            for (const auto& v : am.vertices()) {
                SpaceVec pos = am.position_of(v);
                SpaceVec displ = pos - origin_0;

                // relative radial position
                double rho = arma::norm(displ) - R0; 
                // relative angular position
                double theta = std::atan(displ[0] / displ[1]) / max_theta_0;

                SpaceVec new_pos = (
                    origin_1 + (rho + R1) * SpaceVec({sin(theta * max_theta_1), 
                                                      cos(theta * max_theta_1)})
                );

                am.move_to(v, new_pos);
            }
        }
    };

    return std::make_pair(operation, params);
}

} // namespace OperationCollection
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif