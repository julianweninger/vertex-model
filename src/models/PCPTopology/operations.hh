#ifndef UTOPIA_MODELS_PCPTOPOLOGY_OPERATIONS_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_OPERATIONS_HH

#include <typeinfo>
#include <algorithm>
#include <iterator>

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

        return iters * minimization_params.num_repeat;
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

    if (probability > 0 and fraction > 0) {
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

    if (fraction > 0) {
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

/// An operation to fix boundary vertices in space
OperationBundle build_fix_boundary (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    Operation operation = [] (PCPVertex& vertex_model)
    {
        const auto& am = vertex_model.get_am();
        apply_rule<Update::sync>(
            [am](const auto& vertex) {
                auto state = vertex->state;
                if (am.is_boundary(vertex)) {
                    state.fix_in_space = true;
                }
                else {
                    state.fix_in_space = false;
                }
                return state;
            },
            am.vertices()
        );
    };

    return std::make_pair(operation, params);
}

/// A convergence and extension model
/** The outermost vertices in the vertical direction are moved towards the
 *  horizontal tissue axis and fixed in space for minimization.
 *  Horizontal boundary vertices are free to move and are thought to move
 *  outwards to compensate increased pressure.
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `dH` (double): The length by which the vertical axis is reduced
 */
OperationBundle build_convergence_and_extension (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    using SpaceVec = PCPVertex::SpaceVec;

    OperationParams params(name, cfg, default_minim_params);
    double dH = get_as<double>("dH", cfg);
    Operation operation = [dH] (PCPVertex& vertex_model)
    {
        const auto& am = vertex_model.get_am();

        double min_y = std::numeric_limits<double>::max();
        double max_y = std::numeric_limits<double>::min();

        for (const auto& vertex : am.vertices()) {
            SpaceVec pos = am.position_of(vertex);
            min_y = std::min(min_y, pos[1]);
            max_y = std::max(max_y, pos[1]);
        }
        
        // move the outermost vertices by dH / 2. towards hor. axis
        // fix moved vertices permanently in space 
        apply_rule<Update::sync>(
            [am, dH, min_y, max_y] (const auto& vertex)
            {
                auto state = vertex->state;
                SpaceVec pos = am.position_of(vertex);
                if (fabs(pos[1] - min_y) < 3 * dH) {
                    am.move_by(vertex, SpaceVec({0.,  dH / 2.}));
                    state.fix_in_space = true;
                }
                else if (fabs(pos[1] - max_y) < 3 * dH) {
                    am.move_by(vertex, SpaceVec({0., -dH / 2.}));
                    state.fix_in_space = true;
                }
                
                return state;
            },
            am.vertices()
        );
    };

    return std::make_pair(operation, params);
}

/// The operation to increment preferential area
/** \details The following parameter are extracted from cfg 
 *           (besides those passed to `OperationParams`):
 *               - `progenitor` (double, default: 0.): incremental value for
 *                      cells of type progenitor
 *               - `stddev_progenitor` (double, default: 0.): stddev for cells
 *                      of type progenitor; using normal distribution
 *               - `hair` (double, default: 0.): incremental value for cells
 *                      of type hair
 *               - `stddev_hair` (double, default: 0.): stddev for cells
 *                      of type hair; using normal distribution
 *               - `support` (double, default: 0.): incremental value for cells
 *                      of type support
 *               - `stddev_support` (double, default: 0.): stddev for cells
 *                      of type support; using normal distribution
 *               - `adapt_support` (bool, default: false): If true, the total
 *                      area of hair and support cells remains constant, hence
 *                      support cells compensate area changes from hair cells.
 *               - `relax_domain` (bool, default: false): If true, the domain
 *                      size is adapted to fit the cells as per preferential 
 *                      area.
 * 
 *  \note This affects the current entities properties, it does not overwrite
 *        changes in the past.
 */
OperationBundle build_increment_area (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    double incr_prog(get_as<double>("progenitor", cfg, 0.));
    double stddev_prog(get_as<double>("stddev_progenitor", cfg, 0.));
    double incr_hair(get_as<double>("hair", cfg, 0.));
    double stddev_hair(get_as<double>("stddev_hair", cfg, 0.));
    double incr_support(get_as<double>("support", cfg, 0.));
    double stddev_support(get_as<double>("stddev_support", cfg, 0.));
    bool adapt_support(get_as<bool>("adapt_support", cfg, false));
    bool relax_domain(get_as<bool>("relax_domain", cfg, false));
    
    if (adapt_support and incr_support != 0) {
        throw std::invalid_argument(fmt::format(
            "In cfg {}: If `adapt_support = True`, "
            "then the increment value for the support cells "
            "`support` must be zero, but was {}!", name, incr_support));
    }

    if (adapt_support and relax_domain) {
        throw std::invalid_argument(fmt::format(
            "In cfg {}: `adapt_support` and `relax_domain` cannot both be "
            "true!", name));
    }

    Operation operation = [incr_prog, stddev_prog,
                           incr_hair, stddev_hair,
                           incr_support, stddev_support, adapt_support,
                           relax_domain]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        std::normal_distribution<> dist_prog{incr_prog, stddev_prog};
        std::normal_distribution<> dist_hair{incr_hair, stddev_hair};
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
                                                      stddev_support);
        }
        else {
            dist_support = std::normal_distribution<>(incr_support,
                                                      stddev_support);
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

        if (not relax_domain) {
            return;
        }

        double area = std::accumulate(cells.begin(), cells.end(), 0.,
                            [](const double& val, const auto& cell) {
                                return val + cell->state.area_preferential; });

        PCPVertex::SpaceVec domain_size = 
                vertex_model.get_space()->get_domain_size();
        vertex_model.increase_domain_size(area - domain_size[0]*domain_size[1]);

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
 *  \note This increments the matrix of linetension, hence affects current and
 *        future edges between cells of corresponding type.
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
        contractility(CellType::support, CellType::hair) += incr_hair_supp;
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
        linetension(CellType::support, CellType::hair) += incr_hair_supp;
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
    MinimizationParams minimization_after_increase(
        get_as<Config>("minimization_after_increase", cfg, Config()),
        params.minimization_params);

    auto threshold(get_as<double>("area_threshold", cfg, 0.5));
    
    auto generation_max_id = std::make_shared<std::size_t>(0);
    
    std::uniform_real_distribution<double> prob_distr(0.,1.);

    Operation operation = [num_increases, minimization_after_increase,
                           threshold, generation_max_id,
                           prob_distr{std::move(prob_distr)}, params]
            (PCPVertex& vertex_model) mutable
    {
        using Cell = typename PCPVertex::Cell;

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        AgentContainer<Cell> current_generation(cells.size());
        auto it = std::copy_if (cells.begin(), cells.end(),
                                current_generation.begin(),
                                [generation_max_id, am](const auto& cell){
                                    return (    cell->id() < *generation_max_id
                                            and not am.is_boundary(cell));
                                } );
        current_generation.resize(
            std::distance(current_generation.begin(), it));
        if (current_generation.size() == 0) {
            for (const auto& c : cells) {
                if (am.is_boundary(c)) { continue; }
                *generation_max_id = std::max(*generation_max_id, c->id());
            }
            *generation_max_id = *generation_max_id + 1;
            
            current_generation.clear();
            current_generation.resize(cells.size());
            auto it = std::copy_if (cells.begin(), cells.end(),
                                    current_generation.begin(),
                                    [generation_max_id, am](const auto& cell){
                                        return (cell->id() < *generation_max_id
                                            and not am.is_boundary(cell));
                                    } );
            current_generation.resize(
                std::distance(current_generation.begin(), it));
        }

        std::uniform_int_distribution<> int_dist(
            0, current_generation.size() - 1);
        
        auto cell = current_generation[int_dist(*vertex_model.get_rng())];

        double area_preferential = cell->state.area_preferential;
        if (num_increases > 0) {
            double dA = cell->state.area_preferential / num_increases;
            for (unsigned int i = 0; i < num_increases; i++) {            
                vertex_model.increase_domain_size(dA);
                cell->state.area_preferential += dA;
                vertex_model.minimize_energy(minimization_after_increase);
            }

            if (am.area_of(cell) < threshold * cell->state.area_preferential)
            {
                throw std::runtime_error(fmt::format("Cell division failed! "
                    "Cell did not grow to area larger than threshold. "
                    "For division requested minimal area: {}. \n"
                    "For division preferred area: {}. \n"
                    "Area reached: {}.",
                    threshold * cell->state.area_preferential,
                    cell->state.area_preferential, am.area_of(cell)));
            }

            cell->state.area_preferential = area_preferential;
        }
        else {
            vertex_model.increase_domain_size(area_preferential);
        }

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

/// The operation to relax preferential area
/** \details The following parameter are extracted from cfg 
 *           (besides those passed to `OperationParams`):
 *               - `factor` (double): incremental factor c. Should be in [0, 1].
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
 *  \f$ A^\prime_0  = A_0 + c (A - A_0) \f$
 * 
 *  Operation is applied to all cells with respective \f$A\f$ and \f$A_0\f$.
 */
OperationBundle build_relax_area (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    double factor(get_as<double>("factor", cfg));
    double minimum(get_as<double>("minimum", cfg, 0.));
    double maximum(get_as<double>("maximum", cfg, 0.));


    double minimum_HC(get_as<double>("minimum_hair", cfg, 0.));
    double maximum_HC(get_as<double>("maximum_hair", cfg, 0.));

    double minimum_SC(get_as<double>("minimum_support", cfg, 0.));
    double maximum_SC(get_as<double>("maximum_support", cfg, 0.));

    Operation operation = [factor, minimum, maximum, minimum_HC, maximum_HC,
                           minimum_SC, maximum_SC]
            (PCPVertex& vertex_model)
    {
        using CellType = PCPVertex::CellType;

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        PCPVertex::RuleFuncCell update = [factor,
                                          minimum, maximum, 
                                          minimum_HC, maximum_HC,
                                          minimum_SC, maximum_SC,
                                          am](const auto& cell)
        {
            auto state = cell->state;
            state.area_preferential += factor * (  am.area_of(cell)
                                                 - state.area_preferential);

            // general maximum and minimum
            state.area_preferential = std::max(state.area_preferential,
                                               minimum);
            if (maximum > 0.) {
                state.area_preferential = std::min(state.area_preferential,
                                                   maximum);
            }

            // hair cell maximum and minimum 
            if (state.type == CellType::hair) {
                state.area_preferential = std::max(state.area_preferential,
                                                   minimum_HC);
                if (maximum_HC > 0.) {
                    state.area_preferential = std::min(state.area_preferential,
                                                       maximum_HC);
                }
            }
            // hair cell maximum and minimum 
            else if (state.type == CellType::support) {
                state.area_preferential = std::max(state.area_preferential,
                                                   minimum_SC);
                if (maximum_SC > 0.) {
                    state.area_preferential = std::min(state.area_preferential,
                                                       maximum_SC);
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

        std::function<std::pair<double, double>(double, double)> m_s = 
            [](double mean, double stddev)
        {
            if (mean < 0.) {
                throw std::invalid_argument(fmt::format("Cannot construct "
                    "lognormal distribution with mean {} <= 0!", stddev));
            }
            if (stddev <= 0.) {
                throw std::invalid_argument(fmt::format("Cannot construct "
                    "lognormal distribution with stddev {} <= 0!", stddev));
            }
            double m = log(  std::pow(mean, 2)
                           / sqrt(std::pow(stddev, 2) + std::pow(mean, 2)));
            double s = sqrt(log(std::pow(stddev, 2) / std::pow(mean, 2) + 1));
            return std::make_pair(m, s);
        };

        std::lognormal_distribution<> dist_prog, dist_hair, dist_support;
        if (prog > 0.0 and stddev_prog > 0.) {
            auto [m, s] = m_s(prog, stddev_prog);
            dist_prog = std::lognormal_distribution<>(m, s);
        }
        if (hair > 0.0 and stddev_hair > 0.) {
            auto [m, s] = m_s(hair, stddev_hair);
            dist_hair = std::lognormal_distribution<>(m, s);
        }
        if (support > 0.0 and stddev_support > 0.) {
            auto [m, s] = m_s(support, stddev_support);
            dist_support = std::lognormal_distribution<> (m, s);            
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
                if (stddev > 0.) {
                    return distr(*vertex_model.get_rng());
                }
                else {
                    return mean;
                }
            };

            auto state = cell->state;
            if (state.type == CellType::progenitor and prog > 0.) {
                state.area_preferential = new_area(prog, stddev_prog,
                                                   dist_prog);
            }
            else if (state.type == CellType::hair and hair > 0.) {
                state.area_preferential = new_area(hair, stddev_hair,
                                                   dist_hair);
            }
            else if (state.type == CellType::support and support > 0.) {
                state.area_preferential = new_area(support, stddev_support,
                                                   dist_support);
            }
            return state;
        };
        
        apply_rule<Update::sync>(update, cells);

        if (relax_domain) {
            double area = std::accumulate(cells.begin(), cells.end(), 0.,
                            [](const double& val, const auto& cell) {
                                return val + cell->state.area_preferential; });

            PCPVertex::SpaceVec domain_size = 
                    vertex_model.get_space()->get_domain_size();
            vertex_model.increase_domain_size(
                    area - domain_size[0]*domain_size[1]
            );
        }
        else if (relax_domain_PD_axis) {
            double area = std::accumulate(cells.begin(), cells.end(), 0.,
                            [](const double& val, const auto& cell) {
                                return val + cell->state.area_preferential; });

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

OperationBundle build_set_boundary_parameter (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    auto boundary_cfg = get_as<Config>("boundary_parameter", cfg);

    Operation operation = [boundary_cfg] (PCPVertex& vertex_model)
    {
        vertex_model.set_boundary_parameter(boundary_cfg);
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
    double max_shear = get_as<double>("max_shear", cfg);

    auto [op_fix_bc, __params_fix_bc] = build_fix_boundary("fix_boundary", cfg,
                                             default_minim_params);

    Operation operation = [max_shear, op_fix_bc] (PCPVertex& vertex_model)
    {
        const auto& am = vertex_model.get_am();

        // fix the current boundary cells
        op_fix_bc(vertex_model);

        double min_y = std::numeric_limits<double>::max();
        double max_y = std::numeric_limits<double>::min();

        for (const auto& vertex : am.vertices()) {
            SpaceVec pos = am.position_of(vertex);
            min_y = std::min(min_y, pos[1]);
            max_y = std::max(max_y, pos[1]);
        }

        double gradient = max_shear / (max_y - min_y);
        
        // move the outermost vertices by dH / 2. towards hor. axis
        // fix moved vertices permanently in space 
        apply_rule<Update::sync>(
            [am, gradient, min_y] (const auto& vertex)
            {
                SpaceVec pos = am.position_of(vertex);
                double y = pos[1] - min_y;
                if (vertex->state.fix_in_space) {
                    am.move_by(vertex, SpaceVec({gradient * y, 0.}));
                }
                
                return vertex->state;
            },
            am.vertices());
    };

    return std::make_pair(operation, params);
}


} // namespace OperationCollection
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif