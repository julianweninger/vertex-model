#ifndef UTOPIA_MODELS_PCPTOPOLOGY_OPERATIONS_DIFFERENTIATION_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_OPERATIONS_DIFFERENTIATION_HH

#include <utopia/core/types.hh>

#include "../PCPVertex/utils.hh"
#include "../PCPVertex/PCPVertex.hh"

#include "operations.hh"


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace OperationCollection {


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
            else {
                c_cells[iterator]->state.cell_type = CCellType(type);
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

            cell->state.type = int(c_cell->state.cell_type);
        }
    };

    return std::make_pair(operation, params);
}

/// The operation to differentiate the types of cells into regular domains
/** Along horizontal (width) and vertical (height) axis straight domain walls
 *  will be differentiated, according to cell centers. Lower left domain will be
 *  a rectangle, others will be 2 rectangles combined at 90 degree angle.
 *  Using relative coordinates.
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `relative_widths` (list of double): The widths of vertical domains.
 *             Sum of widhts needs to be 1.
 *      - `relative_heights` (list of double): The widths of horizontal domains.
 *             Sum of widhts needs to be 1.
 *      - `specify_types` (list of int): Customized types in order of appearance.
 *              If not provided, range(N) used.
 */
OperationBundle build_differentiate_domain (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);
    
    auto _xs(get_as<std::vector<double>>("relative_widths", cfg));
    auto _ys(get_as<std::vector<double>>("relative_heights", cfg));
    auto types(get_as<std::vector<double>>("specify_types", cfg, {}));

    std::vector<double> xs{};
    std::vector<double> ys{};

    double sum_x = 0.;
    for (const auto& x : _xs) {
        sum_x += x;
        xs.push_back(sum_x);
    }
    double sum_y = 0.;
    for (const auto& y : _ys) {
        sum_y += y;
        ys.push_back(sum_y);
    }

    if (fabs(sum_x - 1.) + fabs(sum_y - 1.) > 1.e-6) {
        throw std::runtime_error(fmt::format("The sum of `relative_widths` and "
            "`relative_heights`, respectively, needs to be 1, but was {} and {}, resp.",
            sum_x, sum_y));
    }

    Operation operation = [xs, ys, types](PCPVertex& vertex_model)
    {
        using SpaceVec = PCPVertex::SpaceVec;
        
        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();

        double Lx, Ly;
        double min_x, min_y;
        if (vertex_model.get_space()->periodic) {
            SpaceVec domain = vertex_model.get_space()->get_domain_size();
            Lx = domain[0];
            Ly = domain[1];
            min_x = 0.;
            min_y = 0.;
        }
        else {
            min_x = std::numeric_limits<double>::max();
            double max_x = std::numeric_limits<double>::min();
            min_y = std::numeric_limits<double>::max();
            double max_y = std::numeric_limits<double>::min();
            for (const auto& vertex : am.vertices()) {
                SpaceVec pos = am.position_of(vertex);
                min_x = std::min(min_x, pos[0]);
                max_x = std::max(max_x, pos[0]);

                min_y = std::min(min_y, pos[1]);
                max_y = std::max(max_y, pos[1]);
            }
            Lx = (max_x - min_x) - min_x;
            Ly = (max_y - min_y) - max_x;
        }

        for (const auto& cell : cells) {
            SpaceVec rpos = am.barycenter_of(cell) - SpaceVec({min_x, min_y});
            rpos /= SpaceVec({Lx, Ly});

            std::size_t type = 0;
            while(   rpos[0] > xs[std::min(type, xs.size()-1)]
                  or rpos[1] > ys[std::min(type, ys.size()-1)])
            {
                type += 1;
            }
            if (types.size() == 0) {
                cell->state.type = type;
            }
            else if (type < types.size()) {
                cell->state.type = types[type];
            }
            else {
                throw std::runtime_error(fmt::format("Not enough types "
                    "specified in differentiate domain! "
                    "{} types were specified as `specify_types'.", 
                    types.size()
                ));
            }
        }
    };

    return std::make_pair(operation, params);
}


/// The operation to differentiate cells in contact of type 1, surrounded type 0
/** \note mainly for testing purposes, i.e. which mechanism can separate two
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
OperationBundle build_differentiate_cluster (
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
                cell->state.type = 1;
                return cell->state;
            };

            apply_rule<Update::sync>(differentiate, cluster);
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

        std::size_t iterator;
        for (iterator = 0; iterator < cells.size(); iterator++) {
            cells[iterator]->custom_links().nd_cell = nd_cells[iterator];
            
            auto type = cells[iterator]->state.type;
            if (not *prolog and not discard_initialisation) {
                void(); // use initialization of ND Model
            }
            nd_cells[iterator]->state.cell_type = NDCellType(type);
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
            cell->state.type = int(type);
        }
    };

    return std::make_pair(operation, params);
}


/// The operation to differentiate the types of cells randomly
/** Will assign states randomly, given either a distribution or a fraction
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `probability` (double): The probability for a cell to differentiate
 *              to a hair cell; otherwise support cell. Gives a fraction of 
 *              hair cells for large enough number of cells.
 *      - `fraction` (double): Similar to `probability`, but a fraction of the 
 *              cells is selected for HC differentiation. Resulting HC fraction
 *              is fixed. Number of cells is rounded to next smaller integer.
 *      - `overwrite` (bool, default: true): Whether non-zero types will be overwritten
 * 
 *  \note The entities properties do not change during differentiation.
 */
OperationBundle build_differentiate_random (
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    auto _distribution(get_as<std::vector<double>>("distribution", cfg));
    auto _method(get_as<std::string>("method", cfg));
    bool overwrite(get_as<bool>("overwrite", cfg, true));


    std::vector<double> distribution({});
    double sum = 0.;
    for (auto f : _distribution) {
        sum += f;
        distribution.push_back(sum);
    }

    if (distribution.size() < 2) {
        throw std::invalid_argument(fmt::format(
            "At least 2 values have to be provided in `distribution`, but got "
            "{}!",
            distribution.size()
        ));
    }
    if (fabs(1 - sum) > 1.e-8) {
        throw std::invalid_argument(fmt::format(
            "The sum over the `distribution` of cell types needs to be 1, but "
            "was {}!",
            sum
        ));
    }
    

    enum Method {
        Fraction,
        Probability
    } method;
    if (_method == "fraction") {
        method = Method::Fraction;
    }
    else if (_method == "probability") {
        method = Method::Probability;
    }
    else {
        throw std::invalid_argument(fmt::format(
            "The `method` for cell type random distribution has to be one of "
            "the following, but was {}:\n{}",
            _method,
            "- fraction\n"
            "- probability\n"
        ));
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


    Operation operation;

    if (method == Method::Fraction) {
        operation = [distribution, overwrite] (PCPVertex& vertex_model)
        {
            auto cells = vertex_model.get_am().cells();
            if (not overwrite) {
                cells.erase(
                    std::remove_if(
                        cells.begin(), cells.end(), 
                        [](const auto& cell) {
                            return cell->state.type != 0;
                        }
                    ),
                    cells.end()
                );
            }
            std::shuffle(cells.begin(), cells.end(), *vertex_model.get_rng());

            std::size_t type = 0;
            for (std::size_t i = 0; i < cells.size(); i++) {
                while (i >= std::ceil(distribution[type] * cells.size())) {
                    type += 1;
                }
                cells[i]->state.type = type;
            }
        };
    }
    else if (method == Method::Probability) {
        operation = [distribution, overwrite] (PCPVertex& vertex_model)
        {
            std::uniform_real_distribution<double> prob_distr(0., 1.);

            auto cells = vertex_model.get_am().cells();
            if (not overwrite) {
                cells.erase(
                    std::remove_if(
                        cells.begin(), cells.end(), 
                        [](const auto& cell) {
                            return cell->state.type != 0;
                        }
                    ),
                    cells.end()
                );
            }
            for (const auto& cell : cells) {
                double rn = prob_distr(*vertex_model.get_rng());
                std::size_t type = 0;
                while (rn > distribution[type]) {
                    type += 1;
                }
                cell->state.type = type;
            }
        };
    }
    else {
        throw std::runtime_error("Unknown Method!");
    }

    return std::make_pair(operation, params);
}


/// The operation minimizes cell contacts between same type cells
/** Iteratively, for cells of specified type, the type of the cell is swapped 
 *  with the neighbor that has least number of max type neighbours.
 * 
 *  The following parameter are extracted from cfg 
 *  (besides those passed to `OperationParams`):
 *      - `num_steps` (uint): How many steps to iterate the algorithm. In every
 *              step one label is compared with its neighbors.
 *      - `type` (uint): The cell type within which to reduce number of contacts
 *      - `neutral_probability` (double): If swapping is neutral, do with a
 *              probability anyways.
 */
OperationBundle build_minimize_cell_contacts(
        std::string name, const Config& cfg,
        const MinimizationParams& default_minim_params)
{
    OperationParams params(name, cfg, default_minim_params);

    auto num_steps(get_as<std::size_t>("num_steps", cfg));
    auto type(get_as<std::size_t>("type", cfg));
    auto probability(get_as<double>("neutral_probability", cfg));

    Operation operation = [num_steps, type, probability]
    (PCPVertex& vertex_model)
    {
        const auto& am = vertex_model.get_am();
        auto cells = am.cells();

        cells.erase(
            std::remove_if(
                cells.begin(), cells.end(), 
                [type](const auto& c) { return c->state.type != type; }
            ),
            cells.end()
        );

        if (cells.size() == 0) {
            vertex_model.get_logger()->warn("No cells of type {} in tissue. "
                "Aborting minimize contacts of this type.", type);
            return;
        }

        std::uniform_int_distribution<std::size_t> distr(0, cells.size() - 1);
        std::uniform_real_distribution<double> _distr(0, 1.);

        for (std::size_t i = 0; i < num_steps; i++) {
            auto cell = cells[distr(*vertex_model.get_rng())];
            if (cell->state.type != type) { throw; }

            auto nbs = am.neighbors_of(cell);
            std::size_t cnt = std::count_if(
                nbs.begin(), nbs.end(),
                [type](const auto& n) { return n->state.type == type; }
            );

            if (cnt == 0 and probability < 1.e-6) {
                continue;
            }

            std::size_t min_n = cnt+1;
            std::vector<std::shared_ptr<PCPVertex::Cell>> candidates{};
            for (std::size_t n = 0; n < nbs.size(); n++) {
                if (nbs[n]->state.type == type) {
                    continue;
                }
                auto _nbs = am.neighbors_of(nbs[n]);
                std::size_t _cnt = std::count_if(
                    _nbs.begin(), _nbs.end(),
                    [type](const auto& n) { return n->state.type == type; }
                );

                if (_cnt < min_n) {
                    min_n = _cnt;
                    candidates.clear();
                    candidates.push_back(nbs[n]);
                }
                // Center cell is counted in _cnt, hence equal N is better to swap
                else if (_cnt == min_n) {
                    candidates.push_back(nbs[n]);
                }
            }

            if (candidates.size() > 0) {
                std::shuffle(
                    candidates.begin(), candidates.end(),
                    *vertex_model.get_rng()
                );
                // if (   min_n < cnt
                //     or _distr(*vertex_model.get_rng()) < probability)
                // {
                    std::swap(cell->state.type, candidates[0]->state.type);
                    cells.erase(
                        std::remove(cells.begin(), cells.end(), cell),
                        cells.end()
                    );
                    cells.push_back(candidates[0]);
                // }
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