#ifndef UTOPIA_MODELS_PCPVERTEX_OPERATIONS_HH
#define UTOPIA_MODELS_PCPVERTEX_OPERATIONS_HH

#include <typeinfo>

namespace Utopia {
namespace Models {
namespace PCPVertex {

/// Apply a perturbation to the position of vertices
/** Move the x and y position by a random value in [-intensity, intensity]
 *  using a uniform distribution.
 * 
 * TODO write test
 */
void PCPVertex::jiggle_vertices(double intensity)
{
    this->_log->debug("Jiggling the vertices on a length scale of "
                        "{} ..", intensity);
    
    const RuleFuncVertex jiggle = [this, intensity] (const auto& vertex)
    {
        SpaceVec displ = SpaceVec({this->_prob_distr(*this->_rng),
                                   this->_prob_distr(*this->_rng)}) *
                         2. * intensity - intensity;
        this->_am.move_by(vertex, displ);
        return vertex->state;
    };

    apply_rule<Update::sync>(jiggle, _am.vertices());

    this->init_minimisation();
};

/// Helper for differentiation of progenitor cells
/** \param linetension  The symmetric matrix of linetension interactions
 *                      between two cells of same or different type
 *  \param edge_contractility  The symmetric matrix of contractility interactions
 *                             between two cells of same or different type
 *  \param area_preferential    The preferential cell area of the different
 *                              cell types
 */
void PCPVertex::differentiate_hair_cells_hlpr(
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> linetension,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> edge_contractility,
        arma::Col<double>::fixed<CellType::num_cell_types> area_preferential)
{
    // check parameter
    for (int i = 0; i < CellType::num_cell_types; i++) {
        for (int j = i+1; j < CellType::num_cell_types; j++) {
            if (linetension(j, i) == linetension(i, j)) {
                continue;
            }
            
            this->_log->warn("Invalid argument in differentiate_hair_cells."
                "Got non symmetric 'linetension' matrix!");
            throw std::invalid_argument("Non symmetric 'linetension'"
                "matrix");
        }
    }
    for (int i = 0; i < CellType::num_cell_types; i++) {
        for (int j = i+1; j < CellType::num_cell_types; j++) {
            if (edge_contractility(j, i) == edge_contractility(i, j)) {
                continue;
            }
            
            this->_log->warn("Invalid argument in differentiate_hair_cells."
                "Got non symmetric 'edge_contractility' matrix!");
            throw std::invalid_argument("Non symmetric 'edge_contractility'"
                "matrix");
        }
    }

    const RuleFuncCell update_area_pref = [this, area_preferential] (
            const auto& cell)
    {
        auto state = cell->state;
        state.area_preferential = area_preferential(state.type);
        return state;
    };

    _linetension = linetension;
    _edge_contractility = edge_contractility;

    const RuleFuncEdge set_edge_prop = [this] (
            const auto& edge)
    {
        const auto& [adj_cell_a, adj_cell_b] = _am.adjoints_of(edge);
        CellType cell_a_type, cell_b_type;
        if (adj_cell_a) {
            cell_a_type = adj_cell_a->state.type;
        }
        else {
            cell_a_type = CellType::support;
        }
        if (adj_cell_b) {
            cell_b_type = adj_cell_b->state.type;
        }
        else {
            cell_b_type = CellType::support;
        }

        auto state = edge->state;
        state.linetension = this->_linetension(cell_a_type, cell_b_type);
        state.contractility = this->_edge_contractility(cell_a_type,
                                                        cell_b_type);

        return state;
    };

    apply_rule<Update::sync>(update_area_pref, _am.cells());
    apply_rule<Update::sync>(set_edge_prop, _am.edges());
}

/// Differentiates progenitor cells with random hair cell distribution
/** \param fraction     fraction of hair cells. Others are support cells
 *  \param linetension  The symmetric matrix of linetension interactions
 *                      between two cells of same or different type
 *  \param edge_contractility  The symmetric matrix of contractility interactions
 *                             between two cells of same or different type
 *  \param area_preferential    The preferential cell area of the different
 *                              cell types
 */
void PCPVertex::differentiate_hair_cells_random(
        double fraction,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> linetension,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> edge_contractility,
        arma::Col<double>::fixed<CellType::num_cell_types> area_preferential)
{
    this->_log->info("Differentiating progenitor cells to {}% hair cells "
        "and {}% support cells with uniform spatial distribution ...",
        fraction, 1-fraction);

    

    const RuleFuncCell set_type_rand = [this, fraction] (
            const auto& cell)
    {
        auto state = cell->state;
        if (this->_prob_distr(*this->_rng) < fraction) {
            state.type = CellType::hair;
        }
        else {
            state.type = CellType::support;
        }
        return state;
    };

    apply_rule<Update::sync>(set_type_rand, _am.cells());

    differentiate_hair_cells_hlpr(linetension, edge_contractility,
                                  area_preferential);
};

/// Differentiates progenitor cells using the NotchDelta::NotchDelta model
/** \param notch_delta  The pointer to the differentiation model
 *  \param steps        Number of iteration steps the NotchDelta model is run.
 *  \param linetension  The symmetric matrix of linetension interactions
 *                      between two cells of same or different type
 *  \param edge_contractility  The symmetric matrix of contractility interactions
 *                             between two cells of same or different type
 *  \param area_preferential    The preferential cell area of the different
 *                              cell types
 */
template <class NotchDelta>
void PCPVertex::differentiate_hair_cells_NotchDelta(
        std::shared_ptr<NotchDelta> notch_delta, int steps,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> linetension,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> edge_contractility,
        arma::Col<double>::fixed<CellType::num_cell_types> area_preferential)
{
    this->_log->debug("Differentiating progenitor cells to hair "
        "and support cells using the NotchDelta model ...");

    const auto nd_cells = notch_delta->get_cm()->cells();
    const auto cells = _am.cells();

    std::unordered_map<std::shared_ptr<Cell>,
                       std::shared_ptr<typename NotchDelta::Cell>> cell_map;
    if (cells.size() > nd_cells.size()) {
        this->_log->warn("Cells in NotchDelta: {}. Cells in Vertex: {}",
            nd_cells.size(), cells.size());
        throw std::runtime_error("Cannot link cells of NotchDelta and Vertex "
            "models. More cells in Vertex than in NotchDelta model!");
    }
    unsigned int iterator;
    cell_map.reserve(cells.size());
    for (iterator = 0; iterator < cells.size(); iterator++) {
        cell_map.insert({cells[iterator], nd_cells[iterator]});
    }
    for (void(); iterator < nd_cells.size(); iterator++) {
        nd_cells[iterator]->state.cell_type = NotchDelta::CellType::inactive;
    }

    for (const auto& c : cells) {
        auto mapped_cell = cell_map.at(c);
        mapped_cell->custom_links().neighbors.clear();
        for (auto n : _am.neighbors_of(c)) {
            mapped_cell->custom_links().neighbors.push_back(cell_map.at(n));
        }
    }

    notch_delta->prolog();
    for (int i = 0; i < steps; i++) {
        notch_delta->iterate();
    }
    notch_delta->epilog();

    for (auto pair = cell_map.begin(); pair != cell_map.end(); pair++) {
        auto type = pair->second->state.cell_type;
        if (type == NotchDelta::CellType::hair) {
            pair->first->state.type = CellType::hair;
        }
        else if (type == NotchDelta::CellType::support) {
            pair->first->state.type = CellType::support;
        }
        else {
            pair->first->state.type = CellType::progenitor;
        }
    }

    differentiate_hair_cells_hlpr(linetension, edge_contractility,
                                  area_preferential);
};

/// Increase the domain size by a certain area
/** This remaps the domain of size A to size A + dA while keeping the 
 *  relation of Lx to Ly constant.
 * 
 *  Thereby proliferation of cells can be performed in a periodic setup
 *  without changing the parameters of the system. 
 */
void PCPVertex::increase_domain_size(double area)
{
    const auto domain = this->_space->get_domain_size();
    if (-1. * area > domain[0] * domain[1]) {
        throw std::invalid_argument("Cannot decrease the domain size by "
                "an area larger than the domain size. dA = " + 
                std::to_string(area) + " and A = " +
                std::to_string(domain[0] * domain[1]));
    }
    double ratio = domain[0] / domain[1];
    double ly = std::sqrt(domain[1] * domain[1] + area / ratio);
    double lx = ratio * ly;

    this->_space->set_domain_size({lx, ly});
};

/// Stretch the domain size
/** \param stretch   The stretching distance
 *  \param compensate   Whether to compensate the growth of the dissue by 
 *                      increase of preferential area
 *  \param fix_hc_volume   Whether to fix the volume of type CellType::hair
 * 
 *  \return The total change in area
 */
double PCPVertex::stretch_domain(SpaceVec stretch, bool compensate,
        bool fix_hc_volume)
{
    this->_log->debug("stretching domain by ({}, {}). Compensate {}, "
                      "fix hair cell volume {}", stretch[0], stretch[1], 
                      compensate, fix_hc_volume);

    this->_space->set_domain_size(this->_space->get_domain_size() + stretch);
    auto domain = this->_space->get_domain_size();

    if (not compensate) {
        return stretch[0]*domain[1] + stretch[1]*domain[0];
    }

    int num_cells = _am.cells().size();
    if (fix_hc_volume) {
        for (const auto& c : _am.cells()) {
            num_cells -= (c->state.type == CellType::hair);
        }
        if (num_cells == 0) {
            throw std::runtime_error("All cells eliminated!");
        }
    }

    double dA = (stretch[0]*domain[1] + stretch[1]*domain[0]) / num_cells;

    if (not fix_hc_volume) {        
        const RuleFuncCell compensate_dA = [dA](const auto& cell) {
            cell->state.area_preferential += dA;
            return cell->state;
        };
        apply_rule<Update::sync>(compensate_dA, _am.cells());
    }
    else {
        const RuleFuncCell compensate_dA = [dA](const auto& cell) {
            auto state = cell->state;
            if (state.type != CellType::hair) {
                cell->state.area_preferential += dA;
            }
            return state;
        };
        apply_rule<Update::sync>(compensate_dA, _am.cells());
    }
    
    return stretch[0]*domain[1] + stretch[1]*domain[0];
};

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif