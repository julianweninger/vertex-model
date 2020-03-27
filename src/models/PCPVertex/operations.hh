#ifndef UTOPIA_MODELS_SAVANNAHETEROGENEOUS_OPERATIONS_HH
#define UTOPIA_MODELS_SAVANNAHETEROGENEOUS_OPERATIONS_HH

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
template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc, polarity_proteins>::jiggle_vertices(
        double intensity)
{
    this->_log->debug("Jiggling the vertices on a length scale of "
                        "{} ..", intensity);
    for (auto v : _vertices) {
        v->x += 2*intensity * _prob_distr(*this->_rng) - intensity;
        v->y += 2*intensity * _prob_distr(*this->_rng) - intensity;
        correct_periodic_bc<periodic_bc>(v);
    }

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
template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc, polarity_proteins>::differentiate_hair_cells_hlpr(
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
    
    // Set the cell preferential area
    for (auto &c : _cells) {
        c->area_preferential = area_preferential(c->type);
    }

    // Set the surface tension
    for (auto &e : _edges) {
        Cell::CellType cell_a_type, cell_b_type;
        if (e->adj_cell_a.expired()) {
            cell_a_type = Cell::CellType::support; }
        else { cell_a_type = e->adj_cell_a.lock()->type; }
        if (e->adj_cell_b.expired()) {
            cell_b_type = Cell::CellType::support; }
        else { cell_b_type = e->adj_cell_b.lock()->type; }

        e->linetension = linetension(cell_a_type, cell_b_type);
        e->contractility = edge_contractility(cell_a_type, cell_b_type);
    }

    _linetension = linetension;
    _edge_contractility = edge_contractility;

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
template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc, polarity_proteins>::differentiate_hair_cells_random(
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
    
    // Set the cell type
    for (auto &c : _cells) {
        if (_prob_distr(*this->_rng) < fraction) { 
            c->type = CellType::hair; }
        else { 
            c->type = CellType::support; }
    }

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
template <bool periodic_bc, bool polarity_proteins>
template <class NotchDelta>
void PCPVertex<periodic_bc,
               polarity_proteins>::differentiate_hair_cells_NotchDelta(
        std::shared_ptr<NotchDelta> notch_delta, int steps,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> linetension,
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> edge_contractility,
        arma::Col<double>::fixed<CellType::num_cell_types> area_preferential)
{
    this->_log->debug("Differentiating progenitor cells to hair "
        "and support cells using the NotchDelta model ...");

    auto nd_cells = notch_delta->get_cm()->cells();

    std::unordered_map<Cell_ptr, typeof(nd_cells.back())> cell_map;
    if (_cells.size() > nd_cells.size()) {
        this->_log->warn("Cells in NotchDelta: {}. Cells in Vertex: {}",
            nd_cells.size(), _cells.size());
        throw std::runtime_error("Cannot link cells of NotchDelta and Vertex "
            "models. More cells in Vertex than in NotchDelta model!");
    }
    for (int i = 0; i < _cells.size(); i++) {
        cell_map.insert(std::make_pair(_cells[i], nd_cells[i]));
    }
    for (auto c : _cells) {
        auto mapped_cell = cell_map.at(c);
        mapped_cell->custom_links().neighbors.clear();
        for (auto n : neighbors_of(c)) {
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
            pair->first->type =  Cell::CellType::hair;
        }
        else if (type == NotchDelta::CellType::support) {
            pair->first->type =  Cell::CellType::support;
        }
        else {
            pair->first->type =  Cell::CellType::progenitor;
        }
    }

    differentiate_hair_cells_hlpr(linetension, edge_contractility,
                                  area_preferential);
};

/// Perform a cell division on specific cell
/** Divides a specific cell into two identical cells with properties derived
 *  from the common parent cell. 
 *  The division is performed at a given angle through the parent cell's
 *  center. This defines the axis of division that will form a new edge 
 *  between the two new cells.
 * 
 *  \param cell     pointer to the cell that is to be divided
 *  \param division_angle   angle (in rad) at which the cell
 */
template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc, polarity_proteins>::divide_cell(
        Cell_ptr cell, double division_angle)
{
    auto cell_it = std::find(_cells.begin(), _cells.end(), cell);

    if (cell_it == _cells.end()) {
        throw std::invalid_argument("Cannot divide cell at position "
            "({}, {}), because its not a member of cells in vertex model!");
    }

    this->divide_cell(cell_it, division_angle);
};

/// Increase the domain size by a certain area
/** This remaps the domain of size A to size A + dA while keeping the 
 *  relation of Lx to Ly constant.
 * 
 *  Thereby proliferation of cells can be performed in a periodic setup
 *  without changing the parameters of the system. 
 */
template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc, polarity_proteins>::increase_domain_size(
        double area)
{
    if (-1. * area > _Lx * _Ly) {
        throw std::invalid_argument("Cannot decrease the domain size by "
                "an area larger than the domain size. dA = " + 
                std::to_string(area) + " and A = " +
                std::to_string(_Lx * _Ly));
    }
    double ratio = _Lx / double(_Ly);
    _Ly = std::sqrt(_Ly*_Ly + area / ratio);
    _Lx = ratio * _Ly;
};

/// Stretch the domain size
/** \param dx   The stretching distance in x
 *  \param dy   The stretching distance in y
 *  \param compensate   Whether to compensate the growth of the dissue by 
 *                      increase of preferential area
 *  \param fix_hc_volume   Whether to fix the volume of type CellType::hair
 * 
 *  \return The total change in area
 */
template <bool periodic_bc, bool polarity_proteins>
double PCPVertex<periodic_bc, polarity_proteins>::stretch_domain(
        double dx, double dy, bool compensate, bool fix_hc_volume)
{
    this->_log->debug("stretching domain by ({}, {}). Compensate {}, "
                        "fix hair cell volume {}", dx, dy, compensate, 
                        fix_hc_volume);
    _Lx += dx;
    _Ly += dy;

    if (not compensate) {
        return dx * _Ly + dy * _Lx;
    }

    int num_cells = _cells.size();
    if (fix_hc_volume) {
        for (auto c : _cells) {
            num_cells -= c->type == CellType::hair;
        }
        if (num_cells == 0) {
            throw std::runtime_error("All support cells eliminated!");
        }
    }

    double dA = (dx * _Ly + dy * _Lx) / num_cells;

    std::function<void(Cell_ptr&)> compensate_dA = [dA](Cell_ptr& cell) {
        cell->area_preferential += dA;
        return;
    };
    std::function<void(Cell_ptr&)> compensate_dA_non_hc = [dA](
            Cell_ptr& cell)
    {
        if (cell->type != CellType::hair) {
            cell->area_preferential += dA;
        }
        return;
    };

    if (not fix_hc_volume) {
        std::for_each(_cells.begin(), _cells.end(), compensate_dA);
    }
    else {
        std::for_each(_cells.begin(), _cells.end(), compensate_dA_non_hc);
    }

    return dx * _Ly + dy * _Lx;
};

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif