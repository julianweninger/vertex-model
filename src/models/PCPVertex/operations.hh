#ifndef UTOPIA_MODELS_PCPVERTEX_OPERATIONS_HH
#define UTOPIA_MODELS_PCPVERTEX_OPERATIONS_HH

#include <typeinfo>

namespace Utopia {
namespace Models {
namespace PCPVertex {

/// Apply a perturbation to the position of vertices
/** \details Move every vertex randomly on a given lengthscale.
 *           The length scale is \f$ l = I * \sqrt{A_{domain} / num_cells} \f$,
 *           with the intensity I.
 *           The displacement is pulled from a uniform distribution form [-l, l]
 *           for all coordinates.
 */
void PCPVertex::jiggle_vertices(double intensity)
{
    auto num_cells = this->_am.cells().size();
    const SpaceVec domain = this->_space->get_domain_size();
    double domain_area;
    if (this->_space->periodic) {
        domain_area = domain[0]*domain[1];
    }
    else {
        domain_area = 0.;
        for (const auto& c : this->_am.cells()) {
            domain_area += c->state.area_preferential();
        }
    }
    intensity *= sqrt(domain_area / num_cells);

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

    this->init_minimization();
};

/// Increase the domain size by a certain area
/** This remaps the domain of size A to size A + dA while keeping the 
 *  relation of Lx to Ly constant.
 * 
 *  Thereby proliferation of cells can be performed in a periodic setup
 *  without changing the parameters of the system. 
 * 
 *  If deform_plastic, all distances are stretched correspondingly.
 *  I.e. relative positions in space are maintained
 */
void PCPVertex::increase_domain_size(double area,
                                     bool deform_plastic = true)
{
    if (not this->_space->periodic) {
        return;
    }

    const SpaceVec domain = this->_space->get_domain_size();
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

    const SpaceVec new_domain = this->_space->get_domain_size();
    if (deform_plastic) {
        if (this->_space->get_curvature() > 1.e-12) {
            throw std::runtime_error("Increase domain with curvature not "
                "implemented");
        }
        for (const auto& vertex : _am.vertices()) {
            _am.move_to(vertex, _am.position_of(vertex) / domain % new_domain);
        }
    }
};

/// Stretch the domain size
/** \param stretch   The stretching distance
 *  \param compensate   Whether to compensate the growth of the dissue by 
 *                      increase of preferential area
 *  \param fix_hc_area   Whether to fix the area of type CellType::hair
 *  \param fix_sc_area   Whether to fix the area of type CellType::support
 * 
 *  \return The total change in area
 */
double PCPVertex::stretch_domain(SpaceVec stretch, bool compensate,
        bool fix_hc_area, bool fix_sc_area, bool deform_plastic = true)
{
    if (not this->_space->periodic) {
        if (compensate) {
            this->_log->warn("Skipping stretch domain in non-periodic boundary "
                "conditions. The domain is large enough to fit a tissue of any "
                "size.. However, `compensate` is requested, which will also "
                "be skipped!");
        }
        return 0.;
    }

    if (compensate and fix_hc_area and fix_sc_area) {
        throw std::invalid_argument("Cannot compensate area while fixing "
            "hair and support cell area in operation 'increment_domain'!");
    }

    this->_log->debug("Stretching domain by ({}, {}). {}compensating area, "
                      "{}fixing hair cell area, {}fixing support cell area ...",
                      stretch[0], stretch[1],
                      compensate ? "":"not ",
                      fix_hc_area ? "":"not ",
                      fix_sc_area ? "":"not ");
    
    const SpaceVec domain = this->_space->get_domain_size();
    this->_space->set_domain_size(domain + stretch);

    const SpaceVec new_domain = this->_space->get_domain_size();
    const double area_change = (  stretch[0] * new_domain[1]
                                + stretch[1] * new_domain[0]);


    if (deform_plastic) {
        if (this->_space->get_curvature() > 1.e-12) {
            throw std::runtime_error("Increase domain with curvature not "
                "implemented");
        }
        for (const auto& vertex : _am.vertices()) {
            SpaceVec pos = _am.position_of(vertex);
            _am.move_to(vertex, _am.position_of(vertex) / domain % new_domain);
        }
    }


    if (not compensate) {
        return area_change;
    }

    const auto& cells = _am.cells();
    if (fix_hc_area) {
        const auto num_cells = std::count_if(
            cells.begin(), cells.end(), 
            [](const auto& c) {
                return c->state.type != CellType::hair;
            });
        double dA = area_change / num_cells;
        const RuleFuncCell compensate_dA = [dA](const auto& cell) {
            auto state = cell->state;
            if (state.type != CellType::hair) {
                state._area_preferential += dA;
            }
            return state;
        };
        apply_rule<Update::sync>(compensate_dA, cells);    
    }
    else if (fix_sc_area) {
        const auto num_cells = std::count_if(
            cells.begin(), cells.end(), 
            [](const auto& c) {
                return c->state.type != CellType::support;
            });
        double dA = area_change / num_cells;
        const RuleFuncCell compensate_dA = [dA](const auto& cell) {
            auto state = cell->state;
            if (state.type != CellType::support) {
                state._area_preferential += dA;
            }
            return state;
        };
        apply_rule<Update::sync>(compensate_dA, cells);
    }
    else {
        const auto num_cells = cells.size();
        double dA = area_change / num_cells;
        const RuleFuncCell compensate_dA = [dA](const auto& cell) {
            cell->state._area_preferential += dA;
            return cell->state;
        };
        apply_rule<Update::sync>(compensate_dA, cells);
    }
    
    return area_change;
};

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif