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

/// Add skew to the domain's boundary condition
/** Args:
 *      - add_skew (SpaceVec): The additional skew
 *      - absolute (bool): Whether add_skew is absolute values on the length
 *                     of the domain or relative per unit length
 *      - deform_plastic (bool): Whether to move vertices according to
 *                     deformation in a plastic manner
 * 
 *  Note: Only in periodic boundary conditions as it uses the space's skew
 *      (Lees-Edwards) boundary condition
 */
PCPVertex::SpaceVec PCPVertex::skew_domain
(SpaceVec add_skew, bool absolute, bool deform_plastic) 
{
    if (not this->_space->periodic) {
        throw std::runtime_error("Cannot Skipping skew domain in non-periodic "
            "boundary conditions.");
    }

    const SpaceVec skew = _space->get_skew();
    if (not absolute) {
        add_skew = add_skew % _space->get_domain_size();
    }

    this->_log->debug("Skewing domain by ({}, {}){}",
                      add_skew[0], add_skew[1],
                      deform_plastic ? "..":" in a plastic deformation.");
    _space->set_skew(skew + add_skew);
    const SpaceVec new_skew = _space->get_skew();

    if (deform_plastic) {
        const SpaceVec domain = _space->get_domain_size();
        if (fabs(_space->get_curvature()) > 1.e-12) {
            throw std::runtime_error("not implemented skew.");
            double curvature = _space->get_curvature();
            for (const auto& vertex : _am.vertices()) {
                SpaceVec pos = _am.position_of(vertex);
                auto [rho, theta] = _space->transform_radial(pos);
                double skew_theta = add_skew[0] * curvature;
                double max_theta = (domain[0] / 2.) * curvature;

                double r = (rho - 1. / curvature) / domain[1] + 0.5;
                double r_theta = theta / (2 * max_theta);

                theta += skew_theta * r;
                rho += add_skew[1] * (r_theta + max_theta);

                pos = _space->transform_cartesian(rho, theta);
                _am.move_to(vertex, pos);
            }
        }
        else{
            for (const auto& vertex : _am.vertices()) {
                SpaceVec pos = _am.position_of(vertex);
                SpaceVec box_pos = pos - skew % arma::shift(pos / domain, -1);

                SpaceVec new_pos = box_pos + new_skew % arma::shift(box_pos / domain, -1);
                
                _am.move_to(vertex, new_pos);
            }
        }
    }
    
    return new_skew;
}

/// Use boundary conditions to imply a curvature to system
/** In periodic BC, use curved bc
 *  In non-periodic BC, use fixed boundary vertices to imply curvature
 * 
 *  \param add_curvature        The additional curvature
 *  \param deform_plastic       Whether to update vertices positions in a
 *                              plastic deformation update
 */
void PCPVertex::curve_boundary
(double add_curvature, bool deform_plastic)
{
    double curvature;
    if (this->_space->periodic) {
        curvature = _space->get_curvature() + add_curvature;
        _space->set_curvature(curvature);
    }
    else {
        this->fix_boundary(true);
        curvature = _boundary_param.get_curvature() + add_curvature;
        _boundary_param.set_curvature(curvature);
    }

    if (not this->_space->periodic or deform_plastic) {

        // previous curvature
        double k0 = std::max(curvature - add_curvature, 1.e-8);
        double R0 = 1. / k0;

        // new curvature
        double k1 = std::max(curvature, 1.e-8);
        double R1 = 1. / k1;
        
        double max_theta_0, max_theta_1;
        SpaceVec origin_0, origin_1;
        if (_space->periodic) {
            const SpaceVec domain = _space->get_domain_size();

            max_theta_0 = (domain[0] / 2.) * k0;
            max_theta_1 = (domain[0] / 2.) * k1;

            origin_0 = SpaceVec({domain[0] / 2., domain[1] / 2. - R0});
            origin_1 = SpaceVec({domain[0] / 2., domain[1] / 2. - R1});
        }
        else {
            origin_0 = SpaceVec({0., -R0});
            origin_1 = SpaceVec({0., -R1});

            double max_t = std::numeric_limits<double>::min();
            double min_t = std::numeric_limits<double>::max();

            for (const auto& v : _am.vertices()) {
                SpaceVec pos = _am.position_of(v);
                SpaceVec displ = pos - origin_0;

                double theta = std::atan2(displ[0], displ[1]);
                max_t = std::max(max_t, theta);
                min_t = std::min(min_t, theta);
            }

            max_theta_0 = max_t - min_t;
            
            // conserve length of the domain
            max_theta_1 = R0 / R1 * max_theta_0;
        }

        for (const auto& v : _am.vertices()) {
            if (not deform_plastic and not _am.is_boundary(v)) {
                continue;
            }

            SpaceVec pos = _am.position_of(v);
            SpaceVec displ = pos - origin_0;

            // relative radial position
            double rho = arma::norm(displ) - R0; 
            // relative angular position
            double theta = std::atan2(displ[0], displ[1]) / max_theta_0;

            SpaceVec new_pos = (
                origin_1 + (rho + R1) * SpaceVec({sin(theta * max_theta_1), 
                                                  cos(theta * max_theta_1)})
            );

            _am.move_to(v, new_pos);
        }
    }
}


} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif