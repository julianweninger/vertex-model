#ifndef UTOPIA_MODELS_PCPVERTEX_ENERGY_HH
#define UTOPIA_MODELS_PCPVERTEX_ENERGY_HH

namespace Utopia {
namespace Models {
namespace PCPVertex {


/// The energy associated with linetension per edge
/** \f$ E_{ij} = \lambda_{ij} l_{ij} \f$ for the edge connecting vertices
 *  \f$ i \f$ and \f$ j \f$ with length \f$ l_{ij} \f$.
 * 
 *  \param edge     The object
 *  \param beta     If > 0. predicts the energy at this step size along
 *                  direction of update. I.e. vertices are virtually moved in
 *                  direction of update and on that configuration the energy is
 *                  calculated.
 *                  For beta = 0, returns the energy as in current vertex
 *                  position.
 */
double PCPVertex::line_tension_energy (
        const std::shared_ptr<Edge>& edge, double beta) const
{
    if (fabs(edge->state.linetension()) < 1.e-12) {
        return 0.;
    }

    double length;
    if (beta > 1.e-14) {
        const SpaceVec &a = _am.displace_virtual(edge->custom_links().a, beta);
        const SpaceVec &b = _am.displace_virtual(edge->custom_links().b, beta);
        length = this->_space->distance(a, b);
    }
    else {
        length = this->_space->distance(
            _am.position_of(edge->custom_links().a),
            _am.position_of(edge->custom_links().b));
    }
    
    return edge->state.linetension() * length;
};

/// The energy associated with contractility per edge
/** \f$ E_{ij} = \frac{1}{2} \Gamma{ij} l_{ij}^2 \f$ for the edge connecting
 *  vertices \f$ i \f$ and \f$ j \f$ with length \f$ l_{ij} \f$.
 * 
 *  \param edge     The object
 *  \param beta     If > 0. predicts the energy at this step size along
 *                  direction of update. I.e. vertices are virtually moved in
 *                  direction of update and on that configuration the energy is
 *                  calculated.
 *                  For beta = 0, returns the energy as in current vertex
 *                  position.
 */
double PCPVertex::edge_contractility_energy (
        const std::shared_ptr<Edge>& edge, double beta) const
{
    if (fabs(edge->state.contractility()) < 1.e-12) {
        return 0.;
    }

    double length;
    if (beta > 1.e-14) {
        const SpaceVec &a = _am.displace_virtual(edge->custom_links().a, beta);
        const SpaceVec &b = _am.displace_virtual(edge->custom_links().b, beta);
        length = this->_space->distance(a, b);
    }
    else {
        length = this->_space->distance(_am.position_of(edge->custom_links().a),
                                        _am.position_of(edge->custom_links().b));
    }
    
    return 0.5 * edge->state.contractility() * pow(length, 2);
};

/// The energy associated with area elasticity per cell
/** \f$ E_i = K/2 * (A_i - A0_i)**2 \f$ for cell \f$ i \f$ with area \f$ A_i \f$
 *  and preferential area \f$ A0_i \f$.
 * 
 *  \param edge     The object
 *  \param beta     If > 0. predicts the energy at this step size along
 *                  direction of update. I.e. vertices are virtually moved in
 *                  direction of update and on that configuration the energy is
 *                  calculated.
 *                  For beta = 0, returns the energy as in current vertex
 *                  position.
 */
double PCPVertex::area_elasticity_energy (
        const std::shared_ptr<Cell>& cell, double beta) const
{
    // for beta = 0, returns same as area_of(cell)
    double rel_area = (  _am.area_of(cell, beta)
                       / cell->state.area_preferential());
    return 0.5 * _area_elasticity * pow(rel_area - 1., 2);
};

/// The energy associated with contractility per cell
/** \f$ E_i = \Gamma/2 * (p_i - p0_i)**2 \f$ for cell \f$ i \f$ with
 *  shape index \f$ p_i \f$ and preferential shape index \f$ p0_i \f$.
 *  Shape index \f$ p = P / \sqrt(A) \f$ with the cell's perimeter \f$ P \f$,
 *  analogous for the preferential shape index.
 * 
 *  \param edge     The object
 *  \param beta     If > 0. predicts the energy at this step size along
 *                  direction of update. I.e. vertices are virtually moved in
 *                  direction of update and on that configuration the energy is
 *                  calculated.
 *                  For beta = 0, returns the energy as in current vertex
 *                  position.
 */
double PCPVertex::cell_contractility_energy (
        const std::shared_ptr<Cell>& cell, double beta) const
{
    const auto state = cell->state;

    if (fabs(state.contractility) < 1.e-12)
    {
        return 0.;
    }

    if (this->_cell_contractility_impl == Contractile) {
        double perimeter = _am.perimeter_of(cell, beta);
        return (  0.5 * state.contractility
                * std::pow(  perimeter / sqrt(state.area_preferential())
                           - state.shape_index_preferential, 2));
    }
    else if (this->_cell_contractility_impl == Shape_elastic) {
        double shape_index = _am.shape_index_of(cell, beta);
        return (  0.5 * state.contractility
                * std::pow(shape_index - state.shape_index_preferential, 2));
    }
    else {
        throw std::runtime_error(fmt::format(
            "Not Implemented Error: "
            "Cell Contractility Implementation {}!",
            this->_cell_contractility_impl));
    }
};

/// Getter for energy associated with linetension
/** Sums PCPVertex::line_tension_energy for all entities
 * 
 *  \param beta     If > 0. predicts the energy at this step size along
 *                  direction of update. I.e. vertices are virtually moved in
 *                  direction of update and on that configuration the energy is
 *                  calculated.
 *                  For beta = 0, returns the energy as in current vertex
 *                  position.
 */
double PCPVertex::get_energy_linetension(
        const AgentContainer<Edge>& es, double beta) const
{
    double energy = 0.;
    for (const auto& e : es) {
        energy += line_tension_energy(e, beta);
    }
    return energy;
}

/// Getter for energy associated with contractility of junctions
/** Sums PCPVertex::edge_contractility_energy for all entities
 * 
 *  \param beta     If > 0. predicts the energy at this step size along
 *                  direction of update. I.e. vertices are virtually moved in
 *                  direction of update and on that configuration the energy is
 *                  calculated.
 *                  For beta = 0, returns the energy as in current vertex
 *                  position.
 */
double PCPVertex::get_energy_edge_contractility(
        const AgentContainer<Edge>& es, double beta) const
{
    double energy = 0.;
    for (const auto& e : es) {
        energy += edge_contractility_energy(e, beta);
    }
    return energy;
}

/// Getter for energy associated with area elasticity
/** Sums PCPVertex::area_elasticity_energy for all entities
 * 
 *  \param beta     If > 0. predicts the energy at this step size along
 *                  direction of update. I.e. vertices are virtually moved in
 *                  direction of update and on that configuration the energy is
 *                  calculated.
 *                  For beta = 0, returns the energy as in current vertex
 *                  position.
 */
double PCPVertex::get_energy_areaelasticity(
        const AgentContainer<Cell>& cs, double beta) const
{
    double energy = 0.;
    for (const auto& c : cs) {
        energy += area_elasticity_energy(c, beta);
    }
    return energy;
}

/// Getter for energy associated with contractility of cells
/** Sums PCPVertex::cell_contractility_energy for all entities
 * 
 *  \param beta     If > 0. predicts the energy at this step size along
 *                  direction of update. I.e. vertices are virtually moved in
 *                  direction of update and on that configuration the energy is
 *                  calculated.
 *                  For beta = 0, returns the energy as in current vertex
 *                  position.
 */
double PCPVertex::get_energy_cell_contractility(
        const AgentContainer<Cell>& cs, double beta) const
{
    double energy = 0.;
    for (const auto& c : cs) {
        energy += cell_contractility_energy(c, beta);
    }
    return energy;
}

/// Getter for the boundary area energy 
/** Similar to a cell, the domain's area has quadratic energy contribution
 */ 
double PCPVertex::get_boundary_area_energy (double beta) const {
    if (   _space->periodic
        or fabs(_boundary_param.area_elasticity) < 1.e-12)
    {
        return 0.;
    }

    const auto boundary = _am.get_boundary_edges();
    const auto N = _am.cells().size();

    double area = _am.area_of(boundary, beta);

    double rel_area = area / this->_boundary_param.area_preferential(N);
    double area_elasticity = _boundary_param.area_elasticity;

    return 0.5 * area_elasticity * std::pow(rel_area - 1., 2);
}

/// Getter for the boundary shape energy
/** Similar to a cell, the boundaries shape has quadratic energy contribution
 */ 
double PCPVertex::get_boundary_shape_energy(double beta) const {
    if (   _space->periodic
        or fabs(_boundary_param.contractility) < 1.e-12)
    {
        return 0.;
    }
    
    const auto boundary = _am.get_boundary_edges();
    double contractility = _boundary_param.contractility;
    double shape_index;
    double shape_index_pref = _boundary_param.shape_index_preferential;

    if (this->_cell_contractility_impl == Contractile) {
        const auto N = _am.cells().size();
        double perimeter = _am.perimeter_of(boundary, beta);
        shape_index = perimeter / sqrt(_boundary_param.area_preferential(N));
    }
    else if (this->_cell_contractility_impl == Shape_elastic) {
        shape_index = _am.shape_index_of(boundary, beta);
    }
    else {
        throw std::runtime_error(fmt::format(
            "Not Implemented Error: "
            "Cell Contractility Implementation {}!",
            this->_cell_contractility_impl));
    }


    return (  0.5 * contractility
            * std::pow(shape_index - shape_index_pref, 2.));
}

/// Getter for the boundary stripe energy 
/** A quadratic potential for boundary vertices that are outside a stripe
 */ 
double PCPVertex::get_boundary_stripe_energy (double beta) const {
    if (   _space->periodic
        or fabs(_boundary_param.stripe_potential_constant) < 1.e-12)
    {
        return 0.;
    }

    double curvature = std::max(_boundary_param.stripe_curvature, 1.e-10);
    double R = 1. / curvature;

    // the (fixed) center of the circle stripe
    SpaceVec origin = *_boundary_param.stripe_origin - SpaceVec({0., R});

    double inner_radius = (R - _boundary_param.stripe_width / 2.);
    double outer_radius = (R + _boundary_param.stripe_width / 2.);
    
    double energy = 0.;
    // apply to all vertices outside the domain
    for (const auto &v : _am.vertices()) {
        // the position wrt origin
        SpaceVec pos = _am.displace_virtual(v, beta) - origin;

        double radius = arma::norm(pos);
        if (radius < inner_radius) {
            energy += std::pow(radius - inner_radius, 2.);
        }
        else if (radius > outer_radius) {
            energy += std::pow(radius - outer_radius, 2.);
        }
    }

    return 0.5 * _boundary_param.stripe_potential_constant * energy;
}

/// Getter for energy
/** Sums the following energies for provided entities
 *      -# PCPVertex::get_energy_linetension
 *      -# PCPVertex::get_energy_edge_contractility
 *      -# PCPVertex::get_energy_areaelasticity
 *      -# PCPVertex::get_energy_cell_contractility
 * 
 *  \param beta     If > 0. predicts the energy at this step size along
 *                  direction of update. I.e. vertices are virtually moved in
 *                  direction of update and on that configuration the energy is
 *                  calculated.
 *                  For beta = 0, returns the energy as in current vertex
 *                  position.
 */
double PCPVertex::get_energy (
        const AgentContainer<Edge>& es,
        const AgentContainer<Cell>& cs,
        double beta) const
{
    // WARN It is supposed to work on the AgentContainers es and cs.
    //      If functions work on the globally defined edges() and cells(),
    //      handle with care and check the entities_manager transitions.

    if (not std::isfinite(beta)) {
        throw std::runtime_error(fmt::format("Cannot calculate energy "
            "with non finite beta={}", beta));
    }

    return (  get_energy_linetension(es, beta)
            + get_energy_edge_contractility(es, beta)
            + get_energy_areaelasticity(cs, beta)
            + get_energy_cell_contractility(cs, beta)
            + get_boundary_energy(beta));
}

/// Getter for the relative energy change from previous to last step
double PCPVertex::get_rel_energy_change () const
{
    double energy = get_energy();
    
    double energy_change = energy - _energy_previous_step;
    return energy_change / std::max(fabs(energy), 1e-14);
}

/// Getter for the relative energy change from previous to last step
double PCPVertex::get_rel_energy_change (double energy, double energy_0) const
{
    return (energy - energy_0) / std::max(fabs(energy), 1e-14);
}

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif