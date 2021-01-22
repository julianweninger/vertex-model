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
    if (edge->state.linetension == 0.) {
        return 0.;
    }

    double length;
    if (beta > 0) {
        const SpaceVec &a = _am.displace_virtual(edge->custom_links().a, beta);
        const SpaceVec &b = _am.displace_virtual(edge->custom_links().b, beta);
        length = this->_space->distance(a, b);
    }
    else {
        length = this->_space->distance(_am.position_of(edge->custom_links().a),
                                        _am.position_of(edge->custom_links().b));
    }
    
    return edge->state.linetension * length;
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
    if (edge->state.contractility == 0.) {
        return 0.;
    }

    double length;
    if (beta > 0) {
        const SpaceVec &a = _am.displace_virtual(edge->custom_links().a, beta);
        const SpaceVec &b = _am.displace_virtual(edge->custom_links().b, beta);
        length = this->_space->distance(a, b);
    }
    else {
        length = this->_space->distance(_am.position_of(edge->custom_links().a),
                                        _am.position_of(edge->custom_links().b));
    }
    
    return 0.5 * edge->state.contractility * pow(length, 2);
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
                       / cell->state.area_preferential);
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
    double shape_index = _am.shape_index_of(cell, beta);
    return (  0.5 * state.contractility
            * std::pow(shape_index - state.shape_index_preferential, 2));
};

// /// The energy associated with cell-cell polarity
// /** \f$ E = J_1 \sum_i \sigma_i^\alpha \sigma_i^\beta \f$, where
//  *  \f$ i \f$ is naming an edge, \f$ \alpha \f$ and \f$ \beta \f$ naming the
//  *  two neighbouring cells to edge \f$ i \f$. \f$ J_1 \f$ is an interaction
//  *  parameter PCPVertex::_cell_cell_polarity_interaction.
//  */
// double PCPVertex::cell_cell_polarity_energy (
//         Edge_ptr &e) const
// {
//     return _cell_cell_polarity_interaction * e->sigma_a * e->sigma_b;
// };

// /// The energy associated with cell intrinsic exclusion of polarity proteins
// /** Interaction of proteins on neighbouring edges.
//  *  For `cell_polarity_exclusion` > 0 accumulation of opposite sign polarity
//  *  proteins is disfavoured on neighbouring edges
//  *
//  *  \f$ E = - J_2 \sum_{<i, j>} \sigma_i^\alpha \sigma_j^\alpha \f$, where
//  *  edges \f$ i \f$ and \f$ j \f$ are adjacent bonds of cell \f$ \alpha \f$.
//  *  \f$ J_2 \f$ is interaction parameter PCPVertex::_cell_polarity_exclusion.
//  * 
//  *  \param a    edge a, alias \f$ i \f$
//  *  \param b    edge b, that is to be an interaction partner of a,
//  *              for instance a neighbour, alias edge \f$ j \f$
//  *  \param cell The cell \f$ \alpha \f$ within which exclusion is applied
//  */
// double PCPVertex::polarity_exclusion_energy (
//         Edge_ptr &a, Edge_ptr &b, const Cell_ptr &cell) const
// {        
//     double sigma_a = a->get_sigma(cell);
//     double sigma_b = b->get_sigma(cell);
    
//     return - _cell_polarity_exclusion * sigma_a * sigma_b;
// };

// /// The energy of polarity exclusion for the entire cell
// /** Cumulative for the pairwise interactions of edges within cell, 
//  *  see PCPVertex::polarity_exclusion energy
//  */
// double PCPVertex::cell_polarity_exclusion_energy (
//         Cell_ptr &c) const
// {
//     EdgeContainer edges;
//     double energy = 0.;
//     for (auto [e, flip] : c->edges_ordered) {
//         edges.push_back(e);
//     }
//     std::function<double(int, int)> factory = [c, edges, this](
//             int pos_a, int pos_b)
//     {
//         auto a = edges[pos_a];
//         auto b = edges[pos_b];

//         return this->polarity_exclusion_energy(a, b, c);
//     };
//     for (int i = 1; i < edges.size(); i++) {
//         energy += factory(i-1, i);
//     }
//     energy += factory(edges.size()-1, 0);
    
//     return energy;
// };

// /// The energy associated with constraint of zero net polarisation
// /** Within a cell the net polarisation is zero.
//  *  Included via lagrange multiplier I
//  * 
//  *  \f$ E = -\lambda_I^\alpha \sum_i \sigma_i^\alpha \f$
//  * 
//  *  \warning There is no argument available why 
//  *           \f$\dot{\lambda} = \gamma dE/d\lambda\f$
//  *           instead of \f$\dot{\lambda} = -\gamma dE/d\lambda\f$
//  */
// 
// double PCPVertex::lagrange_net_polarisation_energy (
//         Cell_ptr &c) const
// {
//     double cum_sigma = 0.;
//     for (auto [e, flip] : c->edges_ordered) {
//         cum_sigma += e->get_sigma(c);
//     }

//     return - c->lagrange_net_polarisation * cum_sigma;
// };

// /// Energy associated with constraint of constant protein level
// /** Within a cell the concentration of proteins is constant
//  *  Included via lagrange multiplier II     * 
//  * 
//  *  \f$ E = -\lambda_{II}^\alpha (\sum_i (\sigma_i^\alpha)^2 - c^\alpha) \f$
//  */
// double PCPVertex::lagrange_const_concentration_energy (
//         Cell_ptr &c) const
// {
//     double concentration = 0.;
//     for (auto [e, flip] : c->edges_ordered) {
//         concentration += std::pow(e->get_sigma(c), 2);
//     }

//     double lagrange = c->lagrange_const_concentration;
//     return - lagrange * (concentration - c->protein_concentration);
// };

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


// /// Getter for energy associated with cell-cell polarity
// double PCPVertex::get_energy_cell_cell_polarity(
//         const EdgeContainer& es) const
// {
//     double energy = 0.;
//     // for (auto &&e : es) {
//     //     energy += cell_cell_polarity_energy(e);
//     // }
//     return energy;
// }

// /// Getter for energy associated with polarity exclusion
// double PCPVertex::get_energy_polarity_exclusion (
//         const CellContainer& cs) const
// {
//     double energy = 0.;
//     // for (auto &&c : cs) {
//     //     energy += cell_polarity_exclusion_energy(c);
//     // }
//     return energy;
// }

// /// Getter for energy associated with lagrange multiplier I
// /** Langrange multiplier I is the constrain of zero net polarisation within
//  *  a cell.
//  */
// double PCPVertex::get_energy_lagrange_net_polarisation(
//         const CellContainer& cs) const
// {
//     double energy = 0.;
//     // for (auto &&c : cs) {
//     //     energy += lagrange_net_polarisation_energy(c);
//     // }
//     return energy;
// }

// /// Getter for energy associated with lagrange multiplier II
// /** Langrange multiplier II is the constrain of constant level of proteins
//  */
// double PCPVertex::get_energy_lagrange_const_concentration(
//         const CellContainer& cs) const
// {
//     double energy = 0.;
//     // for (auto &&c : cs) {
//     //     energy += lagrange_const_concentration_energy(c);
//     // }
//     return energy;
// }

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
   return (  get_energy_linetension(es, beta)
           + get_energy_edge_contractility(es, beta)
           + get_energy_areaelasticity(cs, beta)
           + get_energy_cell_contractility(cs, beta));
}

/// Getter for the relative energy change from previous to last step
double PCPVertex::get_rel_energy_change () const
{
    double energy = get_energy();
    
    double energy_change = energy - _energy_previous_step;
    return energy_change / (energy + 1e-14);
}

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif