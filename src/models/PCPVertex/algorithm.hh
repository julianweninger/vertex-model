#ifndef UTOPIA_MODELS_SAVANNAHETEROGENEOUS_ALGORITHM_HH
#define UTOPIA_MODELS_SAVANNAHETEROGENEOUS_ALGORITHM_HH

namespace Utopia {
namespace Models {
namespace PCPVertex {

/// Determine step size to reach the line minimum along direction of update
/** See www.acclab.helsinki.fi/~aakurone/atomistiset/lecturenotes/lecture12_2up.pdf
 *  for details of algorithm
 *  
 *  \param dt           The step size of previous step
 *  \param energy_0     The energy at current position
 *  \return {step size to reach line minimum, energy at line minimum}
 */
template <bool periodic_bc>
std::pair<double, double> PCPVertex<periodic_bc>::determine_timestep (
        double dt, const double energy_0) const
{
    // Bracket the minimum
    dt = std::min(2*dt, 2.);

    // check that actually moving towards a minimum
    if (this->get_energy(1e-6) >= energy_0) {
        this->_log->debug("Already in minimum. At step size of 1e-8 the energy "
            "along direction of update increased.");
        
        return std::make_pair(0., energy_0);
    }
    
    // check that direction of update actually is non-zero
    double f_squared = 0.;
    const auto& vertices = _am.vertices();
    for (const auto& v : vertices) {
        f_squared += arma::norm(v->state.f);
    }
    if (f_squared/vertices.size() < _minimisation_precision/100.) {
        this->_log->debug("Linear extrapolation: {}", f_squared);
        for (double dt = 1e-11; dt < 100.; dt *= 2) {
            this->_log->debug("DEBUG At step of {} along direction of "
                "update the energy is changing by {}", dt, 
                (this->get_energy(dt) - energy_0)/dt);
        }
        return std::make_pair(0., energy_0);
    }
    // TODO remove this part if not causing abortions 

    double pos_1 = 0.; // left boundary
    bool calculate_energy_min = true;
    double energy_1 = energy_0;
    bool calculate_energy_2 = true;
    double pos_2 = dt; // right boundary
    double energy_2 = this->get_energy(pos_2);
    double pos_min = dt/2.;
    double energy_min = this->get_energy(pos_min);
    while (true) {
        if (dt/2 > 100.) {
            if (fabs(energy_2 - energy_1) < _minimisation_precision) {
                // the energy function is flat
                return std::make_pair(0., energy_0);
            }
            this->_log->warn("At step of {} along direction of "
                "update the energy is still decreasing by {}", dt/2, 
                energy_2 - energy_1);
            for (double dt = 1e-11; dt < 100.; dt *= 2) {
                this->_log->warn("DEBUG At step of {} along direction of "
                    "update the energy is changing by {}", dt, 
                    (this->get_energy(dt) - energy_0)/dt);
            }
            // NOTE only if energy_2 < energy_1: dt -> 2*dt
            throw std::runtime_error("Unable to find bracket to "
                "local energy minimum!");
        }
        if (dt < 1e-10) {
            // the energy is increasing, hence already at local minimum
            // NOTE only if energy_min > energy_1: dt -> dt / 2
            return std::make_pair(0., energy_0);
        }

        if (energy_2 < energy_1) {
            dt *= 2.;
            energy_min = energy_2;
            energy_2 = this->get_energy(dt);
            continue; // dt was not a right boundary to minimum
        }

        if (energy_min > energy_1) {
            dt = dt/2;
            energy_2 = energy_min;
            energy_min = this->get_energy(dt/2);
            continue; // we know that close to pos_1 there is a minimum
        }

        // Both conditions fulfilled!
        // there is a minimum within interval [0, pos_2]
        pos_2 = dt;
        pos_min = dt / 2.;
        break;
    }

    if (fabs(energy_min - energy_0) < _minimisation_precision) {
        // the found minimum fulfills our condition 
        return std::make_pair(pos_min, energy_min);
    }

    // Find the minimum between brackets with parabolic interpolation
    while(true) {
        double term_3 = (pos_2 - pos_1)*(energy_2 - energy_min);
        double term_4 = (pos_2 - pos_min)*(energy_2 - energy_1);
        double term_1 = (pos_2 - pos_1)*term_3;
        double term_2 = (pos_2 - pos_min)*term_4;
        
        // the minimum of parabola through 1, 2, min
        double pos_4 = pos_2 - 0.5*(term_1 - term_2)/(term_3 - term_4);
        if (pos_4 < pos_1 or pos_4 > pos_2) {
            this->_log->warn("pos_1={}, pos_min={}, pos_2={}",
                                pos_1, pos_min, pos_2);
            this->_log->warn("D_energy_1={}, D_energy_min={}, "
                                "D_energy_2={}. wrt energy_0",
                                energy_0-energy_1, energy_0-energy_min,
                                energy_0-energy_2);
            this->_log->warn("Fitted minimum: x={}", pos_4);
            throw std::runtime_error("Energy minimisation failed! "
                                        "Parabola fit outside brackets");
        }

        double energy_4 = this->get_energy(pos_4);
        if (fabs(energy_4 - energy_min) < _minimisation_precision) {
            dt = pos_4;
            break; // found the minimum with required precision
        }
        
        if (energy_4 - energy_0 > 0.) {
            // there is a closer maximum. Retry with smaller step size
            return determine_timestep(pos_4 / 2., energy_0);
        }

        // ** choose new brackets
        // pos_4 is to the left of minimum and is smaller
        // pos_4 becomes new minimum between 1 and 2->min
        if (pos_4 - pos_min < 0. and energy_4 - energy_min < 0.) {
            pos_2 = pos_min; energy_2 = energy_min;
            pos_min = pos_4; energy_min = energy_4;
        }
        // pos_4 is to the left of minimum but is larger
        // pos_4 becomes new pos 1
        else if (pos_4 - pos_min < 0.) {
            pos_1 = pos_4; energy_1 = energy_4;
        }
        // pos_4 is to the right of minimum and is smaller
        // pos_4 becomes new minimum between 1->min and 2
        else if (energy_4 - energy_min < 0.) {
            pos_1 = pos_min; energy_1 = energy_min;
            pos_min = pos_4; energy_min = energy_4;
        }
        // pos_4 is to the right of minimum but is larger
        // pos_4 becomes new pos 2
        else {
            pos_2 = pos_4; energy_2 = energy_4;
        }
    }

    return std::make_pair(dt, energy_min);
};

/// Initialisation of the energy minisation process
template <bool periodic_bc>
void PCPVertex<periodic_bc>::init_minimisation ()
{
    if (this->_update_scheme != ConjugateGradient) {
        return;
    }

    RuleFuncVertex init = [this](const auto& vertex) {
        auto state = vertex->state;
        state.g = this->_am.position_of(vertex);
        state.h = state.g;
        return state;
    };

    set_gradient();

    apply_rule<Update::sync>(init, _am.vertices());

    _dt = 1e-3;
};

/// Single step in direction of steepest gradient
/** \param adative_step     If true, a line minimisation along steepest gradient
 *                          is performed. Else, update with fixed stepsize.
 *  \return the energy after upate
 */
template <bool periodic_bc>
double PCPVertex<periodic_bc>::steepest_gradient_step (
        bool adaptive_step)
{
    set_gradient();

    double new_energy;
    if (adaptive_step) {
        std::tie(_dt, new_energy) = determine_timestep(_dt, _energy);
        double energy_change = (new_energy - _energy) / new_energy;

        if (energy_change < -1e-14) {
            this->_log->debug("Updating with timestep {} at energy change {}",
                              _dt, energy_change);
        }
        else {
            this->_log->debug("NOT updating with step size {} along direction "
                              "of update at energy change {}", _dt,
                              energy_change);
            return this->_energy;
        }
    }
    else {
        new_energy = this->get_energy(_dt);
    }

    apply_rule<Update::sync>(update_position, _am.vertices());

    return new_energy;
};

/// Single step in direction of conjugate gradient
/** \return the energy after upate
 */
template <bool periodic_bc>
double PCPVertex<periodic_bc>::conjugate_gradient_step ()
{
    const double prev_dt = _dt;

    // line minimisation along direction of update h
    double new_energy;
    std::tie(_dt, new_energy) = determine_timestep(_dt, _energy);
    double energy_change = (new_energy - _energy) / new_energy;

    if (energy_change < -1e-14) {
        this->_log->debug("Updating with timestep {} at energy change {}",
                            _dt, energy_change);
        apply_rule<Update::sync>(update_position, _am.vertices());
    }
    else {
        this->_log->debug("NOT updating with step size {} along direction "
                          "of update at energy change {}", _dt,
                          energy_change);
        return this->_energy;
    }

    // determine the conjugate gradient direction
    set_gradient();
    double gamma = 0.;
    double g_square = 0.;
    for (auto&& v : _am.vertices()) {
        gamma += arma::norm(v->state.f);
        g_square += arma::norm(v->state.g);

        v->state.g = v->state.f;
    }
    gamma /= g_square;
    for (auto&& v : _am.vertices()) {
        v->state.h = v->state.g + gamma * v->state.h;
        v->state.f = v->state.h;
    }

    return new_energy;
};

/// Select the chosen update scheme
/** \return the energy after upate
 */
template <bool periodic_bc>
double PCPVertex<periodic_bc>::perform_update_step(
        UpdateScheme update_scheme)
{
    if (update_scheme == SteepestGradient) {
        return steepest_gradient_step(false);
    }
    else if (update_scheme == SteepestGradientAdaptive) {
        return steepest_gradient_step(true);
    }
    else if (update_scheme == ConjugateGradient) {
        return conjugate_gradient_step();
    }
    else {
        throw std::logic_error("Chosen update scheme not implemented!");
    }
};

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif