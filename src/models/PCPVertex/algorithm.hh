#ifndef UTOPIA_MODELS_PCPVERTEX_ALGORITHM_HH
#define UTOPIA_MODELS_PCPVERTEX_ALGORITHM_HH

namespace Utopia {
namespace Models {
namespace PCPVertex {

/// Determine step size to reach the line minimum along direction of update
/** See www.acclab.helsinki.fi/~aakurone/atomistiset/lecturenotes/lecture12_2up.pdf
 *  for details of algorithm
 *  
 *  \param dt           The step size of previous step
 *  \param energy_0     The energy at current position
 *  \param tolerance    Tolerance to find the minimum
 *  \return {step size to reach line minimum, energy at line minimum}
 */
std::pair<double, double> PCPVertex::determine_timestep (
        const double dt, const double energy_0,
        const double tolerance) const
{
    /// Capture nan energy values
    const std::function<double(double)> get_energy = [this](double beta) {
        double energy = this->get_energy(beta);
        if (isnan(energy)) {
            return std::numeric_limits<double>::max();
        }
        return energy;
    };

    // check that energy is actually decreasing in direction of update
    double tmp_energy = get_energy(1e-8);
    if (tmp_energy >= energy_0) {
        this->_log->trace("Already in minimum. At step size of 1e-8 the energy "
            "along direction of update increased by {}.",
            tmp_energy - energy_0);

        return std::make_pair(0., energy_0);
    }

    // find a minimum that is closer than a local maxima, i.e.
    // E_2 (dt = pos_2) > E_0(dt = 0) > E_min (dt = pos_min)
    // Start search with previous dt
    double pos_2 = std::max(std::min(2*dt, 2.), 1e-8); // right boundary
    double energy_2 = get_energy(pos_2); // energy at right boundary
    double pos_min = pos_2 / 2.; // the minimum
    double energy_min = get_energy(pos_min); // energy at minimum
    while (true) {
        if (pos_min > 3000.) {
            if (fabs(energy_2 - energy_0) < tolerance) {
                // the energy function is flat
                return std::make_pair(0., energy_0);
            }
            this->_log->error("At step of {} along direction of "
                "update the energy is still decreasing by {}", pos_min, 
                energy_2 - energy_0);
            for (double beta = 1e-11; beta < 1000.; beta *= 2) {
                double max_displ = 0;
                for (auto v : _am.vertices()) {
                    max_displ = std::max(max_displ,
                                         arma::norm(v->state.f) * beta);
                }
                this->_log->error("DEBUG At step of {} along direction of "
                    "update the energy is changing by {} at max displacement "
                    "of {}", beta, get_energy(beta) - energy_0,
                    max_displ);
            }
            // NOTE only if energy_2 < energy_1: dt -> 2*dt
            throw std::runtime_error("Unable to find bracket to "
                "local energy minimum!");
        }
        if (pos_2 < 1e-10) {
            // the energy is increasing, hence already at local minimum
            // because only if energy_min > energy_2: pos_2 -> pos_2 / 2
            this->_log->trace("Already in minimum. At step size of {} the "
                "energy along direction of update increased by {}.",
                pos_2, energy_2 - energy_0);
            return std::make_pair(0., energy_0);
        }

        if (energy_2 < energy_0) {
            // pos_2 is not a right boundary, since smaller than E_0
            // Use pos_2 as new minimum
            // and increase interval to find a right boundary
            pos_min = pos_2;
            energy_min = energy_2;

            // increase interval
            pos_2 *= 2.;
            energy_2 = get_energy(pos_2);
            continue;
        }

        if (energy_min > energy_0) {
            // pos_min defines a closer right bracket than pos_2
            // we know that close to pos_0 there is a minimum, keep looking
            pos_2 = pos_min;
            energy_2 = energy_min;

            pos_min = pos_2 / 2.;
            energy_min = get_energy(pos_min);
            continue;
        }

        // Both conditions fulfilled!
        // there is a minimum within interval [0, pos_2]
        break;
    }

    // is minimum already good enough?
    if (fabs(energy_min - energy_0) < tolerance) {
        // the found minimum fulfills our condition
        return std::make_pair(pos_min, energy_min);
    }

    // Improve knowledge of minimum by decreasing the interval size
    double pos_1 = 0.; // the left boundary
    double energy_1 = energy_0; // energy at left boundary

    // Find the minimum between brackets with parabolic interpolation
    while(true) {
        double term_3 = (pos_2 - pos_1)*(energy_2 - energy_min);
        double term_4 = (pos_2 - pos_min)*(energy_2 - energy_1);
        double term_1 = (pos_2 - pos_1)*term_3;
        double term_2 = (pos_2 - pos_min)*term_4;

        // the minimum of parabola through 1, 2, min
        double pos_3 = pos_2 - 0.5*(term_1 - term_2)/(term_3 - term_4);
        if (pos_3 < pos_1 or pos_3 > pos_2) {
            this->_log->error("pos_1={}, pos_min={}, pos_2={}",
                              pos_1, pos_min, pos_2);
            this->_log->error("D_energy_1={}, D_energy_min={}, "
                              "D_energy_2={}. wrt energy_0",
                              energy_0-energy_1, energy_0-energy_min,
                              energy_0-energy_2);
            this->_log->error("Fitted minimum: x={}", pos_3);
            throw std::runtime_error("Energy minimization failed! "
                                     "Parabola fit outside brackets");
        }

        // if the found minimum is precise enough
        double energy_3 = get_energy(pos_3);
        if (fabs(energy_3 - energy_min) < tolerance) {
            // found the minimum with required precision
            if (energy_3 < energy_min) {
                return std::make_pair(pos_3, energy_3);
            }
            else {
                return std::make_pair(pos_min, energy_min);
            }
        }
        
        if (energy_3 > energy_0) {
            // there is a closer maximum. Retry with smaller step size
            return determine_timestep(pos_3 / 2., energy_0, tolerance);
        }

        // Close brackets in on new minimum
        if (energy_3 < energy_min) {
            // pos_3 is to the left of minimum
            if (pos_3 < pos_min) {
                // minimum becomes right bracket
                pos_2 = pos_min; energy_2 = energy_min;
            }
            // else is to the right of minimum
            else {
                // the new left bracket
                pos_1 = pos_min; energy_1 = energy_min;
            }

            // the new minimum
            pos_min = pos_3; energy_min = energy_3;
        }
        // the estimate is not a better minimum, use as new bracket
        else {
            // minimum stays as is

            // pos_3 to the left of minimum
            if (pos_3 < pos_min) {
                // pos_3 is a better left bracket
                pos_1 = pos_3; energy_1 = energy_3;
            }
            // else to the right of minimum
            else {
                // pos_3 is a better right bracket
                pos_2 = pos_3; energy_2 = energy_3;
            }
        }

        // repeat with decreased brackets
    }
};

/// Initialisation of the energy minisation process
void PCPVertex::init_minimization ()
{
    _energy = get_energy();

    set_gradient();

    if (this->_update_scheme != UpdateScheme::ConjugateGradient) {
        return;
    }

    apply_rule<Update::sync>(
        [this](const auto& vertex) {
            auto state = vertex->state;
            state.g = state.f;
            state.h = state.f;
            return state;
        },
        _am.vertices()
    );

    _dt = _default_minimization_params.dt;
};

/// Single step in direction of steepest gradient
/** \param adative_step     If true, a line minimization along steepest gradient
 *                          is performed. Else, update with fixed stepsize.
 *  \return the energy after upate
 */
double PCPVertex::steepest_gradient_step (bool adaptive_step)
{
    set_gradient();

    double new_energy;
    if (adaptive_step) {
        std::tie(_dt, new_energy) = determine_timestep(
            _dt, _energy, _minimization_tolerance);
        double energy_change = get_rel_energy_change(new_energy, _energy);

        if (energy_change < -1e-14) {
            this->_log->trace("Updating with timestep {} at relative energy "
                              "change {}",
                              _dt, energy_change);
        }
        else {
            this->_log->trace("NOT updating with step size {} along direction "
                              "of update at energy change {}", _dt,
                              energy_change);
            return this->_energy;
        }
    }
 
    apply_rule<Update::sync>(update_position, _am.vertices());

    if (not adaptive_step and _distr_temperature.param().stddev() > 0) {
        apply_rule<Update::sync>(update_brownian_motion, _am.vertices());        
    }

    if (not adaptive_step) {
        return get_energy();
    }
    else {
        return new_energy;
    }
};

/// Single step in direction of conjugate gradient
/** \return the energy after upate
 */
double PCPVertex::conjugate_gradient_step ()
{
    // line minimization along direction of update h
    double new_energy;
    std::tie(_dt, new_energy) = determine_timestep(
        _dt, _energy, _minimization_tolerance);
    double energy_change = get_rel_energy_change(new_energy, _energy);

    if (energy_change < -1e-14) {
        this->_log->trace("Updating with timestep {} at energy change {}",
                            _dt, energy_change);
        
        apply_rule<Update::sync>(update_position, _am.vertices());
    }
    else {
        this->_log->trace("NOT updating with step size {} along direction "
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
double PCPVertex::perform_update_step(
        UpdateScheme update_scheme)
{
    if (update_scheme == UpdateScheme::SteepestGradient) {
        return steepest_gradient_step(false);
    }
    else if (update_scheme == UpdateScheme::SteepestGradientAdaptive) {
        return steepest_gradient_step(true);
    }
    else if (update_scheme == UpdateScheme::ConjugateGradient) {
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