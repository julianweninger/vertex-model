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
 *  \return {whether successful, step size to reach line minimum,
 *           energy at line minimum}
 */
std::tuple<bool, double, double> PCPVertex::determine_timestep (
        const double dt, const double energy_0,
        const double tolerance) const
{
    // whether an error occurred - infinite energy
    auto success = std::make_shared<bool>(true);

    // Capture nan energy values
    // nan energy is mapped to infinity
    // NOTE arise from bad prediction of area
    const std::function<double(double, std::shared_ptr<bool>)> get_energy =
    [this](double beta, std::shared_ptr<bool> success)
    {
        double energy = this->get_energy(beta);
        if (std::isnan(energy)) {
            *success = false;
            return std::numeric_limits<double>::max();
        }
        return energy;
    };

    // check that energy is actually decreasing in direction of update
    double tmp_energy = get_energy(1e-4, success);
    if (tmp_energy >= energy_0) {
        this->_log->trace("Already in minimum. At step size of 1e-4 the energy "
            "along direction of update increased by {}.",
            tmp_energy - energy_0);

        return std::make_tuple(*success, 0., energy_0);
    }

    // find a minimum that is closer than a local maxima, i.e.
    // E_2 (dt = pos_2) > E_0(dt = 0) > E_min (dt = pos_min)
    // Start search with previous dt
    double pos_2 = std::max(std::min(2*dt, 2.), 1e-8); // right boundary
    double energy_2 = get_energy(pos_2, success); // energy at right boundary
    double pos_min = pos_2 / 2.; // the minimum
    double energy_min = get_energy(pos_min, success); // energy at minimum
    while (true) {
        if (pos_min > 3000.) {
            if (fabs(energy_2 - energy_0) < tolerance) {
                // the energy function is flat
                return std::make_tuple(true, 0., energy_0);
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
                    "of {}", beta, get_energy(beta, success) - energy_0,
                    max_displ);
            }
            // NOTE only if energy_2 < energy_1: dt -> 2*dt
            throw std::runtime_error("Unable to find bracket to "
                "local energy minimum!");
        }
        if (pos_2 < 1e-8) {
            // the energy is increasing, hence already at local minimum
            // because only if energy_min > energy_2: pos_2 -> pos_2 / 2
            this->_log->trace("Already in minimum. At step size of {} the "
                "energy along direction of update increased by {}.",
                pos_2, energy_2 - energy_0);
            return std::make_tuple(true, 0., energy_0);
        }

        if (energy_2 < energy_0) {
            // pos_2 is not a right boundary, since smaller than E_0
            // Use pos_2 as new minimum
            // and increase interval to find a right boundary
            pos_min = pos_2;
            energy_min = energy_2;

            // increase interval
            pos_2 *= 2.;
            energy_2 = get_energy(pos_2, success);
            continue;
        }

        if (energy_min > energy_0) {
            // pos_min defines a closer right bracket than pos_2
            // we know that close to pos_0 there is a minimum, keep looking
            pos_2 = pos_min;
            energy_2 = energy_min;

            pos_min = pos_2 / 2.;
            energy_min = get_energy(pos_min, success);
            continue;
        }

        // Both conditions fulfilled!
        // there is a minimum within interval [0, pos_2]
        break;
    }

    // check that pos_2 is a true right boundary
    if (std::isnan(energy_2) or energy_2 > 1.e16) {
        if (std::isnan(energy_min) or energy_min > 1.e15) {
            return std::make_tuple(false, std::nan("1"), std::nan("1"));
        }

        *success = true;

        // choose a closer right boundary
        pos_2 = 0.5 * (pos_2 + pos_min);
        energy_2 = get_energy(pos_2, success);

        // failed to find right bracket
        if (energy_2 < energy_0) {
            return std::make_tuple(false, std::nan("1"), std::nan("1"));
        }

        // continue with new right bracket
    }

    // is minimum already good enough?
    if (fabs(energy_min - energy_0) < tolerance) {
        // the found minimum fulfills our condition
        return std::make_tuple(*success, pos_min, energy_min);
    }
    

    // Improve knowledge of minimum by decreasing the interval size
    double pos_1 = 0.; // the left boundary
    double energy_1 = energy_0; // energy at left boundary

    // Find the minimum between brackets with parabolic interpolation
    while(*success) {
        double term_3 = (pos_2 - pos_1)*(energy_2 - energy_min);
        double term_4 = (pos_2 - pos_min)*(energy_2 - energy_1);
        double term_1 = (pos_2 - pos_1)*term_3;
        double term_2 = (pos_2 - pos_min)*term_4;

        // the minimum of parabola through 1, 2, min
        double pos_3 = pos_2 - 0.5*(term_1 - term_2)/(term_3 - term_4);
        if (not std::isfinite(pos_3) or pos_3 < pos_1 or pos_3 > pos_2) {
            this->_log->error("pos_1={}, pos_min={}, pos_2={}",
                              pos_1, pos_min, pos_2);
            this->_log->error("D_energy_1={}, D_energy_min={}, "
                              "D_energy_2={}. wrt energy_0",
                              energy_0-energy_1, energy_0-energy_min,
                              energy_0-energy_2);
            this->_log->error("Fitted minimum: x={}", pos_3);

            // catch with steepest gradient step
            return std::make_tuple(false, std::nan("1"), std::nan("1"));
        }

        // if the found minimum is precise enough
        double energy_3 = get_energy(pos_3, success);
        if (fabs(energy_3 - energy_min) < tolerance) {
            // found the minimum with required precision
            if (energy_3 < energy_min) {
                return std::make_tuple(*success, pos_3, energy_3);
            }
            else {
                return std::make_tuple(*success, pos_min, energy_min);
            }
        }
        
        if (energy_3 > energy_0) {
            if (not *success) {
                return std::make_tuple(false, std::nan("1"), std::nan("1"));
            }

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

    // only reaches here when success == false
    return std::make_tuple(false, std::nan("1"), std::nan("1"));
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
        [](const auto& vertex) {
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
        bool success;
        std::tie(success, _dt, new_energy) = determine_timestep(
            _dt, _energy, _minimization_tolerance);
        double energy_change = get_rel_energy_change(new_energy, _energy);

        if (not success) {            
            this->init_minimization();
            // NOTE this resets _dt

            this->_log->warn("Failed to determine timestep with nan energy. "
                "Performing a steepest gradient step with dt = {}!", _dt);
            
            return steepest_gradient_step(false);
        }

        if (energy_change < -1e-14) {
            this->_log->trace("Updating with timestep {} at relative energy "
                              "change {}",
                              _dt, energy_change);
        }
        else {
            this->_log->trace("NOT updating with step size {} along direction "
                              "of update at energy change {}", _dt,
                              energy_change);
                              
            this->init_minimization();
            // NOTE this resets _dt

            return steepest_gradient_step(false);
        }
    }
 
    apply_rule<Update::sync>(update_position, _am.vertices());

    if (not adaptive_step) {
        if (_distr_temperature.param().stddev() > 1.e-12) {
            apply_rule<Update::sync>(update_brownian_motion, _am.vertices());
        }
        if (std::get<1>(_linetension_fluctuations) > 1.e-12) {
            apply_rule<Update::sync>(update_linetension_ornstein, _am.edges());
        }
        if (std::get<1>(_contractility_activity) > 1.e-12) {
            apply_rule<Update::sync>(update_edge_contractility, _am.edges());
        }
        if (std::get<1>(_area_fluctuations) > 1.e-12) {
            apply_rule<Update::sync>(update_area_preferential_ornstein,
                                     _am.cells());
        }
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
    bool success;
    std::tie(success, _dt, new_energy) = determine_timestep(
        _dt, _energy, _minimization_tolerance);
    double energy_change = get_rel_energy_change(new_energy, _energy);

    if (not success) {
        this->init_minimization();
        // NOTE this resets dt

        this->_log->debug("Failed to determine timestep with nan energy. "
            "Performing a steepest gradient step with dt = {}!", _dt);

        new_energy = steepest_gradient_step(false);
        this->init_minimization();
        
        return new_energy;
    }

    if (energy_change < -1e-14) {
        this->_log->trace("Updating with timestep {} at energy change {}",
                            _dt, energy_change);
        
        apply_rule<Update::sync>(update_position, _am.vertices());
    }
    else {
        this->init_minimization();
        // NOTE this resets dt

        this->_log->debug("Failed to determine timestep with dt = 0. or "
            "dE = {} >= 0."
            "Performing a steepest gradient step with dt = {}!",
            energy_change, _dt);

        new_energy = steepest_gradient_step(false);
        this->init_minimization();
        
        return new_energy;
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
    if (this->_apical_contractility) {
        const auto [apical,
                    basal] = this->_am.get_apical_basal_boundary_edges();

        for (const auto& [e, flip] : basal) {
            e->state.contractility_on = false;
        }
    }


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


/// Perform all topological transitions
/** Performs the topological transitions in this order: 
 *      -# T2 cell extrusion
 *      -# T1 edge remodelling
 * 
 *  \returns whether any transition was successfull.
 */
bool PCPVertex::perform_transitions(bool enabled = true)
{
    if (not enabled) {
        return false;
    }

    bool transition_occurred = false;

    // T2 transitions -- cell extrusion
    _num_T2s = 0;
    if (_enable_T2_transitions) {
        auto cells = _am.cells();
        std::shuffle(cells.begin(), cells.end(), *this->_rng);
        for (int i = cells.size() - 1; i >= 0; i--) {
            double area = _am.area_of(cells[i]);
            if (area < _T2_threshold) {
                this->_log->debug("Removing cell {} in T2 transition "
                                  "in step {}..", cells[i]->id(), this->_time);
                bool T2 = _am.remove_cell_T2(cells[i]);
                transition_occurred = (transition_occurred or T2);
                _num_T2s += T2;
            }
        }
    }

    // T1 transition -- neighborhood change
    _num_T1s = 0;
    _num_T1s_attempted = 0;
    if (_enable_T1_transitions) {
        auto edges = _am.edges();
        std::shuffle(edges.begin(), edges.end(), *this->_rng);
        for (int i = edges.size() - 1; i >= 0; i--) {
            auto& edge = edges[i];
            double length = _am.length_of(edge);
            if (length < _T1_threshold)
            {
                if (_prob_distr(*this->_rng) > _T1_probability) {
                    this->_log->trace("Skipping T1 transition on edge {} of "
                                      "length {} < {} because of freq {} < 1.",
                                      edge->id(), length, _T1_threshold,
                                      _T1_probability);
                    continue;
                }
                if (    edge->state.last_T1_attempt > 0
                    and (  this->_time - edge->state.last_T1_attempt
                         < _T1_timeout)) 
                {
                    this->_log->trace("Skipping T1 transition on edge {} of "
                                      "length {} < {} because of lifetime "
                                      "{} < {} (timeout).",
                                      edge->id(), length, _T1_threshold,
                                      this->_time - edge->state.last_T1_attempt,
                                      _T1_timeout);
                    continue;
                }

                this->_log->debug("Removing edge {} in T1 transition "
                                  "in step {}..", edge->id(), this->_time);
                bool T1 = _am.remove_edge_T1(edge,
                        _linetension, _edge_contractility,
                        [this](const AgentContainer<Edge>& es,
                                const AgentContainer<Cell>& cs) { 
                                    return this->get_energy(es, cs, 0.); },
                        _T1_separation, _T1_barrier, _prob_distr(*this->_rng));
                
                if (not T1) {
                    edge->state.last_T1_attempt = this->_time;
                }
                // else: edge was removed

                transition_occurred = (transition_occurred or T1);
                _num_T1s += T1;
                _num_T1s_attempted++;
            }
        }
    }

    if (_num_T2s > 0 or _num_T1s > 0 or _num_T1s_attempted > 0) {
        this->_log->info("Removed {} cell{} and {} edge{} ({} attempted) "
                            "in step {}",
                        _num_T2s, _num_T2s != 1 ? "s":"", 
                        _num_T1s, _num_T1s != 1 ? "s":"",
                        _num_T1s_attempted, this->_time);
    }
    _num_T1s_total += _num_T1s;
    _num_T1s_attempted_total += _num_T1s_attempted;
    _num_T2s_total += _num_T2s;

    return transition_occurred;
}

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif