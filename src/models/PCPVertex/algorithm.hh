#ifndef UTOPIA_MODELS_SAVANNAHETEROGENEOUS_ALGORITHM_HH
#define UTOPIA_MODELS_SAVANNAHETEROGENEOUS_ALGORITHM_HH

namespace Utopia {
namespace Models {
namespace PCPVertex {

/// The step size to reach the line minimum along steepest gradient
/** See www.acclab.helsinki.fi/~aakurone/atomistiset/lecturenotes/lecture12_2up.pdf
 */
template <bool periodic_bc, bool polarity_proteins>
std::pair<double, double> PCPVertex<periodic_bc,
                                    polarity_proteins>::determine_timestep (
        const double energy_0) const
{
    // Bracket the minimum
    double dt = std::min(std::max(_dt, 1e-5), 1e-2);

    double pos_1 = 0.; // left boundary
    bool calculate_energy_min = true;
    double energy_1 = energy_0;
    bool calculate_energy_2 = true;
    double pos_2 = dt; // right boundary
    double energy_2 = this->get_energy(_edges, _cells, pos_2);
    double pos_min = dt/2.;
    double energy_min = this->get_energy(_edges, _cells, pos_min);
    while (true) {
        if (dt > 3.) {
            if (fabs(energy_2 - energy_1) < _minimisation_precision) {
                // the energy function is flat
                return std::make_pair(0., energy_0);
            }
            this->_log->warn("At step of {} along direction of "
                "update the energy is still decreasing by {}", dt/2, 
                energy_2 - energy_1);
            for (double dt = 1e-11; dt < 1000.; dt *= 2) {
                this->_log->warn("At step of {} along direction of "
                    "update the energy is changing by {}", dt, 
                    this->get_energy(_edges, _cells, dt) - energy_1);
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
            energy_2 = this->get_energy(_edges, _cells, dt);
            continue; // dt was not a right boundary to minimum
        }

        if (energy_min > energy_1) {
            dt = dt/2;
            energy_2 = energy_min;
            energy_min = this->get_energy(_edges, _cells, dt/2);
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

        double energy_4 = this->get_energy(_edges, _cells, pos_4);
        if (fabs(energy_4 - energy_min) < _minimisation_precision) {
            dt = pos_4;
            break; // found the minimum with required precision
        }
        
        if (energy_4 - energy_0 > 0.) {
            throw std::runtime_error("Energy minimisation failed! "
                                        "Found closer local maximum, hence "
                                        "overestimated size of brackets.");
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

template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc, polarity_proteins>::init_minimisation ()
{
    set_gradient();

    for (auto v : _vertices) {
        v->gx = v->x; v->gy = v->y;
        v->hx = v->x; v->hy = v->y;
    }
};

template <bool periodic_bc, bool polarity_proteins>
void PCPVertex<periodic_bc, polarity_proteins>::conjugant_gradient_step ()
{
    // line minimisation along direction of update h
    double new_energy;
    std::tie(_dt, new_energy) = determine_timestep(_energy_previous_step);
    _energy_change = (new_energy - _energy_previous_step) / new_energy;

    if (_energy_change < -1e-14) {
        this->_log->debug("Updating with timestep {} at energy change {}",
                            _dt, _energy_change);
        std::for_each(_vertices.begin(), _vertices.end(),
                        update_position);
    }
    else {
        this->_log->debug("NOT updating with step size {} along direction "
                            "of update at energy change {}", _dt,
                            _energy_change);
        return;
    }        

    // determine the conjugate gradient direction
    set_gradient();
    double gamma = 0.;
    double g_square = 0.;
    for (auto v : _vertices) {
        gamma += std::pow(v->fx, 2) + std::pow(v->fy, 2);
        g_square += std::pow(v->gx, 2) + std::pow(v->gy, 2);
        
        v->gx = v->fx; v->gy = v->fy;
    }
    gamma /= g_square;
    for (auto v : _vertices) {
        v->hx = v->gx + gamma * v->hx;
        v->hy = v->gy + gamma * v->hy;

        v->fx = v->hx; v->fy = v->hy; 
    }
};

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif