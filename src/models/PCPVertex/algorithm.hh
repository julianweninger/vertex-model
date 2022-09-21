#ifndef UTOPIA_MODELS_PCPVERTEX_ALGORITHM_HH
#define UTOPIA_MODELS_PCPVERTEX_ALGORITHM_HH

namespace Utopia {
namespace Models {
namespace PCPVertex {

/// Initialisation of the energy minisation process
void PCPVertex::init_minimization ()
{
    _energy_buffer.clear();
    _energy_buffer.push_back(get_energy());
};

/// Single step in direction of steepest gradient
/** \param adative_step     If true, a line minimization along steepest gradient
 *                          is performed. Else, update with fixed stepsize.
 *  \return the energy after upate
 */
double PCPVertex::steepest_gradient_step ()
{
    compute_forces();

    apply_rule<Update::sync>(update_position, _am.vertices());

    return get_energy();
};

/// Select the chosen update scheme
/** \return the energy after update
 */
double PCPVertex::perform_update_step() {    
    return steepest_gradient_step();
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
                bool T1 = _am.remove_edge_T1(
                    edge,
                    [this](const AgentContainer<Vertex>& vs,
                           const AgentContainer<Edge>& es,
                           const AgentContainer<Cell>& cs)
                    {
                        return this->get_energy(vs, es, cs); 
                    },
                    _T1_separation,
                    _T1_barrier,
                    _prob_distr(*this->_rng),
                    this->_time
                );
                
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
    _T1_frequency_acc += _num_T1s / static_cast<double>(_am.edges().size());

    _num_T1s_attempted_total += _num_T1s_attempted;
    _T1_attempt_frequency_acc += (  _num_T1s_attempted
                                 / static_cast<double>(_am.edges().size()));
                                 
    _num_T2s_total += _num_T2s;
    _T2_frequency_acc += _num_T2s / static_cast<double>(_am.cells().size());

    return transition_occurred;
}

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif