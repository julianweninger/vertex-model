#ifndef UTOPIA_MODELS_NOTCHDELTA_HH
#define UTOPIA_MODELS_NOTCHDELTA_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/cell_manager.hh>
#include <utopia/core/apply.hh>
#include <utopia/core/select.hh>

#include "Differentiation.hh"

namespace Utopia::Models {
namespace NotchDelta {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/// The type of a cell's state
struct CellState : Differentiation::CellState {
    /// the level of atoh1
    double atoh1;

    /// Construct the cell state from a configuration
    CellState()
    :
        Differentiation::CellState(),
        atoh1(0.)
    { }
};

/// Specialize the CellTraits type helper for this model
using CellTraits = Differentiation::CellTraits<CellState>;


// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// The NotchDelta Model
/** A Differentiation model using T1 transitions and atoh1 lateral inhibition
 *  to model spatial patterning. 
 */
class NotchDelta:
    public Differentiation::Differentiation<CellTraits>
{
// using typenames from Differentiation

private:
    // -- Members -------------------------------------------------------------
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _log, _space
    // _cm, _prob_distr
    // ... but you should definitely check out the documentation ;)

    /// The rate of progenitor to hair cell transition
    /** The first entry is for atoh1 levels above threshold, the latter entries
     *  are inverse-linearly mapped to atoh1 levels below threshold
     */
    std::vector<double> _rate_ph;

    /// The rate of progenitor to support cell transition
    double _rate_ps;

    /// The rate of atoh1 suppression
    /** The entries are for 0 ... N neighbors of type hair
     */
    std::vector<double> _rate_atoh1;
    
    /// The threshold level of atoh1
    /** Above threshold level progenitor cells differentiate fastest to hair
     *  cell
     */
    double _atoh1_threshold;

    /// The rate of T1 intercalation at hair cell - hair cell contact
    /** The repulsion of hair - hair cell types induces T1 transitions -- 
     *  neighborhood exchange at cell-cell junction.
     *  This is the rate at which a hair cell swaps state with a non-hair cell,
     *  if it is in contact with another hair cell
     */
    double _rate_swap;

    /// The exponential distribution used for atoh1 production
    std::exponential_distribution<double> _exp_distr;

public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the NotchDelta model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel, class... WriterArgs>
    NotchDelta (const std::string name, ParentModel &parent_model,
                const DataIO::Config& custom_cfg = {},
                std::tuple<WriterArgs...> &&writer_args = {})
    :
        // Initialize first via base model
        Differentiation::Differentiation<CellTraits>(name, parent_model,
                                                     custom_cfg, writer_args),

        // Initialize model parameters
        _rate_ph(setup_as_vector(
            get_as<Config>("rate_ph", this->_cfg), "rate_ph")),
        _rate_ps(get_as<double>("rate_ps", this->_cfg)),
        _rate_atoh1(setup_as_vector(
            get_as<Config>("atoh1_suppression", this->_cfg),
            "atoh1_suppression")),
        _atoh1_threshold(get_as<double>("atoh1_threshold", this->_cfg)),
        _rate_swap(get_as<double>("rate_swap", this->_cfg)),
        
        _exp_distr(get_as<double>("lambda", this->_cfg))
    {
        if (_atoh1_threshold <= 0.) {
            throw std::invalid_argument(fmt::format("Value of "
                "'atho1_threshold' must be larger than 0, but was {}.",
                _atoh1_threshold));
        }

        this->_log->debug("{} model fully set up.", this->_name);
    }


private:
    // .. Setup functions .....................................................
    /// Extract a collection of rates as a vector from a config node
    /** Instead of defining `some_rate_vector: [0., 0., 0.]`, define it as
     *  `some_rate_vector: {rate_0: 0., rate_1: 0., rate_2: 0.}`.
     */
    std::vector<double> setup_as_vector(const Config& cfg,
            const std::string&& name)
    {
        std::vector<double> vec;
        for (std::size_t i = 0; cfg["rate_"+std::to_string(i)]; i++) {
            vec.push_back(get_as<double>("rate_" + std::to_string(i), cfg));
        }

        if (vec.size() == 0) {
            throw std::invalid_argument(fmt::format("Provide a minimum of 1 "
                "rate to initilize the rate vector {}. Expected entries of "
                "type 'rate_i' with i consecutive uints in [0, N].", name));
        }

        this->_log->debug("Set up rates vector {} with {} entr{}", name, 
                          vec.size(), vec.size() != 1 ? "ies" : "y");

        return vec;
    }
    
    // .. Helper functions ....................................................


    // .. Rule functions ......................................................
    /// Accumulate atoh1 in all cells
    /** Use an exponential distribution with mean 1 / lambda
     */
    const RuleFunc accumulate_atoh1 = [this](const auto& cell) {
        auto state = cell->state;
        state.atoh1 += _exp_distr(*this->_rng);
        return state;
    };
   
    /// The differentiation rule for progenitor cells
    /** A progenitor cell differentiates
     *      - to hair cell with probability depending on atoh1 level
     *      - to support cell with probability NotchDelta::_rate_ps
     */
    const RuleFunc transition = [this](const auto& cell){
        auto state = cell->state;

        if (state.cell_type == CellType::progenitor) {
            int mapping = ceil((1 - state.atoh1/_atoh1_threshold) * 
                               (_rate_ph.size() - 1));
            double rate_ph = _rate_ph[std::max(mapping, 0)];
            if (this->_prob_distr(*this->_rng) < rate_ph) {
                state.cell_type = CellType::hair;
            }
            else if (this->_prob_distr(*this->_rng) < _rate_ps) {
                state.cell_type = CellType::support;
            }
        }
        
        return state;
    };

    /// The suppression of atoh1 rule
    const RuleFunc suppress_atoh1 = [this](
            const auto& cell)
    {
        auto state = cell->state;

        if (state.cell_type == CellType::inactive) { return state; }

        // number of neighboring hair cells
        const auto& neighbors = cell->custom_links().neighbors;
        std::size_t nbs_hair = std::count_if(neighbors.begin(), neighbors.end(),
            [](const auto& nb) {
                return nb->state.cell_type == CellType::hair; });

        // the rate of atoh1 change
        auto mapping = std::min(nbs_hair, _rate_atoh1.size() - 1);
        state.atoh1 /= _rate_atoh1[mapping];

        return state;
    };

    /// The T1 transition rule
    /** Tag cells using Differentiation::tag_rosettes.
     *  If tagges, it swaps state with a random non-hair cell typed neighbor
     * 
     *  \note This is an asynchronous rule! It cannot be applied synchronously.
     */
    const RuleFunc T1_transition = [this](const auto& cell){
        auto state = tag_rosettes(cell);

        if (state.cell_type != CellType::hair or state.is_rosette)
        {
            return state;
        }

        // a list of non hair neighbors
        auto neighbors = cell->custom_links().neighbors;
        neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(),
                [cell](auto& n) {
                    if (n->state.cell_type == CellType::hair) {
                        return true;
                    }
                    for (auto& nn : n->custom_links().neighbors) {
                        if (    nn != cell
                            and nn->state.cell_type == CellType::hair) {
                            return true;
                        }
                    }
                    return false;
                }),
            neighbors.end());

        if (neighbors.empty()) {
            neighbors = cell->custom_links().neighbors;
            neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(),
                    [](auto& n) {
                        return (n->state.cell_type == CellType::hair);
                    }),
                neighbors.end());
        }

        if (not neighbors.empty()) {
            // select a random neighbor
            std::shuffle(neighbors.begin(), neighbors.end(), *this->_rng);
            auto n = neighbors.back(); 

            // swap the state with neighbor
            std::swap(state, n->state);
        }
        
        return state;
    };

public:
    // .. Helper functions ....................................................

    // -- Public Interface ----------------------------------------------------
    
    // .. Simulation Control ..................................................

    /// Iterate a single step
    /** \details Performs the following rules
     *      -# NotchDelta::accumulate_atoh1
     *      -# NotchDelta::suppress_atoh1
     *      -# NotchDelta::transition
     *      -# NotchDelta::T1_transition
     */
    void perform_step () {
        apply_rule<Update::sync>(accumulate_atoh1, _cm.cells());
        apply_rule<Update::sync>(suppress_atoh1, _cm.cells());
        apply_rule<Update::sync>(transition, _cm.cells());

        // pre-select only the fraction of non-rosette hair cells.
        // NOTE because the selection is fixed, from pair of hair cells in
        //      contact only one can migrate over the distance of a single cell.
        auto T1_cells = _cm.cells();
        T1_cells.erase(std::remove_if(T1_cells.begin(), T1_cells.end(),
                [this, _prob_distr{std::move(_prob_distr)}]
                (const auto& cell) mutable {
                    return (   (cell->state.cell_type != CellType::hair)
                            or (_prob_distr(*this->_rng) > this->_rate_swap)
                            or this->tag_rosettes(cell).is_rosette);
                }),
            T1_cells.end());
        apply_rule<Update::async>(T1_transition, T1_cells, *this->_rng);
    }


    // .. Getters and setters .................................................
    // use interface of Base!

    /// Getter for density of atoh1
    std::vector<double> get_protein_densities () const {
        const auto& cells = _cm.cells();
        double atoh1 = 0.;
        for (const auto& cell : cells) {
            atoh1 += cell->state.atoh1;
        }
        std::size_t num_cells = std::count_if(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.cell_type != CellType::inactive; });
        atoh1 /= static_cast<double>(num_cells);
        return {atoh1};
    }
};

} // namespace NotchDelta
} // namespace Utopia::Models

#endif // UTOPIA_MODELS_NOTCHDELTA_HH
