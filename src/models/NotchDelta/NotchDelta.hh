#ifndef UTOPIA_MODELS_NOTCHDELTA_HH
#define UTOPIA_MODELS_NOTCHDELTA_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/cell_manager.hh>
#include <utopia/core/apply.hh>


namespace Utopia::Models::NotchDelta {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/// The type of a cell's state
struct CellState {
    /// The state
    enum StateType {
        progenitor,
        hair,
        support
    } cell_type;

    /// The level of atoh1;
    double atoh1;

    bool has_hair_neighbor;

    /// Construct the cell state from a configuration
    CellState(const DataIO::Config& cfg)
    :
        cell_type(progenitor),
        atoh1(0.),
        has_hair_neighbor(false)
    {}
};


/// Specialize the CellTraits type helper for this model
using CellTraits = Utopia::CellTraits<CellState, Update::manual>;


/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<>;


// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// The NotchDelta Model; a good start for a CA-based model
/** TODO Add your model description here.
 *  This model's only right to exist is to be a template for new models.
 *  That means its functionality is based on nonsense but it shows how
 *  actually useful functionality could be implemented.
 */
class NotchDelta:
    public Model<NotchDelta, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<NotchDelta, ModelTypes>;

    /// Data type of the group to write model data to, holding datasets
    using DataGroup = typename Base::DataGroup;

    /// Data type for a dataset
    using DataSet = typename Base::DataSet;

    /// Type of the CellManager to use
    using CellManager = Utopia::CellManager<CellTraits, NotchDelta>;
    // NOTE that it requires the model's type as second template argument

    /// Extract the type of the rule function from the CellManager
    /** This is a function that receives a reference to a cell and returns the 
      * new cell state. For more details, check out \ref Utopia::CellManager
      *
      * \note   Whether the cell state gets applied directly or
      *         requires a synchronous update depends on the update mode
      *         specified in the cell traits.
      */
    using RuleFunc = typename CellManager::RuleFunc;

    using CellType = CellState::StateType;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _log, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The cell manager
    CellManager _cm;

    /// The rate of progenitor to hair cell transition
    /** The first entry is for atoh1 levels above threshold, the latter entries
     *  are linearly mapped to atoh1 levels below threshold
     */
    std::vector<double> _rate_ph;

    /// The rate of progenitor to support cell transition
    double _rate_ps;

    /// The rate of atoh1 accumulation
    /** This is the mean rate.
     *  Values are uniformly distributed between [0, 2*rate] within each cell
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


    /// A re-usable uniform real distribution to evaluate probabilities
    std::uniform_real_distribution<double> _prob_distr;

    // .. Temporary objects ...................................................


    // .. Datasets ............................................................
    /// A dataset for storing all cells' state
    std::shared_ptr<DataSet> _dset_state;

    /// A dataset for storing all cells' atoh1 level
    std::shared_ptr<DataSet> _dset_atoh1;


public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the NotchDelta model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel>
    NotchDelta (const std::string name, ParentModel &parent)
    :
        // Initialize first via base model
        Base(name, parent),

        // Now initialize the cell manager
        _cm(*this),

        // Initialize model parameters
        _rate_ph(get_as<std::vector<double>>("rate_ph", this->_cfg)),
        _rate_ps(get_as<double>("rate_ps", this->_cfg)),
        _rate_atoh1(get_as<std::vector<double>>("rate_atoh1", this->_cfg)),
        _atoh1_threshold(get_as<double>("atoh1_threshold", this->_cfg)),
        _rate_swap(get_as<double>("rate_swap", this->_cfg)),
        
        _prob_distr(0., 1.),

        // Datasets
        _dset_state(this->create_cm_dset("cell_type", _cm)),
        _dset_atoh1(this->create_cm_dset("atoh1_level", _cm))
    {
        this->_log->debug("{} model fully set up.", this->_name);
    }


private:
    // .. Setup functions .....................................................
    
    // .. Helper functions ....................................................

    // .. Rule functions ......................................................
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
            if (_prob_distr(*this->_rng) < _rate_ph[std::max(mapping, 0)]) {
                state.cell_type = CellType::hair;
            }
            else if (_prob_distr(*this->_rng) < _rate_ps) {
                state.cell_type = CellType::support;
            }
        }
        
        return state;
    };

    /// A direct differentiation rule
    /** TODO
     */
    const RuleFunc transition_direct = [this](const auto& cell){
        auto state = cell->state;

        if (state.cell_type == CellType::progenitor) {
            if (_prob_distr(*this->_rng) < _rate_ph[0]) {
                state.cell_type = CellType::hair;
            }
            else if (_prob_distr(*this->_rng) < _rate_ps) {
                state.cell_type = CellType::support;
            }
        }
        
        return state;
    };

    /// The accumulation of atoh1 rule
    /** Atoh1 accumulates with mean per cell rate NotchDelta::_rate_atoh1.
     *  Values are uniformly distributed between [0, 2*rate] among cells
     */
    const RuleFunc accumulate_atoh1 = [this](const auto& cell) {
        auto state = cell->state;

        // number of neighboring hair cells
        int ns_hair = 0;
        for (const auto& n : this->_cm.neighbors_of(cell)) {
            if (n->state.cell_type == CellType::hair) { ns_hair++; }
        }

        // the rate of atoh1 change
        int mapping = std::min(ns_hair, int(_rate_atoh1.size() - 1));
        double d_atoh1 = 2 * _prob_distr(*this->_rng) * _rate_atoh1[mapping];
        state.atoh1 = std::max(state.atoh1 + d_atoh1, 0.);

        return state;
    };
    
    /// The preparation stage for T1_transitions
    /** Tagges all hair cells, that have at least one hair cell neighbor
     */
    const RuleFunc T1_transition_tag = [this](const auto& cell){
        auto state = cell->state;

        if (state.cell_type != CellType::hair) {
            return state;
        }

        for (auto&& n :  this->_cm.neighbors_of(cell)) {
            if (n->state.cell_type == CellType::hair) {
                state.has_hair_neighbor = true;
                return state;
            }
        }
        state.has_hair_neighbor = false;
        
        return state;
    };

    /// The T1 intercalation rule
    /** If cell is tagged from NotchDelta::T1_transition_tag, then it swaps
     *  state with a non-hair cell typed neighbor with probability
     *  NotchDelta::_rate_swap 
     * 
     *  \note This is an asynchronous rule! It cannot be applied synchronously.
     */
    const RuleFunc T1_transition = [this](auto& cell){
        auto state = cell->state;

        if (not state.has_hair_neighbor or 
                _prob_distr(*this->_rng) > _rate_swap)
        {
            return state;
        }

        auto neighbors = this->_cm.neighbors_of(cell);
        neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(),
                            [](auto n) { 
                                return n->state.cell_type == CellType::hair;
                            }),
                        neighbors.end());
        std::shuffle(neighbors.begin(), neighbors.end(), *this->_rng);

        if (not neighbors.empty()) {
            std::swap(state, neighbors.back()->state);
        }
        
        return state;
    };

public:
    // -- Public Interface ----------------------------------------------------
    // .. Simulation Control ..................................................

    /// Iterate a single step
    void perform_step () {
        apply_rule<Update::sync>(accumulate_atoh1, _cm.cells());
        apply_rule<Update::sync>(transition, _cm.cells());

        apply_rule<Update::sync>(T1_transition_tag, _cm.cells());
        apply_rule<Update::async>(T1_transition, _cm.cells(), *this->_rng);
    }

    /// Monitor model information
    void monitor () { }


    /// Write data
    void write_data () {
        // Write out the some_state of all cells
        _dset_state->write(_cm.cells().begin(), _cm.cells().end(),
            [](const auto& cell) {
                return int(cell->state.cell_type);
        });

        // Write out the some_trait of all cells
        _dset_atoh1->write(_cm.cells().begin(), _cm.cells().end(),
            [](const auto& cell) {
                return cell->state.atoh1;
        });
    }


    // .. Getters and setters .................................................
    // Add getters and setters here to interface with other models

};

} // namespace Utopia::Models::NotchDelta

#endif // UTOPIA_MODELS_COPYMEGRID_HH
