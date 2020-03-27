#ifndef UTOPIA_MODELS_NOTCHDELTA_HH
#define UTOPIA_MODELS_NOTCHDELTA_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/cell_manager.hh>
#include <utopia/core/apply.hh>

#include <utopia/models/Environment/Environment.hh>


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

    bool has_hair_neighbor;

    /// Construct the cell state from a configuration
    CellState(const DataIO::Config& cfg)
    :
        cell_type(progenitor),
        has_hair_neighbor(false)
    {}
};

/// State of the Environment model
struct EnvCellState : Environment::BaseEnvCellState {
    /// the level of atoh1
    double atoh1;

    /// Constructor a uniform background
    EnvCellState()
    :
        atoh1(0.)
    { }

    /// Constructor a uniform background
    EnvCellState(const DataIO::Config& cfg)
    :
        atoh1(get_as<double>("atoh1", cfg, 0.))
    { }

    ~EnvCellState() = default;

    /// Getter
    double get_env(const std::string& key) const {
        if (key == "atoh1") {
            return atoh1;
        }
        else {
            throw std::invalid_argument("No parameter '"+ key +
                                        "' available in EnvCellState!");
        }
    }

    /// Setter
    void set_env(const std::string& key, const double& value) {        
        if (key == "atoh1") {
            atoh1 = value;
        }
        else {
            throw std::invalid_argument("No parameter '"+ key +
                                        "' available in EnvCellState!");
        }
    }
};

using EnvModel = Environment::Environment<Environment::DummyEnvParam,
                                          EnvCellState>;
using EnvCell = EnvModel::CellManager::Cell;

/// The type of the link container of cells in the Environment model
template<class EntityContainerType>
struct EnvLinks {
    /// Link to the associated cell in Environment model
    std::shared_ptr<EnvCell> env;

    EntityContainerType neighbors;
};


/// Specialize the CellTraits type helper for this model
using CellTraits = Utopia::CellTraits<CellState, Update::manual, false,
                                      EmptyTag, EnvLinks>;


/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<DefaultRNG, WriteMode::managed>;


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

    using Cell = CellManager::Cell;

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

    /// The Environment model
    EnvModel _envm;

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
    bool _progenitors_depleted;

    bool _end_simulation;


public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the NotchDelta model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel, typename... Taskargs>
    NotchDelta (const std::string name, ParentModel &parent,
                Taskargs&&... taskargs)
    :
        // Initialize first via base model
        Base(name, parent, std::forward<Taskargs>(taskargs)...),

        // Now initialize the cell manager
        _cm(*this),
        _envm("Environment", *this, _cm),

        // Initialize model parameters
        _rate_ph(),
        _rate_ps(get_as<double>("rate_ps", this->_cfg)),
        _rate_atoh1(),
        _atoh1_threshold(get_as<double>("atoh1_threshold", this->_cfg)),
        _rate_swap(get_as<double>("rate_swap", this->_cfg)),
        
        _prob_distr(0., 1.),
        _progenitors_depleted(false),
        _end_simulation(false)
    {
        if (get_as<std::string>("mode", _cm.cfg()["neighborhood"]) != "empty") {
            this->_log->info("Setting up costum neighborhood from cell "
                "manager ..");
            for (auto c : _cm.cells()) {
                auto neighbors = this->get_cm()->neighbors_of(c);
                c->custom_links().neighbors.insert(
                    c->custom_links().neighbors.begin(),
                    neighbors.begin(),
                    neighbors.end());
            }
        }
        else {
            this->_log->info("No neighborhood set up from cell manager. "
                "Remember to define costum neighborhood for cells");
        }

        if (not this->_cfg["rate_ph"]) {
            throw std::invalid_argument("Missing cfg entry: Expected dict with "
                "key 'rate_ph'.");
        }
        _rate_ph.clear();
        for (int it = 0; true; it++) {
            if (not this->_cfg["rate_ph"]["rate_"+std::to_string(it)]) {
                break;
            }
            _rate_ph.push_back(get_as<double>("rate_"+std::to_string(it),
                                              this->_cfg["rate_ph"]));
        }
        if (_rate_ph.size() < 2) {
            throw std::invalid_argument("Missing cfg entry: Expected at least "
                "2 entries in dict 'rate_ph'!");
        }

        if (not this->_cfg["atoh1_suppression"]) {
            throw std::invalid_argument("Missing cfg entry: Expected dict with "
                "key 'atoh1_suppression'.");
        }
        _rate_atoh1.clear();
        for (int it = 0; true; it++) {
            if (not this->_cfg["atoh1_suppression"]["rate_"+std::to_string(it)])
            {
                break;
            }
            _rate_atoh1.push_back(get_as<double>("rate_"+std::to_string(it),
                                                 this->_cfg["atoh1_suppression"]));
        }
        if (_rate_atoh1.size() < 2) {
            throw std::invalid_argument("Missing cfg entry: Expected at least "
                "2 entries in dict 'rate_atoh1'!");
        }

        if (_atoh1_threshold <= 0.) {
            _atoh1_threshold = 1e-10;
        }

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
        auto atoh1 = cell->custom_links().env->state.atoh1;

        if (state.cell_type == CellType::progenitor) {
            int mapping;
            mapping = ceil((1 - atoh1/_atoh1_threshold) * 
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

    /// The suppression of atoh1 rule
    const RuleFunc suppress_atoh1 = [this](const auto& cell) {
        auto state = cell->state;
        auto env_state = cell->custom_links().env->state;

        // number of neighboring hair cells
        int ns_hair = 0;
        for (const auto& n : cell->custom_links().neighbors) {
            if (n->state.cell_type == CellType::hair) { ns_hair++; }
        }

        if (ns_hair == 0.) {
            return state;
        }

        // the rate of atoh1 change
        int mapping = std::min(ns_hair, int(_rate_atoh1.size() - 1));
        env_state.atoh1 /= _rate_atoh1[mapping];

        cell->custom_links().env->state = env_state;

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

        for (auto&& n : cell->custom_links().neighbors) {
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

        auto neighbors = cell->custom_links().neighbors;
        neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(),
                            [](auto n) { 
                                return n->state.cell_type == CellType::hair;
                            }),
                        neighbors.end());
        std::shuffle(neighbors.begin(), neighbors.end(), *this->_rng);

        if (not neighbors.empty()) {
            auto n = neighbors.back();
            std::swap(state, n->state);
            std::swap(cell->custom_links().env->state,
                      n->custom_links().env->state);
        }
        
        return state;
    };

public:
    // -- Public Interface ----------------------------------------------------
    
    // .. Simulation Control ..................................................

    /// Iterate a single step
    void perform_step () {
        if (_progenitors_depleted) {
            if (not _end_simulation) {
                this->_log->warn("All progenitors differentiated. Ending "
                        "simulation at time {}!", this->_time);
                _end_simulation = true;
            }
            return;
        }

        _envm.iterate();

        apply_rule<Update::sync>(suppress_atoh1, _cm.cells());
        apply_rule<Update::sync>(transition, _cm.cells());

        apply_rule<Update::sync>(T1_transition_tag, _cm.cells());
        apply_rule<Update::async>(T1_transition, _cm.cells(), *this->_rng);
    }

    /// Monitor model information
    void monitor () {
        auto densities = this->get_densities();
        if (densities[CellType::progenitor] == 0.) {
            _progenitors_depleted = true;
        }

        this->_monitor.set_entry("density_progenitor",
                                 densities[CellType::progenitor]);
        this->_monitor.set_entry("density_hair",
                                 densities[CellType::hair]);
        this->_monitor.set_entry("density_support",
                                 densities[CellType::support]);
    }

    void prolog () {
        _envm.prolog();
        return this->__prolog();
    }

    void epilog () {
        _envm.epilog();
        return this->__epilog();
    }


    // .. Getters and setters .................................................
    // Add getters and setters here to interface with other models
    /// Getter for density of different cell types
    std::vector<double> get_densities () const {
        std::vector<int> count(3, 0.);
        for (auto c : _cm.cells()) {
            int type = c->state.cell_type;
            count[int(c->state.cell_type)]++;
        }
        double num_cells = _cm.cells().size();
        return {count[0]/num_cells, count[1]/num_cells, count[2]/num_cells}; 
    }

    int get_hh_contacts() const {
        apply_rule<Update::sync>(T1_transition_tag, _cm.cells());

        int cnt = 0;
        for (auto c : _cm.cells()) {
            cnt += c->state.has_hair_neighbor;
        }

        return cnt / 2;
    }

    auto get_cm () const {
        return std::make_shared<CellManager>(this->_cm);
    }

    bool simulation_ended () const {
        return _end_simulation;
    }
};

} // namespace Utopia::Models::NotchDelta

#endif // UTOPIA_MODELS_COPYMEGRID_HH
