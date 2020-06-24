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

#include <utopia/models/Environment/Environment.hh>


namespace Utopia::Models::NotchDelta {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/// The type of a cell's state
struct CellState {
    /// The type of a cell
    enum StateType {
        progenitor,
        hair,
        support,
        inactive,
        num_states
    } cell_type;

    /// Whether has neighbor of type hair
    bool has_hair_neighbor;

    /// An ID denoting to which cluster this cell belongs
    std::size_t cluster_id;

    /// Construct the cell state from a configuration
    CellState()
    :
        cell_type(progenitor),
        has_hair_neighbor(false),
        cluster_id(0)
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
template<class CellContainer>
struct EnvLinks {
    /// Link to the associated cell in Environment model
    std::shared_ptr<EnvCell> env;

    /// The custom neighborhood
    CellContainer neighbors;
};


/// Specialize the CellTraits type helper for this model
using CellTraits = Utopia::CellTraits<CellState, Update::manual, true,
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

    /// Type of a cell
    using Cell = CellManager::Cell;

    /// The state type of a cell
    using CellType = CellState::StateType;

    /// Extract the type of the rule function from the CellManager
    /** This is a function that receives a reference to a cell and returns the 
      * new cell state. For more details, check out \ref Utopia::CellManager
      *
      * \note   Whether the cell state gets applied directly or
      *         requires a synchronous update depends on the update mode
      *         specified in the cell traits.
      */
    using RuleFunc = typename CellManager::RuleFunc;

    /// The type of a function returning the neighborhood of a cell
    using NBFuncCell = std::function<CellContainer<Cell>(
                                const std::shared_ptr<Cell>&)>;



private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _log, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The cell manager
    /** \note the neighborhood defined in the cell manager is copied to the 
     *        custom neighborhood in `custom_links` and should not be accessed
     *        in the cell manager!
     */
    CellManager _cm;

    /// The Environment model
    EnvModel _envm;

    /// The rate of progenitor to hair cell transition
    /** The first entry is for atoh1 levels above threshold, the latter entries
     *  are inverse-linearly mapped to atoh1 levels below threshold
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
    /// The incremental cluster tag
    unsigned int _cluster_id_cnt;

    /// The time cluster tags were updated last
    Time _cluster_id_time;

    /// A temporary container for use in cluster identification
    CellContainer<Cell> _cluster_members;

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
        Base(name, parent_model, custom_cfg, writer_args),

        // Now initialize the cell manager
        _cm(*this),
        _envm("Environment", *this, _cm),

        // Initialize model parameters
        _rate_ph(setup_as_vector(
            get_as<Config>("rate_ph", this->_cfg), "rate_ph")),
        _rate_ps(get_as<double>("rate_ps", this->_cfg)),
        _rate_atoh1(setup_as_vector(
            get_as<Config>("atoh1_suppression", this->_cfg),
            "atoh1_suppression")),
        _atoh1_threshold(get_as<double>("atoh1_threshold", this->_cfg)),
        _rate_swap(get_as<double>("rate_swap", this->_cfg)),
        
        _prob_distr(0., 1.),
        _cluster_id_cnt(),
        _cluster_id_time(1e9),
        _cluster_members()
    {
        // copy the cm neighborhood to custom links
        if (get_as<std::string>("mode", _cm.cfg()["neighborhood"]) != "empty") {
            for (auto c : _cm.cells()) {
                c->custom_links().neighbors = _cm.neighbors_of(c);
            }
        }

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

    /// Identify each cluster of hair cells
    RuleFunc _identify_cluster = [this](const auto& cell){
        if (cell->state.cluster_id != 0 or
            cell->state.cell_type != CellType::hair)
        {
            // already labelled, nothing to do. Return current state
            return cell->state;
        }
        // else: need to label this cell

        // Increment the cluster ID counter and label the given cell
        _cluster_id_cnt++;
        cell->state.cluster_id = _cluster_id_cnt;

        // Use existing cluster member container, clear it, add current cell
        auto& cluster = _cluster_members;
        cluster.clear();
        cluster.push_back(cell);

        // Perform the percolation
        for (unsigned int i = 0; i < cluster.size(); ++i) {
            // Iterate over all potential cluster members c, i.e. all
            // neighbors of cell cluster[i] that is already in the cluster
            for (const auto& nb : cluster[i]->custom_links().neighbors) {
                // If it is a hair cell that is not yet in the cluster, add it.
                if (    nb->state.cluster_id == 0
                    and nb->state.cell_type == CellType::hair)
                {
                    nb->state.cluster_id = _cluster_id_cnt;
                    cluster.push_back(nb);
                    // This extends the outer for-loop...
                }
            }
        }

        return cell->state;
    };

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
            int mapping = ceil((1 - atoh1/_atoh1_threshold) * 
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
    const RuleFunc suppress_atoh1 = [this](
            const auto& cell)
    {
        auto state = cell->state;
        auto env_state = cell->custom_links().env->state;

        if (state.cell_type == CellType::inactive) { return state; }

        // number of neighboring hair cells
        std::size_t nbs_hair = std::count_if(
            cell->custom_links().neighbors.begin(),
            cell->custom_links().neighbors.end(),
            [](const auto& nb) {
                return nb->state.cell_type == CellType::hair; });

        // the rate of atoh1 change
        auto mapping = std::min(nbs_hair, _rate_atoh1.size() - 1);
        env_state.atoh1 /= _rate_atoh1[mapping];

        cell->custom_links().env->state = env_state;

        return state;
    };
    
    /// The preparation stage for T1_transitions
    /** Tagges a hair cell, if has at least one hair cell neighbor
     */
    const RuleFunc T1_transition_tag = [this](const auto& cell)
    {
        auto state = cell->state;
        state.has_hair_neighbor = false;

        if (state.cell_type != CellType::hair) {
            return state;
        }

        for (auto&& n : cell->custom_links().neighbors) {
            if (n->state.cell_type == CellType::hair) {
                state.has_hair_neighbor = true;
                break;
            }
        }

        return state;
    };

    /// The T1 transition rule
    /** Tag cells using NotchDelta::T1_transition_tag.
     *  If tagges, it swaps state with a random non-hair cell typed neighbor
     * 
     *  \note This is an asynchronous rule! It cannot be applied synchronously.
     *  \note This rule is deterministic
     */
    const RuleFunc T1_transition = [this](
            auto& cell){
        auto state = cell->state;

        state = T1_transition_tag(cell);

        if (not state.has_hair_neighbor)
        {
            return state;
        }

        // remove the hair cells from neighbors
        auto neighbors = cell->custom_links().neighbors;
        neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(),
                [](auto n) {
                    return (n->state.cell_type == CellType::hair);
                }),
            neighbors.end());

        if (not neighbors.empty()) {
            // select a random neighbor
            std::shuffle(neighbors.begin(), neighbors.end(), *this->_rng);
            auto n = neighbors.back(); 

            // swap the state and the env->state (atoh1, etc.)
            std::swap(state, n->state);
            std::swap(cell->custom_links().env->state,
                      n->custom_links().env->state);
            // NOTE keep the geometric properties like neighbors
        }
        
        return state;
    };

public:
    // .. Helper functions ....................................................
    /// Identify clusters
    /** This function identifies clusters and updates the cell
     *  specific cluster_id as well as the member variable 
     *  cluster_id_cnt that counts the number of ids
     * 
     *  \note This function tracks the last time it was applied and is hence not
     *        applied twice at the same timepoint
     */
    void identify_clusters(){
        if (_cluster_id_time == this->_time) {
            // already applied in this step
            return;
        }

        this->_log->debug("Identifying cluster ids");

        // reset cluster counter
        _cluster_id_cnt = 0;
        apply_rule<Update::sync>(
            [](const auto& cell) {
                cell->state.cluster_id = 0;
                return cell->state; },
            _cm.cells() );
        
        apply_rule<Update::async, Shuffle::off>(_identify_cluster, 
                                                _cm.cells());

        _cluster_id_time = this->_time;
    }

    // -- Public Interface ----------------------------------------------------
    
    // .. Simulation Control ..................................................

    /// Iterate a single step
    /** \details Performs the following rules
     *      -# NotchDelta::suppress_atoh1
     *      -# NotchDelta::transition
     *      -# NotchDelta::T1_transition
     */
    void perform_step () {
        _envm.iterate();

        apply_rule<Update::sync>(suppress_atoh1, _cm.cells());
        apply_rule<Update::sync>(transition, _cm.cells());

        apply_rule<Update::async>(
            T1_transition,
            select_entities<SelectionMode::probability>(_cm, _rate_swap),
            *this->_rng);
    }

    /// Monitor model information
    void monitor () {
        auto densities = this->get_densities();

        this->_monitor.set_entry("density_progenitor",
                                 densities[CellType::progenitor]);
        this->_monitor.set_entry("density_hair",
                                 densities[CellType::hair]);
        this->_monitor.set_entry("density_support",
                                 densities[CellType::support]);
    }

    /// The custom prolog
    void prolog () {
        _envm.prolog();
        return this->__prolog();
    }

    //// The custom epilog
    void epilog () {
        _envm.epilog();
        return this->__epilog();
    }


    // .. Getters and setters .................................................
    // Add getters and setters here to interface with other models
    /// Getter for density of different cell types
    std::vector<double> get_densities () const {
        std::vector<int> count(CellType::num_states, 0.);
        for (auto c : _cm.cells()) {
            count[int(c->state.cell_type)]++;
        }
        double num_cells = _cm.cells().size() - count[int(CellType::inactive)];
        return {count[0]/num_cells, count[1]/num_cells, count[2]/num_cells}; 
    }

    /// Get the number of hair cells that have no contact to another hair cell
    unsigned int get_num_rosettes() const {
        apply_rule<Update::sync>(T1_transition_tag, _cm.cells());

        unsigned int cnt = 0;
        for (auto c : _cm.cells()) {
            if (c->state.cell_type == CellType::hair) {
                cnt += (not c->state.has_hair_neighbor);
            }
        }

        return cnt;
    }

    /// Getter for the cell manager
    const auto& get_cm () const {
        return this->_cm;
    }
};

} // namespace Utopia::Models::NotchDelta

#endif // UTOPIA_MODELS_COPYMEGRID_HH
