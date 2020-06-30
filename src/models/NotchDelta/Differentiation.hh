#ifndef UTOPIA_MODELS_DIFFERENTIATION_HH
#define UTOPIA_MODELS_DIFFERENTIATION_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/cell_manager.hh>
#include <utopia/core/apply.hh>
#include <utopia/core/select.hh>

namespace Utopia::Models::Differentiation {

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

    /// Whether forms a rosette
    /** A rosette is a hair cell with no hair cell neighbors
     */
    bool is_rosette;

    /// An ID denoting to which cluster this cell belongs
    std::size_t cluster_id;

    /// Construct the cell state from a configuration
    CellState()
    :
        cell_type(progenitor),
        is_rosette(false),
        cluster_id(0)
    { }
};


/// The type of the link container
template<class CellContainer>
struct CustomLinks {
    /// The custom neighborhood
    CellContainer neighbors;
};

/// Specialize the CellTraits type helper for this model
template<typename CellState>
using CellTraits = Utopia::CellTraits<CellState, Update::manual, true,
                                      EmptyTag, CustomLinks>;

/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<DefaultRNG, WriteMode::managed>;


// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// The base of Differentiation models; a good start for a Differentiation model
/** This model provides a base for models that want to describe the
 *  differentiation of progenitor cells to hair and support cells.
 */
template<typename CellTraits>
class Differentiation:
    public Model<Differentiation<CellTraits>, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<Differentiation<CellTraits>, ModelTypes>;

    /// The type of Time
    using Time = typename Base::Time;

    /// The type of a configuration
    using Config = typename Base::Config;

    /// Type of the CellManager to use
    using CellManager = Utopia::CellManager<CellTraits,
                                            Differentiation<CellTraits>>;

    /// Type of a cell
    using Cell = typename CellManager::Cell;

    /// The state of a cell
    using CellState = typename Cell::State;

    /// The cell's type (hair, support, etc.)
    using CellType = typename CellState::StateType;

    /// Extract the type of the rule function from the CellManager
    /** This is a function that receives a reference to a cell and returns the 
      * new cell state. For more details, check out \ref Utopia::CellManager
      *
      * \note   Whether the cell state gets applied directly or
      *         requires a synchronous update depends on the update mode
      *         specified in the cell traits.
      */
    using RuleFunc = typename CellManager::RuleFunc;



protected:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _log, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The cell manager
    /** \note the neighborhood defined in the cell manager is copied to the 
     *        custom neighborhood in `custom_links` and should not be accessed
     *        in the cell manager!
     */
    CellManager _cm;

    /// A re-usable uniform real distribution to evaluate probabilities
    std::uniform_real_distribution<double> _prob_distr;

private:
    // .. Temporary objects ...................................................
    /// The incremental cluster tag
    unsigned int _cluster_id_cnt;

    /// The time cluster tags were updated last
    Time _cluster_id_time;

    /// A temporary container for use in cluster identification
    CellContainer<Cell> _cluster_members;

public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the Differentiation model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel, class... WriterArgs>
    Differentiation (const std::string name, ParentModel &parent_model,
                const DataIO::Config& custom_cfg = {},
                std::tuple<WriterArgs...> writer_args = {})
    :
        // Initialize first via base model
        Base(name, parent_model, custom_cfg, writer_args),

        // Now initialize the cell manager
        _cm(*this),
        
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

        this->_log->debug("Differentiation model base constructor finished.");
    }


private:
    // .. Helper functions ....................................................
    /// Identify each cluster of hair cells
    RuleFunc _identify_cluster = [this](const auto& cell) {
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

protected:
    // .. Helper functions ....................................................
    /// Tag those cells that form a rosette
    /** I.e. hair cells that have no hair cell neighbors
     */
    const RuleFunc tag_rosettes = [this](const auto& cell)
    {
        auto state = cell->state;

        if (state.cell_type != CellType::hair) {
            state.is_rosette = false;
            return state;
        }

        state.is_rosette = true;
        for (auto&& n : cell->custom_links().neighbors) {
            if (n->state.cell_type == CellType::hair) {
                state.is_rosette = false;
                break;
            }
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
    void identify_clusters() {
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
    virtual void perform_step () { };

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

    /// Getter for density of protein densities
    virtual std::vector<double> get_protein_densities () const {
        return {};
    }

    /// Get the number of rosettes
    /** I.e. hair cells that have no contact to another hair cell
     */
    unsigned int get_num_rosettes() const {
        apply_rule<Update::sync>(tag_rosettes, _cm.cells());

        unsigned int cnt = 0;
        for (auto c : _cm.cells()) {
            cnt += c->state.is_rosette;
        }

        return cnt;
    }

    /// Getter for the cell manager
    const auto& get_cm () const {
        return this->_cm;
    }
};

} // namespace Utopia::Models::Differentiation

#endif // UTOPIA_MODELS_COPYMEGRID_HH
