#ifndef UTOPIA_MODELS_COLLIER_HH
#define UTOPIA_MODELS_COLLIER_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/cell_manager.hh>
#include <utopia/core/apply.hh>
#include <utopia/core/select.hh>

#include "../NotchDelta/Differentiation.hh"

namespace Utopia::Models {
namespace Collier {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/// The type of a cell's state
struct CellState : Differentiation::CellState {
    /// The level of Notch
    /** The receptor
     */
    double notch;

    /// The level of Delta
    /** The ligand
     */
    double delta;

    /// The level of NICD
    /** The repressor: Notch intracellular domain, a cleavage product of Notch
     */
    double nicd;

    /// Construct the cell state from a configuration
    /** Initialize with very low concentration
     */
    template<class RNG>
    CellState(const DataIO::Config& cfg, const std::shared_ptr<RNG>& rng)
    :
        Differentiation::CellState()
    {
        notch = std::max(
                    std::normal_distribution<double>(
                        get_as<double>("notch", cfg),
                        get_as<double>("notch_sigma", cfg)
                    )(*rng),
                    1e-9);
        delta = std::max(
                    std::normal_distribution<double>(
                        get_as<double>("delta", cfg),
                        get_as<double>("delta_sigma", cfg)
                    )(*rng),
                    1e-9);
        nicd = std::max(
                std::normal_distribution<double>(
                        get_as<double>("nicd", cfg),
                        get_as<double>("nicd_sigma", cfg)
                    )(*rng),
                    1e-9);
    }
};

/// Specialize the CellTraits type helper for this model
using CellTraits = Utopia::CellTraits<CellState, Update::manual, false,
                                      EmptyTag, Differentiation::CustomLinks>;


// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// The Collier Model
/** A Differentiation model using 3 coupled differential equations for the
 *  concentrations of Notch, Delta and Repressor (NICD). These equations
 *  describe intracellular inhibition of Delta by the Repressor,
 *  trans-annihilation of Notch and Delta (cleavage and endocytosis) while
 *  up-regulating the production of the Repressor (due to the cleavage product),
 *  and cis-inactivation of Notch and Delta, 
 */
class Collier:
    public Differentiation::Differentiation<CellTraits>
{
// using typenames from Differentiation

private:
    // -- Members -------------------------------------------------------------
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _log, _space
    // _cm, _prob_distr
    // ... but you should definitely check out the documentation ;)

    /// The temporal resolution
    /** The time is normalised to the timescale of NICD degradation.
     *  The temporal resolution defines the number of steps per characteristic
     *  time of NICD (_degradation_nicd).
     */
    unsigned int _time_resolution;

    /// The nondimensional degradation rate of NICD
    /** \note This defines the timescale of the model
     */
    double _degradation_nicd;

    /// The nondimensional degradation rate of Notch
    /** \note This defines a timescale for Notch
     */
    double _degradation_notch;

    /// The nondimensional degradation rate of Delta
    /** \note This defines a timescale for Delta
     */
    double _degradation_delta;

    /// The nondimensional production rate of Notch
    double _production_notch;

    /// The nondimensional production rate of Delta
    double _production_delta;

    /// The nondimensional production rate of NICD
    double _production_nicd;
    
    /// The nondimensional strength of trans interaction
    double _interaction_strength_trans;

    /// The nondimensional strength of cis interaction
    double _interaction_strength_cis;

    /// The cooperativity of Delta inhibition
    double _cooperativity_delta;

    /// The cooperativity of NICD activity
    double _cooperativity_nicd;
    
    /// A factor for calculating _notch_threshold
    /** \f$ \theta = min + f * (max - min) \f$
     *  with \f$ \theta \f$ the threshold and min, max the minimum and maximum
     *  of notch levels accross cells.
     */
    double _notch_threshold_factor;

    // .. Temporary objects ...................................................
    /// The minimum notch level accross cells
    double _minimum_notch;

    /// The maximum notch level accross cells
    double _maximum_notch;
    
    /** If a cell's notch concentration falls below threshold it is considered
     *  a hair cell; otherwise it is considered a support cell
     */
    double _notch_threshold;


public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the Collier model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel, class... WriterArgs>
    Collier (const std::string name, ParentModel &parent_model,
                const DataIO::Config& custom_cfg = {},
                std::tuple<WriterArgs...> &&writer_args = {})
    :
        // Initialize first via base model
        Differentiation::Differentiation<CellTraits>(name, parent_model,
                                                     custom_cfg, writer_args),

        // Initialize model parameters
        _time_resolution(
            get_as<unsigned int>("time_resolution", this->_cfg)),
        _degradation_nicd(1. / static_cast<double>(_time_resolution)),

        _degradation_notch(
            get_as<double>("degradation_notch", this->_cfg)
            * _degradation_nicd),
        _degradation_delta(
            get_as<double>("degradation_delta", this->_cfg)
            * _degradation_nicd),
        
        _production_notch(
            get_as<double>("production_notch", this->_cfg)),
        _production_delta(
            get_as<double>("production_delta", this->_cfg)),
        _production_nicd(
            get_as<double>("production_nicd", this->_cfg)),
        
        _interaction_strength_trans(
            get_as<double>("interaction_strength_trans", this->_cfg)),
        _interaction_strength_cis(
            get_as<double>("interaction_strength_cis", this->_cfg)),
        
        _cooperativity_delta(
            get_as<double>("cooperativity_delta", this->_cfg)),
        _cooperativity_nicd(
            get_as<double>("cooperativity_nicd", this->_cfg)),
         
        _notch_threshold_factor(
            get_as<double>("notch_threshold_factor", this->_cfg))
        
    {
        if (_degradation_notch < 0) {
            throw std::invalid_argument(fmt::format("Model parameter "
                "'degradation_nicd' is a rate and hence cannot be "
                "negative, but was specified {} < 0!",
                _degradation_notch)); }
        if (_degradation_delta < 0) {
            throw std::invalid_argument(fmt::format("Model parameter "
                "'degradation_delta' is a rate and hence cannot be "
                "negative, but was specified {} < 0!",
                _degradation_delta)); }
        if (_production_notch < 0) {
            throw std::invalid_argument(fmt::format("Model parameter "
                "'production_notch' is a rate and hence cannot be "
                "negative, but was specified {} < 0!",
                _production_notch)); }
        if (_production_delta < 0) {
            throw std::invalid_argument(fmt::format("Model parameter "
                "'production_delta' is a rate and hence cannot be "
                "negative, but was specified {} < 0!",
                _production_delta)); }
        if (_production_nicd < 0) {
            throw std::invalid_argument(fmt::format("Model parameter "
                "'production_nicd' is a rate and hence cannot be "
                "negative, but was specified {} < 0!",
                _production_nicd)); }
        if (_interaction_strength_trans < 0) {
            throw std::invalid_argument(fmt::format("Model parameter "
                "'interaction_strength_trans' is a rate and hence cannot be "
                "negative, but was specified {} < 0!",
                _interaction_strength_trans)); }
        if (_interaction_strength_cis < 0) {
            throw std::invalid_argument(fmt::format("Model parameter "
                "'interaction_strength_cis' is a rate and hence cannot be "
                "negative, but was specified {} < 0!",
                _interaction_strength_cis)); }
        if (_cooperativity_delta < 0) {
            throw std::invalid_argument(fmt::format("Model parameter "
                "'cooperativity_delta' is a rate and hence cannot be "
                "negative, but was specified {} < 0!",
                _cooperativity_delta)); }
        if (_cooperativity_nicd < 0) {
            throw std::invalid_argument(fmt::format("Model parameter "
                "'cooperativity_nicd' is a rate and hence cannot be "
                "negative, but was specified {} < 0!",
                _cooperativity_nicd)); }

        if (_time_resolution < 3) {
            this->_log->warn("Time resolution was chosen {}. However, it is "
                "recommended to be >= 3. It defines the number of steps to "
                "represent the characteristic time of the model, as defined by "
                "the decay of the repressor protein NICD");
        }

        apply_rule<Update::sync>(Utopia::ExecPolicy::par,
                                 determine_cell_type, _cm.cells());
        apply_rule<Update::sync>(
            Utopia::ExecPolicy::par,
            [](const auto& cell) {
                auto state = cell->state;
                if (state.cell_type == CellType::inactive) {
                    state.notch = 0.;
                    state.delta = 0.;
                    state.nicd = 0.;
                }
                return state;
            },
            _cm.cells()
        );

        this->_log->debug("{} model fully set up.", this->_name);
    }


private:
    // .. Setup functions .....................................................
    
    // .. Helper functions ....................................................

    // .. Rule functions ......................................................
    /// Update the cell's protein concentrations
    /** TODO put differential equations here
     */
    const RuleFunc update = [this](const auto& cell) {
        auto state = cell->state;

        // protein levels
        const double notch = state.notch;
        const double delta = state.delta;
        const double nicd = state.nicd;

        // mean protein levels in neighbors
        double notch_nbs = 0.;
        double delta_nbs = 0.;
        for (const auto& neighbor : cell->custom_links().neighbors) {
            notch_nbs += neighbor->state.notch;
            delta_nbs += neighbor->state.delta;
        }
        notch_nbs /= static_cast<double>(cell->custom_links().neighbors.size());
        delta_nbs /= static_cast<double>(cell->custom_links().neighbors.size());

        // interactions
        double trans_annihilation = _interaction_strength_trans * notch *
                                    delta_nbs;
        double cis_inactivation = _interaction_strength_cis * notch * delta;

        state.notch += _degradation_notch * (
                            _production_notch -
                            trans_annihilation -
                            cis_inactivation -
                            notch
                       );
        state.delta += _degradation_delta * (
                            _production_delta
                              / (1. + std::pow(nicd, _cooperativity_delta)) -
                            _interaction_strength_trans * notch_nbs * delta -
                            cis_inactivation -
                            delta
                       );
        state.nicd += _degradation_nicd * (
                            _production_nicd *
                              std::pow(trans_annihilation, _cooperativity_nicd) /
                              (1. + std::pow(trans_annihilation,
                                             _cooperativity_nicd)) -
                            nicd
                      );

        return state;
    };
    
    /// Determine the type of this cell
    /** If notch concentration falls below _notch_threshold, the cell is
     *  considered a hair cell and a support cell otherwise.
     */
    const RuleFunc determine_cell_type = [this](const auto& cell) {
        auto state = cell->state;
        if (state.notch < _notch_threshold) {
            state.cell_type = CellType::hair;
        }
        else {
            state.cell_type = CellType::support;
        }
        return state;
    };

public:
    // .. Helper functions ....................................................

    // -- Public Interface ----------------------------------------------------
    
    // .. Simulation Control ..................................................

    /// Iterate a single step
    /** \details Performs the following rules
     *      -# Collier::update
     *      -# Collier::determine_cell_type
     */
    void perform_step () {
        apply_rule<Update::sync>(Utopia::ExecPolicy::par, update, _cm.cells());

        _minimum_notch = 1.e8;
        _maximum_notch = 0.;
        for (const auto& cell : _cm.cells()) {
            double notch = cell->state.notch;
            _minimum_notch = std::min(_minimum_notch, notch);
            _maximum_notch = std::max(_maximum_notch, notch);
        }
        _notch_threshold = _minimum_notch +
                           _notch_threshold_factor * 
                                (_maximum_notch - _minimum_notch);

        apply_rule<Update::sync>(Utopia::ExecPolicy::par, 
                                 determine_cell_type, _cm.cells());
    }


    // .. Getters and setters .................................................
    // use interface of Base!

    /// Getter for density of protein densities
    /** In this order: notch, delta, nicd
     */
    std::vector<double> get_protein_densities () const {
        const auto& cells = _cm.cells();
        double notch = 0., delta = 0., nicd = 0.;
        for (const auto& cell : cells) {
            const auto state = cell->state;
            notch += state.notch;
            delta += state.delta;
            nicd += state.nicd;
        }
        std::size_t num_cells = std::count_if(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.cell_type != CellType::inactive; });
        notch /= static_cast<double>(num_cells);
        delta /= static_cast<double>(num_cells);
        nicd /= static_cast<double>(num_cells);
        return {notch, delta, nicd}; 
    }
};

} // namespace Collier
} // namespace Utopia::Models

#endif // UTOPIA_MODELS_COLLIER_HH
