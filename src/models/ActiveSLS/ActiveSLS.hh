#ifndef UTOPIA_MODELS_ACTIVESLS_HH
#define UTOPIA_MODELS_ACTIVESLS_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>


namespace Utopia {
namespace Models {
namespace ActiveSLS {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<>;


// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// The ActiveSLS Model; the bare-basics a model needs
/** TODO Add your class description here.
 *  ...
 */
class ActiveSLS:
    public Model<ActiveSLS, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<ActiveSLS, ModelTypes>;

    /// Data type of the group to write model data to, holding datasets
    using DataGroup = typename Base::DataGroup;

    /// Data type for a dataset
    using DataSet = typename Base::DataSet;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space

    // -- Members -------------------------------------------------------------
    double _length;
    double _activity;
    
    double _target_length;

    double _tau_activation;
    double _activation_rate;
    double _activation_rate_slope;
    double _threshold_tension;

    double _dissipation_tau;

    double _external_tension;

    double _maxwell_stiffness;
    double _maxwell_tau_viscous;
    double _spring_constant;
    double _spring_restlength;
    
    double _active_tension;
    double _background_activity;



    // .. Temporary objects ...................................................


    // .. Datasets ............................................................
    std::shared_ptr<DataSet> _dset_junction;    /// Properties of the junction
    std::shared_ptr<DataSet> _dset_tension;     /// Tensions
    std::shared_ptr<DataSet> _dset_activity;    /// Dataset of activity


public:
    // -- Model Setup ---------------------------------------------------------

    /// Construct the ActiveSLS model
    /** \param name             Name of this model instance; is used to extract
     *                          the configuration from the parent model and
     *                          set up a HDFGroup for this instance
     *  \param parent_model     The parent model this model instance resides in
     *  \param custom_cfg       A custom configuration to use instead of the
     *                          one extracted from the parent model using the
     *                          instance name
     */
    template<class ParentModel>
    ActiveSLS (
        const std::string& name,
        ParentModel& parent_model,
        const DataIO::Config& custom_cfg = {}
    )
    :
        Base(name, parent_model, custom_cfg),
        _length(get_as<double>("length", this->_cfg)),
        _activity(get_as<double>("activity", this->_cfg)),
        _target_length(get_as<double>("target_length", this->_cfg)),
        _tau_activation(get_as<double>("tau_activation", this->_cfg)),
        _activation_rate(get_as<double>("activation_rate", this->_cfg)),
        _activation_rate_slope(get_as<double>("activation_rate_slope", this->_cfg)),
        _threshold_tension(get_as<double>("threshold_tension", this->_cfg)),
        _dissipation_tau(get_as<double>("dissipation_tau", this->_cfg)),
        _external_tension(get_as<double>("external_tension", this->_cfg)),
        _maxwell_stiffness(get_as<double>("maxwell_stiffness", this->_cfg)),
        _maxwell_tau_viscous(get_as<double>("maxwell_tau_viscous", this->_cfg)),
        _spring_constant(get_as<double>("spring_constant", this->_cfg)),
        _spring_restlength(get_as<double>("spring_restlength", this->_cfg)),
        _active_tension(get_as<double>("active_tension", this->_cfg)),
        _background_activity(get_as<double>("background_activity", this->_cfg)),

        _dset_junction(this->create_dset("junction", {2})),
        _dset_tension(this->create_dset("tension", {4})),
        _dset_activity(this->create_dset("activity", {}))
    {
        // Add attributes to junction dataset that provide coordinates
        _dset_junction->add_attribute("dim_name__1", "property");
        _dset_junction->add_attribute("coords_mode__property", "values");
        _dset_junction->add_attribute("coords__property",
            std::vector<std::string>{
                "length", "target_length"
            });
        // Add attributes to density dataset that provide coordinates
        _dset_tension->add_attribute("dim_name__1", "kind");
        _dset_tension->add_attribute("coords_mode__kind", "values");
        _dset_tension->add_attribute("coords__kind",
            std::vector<std::string>{
                "maxwell_elastic", "active", "elastic", "total"
            });
        this->_log->debug("Added coordinates to densities dataset.");
    }


private:
    // .. Setup functions .....................................................

    // .. Helper functions ....................................................
    std::function<double(double)> activation_rate = [this](double tension) {
        return (
              this->_activation_rate
            + exp(-_activation_rate_slope * (tension - _threshold_tension))
        );
    };

public:
    // -- Public Interface ----------------------------------------------------
    
    inline double tension__maxwell_elastic() {
        return _maxwell_stiffness * (_length - _target_length);
    }

    inline double tension__active() {
        return _active_tension * (_activity - _background_activity);
    }

    inline double tension__elastic() {
        return _spring_constant * (_length - _spring_restlength);
    }
    
    double tension() {        
        return (
              tension__maxwell_elastic() 
            + tension__active()
            + tension__elastic()
        );
    }


    // .. Simulation Control ..................................................

    /// Iterate a single step
    /** \details Here you can add a detailed description what exactly happens 
      *         in a single iteration step
      */
    void perform_step () {
        double T = tension();

        _activity += (
              1. 
            - _activity * this->activation_rate(tension())
        ) / _tau_activation;

        _target_length += (_length - _target_length) / _maxwell_tau_viscous;
        _length += (_external_tension - T) / _dissipation_tau;

        this->_log->debug("tension: {}", T);
        this->_log->debug("length: {}", _length);
        this->_log->debug("target length: {}", _target_length);
    }

    /// The prolog
    void prolog () {
        this->__prolog();

        double activity_equilibrium = 1. / activation_rate(_external_tension);
        double active_tension_eq = _active_tension * (  activity_equilibrium
                                                      - _background_activity);

        _dset_tension->add_attribute(
            "maximum_barrier_height", 
             active_tension_eq - _external_tension
        );
    }


    /// Monitor model information
    void monitor () {
        this->_monitor.set_entry("length", _length);
        this->_monitor.set_entry(
            "maxwell_length_difference",
            _length - _target_length
        );
        this->_monitor.set_entry(
            "excess_tension",
            tension() - _external_tension
        );
        this->_monitor.set_entry("activity", _activity);
    }


    /// Write data
    void write_data () {
        _dset_junction->write(std::vector<double>({
            _length,
            _target_length
        }));
        _dset_tension->write(std::vector<double>({
            tension__maxwell_elastic(),
            tension__active(),
            tension__elastic(),
            tension()
        }));
        _dset_activity->write(_activity);
    }


    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

};

} // namespace ActiveSLS
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_ACTIVESLS_HH
