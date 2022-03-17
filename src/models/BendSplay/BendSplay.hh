#ifndef UTOPIA_MODELS_BENDSPLAY_HH
#define UTOPIA_MODELS_BENDSPLAY_HH

// standard library includes
#include <random>
#include <armadillo>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/cell_manager.hh>
#include <utopia/core/apply.hh>

#include "utils.hh"


namespace Utopia::Models::BendSplay {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// The type of a cell's state
struct CellState {
    double angle;

    /// Construct the cell state from a configuration
    CellState()
    :
        angle(0.)
    {}
};


/// Specialize the CellTraits type helper for this model
/** Specifies the type of each cells' state as first template argument and the 
  * update mode as second.
  *
  * See \ref Utopia::CellTraits for more information.
  */
using CellTraits = Utopia::CellTraits<CellState, Update::sync, true>;


/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<DefaultRNG, WriteMode::managed>;


// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// The BendSplay Model; a good start for a CA-based model
/** TODO Add your model description here.
 *  This model's only right to exist is to be a template for new models.
 *  That means its functionality is based on nonsense but it shows how
 *  actually useful functionality could be implemented.
 */
class BendSplay:
    public Model<BendSplay, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<BendSplay, ModelTypes>;

    /// Data type of the group to write model data to, holding datasets
    using DataGroup = typename Base::DataGroup;

    /// Data type for a dataset
    using DataSet = typename Base::DataSet;

    /// Type of the CellManager to use
    using CellManager = Utopia::CellManager<CellTraits, BendSplay>;

    using GridMat = arma::mat;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _log, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The cell manager
    // CellManager _cm;

    std::size_t _N_rho;
    std::size_t _N_phi;

    /// The angle of the polarity vectors
    /** The angle is measured wrt the S-I tissue axis, i.e. 0 corresponds
     *  to the normal vector (the radial vector in curved tissues) of the long
     *  axis (P-D).
     */
    GridMat _angles;

    /// The inner and the outer radius of the tissue
    std::pair<double, double> _radii;

    /// The opening angle of the circle section in rad
    double _opening_angle;

    /// The radial coordinates
    GridMat _coords_rho;

    /// The azimuthal coordinates
    GridMat _coords_phi;

    /// The numerical step size for Euler-Scheme
    double _dt;

    /// The lattice point distance in radial direction
    double _drho;

    /// The lattice point distance in azimuthal direction
    double _dphi;

    /// The elastic modulus
    double _K;

    /// The difference of splay vs bending modulus
    double _dK;

    /// A re-usable uniform real distribution to evaluate probabilities
    std::uniform_real_distribution<double> _prob_distr;

    std::normal_distribution<double> _angle_distr;

    // .. Temporary objects ...................................................


    // .. Datasets ............................................................
    std::shared_ptr<DataSet> _dset_angles;
    std::shared_ptr<DataSet> _dset_coords_x;
    std::shared_ptr<DataSet> _dset_coords_y;

    std::shared_ptr<DataSet> _dset_energy;
    std::shared_ptr<DataSet> _dset_orientation_mean;
    std::shared_ptr<DataSet> _dset_orientation_std;


public:
    // -- Model Setup ---------------------------------------------------------

    /// Construct the BendSplay model
    /** \param name             Name of this model instance; is used to extract
     *                          the configuration from the parent model and
     *                          set up a HDFGroup for this instance
     *  \param parent_model     The parent model this model instance resides in
     *  \param custom_cfg       A custom configuration to use instead of the
     *                          one extracted from the parent model using the
     *                          instance name
     */
    template<class ParentModel, typename... WriterArgs>
    BendSplay (
        const std::string& name,
        ParentModel& parent_model,
        const DataIO::Config& custom_cfg = {},
        std::tuple<WriterArgs...> &&writer_args = {}
    )
    :
        // Initialize first via base model
        Base(name, parent_model, custom_cfg, writer_args),

        // Now initialize the cell manager
        // _cm(*this),
        
        _N_rho(get_as<std::size_t>("N_rho", this->_cfg)),
        _N_phi(get_as<std::size_t>("N_phi", this->_cfg)),
        _angles(_N_rho, _N_phi, arma::fill::zeros),
        _radii(std::make_pair(
            get_as<double>("inner_radius", this->_cfg),
            get_as<double>("outer_radius", this->_cfg)
        )),
        _opening_angle(get_as<double>("opening_angle", this->_cfg)),
        _coords_rho(_N_rho, _N_phi, arma::fill::zeros),
        _coords_phi(_N_rho, _N_phi, arma::fill::zeros),
        _dt(get_as<double>("dt", this->_cfg)),
        _drho(  (std::get<1>(_radii) - std::get<0>(_radii))
              / static_cast<double>(_N_rho)),
        _dphi(_opening_angle / static_cast<double>(_N_phi)),
        _K(get_as<double>("K", this->_cfg)),
        _dK(get_as<double>("dK", this->_cfg)),


        // Initialize the uniform real distribution to range [0., 1.]
        _prob_distr(0., 1.),
        _angle_distr(get_as<double>("initial_orientation", this->_cfg),
                     get_as<double>("initial_orientation_stddev", this->_cfg)),

        _dset_angles(this->create_grid_dataset("angles")),
        _dset_coords_x(this->create_grid_dataset("coords_x")),
        _dset_coords_y(this->create_grid_dataset("coords_y")),

        _dset_energy(this->create_dset("energy", {})),
        _dset_orientation_mean(this->create_dset("orientation_mean", {})),
        _dset_orientation_std(this->create_dset("orientation_std", {}))
    {

        _angles.transform([this](double) {
            return this->_angle_distr(*this->get_rng()) + M_PI;
        });
        _angles.transform([](double&  val) { return constrain_angle(val); });

        if (get_as<bool>("init_symmetric", this->_cfg)) {
            for (std::size_t col = 0; col < _N_phi / 2; col++) {
                _angles.col(col) *= -1;
            }
            
            _angles.transform([](double&  val) { 
                return constrain_angle(val); 
            });
        }

        _coords_rho += (  _coords_rho.each_col()
                        + arma::linspace(std::get<0>(_radii) + _drho / 2.,
                                         std::get<1>(_radii) - _drho / 2.,
                                         _N_rho));
        
        if (_opening_angle > 2 * M_PI) {
            throw std::invalid_argument(fmt::format("Opening angle ({}) cannot "
                "be larger than 2Pi ({})", _opening_angle, 2 * M_PI));
        }

        _coords_phi += (  _coords_phi.each_row()
                        + arma::linspace(- _opening_angle / 2. + _dphi / 2.,
                                           _opening_angle / 2. - _dphi / 2.,
                                         _N_phi).t());

        // this->_log->debug("Grid initialised with angle {} "
        //                   "(wrt radial direction) and stddev {}.",
        //                   get_as<double>("angle", this->_cfg),
        //                   stddev > 1.e-12 ? stddev : 0.);
        

        // Initialization should be finished here.
        this->_log->debug("{} model fully set up.", this->_name);

        std::cout << free_energy() << std::endl;
    }


private:
    // .. Setup functions .....................................................
    std::shared_ptr<DataSet>
        create_grid_dataset(const std::string name,
                            const std::size_t compression_level=1,
                            const std::vector<hsize_t> chunksize={})
    {
        _log->debug("Creating dataset '{}' in group '{}' ...",
                    name, this->get_hdfgrp()->get_path());
        
        const hsize_t num_write_ops = this->get_remaining_num_writes();
        auto write_shape = std::vector<hsize_t>();

        const auto dset = this->get_hdfgrp()->open_dataset(
            name, {num_write_ops, _N_rho * _N_phi},
            chunksize, compression_level
        );

        _log->debug("Successfully created dataset '{}'.", name);

        dset->add_attribute("dim_name__0", "time");
        dset->add_attribute("coords_mode__time", "start_and_step");
        dset->add_attribute("coords__time",
                            std::vector<std::size_t>({_write_start,
                                                      _write_every}));
        dset->add_attribute("dim_name__1", "ids");
        dset->add_attribute("coords_mode__ids", "trivial");

        dset->add_attribute("N_rho", _N_rho);
        dset->add_attribute("N_phi", _N_phi);
        
        _log->debug("Added time dimension labels and coordinates to "
                    "dataset '{}'.", name);

        return dset;
    }

    // .. Helper functions ....................................................
    /// Calculate the (midpoint) derivative of x for phi
    /** using von Neumann zero boundary condition
     */
    GridMat d_dphi(const GridMat& X) const {
        // derivatives with periodic boundary conditions
        GridMat dphi_p = arma::shift(X, -1, 1) - X; // right difference
        GridMat dphi_m = X - arma::shift(X, 1, 1);  // left difference

        // von Neumann zero-BC
        dphi_m.col(0)        = arma::zeros(_N_rho);
        dphi_p.col(_N_phi-1) = arma::zeros(_N_rho);

        // round for 2 Pi
        dphi_p -= arma::round(dphi_p / (2 * M_PI)) * (2 * M_PI);
        dphi_m -= arma::round(dphi_m / (2 * M_PI)) * (2 * M_PI);

        // mid-point derivative
        return (dphi_p + dphi_m) / (2 * _dphi);
    }


    /// Calculate the (midpoint) second derivative of x for phi
    /** using von Neumann zero boundary condition
     */
    GridMat d2_dphi2(const GridMat& X) const {
        // derivatives with periodic boundary conditions
        GridMat dphi_p = arma::shift(X, -1, 1) - X; // right difference
        GridMat dphi_m = X - arma::shift(X, 1, 1);  // left difference

        // von Neumann zero-BC
        dphi_m.col(0)        = arma::zeros(_N_rho);
        dphi_p.col(_N_phi-1) = arma::zeros(_N_rho);

        // round for 2 Pi
        dphi_p -= arma::round(dphi_p / (2 * M_PI)) * (2 * M_PI);
        dphi_m -= arma::round(dphi_m / (2 * M_PI)) * (2 * M_PI);

        // mid-point derivative
        return (dphi_p - dphi_m) / std::pow(_dphi, 2);
    }

    /// Calculate the (midpoint) derivative of x for rho
    /** using periodic boundary condition
     */
    GridMat d_drho(const GridMat& X) const {
        // derivatives with periodic boundary conditions
        GridMat drho_p = arma::shift(X, -1, 0) - X; // upper difference
        GridMat drho_m = X - arma::shift(X, 1, 0); // lower difference

        // round for 2 Pi
        drho_p -= arma::round(drho_p / (2 * M_PI)) * (2 * M_PI);
        drho_m -= arma::round(drho_m / (2 * M_PI)) * (2 * M_PI);

        // mid-point derivative
        return (drho_p + drho_m) / (2 * _drho);
    }


    /// Calculate the (midpoint) second derivative of x for rho
    /** using periodic boundary condition
     */
    GridMat d2_drho2(const GridMat& X) const {
        // derivatives with periodic boundary conditions
        GridMat drho_p = arma::shift(X, -1, 0) - X; // upper difference
        GridMat drho_m = X - arma::shift(X, 1, 0); // lower difference

        // round for 2 Pi
        drho_p -= arma::round(drho_p / (2 * M_PI)) * (2 * M_PI);
        drho_m -= arma::round(drho_m / (2 * M_PI)) * (2 * M_PI);

        // mid-point derivative
        return (drho_p - drho_m) / std::pow(_drho, 2);
    }


    GridMat delta_free_energy(const GridMat& X) const {
        return -1. * (
            _K * (
                  d2_dphi2(X) / _coords_rho
                + _coords_rho % d2_drho2(X)
                + d_drho(X)
            )
            + _dK * (
                - 0.5 * arma::sin(2 * _angles) % arma::pow((1. + d_dphi(_angles)), 2) / _coords_rho
                + arma::sin(2 * _angles) % (1 + d_dphi(_angles)) % d_dphi(_angles) / _coords_rho
                + arma::pow(arma::sin(_angles), 2) % d2_dphi2(_angles) / _coords_rho
                - 0.5 * _coords_rho % arma::sin(2 * _angles) % arma::pow(d_dphi(_angles), 2)
                + _coords_rho % arma::pow(arma::cos(_angles), 2) % d2_drho2(_angles)
                + arma::pow(arma::cos(_angles), 2) % d_drho(_angles)
                + arma::cos(2 * _angles) % d_dphi(_angles) % d_drho(_angles)
                + arma::sin(2 * _angles) % d_drho(d_dphi(_angles))
            )
        );
    }


    // .. Rule functions ......................................................


public:
    // -- Public Interface ----------------------------------------------------
    // .. Simulation Control ..................................................

    /// Iterate a single step
    /** \details Here you can add a detailed description what exactly happens
      *         in a single iteration step
      */
    void perform_step () {
        _angles -= delta_free_energy(_angles) * _dt;
        
        _angles.transform([](double&  val) { return constrain_angle(val); });
    }


    /// Monitor model information
    /** \details Here, functions and values can be supplied to the monitor that
     *           are then available to the frontend. The monitor() function is
     *           _only_ called if a certain emit interval has passed; thus, the
     *           performance hit is small.
     *           With this information, you can then define stop conditions on
     *           frontend side, that can stop a simulation once a certain set
     *           of conditions is fulfilled.
     */
    void monitor () {
        this->_monitor.set_entry("energy", free_energy());
        this->_monitor.set_entry("orientaion", circular_mean(_angles));
        this->_monitor.set_entry("max_change", 
                                 arma::abs(delta_free_energy(_angles)).max());
    }


    // .. Getters and setters .................................................
    // Add getters and setters here to interface with other models
    double free_energy() const {
        return arma::accu(
            _K / 2 * (
                  arma::pow(1 + d_dphi(_angles), 2) / _coords_rho
                + _coords_rho % arma::pow(d_drho(_angles), 2)
            )
            + _dK / 2 * (
                  arma::pow(arma::sin(_angles), 2) % arma::pow(1 + d_dphi(_angles), 2) / _coords_rho
                + arma::sin(2 * _angles) % d_drho(_angles) % (1  + d_dphi(_angles))
                + _coords_rho % arma::pow(arma::cos(_angles) % d_drho(_angles), 2)
            )
        );
    }

    const arma::mat& get_orientations() const {
        return _angles;
    }

    arma::mat get_coords_x() const {
        return _coords_rho % arma::sin(_coords_phi);
    }

    arma::mat get_coords_y() const {
        return _coords_rho % arma::cos(_coords_phi);
    }
};

} // namespace Utopia::Models::BendSplay

#endif // UTOPIA_MODELS_BENDSPLAY_HH
