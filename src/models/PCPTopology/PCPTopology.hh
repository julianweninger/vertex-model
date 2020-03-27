#ifndef UTOPIA_MODELS_PCPTOPOLOGY_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>

#include "../PCPVertex/geometry.hh"
#include "../PCPVertex/PCPVertex.hh"
#include "../PCPVertex/initialisation.hh"
#include "../PCPVertex/energy.hh"
#include "../PCPVertex/algorithm.hh"
#include "../PCPVertex/transitions.hh"
#include "../PCPVertex/operations.hh"
#include "../PCPVertex/PCPVertex_write_tasks.hh"

#include "../NotchDelta/NotchDelta.hh"
#include "../NotchDelta/NotchDelta_write_tasks.hh"

#include "PCPTopology_write_tasks.hh"

#include <utopia/models/Environment/Environment.hh>


namespace Utopia {
namespace Models {
namespace PCPVertex {
// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// Parameter of the Environment model
struct EnvParam : Utopia::Models::Environment::BaseEnvParam
{
    /// The preferential area of hair cells
    /** Since the volume is experimentally constant, this maps the detachment
     *  of HCs from the base. Assuming uniform height of the HCs. The SCs 
     *  automatically fill the volume.
     */
    double area_preferential_hair;
    
    /// The preferential area of support cells
    /** Since the volume is experimentally constant, this maps the detachment
     *  of HCs from the base. Assuming uniform height of the HCs. The SCs 
     *  automatically fill the volume.
     */
    double area_preferential_support;

    EnvParam(const Utopia::DataIO::Config& cfg)
    :
        area_preferential_hair(Utopia::get_as<double>("area_preferential_hair",
                                                      cfg)),
        area_preferential_support(Utopia::get_as<double>("area_preferential"
                                                         "_support", cfg))
    { }

    ~EnvParam() = default;

    /// Getter
    double get_env(const std::string& key) const override {
        if (key == "area_preferential_hair") {
            return area_preferential_hair;
        }
        else if (key == "area_preferential_support") {
            return area_preferential_support;
        }
        throw std::invalid_argument("No access method for key '" + key
                                    + "' in EnvParam!");
    }

    /// Setter
    void set_env(const std::string& key,
                const double& value) override
    {
        if (key == "area_preferential_hair") {
            area_preferential_hair = value;
        }
        else if (key == "area_preferential_support") {
            area_preferential_support = value;
        }
        else {
            throw std::invalid_argument("No setter method for key '" + key
                                        + "' in EnvParam!");
        }
    }
};

using EnvCellState = Environment::DummyEnvCellState;
using EnvModel = Environment::Environment<EnvParam, EnvCellState>;

/// Type helper to define types used by the model
using PCPTopologyModelTypes = Utopia::ModelTypes<DefaultRNG,
                                                 WriteMode::managed>;

// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// The PCPTopology Model; the bare-basics a model needs
template <bool periodic_bc, bool polarity_proteins>
class PCPTopology:
    public Model<PCPTopology<periodic_bc, polarity_proteins>,
                 PCPTopologyModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<PCPTopology<periodic_bc, polarity_proteins>,
                                   PCPTopologyModelTypes>;

    /// The types of a cell
    using CellType = typename PCPVertex<periodic_bc,
                                        polarity_proteins>::CellType;

private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The Vertex model
    PCPVertex<periodic_bc, polarity_proteins> _vertex_model;

    /// A tolerance value for equilibrium
    double _equilibration_tolerance;

    /// Number of steps performed in VertexModel per iteration
    int _num_equilibration_steps;

    /// Number of max iterations performed in VertexModel before aborting
    int _max_equilibration_iterations;

    /// How often to shake the system during equilibration
    /** Equilibration processes might get stuck in local minima, hence the 
     *  equilibrated cellular arrangement is perturbed and equilibration is
     *  repeated.
     */
    int _num_jiggle_per_equilibration;
    
    /// The intensity of jiggle perturbation
    /** Use a fraction of the typical length scale of juncitons 
     */
    double _jiggle_intensity;

    /// The tolerance value for equilibrium at jiggling perturbations
    /** After the last jiggling perturbation equilibration is performed with
     *  PCPTopology::_equilibration_tolerance, but equilibrations after jiggling
     *  might be performed at lower precision
     */
    double _jiggle_equilibration_tolerance;

    /// The frequency of cell divisions
    double _cell_divisions_per_step;

    std::pair<double, double> _tissue_stretch_speed;

    /// A model where the parameters are changed over time
    /** See PCPTopology::EnvParam for available parameter
     */
    EnvModel _envm;

    /// The parameter of preferential area per cell type
    /** This parameter is updated in PCPTopology::_envm
     */
    arma::Col<double>::fixed<CellType::num_cell_types> _area_preferential;

    /// Whether to fix hair cell volume
    /** \note Hair cell volume and hair cell apical area not necessarily 
     *        behave the same way. Here volume is kept const.
     */
    bool _fix_hair_cell_volume;
    
    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::uniform_real_distribution<double> _prob_distr;

    // .. Temporary objects ...................................................
    bool _equilibrated;

public:
    // -- Model Setup ---------------------------------------------------------
    /// Construct the PCPTopology model
    /** \param name     Name of this model instance
     *  \param parent   The parent model this model instance resides in
     */
    template<class ParentModel, typename... Taskargs>
    PCPTopology (const std::string name, ParentModel& parent, 
                 Taskargs&&... taskargs)
    :
        // Initialize first via base model
        Base(name, parent, std::forward<Taskargs>(taskargs)...),
        
        // construct the vertex model with an external maximum time stamp
        _vertex_model("PCPVertex", *this,
                    DataIO::time_energy_adaptor, DataIO::energy_adaptor,
                    DataIO::areaelasticity_adaptor,
                    DataIO::linetension_adaptor,
                    DataIO::contractility_adaptor,
                    DataIO::cell_cell_polarity_adaptor,
                    DataIO::polarity_exclusion_adaptor,
                    DataIO::lagrange_net_polarisation_adaptor,
                    DataIO::lagrange_const_concentration_adaptor,
                    DataIO::vertex_position_adaptor,  
                    DataIO::cell_position_adaptor<periodic_bc>,
                    DataIO::edge_link_adaptor),
        
        // the parameter
        _equilibration_tolerance(0.),
        _num_equilibration_steps(0),
        _max_equilibration_iterations(0),
        _num_jiggle_per_equilibration(0),
        _jiggle_intensity(0.),
        _jiggle_equilibration_tolerance(0.),
        _cell_divisions_per_step(get_as<double>("cell_divisions_per_step",
                                                this->_cfg)),
        _tissue_stretch_speed(get_as<std::pair<double, double>>(
                                    "tissue_stretch_speed", this->_cfg)),
        _envm("Environment", *this),
        _area_preferential(),
        _fix_hair_cell_volume(get_as<bool>("fix_hair_cell_volume", this->_cfg)),
        _prob_distr(0.,1.),
        _equilibrated(false)
    {
        if (not this->_cfg["equilibration"]) {
            throw std::invalid_argument("No cfg entry 'equilibration' "
                "available for model " + this->_name + "!");
        }
        _equilibration_tolerance = get_as<double>("tolerance",
                this->_cfg["equilibration"]);
        _num_equilibration_steps = get_as<int>("num_steps",
                this->_cfg["equilibration"]);
        _max_equilibration_iterations = get_as<int>("num_iterations",
                this->_cfg["equilibration"]);
        _num_jiggle_per_equilibration = get_as<int>("num_jiggle",
                this->_cfg["equilibration"]);
        _jiggle_intensity = get_as<double>("jiggle_intensity",
                this->_cfg["equilibration"], 0.);
        _jiggle_equilibration_tolerance = get_as<double>("jiggle_tolerance",
                this->_cfg["equilibration"], _equilibration_tolerance);
        
        if (_jiggle_equilibration_tolerance < _equilibration_tolerance) {
            throw std::invalid_argument("In model " + this->_name + ": "
                "'equilibration' expected jiggle_tolerance to be larger "
                "than tolerance, but " + 
                std::to_string(_jiggle_equilibration_tolerance) + " < " +
                std::to_string(_equilibration_tolerance) + "!");
        }

        _envm.track_parameters({"area_preferential_hair",
                                "area_preferential_support"});
        this->_log->info("Model set up.");
    }


private:
    // .. Setup functions .....................................................

    // .. Helper functions ....................................................
    /// Equilibrates the vertex model
    /** Iterate the vertex model until it reaches an equilibrium state
     *  (see PCPVertex<periodic_bc>::equilibrium_state_reached(double tolerance) const)
     *  with tolerance = _equilibrium_tolerance
     * 
     *  More precisely, the _vertex_model is iterated for 
     *  _num_equilibration_steps, then the equilibrium condition (see above)
     *  is checked:
     * 
     *      If true, iteration is stopped and update_object_containers() called
     * 
     *      If false, _vertex_model is iterated again for _num_equilibration_steps 
     *          this is repeated for a maximum of _max_equilibration_iterations
     */
    void equilibrate_vertex_model() {
        if (_equilibrated) {
            return;
        }

        double tolerance = _jiggle_equilibration_tolerance;
        int time_0 = _vertex_model.get_time();
        bool repeat = false;
        for (int it_jiggle = 0; it_jiggle <= _num_jiggle_per_equilibration;
             it_jiggle++)
        {
            if (it_jiggle == _num_jiggle_per_equilibration) {
                tolerance = _equilibration_tolerance;
            }
            
            auto [Lx, Ly] = _vertex_model.get_domain_size();
            auto num_cells = _vertex_model.get_cells().size();
            double intensity = _jiggle_intensity * sqrt(Lx*Ly / num_cells);
            // NOTE the sqrt(x) defines a typical lengthscale under the 
            //      assumption of isotropic cells
            
            int time_start = _vertex_model.get_time(); 
            this->_log->debug("  Jiggling the vertices on a length scale of "
                "{}. Then equilibrating the vertex model ..", intensity);
            
            _vertex_model.jiggle_vertices(intensity);
            _vertex_model.set_minimisation_precision(tolerance);
            _equilibrated = false;
            while (not _equilibrated) {
                _vertex_model.iterate();

                double energy_change;
                std::tie(_equilibrated,
                    energy_change) = _vertex_model.equilibrium_state_reached();
                
                if (not _equilibrated
                    and _vertex_model.get_time() - time_start >= _num_equilibration_steps)
                {
                    this->_log->warn("ERROR Equilibrium not reached within {} "
                        "steps at a tolerance of {}! Energy change in last "
                        "step was {}.", 
                        _num_equilibration_steps,  _equilibration_tolerance,
                        energy_change);
                    if (not repeat) {
                        it_jiggle--;
                        repeat = true;
                        break;
                    }
                    #ifdef NDEBUG
                    this->_log->warn("Running model in release mode. Some known "
                        "exceptions are only evaluated in debug mode, the author "
                        "recommends to build the model in debug mode!");
                    #endif
                    throw std::runtime_error("Equilibrium not reached!");
                }
                if (not _equilibrated
                    and (_vertex_model.get_time() - time_start) % 10 == 0)
                {
                    _vertex_model.init_minimisation();
                }

                if (stop_now.load()) {
                    this->_log->warn("Was told to stop. Not iterating "
                        "vertex model further ...");
                    throw GotSignal(received_signum.load());
                }
            }
            if (_equilibrated) { repeat = false; }
        }

        this->_log->debug("Vertex model equilibrated within {} steps", 
                        _vertex_model.get_time() - time_0);
        return;
    }

    /// Perform a cell division on a random cell
    /** Divides a random cell into two daughter cells.
     *  The cell is first expanded to the double of its preferential area.
     *  The cell is then divided at an axis through it's center at a random
     *  angle, which creates an edge between the two daughter cells.
     * 
     *  \param threshold    Fraction of the area_preferential at which cell
     *                      is not divided
     */
    void divide_random_cell(double threshold = 0.5) {
        auto cells = _vertex_model.get_cells();
        std::uniform_int_distribution<> int_dist(0, cells.size() - 1);
        
        auto c = cells[int_dist(*this->_rng)].lock();
        auto [Lx, Ly] = _vertex_model.get_domain_size();

        if constexpr (not periodic_bc) {
            for (auto [e, flip] : c->edges_ordered) {
                if (e->adj_cell_a.expired() or e->adj_cell_b.expired()) {
                    this->_log->warn("Cannot divide randomly chosen cell, "
                                      "because it is a boundary cell. Division "
                                      "od boundary cells in non-periodic "
                                      "boundary conditions is not implemented."
                                      "CONTINUING WITHOUT DIVISION.");
                    return;
                }
            }
        }

        _vertex_model.increase_domain_size(c->area_preferential);
        c->area_preferential *= 2;
        _equilibrated = false;

        equilibrate_vertex_model();

        if (c->template area_abs<periodic_bc>(Lx, Ly) < threshold * c->area_preferential) {
            this->_log->warn("Could not divide cell, because it would "
                "not grow to sufficient area. For division requested area: "
                "75\% of {}. Area reached: {}. !!ABORTING!!",
                c->area_preferential, c->template area_abs<periodic_bc>(Lx, Ly));
            throw std::runtime_error("Cell division not possible!");
        }

        c->area_preferential /= 2;
        _vertex_model.divide_cell(c, _prob_distr(*this->_rng) * PI);
        _equilibrated = false;

        equilibrate_vertex_model();
        return;
    }

    /// Perform N cell divisions
    /** \param num_cell_divisions number of cell divisions to be performed
     * 
     */
   void perform_cell_divisions(double num_cell_divisions) {
        int i;
        for (i = 1; i <= num_cell_divisions; ++i) {
            divide_random_cell();
        }
        i -= 1;
        if (num_cell_divisions - i > 0 and 
                _prob_distr(*this->_rng) < num_cell_divisions - i)
        {
            divide_random_cell();
        }
    }
    /// Perform stretch tissue
    /** \param num_cell_divisions number of cell divisions to be performed
     * 
     */
   void stretch_domain (std::pair<double, double> stretch_speed) {
        if (std::get<0>(stretch_speed) == 0 and
            std::get<1>(stretch_speed) == 0)
        {
            return;
        }
        double dA = _vertex_model.stretch_domain(
                            std::get<0>(stretch_speed),
                            std::get<1>(stretch_speed),
                            true, _fix_hair_cell_volume);
        if (_fix_hair_cell_volume) {
            int num_cells = _vertex_model.get_cells().size();
            for (auto c : _vertex_model.get_cells()) {
                num_cells -= c.lock()->type == CellType::hair;
            }
            for (int i = 0; i < CellType::num_cell_types; i++) {
                if (i == CellType::hair) { continue; }                    
                _area_preferential(i) += dA / num_cells;
            }
        }
        else {
            int num_cells = _vertex_model.get_cells().size();
            for (int i = 0; i < CellType::num_cell_types; i++) {                
                _area_preferential(i) += dA / num_cells;
            }                
        }
        _envm.set_parameter("area_preferential_hair",
                            _area_preferential(CellType::hair));
        _envm.set_parameter("area_preferential_support",
                            _area_preferential(CellType::support));

        _equilibrated = false;
        return;
   }

    void differentiate_cells () {
        if (not this->_cfg["differentiation"]
            or not get_as<bool>("active", this->_cfg["differentiation"],
                                     true))
        {

            this->_log->debug("No differentiation requested. Continuing.");

            double area_preferential = get_as<double>("area_preferential",
                        this->_cfg["PCPVertex"]);            
            for (int i = 0; i < CellType::num_cell_types; i++) {
                _area_preferential(i) = area_preferential;
            }

            return;
        }
        this->_log->debug("Preparing differentiating progenitor cells to hair- "
            "and support-cells ...");

        static_assert(CellType::num_cell_types == 3, "Initialisation of "
            "interaction matrices `linetension` and `area_preferential` only "
            "defined for 3 cell types.");
        
        this->_log->debug("Extracting area-preferential (expecting {} "
                          "entries) ..", CellType::num_cell_types);
        if (not this->_cfg["differentiation"]["area_preferential"]) {
            throw std::invalid_argument("Expected cfg dict 'area_preferential' "
                "not available in proliferation!");
        }

        _area_preferential(CellType::progenitor) = get_as<double>(
            "progenitor", this->_cfg["differentiation"]["area_preferential"]);
        _area_preferential(CellType::hair) = get_as<double>(
            "hair", this->_cfg["differentiation"]["area_preferential"]);
        _area_preferential(CellType::support) = get_as<double>(
            "support", this->_cfg["differentiation"]["area_preferential"]);

        this->_log->debug("Extracting linetension (expecting {}! entries, "
                          "i.e. the upper diagonal matrix of a {}x{} matrix) ..",
                          CellType::num_cell_types, CellType::num_cell_types,
                          CellType::num_cell_types);
        if (not this->_cfg["differentiation"]["linetension"]) {
            throw std::invalid_argument("Expected cfg dict 'linetension' "
                "not available in proliferation!");
        }
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> linetension;
        linetension(CellType::progenitor, CellType::progenitor) = get_as<double>(
            "progenitor_progenitor", this->_cfg["differentiation"]["linetension"]);
        linetension(CellType::progenitor, CellType::hair) = get_as<double>(
            "progenitor_hair", this->_cfg["differentiation"]["linetension"]);
        linetension(CellType::progenitor, CellType::support) = get_as<double>(
            "progenitor_support", this->_cfg["differentiation"]["linetension"]);
        linetension(CellType::hair, CellType::hair) = get_as<double>(
            "hair_hair", this->_cfg["differentiation"]["linetension"]);
        linetension(CellType::hair, CellType::support) = get_as<double>(
            "hair_support", this->_cfg["differentiation"]["linetension"]);
        linetension(CellType::support, CellType::support) = get_as<double>(
            "support_support", this->_cfg["differentiation"]["linetension"]);
        for (int i = 0; i < CellType::num_cell_types; i++) {
            for (int j = i+1; j < CellType::num_cell_types; j++) {
                linetension(j, i) = linetension(i, j);
            }
        }
        arma::Mat<double>::fixed<CellType::num_cell_types,
                                 CellType::num_cell_types> contractility;
        contractility(CellType::progenitor, CellType::progenitor) = get_as<double>(
            "progenitor_progenitor", this->_cfg["differentiation"]["contractility"]);
        contractility(CellType::progenitor, CellType::hair) = get_as<double>(
            "progenitor_hair", this->_cfg["differentiation"]["contractility"]);
        contractility(CellType::progenitor, CellType::support) = get_as<double>(
            "progenitor_support", this->_cfg["differentiation"]["contractility"]);
        contractility(CellType::hair, CellType::hair) = get_as<double>(
            "hair_hair", this->_cfg["differentiation"]["contractility"]);
        contractility(CellType::hair, CellType::support) = get_as<double>(
            "hair_support", this->_cfg["differentiation"]["contractility"]);
        contractility(CellType::support, CellType::support) = get_as<double>(
            "support_support", this->_cfg["differentiation"]["contractility"]);
        for (int i = 0; i < CellType::num_cell_types; i++) {
            for (int j = i+1; j < CellType::num_cell_types; j++) {
                contractility(j, i) = contractility(i, j);
            }
        }
        

        auto method = get_as<std::string>("method",
                                          this->_cfg["differentiation"]);
        if (method == "random") {
            double hair_cell_fraction = get_as<double>("hair_cell_fraction",
                                            this->_cfg["differentiation"]);
            _vertex_model.differentiate_hair_cells_random(hair_cell_fraction,
                linetension, contractility, _area_preferential);
        }
        else if (method == "NotchDelta") {
            auto notch_delta = NotchDelta::NotchDelta("NotchDelta", *this,
                    NotchDelta::DataIO::density_time,
                    NotchDelta::DataIO::density_progenitor,
                    NotchDelta::DataIO::density_hair,
                    NotchDelta::DataIO::density_support,
                    NotchDelta::DataIO::density_ratio_hair_support,
                    NotchDelta::DataIO::number_hair_hair_contacts);
            auto steps = get_as<int>("notch_delta_steps",
                                     this->_cfg["differentiation"]);
            _vertex_model.differentiate_hair_cells_NotchDelta(
                std::make_shared<NotchDelta::NotchDelta>(notch_delta), steps,
                linetension, contractility, _area_preferential);
        }
        else {
            throw KeyError("method", this->_cfg["differentiation"], "Method "
                "for differentiation can be: "
                "random, " "NotchDelta");
        }

        auto cells = _vertex_model.get_cells();
        int num_hc = 0;
        for (auto c : cells) {
            if (c.lock()->type == CellType::hair) {
                num_hc++;
            }
        }
        
        _equilibrated = false;
        this->_log->info("Differentiated progenitor cells; "
            "there are {} hair cells out of {} cells ({}%)", num_hc,
            cells.size(), double(num_hc)/_vertex_model.get_cells().size());
        
        return;
    }

    void update_area_preferential () {
        double new_value = _envm.get_parameter("area_preferential_hair");
        if (_area_preferential(CellType::hair) == new_value) {
            return;
        }

        // area increase per hair cell
        _area_preferential(CellType::hair) = new_value;
        double dA = new_value - _area_preferential(CellType::hair);
        
        auto cells = _vertex_model.get_cells();
        int num_hcs = 0;
        for (auto c : cells) {
            num_hcs += (c.lock()->type == CellType::hair);
        }

        dA *= num_hcs; // the total increase of area by all hcs
        if (num_hcs < cells.size()) {
            // area change compensation per non-hair cell
            dA /= cells.size() - num_hcs;
        }
        else { dA = 0.; }

        // update the parameter
        for (int i = 0; i < CellType::num_cell_types; i++) {
            if (i == CellType::hair) { continue; }
            _area_preferential(i) -= dA;
        }

        for (auto c_weak : cells) {
            auto c = c_weak.lock();
            c->area_preferential = _area_preferential(c->type);
        }
        
        this->_log->debug("Updated area preferential from Environment model "
            "to ({}, {}, {}.", _area_preferential(0), _area_preferential(1),
            _area_preferential(2));

        _envm.set_parameter("area_preferential_support",
                            _area_preferential(CellType::support));
        
        _equilibrated = false;
        return;
    }

public:
    // -- Public Interface ----------------------------------------------------
    // .. Simulation Control ..................................................
    /// Iterate a single step
    void perform_step () {
        _envm.iterate();

        // deformations
        this->update_area_preferential();
        this->stretch_domain(_tissue_stretch_speed);
        this->equilibrate_vertex_model();

        this->perform_cell_divisions(_cell_divisions_per_step);
    }

    /// Monitor model information
    void monitor () {        
        this->_monitor.set_entry("num cells", _vertex_model.get_cells().size());
    }

    /// The prolog
    /** Performs the following tasks:
     *      1. call prolog of vertex model
     *      2. equilibrate vertex model
     *      3. default prolog tasks
     */
    void prolog () {
        _vertex_model.prolog();

        if (this->_cfg["proliferation"]
            and get_as<bool>("active", this->_cfg["proliferation"], true))
        {
            auto num_divisions = get_as<int>("num_cell_divisions",
                                             this->_cfg["proliferation"]);
            int emit_interval = get_as<int>("emit_interval",
                                            this->_cfg["proliferation"], 1);
                                            
            this->_log->info("Performing {} consecutive cell divisions ...",
                             num_divisions);

            this->equilibrate_vertex_model();
            for (int i = 0; i < num_divisions; i++) {
                if (i % emit_interval == 0) {
                    this->_log->info("   Performing cell division {} of {} "
                                      "...", i + 1, num_divisions);
                }
                else {
                    this->_log->debug("   Performing cell division {} of {} "
                                      "...", i + 1, num_divisions);
                }
                this->perform_cell_divisions(1.0);
            }

            this->_log->info("Model initialised with proliferated vertex model. "
                             "There are {} cells on equilibrated tissue.",
                             _vertex_model.get_cells().size());
        }
        
        differentiate_cells();
        
        _envm.prolog();
        update_area_preferential();

        this->equilibrate_vertex_model();

        return this->__prolog();
    }

    /// The epilog
    /** Performs the following tasks:
     *      1. (optional) Equilibrate the vertex model with changed noise level
     *      2. default epilog tasks
     */
    void epilog () {
        _vertex_model.epilog();
        _envm.epilog();
        
        auto [Lx, Ly] = _vertex_model.get_domain_size();
        this->_log->info("Domain size is {} x {}.", Lx, Ly);

        if (not this->_cfg["epilog"]) {
            return this->__epilog();
        }

        auto epilog_cfg = this->_cfg["epilog"];
        if (epilog_cfg["set_noise_const"]) {
            _vertex_model.set_noise_const(get_as<double>("set_noise_const",
                                                          epilog_cfg));
        }
        if (epilog_cfg["set_noise_linear"]) {
            _vertex_model.set_noise_linear(get_as<double>("set_noise_linear",
                                                          epilog_cfg));
        }

        int time_start = _vertex_model.get_time();
        int num_steps = get_as<int>("num_epilog_steps", epilog_cfg);

        this->_log->info("Iterating vertex model from time {} to {}", 
                         time_start, time_start + num_steps);

        for (int i = 0; i < num_steps; ++i) {
            _vertex_model.iterate();

            if (stop_now.load()) {
                this->_log->warn("Was told to stop. Not iterating vertex "
                    "model further ...");
                throw GotSignal(received_signum.load());
            }
        }

        return this->__epilog();
    }

    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

    /// Getter for vertices
    std::vector<std::weak_ptr<Vertex>> get_vertices () {
        return _vertex_model.get_vertices();
    }

    /// Getter for edges
    std::vector<std::weak_ptr<Edge>> get_edges () {
        return _vertex_model.get_edges();
    }

    /// Getter for cells
    std::vector<std::weak_ptr<Cell>> get_cells () {
        return _vertex_model.get_cells();
    }

    auto get_domain_size() const {
        return _vertex_model.get_domain_size();
    }

};

} // namespace PCPTopology
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPTOPOLOGY_HH
