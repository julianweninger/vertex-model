#ifndef UTOPIA_MODELS_PCPTOPOLOGY_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_HH

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>

#include "../PCPVertex/PCPVertex.hh"
#include "../PCPVertex/space.hh"
#include "../PCPVertex/energy.hh"
#include "../PCPVertex/algorithm.hh"
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
using PCPTopologyModelTypes = Utopia::ModelTypes<DefaultRNG, WriteMode::managed,
                                                 Space::CustomSpace<2>>;

// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// The PCPTopology Model; the bare-basics a model needs
class PCPTopology:
    public Model<PCPTopology, PCPTopologyModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<PCPTopology, PCPTopologyModelTypes>;

    /// Type of the config
    using typename Base::Config;

    /// The types of a cell (support, hair)
    using CellType = typename PCPVertex::CellType;

    /// The type of coordinates and vectors in space
    using SpaceVec = typename PCPVertex::SpaceVec;
                                    

private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _space
    // ... but you should definitely check out the documentation ;)

    // -- Members -------------------------------------------------------------
    /// The Vertex model
    PCPVertex _vertex_model;

    /// A tolerance value for equilibrium
    double _equilibration_tolerance;

    /// Number of steps performed in VertexModel per iteration
    unsigned int _num_equilibration_steps;

    /// Number of max iterations performed in VertexModel before aborting
    unsigned int _max_equilibration_iterations;

    /// How often to shake the system during equilibration
    /** Equilibration processes might get stuck in local minima, hence the 
     *  equilibrated cellular arrangement is perturbed and equilibration is
     *  repeated.
     */
    unsigned int _num_jiggle_per_equilibration;
    
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

    /// A model where the parameters are changed over time
    /** See PCPTopology::EnvParam for available parameter
     */
    EnvModel _envm;

    /// The parameter of preferential area per cell type
    /** This parameter is updated in PCPTopology::_envm
     */
    arma::Col<double>::fixed<CellType::num_cell_types> _area_preferential;
    
    /// Config of operation proliferation
    Config _cfg_proliferation;

    /// Config of operation differentaition
    Config _cfg_differentiation;

    /// Config of operation stretch domain
    Config _cfg_stretch_domain;
    
    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::uniform_real_distribution<double> _prob_distr;

    // .. Temporary objects ...................................................

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
                    DataIO::cell_position_adaptor<SpaceVec>,
                    DataIO::edge_link_adaptor),
        
        // the parameter
        _equilibration_tolerance(0.),
        _num_equilibration_steps(0),
        _max_equilibration_iterations(0),
        _num_jiggle_per_equilibration(0),
        _jiggle_intensity(0.),
        _jiggle_equilibration_tolerance(0.),
        _envm("Environment", *this),
        _area_preferential(),
        _prob_distr(0.,1.)
    {
        this->_space = _vertex_model.get_space();

        if (not this->_cfg["equilibration"]) {
            throw std::invalid_argument("No cfg entry 'equilibration' "
                "available for model " + this->_name + "!");
        }
        _equilibration_tolerance = get_as<double>("tolerance",
                this->_cfg["equilibration"]);
        _num_equilibration_steps = get_as<unsigned int>("num_steps",
                this->_cfg["equilibration"]);
        _max_equilibration_iterations = get_as<unsigned int>("num_iterations",
                this->_cfg["equilibration"]);
        _num_jiggle_per_equilibration = get_as<unsigned int>("num_jiggle",
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
                                
        double area_preferential = get_as<double>("area_preferential",
                                        this->_cfg["PCPVertex"]["agent_manager"]
                                        ["cell_manager"]["agent_params"]);            
        for (int i = 0; i < CellType::num_cell_types; i++) {
            _area_preferential(i) = area_preferential;
        }

        // copy operation configs
        if (not this->_cfg["proliferation"]) {
            throw KeyError("proliferation", this->_cfg);
        }
        _cfg_proliferation = this->_cfg["proliferation"];

        if (not this->_cfg["differentiation"]) {
            throw KeyError("differentiation", this->_cfg);
        }
        _cfg_differentiation = this->_cfg["differentiation"];
                
        if (not this->_cfg["stretch_domain"]) {
            throw KeyError("stretch_domain", this->_cfg);
        }
        _cfg_stretch_domain = this->_cfg["stretch_domain"];

        this->_log->info("Model set up.");
    }


private:
    // .. Setup functions .....................................................

    // .. Helper functions ....................................................
    /// Equilibrates the vertex model
    /** Iterate the vertex model until it reaches an equilibrium state
     *  (see PCPVertex::equilibrium_state_reached(double tolerance) const)
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
        double tolerance = _jiggle_equilibration_tolerance;
        auto time_0 = _vertex_model.get_time();
        bool repeat = false;
        for (unsigned int it_jiggle = 0;
             it_jiggle <= _num_jiggle_per_equilibration; it_jiggle++)
        {
            if (it_jiggle == _num_jiggle_per_equilibration) {
                tolerance = _equilibration_tolerance;
            }
            
            auto num_cells = _vertex_model.get_am().cells().size();
            const auto domain = _vertex_model.get_space()->get_domain_size();
            double intensity = _jiggle_intensity * sqrt(domain[0]*domain[1] / 
                                                        num_cells);
            // NOTE the sqrt(x) defines a typical lengthscale under the 
            //      assumption of isotropic cells
            
            int time_start = _vertex_model.get_time();
            
            _vertex_model.jiggle_vertices(intensity);
            _vertex_model.set_minimisation_precision(tolerance);
            bool equilibrated = false;
            while (not equilibrated) {
                _vertex_model.iterate();

                double energy_change;
                std::tie(equilibrated,
                    energy_change) = _vertex_model.equilibrium_state_reached();
                
                if (not equilibrated
                    and _vertex_model.get_time() - time_start >= 
                            _num_equilibration_steps)
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
                if (not equilibrated
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
            if (equilibrated) { repeat = false; }
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
        const auto& am = _vertex_model.get_am();
        const auto& cells = am.cells();
        std::uniform_int_distribution<> int_dist(0, cells.size() - 1);
        
        auto c = cells[int_dist(*this->_rng)];
        const auto& domain = _vertex_model.get_space()->get_domain_size();

        if (not this->_space->periodic) {
            for (auto [e, flip] : c->custom_links().edges) {
                const auto [adj_cell_a, adj_cell_b] = am.adjoints_of(e);
                if (not adj_cell_a or not adj_cell_b) {
                    this->_log->warn("Cannot divide randomly chosen cell, "
                                      "because it is a boundary cell. Division "
                                      "od boundary cells in non-periodic "
                                      "boundary conditions is not implemented."
                                      "CONTINUING WITHOUT DIVISION.");
                    return;
                }
            }
        }

        _vertex_model.increase_domain_size(c->state.area_preferential);
        c->state.area_preferential *= 2;

        equilibrate_vertex_model();

        if (am.area_of(c) < threshold * c->state.area_preferential) {
            this->_log->warn("Could not divide cell, because it would "
                "not grow to sufficient area. For division requested area: "
                "75\% of {}. Area reached: {}. !!ABORTING!!",
                c->state.area_preferential, am.area_of(c));
            throw std::runtime_error("Cell division not possible!");
        }

        c->state.area_preferential /= 2;
        _vertex_model.divide_cell(c, _prob_distr(*this->_rng) * PI);

        return;
    }
    
    /// Perform stretch domain
    void stretch_domain () {
        auto stretch_speed = get_as_SpaceVec<2>(
                "stretch_speed", this->_cfg["stretch_domain"]);
        
        if (stretch_speed[0] == 0 and stretch_speed[1] == 0) {
            return;
        }
        auto fix_hair_cell_volume = get_as<bool>("fix_hair_cell_volume",
            this->_cfg["stretch_domain"], false);
        double dA = _vertex_model.stretch_domain(stretch_speed, true,
                                                 fix_hair_cell_volume);
        
        if (fix_hair_cell_volume) {
            int num_cells = _vertex_model.get_am().cells().size();
            for (auto c : _vertex_model.get_am().cells()) {
                num_cells -= (c->state.type == CellType::hair);
            }
            for (int i = 0; i < CellType::num_cell_types; i++) {
                if (i == CellType::hair) { continue; }                    
                _area_preferential(i) += dA / num_cells;
            }
        }
        else {
            int num_cells = _vertex_model.get_am().cells().size();
            for (int i = 0; i < CellType::num_cell_types; i++) {                
                _area_preferential(i) += dA / num_cells;
            }                
        }
        _envm.set_parameter("area_preferential_hair",
                            _area_preferential(CellType::hair));
        _envm.set_parameter("area_preferential_support",
                            _area_preferential(CellType::support));
        return;
   }

    /// Differentiate cells
    void differentiate_cells () {
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

        // initialise area preferential from Vertex model
        double area_preferential = get_as<double>("area_preferential",
                                    this->_cfg["PCPVertex"]["agent_manager"]
                                            ["cell_manager"]["agent_params"]);          
        for (int i = 0; i < CellType::num_cell_types; i++) {
            _area_preferential(i) = area_preferential;
        }

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

        int num_hc = 0;
        for (const auto& c : _vertex_model.get_am().cells()) {
            if (c->state.type == CellType::hair) {
                num_hc++;
            }
        }
        
        this->_log->info("Differentiated progenitor cells; "
            "there are {} hair cells out of {} cells ({}%)", num_hc,
            _vertex_model.get_am().cells().size(),
            double(num_hc)/_vertex_model.get_am().cells().size());
        
        return;
    }

    /// Update the parameters for area preferential from environment model
    bool update_area_preferential () {
        double new_value = _envm.get_parameter("area_preferential_hair");
        if (_area_preferential(CellType::hair) == new_value) {
            return false;
        }

        // area increase per hair cell
        _area_preferential(CellType::hair) = new_value;
        double dA = new_value - _area_preferential(CellType::hair);
        
        const auto& cells = _vertex_model.get_am().cells();
        unsigned int num_hcs = 0;
        for (const auto& c : cells) {
            num_hcs += (c->state.type == CellType::hair);
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

        for (auto& c : cells) {
            c->state.area_preferential = _area_preferential(c->state.type);
        }
        
        this->_log->debug("Updated area preferential from Environment model "
            "to ({}, {}, {}.", _area_preferential(0), _area_preferential(1),
            _area_preferential(2));

        _envm.set_parameter("area_preferential_support",
                            _area_preferential(CellType::support));
        
        return true;
    }

public:
    // -- Public Interface ----------------------------------------------------
    /// Perform an operation
    /** \param operation    The operation to perform
     *  \param name         The name of the operation
     *  \param iterates     How often to apply the operation
     *  \param emit_interval    How often to emit information on the progress
     *                          If 0, no emit
     * 
     *  After every iteration the vertex model is equilibrated
     */
    bool perform_operation(std::function<void()> operation, std::string name,
                           int iterates, int emit_interval)
    {
        if (iterates == 0) {
            return false;
        }

        if (emit_interval > 0) {
            this->_log->info("Performing operation '{}' with {} iterates ..",
                             name, iterates);
        }
        else {
            this->_log->debug("Performing operation '{}' with {} iterates ..",
                              name, iterates);
        }

        if (emit_interval == 0) { emit_interval = iterates + 1; }

        for (int i = 0; i < iterates; i++) {
            operation();            
            this->equilibrate_vertex_model();

            if ((i+1) % emit_interval == 0) {
                this->_log->info("   Performed iterate {} of {} on operation "
                                 "'{}'.", i+1, iterates, name);
            }
            else {
                this->_log->debug("   Performed iterate {} of {} on operation "
                                  "'{}'.", i+1, iterates, name);
            }
        }
        return true;
    }

    /// Decider whether to perform an operation given a config
    /** Decides whether to perform the information dependent on the config at 
     *  this->_cfg[name].
     * 
     *  \param operation    The operation to perform
     *  \param name         (optional) The name of the operation
     *  \param prolog       (optional) Whether called during prolog
     *  \param epilog       (optional) Whether called during epilog
     * 
     *  \return whether operation performed
     */
    bool perform_operation(std::function<void()> operation, std::string name,
                           Utopia::DataIO::Config cfg, 
                           bool prolog=false, bool epilog=false)
    {
        int num_steps;
        if (not get_as<bool>("active", cfg, true)) {
           return false;
        }
        else if (not cfg["times"] and not prolog and not epilog) { 
            num_steps = 1;
        }
        else if (prolog) {
            num_steps = get_as<int>("prolog", cfg["times"], 0);
        }
        else if (epilog) {
            num_steps = get_as<int>("epilog", cfg["times"], 0);
        }
        else {
            auto time = this->get_time();
            if (get_as<unsigned int>("begin", cfg["times"], 0) > time or 
                get_as<unsigned int>("end", cfg["times"],
                                     this->get_time_max()) < time)
            {
                return false;
            }
            double probability = get_as<double>("probability", cfg["times"], 1);
            if (probability != 1 and _prob_distr(*this->_rng) > probability) {
                return false;
            }
            
            num_steps = get_as<int>("iterates", cfg["times"], 1);
        }

        int emit_interval = get_as<int>("emit_interval", cfg, 0);
        
        return perform_operation(operation, name, num_steps, emit_interval);
    }

    // .. Simulation Control ..................................................
    /// Iterate a single step
    void perform_step () {
        this->_envm.iterate();
        bool update = this->update_area_preferential();
        if (update) {
            perform_operation([] () {}, "hair cell growth", 1, 0);
        }

        // deformations
        perform_operation(
            [this] () { return this->divide_random_cell(); },
            "proliferation", _cfg_proliferation);
        perform_operation(
            [this] () { return this->differentiate_cells(); },
            "differentiation", _cfg_differentiation);

        perform_operation(
            [this] () { return this->stretch_domain(); },
            "stretch domain", _cfg_stretch_domain);
    }

    /// Monitor model information
    void monitor () {        
        this->_monitor.set_entry("num cells",
                                 _vertex_model.get_am().cells().size());
    }

    /// The prolog
    /** Performs the following tasks:
     *      1. call prolog of vertex model
     *      2. equilibrate vertex model
     *      3. default prolog tasks
     */
    void prolog () {
        perform_operation(
            [this] () { return _vertex_model.prolog(); },
            "initialise cells", 1, 0);

        perform_operation(
            [this] () { return this->divide_random_cell(); },
            "proliferation", _cfg_proliferation, true);
        this->_log->debug("Model initialised with proliferated vertex model. "
                          "There are {} cells on equilibrated tissue.",
                          _vertex_model.get_am().cells().size());
        
        perform_operation(
            [this] () { return this->differentiate_cells(); },
            "differentiation", _cfg_differentiation, true);            
        
        _envm.track_parameters({"area_preferential_hair",
                                "area_preferential_support"});
        this->_envm.prolog();
        bool update = this->update_area_preferential();
        if (update) {
            perform_operation([] () {}, "hair cell growth", 1, 0);
        }
        
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
        
        const auto domain = _vertex_model.get_space()->get_domain_size();
        this->_log->info("Domain size is {} x {}.", domain[0], domain[1]);

        if (not this->_cfg["epilog"]) {
            return this->__epilog();
        }

        auto epilog_cfg = this->_cfg["epilog"];

        auto num_steps = get_as<int>("num_equilibrations", epilog_cfg);
        this->_log->info("Equilibrating vertex model another {} times ..", 
                         num_steps);

        for (int i = 0; i < num_steps; ++i) {
            this->equilibrate_vertex_model();
        }

        return this->__epilog();
    }

    // Getters and setters ....................................................
    // Add getters and setters here to interface with other model

    /// Getter for vertices
    const auto& get_am () const {
        return _vertex_model.get_am();
    }
};

} // namespace PCPTopology
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPTOPOLOGY_HH
