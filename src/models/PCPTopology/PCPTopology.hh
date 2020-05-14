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

    /// The type of a cell
    using Cell = typename PCPVertex::Cell;

    /// The types a cell can take (support, hair)
    using CellType = typename PCPVertex::CellType;

    /// The type of a rule function acting on vertices of the agent manager
    using RuleFuncVertex = typename PCPVertex::RuleFuncVertex;
    
    /// The type of a rule function acting on edges of the agent manager
    using RuleFuncEdge = typename PCPVertex::RuleFuncEdge;
    
    /// The type of a rule function acting on cells of the agent manager
    using RuleFuncCell = typename PCPVertex::RuleFuncCell;

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

    /// The parameter for proliferation
    struct ParamsProliferation {
        /// The number of steps to double the cell's area  
        unsigned int num_increases;

        /// The minimal area of a cell to divide
        /** \details The fraction to which a cell has to reach the target area
         *           so that it can be divided
         */
        double threshold;

        Utopia::IndexType generation_max_id;

        ParamsProliferation (const Config& cfg)
        :
            num_increases(get_as<unsigned int>("num_increases", cfg)),
            threshold(get_as<double>("area_threshold", cfg, 0.5)),
            generation_max_id(0) 
        {
            if (num_increases == 0) {
                throw std::invalid_argument("Parameter 'num_increases' in "
                    "proliferation must be larger than 0. It specifies in how "
                    "many steps the area of the cell should be increased.");
            }
        }
    } _params_proliferation;

    /// Config of operation proliferation
    Config _cfg_proliferation;

    /// Config of operation differentaition
    Config _cfg_differentiation;

    /// Config of operation stretch domain
    Config _cfg_stretch_domain;
    
    /// The parameters for increment of cell area
    struct ParamsIncrementArea {
        /// The value of increment for area of hair cells
        double hair;

        /// The value of increment for area of support cells
        double support;
        
        /// The variance of increment for area of hair cells
        /** \details using a normal distribution
         */
        double var_hair;

        /// The variance of increment for area of support cells
        /** \details using a normal distribution
         */
        double var_support;

        /// Whether to compensate area change in support cells
        /** \details If true N_hc * A^0_hc + N_sc * A^0_sc = const
         */
        bool adapt_support;

        /// The name of this operation
        const std::string name;

        /// The original configuration
        const Config cfg;

        ParamsIncrementArea (const Config& cfg)
        :
            hair(get_as<double>("hair", cfg)),
            support(get_as<double>("support", cfg)),
            var_hair(get_as<double>("var_hair", cfg, 0.)),
            var_support(get_as<double>("var_support", cfg, 0.)),
            adapt_support(get_as<bool>("adapt_support", cfg)),
            name("increment_area"),
            cfg(cfg)
        {
            if (adapt_support and support != 0) {
                throw std::invalid_argument(fmt::format(
                    "In cfg {}: If `adapt_support = True`, "
                    "then the increment value for the support cells "
                    "`support` must be zero, but was {}!", name, support));
            }
        }
    } _params_increment_area;
    
    /// The parameters for increment of cell shape-index
    struct ParamsIncrementShapeindex {
        /// The value of increment for shape-index of hair cells
        double hair;

        /// The value of increment for shape-index of support cells
        double support;

        /// The name of this operation
        const std::string name;

        /// The original configuration
        const Config cfg;

        ParamsIncrementShapeindex (const Config& cfg)
        :
            hair(get_as<double>("hair", cfg, 0.)),
            support(get_as<double>("support", cfg, 0.)),
            name("increment_shape_index"),
            cfg(cfg)
        {
            
        }
    } _params_increment_shape_index;
    
    /// The parameters for increment of edge linetension
    struct ParamsIncrementLinetension {
        /// The value of increment for hair-hair junctions
        double hair_hair;

        /// The value of increment for hair-support junctions
        double hair_support;

        /// The value of increment for support-support junctions
        double support_support;

        /// The name of this operation
        const std::string name;

        /// The original configuration
        const Config cfg;

        ParamsIncrementLinetension (const Config& cfg)
        :
            hair_hair(get_as<double>("hair_hair", cfg, 0.)),
            hair_support(get_as<double>("hair_support", cfg, 0.)),
            support_support(get_as<double>("support_support", cfg, 0.)),
            name("increment_linetension"),
            cfg(cfg)
        {

        }
    } _params_increment_linetension;
    
    /// The parameters for increment of edge contractility
    struct ParamsIncrementContractility {
        /// The value of increment for hair-hair junctions
        double hair_hair;

        /// The value of increment for hair-support junctions
        double hair_support;

        /// The value of increment for support-support junctions
        double support_support;

        /// The name of this operation
        const std::string name;

        /// The original configuration
        const Config cfg;

        ParamsIncrementContractility (const Config& cfg)
        :
            hair_hair(get_as<double>("hair_hair", cfg, 0.)),
            hair_support(get_as<double>("hair_support", cfg, 0.)),
            support_support(get_as<double>("support_support", cfg, 0.)),
            name("increment_contractility"),
            cfg(cfg)
        {

        }
    } _params_increment_contractility;
    
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
                    DataIO::vertices_adaptor,  
                    DataIO::cells_adaptor<SpaceVec, CellType>,
                    DataIO::edges_adaptor),
        
        // the parameter
        _equilibration_tolerance(0.),
        _num_equilibration_steps(0),
        _max_equilibration_iterations(0),
        _num_jiggle_per_equilibration(0),
        _jiggle_intensity(0.),
        _jiggle_equilibration_tolerance(0.),
        _params_proliferation(this->_cfg["proliferation"]),
        _params_increment_area(extract_cfg("increment_area", this->_cfg)),
        _params_increment_shape_index(extract_cfg("increment_shape_index",
                                                  this->_cfg)),
        _params_increment_linetension(extract_cfg("increment_linetension",
                                                  this->_cfg)),
        _params_increment_contractility(extract_cfg("increment_contractility",
                                                    this->_cfg)),
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
    /// Extract the configuration of an operation
    Config extract_cfg (std::string name, const Config& cfg) {
        if (not cfg[name]) {
            throw KeyError(name, cfg);
        }

        return cfg[name];
    }

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
                    this->_log->error("ERROR Equilibrium not reached within {} "
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
                    this->_log->error("Running model in release mode. Some known "
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
    void divide_random_cell() {
        const auto& am = _vertex_model.get_am();
        const auto& cells = am.cells();

        Utopia::IndexType generation_max_id = _params_proliferation.generation_max_id;
        AgentContainer<Cell> cells_of_current_generation(cells.size());
        auto it = std::copy_if (cells.begin(), cells.end(),
                                cells_of_current_generation.begin(),
                                [generation_max_id](const auto& cell){
                                    return cell->id() <
                                        generation_max_id; } );
        cells_of_current_generation.resize(
            std::distance(cells_of_current_generation.begin(), it));
        if (cells_of_current_generation.size() == 0) {
            for (const auto& c : cells) {
                generation_max_id = std::max(generation_max_id, c->id());
            }
            _params_proliferation.generation_max_id = generation_max_id + 1;
            cells_of_current_generation = cells;
        }

        std::uniform_int_distribution<> int_dist(
            0, cells_of_current_generation.size() - 1);
        
        auto cell = cells_of_current_generation[int_dist(*this->_rng)];
        if (not this->_space->periodic) {
            for (auto [e, flip] : cell->custom_links().edges) {
                const auto [adj_cell_a, adj_cell_b] = am.adjoints_of(e);
                if (not adj_cell_a or not adj_cell_b) {
                    this->_log->error("Cannot divide randomly chosen cell, "
                                      "because it is a boundary cell. Division "
                                      "od boundary cells in non-periodic "
                                      "boundary conditions is not implemented."
                                      "CONTINUING WITHOUT DIVISION.");
                    return;
                }
            }
        }

        double area_preferential = cell->state.area_preferential;
        double dA = cell->state.area_preferential /
                    _params_proliferation.num_increases;
        for (unsigned int i = 0; i < _params_proliferation.num_increases; i++) {            
            _vertex_model.increase_domain_size(dA);
            cell->state.area_preferential += dA;
            equilibrate_vertex_model();
        }

        if (am.area_of(cell) < _params_proliferation.threshold *
                            cell->state.area_preferential)
        {
            this->_log->error("Could not divide cell, because it would "
                "not grow to sufficient area. For division requested area: "
                "{}\% of {}. Area reached: {}. !!ABORTING!!",
                _params_proliferation.threshold * 100,
                cell->state.area_preferential, am.area_of(cell));
            throw std::runtime_error("Cell division not possible!");
        }

        cell->state.area_preferential = area_preferential;
        _vertex_model.divide_cell(cell, _prob_distr(*this->_rng) * PI);

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
        
        arma::Col<double>::fixed<CellType::num_cell_types> area_change;
        int num_cells = _vertex_model.get_am().cells().size();
        if (fix_hair_cell_volume) {
            for (auto c : _vertex_model.get_am().cells()) {
                num_cells -= (c->state.type == CellType::hair);
            }
        }
        
        for (int i = 0; i < CellType::num_cell_types; i++) {
            if (fix_hair_cell_volume and i == CellType::hair) {
                area_change[i] = 0;
            }
            else {
                area_change[i] = dA / num_cells;
            }
        }

        PCPVertex::RuleFuncCell update = [this, area_change](const auto& cell) {
            auto state = cell->state;
            state.area_preferential -= area_change[state.type];
            return state;
        };

        apply_rule<Update::sync>(update, _vertex_model.get_am().cells());
   }

    /// Differentiate cells
    void differentiate_cells () {
        this->_log->debug("Preparing differentiating progenitor cells to hair- "
            "and support-cells ...");

        auto method = get_as<std::string>("method",
                                          this->_cfg["differentiation"]);
        if (method == "random") {
            double hair_cell_fraction = get_as<double>("hair_cell_fraction",
                                            this->_cfg["differentiation"]);
            _vertex_model.differentiate_hair_cells_random(hair_cell_fraction);
        }
        else if (method == "NotchDelta") {
            auto notch_delta = std::make_shared<NotchDelta::NotchDelta>(
                    "NotchDelta", *this,
                    NotchDelta::DataIO::density_time,
                    NotchDelta::DataIO::density_progenitor,
                    NotchDelta::DataIO::density_hair,
                    NotchDelta::DataIO::density_support,
                    NotchDelta::DataIO::density_ratio_hair_support,
                    NotchDelta::DataIO::number_hair_hair_contacts);
            
            auto steps = get_as<int>("notch_delta_steps",
                                     this->_cfg["differentiation"]);
            _vertex_model.differentiate_hair_cells_NotchDelta(notch_delta,
                                                              steps);
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

    void increment_area () {
        const auto& am = _vertex_model.get_am();

        double hair = _params_increment_area.hair;
        
        double support;
        if (not _params_increment_area.adapt_support) {
            support = _params_increment_area.support;
        }
        else {
            unsigned int num_hcs = 0;
            for (const auto& c : am.cells()) {
                if (c->state.type == CellType::hair) { num_hcs++; }
            }
            if (num_hcs == am.cells().size()) {
                support = 0;
            }
            else {
                support = -1. * hair * num_hcs / (am.cells().size() - num_hcs);
            }
        }

        std::normal_distribution<> dist_hair{hair,
                                             _params_increment_area.var_hair};
        std::normal_distribution<> dist_support{support,
                                            _params_increment_area.var_support};

        RuleFuncCell update = [hair, support, this,
                               dist_hair{std::move(dist_hair)},
                               dist_support{std::move(dist_support)}]
                (const auto& cell) mutable
        {
            auto state = cell->state;
            if (state.type == CellType::support) {
                state.area_preferential += dist_support(*this->_rng);
            }
            else if (state.type == CellType::hair) {
                state.area_preferential += dist_hair(*this->_rng);
            }
            return state;
        };
        
        apply_rule<Update::sync>(update, am.cells());
    }

    void increment_shape_index () {
        double hair = _params_increment_shape_index.hair;
        double support = _params_increment_shape_index.support;

        RuleFuncCell update = [hair, support](const auto& cell) {
            auto state = cell->state;
            if (state.type == CellType::support) {
                state.shape_index_preferential += support;
            }
            else if (state.type == CellType::hair) {
                state.shape_index_preferential += hair;
            }
            return state;
        };

        apply_rule<Update::sync>(update, _vertex_model.get_am().cells());
    }

    void increment_linetension () {
        auto linetension = _vertex_model.get_linetension();

        double hair_hair = _params_increment_linetension.hair_hair;
        double hair_support = _params_increment_linetension.hair_support;
        double support_support = _params_increment_linetension.support_support;

        linetension(CellType::hair, CellType::hair) += hair_hair;
        linetension(CellType::support, CellType::support) += support_support;
        linetension(CellType::hair, CellType::support) += hair_support;
        linetension(CellType::support, CellType::hair) = linetension(
                                                            CellType::hair,
                                                            CellType::support);

        // set the contractility and update the edge properties
        _vertex_model.set_linetension(linetension, true);
    }

    void increment_contractility () {
        auto contractility = _vertex_model.get_edge_contractility();

        double hair_hair = _params_increment_contractility.hair_hair;
        double hair_support = _params_increment_contractility.hair_support;
        double support_support = _params_increment_contractility.support_support;

        contractility(CellType::hair, CellType::hair) += hair_hair;
        contractility(CellType::support, CellType::support) += support_support;
        contractility(CellType::hair, CellType::support) += hair_support;
        contractility(CellType::support, CellType::hair) = contractility(
                                                            CellType::hair,
                                                            CellType::support);

        // set the contractility and update the edge properties
        _vertex_model.set_edge_contractility(contractility, true);
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
    /** \param operation    The operation to perform
     *  \param name         The name of the operation
     *  \param cfg          The config from which to extract the times of
     *                      application.
     *  \param prolog       (optional) Whether called during prolog
     *  \param epilog       (optional) Whether called during epilog
     * 
     *  The cfg can have entry `active`: bool, and a `times`: dict:
     *      * `begin` (uint, default: 0): the first step of application
     *      * `end` (uint, default: max_steps): the last step of
     *        application
     *      * `iterates` (uint, default: 1): The number of applications per step
     *      * `probability` (double, [0, 1], default: 1): The probability that 
     *         the operations is performed in this step. If evaluated true, all
     *         iterates are applied
     *      * `emit_interval` (uint, default: 0): The interval to emit
     *        information on progress in `info` level, otherwise progress in
     *        `debug` level
     * 
     *  \return whether operation was performed
     */
    bool perform_operation(std::function<void()> operation, std::string name,
                           const Utopia::DataIO::Config& cfg, 
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
            
            num_steps = get_as<unsigned int>("iterates", cfg["times"], 1);
        }

        int emit_interval = get_as<unsigned int>("emit_interval",
                                                 cfg["times"], 0);
        
        return perform_operation(operation, name, num_steps, emit_interval);
    }

    // .. Simulation Control ..................................................
    /// Iterate a single step
    void perform_step () {
        // increment entity parameter
        perform_operation(
            [this] () { return this->increment_area(); },
            "increment_area", _params_increment_area.cfg);
        perform_operation(
            [this] () { return this->increment_shape_index(); },
            "increment_shape_index", _params_increment_shape_index.cfg);
        perform_operation(
            [this] () { return this->increment_linetension(); },
            "increment_linetension", _params_increment_linetension.cfg);
        perform_operation(
            [this] () { return this->increment_contractility(); },
            "increment_contractility", _params_increment_contractility.cfg);

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

        // increment entity parameter
        perform_operation(
            [this] () { return this->increment_area(); },
            "increment_area", _params_increment_area.cfg,
            true);
        perform_operation(
            [this] () { return this->increment_shape_index(); },
            "increment_shape_index", _params_increment_shape_index.cfg,
            true);
        perform_operation(
            [this] () { return this->increment_linetension(); },
            "increment_linetension", _params_increment_linetension.cfg,
            true);
        perform_operation(
            [this] () { return this->increment_contractility(); },
            "increment_contractility", _params_increment_contractility.cfg,
            true);
        
        return this->__prolog();
    }

    /// The epilog
    /** Performs the following tasks:
     *      1. (optional) Equilibrate the vertex model with changed noise level
     *      2. default epilog tasks
     */
    void epilog () {
        _vertex_model.epilog();
        
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
