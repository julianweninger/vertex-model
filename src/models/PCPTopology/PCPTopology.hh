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
#include "../PCPVertex/transitions.hh"
#include "../PCPVertex/PCPVertex_write_tasks.hh"
#include "PCPTopology_write_tasks.hh"


namespace Utopia {
namespace Models {
namespace PCPVertex {
// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

    /// The frequency of cell divisions
    double _cell_divisions_per_step;

    std::pair<double, double> _tissue_stretch_speed;
    
    /// A [0,1]-range uniform distribution used for evaluating probabilities
    std::uniform_real_distribution<double> _prob_distr;

    // .. Temporary objects ...................................................


    // .. Datasets ............................................................
    // NOTE They should be named '_dset_<name>', where <name> is the
    //      dataset's actual name as set in its constructor. Ideally, do not
    //      hide them inside a struct ...
    // std::shared_ptr<DataSet> _dset_my_var;


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

        _num_equilibration_steps(get_as<int>("num_equilibration_steps",
                                             this->_cfg)),
        _max_equilibration_iterations(get_as<int>("max_equilibration_iterations",
                                                  this->_cfg)),
        
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

        _equilibration_tolerance(get_as<double>("equilibration_tolerance",
                                                this->_cfg)),
        _cell_divisions_per_step(get_as<double>("cell_divisions_per_step",
                                                this->_cfg)),
        _tissue_stretch_speed(get_as<std::pair<double, double>>(
                                    "tissue_stretch_speed", this->_cfg)),
        _prob_distr(0.,1.)
    {
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
        bool equilibrated = false;
        int time_start = _vertex_model.get_time();
        while (not equilibrated) {
            this->_log->debug("Iterating vertex model for equilibration from "
                "time {} to {}", _vertex_model.get_time(),
                _vertex_model.get_time() + _num_equilibration_steps);
            for (int i = 0; i < _num_equilibration_steps; ++i) {
                _vertex_model.iterate();

                if (stop_now.load()) {
                    this->_log->warn("Was told to stop. Not iterating vertex "
                        "model further ...");
                    throw GotSignal(received_signum.load());
                }
            }
            equilibrated = _vertex_model.equilibrium_state_reached(
                                            _equilibration_tolerance);
            int max_steps = _num_equilibration_steps * _max_equilibration_iterations;
            if (_vertex_model.get_time() - time_start >= max_steps and 
                not equilibrated)
            {
                throw std::runtime_error("Equilibrium not reached within " +
                        std::to_string(_max_equilibration_iterations) + 
                        " iterations of " + 
                        std::to_string(_num_equilibration_steps) + " steps each "
                        "at a tolerance of " +
                        std::to_string(_equilibration_tolerance) + "! "
                        "The relative change in energy in last step was " +
                        std::to_string(_vertex_model.get_rel_energy_change()) + 
                        ".");
            }
        }

        this->_log->debug("Vertex model equilibrated within {} steps", 
                         _vertex_model.get_time() - time_start);
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

        _vertex_model.increase_domain_size(c->area_preferential);
        c->area_preferential *= 2;

        equilibrate_vertex_model();

        if (c->area_abs(Lx, Ly) < threshold * c->area_preferential) {
            this->_log->warn("Could not divide cell, because it would "
                "not grow to sufficient area. For division requested area: "
                "75\% of {}. Area reached: {}. !!ABORTING!!",
                c->area_preferential, c->area_abs(Lx, Ly));
            throw std::runtime_error("Cell division not possible!");
        }

        c->area_preferential /= 2;

        _vertex_model.divide_cell(c, _prob_distr(*this->_rng) * PI);

        equilibrate_vertex_model();
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

public:
    // -- Public Interface ----------------------------------------------------
    // .. Simulation Control ..................................................
    /// Iterate a single step
    void perform_step () {
        this->perform_cell_divisions(_cell_divisions_per_step);

        if (std::get<0>(_tissue_stretch_speed) != 0 or 
            std::get<1>(_tissue_stretch_speed) != 0)
        {
            _vertex_model.stretch_domain(std::get<0>(_tissue_stretch_speed),
                                         std::get<1>(_tissue_stretch_speed),
                                         true);
            this->equilibrate_vertex_model();
        }
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

        if (not this->_cfg["proliferation"]) {
            this->_log->info("Equilibrating vertex model ...");        
            this->equilibrate_vertex_model();

            this->_log->info("Model initialised with equilibrated vertex "
                             "model.");
        }
        else {
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

        return this->__prolog();
    }

    /// The epilog
    /** Performs the following tasks:
     *      1. (optional) Equilibrate the vertex model with changed noise level
     *      2. default epilog tasks
     */
    void epilog () {
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

    /// Histogram with number of cells with a specific number of neighbours
    /** bins are range(1, 10), the number of neighbouring cells, where
     *  - last bin for 9 or more neighbours
     */
    std::array<int, 9> get_cell_neighbourhood_histogram () {
        auto cells = _vertex_model.get_cells();
        std::array<int, 9> histogram = {0};
        for (auto c : cells) {
            int neighbours = c.lock()->edges_ordered.size();
            neighbours = std::min(neighbours, 9);
            histogram[neighbours - 1]++;
        }

        return histogram;
    }
    
    /// Average cell size for cells with specific number of neighbours
    /** bins are range(0, num_bins), the number of neighbouring cells, where
     *  - first bin for global average 
     *  - last bin for 9 or more neighbours
     */
    std::array<double, 10> get_cell_size_average () {
        auto cells = _vertex_model.get_cells();
        std::array<int, 10> histogram = {0};
        std::array<double, 10> area = {0};
        histogram[0] = cells.size();
        for (auto c_weak : cells) {
            auto c = c_weak.lock();
            int neighbours = c->edges_ordered.size();
            neighbours = std::min(neighbours, 9);
            histogram[neighbours]++;

            area[neighbours] += c->template cell_area<periodic_bc>();
        }
        for (int i = 1; i < 10; i++) {
            if (histogram[i] == 0) { continue; }
            area[0] += area[i];
            area[i] /= double(histogram[i]);
        }
        area[0] /= double(histogram[0]);

        return area;
    }

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

};

} // namespace PCPTopology
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_PCPTOPOLOGY_HH
