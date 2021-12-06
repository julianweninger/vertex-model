#ifndef UTOPIA_MODELS_PLANARCELLPOLARITY_HH
#define UTOPIA_MODELS_PLANARCELLPOLARITY_HH
// TODO Adjust above include guard (and at bottom of file)

// standard library includes
#include <random>

// third-party library includes

// Utopia-related includes
#include <utopia/core/types.hh>
#include <utopia/core/model.hh>
#include <utopia/core/agent_manager.hh>
#include <utopia/core/apply.hh>

#include "../PCPVertex/entities.hh"
#include "../PCPVertex/entities_manager.hh"
#include "../PCPVertex/space.hh"


namespace Utopia::Models::PlanarCellPolarity {

// ++ Type definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/// Type helper to define types used by the model
using ModelTypes = Utopia::ModelTypes<DefaultRNG, WriteMode::managed>;



// ++ Model definition ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// The PlanarCellPolarity Model; a good start for a CA-based model
/** TODO Add your model description here.
 *  This model's only right to exist is to be a template for new models.
 *  That means its functionality is based on nonsense but it shows how
 *  actually useful functionality could be implemented.
 */
template<typename CellManager>
class PlanarCellPolarity:
    public Model<PlanarCellPolarity<CellManager>, ModelTypes>
{
public:
    /// The type of the Model base class of this derived class
    using Base = Model<PlanarCellPolarity<CellManager>, ModelTypes>;

    /// The type of Time
    using Time = typename Base::Time;

    /// The type of a configuration
    using Config = typename Base::Config;
    
    /// Type of an edge
    using Edge = typename CellManager::Edge;

    /// Type of a cell
    using Cell = typename CellManager::Cell;

    /// The state of a cell
    using CellState = typename Cell::State;

    /// Extract the type of the rule function from the CellManager
    using RuleFuncEdge = typename CellManager::RuleFuncEdge;

    /// Extract the type of the rule function from the CellManager
    using RuleFuncCell = typename CellManager::RuleFuncCell;

    using SpaceVec = typename CellManager::Space::SpaceVec;

    using ProteinVec = typename arma::Col<double>;


private:
    // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor, _log, _space
    // ... but you should definitely check out the documentation ;)


    // -- Members -------------------------------------------------------------
    /// The cell manager
    const CellManager& _cm;

    /** The type of boundary condition
     *  Possible choices:
     *      - Dirichlet: On the outside boundary of the tissue the polarity
     *              proteins have value zero.
     *      - Neumann: The polarity proteins on the outside of a boundary bond
     *              have value -sigma, where sigma is the protein level on the
     *              inside of the same bond.
     */
    enum BoundaryType {
        Dirichlet,  // sigma_b = 0 on boundary edges
        Neumann     // sigma_b = -sigma_a on boundary edges
    } _boundary;

    /// The polarity proteins
    /** Every edge is assigned 2 polarity proteins, on left and right side.
     *  NOTE use [a, b] = _cm.template adjoints_of<true>(edge) such that a and b
     *       are the cells corresponding to the proteins sigma_a and sigma_b
     *       respectively.
     */
    std::map<std::shared_ptr<Edge>,
             std::pair<ProteinVec, ProteinVec>> polarity_proteins;

    /// The polarity proteins
    std::map<std::weak_ptr<Edge>,
             std::pair<ProteinVec, ProteinVec>,
             std::owner_less<std::weak_ptr<Edge>>> protein_fluctuations;

    /// timestep scaling
    double _dt;

    /// Number of proteins per edge
    /** Proteins are equally spaced and fixed wrt edge */
    std::size_t _N;

    /// Interaction parameter of cell-cell polarity interaction
    double _cell_cell_polarity_interaction;

    /// Interaction parameter of cell-internal polarity interaction
    double _cell_polarity_exclusion;

    /// A constant used in penalty method on net cell polarisation
    double _penalty_net_polarisation;

    /// A constant used in penalty method on constant cell protein levels
    double _penalty_const_proteins;

    bool _use_lagrange_multipliers;

    /// Ornstein-Uhlenbeck fluctuations: timescale tau and amplitude A
    std::pair<double, double> _fluctuation;

    /// A re-usable uniform real distribution to evaluate probabilities
    std::uniform_real_distribution<double> _prob_distr;
    
    /// A [0,1]-range normal distribution
    std::normal_distribution<double> _normal_distr;


    // .. Temporary objects ...................................................
    std::map<std::weak_ptr<Edge>, ProteinVec,
             std::owner_less<std::weak_ptr<Edge>>> gradient_sigma_a;
    std::map<std::weak_ptr<Edge>, ProteinVec,
             std::owner_less<std::weak_ptr<Edge>>> gradient_sigma_b;

    std::map<std::weak_ptr<Cell>, double,
             std::owner_less<std::weak_ptr<Cell>>> lagrange_net_polarisation;
    std::map<std::weak_ptr<Cell>, double,
             std::owner_less<std::weak_ptr<Cell>>> lagrange_const_proteins;

    std::size_t cells_id_counter;
    std::size_t edges_id_counter;
    std::size_t vertices_id_counter;

    double _energy;
    double _energy_previous;

public:
    // -- Model Setup ---------------------------------------------------------

    /// Construct the PlanarCellPolarity model
    /** \param name             Name of this model instance; is used to extract
     *                          the configuration from the parent model and
     *                          set up a HDFGroup for this instance
     *  \param parent_model     The parent model this model instance resides in
     *  \param custom_cfg       A custom configuration to use instead of the
     *                          one extracted from the parent model using the
     *                          instance name
     */
    template<class ParentModel, typename... WriterArgs>
    PlanarCellPolarity (
        const std::string& name,
        ParentModel& parent_model,
        const CellManager& cm,
        const DataIO::Config& custom_cfg = {},
        std::tuple<WriterArgs...> &&writer_args = {}
    )
    :
        // Initialize first via base model
        Base(name, parent_model, custom_cfg, writer_args),

        // Now initialize the cell manager
        _cm(cm),

        _boundary(setup_boundary(this->_cfg)),

        polarity_proteins{},

        // Initialize model parameters

        // model dynamics parameter
        _dt(get_as<double>("dt", this->_cfg)),
        _N(get_as<std::size_t>("proteins_per_edge", this->_cfg)),
        _cell_cell_polarity_interaction(get_as<double>(
            "cell_cell_polarity_interaction",this->_cfg)),
        _cell_polarity_exclusion(get_as<double>(
            "cell_polarity_exclusion", this->_cfg)),
        _penalty_net_polarisation(get_as<double>(
            "penalty_net_polarisation", this->_cfg)),
        _penalty_const_proteins(get_as<double>(
            "penalty_const_proteins", this->_cfg)),
        _use_lagrange_multipliers(get_as<bool>(
            "use_lagrange_multipliers", this->_cfg)),
        _fluctuation(
            get_as<double>("fluctuation_tau", this->_cfg),
            get_as<double>("fluctuation", this->_cfg)),
        // when adding more parameters: add also to `update_parameters(cfg)`

        // Initialize the uniform real distribution to range [0., 1.]
        _prob_distr(0., 1.),
        _normal_distr(0.,1.),
        gradient_sigma_a{},
        gradient_sigma_b{},
        lagrange_net_polarisation{},
        lagrange_const_proteins{},
        cells_id_counter(0),
        edges_id_counter(0),
        vertices_id_counter(0),
        _energy(0.),
        _energy_previous(0.)
    {
        // Initialization should be finished here.
        this->_log->debug("{} model fully set up.", this->_name);
    }

    virtual ~PlanarCellPolarity() {}


private:
    // .. Setup functions .....................................................
    /// Extract boundary condition from cfg
    /** Only in non-periodic boundary conditions.
     *  Possible choices:
     *      - Dirichlet: On the outside boundary of the tissue the polarity
     *              proteins have value zero.
     *      - Neumann: The polarity proteins on the outside of a boundary bond
     *              have value -sigma, where sigma is the protein level on the
     *              inside of the same bond.
     */
    BoundaryType setup_boundary(const Config& cfg) const {
        if (_cm.get_space()->periodic) {
            return BoundaryType::Dirichlet;
        }

        auto type = get_as<std::string>("boundary_type", cfg);

        if (type == "Dirichlet") {
            return BoundaryType::Dirichlet;
        }
        else if (type == "Neumann") {
            return BoundaryType::Neumann;
        }
        else {
            throw std::runtime_error(fmt::format("Unknown boundary_type {}! "
                "Please choose one of the following: {}.",
                type, "zero, mirror"));
        }
    }

    void initialize_cells() {
        this->_log->debug("Initializing PCP cells ...");

        for (const auto& edge : _cm.edges()) {
            ProteinVec null = ProteinVec(_N, arma::fill::zeros);
            polarity_proteins[edge] = std::make_pair(null, null);
            protein_fluctuations[edge] = std::make_pair(null, null);
            gradient_sigma_a[edge] = null;
            gradient_sigma_b[edge] = null;
        }
        edges_id_counter = _cm.edges_id_counter();
        for (const auto& cell : _cm.cells()) {
            lagrange_net_polarisation[cell] = 0.;
            lagrange_const_proteins[cell] = 0.;
        }
        cells_id_counter = _cm.cells_id_counter();

        auto method = get_as<std::string>("initialization", this->_cfg);
        if (method == "random") {
            initialize_cells_random();
        }
        else if (method == "aligned") {
            initialize_cells_aligned(
                get_as<double>("alignment_orientation", this->_cfg),
                get_as<double>("alignment_orientation_stddev", this->_cfg)
            );
        }
        else if (method == "aligned_antisymmetric") {
            initialize_cells_aligned_antisymmetric(
                get_as<double>("alignment_orientation", this->_cfg),
                get_as<double>("alignment_orientation_stddev", this->_cfg)
            );
        }
        else if (method == "uniform_angle_distribution") {
            initialize_cells_uniform_angle_distribution();
        }
        else {
            this->_log->error("Unexpected method for cell initialization "
                "(initialization='{}'). Must be one of: "
                    "random, "
                    "aligned.", method);
            throw std::runtime_error("Unexpected method for cell "
                "initialization!");
        }


        this->_log->info("Initialized all cells with constant protein level "
                         "and no net polarization.");
    }

    void initialize_cells_random(double concentration = 1.) {
        this->_log->debug("Initializing random PCP protein levels ..");
        for (const auto& cell : _cm.cells()) {
            const auto& edges = cell->custom_links().edges;
            const auto num_nbs = edges.size();

            // the last rn is dependent
            std::vector<double> rns(num_nbs - 1);
            std::generate(rns.begin(), rns.end(), 
                          [this](){
                              return this->_prob_distr(*this->_rng) - 0.5;
                          });
            
            // the sum of all rns must be 0
            // No net polarity
            rns.push_back(- std::accumulate(rns.begin(), rns.end(), 0.));
            
            // Normalize the square sum, i.e. protein concentration
            double norm_2 = std::accumulate(rns.begin(), rns.end(), 0.,
                                            [](const double& val, double r) {
                                                return val + std::pow(r, 2);
                                            });
            double norm = sqrt(concentration / norm_2);

            // the norm of the rns must be protein_level
            std::transform(rns.begin(), rns.end(), rns.begin(),
                           [norm](auto r) {
                               return r * norm;
                           });

            for (std::size_t it = 0; it < num_nbs; it++) {
                const auto& [edge, flip] = edges[it];
                auto [sigma_a, sigma_b] = polarity_proteins[edge];
                if (not flip) { 
                    sigma_a = rns[it];
                }
                else {
                    sigma_b = rns[it];
                }
                polarity_proteins[edge] = std::make_pair(sigma_a, sigma_b);
            }
        }
    }

    void initialize_cell_with_orientation(const std::shared_ptr<Cell>& cell,
                                          double angle, double concentration)
    {
        const auto space = _cm.get_space();

        SpaceVec axis({cos(angle), sin(angle)});

        SpaceVec cell_center = _cm.barycenter_of(cell);

        const auto& edges = cell->custom_links().edges;
        std::vector<double> orientations{};
        orientations.reserve(edges.size());
        for (const auto& [edge, flip] : edges) {
            SpaceVec a = _cm.position_of(edge->custom_links().a);
            SpaceVec e_vec = _cm.displacement(edge);
            double dN = 1. / static_cast<double>(_N);
            for (std::size_t i = 0; i < _N; i++) {
                SpaceVec pos = a + (i + 0.5) * dN * e_vec;
                SpaceVec orientation = space->displacement(cell_center, pos);
                orientation /= arma::norm(orientation);
                
                orientations.push_back(arma::dot(axis, orientation));
            }
        }

        double net = std::accumulate(orientations.begin(),
                                     orientations.end(), 0.);
        net /= orientations.size();

        for (std::size_t i = 0; i < orientations.size(); i++) {
            orientations[i] -= net;
        }
        // NOTE net polarity is zero
        
        // Normalize the square sum, i.e. protein concentration
        double norm_2 = std::accumulate(orientations.begin(),
                                        orientations.end(),
                                        0.,
                                        [](const double& val, double r) {
                                            return val + std::pow(r, 2);
                                        });
        double norm = sqrt(concentration / norm_2);

        // the norm of the polairty must be protein_level
        std::transform(orientations.begin(), orientations.end(),
                        orientations.begin(),
                        [norm](auto val) {
                            return val * norm;
                        });

        for (std::size_t it = 0; it < edges.size(); it++) {
            const auto& [edge, flip] = edges[it];
            
            ProteinVec s(_N);
            for (std::size_t i = 0; i < _N; i++) {
                s[i] = orientations[it * _N + i];
            }

            ProteinVec sigma_a, sigma_b;
            std::tie(sigma_a, sigma_b) = polarity_proteins[edge];
            if (not flip) {
                sigma_a = s;
            }
            else {
                sigma_b = s;
            }
            polarity_proteins[edge] = std::make_pair(sigma_a, sigma_b);
        }
    }

    /// Initialize cells with aligned polarity
    /** Distribute polarity proteins to result in a given orientation.
     *  Orientation is drawn from normal distribution with mean orientation 
     *  and stddev.
     */
    void initialize_cells_aligned(double orientation, double stddev,
                                  double concentration = 1.)
    {
        this->_log->debug("Initializing aligned PCP protein levels ..");

        std::normal_distribution<double> distr(orientation, stddev);

        for (const auto& cell : _cm.cells()) {
            double angle = distr(*this->_rng);
            initialize_cell_with_orientation(cell, angle, concentration);
        }
    }

    /// Initialize cells with aligned polarity
    /** Distribute polarity proteins to result in a given orientation.
     *  Orientation is drawn from normal distribution with mean orientation 
     *  and stddev.
     *  The global orientations are anti-symmetric about a vertical axis through
     *  the barycenter of the tissue
     */
    void initialize_cells_aligned_antisymmetric(double orientation,
                                                double stddev,
                                                double concentration = 1.)
    {
        this->_log->debug("Initializing aligned PCP protein levels with "
                          "orientational antisymmetry about the vertical "
                          "center.");

        std::normal_distribution<double> distr(orientation, stddev);

        const auto boundary_edges = _cm.get_boundary_edges();
        SpaceVec barycenter = _cm.barycenter_of(boundary_edges);

        for (const auto& cell : _cm.cells()) {
            double angle = distr(*this->_rng);
            SpaceVec pos = _cm.barycenter_of(cell);
            if (pos[0] < barycenter[0]) {
                angle = M_PI - angle;
            }
            initialize_cell_with_orientation(cell, angle, concentration);
        }
    }

    /// Initialize cells with aligned polarity
    /** Distribute polarity proteins to result in a given orientation.
     *  Orientation is drawn from uniform distribution [0, 2*Pi]
     */
    void initialize_cells_uniform_angle_distribution(double concentration = 1.)
    {
        this->_log->debug("Initializing uniformly distributed angles");

        std::uniform_real_distribution<double> distr(0., 2 * M_PI);

        for (const auto& cell : _cm.cells()) {
            double angle = distr(*this->_rng);
            initialize_cell_with_orientation(cell, angle, concentration);
        }
    }
    

    // .. Helper functions ....................................................

    // .. Rule functions ......................................................
    void reset_gradient_polarity() {
        // reset the gradient
        gradient_sigma_a.clear();
        gradient_sigma_b.clear();
        for (const auto& edge : _cm.edges()) {
            gradient_sigma_a[edge] = ProteinVec(_N, arma::fill::zeros);
            gradient_sigma_b[edge] = ProteinVec(_N, arma::fill::zeros);
        }

        // update links of cells
        bool update = false;
        if (cells_id_counter != _cm.cells_id_counter()) {
            update = true;

            // erase old links
            for (auto iter = lagrange_net_polarisation.begin();
                iter != lagrange_net_polarisation.end(); 
                // None
            )
            {
                if ((iter->first).expired()) {
                    iter = lagrange_net_polarisation.erase(iter);
                } 
                else {
                    ++iter;
                }
            }
            for (auto iter = lagrange_const_proteins.begin();
                iter != lagrange_const_proteins.end(); 
                // None
            )
            {
                if ((iter->first).expired()) {
                    iter = lagrange_const_proteins.erase(iter);
                } 
                else {
                    ++iter;
                }
            }

            // add new links
            for (const auto& cell : _cm.cells()) {
                if (   lagrange_net_polarisation.find(cell)
                    == lagrange_net_polarisation.end())
                {
                    lagrange_net_polarisation[cell] = 0.;
                    lagrange_const_proteins[cell] = 0.;
                }
            }
            cells_id_counter = _cm.cells_id_counter();
        }
        // update links of edges
        if (edges_id_counter != _cm.edges_id_counter()) {
            update = true;

            // conserve the polarity proteins of edges that have been split
            // during cell division, i.e. those for which one vertex remains

            // erase old links
            std::map<std::shared_ptr<Edge>, 
                     std::pair<std::pair<ProteinVec, ProteinVec>,
                               std::pair<ProteinVec, ProteinVec>>> removals;
            for (auto iter = polarity_proteins.begin();
                iter != polarity_proteins.end(); 
                // None
            )
            {
                if (iter->first->state.remove) {
                    std::shared_ptr<Edge> edge = iter->first;
                    removals[edge] = std::make_pair(
                        polarity_proteins[edge],
                        protein_fluctuations[edge]
                    );
                    iter = polarity_proteins.erase(iter);
                } 
                else {
                    ++iter;
                }
            }
            for (auto iter = protein_fluctuations.begin();
                iter != protein_fluctuations.end(); 
                // None
            )
            {
                if (   (iter->first).expired()
                    or (iter->first).lock()->state.remove) {
                    iter = protein_fluctuations.erase(iter);
                } 
                else {
                    ++iter;
                }
            }

            for (const auto& edge : _cm.edges()) {
                if (polarity_proteins.find(edge) == polarity_proteins.end()) {
                    auto find = std::find_if(
                        removals.begin(), removals.end(),
                        [edge](const auto& e_sigma){
                            auto a = edge->custom_links().a;
                            auto b = edge->custom_links().b;
                            auto _e = std::get<0>(e_sigma);
                            auto _a = _e->custom_links().a;
                            auto _b = _e->custom_links().b;
                            if (a == _a or b == _b) {
                                return true;
                            }
                            return false;
                        });
                    auto find_flip = std::find_if(
                        removals.begin(), removals.end(),
                        [edge](const auto& e_sigma){
                            auto a = edge->custom_links().a;
                            auto b = edge->custom_links().b;
                            auto _e = std::get<0>(e_sigma);
                            auto _a = _e->custom_links().a;
                            auto _b = _e->custom_links().b;
                            if (a == _b or b == _a) {
                                return true;
                            }
                            return false;
                        });
                    if (find != removals.end()) {
                        polarity_proteins[edge] = std::get<0>(find->second);
                        protein_fluctuations[edge] = std::get<1>(find->second);
                        removals.erase(find);
                    }
                    else if (find_flip != removals.end()) {
                        auto [sa, sb] = std::get<0>(find_flip->second);
                        auto [_sa, _sb] = std::get<1>(find_flip->second);
                        polarity_proteins[edge] = std::make_pair(sb, sa);
                        protein_fluctuations[edge] = std::make_pair(_sb, _sa);
                        removals.erase(find_flip);
                    }
                    else {
                        polarity_proteins[edge] = std::make_pair(0., 0.);
                        protein_fluctuations[edge] = std::make_pair(0., 0.);
                    }
                }
            }

            edges_id_counter = _cm.edges_id_counter();
        }

        if (update) {
            _energy = this->total_energy();
        }
        _energy_previous = _energy;
    }    

    /** The update of polarity protein levels
     * 
     *  Change polarity level proportional to the gradient of energy (force)
     * 
     *  @param e    The pointer to the edge to update
     */
    void update_polarity () {
        for (const auto& edge : _cm.edges()) {
            auto [sigma_a, sigma_b] = polarity_proteins[edge];

            auto [adj_a, adj_b] = _cm.template adjoints_of<true>(edge);
            if (_boundary == BoundaryType::Dirichlet) {
                // exponential decay
                if (adj_a == nullptr) {
                    gradient_sigma_a[edge] = sigma_a;
                }
                if (adj_b == nullptr) {
                    gradient_sigma_b[edge] = sigma_b;
                }
            }

            // steepest gradient
            sigma_a -= gradient_sigma_a[edge] * this->_dt;
            sigma_b -= gradient_sigma_b[edge] * this->_dt;

            if (_boundary == BoundaryType::Neumann) {
                if (not adj_a) {
                    sigma_a = -sigma_b;
                }
                if (not adj_b) {
                    sigma_b = -sigma_a;
                }
            }

            polarity_proteins[edge] = std::make_pair(sigma_a, sigma_b);


            // Ornstein-Uhlenbeck fluctuations
            auto [_sa, _sb] = protein_fluctuations[edge];

            auto [tau, dp] = _fluctuation;

            ProteinVec rand_a(_N), rand_b(_N);

            for (std::size_t i = 0; i < _N; i++) {
                rand_a[i] = (  dp * sqrt(2. * _dt / tau)
                             * _normal_distr(*this->_rng));
                rand_b[i] = (  dp * sqrt(2. * _dt / tau)
                             * _normal_distr(*this->_rng));
            }

            if (adj_a or _boundary == BoundaryType::Neumann) {
                _sa += rand_a - _dt / tau * _sa;
            }
            else {
                _sa -= _dt / tau * _sa;
            }
            if (adj_b or _boundary == BoundaryType::Neumann) {
                _sb += rand_b - _dt / tau * _sb;
            }
            else {
                _sb -= _dt / tau * _sb;
            }

            protein_fluctuations[edge] = std::make_pair(_sa, _sb);
        }
    };

    void update_lagrange_multiplier() {
        if (not _use_lagrange_multipliers) {
            return;
        }

        for (const auto& cell : _cm.cells()) {
            double cum_sigma = 0.;
            double cum_sigma_2 = 0.;
            for (const auto& e_pair : cell->custom_links().edges) {
                ProteinVec sigma = std::get<0>(get_polarity_proteins(e_pair));
                cum_sigma += arma::accu(sigma);
                cum_sigma_2 += std::pow(arma::norm(sigma), 2);
            }   

            lagrange_net_polarisation[cell] += cum_sigma * _dt;
            lagrange_const_proteins[cell] += (cum_sigma_2 - 1.) * _dt;
        }
    }


    void set_gradient_cell_cell_interaction(const AgentContainer<Edge>& edges) {
        for (const auto& edge : edges) {
            const auto [sigma_a, sigma_b] = get_polarity_proteins(edge);

            gradient_sigma_a[edge] += _cell_cell_polarity_interaction * sigma_b;
            gradient_sigma_b[edge] += _cell_cell_polarity_interaction * sigma_a;
        }
    }

    void set_gradient_polarity_exclusion (const AgentContainer<Cell>& cells) {
        for (const auto& cell : cells) {
            const auto& edges = cell->custom_links().edges;
            std::size_t N = edges.size() * _N;

            // create a vector that contains all sigma in correct order
            ProteinVec sigma(N);
            for (std::size_t it = 0; it < edges.size(); it++) {
                auto [edge, flip] = edges[it];
                ProteinVec _s = std::get<0>(get_polarity_proteins(edge, flip));
                for (std::size_t pi = 0; pi < _N; pi++) {
                    sigma[it * _N + pi] = _s[pi];
                }
            }

            // shift in circular manner
            ProteinVec sigma_p1 = arma::shift(sigma,  1);
            ProteinVec sigma_m1 = arma::shift(sigma, -1);

            // elements interacts with neighbors only
            ProteinVec grad = -_cell_polarity_exclusion * (sigma_p1 + sigma_m1);

            // assign gradients
            for (std::size_t it = 0; it < edges.size(); it++) {
                auto [edge, flip] = edges[it];
                ProteinVec _s = std::get<0>(get_polarity_proteins(edge, flip));
                for (std::size_t pi = 0; pi < _N; pi++) {
                    if (not flip) {
                        gradient_sigma_a[edge] += grad[it * _N + pi];
                    }
                    else {
                        gradient_sigma_b[edge] += grad[it * _N + pi];
                    }
                }
            }
        } 
    }


    void set_gradient_lagrange () {
        if (_use_lagrange_multipliers) {
            for (const auto& edge : _cm.edges()) {
                auto [adj_a, adj_b] = _cm.template adjoints_of<true>(edge);
                auto [sigma_a, sigma_b] = get_polarity_proteins(edge);

                gradient_sigma_a[edge] += lagrange_net_polarisation[adj_a];
                gradient_sigma_b[edge] += lagrange_net_polarisation[adj_b];

                gradient_sigma_a[edge] += (  lagrange_const_proteins[adj_a]
                                           * sigma_a * 2.);
                gradient_sigma_b[edge] += (  lagrange_const_proteins[adj_b]
                                           * sigma_b * 2.);
            }
        }

        // penalty method
        for (const auto& cell : _cm.cells()) {
            double cum_sigma = 0.;
            double cum_sigma_2 = 0.;
            for (const auto& e_pair : cell->custom_links().edges) {
                auto sigma = std::get<0>(get_polarity_proteins(e_pair));
                cum_sigma += arma::accu(sigma);
                cum_sigma_2 += std::pow(arma::norm(sigma), 2);
            }
            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto sigma = std::get<0>(get_polarity_proteins(edge, flip));
                if (not flip) {
                    gradient_sigma_a[edge] += (
                        _penalty_net_polarisation * cum_sigma);
                    gradient_sigma_a[edge] += (
                          _penalty_const_proteins
                        * (cum_sigma_2 - 1.) * 2. * sigma);
                }
                else {
                    gradient_sigma_b[edge] += (
                        _penalty_net_polarisation * cum_sigma);
                    gradient_sigma_b[edge] += (
                          _penalty_const_proteins
                        * (cum_sigma_2 - 1.) * 2. * sigma);
                }
            }
        }
    }

public:
    // -- Public Interface ----------------------------------------------------
    // .. Simulation Control ..................................................

    /// Iterate a single step
    void perform_step () {
        reset_gradient_polarity();

        set_gradient_cell_cell_interaction(_cm.edges());
        set_gradient_polarity_exclusion(_cm.cells());

        set_gradient_lagrange();

        update_lagrange_multiplier();
        update_polarity();

        _energy = total_energy();
        
        this->_log->debug("Energy changed by: {}",
                          (_energy - _energy_previous) / fabs(_energy));
    }


    /// Monitor model information
    void monitor () {
        this->_monitor.set_entry("time", this->get_time());
        this->_monitor.set_entry("Energy", _energy);
        this->_monitor.set_entry("dE", 
                                 (_energy - _energy_previous) / fabs(_energy));
    }


    /// Write data
    void write_data () {
        
    }
    
    void prolog () {
        initialize_cells();
        reset_gradient_polarity();

        return this->__prolog();
    }

    void epilog () {
        double cum_sigma_max = 0.;
        double cum_sigma_2_max = 0.;
        for (const auto& cell : _cm.cells()) {
            double cum_sigma = 0.;
            double cum_sigma_2 = 0.;
            for (const auto& e_pair : cell->custom_links().edges) {
                ProteinVec sigma = std::get<0>(get_polarity_proteins(e_pair));
                cum_sigma += arma::accu(sigma);
                cum_sigma_2 += std::pow(arma::norm(sigma), 2);
            }
            cum_sigma_max = std::max(cum_sigma_max, fabs(cum_sigma));
            cum_sigma_2_max = std::max(cum_sigma_2_max, fabs(cum_sigma_2 - 1.));
        }
        
        this->_log->warn("The max difference to no cell polarisation is {}.",
                         cum_sigma_max);
        this->_log->warn("The max difference to cell protein level 1. is {}.",
                         cum_sigma_2_max);

        return this->__epilog();
    }


    // .. Getters and setters .................................................
    // Add getters and setters here to interface with other models
    std::pair<ProteinVec, ProteinVec> get_polarity_proteins
    (const std::shared_ptr<Edge>& edge, bool flip = false) const
    {
        if (edges_id_counter != _cm.edges_id_counter()) {
            if (polarity_proteins.find(edge) == polarity_proteins.end()) {
                ProteinVec null(_N, arma::fill::zeros);
                return std::make_pair(null, null);
            }
        }

        ProteinVec a, b, _a, _b;
        std::tie( a,  b) = polarity_proteins.at(edge);
        std::tie(_a, _b) = protein_fluctuations.at(edge);
        a += _a; b += _b;
        if (flip) {
            std::swap(a, b);
        }
        return std::make_pair(a, b);
    }

    std::pair<ProteinVec, ProteinVec> get_polarity_proteins
    (const std::pair<std::shared_ptr<Edge>, bool>& e_pair) const 
    {
        const auto& [edge, flip] = e_pair;
        return get_polarity_proteins(edge, flip);
    }

    /// Update the model parameters from config
    void update_parameters(const Config& cfg) {
        if (cfg.size() == 0) {
            return;
        }

        this->_log->debug("Updating parameters from config: {}",
                          to_string(cfg));

        _dt = get_as<double>("dt", cfg, _dt);
        _cell_cell_polarity_interaction = get_as<double>(
            "cell_cell_polarity_interaction", cfg, 
            _cell_cell_polarity_interaction);
        _cell_polarity_exclusion = get_as<double>(
            "cell_polarity_exclusion", cfg, 
            _cell_polarity_exclusion);
        _penalty_net_polarisation = get_as<double>(
            "penalty_net_polarisation", cfg, 
            _penalty_net_polarisation);
        _penalty_const_proteins = get_as<double>(
            "penalty_const_proteins", cfg, 
            _penalty_const_proteins);
        _use_lagrange_multipliers = get_as<bool>(
            "use_lagrange_multipliers", cfg, 
            _use_lagrange_multipliers);
        _fluctuation = std::make_pair(
            get_as<double>("fluctuation_tau", cfg,
               std::get<0>(_fluctuation)),
            get_as<double>("fluctuation", cfg,
               std::get<1>(_fluctuation)));
    }
    
    // .. Energy terms ........................................................
protected:
    double energy_cell_cell_polarity (const std::shared_ptr<Edge>& edge) const {
        ProteinVec sigma_a, sigma_b;
        std::tie(sigma_a, sigma_b) = get_polarity_proteins(edge);

        return _cell_cell_polarity_interaction * arma::accu(sigma_a % sigma_b);
    }

    double energy_cell_cell_polarity (const AgentContainer<Edge>& edges) const {
        double energy = 0.;
        for (const auto& edge : edges) {
            energy += energy_cell_cell_polarity(edge);
        }

        return energy;
    }

    double energy_polarity_exclusion(const std::shared_ptr<Cell>& cell) const {
        const auto& edges = cell->custom_links().edges;
        std::size_t N = edges.size() * _N;
    
        // create a vector that contains all sigma in correct order
        ProteinVec sigma(N);
        for (std::size_t it = 0; it < edges.size(); it++) {
            auto [edge, flip] = edges[it];
            ProteinVec _s = std::get<0>(get_polarity_proteins(edge, flip));
            for (std::size_t pi = 0; pi < _N; pi++) {
                sigma[it * _N + pi] = _s[pi];
            }
        }

        // shift in circular manner
        ProteinVec sigma_p1 = arma::shift(sigma,  1);

        return - _cell_polarity_exclusion * arma::accu(sigma % sigma_p1);
    }

    double energy_polarity_exclusion(const AgentContainer<Cell>& cells) const {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += energy_polarity_exclusion(cell);
        }
        return energy;
    }


    double energy_lagrange_net_polarisation
    (const std::shared_ptr<Cell>& cell) const
    {
        if (   lagrange_net_polarisation.find(cell)
            == lagrange_net_polarisation.end())
        {
            return 0.;
        }

        double cum_sigma = 0.;
        for (auto e_pair : cell->custom_links().edges) {
            cum_sigma += arma::accu(std::get<0>(get_polarity_proteins(e_pair)));
        }

        return (  lagrange_net_polarisation.at(cell) * cum_sigma
                + 0.5 * _penalty_net_polarisation * std::pow(cum_sigma, 2));
    }
    /// Getter for energy associated with lagrange multiplier I
    /** Langrange multiplier I is the constrain of zero net polarisation within
     *  a cell.
     */
    double energy_lagrange_net_polarisation
    (const AgentContainer<Cell>& cells) const
    {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += energy_lagrange_net_polarisation(cell);
        }
        return energy;
    }
    
    double energy_lagrange_const_proteins
    (const std::shared_ptr<Cell>& cell) const
    {
        if (   lagrange_const_proteins.find(cell)
            == lagrange_const_proteins.end())
        {
            return 0.;
        }

        double concentration = 0.;
        for (const auto& e_pair : cell->custom_links().edges) {
            ProteinVec sigma = std::get<0>(get_polarity_proteins(e_pair));
            concentration += std::pow(arma::norm(sigma), 2);
        }

        return (  lagrange_const_proteins.at(cell) * (concentration - 1.)
                + (  0.5 * _penalty_const_proteins 
                   * std::pow((concentration - 1.), 2)));
    }

    /// Getter for energy associated with lagrange multiplier II
    /** Langrange multiplier II is the constrain of constant level of proteins
     */
    double energy_lagrange_const_concentration
    (const AgentContainer<Cell>& cells) const
    {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += energy_lagrange_const_proteins(cell);
        }
        return energy;
    }

    double total_energy (const AgentContainer<Cell>& cells, 
                         const AgentContainer<Edge>& edges) const {
        return (  energy_cell_cell_polarity(edges)
                + energy_polarity_exclusion(cells)
                + energy_lagrange_net_polarisation(cells)
                + energy_lagrange_const_concentration(cells));
    }


public:
    double energy_cell_cell_polarity () const {
        return energy_cell_cell_polarity(_cm.edges());
    }

    double energy_polarity_exclusion () const {
        return energy_polarity_exclusion(_cm.cells());
    }

    double energy_lagrange_net_polarisation () const {
        return energy_lagrange_net_polarisation(_cm.cells());
    }

    double energy_lagrange_const_concentration () const {
        return energy_lagrange_const_concentration(_cm.cells());
    }

    double total_energy () const {
        return total_energy(_cm.cells(), _cm.edges());
    }

    std::size_t proteins_per_edge () const {
        return _N;
    }

    /// The polarity of a cell
    /** Polarity proteins accumulated with a weight that is the vector
     *  between the position of the protein and the cell center
     */
    SpaceVec pcp_polarity (const std::shared_ptr<Cell>& cell) const {        
        const auto& space = _cm.get_space();

        SpaceVec center = _cm.barycenter_of(cell);
        SpaceVec polarity({0., 0.});
        for (auto [edge, flip] : cell->custom_links().edges) {
            ProteinVec sigma = std::get<0>(get_polarity_proteins(edge, flip));
            std::size_t N = sigma.n_elem;

            SpaceVec vertex = _cm.position_of(edge->custom_links().a);
            SpaceVec e_vec = _cm.displacement(edge);

            // iterate the proteins along edge
            double dN = 1. / static_cast<double>(N);
            for (std::size_t i = 0; i < N; i++) {
                double s = sigma.at(i);

                SpaceVec pos = vertex + (i + 0.5) * dN * e_vec;
                SpaceVec displ = space->displacement(center, pos);
                
                polarity += s * displ / arma::norm(displ);
            }
        }

        return polarity;
    }
};

} // namespace Utopia::Models::PlanarCellPolarity

#endif // UTOPIA_MODELS_PLANARCELLPOLARITY_HH
