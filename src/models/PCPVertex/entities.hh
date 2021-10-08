#ifndef UTOPIA_MODELS_PCPVERTEX_ENTITIES_HH
#define UTOPIA_MODELS_PCPVERTEX_ENTITIES_HH

#include "utils.hh"

namespace Utopia::Models::PCPVertex {

/// The Vertex
struct VertexState {
    using SpaceVec = Utopia::SpaceVecType<2>;

    /// The steepest gradient slope
    /** \f$ f = - \nabla V (r_0) \f$ with \f$r_0\f$ the position of this vertex 
     */
    SpaceVec f;

    /// Conjugate gradient
    SpaceVec g;

    /// Direction of update (conjugate gradient)
    SpaceVec h;

    /// The position of the vertex when moved by beta along direction of update.
    /** The first value (double) is beta, the latter (SpaceVec) the virtual
     *  position.
     */
    std::pair<double, SpaceVec> virtual_pos;

    /// Whether to remove
    bool remove;

    /// Whether to fix the position of this vertex, i.e. boundary condition
    bool fix_in_space;

    /// Constructor
    VertexState ()
    :  
        f(), g(), h(),
        remove(false),
        fix_in_space(false)
    {
        f.zeros();
        g.zeros();
        h.zeros();
    }
};


/// The Edge object
/** It is defined as the shortest connection between two vertices a and b
 * 
 *  Every edge has two adjoint cells, one to either side, except boundary edges
 *  have only one.
 */
struct EdgeState {
    /// The linetension parameter
    double _linetension;

    /// Fluctuations to the linetension parameter (Ornstein-Uhlenbeck process)
    double _linetension_fluctuation;

    /// The linetension property to this edge
    double linetension() const { 
        return _linetension + _linetension_fluctuation;
    };

    /// The contractility parameter
    double _contractility;

    double contractility() const {
        if (not contractility_on) {
            return 0.;
        }
        else {
            return _contractility;
        }
    }

    bool contractility_on;
    /// The time T1 transition was last attempted
    /** 0 if never attempted */
    std::size_t last_T1_attempt;

    /// Whether this object is to be removed 
    bool remove;


    /// Constructor
    /** \param linetension   The linetension property
     *  \param contractility The contractility parameter
     */
    EdgeState (const Utopia::DataIO::Config& cfg)
    :
        _linetension(get_as<double>("linetension", cfg)),
        _linetension_fluctuation(0.),
        _contractility(get_as<double>("contractility", cfg)),
        contractility_on(true),
        last_T1_attempt(0),
        remove(false)
    { }

    Utopia::DataIO::Config get_cfg() const {
        Utopia::DataIO::Config cfg;

        cfg["linetension"] = _linetension;
        cfg["contractility"] = _contractility;

        return cfg;
    }
};


struct RotationCellState {
    double torque;

    double tracked_rotation;

    double tracking_persistence;

    RotationCellState(const Utopia::DataIO::Config& cfg)
    :
        torque(get_as<double>("torque", cfg)),
        tracked_rotation(0.),
        tracking_persistence(1. / get_as<double>("tracking", cfg, 
                                        std::numeric_limits<double>::max()))
    { }

    template <typename Cell, typename AgentManager>
    double angular_velocity(const Cell& cell, const AgentManager& am) const
    {
        using SpaceVec = typename AgentManager::SpaceVec;

        const SpaceVec cell_center = am.barycenter_of(cell);

        double angular_velocity = 0.;
        for (const auto& vertex : cell->custom_links().vertices) {
            SpaceVec pos = am.position_of(vertex);
            SpaceVec displ = am.get_space()->displacement(cell_center, pos);
            double length = arma::norm(displ);

            SpaceVec vel = vertex->state.f;

            double cross_prod = (displ[0] * vel[1] - displ[1] * vel[0]);
            angular_velocity += cross_prod / std::pow(length, 2.);
        }

        return angular_velocity / cell->custom_links().vertices.size();
    }
};


/// The Cell defined by its id, its vertices and its area
struct CellState {
    using SpaceVec = Utopia::SpaceVecType<2>;

    /// The type of a cell
    enum CellType {
        progenitor,
        hair,
        support,
        num_cell_types,
    } type;

    /// The preferential area of the cell
    double _area_preferential;

    /// Fluctuation of the A0 parameter (Ornstein-Uhlenbeck process)
    double _area_preferential_fluctuations;

    double area_preferential () const {
        return _area_preferential + _area_preferential_fluctuations;
    }

    /// The reference shape index
    double shape_index_preferential;

    /// The contractility of the cell
    double contractility;

    /// An intrinsic polarity angle wrt x axis
    double _polarity;

    double _polarity_fluctuations;

    double polarity(double rotate = 0.) const {
        return std::fmod(  _polarity + _polarity_fluctuations
                         + rotate + 2 * M_PI, 2 * M_PI);
    }

    SpaceVec polarity_vec(double rotate = 0.) {
        double pol = this->polarity(rotate);
        return SpaceVec({cos(pol), sin(pol)});
    }

    double polarity_torque;

    bool fix_polarity;

    /// Whether this object is to be removed 
    bool remove;

    /// Create a config from the cell's properties
    DataIO::Config create_cfg_from_props () const {
        DataIO::Config cfg;
        cfg["area_preferential"] = _area_preferential;
        cfg["shape_index_preferential"] = shape_index_preferential;
        cfg["contractility"] = contractility;
        cfg["polarity"] = this->polarity();
        
        if (type == CellType::progenitor) {
            cfg["cell_type"] = "progenitor";
        }
        else if (type == CellType::hair) {
            cfg["cell_type"] = "hair";
        }
        else if (type == CellType::support) {
            cfg["cell_type"] = "support";
        }
        else { 
            cfg["cell_type"] = "not implemented within division";
        }

        return cfg;
    }

    /// Constructor of a cell
    /** Construct a cell from a configuration
     * 
     *  \param cfg      The configuration
     *      -  `area_preferential` (double): The preferential size of this cell
     *      -  `area_preferential_var` (double, default: 0): The variation of
     *              the preferential area in a normal distribution
     *      -  `shape_index_preferential` (double):    The preferential shape 
     *              index of this cell. 
     *              p0 = Perimeter0 / sqrt(area preferential)
     *      -  `contractility` (double): The contractility of the cell
     *              associated with contractility of the actin-myosin ring
     *      -  `cell_type` (str): The type of cell. See CellState::setup_type.
     * 
     *  \param rng      A random number generator
     */
    template<class RNGType>
    CellState (const Utopia::DataIO::Config& cfg,
               const std::shared_ptr<RNGType>& rng)
    :
        type(setup_type(cfg)),
        _area_preferential(get_as<double>("area_preferential", cfg)),
        _area_preferential_fluctuations(0.),
        shape_index_preferential(get_as<double>("shape_index_preferential",
                                 cfg)),
        contractility(get_as<double>("contractility", cfg)),
        _polarity(get_as<double>("polarity", cfg, -M_PI_2)),
        _polarity_fluctuations(0.),
        polarity_torque(0.),
        fix_polarity(false),
        remove(false)
    {
        double area_preferential_var = get_as<double>("area_preferential_var",
                                                      cfg, 0.);
        if (area_preferential_var > 1.e-12) {
            auto dist = get_lognormal_distribution(_area_preferential,
                                                   area_preferential_var);
            _area_preferential_fluctuations = dist(*rng) - _area_preferential;
        }

        // give it another try, but negative area_preferential wont work ..
        if (area_preferential() <= 0.)
        {
            throw std::invalid_argument(fmt::format("Cannot construct a cell "
                "with negative or 0 preferential area! Received preferential "
                "area: {}", area_preferential()));
        }
        if (shape_index_preferential <= 0.) {
            throw std::invalid_argument(fmt::format("Cannot construct a cell "
                "with negative or 0 preferential shape_index! Received "
                "preferential shape_index: {}", shape_index_preferential));
        }
    }

private:
    /// Setup the type of the cell from config
    /** \param cfg      The configuration containing
     *      - `cell_type` (str): Can be one of
     *             - `progenitor`
     *             - `hair`
     *             - `support`
     */
    CellType setup_type (const Utopia::DataIO::Config& cfg) {
        auto cell_type = get_as<std::string>("cell_type", cfg);

        if (cell_type == "progenitor") {
            return CellType::progenitor;
        }
        else if (cell_type == "hair") {
            return CellType::hair;
        }
        else if (cell_type == "support") {
            return CellType::support;
        }
        else {
            throw KeyError(cell_type, cfg, "Cell type can be 'progenitor', "
                           "'hair', or 'support'.");
        }
    }
};

} // namespace Utopia::Models::PCPVertex

#endif