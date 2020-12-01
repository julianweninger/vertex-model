#ifndef UTOPIA_MODELS_PCPVERTEX_ENTITIES_HH
#define UTOPIA_MODELS_PCPVERTEX_ENTITIES_HH

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

    /// Constructor
    VertexState ()
    :  
        f(), g(), h(),
        remove(false)
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
    /// The linetension parameter property to this edge
    double linetension;

    /// The contractility parameter
    double contractility;

    /// The time T1 transition was last attempted
    /** 0 if never attempted */
    std::size_t last_T1_attempt;

    /// Polartity protein level on side to cell a
    double sigma_a;

    /// Polartity protein level on side to cell b
    double sigma_b;

    /// Steepest descent in protein level (side a)
    double d_sigma_a;

    /// Steepest descent in protein level (side b)
    double d_sigma_b;

    /// Whether this object is to be removed 
    bool remove;

    /// Constructor
    /** \param linetension   The linetension property
     *  \param contractility The contractility parameter
     */
    EdgeState (const Utopia::DataIO::Config& cfg)
    :
        linetension(get_as<double>("linetension", cfg)),
        contractility(get_as<double>("contractility", cfg)),
        last_T1_attempt(0),
        sigma_a(0.), sigma_b(0.),
        d_sigma_a(0.), d_sigma_b(0.),
        remove(false)
    { }
};



/// The Cell defined by its id, its vertices and its area
struct CellState {
    /// The type of a cell
    enum CellType {
        progenitor,
        hair,
        support,
        num_cell_types,
    } type;

    /// The preferential area of the cell
    /** In absolute coordinates
     */
    double area_preferential;

    /// The variance of the area preferential
    double area_preferential_var;

    /// The reference shape index
    double shape_index_preferential;

    /// The preferential perimeter of the cell
    /** \f$ p_0 = shape_index_preferential * \sqrt{area_preferential} \f$
     */
    double perimeter_preferential () const {
        return shape_index_preferential * sqrt(area_preferential);
    }

    /// The contractility of the cell
    double contractility;

    /// Lagrange multiplier I
    /** Constraint of zero net polarisation
     */
    double lagrange_net_polarisation;

    /// Lagrange multiplier II
    /** Constraint of constant protein level
     */
    double lagrange_const_concentration;

    /// The initial protein concentration
    double protein_concentration;

    /// Whether this object is to be removed 
    bool remove;

    /// Create a config from the cell's properties
    DataIO::Config create_cfg_from_props () const {
        DataIO::Config cfg;
        cfg["area_preferential"] = area_preferential;
        cfg["area_preferential_var"] = area_preferential_var;
        cfg["shape_index_preferential"] = shape_index_preferential;
        cfg["contractility"] = contractility;
        cfg["protein_concentration"] = protein_concentration;
        
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
     *      - `protein_concentration` (double, unused): Currently unused
     *      -  `cell_type` (str): The type of cell. See CellState::setup_type.
     * 
     *  \param rng      A random number generator
     */
    template<class RNGType>
    CellState (const Utopia::DataIO::Config& cfg,
               const std::shared_ptr<RNGType>& rng)
    :
        type(setup_type(cfg)),
        area_preferential(get_as<double>("area_preferential", cfg)),
        area_preferential_var(get_as<double>("area_preferential_var", cfg, 0.)),
        shape_index_preferential(get_as<double>("shape_index_preferential",
                                 cfg)),
        contractility(get_as<double>("contractility", cfg)),
        lagrange_net_polarisation(0.),
        lagrange_const_concentration(0.),
        protein_concentration(get_as<double>("protein_concentration", cfg, 0.)),
        remove(false)
    {
        if (area_preferential_var > 0) {
            std::normal_distribution<> dist{area_preferential,
                                            area_preferential_var};
            area_preferential = dist(*rng);

            // give it another try
            if (area_preferential <= 0.) {
                area_preferential = dist(*rng);
            }
        }

        // give it another try, but negative area_preferential wont work ..
        if (area_preferential <= 0.)
        {
            throw std::invalid_argument(fmt::format("Cannot construct a cell "
                "with negative or 0 preferential area! Received preferential "
                "area: {}", area_preferential));
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