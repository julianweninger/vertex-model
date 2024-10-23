#ifndef UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION2DPLUSTRIANGLE_HH
#define UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION2DPLUSTRIANGLE_HH


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction2Dplus {

/// @brief The base class for Triangle Approximation of 2D plus WF class
/// @tparam Model 
/** The base class defines the background tissue height, as well as height and 
 *  height derivative cell parameter names.
 * 
 *  Parameters:
 *      - `tissue_height` (double): height of cells with basal contact
 * 
 */
template <typename Model>
class TriangleApproximationBase : public WorkFunction::WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunction::WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// @brief The height of the tissue
    double _tissue_height;

    /// @brief  The name of the cell's parameter referring to cell's height
    const std::string _height;

    /// @brief  The name of the cell's parameter referring to cell's height
    /** The derivate dW/dH. 
     *  Note that the resulting force is F_z = -1/\gamma dW/dH 
     */
    const std::string _height_derivative;

public:
    /// @brief Check whether cell has a basal contact
    /// @param height The height of considered cell
    /// @return Whether cell has a basal contact
    bool has_basal_contact(const std::shared_ptr<Cell>& cell) const {
        #if UTOPIA_DEBUG
            if (not cell) {
                throw std::runtime_error(
                    "Checking for Basal contact on a boundary cell in "
                    "PCPVertex::WorkFunction::VolumeElasticityTriangle!"
                );
            }
        #endif

        return cell->state.type != 1;
    }

    /// @brief  A factor that rescales depth of triangles to match cell's area
    double triangle_factor (const std::shared_ptr<Cell>& cell) const {
        #if UTOPIA_DEBUG
            if (not cell) {
                throw std::runtime_error(
                    "Checking for triangle_factor on a boundary cell in "
                    "PCPVertex::WorkFunction::VolumeElasticityTriangle!"
                );
            }
        #endif

        double l2 = 0.;
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            auto [__c, nb] = this->_am.adjoints_of(edge);
            if (__c != cell) { std::swap(__c, nb); }

            // Exclude HCs and boundary
            if (not nb or nb->state.type == 1) {
                continue;
            }

            l2 += std::pow(this->_am.length_of(edge), 2);
        }

        if (l2 < 1.e-8) {
            return 0.;
        }

        return 2 * this->_am.area_of(cell) / l2;
    }

    /// @brief  The height of a cell
    /// @param cell 
    /// @return 
    inline double height_of (const std::shared_ptr<Cell>& cell) const {
        #if UTOPIA_DEBUG
            if (not cell) {
                throw std::runtime_error(
                    "Checking for height on a boundary cell in "
                    "PCPVertex::WorkFunction::VolumeElasticityTriangle!"
                );
            }
        #endif

        return cell->state.get_parameter(_height);
    }
    
    TriangleApproximationBase (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _tissue_height(get_as<double>("tissue_height", cfg)),
        _height("TriangleApproximation__height"),
        _height_derivative(_height + "_derivative")
    {}

    /// @brief The instruction how to calculate forces acting on vertices
    ///        from this term, i.e the derivative of the energy.
    void compute_and_set_forces () override { }

    /// @brief  The total energy for a subset of entities.
    /** @param vertices  The container of vertices
     *  @param edges     The container of edges
     *  @param cells     The container of cells
     *  
     *  @return energy (double)
    **/ 
    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const override
    {
        return 0.;
    }

    /// @brief  Update the WF's parameters
    /// @param cfg      The update configuration
    /** Parameters:
     *      - `tissue_height` (double): Updates the height of the tissue
     *      - `set_height` (double): Sets the height of cells without basal 
     *              contact
     *  
     */
    void update_parameters (const DataIO::Config& cfg) override {
        double tmp = get_as<double>("tissue_height", cfg, _tissue_height);
        if (fabs(tmp - _tissue_height) > 1.e-8) {
            _tissue_height = tmp;
            this->_am.get_logger()->debug("Setting height of cells excluding "
                "type 1 to {}", _tissue_height);
            for (const auto& cell : this->_am.cells()) {
                if (has_basal_contact(cell)) {
                    cell->state.update_parameter(_height, _tissue_height);

                }
            }
        }
        if (cfg["set_height"]) {
            auto height = get_as<double>("set_height", cfg);
            this->_am.get_logger()->debug("Setting height of cells with type 1 "
                "to {}", height);
            for (const auto& cell : this->_am.cells()) {
                if (not has_basal_contact(cell)) {
                    cell->state.update_parameter(_height, height);
                }
            }
        }
    }

};

/// The update function to the triangle Approximation
/** Performs a steepest gradient update of height variables in cells without
 *  basal contact
 * 
 *  Parameters:
 *      - `gamma` (double): Damping factor relative to time step size of 
 *              vertex-model.
 *      - `tissue_height` (double): height of cells with basal contact
 *      - `maximum_height` (double): Upper limit for cell height
 *      - `minimum_height` (double): Lower limit for cell height
 *      - `set_height` (double, default: tissue_height): Initial value
 *              of cell height for cells without basal contact
 *      - Parameters of TriangleApproximationBase
 */
template <typename Model>
class TriangleApproximationUpdate : public TriangleApproximationBase<Model>
{
public:
    using Base = TriangleApproximationBase<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// @brief  Damping factor for hight development additional to dt
    double _gamma;

    /// @brief Maximum value for cell height
    double _maximum_height;

    /// @brief Minimum value for cell height
    double _minimum_height;

public:
    TriangleApproximationUpdate (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _gamma(get_as<double>("gamma", cfg)),
        _maximum_height(get_as<double>("maximum_height", cfg)),
        _minimum_height(get_as<double>("minimum_height", cfg))
    {
        auto height = get_as<double>("set_height", cfg, this->_tissue_height);
        for (const auto& cell : this->_am.cells()) {
            if (cell->state.type != 1) {
                cell->state.register_parameter(
                        this->_height,this->_tissue_height);
            }
            else {
                cell->state.register_parameter(this->_height, height);
            }
            cell->state.register_parameter(this->_height_derivative, 0.);
            cell->state.register_parameter(
                    this->_height_derivative + "_monitor", 0.);
        }
    }

    ~TriangleApproximationUpdate () {
        for (const auto& cell : this->_am.cells()) {
            cell->state.unregister_parameter(this->_height);
            cell->state.unregister_parameter(this->_height_derivative);
            cell->state.unregister_parameter(
                    this->_height_derivative + "_monitor");
        }
    }

    /// @brief The instruction how to calculate forces acting on vertices
    ///        from this term, i.e the derivative of the energy.
    void compute_and_set_forces () override { }

    /// @brief  The total energy for a subset of entities.
    /** @param vertices  The container of vertices
     *  @param edges     The container of edges
     *  @param cells     The container of cells
     *  
     *  @return energy (double)
    **/ 
    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const override
    {
        return 0.;
    }
    
    /// Performs the steepest gradient update of \f$ H_\alpha \f$ variables
    void update (double dt) override {
        for (const auto& cell : this->_am.cells()) {
            if (this->has_basal_contact(cell)) {
                continue;
            }

            double H = cell->state.get_parameter(this->_height);
            double dW_dH = cell->state.get_parameter(this->_height_derivative);
            H -= dt * _gamma * dW_dH;
            H = std::max( std::min(H, _maximum_height), _minimum_height );
            cell->state.update_parameter(this->_height, H);
            cell->state.update_parameter(
                this->_height_derivative + "_monitor", dW_dH
            );
            cell->state.update_parameter(this->_height_derivative, 0.);
        }
    }

    /// Updates parameter values
    void update_parameters (const DataIO::Config& cfg) override {
        Base::update_parameters(cfg);
        _gamma = get_as<double>("gamma", cfg, _gamma);
        _maximum_height = get_as<double>("maximum_height",cfg,_maximum_height);
        _minimum_height = get_as<double>("minimum_height",cfg,_minimum_height);
    }

    std::vector<std::string> write_task_cell_properties_names () const override
    {
        return std::vector<std::string>({"height"});
    }

    std::vector<std::vector<double>> write_cell_properties () const override {
        std::vector<double> heights({});
        for (const auto& cell : this->_am.cells()) {
            heights.push_back(cell->state.get_parameter(this->_height));
        }
        return std::vector<std::vector<double>>({ heights });
    }
};

/// The update function to the triangle Approximation using a graded damping factor
/** Performs a steepest gradient update of height variables in cells without
 *  basal contact. However, depending on x position of cells height variable
 *  develops faster or slower.
 * 
 *  Parameters:
 *      - `gamma` (double): Damping factor relative to time step size of 
 *              vertex-model.
 *      - `gamma_fold_decrease` (double, > 0): Damping factor fold decrease at right-
 *          most cells. A fold decrease of 2 indicates, that gamma is halved.
 *          A fold decrease of 1 indicates a constant gamma. A fold decrease of
 *          0.5 indicates, that gamma doubles. Note that characteristic time
 *          scales with 1/gamma.
 *      - `tissue_height` (double): height of cells with basal contact
 *      - `maximum_height` (double): Upper limit for cell height
 *      - `minimum_height` (double): Lower limit for cell height
 *      - `set_height` (double, default: tissue_height): Initial value
 *              of cell height for cells without basal contact
 */
template <typename Model>
class TriangleApproximationUpdateGraded : public TriangleApproximationUpdate<Model>
{
public:
    using Base = TriangleApproximationUpdate<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// @brief  Fold change in damping factor gamma between left and right
    /** A fold change  */
    double _gamma_fold_change;

    /// The origin and length (horizontal) of the domain
    std::pair<SpaceVec, SpaceVec> _domain_info;

    /// Caluclate the origin and length of the domain
    /** Based on the vertex positions and periodicity of space.
     *  It deals with open, semi-periodic and periodic BC.
     */
    std::pair<SpaceVec, SpaceVec> get_domain_info () const {
        // Establish domain size
        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::min();
        double y_min = std::numeric_limits<double>::max();
        double y_max = std::numeric_limits<double>::min();

        const auto& space = this->_am.get_space();
        SpaceVec domain = space->get_domain_size();
        SpaceVec ref = this->_am.position_of(this->_am.vertices()[0]);
        for (const auto& v : this->_am.vertices()) {
            SpaceVec pos = space->displacement(ref, this->_am.position_of(v));
            x_min = std::min(x_min, pos[0]);
            x_max = std::max(x_max, pos[0]);
            y_min = std::min(y_min, pos[1]);
            y_max = std::max(y_max, pos[1]);
        }
        double Lx = (x_max - x_min);
        double Ly = (y_max - y_min);

        // is truly periodic space
        if (space->periodic and Lx > 0.9 * domain[0]) {
            return std::make_pair(domain / 2., domain);
        }

        return std::make_pair(
            // origin at center
            SpaceVec({x_min + ref[0] + Lx/2., y_min + ref[1] + Ly/2.}),
            SpaceVec({Lx, Ly})  // extent
        );
    }

    /// @brief  Get the damping factor at cell's position
    /// @param cell 
    double get_gamma (const std::shared_ptr<Cell>& cell) const { 
        SpaceVec origin = std::get<0>(_domain_info);
        SpaceVec extent = std::get<1>(_domain_info);

        SpaceVec pos = this->_am.barycenter_of(cell);
        pos = this->_am.get_space()->displacement(origin, pos);
        SpaceVec rpos = pos / extent;

        // NOTE rpos in [-0.5, 0.5] for x and y
        return this->_gamma / ((rpos[0] + 0.5)*_gamma_fold_change + 1);
    }

public:
    TriangleApproximationUpdateGraded (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _gamma_fold_change(get_as<double>("gamma_fold_decrease", cfg)-1),
        _domain_info(get_domain_info())
    {
        if (_gamma_fold_change + 1 < 1.e-8) {
            throw std::runtime_error(fmt::format(
                "TriangleApproximationUpdateGraded: "
                "Fold decrease must be larger than 0., but was {}!",
                _gamma_fold_change + 1.
            ));
        }
    }
    
    /// Performs the update of \f$ H_\alpha \f$, where the damping factor 
    /// depends on the cell's position
    void update (double dt) override {
        _domain_info = get_domain_info();

        for (const auto& cell : this->_am.cells()) {
            if (this->has_basal_contact(cell)) {
                continue;
            }

            double gamma = get_gamma(cell);

            double H = cell->state.get_parameter(this->_height);
            double dW_dH = cell->state.get_parameter(this->_height_derivative);
            H -= dt * gamma * dW_dH;
            H = std::max(std::min(H, this->_maximum_height),
                         this->_minimum_height);
            cell->state.update_parameter(this->_height, H);
            cell->state.update_parameter(
                this->_height_derivative + "_monitor", dW_dH
            );
            cell->state.update_parameter(this->_height_derivative, 0.);
        }
    }
};

/// Setter for height variables in Triangle Approximation
/** Parameters:
 *      - `target_volume` (double): Reference volume from which set height
 *              is calculated as volume/area.
 *      - `area` (double): The mean set area. Set height is calculated as 
 *              target_volume / area.
 *      - `area_gradient` (double): The difference between x=-0.5 and x=0.5,
 *              mean area (`area`) is set for x=0.
 *      - `tissue_height` (double): height of cells with basal contact
 *      - `maximum_height` (double): Upper limit for cell height
 *      - `minimum_height` (double): Lower limit for cell height
 */
template <typename Model>
class TriangleApproximationSetter : public TriangleApproximationBase<Model>
{
public:
    using Base = TriangleApproximationBase<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// @brief Cell target volume
    double _V0;

    /// @brief Mean value of cell area
    double _area;

    /// @brief Gradient of cell area, where left most is value - grad / 2.
    double _area_gradient;

    /// @brief Maximum value for cell height
    double _maximum_height;

    /// @brief Minimum value for cell height
    double _minimum_height;

    /// The origin and length (horizontal) of the domain
    std::pair<SpaceVec, SpaceVec> _domain_info;

    /// Caluclate the origin and length of the domain
    /** Based on the vertex positions and periodicity of space.
     *  It deals with open, semi-periodic and periodic BC.
     */
    std::pair<SpaceVec, SpaceVec> get_domain_info () const {
        // Establish domain size
        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::min();
        double y_min = std::numeric_limits<double>::max();
        double y_max = std::numeric_limits<double>::min();

        const auto& space = this->_am.get_space();
        SpaceVec domain = space->get_domain_size();
        SpaceVec ref = this->_am.position_of(this->_am.vertices()[0]);
        for (const auto& v : this->_am.vertices()) {
            SpaceVec pos = space->displacement(ref, this->_am.position_of(v));
            x_min = std::min(x_min, pos[0]);
            x_max = std::max(x_max, pos[0]);
            y_min = std::min(y_min, pos[1]);
            y_max = std::max(y_max, pos[1]);
        }
        double Lx = (x_max - x_min);
        double Ly = (y_max - y_min);

        // is truly periodic space
        if (space->periodic and Lx > 0.9 * domain[0]) {
            return std::make_pair(domain / 2., domain);
        }

        return std::make_pair(
            // origin at center
            SpaceVec({x_min + ref[0] + Lx/2., y_min + ref[1] + Ly/2.}),
            SpaceVec({Lx, Ly})  // extent
        );
    }

    /// @brief  Update the height parameter of a cell
    /// @param cell 
    void update_height (const std::shared_ptr<Cell>& cell) const { 
        SpaceVec origin = std::get<0>(_domain_info);
        SpaceVec extent = std::get<1>(_domain_info);

        SpaceVec pos = this->_am.barycenter_of(cell);
        pos = this->_am.get_space()->displacement(origin, pos);
        SpaceVec rpos = pos / extent;

        // NOTE rpos in [-0.5, 0.5] for x and y
        double area = _area + rpos[0] * _area_gradient;
        double height = _V0 / area;        
        height = std::min(std::max(height, _minimum_height), _maximum_height);

        cell->state.update_parameter(this->_height, height);
    }

public:
    TriangleApproximationSetter (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _V0(get_as<double>("target_volume", cfg)),
        _area(get_as<double>("area", cfg)),
        _area_gradient(get_as<double>("area_gradient", cfg)),
        _maximum_height(get_as<double>("maximum_height", cfg)),
        _minimum_height(get_as<double>("minimum_height", cfg)),
        _domain_info(get_domain_info())
    {
        for (const auto& cell : this->_am.cells()) {
            if (cell->state.type != 1) {
                cell->state.register_parameter(
                        this->_height, this->_tissue_height);
            }
            else {
                cell->state.register_parameter(this->_height, 0.);
                update_height(cell);
            }
            cell->state.register_parameter(this->_height_derivative, 0.);
            cell->state.register_parameter(
                    this->_height_derivative + "_monitor", 0.);
        }
    }

    ~TriangleApproximationSetter () {
        for (const auto& cell : this->_am.cells()) {
            cell->state.unregister_parameter(this->_height);
            cell->state.unregister_parameter(this->_height_derivative);
            cell->state.unregister_parameter(
                    this->_height_derivative + "_monitor");
        }
    }

    /// @brief The instruction how to calculate forces acting on vertices
    ///        from this term, i.e the derivative of the energy.
    void compute_and_set_forces () override { }

    /// @brief  The total energy for a subset of entities.
    /** @param vertices  The container of vertices
     *  @param edges     The container of edges
     *  @param cells     The container of cells
     *  
     *  @return energy (double)
    **/ 
    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const override
    {
        return 0.;
    }
    
    /// Performs the update of \f$ H_\alpha \f$
    void update (double dt) override {
        _domain_info = get_domain_info();

        for (const auto& cell : this->_am.cells()) {
            if (this->has_basal_contact(cell)) {
                continue;
            }
            
            update_height(cell);
            cell->state.update_parameter(
                this->_height_derivative + "_monitor",
                cell->state.get_parameter(this->_height_derivative)
            );
            cell->state.update_parameter(this->_height_derivative, 0.);
        }
    }

    /// Updates parameters
    /** Parameters:
     *      - `target_volume` (double, optional):
     *              Reference volume from which set height is calculated as 
     *              volume/area.
     *      - `area` (double, optional):
     *              The mean set area. Set height is calculated as 
     *              target_volume / area.
     *      - `increment_area` (double, optional):
     *              The incremental to mean set area. 
     *              Set height is calculated as target_volume / area.
     *      - `area_gradient` (double, optional):
     *              The difference between x=-0.5 and x=0.5, mean area (`area`) 
     *              is set for x=0.
     *      - `increment_area_gradient` (double): The incremental to gradient 
     *              area. 
     *      - `tissue_height` (double, optional):
     *              height of cells with basal contact
     *      - `maximum_height` (double, optional):
     *              Upper limit for cell height
     *      - `minimum_height` (double, optional):
     *              Lower limit for cell height
     */
    void update_parameters (const DataIO::Config& cfg) override {
        Base::update_parameters(cfg);
        if (cfg["increment_area"]) {
            _area += get_as<double>("increment_area", cfg);
        }
        else {
            _area = get_as<double>("area", cfg, _area);
        }
        if (cfg["increment_area_gradient"]) {
            _area_gradient += get_as<double>("increment_area_gradient",cfg);
        }
        else {
            _area_gradient = get_as<double>(
                    "area_gradient", cfg, _area_gradient);
        }

        _V0 = get_as<double>("target_volume", cfg, _V0);
        _maximum_height = get_as<double>("maximum_height",cfg,_maximum_height);
        _minimum_height = get_as<double>("minimum_height",cfg,_minimum_height);
    }

    std::vector<std::string> write_task_cell_properties_names () const override
    {
        return std::vector<std::string>({"height", "height_derivative"});
    }

    std::vector<std::vector<double>> write_cell_properties () const override {
        std::vector<double> heights({});
        std::vector<double> dh({});
        for (const auto& cell : this->_am.cells()) {
            heights.push_back(cell->state.get_parameter(this->_height));
            dh.push_back(
                cell->state.get_parameter(this->_height_derivative+"_monitor"));
        }
        return std::vector<std::vector<double>>({ heights, dh });
    }
};




/// @brief The volume elasticity of a cell base class
/** \f$ E_\alpha = k/2 (V_\alpha / V^{(0)} - 1)^2\f$, an elastic penalty on 
 *  cell volume  \f$ V_\alpha \f$.
 *
 *  Requires an implementation of preferential_volume(cell).
 *  
 *  Cell volume is calculated from a columnar model, using cell apical area and 
 *  a cell height variable H. In this columnar model, cells of type 1 can detach
 *  basally, with \f$ H < H_{tissue} \f$. This basal space is filled up by 
 *  neighboring cells of other types, having a columnar core and basal 
 *  extensions.
 *  
 *  The variables \f$ H_\alpha \f$ develop according to the derivative 
 *  \f$ \frac{dH_\alpha}{dt} = - \gamma \frac{dE}{dH_\alpha} \f$.
 * 
 *  Parameters:
 *      - `elastic_modulus`: (double) the elastic modulus \f$ k \f$
 *      - `preferential_volume`: (double) the preferential volume 
 *              \f$ V^{(0)} \f$
 *      - `tissue_height`: (double) The height of cells of type other than 1, 
 *              also used for initialisation for \f$ H_\alpha \f$.
 *      - `minimum_height`: (double) The minimum value for \f$ H_\alpha \f$. 
 *      - `height_parameter_name`: (string) Name used for the height variable in 
 *              properties of the cell. Can be accessed by other terms. 
 *              The derivative \f$ H_\alpha^\prime \f$ is saved as 
 *              <height_parameter_name>_derivative. Additive contributions from
 *              terms are considered during update.
 *      - `gamma`: (double) The damping factor for H temporal development, 
 *              additional to numeric step size
 */
template <typename Model>
class VolumeElasticityTriangleBase : public TriangleApproximationBase<Model>
{
public:
    using Base = TriangleApproximationBase<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// @brief  Elastic constant associated with volume elasticity
    double _elastic_modulus;

    virtual double preferential_volume(
            const std::shared_ptr<Cell>& cell) const=0;

    /// @brief The volume of a cell
    /** If cell does not have a basal contact, it is considered columnar.
     *  Else, it is columnar and occupies an equal share of those neighbours
     *  that do not have a basal contact.
     *  
     *  @param cell Pointer to cell
     *  @return Volume of cell
     */
    double volume_of (const std::shared_ptr<Cell>& cell) const {
        #if UTOPIA_DEBUG
            if (not cell) {
                throw std::runtime_error(
                    "Checking for Volume on a boundary cell in "
                    "PCPVertex::WorkFunction::VolumeElasticity!"
                );
            }
        #endif

        const double height = this->height_of(cell);
        const double area = this->_am.area_of(cell);

        double volume = area * height;

        // It is columnar itself, no neighbour contributions
        if (not this->has_basal_contact(cell)) {
            return volume;
        }

        // check neighbours for basal detachment
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            auto [c, n] = this->_am.adjoints_of(edge);
            if (c != cell) { std::swap(c, n); }

            // if neighbour has basal contact, does not contribute
            if (n and not this->has_basal_contact(n)) { 
                volume +=   std::pow(this->_am.length_of(edge), 2)
                          * (this->_tissue_height - this->height_of(n))
                          * this->triangle_factor(n) / 2;
            }
        }

        return volume;
    }

public:
    VolumeElasticityTriangleBase (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _elastic_modulus(get_as<double>("elastic_modulus", cfg))
    { }

    void compute_and_set_forces () override {
        // calculate some quantities
        std::unordered_map<std::shared_ptr<Cell>, double> volumes;
        std::unordered_map<std::shared_ptr<Cell>, double> pressures;
        for (const auto& cell : this->_am.cells()) {
            double V  = volume_of(cell);
            double V0 = preferential_volume(cell);

            volumes  [cell] = V;
            pressures[cell] = _elastic_modulus * (V - V0) / std::pow(V0, 2);
        }
        std::unordered_map<std::shared_ptr<Edge>, double> lengths;
        for (const auto& edge : this->_am.edges()) {
            lengths[edge] = this->_am.length_of(edge);
        }

        // iterate all cells
        for (const auto& cell : this->_am.cells()) {
            double dW_dH = 0.;     // the derivative to height
            double dW_dV = pressures[cell];

            double H = cell->state.get_parameter(this->_height);
            double volume = volumes[cell];
            double f = 0.;      // the triangle factor
            double area = 0.;

            // derivative w.r.t. area
            double dWc_dA = dW_dV * H;
            double dWn_dA = 0.;     // accumulate from neighbors
            double L2 = 0.;
            double dW_df = 0.;
            if (not this->has_basal_contact(cell)) {
                f = this->triangle_factor(cell);
                area += this->_am.area_of(cell);
                dW_dH += dW_dV * area;

                // and neighbour derivatives
                for (const auto& [edge, flip] : cell->custom_links().edges) {
                    auto [cell__, n] = this->_am.adjoints_of(edge);
                    if (cell__ != cell) { std::swap(cell__, n); }

                    // account for A in nbs volume (within factor)
                    if (n and this->has_basal_contact(n)) {
                        double p = pressures[n];
                        double l2 = std::pow(lengths[edge], 2);
                        dWn_dA += p * l2;
                        L2 += l2;
                        dW_df += 0.5 * p * l2 * (this->_tissue_height - H);
                    }
                }
                if (L2 > 1.e-8) {
                    dWn_dA = (this->_tissue_height - H) * dWn_dA / L2;
                }
                // A fictional neighbour than wants to shrink to 0 volume
                else {
                    dW_dH -= _elastic_modulus * std::pow(area, 2)
                             * (this->_tissue_height - H)
                             / std::pow(0.01 * preferential_volume(cell), 2);
                }
            }

            // apply forces via edges
            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto a = edge->custom_links().a;
                auto b = edge->custom_links().b;
                if (flip) { std::swap(a, b); }

                SpaceVec displ = this->_am.displacement(a, b);
                SpaceVec normal({displ[1], -displ[0]});
                double l = lengths[edge];

                auto [c, n] = this->_am.adjoints_of(edge);
                if (c != cell) { std::swap(c, n); }


                // consider neighbor contributions
                if (    not this->has_basal_contact(cell) 
                    and n 
                    and this->has_basal_contact(n))
                {
                    double dW_dVn = pressures[n];
                    // The own and neighbor contributions to tension
                    double dW_dl = (
                          dW_dVn * l * f * (this->_tissue_height-H)
                        - dW_df * 4 * area / std::pow(L2, 2) * l
                    );

                    a->state.add_force(+ dW_dl * displ / l);
                    b->state.add_force(- dW_dl * displ / l);

                    dW_dH -= dW_dVn * std::pow(l, 2) * f / 2;
                }


                // apply force from pressure
                SpaceVec force = -(dWc_dA + dWn_dA) * normal / 2.;
                a->state.add_force(force);
                b->state.add_force(force);
            }

            cell->state.update_parameter(
                this->_height_derivative, 
                cell->state.get_parameter(this->_height_derivative) + dW_dH
            );
        }
    }

    /// @brief The derivative of W to area of this cell
    /// @param cell 
    /// @return apical pressure
    double compute_pressure(const std::shared_ptr<Cell>& cell) const override {
        double V0 = preferential_volume(cell);
        return _elastic_modulus * (volume_of(cell) - V0) / std::pow(V0, 2);
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const override {
        double V = volume_of(cell);
        double V0 = preferential_volume(cell);
        return 0.5 * _elastic_modulus * std::pow(V / V0 - 1., 2);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        const AgentContainer<Cell>& cells
    ) const override
    {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += compute_energy(cell);
        }
        return energy;
    }

    void update_parameters (const DataIO::Config& cfg) override {
        Base::update_parameters(cfg);

        _elastic_modulus = get_as<double>("elastic_modulus", cfg, 
                                          _elastic_modulus);
    }

    bool test_constraints (const std::shared_ptr<spdlog::logger>& logger) 
    const override 
    {
        bool PASS_TEST = Base::test_constraints(logger);
        logger->debug("Testing constraints of WF term {} ...", this->_name);


        logger->debug("   Testing height constraint ...");
        bool test_height = true;
        for (const auto& cell : this->_am.cells()) {
            double height = cell->state.get_parameter(this->_height);
            if (height > this->_tissue_height + 1.e-8) {
                test_height = false;
                logger->error("Cell {} (height: {}) exceeds the tissue height "
                              "({})!",
                              cell->id(), height, this->_tissue_height);
            }
            if (height < -1.e-8) {
                test_height = false;
                logger->error("Cell {} has negative hight (height: {} < 0)!",
                              cell->id(), height);
            }
        }
        if (test_height) {
            logger->debug("   Height constraint is fulfilled.");
        }
        else {
            PASS_TEST = false;
            logger->error("  Test of height constraint failed! ");
        }


        logger->debug("   Testing volume constraint ...");
        double tot_volume = 0.;
        for (const auto& cell : this->_am.cells()) {
            tot_volume += volume_of(cell);
        }
        SpaceVec domain = this->_am.get_space()->get_domain_size();
        double area = domain[0] * domain[1];
        bool test_volume = fabs(tot_volume - area * this->_tissue_height) 
                           / tot_volume < 1.e-6;
        if (test_volume) {
            logger->debug("   Volume constraint is fulfilled.");
        }
        else {
            PASS_TEST = false;
            logger->error("  Test of volume constraint failed! "
                          "Total volume was {}, expected volume was {}. "
                          "Relative error {} exceeds {}.",
                          tot_volume, area * this->_tissue_height, 
                          fabs(tot_volume-area*this->_tissue_height)/tot_volume,
                          1.e-6);
        }


        if (PASS_TEST) {
            logger->debug("All tests of constraints of WF term {} passed.",
                          this->_name);
        }
        else {
            logger->error("Test of constraints of WF term {} FAILED!",
                          this->_name);
        }
        return PASS_TEST;
    }

    std::vector<std::string> write_task_cell_properties_names () const override {
        return std::vector<std::string>({"volume"});
    }

    std::vector<std::vector<double>> write_cell_properties () const override {
        std::vector<double> volumes({});

        for (const auto& cell : this->_am.cells()) {
            volumes.push_back(volume_of(cell));
        }
        return std::vector<std::vector<double>>({ volumes });
    }
    

    std::vector<std::string> write_task_cell_energies_names () const override {
        return std::vector<std::string>({
            "energy",
            "pressure",
        });
    }

    std::vector<std::vector<double>> write_cell_energies () const override {
        std::vector<double> energies({});
        std::vector<double> pressures({});
        
        for (const auto& cell : this->_am.cells()) {
            energies.push_back(compute_energy(cell));
            pressures.push_back(compute_pressure(cell));
        }

        return std::vector<std::vector<double>>({
            energies,
            pressures,
        });
    }
};


/// @brief The volume elasticity of a cell with heterotypic target volume
/** Specialisation of the VolumeElasticityTriangleBase class.
 *  Parameters:
 *      - `elastic_modulus`: (double) the elastic modulus \f$ k \f$
 *      - `preferential_volume`: (list(double))the preferential volume 
 *              \f$ V^{(0)} \f$ per cell type.
 *      - `tissue_height`: (double) The height of cells of type other than 1, 
 *              also used for initialisation for \f$ H_\alpha \f$.
 *      - `critical_height`: (double) The crossover height at which HCs loose 
 *              basal contact.
 *      - `minimum_height`: (double) The minimum value for \f$ H_\alpha \f$. 
 *      - `height_parameter_name`: (string) Name used for the height variable in 
 *              properties of the cell. Can be accessed by other terms. 
 *              The derivative \f$ H_\alpha^\prime \f$ is saved as 
 *              <height_parameter_name>_derivative. Additive contributions from
 *              terms are considered during update.
 *      - `gamma`: (double) The damping factor for H temporal development, 
 *              additional to numeric step size
 */
template <typename Model>
class VolumeElasticityTriangle : public VolumeElasticityTriangleBase<Model>
{
public:
    using Base = VolumeElasticityTriangleBase<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// @brief The target volume for a cell at which there is no pressure
    double _preferential_volume;

    inline double preferential_volume
            ([[maybe_unused]] const std::shared_ptr<Cell>& cell) const override
    {
        return _preferential_volume;
    }

public:
    VolumeElasticityTriangle (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _preferential_volume(get_as<double>("preferential_volume", cfg))
    { }

    void update_parameters (const DataIO::Config& cfg) override {
        _preferential_volume = get_as<double>("preferential_volume", cfg, 
                                              _preferential_volume);
        
        Base::update_parameters(cfg);
    }
};

/// @brief The volume elasticity of cells, dependent on cell type
/** \f$ E_\alpha = k/2 (V_\alpha / V^{(0)} - 1)\f$, an elastic penalty on 
 *  cell volume  \f$ V_\alpha \f$.
 *  
 *  Cell volume is calculated from a columnar model, using cell apical area and 
 *  a cell height variable H. In this columnar model, cells of type 1 can detach
 *  basally, with \f$ H < H_{tissue} \f$. This basal space is filled up by 
 *  neighboring cells of other types, having a columnar core and basal 
 *  extensions.
 *  
 *  The variables \f$ H_\alpha \f$ develop according to the derivative 
 *  \f$ \frac{dH_\alpha}{dt} = - \gamma \frac{dE}{dH_\alpha} \f$.
 * 
 *  Parameters:
 *      - `elastic_modulus`: the elastic modulus \f$ k \f$
 *      - `preferential_volume`: the preferential volume \f$ V^{(0)} \f$ 
 *              provided as a vector of doubles, with incrementing cell type.
 *      - `tissue_height`: The height of cells of type other than 1, also used
 *              for initialisation for \f$ H_\alpha \f$.
 *      - `minimum_height`: The minimum value for \f$ H_\alpha \f$. 
 *      - `height_parameter_name`: Name used for the height variable in 
 *              properties of the cell. Can be accessed by other terms. 
 *              The derivative \f$ H_\alpha^\prime \f$ is saved as 
 *              <height_parameter_name>_derivative. Additive contributions from
 *              terms are considered during update.
 *      - `gamma`: The damping factor for H temporal development, additional 
 *              to numeric step size
 */
template <typename Model>
class VolumeElasticityTriangleHeterotypic 
: 
    public VolumeElasticityTriangleBase<Model>
{
public:
    using Base = VolumeElasticityTriangleBase<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// @brief The target volume for cells dependent on their cell-type
    std::vector<double> _preferential_volume;

    double preferential_volume (const std::shared_ptr<Cell>& cell)
    const override
    {
        if (cell->state.type > _preferential_volume.size())
        {
            throw std::runtime_error(fmt::format(
                "In WF-term VolumeElasticityTriangleHeterotypic, no "
                "preferential volume registered for cells of type {}. "
                "Parameters only for {} types registered.",
                cell->state.type,
                _preferential_volume.size()
            ));
        }
        return _preferential_volume[cell->state.type];
    }

public:
    VolumeElasticityTriangleHeterotypic (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _preferential_volume(
            get_as<std::vector<double>>("preferential_volume", cfg))
    { }

    void update_parameters (const DataIO::Config& cfg) override {
        _preferential_volume = get_as<std::vector<double>>(
            "preferential_volume", 
            cfg,
            _preferential_volume
        );
        
        Base::update_parameters(cfg);
    }
};


/// @brief The surface elasticity of a cell
/** \f$ E_\alpha = k/2 (S_\alpha / S^{(0)} - 1)^2\f$, an elastic penalty on 
 *  cell volume  \f$ V_\alpha \f$.
 *
 *  Requires an implementation of preferential_volume(cell).
 *  
 *  Cell volume is calculated from a columnar model, using cell apical area and 
 *  a cell height variable H. In this columnar model, cells of type 1 can detach
 *  basally, with \f$ H < H_{tissue} \f$. This basal space is filled up by 
 *  neighboring cells of other types, having a columnar core and basal 
 *  extensions.
 *  
 *  The variables \f$ H_\alpha \f$ develop according to the derivative 
 *  \f$ \frac{dH_\alpha}{dt} = - \gamma \frac{dE}{dH_\alpha} \f$.
 * 
 *  Parameters:
 *      - `elastic_modulus`: (double) the elastic modulus \f$ k \f$
 *      - `preferential_shape`: (double) the target shape
 *              \f$ S^{(0)} \f$
 *      - `tissue_height`: (double) The height of cells of type other than 1, 
 *              also used for initialisation for \f$ H_\alpha \f$.
 *      - `critical_height`: (double) The crossover height at which HCs loose 
 *              basal contact.
 *      - `minimum_height`: (double) The minimum value for \f$ H_\alpha \f$. 
 *      - `VolumeElasticity_term`: (string) Name of the workfunction registering
 *              the height parameter.
 *      - `height_parameter_name`: (string) Name used for the height variable in 
 *              properties of the cell. Needs to be same as registered in 
 *              `VolumeElasticityTerm`. `VolumeElasticityTerm` needs to update
 *              height parameter from gradient.
 */
template <typename Model>
class SurfaceElasticity : public TriangleApproximationBase<Model>
{
public:
    using Base = TriangleApproximationBase<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// @brief  Elastic constant associated with volume elasticity
    double _elastic_modulus;

    /// @brief The target volume for a cell at which there is no pressure
    double _preferential_surface;

    /// @brief The volume of a cell
    /** If cell does not have a basal contact, it is considered columnar.
     *  Else, it is columnar and occupies an equal share of those neighbours
     *  that do not have a basal contact.
     *  
     *  @param cell Pointer to cell
     *  @return Volume of cell
     */
    double surface_of (const std::shared_ptr<Cell>& cell) const {
        #if UTOPIA_DEBUG
            if (not cell) {
                throw std::runtime_error(
                    "Checking for Volume on a boundary cell in "
                    "PCPVertex::WorkFunction::VolumeElasticity!"
                );
            }
        #endif

        
        // apical and basal core area
        double surface = 2 * this->_am.area_of(cell);

        // Columnar cells, -> lateral contribution uniform
        if (not this->has_basal_contact(cell)) {
            double height = cell->state.get_parameter(this->_height);
            surface += this->_am.perimeter_of(cell) * height;
            
            return surface;
        }
        // Else might have neighbor contributions
        
        // check neighbours for basal detachment
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            auto [c, n] = this->_am.adjoints_of(edge);
            if (c != cell) { std::swap(c, n); }
            double length = this->_am.length_of(edge);

            // if neighbour has basal contact, does not contribute
            if (n and not this->has_basal_contact(n)) {
                double height = this->height_of(n);

                // apical contribution
                surface += length * height;

                // the new lateral faces 
                double f = this->triangle_factor(n);
                surface += length * std::sqrt(1 + 4. * std::pow(f, 2))
                            * (this->_tissue_height - height);

                // the new medial & basal faces
                surface += std::pow(length, 2) * f;
            }
            else {
                surface += length * this->_tissue_height;
            }
        }

        return surface;
    }

public:
    SurfaceElasticity (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _elastic_modulus(get_as<double>("elastic_modulus", cfg)),
        _preferential_surface(get_as<double>("preferential_surface", cfg))
    { }

    void compute_and_set_forces () override {
        std::unordered_map<std::shared_ptr<Cell>, double> surfaces;
        std::unordered_map<std::shared_ptr<Cell>, double> tensions;
        for (const auto& cell : this->_am.cells()) {
            double S  = surface_of(cell);
            double S0 = _preferential_surface;

            surfaces[cell] = S;
            tensions[cell] = _elastic_modulus * (S - S0) / std::pow(S0, 2);
        }
        std::unordered_map<std::shared_ptr<Edge>, double> lengths;
        for (const auto& edge : this->_am.edges()) {
            lengths[edge] = this->_am.length_of(edge);
        }

        for (const auto& cell : this->_am.cells()) {
            double dW_dH = 0.;
            double dW_dS = tensions[cell];
            
            double H = cell->state.get_parameter(this->_height);
            double f = 0.;  // the triangle factor
            double area = 0.;

            double dWc_dA = dW_dS * 2.;
            double dWn_dA = 0.;     // accumulate from neighbors
            double L2 = 0.;
            double dW_df = 0.;
            // NOTE non-columnar neighbour contributions considered below only
            if (not this->has_basal_contact(cell)) {
                f = this->triangle_factor(cell);
                area = this->_am.area_of(cell);

                dW_dH += dW_dS * this->_am.perimeter_of(cell);

                // neighbour derivatives
                for (const auto& [edge, flip] : cell->custom_links().edges) {
                    auto [cell__, n] = this->_am.adjoints_of(edge);
                    if (cell__ != cell) { std::swap(cell__, n); }

                    // account for A in nbs volume (within factor)
                    if (n and this->has_basal_contact(n)) {
                        double t = tensions[n];
                        double l2 = std::pow(lengths[edge], 2);
                        dW_df += t * (
                            l2 
                            + 4 * lengths[edge] * (this->_tissue_height - H) * f
                              / std::sqrt(1 + 4. * std::pow(f, 2))
                        );
                        L2 += l2;
                    }
                }
                if (L2 > 1.e-8) {
                    dWn_dA = 2 * dW_df / L2;
                }
            }

            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto a = edge->custom_links().a;
                auto b = edge->custom_links().b;
                if (flip) { std::swap(a, b); }

                auto [cell__, n] = this->_am.adjoints_of(edge);
                if (cell__ != cell) { std::swap(cell__, n); }

                SpaceVec displ = this->_am.displacement(a, b);
                SpaceVec normal({displ[1], -displ[0]});
                double l = lengths[edge];

                // calculate the basal triangle contributions
                double dWc_dl = 0.;     // the upper lateral contribution
                double dWn_dl = 0.;     // the lower lateral
                if (not this->has_basal_contact(cell)) {
                    dWc_dl = dW_dS * H;
                    
                    if (n and this->has_basal_contact(n)) {
                        // derivative to S of neighbor n
                        double dW_dSn = tensions[n];

                        // Derivative to l of neighbor's S3 contribution
                        dWn_dl += dW_dSn * (
                            // lateral faces
                            (this->_tissue_height - H) 
                            * std::sqrt(1 + 4. * std::pow(f, 2))
                            // medial & basal faces
                            + 2 * l * f
                        );
                        
                        // factor derivative
                        dWn_dl -= dW_df * 4 * area / std::pow(L2, 2) * l;
                        
                        // derivative S3 to H for neighbor n
                        dW_dH -= dW_dSn * l * std::sqrt(1 + 4. * std::pow(f,2));
                    }
                }
                else {
                    if (not n or this->has_basal_contact(n)) {
                        dWc_dl = dW_dS * this->_tissue_height;
                    }
                    else {
                        dWc_dl = dW_dS * this->height_of(n);
                        n->state.update_parameter(
                            this->_height_derivative,
                            n->state.get_parameter(this->_height_derivative)
                            + dW_dS * l
                        );
                    }
                }

                // calculate forces from perimeter contribution
                a->state.add_force(+ (dWc_dl + dWn_dl) * displ / l);
                b->state.add_force(- (dWc_dl + dWn_dl) * displ / l);

                // calculate forces from apical & basal contribution
                SpaceVec force = -(dWc_dA + dWn_dA) * normal / 2.;
                a->state.add_force(force);
                b->state.add_force(force);
            }

            dW_dH += cell->state.get_parameter(this->_height_derivative);
            cell->state.update_parameter(this->_height_derivative, dW_dH);
        }
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const override {
        double S = surface_of(cell);
        return 0.5*_elastic_modulus * std::pow(S/_preferential_surface - 1., 2);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        const AgentContainer<Cell>& cells
    ) const override
    {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += compute_energy(cell);
        }
        return energy;
    }

    void update_parameters (const DataIO::Config& cfg) override {
        Base::update_parameters(cfg);

        _elastic_modulus = get_as<double>("elastic_modulus", cfg, 
                                          _elastic_modulus);
        _preferential_surface = get_as<double>("preferential_surface", cfg, 
                                               _preferential_surface);
    }

    std::vector<std::string> write_task_cell_properties_names () const override
    {
        return std::vector<std::string>({"surface"});
    }

    std::vector<std::vector<double>> write_cell_properties () const override {
        std::vector<double> surfaces({});

        for (const auto& cell : this->_am.cells()) {
            surfaces.push_back(surface_of(cell));
        }
        return std::vector<std::vector<double>>({
            surfaces
        });
    }
    

    std::vector<std::string> write_task_cell_energies_names () const override {
        return std::vector<std::string>({
            "energy",
        });
    }

    std::vector<std::vector<double>> write_cell_energies () const override {
        std::vector<double> energies({});
        
        for (const auto& cell : this->_am.cells()) {
            energies.push_back(compute_energy(cell));
        }

        return std::vector<std::vector<double>>({
            energies,
        });
    }
};

} // namespace WorkFunction2DplusTriangle
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
