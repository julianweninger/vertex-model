#ifndef UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION2DPLUSTRIANGLE_HH
#define UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION2DPLUSTRIANGLE_HH


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction2Dplus {

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
class VolumeElasticityTriangleBase 
: 
    public WorkFunction::WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunction::WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// @brief  Elastic constant associated with volume elasticity
    double _elastic_modulus;

    virtual double preferential_volume(
            const std::shared_ptr<Cell>& cell) const=0;

    /// @brief The height of the tissue, defining a maximum height for every cell
    double _tissue_height;

    /// @brief Minimum value for cell height
    double _minimum_height;

    /// The crossover hight for loss of basal contact
    double _critical_height;

    /// @brief  The name of the cell's parameter referring to cell's height
    const std::string _height;

    /// @brief  The name of the cell's parameter referring to cell's height
    const std::string _height_derivative;

    /// @brief  Damping factor for hight development additional to dt
    double _gamma;
    
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
            l2 += std::pow(this->_am.length_of(edge), 2);
        }

        return 2 * this->_am.area_of(cell) / l2;
    }

    /// @brief A factor for crossover from attached to detached state
    /** A crossover between full basal attachment to basal detachment occurs
     *  at `_critical_height`. It is implemented as a sigmoidal factor, where
     *  depth of triangles is rescales with `distribution_factor`.
     */
    double distribution_factor (const std::shared_ptr<Cell>& cell) const {
        #if UTOPIA_DEBUG
            if (not cell) {
                throw std::runtime_error(
                    "Checking for distribution_factor on a boundary cell in "
                    "PCPVertex::WorkFunction::VolumeElasticityTriangle!"
                );
            }
        #endif

        double H = this->height_of(cell);
        double k = 5. / (_tissue_height - _critical_height);
        return 1. - 1. / (1. + exp(- k * (H - _critical_height)));
    }

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

        const double height = height_of(cell);
        const double area = this->_am.area_of(cell);

        double volume = area * height;

        // It is columnar itself, no neighbour contributions
        // but consider triangles not contributed by neighbors
        if (not has_basal_contact(cell)) {
            double f = triangle_factor(cell);
            double d = distribution_factor(cell);

            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto [c, n] = this->_am.adjoints_of(edge);
                if (c != cell) { std::swap(c, n); }

                // the triangle volume
                double V3 = std::pow(this->_am.length_of(edge), 2)
                            * (_tissue_height - height)
                            * f / 2;

                // the distribution
                if (n and has_basal_contact(n)) {
                    volume += (1 - d) * V3;
                }
                else {
                    volume += V3;
                }
            }
            return volume;
        }

        // check neighbours for basal detachment
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            auto [c, n] = this->_am.adjoints_of(edge);
            if (c != cell) { std::swap(c, n); }

            // if neighbour has basal contact, does not contribute
            if (n and not has_basal_contact(n)) { 
                volume += distribution_factor(n)
                        * std::pow(this->_am.length_of(edge), 2)
                        * (_tissue_height - height_of(n))
                        * triangle_factor(n) / 2;
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
        _elastic_modulus(get_as<double>("elastic_modulus", cfg)),
        _tissue_height(get_as<double>("tissue_height", cfg)),
        _critical_height(get_as<double>("critical_height", cfg)),
        _minimum_height(get_as<double>("minimum_height", cfg)),
        _height(name + "__" + get_as<std::string>("height_parameter_name", cfg)),
        _height_derivative(_height + "_derivative"),
        _gamma(get_as<double>("gamma", cfg))
    {
        auto height = get_as<double>("set_height", cfg, _tissue_height);
        for (const auto& cell : this->_am.cells()) {
            if (cell->state.type != 1) {
                cell->state.register_parameter(_height, _tissue_height);
            }
            else {
                cell->state.register_parameter(_height, height);
            }
            cell->state.register_parameter(_height_derivative, 0.);
            cell->state.register_parameter(_height_derivative + "_monitor", 0.);
        }
    }

    ~VolumeElasticityTriangleBase () {
        for (const auto& cell : this->_am.cells()) {
            cell->state.unregister_parameter(_height);
            cell->state.unregister_parameter(_height_derivative);
            cell->state.unregister_parameter(_height_derivative + "_monitor");
        }
    }

    void compute_and_set_forces () override {
        for (const auto& cell : this->_am.cells()) {
            double pressure = compute_pressure(cell);

            const double H = cell->state.get_parameter(_height);
            double volume = 0.;
            double V0 = preferential_volume(cell);
            double dH = 0.;     // the derivative to height
            double f = 0.;      // the triangle factor
            double d = 0.;      // the distribution factor
            double dW_dV = 0.;
            if (not has_basal_contact(cell)) {
                volume = volume_of(cell);
                f = triangle_factor(cell);
                d = distribution_factor(cell);

                dW_dV = _elastic_modulus * (volume - V0) / std::pow(V0, 2);
                dH += dW_dV * this->_am.area_of(cell);
            }

            // apply forces via edges
            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto [c, n] = this->_am.adjoints_of(edge);
                if (c != cell) { std::swap(c, n); }

                const auto& a = edge->custom_links().a;
                const auto& b = edge->custom_links().b;

                SpaceVec displ = this->_am.displacement(edge);
                double l = arma::norm(displ);
                
                SpaceVec normal({displ[1], -displ[0]});
                if (flip) {
                    normal *= -1;
                }

                // calculate force from pressure orthogonal to junction
                SpaceVec force = - pressure / 2. * normal;
                a->state.add_force(force);
                b->state.add_force(force);

                // consider the basal triangle volumes
                if (not has_basal_contact(cell)) {
                    if (n and has_basal_contact(n)) {
                        double V0n = preferential_volume(n);
                        double Vn = volume_of(n);
                        double dW_dVn = _elastic_modulus * (Vn - V0n)
                                                         / std::pow(V0n, 2);

                        // The own and neighbor contributions to tension
                        double Tn = (dW_dV * (1-d) + dW_dVn * d)
                                    * l * f * (_tissue_height - H);

                        a->state.add_force(+ Tn * displ / l);
                        b->state.add_force(- Tn * displ / l);

                        dH -= (dW_dV * (1-d) + dW_dVn * d)
                              * std::pow(l, 2) * f / 2;
                    }
                    else {
                        // The contribution of perimeter and interfaces
                        double T = dW_dV * l * f * (_tissue_height - H);

                        a->state.add_force(+ T * displ / l);
                        b->state.add_force(- T * displ / l);

                        dH -= dW_dV * std::pow(l, 2) * f / 2;
                    }
                }
            }

            cell->state.update_parameter(
                _height_derivative, 
                cell->state.get_parameter(_height_derivative) + dH
            );
        }
    }

    /// @brief The derivative of W to area of this cell
    /// @param cell 
    /// @return apical pressure
    double compute_pressure(const std::shared_ptr<Cell>& cell) const override {
        const double H = cell->state.get_parameter(_height);
        double V0 = preferential_volume(cell);
        double P = _elastic_modulus * H
                   * (volume_of(cell) - V0) / std::pow(V0, 2);

        return P;
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
    
    /// Performs the update of \f$ H_\alpha \f$
    void update (double dt) override {
        for (const auto& cell : this->_am.cells()) {
            if (has_basal_contact(cell)) {
                continue;
            }

            double H = cell->state.get_parameter(_height);
            double dH = cell->state.get_parameter(_height_derivative);
            cell->state.update_parameter(
                _height,
                std::max(
                    std::min(H - dt * _gamma * dH, _tissue_height), 
                    _minimum_height
                )
            );
            cell->state.update_parameter(
                _height_derivative + "_monitor", 
                cell->state.get_parameter(_height_derivative)
            );
            cell->state.update_parameter(_height_derivative, 0.);
        }
    }

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

        _critical_height = get_as<double>("critical_height", cfg,
                                          _critical_height);
        _minimum_height = get_as<double>("minimum_height",cfg,_minimum_height);
        _elastic_modulus = get_as<double>("elastic_modulus", cfg, 
                                          _elastic_modulus);

        _gamma = get_as<double>("gamma", cfg, _gamma);
    }

    bool test_constraints (const std::shared_ptr<spdlog::logger>& logger) 
    const override 
    {
        bool PASS_TEST = Base::test_constraints(logger);
        logger->debug("Testing constraints of WF term {} ...", this->_name);


        logger->debug("   Testing height constraint ...");
        bool test_height = true;
        for (const auto& cell : this->_am.cells()) {
            double height = cell->state.get_parameter(_height);
            if (height > _tissue_height + 1.e-8) {
                test_height = false;
                logger->error("Cell {} (height: {}) exceeds the tissue height "
                              "({})!",
                              cell->id(), height, _tissue_height);
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
        bool test_volume = fabs(tot_volume - area * _tissue_height) / tot_volume < 1.e-6;
        if (test_volume) {
            logger->debug("   Volume constraint is fulfilled.");
        }
        else {
            PASS_TEST = false;
            logger->error("  Test of volume constraint failed! "
                          "Total volume was {}, expected volume was {}. "
                          "Relative error {} exceeds {}.",
                          tot_volume, area * _tissue_height, 
                          fabs(tot_volume - area * _tissue_height) / tot_volume,
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
        return std::vector<std::string>({"height", "volume"});
    }

    std::vector<std::vector<double>> write_cell_properties () const override {
        std::vector<double> heights({});
        std::vector<double> volumes({});

        for (const auto& cell : this->_am.cells()) {
            heights.push_back(cell->state.get_parameter(_height));
            volumes.push_back(volume_of(cell));
        }
        return std::vector<std::vector<double>>({
            heights, volumes
        });
    }
    

    std::vector<std::string> write_task_cell_energies_names () const override {
        return std::vector<std::string>({
            "energy",
            "pressure",
            "extrusiveness"
        });
    }

    std::vector<std::vector<double>> write_cell_energies () const override {
        std::vector<double> energies({});
        std::vector<double> pressures({});
        std::vector<double> dHs({});
        
        for (const auto& cell : this->_am.cells()) {
            energies.push_back(compute_energy(cell));
            pressures.push_back(compute_pressure(cell));
            dHs.push_back(-1. * cell->state.get_parameter(_height_derivative + "_monitor"));
        }

        return std::vector<std::vector<double>>({
            energies,
            pressures,
            dHs
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
class VolumeElasticityTriangleHeterotypic : public VolumeElasticityTriangleBase<Model>
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
class SurfaceElasticity : public WorkFunction::WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunction::WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    /// @brief  Elastic constant associated with volume elasticity
    double _elastic_modulus;

    /// @brief The target volume for a cell at which there is no pressure
    double _preferential_surface;

    /// @brief The height of the tissue, defining a maximum height for every cell
    double _tissue_height;

    /// The crossover hight for loss of basal contact
    double _critical_height;

    /// @brief  The name of the cell's parameter referring to cell's height
    const std::string _height;

    /// @brief  The name of the cell's parameter referring to cell's height
    const std::string _height_derivative;

    /// @brief Check whether cell has a basal contact
    /// @param height The height of considered cell
    /// @return Whether cell has a basal contact
    bool has_basal_contact(const std::shared_ptr<Cell>& cell) const {
        #if UTOPIA_DEBUG
            if (not cell) {
                throw std::runtime_error(
                    "Checking for Basal contact on a boundary cell in "
                    "PCPVertex::WorkFunction::SurfaceElasticity!"
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
                    "PCPVertex::WorkFunction::SurfaceElasticity!"
                );
            }
        #endif

        double l2 = 0.;
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            l2 += std::pow(this->_am.length_of(edge), 2);
        }

        return 2 * this->_am.area_of(cell) / l2;
    }

    /// @brief A factor for crossover from attached to detached state
    /** A crossover between full basal attachment to basal detachment occurs
     *  at `_critical_height`. It is implemented as a sigmoidal factor, where
     *  depth of triangles is rescales with `distribution_factor`.
     */
    double distribution_factor (const std::shared_ptr<Cell>& cell) const {
        #if UTOPIA_DEBUG
            if (not cell) {
                throw std::runtime_error(
                    "Checking for distribution_factor on a boundary cell in "
                    "PCPVertex::WorkFunction::VolumeElasticityTriangle!"
                );
            }
        #endif

        // double k = 5. / (_tissue_height - _critical_height);
        // return 1. - 1. / (1. + exp(- k * (0.98 - _critical_height)));

        double H = this->height_of(cell);
        double k = 5. / (_tissue_height - _critical_height);
        return 1. - 1. / (1. + exp(- k * (H - _critical_height)));
    }

    inline double height_of (const std::shared_ptr<Cell>& cell) const {
        #if UTOPIA_DEBUG
            if (not cell) {
                throw std::runtime_error(
                    "Checking for height on a boundary cell in "
                    "PCPVertex::WorkFunction::SurfaceElasticity!"
                );
            }
        #endif

        return cell->state.get_parameter(_height);
    }

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

        const double height = cell->state.get_parameter(_height);

        double surface = this->_am.perimeter_of(cell) * height;
        surface += 2 * this->_am.area_of(cell);

        // It is columnar itself, so no neighbour contributions
        if (not has_basal_contact(cell)) {
            double f = triangle_factor(cell);
            double d = distribution_factor(cell);

            // consider triangles that are not occupied by neighbors
            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto [c, n] = this->_am.adjoints_of(edge);
                if (c != cell) { std::swap(c, n); }

                double length = this->_am.length_of(edge);
                if (n and has_basal_contact(n)) {
                    // the new lateral faces
                    double S3 = length * std::sqrt(1 + 4. * std::pow(f*d, 2))
                                * (_tissue_height - height);

                    // additional rescaling subtracting approximated
                    // triangle-triangle interface
                    surface += (1. - d) * S3;
                }
                // no neighbor to occupy triangle
                else {
                    // the new lateral faces, including the outer extension
                    surface += length * (std::sqrt(1 + 4. * std::pow(f,2)) + 1)
                               * (_tissue_height - height);
                }
            }
            
            return surface;
        }

        // check neighbours for basal detachment
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            auto [c, n] = this->_am.adjoints_of(edge);
            if (c != cell) { std::swap(c, n); }

            // if neighbour has basal contact, does not contribute
            if (n and not has_basal_contact(n)) {
                double length = this->_am.length_of(edge);
                double f = triangle_factor(n);
                double d = distribution_factor(n);

                // the new lateral faces 
                // minus what was assumed from columnar
                double S3 = length * (std::sqrt(1 + 4. * std::pow(f*d, 2)) - 1)
                            * (_tissue_height - height_of(n));

                surface += S3;

                // the new medial & basal faces
                surface += std::pow(length, 2) * f * d;
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
        _preferential_surface(get_as<double>("preferential_surface", cfg)),
        _tissue_height(get_as<double>("tissue_height", cfg)),
        _critical_height(get_as<double>("critical_height", cfg)),
        _height(get_as<std::string>("VolumeElasticity_term", cfg) 
                + "__" + get_as<std::string>("height_parameter_name", cfg)),
        _height_derivative(_height + "_derivative")
    { }

    void compute_and_set_forces () override {
        for (const auto& cell : this->_am.cells()) {
            const double H = cell->state.get_parameter(_height);
            double surface = surface_of(cell);
            double area = 0.;
            double dH = 0.;
            double f = 1.;  // the triangle factor
            double d = 1.;  // the distribution factor

            // derivative to S
            double dW_dS = (
                _elastic_modulus
                * (surface - _preferential_surface) 
                / std::pow(_preferential_surface, 2)
            );

            // NOTE non-columnar neighbour contributions considered below only
            if (not has_basal_contact(cell)) {
                area = this->_am.area_of(cell);

                dH += dW_dS * this->_am.perimeter_of(cell);
                
                f = triangle_factor(cell);
                d = distribution_factor(cell);
                // NOTE factors are assumed constant
                //      i.e. df/dl = 0 and df/dA = 0
                //           dd/dl = 0 and dd/dA = 0
            }

            // derivative perimeter contribution
            double T = dW_dS * H;

            // derivative apical & basal contributions
            double pressure = dW_dS * 2;


            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto [c, n] = this->_am.adjoints_of(edge);
                if (c != cell) { std::swap(c, n); }

                SpaceVec displ = this->_am.displacement(edge);
                double l = arma::norm(displ);

                const auto& a = edge->custom_links().a;
                const auto& b = edge->custom_links().b;

                // calculate forces from perimeter contribution
                a->state.add_force(+ T * displ / l);
                b->state.add_force(- T * displ / l);

                // calculate forces from apical & basal contribution
                SpaceVec normal({displ[1], -displ[0]});
                if (flip) {
                    normal *= -1;
                }                
                a->state.add_force(- pressure/2. * normal);
                b->state.add_force(- pressure/2. * normal);

                // calculate the basal triangle contributions
                if (not has_basal_contact(cell)) {
                    if (n and has_basal_contact(n)) {
                        // Derivative to l of cell's S3 contribution
                        double Tn = dW_dS * (1. - d) * (
                            // lateral faces of triangle
                            (_tissue_height - H) 
                            * std::sqrt(1 + 4. * std::pow(f * d, 2))
                        );

                        // derivative to S of neighbor n
                        double dW_dS_n = (
                            _elastic_modulus
                            * (surface_of(n) - _preferential_surface)
                            / std::pow(_preferential_surface, 2)
                        );

                        // Derivative to l of neighbor's S3 contribution
                        Tn += dW_dS_n * (
                            // lateral faces
                            // minus what was assumed from columnar
                            (_tissue_height - H) 
                            * (std::sqrt(1 + 4. * std::pow(f * d, 2)) - 1)
                            // medial & basal faces
                            + 2 * l * f * d
                        );

                        // Force from neighbor's S3 contribution
                        a->state.add_force(+ Tn * displ / l);
                        b->state.add_force(- Tn * displ / l);
                        
                        // derivative S3 to H for cell
                        dH -= dW_dS * l * (1 - d)
                              * std::sqrt(1 + 4. * std::pow(f * d, 2));

                        // derivative S3 to H for cell
                        // the change in distribution factor contribution
                        double k = 5. / (_tissue_height - _critical_height);
                        double expo = exp(-k * (H - _critical_height));
                        dH += dW_dS * l * (_tissue_height - H)
                              * std::sqrt(1 + 4. * std::pow(f * d, 2))
                              * std::pow(1./(1+expo), 2) * k * expo;
                        
                        // derivative S3 to H for neighbor n
                        dH -= dW_dS_n * l 
                              * (std::sqrt(1 + 4. * std::pow(f * d, 2)) - 1);
                        dH -= dW_dS_n * std::pow(l, 2) * f
                              * std::pow(1./(1+expo), 2) * k * expo;
                    }
                    else {
                        // derivative to l of unoccupied triangles
                        double T = dW_dS * (
                            (_tissue_height - H) 
                            * (std::sqrt(1 + 4. * std::pow(f, 2)) + 1)
                        );

                        a->state.add_force(+ T * displ / l);
                        b->state.add_force(- T * displ / l);

                        // derivative to H of unoccupied triangles
                        dH -= dW_dS * (std::sqrt(1 + 4 * std::pow(f, 2)) + 1)*l;
                    }
                }
            }

            dH += cell->state.get_parameter(_height_derivative);
            cell->state.update_parameter(_height_derivative, dH);
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
        _elastic_modulus = get_as<double>("elastic_modulus", cfg, 
                                          _elastic_modulus);
        _preferential_surface = get_as<double>("preferential_surface", cfg, 
                                               _preferential_surface);
        
        // update tissue height within all cells
        double tmp = get_as<double>("tissue_height", cfg, _tissue_height);
        if (fabs(tmp - _tissue_height) > 1.e-8) {
            _tissue_height = tmp;
            this->_am.get_logger()->debug("Setting height of cells excluding "
                "type 1 to {}", _tissue_height);
            for (const auto& cell : this->_am.cells()) {
                if (cell->state.type != 1) {
                    cell->state.update_parameter(_height, _tissue_height);

                }
            }
        }

        _critical_height = get_as<double>("critical_height", cfg,
                                          _critical_height);
    }

    std::vector<std::string> write_task_cell_properties_names () const override {
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
