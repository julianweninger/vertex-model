#ifndef UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION2DPLUSTRIANGLE_HH
#define UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION2DPLUSTRIANGLE_HH


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction2Dplus {

/// @brief The volume elasticity of a cell base class
/** \f$ E_\alpha = k/2 (V_\alpha / V^{(0)} - 1)\f$, an elastic penalty on 
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
 *      - `elastic_modulus`: the elastic modulus \f$ k \f$
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
class VolumeElasticityTriangleBase : public WorkFunction::WorkFunctionTerm<Model>
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

    virtual double preferential_volume(const std::shared_ptr<Cell>& cell) const=0;

    /// @brief The height of the tissue, defining a maximum height for every cell
    double _tissue_height;

    /// @brief Minimum value for cell height
    double _minimum_height;

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
                    "PCPVertex::WorkFunction::VolumeElasticity!"
                );
            }
        #endif

        return cell->state.type != 1;
    }

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

    inline double height_of (const std::shared_ptr<Cell>& cell) const {
        #if UTOPIA_DEBUG
            if (not cell) {
                throw std::runtime_error(
                    "Checking for height on a boundary cell in "
                    "PCPVertex::WorkFunction::VolumeElasticity!"
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

        // It is columnar itself, so no neighbour contributions
        if (not has_basal_contact(cell)) {
            return volume;
        }

        // check neighbours for basal detachment
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            auto [c, n] = this->_am.adjoints_of(edge);
            if (c != cell) { std::swap(c, n); }

            // if neighbour has basal contact, does not contribute
            if (n and not has_basal_contact(n)) {
                volume += std::pow(this->_am.length_of(edge), 2)
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
        _minimum_height(get_as<double>("minimum_height", cfg)),
        _height(name + "_" + get_as<std::string>("height_parameter_name", cfg)),
        _height_derivative(_height + "_derivative"),
        _gamma(get_as<double>("gamma", cfg))
    {
        for (const auto& cell : this->_am.cells()) {
            cell->state.register_parameter(_height, _tissue_height);
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
            double V0 = preferential_volume(cell);
            double dH = 0.;     // the derivative to height
            double f = 0.;      // the triangle factor
            if (not has_basal_contact(cell)) {
                dH += _elastic_modulus
                      * (volume_of(cell) - V0) / std::pow(V0, 2)
                      * this->_am.area_of(cell);
                f += triangle_factor(cell);
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


                // consider the neighbours of extruded cell
                if (not has_basal_contact(cell) and n and has_basal_contact(n)) {
                    double V0n = preferential_volume(n);
                    dH -= _elastic_modulus
                          * (volume_of(n) - V0n) / std::pow(V0n, 2)
                          * std::pow(l, 2)
                          * f / 2;

                    // The contribution of perimeter and interfaces
                    double Tn = (
                          _elastic_modulus
                        * (volume_of(n) - V0n) / std::pow(V0n, 2)
                        * l * f * (_tissue_height - H)
                    );

                    a->state.add_force(+ Tn * displ / l);
                    b->state.add_force(- Tn * displ / l);
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

        _minimum_height = get_as<double>("minimum_height",cfg,_minimum_height);
        _elastic_modulus = get_as<double>("elastic_modulus", cfg, 
                                          _elastic_modulus);

        _gamma = get_as<double>("gamma", cfg, _gamma);
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


/// @brief The volume elasticity of a cell
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
            return surface;
        }

        // check neighbours for basal detachment
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            auto [c, n] = this->_am.adjoints_of(edge);
            if (c != cell) { std::swap(c, n); }

            // if neighbour has basal contact, does not contribute
            if (n and not has_basal_contact(n)) {
                double f2 = std::pow(triangle_factor(n),2);
                surface += this->_am.length_of(edge) 
                           * (std::sqrt(1 + 4. * f2) - 1)
                           * (_tissue_height - height_of(n));
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
        _height(get_as<std::string>("VolumeElasticity_term", cfg) 
                + "_" + get_as<std::string>("height_parameter_name", cfg)),
        _height_derivative(_height + "_derivative")
    { }

    void compute_and_set_forces () override {
        for (const auto& cell : this->_am.cells()) {
            const double H = cell->state.get_parameter(_height);
            double surface = surface_of(cell);
            double dH = 0.;
            double f2 = 1.;  // the triangle factor
            if (not has_basal_contact(cell)) {
                dH += _elastic_modulus
                      * (surface - _preferential_surface)
                      / std::pow(_preferential_surface, 2)
                      * this->_am.perimeter_of(cell);
                
                f2 = std::pow(triangle_factor(cell), 2);
                // NOTE triangle factor is assumed constant
            }


            // calculate tension force
            double T = (
                _elastic_modulus * H
                * (surface - _preferential_surface)
                / std::pow(_preferential_surface, 2)
            );
                // calculate pressure force
            double pressure = (
                _elastic_modulus
                * (surface - _preferential_surface) 
                / std::pow(_preferential_surface, 2)
            );
            

            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto [c, n] = this->_am.adjoints_of(edge);
                if (c != cell) { std::swap(c, n); }

                SpaceVec displ = this->_am.displacement(edge);
                double l = arma::norm(displ);

                const auto& a = edge->custom_links().a;
                const auto& b = edge->custom_links().b;

                // calculate forces from tension
                a->state.add_force(+ T * displ / l);
                b->state.add_force(- T * displ / l);

                // calculate forces from pressure
                SpaceVec normal({displ[1], -displ[0]});
                if (flip) {
                    normal *= -1;
                }                
                edge->custom_links().a->state.add_force(- pressure * normal);
                edge->custom_links().b->state.add_force(- pressure * normal);

                // calculate the contributions of 2D+ deformations 
                if (not has_basal_contact(cell) and n and has_basal_contact(n)) {
                    // derivative of H
                    dH -= (
                        _elastic_modulus
                        * (surface_of(n) - _preferential_surface)
                        / std::pow(_preferential_surface, 2)
                        * (std::sqrt(1 + 4. * f2) - 1)
                        * l
                    );

                    // The contribution of perimeter and interfaces
                    double Tn = (
                        _elastic_modulus * (_tissue_height - H)
                        * (surface_of(n) - _preferential_surface)
                        / std::pow(_preferential_surface, 2)
                        * (std::sqrt(1 + 4. * f2) - 1)
                    );

                    a->state.add_force(+ Tn * displ / l);
                    b->state.add_force(- Tn * displ / l);
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
        _preferential_surface = get_as<double>("preferential_surface", cfg, 
                                               _preferential_surface);
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

        _elastic_modulus = get_as<double>("elastic_modulus", cfg, 
                                          _elastic_modulus);
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
