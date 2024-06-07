#ifndef UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION2DPLUSTRIANGLE_HH
#define UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION2DPLUSTRIANGLE_HH


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction2Dplus {

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
class VolumeElasticityTriangle : public WorkFunction::WorkFunctionTerm<Model>
{
protected:
    using Base = WorkFunction::WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

    /// @brief  Elastic constant associated with volume elasticity
    double _elastic_modulus;

    /// @brief The target volume for a cell at which there is no pressure
    double _preferential_volume;

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

    /// @brief The length ("perimeter") of basal contacts
    /// excluding neighbours with basal contacts
    /// @param cell 
    /// @return The number of neighbors that have basal contact
    double normalization_basal_contact (const std::shared_ptr<Cell>& cell) const {
        double length = 0;
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            auto [c, n] = this->_am.adjoints_of(edge);
            if (c != cell) { std::swap(c, n); }

            if (n and has_basal_contact(n)) {
                length += std::pow(this->_am.length_of(edge), 2);
            }
        }
        return length;
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

        const double height = cell->state.get_parameter(_height);
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
                volume += this->_am.area_of(n)
                        * (_tissue_height - n->state.get_parameter(_height))
                        * std::pow(this->_am.length_of(edge), 2)
                        / normalization_basal_contact(n);
            }
        }

        return volume;
    }

public:
    VolumeElasticityTriangle (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _elastic_modulus(get_as<double>("elastic_modulus", cfg)),
        _preferential_volume(get_as<double>("preferential_volume", cfg)),
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

    ~VolumeElasticityTriangle () {
        for (const auto& cell : this->_am.cells()) {
            cell->state.unregister_parameter(_height);
            cell->state.unregister_parameter(_height_derivative);
            cell->state.unregister_parameter(_height_derivative + "_monitor");
        }
    }

    void compute_and_set_forces () final {
        for (const auto& cell : this->_am.cells()) {
            double pressure = compute_pressure(cell);

            const double H = cell->state.get_parameter(_height);
            double area = 0.;
            double dH = 0.;
            double L = 0;
            if (cell->state.type == 1) {
                area += this->_am.area_of(cell);
                dH += _elastic_modulus
                      * (volume_of(cell) - _preferential_volume)
                      / std::pow(_preferential_volume, 2)
                      * area;
                L = normalization_basal_contact(cell);
                if (L < 1.e-8) {
                    throw std::runtime_error("Cell has no basal length!");
                }
            }

            for (const auto& [edge, flip] : cell->custom_links().edges) {
                SpaceVec displ = this->_am.displacement(edge);
                
                SpaceVec normal({displ[1], -displ[0]});
                if (flip) {
                    normal *= -1;
                }

                auto [c, n] = this->_am.adjoints_of(edge);
                if (c != cell) { std::swap(c, n); }

                if (cell->state.type == 1 and n and has_basal_contact(n)) {
                    dH -= _elastic_modulus
                          * (volume_of(n) - _preferential_volume)
                          / std::pow(_preferential_volume, 2)
                          * area
                          * std::pow(arma::norm(displ), 2) / L;
                }

                // calculate force
                SpaceVec force = - pressure / 2. * normal;
                
                edge->custom_links().a->state.add_force(force);
                edge->custom_links().b->state.add_force(force);
            }

            // The contribution of perimeter and interfaces
            if (cell->state.type == 1) {
                for (const auto& [edge, flip] : cell->custom_links().edges) {
                    auto [c, n] = this->_am.adjoints_of(edge);
                    if (c != cell) { std::swap(c, n); }
                    if (not n) { continue; }

                    SpaceVec displ = this->_am.displacement(edge);
                    double l = arma::norm(displ);

                    double T = (
                        _elastic_modulus
                        * 2 * l * (1. - l / L) / L
                        * area * (_tissue_height - H)
                        * (volume_of(n) - _preferential_volume)
                        / std::pow(_preferential_volume, 2)
                    );

                    const auto& a = edge->custom_links().a;
                    const auto& b = edge->custom_links().b;

                    a->state.add_force(+ T * displ / l);
                    b->state.add_force(- T * displ / l);
                }
            }

            if (cell->state.type == 1) {
                dH += cell->state.get_parameter(_height_derivative);
                cell->state.update_parameter(_height_derivative, dH);
            }
        }
    }

    /// @brief The derivative of W to area of this cell
    /// @param cell 
    /// @return apical pressure
    double compute_pressure(const std::shared_ptr<Cell>& cell) const final {
        const double H = cell->state.get_parameter(_height);
        double P = _elastic_modulus * H
                   * (volume_of(cell) - _preferential_volume)
                   / std::pow(_preferential_volume, 2);
        
        // consider how changes in A impact Volume of neighbor
        if (not has_basal_contact(cell)) {
            const double L = normalization_basal_contact(cell);
            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto [c, n] = this->_am.adjoints_of(edge);
                if (c != cell) { std::swap(c, n); }

                if (n and has_basal_contact(n)) {
                    P += _elastic_modulus 
                         * (volume_of(n) - _preferential_volume)
                         / std::pow(_preferential_volume, 2)
                         * (_tissue_height - H)
                         * std::pow(this->_am.length_of(edge), 2) / L;
                }
            }

        }

        return P;
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        double V = volume_of(cell);
        return 0.5*_elastic_modulus * std::pow(V/_preferential_volume - 1., 2);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += compute_energy(cell);
        }
        return energy;
    }
    
    /// Performs the update of \f$ H_\alpha \f$
    void update (double dt) final {
        for (const auto& cell : this->_am.cells()) {
            if (cell->state.type != 1) {
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
        _preferential_volume = get_as<double>("preferential_volume", cfg, 
                                              _preferential_volume);
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
        if (cfg["set_height"]) {
            auto height = get_as<double>("set_height", cfg);
            this->_am.get_logger()->debug("Setting height of cells with type 1 "
                "to {}", height);
            for (const auto& cell : this->_am.cells()) {
                cell->state.update_parameter(_height, height);
            }
        }

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
        bool test_volume = fabs(tot_volume - area * _tissue_height) < 1.e-8;
        if (test_volume) {
            logger->debug("   Volume constraint is fulfilled.");
        }
        else {
            PASS_TEST = false;
            logger->error("  Test of volume constraint failed! "
                          "Total volume was {}, expected volume was {}. "
                          "Difference exceeds {}.",
                          tot_volume, area * _tissue_height, 1.e-8);
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

    std::vector<std::string> write_task_cell_properties_names () const final {
        return std::vector<std::string>({"height", "volume"});
    }

    std::vector<std::vector<double>> write_cell_properties () const final {
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
    

    std::vector<std::string> write_task_cell_energies_names () const final {
        return std::vector<std::string>({
            "energy",
            "pressure",
            "extrusiveness"
        });
    }

    std::vector<std::vector<double>> write_cell_energies () const final {
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

template <typename Model>
class SurfaceElasticity : public WorkFunction::WorkFunctionTerm<Model>
{
protected:
    using Base = WorkFunction::WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

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

    const std::string _triangle_factor;
    
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

    /// @brief The length ("perimeter") of basal contacts
    /// excluding neighbours with basal contacts
    /// @param cell 
    /// @return The number of neighbors that have basal contact
    double normalization_basal_contact (const std::shared_ptr<Cell>& cell) const {
        double length = 0;
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            auto [c, n] = this->_am.adjoints_of(edge);
            if (c != cell) { std::swap(c, n); }

            if (n and has_basal_contact(n)) {
                length += std::pow(this->_am.length_of(edge), 2);
            }
        }
        return length;
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
                double f2 = std::pow(n->state.get_parameter(_triangle_factor),2);
                surface += this->_am.length_of(edge) 
                           * (std::sqrt(1 + 3 * f2) - 1)
                           * (_tissue_height - n->state.get_parameter(_height));
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
        _height_derivative(_height + "_derivative"),
        _triangle_factor(name + "_triangle_factor")
    {
        for (const auto& cell : this->_am.cells()) {
            cell->state.register_parameter(_triangle_factor, 1.);
        }
    }

    ~SurfaceElasticity () {
        for (const auto& cell : this->_am.cells()) {
            cell->state.unregister_parameter(_triangle_factor);
        }
    }

    void compute_and_set_forces () final {
        for (const auto& cell : this->_am.cells()) {
            const double H = cell->state.get_parameter(_height);
            double surface = surface_of(cell);
            double dH = 0.;
            if (not has_basal_contact(cell)) {
                dH += _elastic_modulus
                      * (surface - _preferential_surface)
                      / std::pow(_preferential_surface, 2)
                      * this->_am.perimeter_of(cell);
            }


            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto [c, n] = this->_am.adjoints_of(edge);
                if (c != cell) { std::swap(c, n); }

                SpaceVec displ = this->_am.displacement(edge);
                double l = arma::norm(displ);

                // calculate tension force
                double T = (
                    _elastic_modulus * H
                    * (surface - _preferential_surface)
                    / std::pow(_preferential_surface, 2)
                );

                const auto& a = edge->custom_links().a;
                const auto& b = edge->custom_links().b;

                a->state.add_force(+ T * displ / l);
                b->state.add_force(- T * displ / l);

                // calculate pressure force
                SpaceVec normal({displ[1], -displ[0]});
                if (flip) {
                    normal *= -1;
                }
                SpaceVec force = - _elastic_modulus
                                 * (surface - _preferential_surface) 
                                 / std::pow(_preferential_surface, 2)
                                 * normal;
                
                edge->custom_links().a->state.add_force(force);
                edge->custom_links().b->state.add_force(force);

                if (not has_basal_contact(cell) and n and has_basal_contact(n)) {
                    double f2 = std::pow(
                        n->state.get_parameter(_triangle_factor), 2
                    );
                    dH -= _elastic_modulus
                          * (surface_of(n) - _preferential_surface)
                          / std::pow(_preferential_surface, 2)
                          * (std::sqrt(1 + 3 * f2) - 1)
                          * l;
                }
            }

            // The contribution of perimeter and interfaces
            if (not has_basal_contact(cell)) {
                for (const auto& [edge, flip] : cell->custom_links().edges) {
                    auto [c, n] = this->_am.adjoints_of(edge);
                    if (c != cell) { std::swap(c, n); }
                    if (not n) { continue; }

                    SpaceVec displ = this->_am.displacement(edge);
                    double l = arma::norm(displ);

                    double f2 = std::pow(
                        n->state.get_parameter(_triangle_factor), 2
                    );
                    double T = (
                        _elastic_modulus * (_tissue_height - H)
                        * (surface_of(n) - _preferential_surface)
                        / std::pow(_preferential_surface, 2)
                        * (std::sqrt(1 + 3 * f2) - 1)
                    );

                    const auto& a = edge->custom_links().a;
                    const auto& b = edge->custom_links().b;

                    a->state.add_force(+ T * displ / l);
                    b->state.add_force(- T * displ / l);
                }

                dH += cell->state.get_parameter(_height_derivative);
                cell->state.update_parameter(_height_derivative, dH);
            }

            // NOTE triangle factor is assumed constant
        }
    }

    double compute_energy(const std::shared_ptr<Cell>& cell) const final {
        double S = surface_of(cell);
        return 0.5*_elastic_modulus * std::pow(S/_preferential_surface - 1., 2);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& cell : cells) {
            energy += compute_energy(cell);
        }
        return energy;
    }
    
    /// Performs the update of \f$ H_\alpha \f$
    void update (double dt) final {
        for (const auto& cell : this->_am.cells()) {
            if (cell->state.type != 1) {
                continue;
            }
            cell->state.update_parameter(
                _triangle_factor,
                triangle_factor(cell)
            );
        }
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

    std::vector<std::string> write_task_cell_properties_names () const final {
        return std::vector<std::string>({"surface"});
    }

    std::vector<std::vector<double>> write_cell_properties () const final {
        std::vector<double> surfaces({});

        for (const auto& cell : this->_am.cells()) {
            surfaces.push_back(surface_of(cell));
        }
        return std::vector<std::vector<double>>({
            surfaces
        });
    }
    

    std::vector<std::string> write_task_cell_energies_names () const final {
        return std::vector<std::string>({
            "energy",
        });
    }

    std::vector<std::vector<double>> write_cell_energies () const final {
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
