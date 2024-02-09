#ifndef UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION2DPLUS_HH
#define UTOPIA_MODELS_PCPVERTEX_WORKFUNCTION2DPLUS_HH


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction2Dplus {

/// @brief The volume elasticity of a cell
/** TODO */
template <typename Model>
class VolumeElasticity : public WorkFunction::WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunction::WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;
    
    /// @brief Check whether cell has a basal contact
    /// @param height The height of considered cell
    /// @return Whether cell has a basal contact
    bool has_basal_contact(const std::shared_ptr<Cell>& cell) const {
        return cell->state.type != 1;
    }

    inline double height_of (const std::shared_ptr<Cell>& cell) const {
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
        const double height = cell->state.get_parameter(_height);
        const double area = this->_am.area_of(cell);

        // It is columnar itself, so no neighbour contributions
        if (not has_basal_contact(cell)) {
            return area * height;
        }

        // check neighbours for basal detachment
        double volume = area * height;
        for (const auto& n : this->_am.neighbors_of(cell)) {
            // if neighbour has basal contact, does not contribute
            if (has_basal_contact(n)) {
                continue;
            }
            volume += this->_am.area_of(n)
                      * (_tissue_height - n->state.get_parameter(_height)) 
                      / num_basal_neighbors(n);
        }

        return volume;
    }

protected:
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

    /// @brief The number of neighbors that have basal contact
    /// @param cell 
    /// @return The number of neighbors that have basal contact
    std::size_t num_basal_neighbors (const std::shared_ptr<Cell>& cell) const {
        std::size_t cnt = 0;
        for (const auto& n : this->_am.neighbors_of(cell)) {
            cnt += has_basal_contact(n);
        }
        return cnt;
    }

public:
    VolumeElasticity (
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

    ~VolumeElasticity () {
        for (const auto& cell : this->_am.cells()) {
            cell->state.unregister_parameter(_height);
            cell->state.unregister_parameter(_height_derivative);
            cell->state.unregister_parameter(_height_derivative + "_monitor");
        }
    }

    void compute_and_set_forces () final {
        for (const auto& cell : this->_am.cells()) {
            double pressure = compute_pressure(cell);

            double area = 0.;
            double dH = 0.;
            std::size_t N_div = 0;
            if (cell->state.type == 1) {
                area += this->_am.area_of(cell);
                dH += _elastic_modulus
                      * (volume_of(cell) - _preferential_volume)
                      / std::pow(_preferential_volume, 2)
                      * area;
                N_div = num_basal_neighbors(cell);
                if (N_div == 0) {
                    dH -= 1.e4;
                }
            }

            for (const auto& [edge, flip] : cell->custom_links().edges) {
                SpaceVec displ = this->_am.displacement(edge);
                auto [c, n] = this->_am.template adjoints_of<true>(edge);
                
                SpaceVec normal({displ[1], -displ[0]});
                if (flip) {
                    normal *= -1;
                    std::swap(c, n);
                }

                if (cell->state.type == 1 and has_basal_contact(n)) {
                    dH -= _elastic_modulus
                          * (volume_of(n) - _preferential_volume)
                          / std::pow(_preferential_volume, 2)
                          * area / N_div;
                }

                // calculate force
                SpaceVec force = - pressure / 2. * normal;
                
                edge->custom_links().a->state.add_force(force);
                edge->custom_links().b->state.add_force(force);
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
            std::size_t N_div = num_basal_neighbors(cell);
            for (const auto& [edge, flip] : cell->custom_links().edges) {
                auto [c, n] = this->_am.template adjoints_of<true>(edge);                
                if (flip) {
                    std::swap(c, n);
                }
                if (has_basal_contact(n)) {
                    P += _elastic_modulus 
                         * (volume_of(n) - _preferential_volume) 
                         / std::pow(_preferential_volume, 2) 
                         * (_tissue_height - H) / N_div;
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

    bool test_constraints (const std::shared_ptr<spdlog::logger>& logger) const override {
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
        const std::size_t N = this->_am.cells().size();
        bool test_volume = fabs(tot_volume - N * _preferential_volume) < 1.e-8;
        if (test_volume) {
            logger->debug("   Volume constraint is fulfilled.");
        }
        else {
            PASS_TEST = false;
            logger->error("  Test of volume constraint failed! "
                          "Total volume was {}, expected volume was {}. "
                          "Difference exceeds {}.",
                          tot_volume, N * _preferential_volume, 1.e-8);
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


/// @brief The linetension term
/** \f$ E_{i,j} = k l_{i,j}\f$, a term linear in edge length \f$ l \f$.
 * 
 *  Parameters:
 *      - `linetension`: the tension \f$ k \f$
 */
template <typename Model>
class SurfaceTension : public WorkFunction::WorkFunctionTerm<Model>
{
    using Base = WorkFunction::WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

protected:
    double _surface_tension;

    const std::string _cell_height;

    const std::string _cell_height_derivative;

    double apical_height_of (const std::shared_ptr<Edge>& edge) const {
        const auto& [ca, cb] = this->_am.adjoints_of(edge);
        return std::min(
            ca->state.get_parameter(_cell_height),
            cb->state.get_parameter(_cell_height)
        );
    }

public:
    SurfaceTension (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _surface_tension(get_as<double>("surface_tension", cfg)),
        _cell_height(get_as<std::string>("VolumeElasticity_term", cfg) 
                     + "_"
                     + get_as<std::string>("height_parameter_name", cfg)),
        _cell_height_derivative(_cell_height + "_derivative")
    { }

    void compute_and_set_forces () final {
        for (const auto& edge : this->_am.edges()) {
            SpaceVec displ = this->_am.displacement(edge);
            SpaceVec director = displ / arma::norm(displ);
            double H = apical_height_of(edge);

            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;

            a->state.add_force(+ _surface_tension * director * H);
            b->state.add_force(- _surface_tension * director * H);

            auto [ca, cb] = this->_am.adjoints_of(edge);
            if (cb->state.type == 1) { std::swap(ca, cb); }
            if (ca->state.type == 1) {
                auto dH = ca->state.get_parameter(_cell_height_derivative);
                dH += _surface_tension * arma::norm(displ);
                ca->state.update_parameter(_cell_height_derivative, dH);
            }
        }
    }

    SpaceVec compute_force(const std::shared_ptr<Vertex>& vertex) const final {
        SpaceVec force({0., 0.});
        for (const auto& edge : this->_am.adjoint_edges_of(vertex)) {
            SpaceVec displ = this->_am.displacement(edge);
            SpaceVec director = displ / arma::norm(displ);
            double H = apical_height_of(edge);

            if (edge->custom_links().a == vertex) {
                force += _surface_tension * director * H;
            }
            else {
                force -= _surface_tension * director * H;
            }
        }
        return force;
    }

    double compute_energy(const std::shared_ptr<Edge>& edge) const final {
        return _surface_tension * this->_am.length_of(edge)
               * apical_height_of(edge);
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& edge : edges) {
            energy += compute_energy(edge);
        }
        return energy;
    }

    void update_parameters (const DataIO::Config& cfg) final {
        _surface_tension  = get_as<double>("surface_tension", cfg, 
                                           _surface_tension);
    }

    std::vector<std::string> write_task_edge_properties_names () const final {
        return std::vector<std::string>({
            "apical_height",
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const final {
        std::vector<double> apical_heights({});
        
        for (const auto& edge : this->_am.edges()) {
            apical_heights.push_back(apical_height_of(edge));
        }

        return std::vector<std::vector<double>>({
            apical_heights
        });
    }

    std::vector<std::string> write_task_edge_energies_names () const final {
        return std::vector<std::string>({
            "energy",
        });
    }

    std::vector<std::vector<double>> write_edge_energies () const final {
        std::vector<double> energies({});
        
        for (const auto& edge : this->_am.edges()) {
            energies.push_back(compute_energy(edge));
        }

        return std::vector<std::vector<double>>({
            energies,
        });
    }
};


} // namespace WorkFunction2Dplus
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
