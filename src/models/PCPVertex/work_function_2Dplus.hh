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
    inline bool has_basal_contact(const double& height) const {
        return fabs(height - _tissue_height) < 1.e-8;
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
        const double& height = cell->state.get_parameter(_height);
        const double area = this->_am.area_of(cell);

        // It is columnar itself, so no neighbour contributions
        if (not has_basal_contact(height)) {
            return area * height;
        }

        // check neighbours for basal detachment
        double volume = area * height;
        for (const auto& n : this->_am.neighbors_of(cell)) {
            const double& h = n->state.get_parameter(_height);
            // if neighbour has basal contact, does not contribute
            if (has_basal_contact(h)) {
                continue;
            }
            const double a = this->_am.area_of(n);
            volume += a * (1.-h) / this->_am.neighbors_of(n).size();
        }

        return volume;
    }

protected:
    /// @brief The target volume for a cell at which there is no pressure
    double _preferential_volume;

    /// @brief The height of the tissue, defining a maximum height for every cell
    double _tissue_height;

    /// @brief  The name of the cell's parameter referring to cell's height
    const std::string _height;

    /// @brief  Elastic constant associated with volume elasticity
    double _elastic_modulus;

public:
    VolumeElasticity (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _preferential_volume(get_as<double>("preferential_volume", cfg)),
        _tissue_height(get_as<double>("tissue_height", cfg)),
        _height(name + "_height"),
        _elastic_modulus(get_as<double>("elastic_modulus", cfg))
    {
        for (const auto& cell : this->_am.cells()) {
            cell->state.register_parameter(_height, _tissue_height);
        }
    }

    ~VolumeElasticity () {
        for (const auto& cell : this->_am.cells()) {
            cell->state.unregister_parameter(_height);
        }
    }

    void compute_and_set_forces () final {
        for (const auto& cell : this->_am.cells()) {
            const double& H = cell->state.get_parameter(_height);

            double pressure = compute_pressure(cell);
            for (const auto& [edge, flip] : cell->custom_links().edges) {
                SpaceVec displ = this->_am.displacement(edge);
                auto [c, n] = this->_am.template adjoints_of<true>(edge);
                
                SpaceVec normal({displ[1], -displ[0]});
                if (flip) {
                    normal *= -1;
                    std::swap(c, n);
                }

                double factor = H;
                double h = n->state.get_parameter(_height);
                if (not has_basal_contact(h)) {
                    double p = compute_pressure(n);
                    factor += p * (_tissue_height - h)
                              / this->_am.neighbors_of(n).size();
                }
                SpaceVec force = - pressure * factor * normal / 2.;
                
                edge->custom_links().a->state.add_force(force);
                edge->custom_links().b->state.add_force(force);
            }
        }
    }

    double compute_pressure(const std::shared_ptr<Cell>& cell) const final {
        return _elastic_modulus * (volume_of(cell) - 1.)
                / std::pow(_preferential_volume, 2);
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

    void update_parameters (const DataIO::Config& cfg) override {
        _elastic_modulus = get_as<double>("elastic_modulus", cfg, 
                                          _elastic_modulus);
        _tissue_height = get_as<double>("tissue_height", cfg, _tissue_height);
        _preferential_volume = get_as<double>("preferential_volume", cfg, 
                                              _preferential_volume);
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
            "pressure"
        });
    }

    std::vector<std::vector<double>> write_cell_energies () const final {
        std::vector<double> energies({});
        std::vector<double> pressures({});
        
        for (const auto& cell : this->_am.cells()) {
            energies.push_back(compute_energy(cell));
            pressures.push_back(compute_pressure(cell));
        }

        return std::vector<std::vector<double>>({
            energies,
            pressures
        });
    }
};

} // namespace WorkFunction2Dplus
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
