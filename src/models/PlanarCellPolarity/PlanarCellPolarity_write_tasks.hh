#ifndef UTOPIA_MODELS_PLANARCELLPOLARITY_WRITETASKS_HH
#define UTOPIA_MODELS_PLANARCELLPOLARITY_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"
#include <complex>
#include <cmath>

using namespace Utopia::DataIO;

namespace Utopia::Models::PlanarCellPolarity::DataIO{

/** Available datatree:
 *      - Energy
 *          - Energy_time
 *          - Energy_total
 *          - Energy_cell_cell_polarity
 *          - Energy_polarity_exclusion
 *          - Energy_lagr_net_polaris
 *          - Energy_lagr_const_concentr
 *      - Polarity
 */

/// Datamanager adaptor for timepoints
auto time_energy_adaptor = std::make_tuple(
    // name of the task
    "Energy_time",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_time());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Time");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
    }
); // end time_energy_adaptor

/// Datamanager adaptor for total energy
auto energy_adaptor = std::make_tuple(

    // name of the task
    "Energy",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        std::vector<double> energy{};
        energy.reserve(5);

        energy.push_back(model.total_energy());
        energy.push_back(model.energy_cell_cell_polarity());
        energy.push_back(model.energy_polarity_exclusion());
        energy.push_back(model.energy_lagrange_net_polarisation());
        energy.push_back(model.energy_lagrange_const_concentration());

        dataset->write(energy);
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Energy", {H5S_UNLIMITED, 5});
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
        
        hdfdataset->add_attribute("dim_name__1", "energy");
        hdfdataset->add_attribute("coords__energy", 
            std::vector<std::string>({
                "total",
                "cell_cell_interaction",
                "polarity_exclusion",
                "lagrange_net_polarity",
                "lagrange_concentration"
            }));
    }
); // end energy_adaptor

template <typename SpaceVec, typename ProteinVec>
auto pcp_cells_adaptor = std::make_tuple(
    // name of the task
    "PCP_Cells",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("PlanarCellPolarity")->open_group("Cells");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& space = model.get_space();
        const auto& am = model.get_am();
        const auto& cells = am.cells();

        std::vector<SpaceVec> polarities;
        for (const auto& cell : cells) {
            SpaceVec center = am.barycenter_of(cell);
            SpaceVec polarity({0., 0.});
            for (auto [edge, flip] : cell->custom_links().edges) {
                ProteinVec sigma = std::get<0>(
                    model.get_polarity_proteins(edge, flip));
                std::size_t N = sigma.n_elem;

                auto a = edge->custom_links().a;
                SpaceVec _a = am.position_of(a);
                SpaceVec e_vec = am.displacement(edge);

                // iterate the proteins along edge
                double dN = 1. / static_cast<double>(N);
                for (std::size_t i = 0; i < N; i++) {
                    double s = sigma.at(i);

                    SpaceVec block_pos = _a + (i + 0.5) * dN * e_vec;
                    SpaceVec displ = space->displacement(center, block_pos);
                    
                    polarity += s * displ;
                }
            }
            polarities.push_back(polarity);
        }
        dataset->write(polarities.begin(), polarities.end(),
                       [](auto&& pol) { return pol[0]; });
        dataset->write(polarities.begin(), polarities.end(),
                       [](auto&& pol) { return pol[1]; });

        dataset->write(cells.begin(), cells.end(),
            [model](const auto& cell) {
                double cum_sigma = 0.;
                for (const auto& [edge, flip] : cell->custom_links().edges) {
                    cum_sigma += arma::accu(
                        std::get<0>(model.get_polarity_proteins(edge, flip)));
                }

                return cum_sigma;
            });
        dataset->write(cells.begin(), cells.end(),
            [model](const auto& cell) {
                double cum_sigma_2 = 0.;
                for (const auto& [edge, flip] : cell->custom_links().edges) {
                    ProteinVec sigma = std::get<0>(
                            model.get_polarity_proteins(edge, flip));
                    cum_sigma_2 += std::pow(arma::norm(sigma), 2);
                }

                return cum_sigma_2;
            });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {4, m.get_am().cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "property");
        hdfdataset->add_attribute("coords__property", 
                std::vector<std::string>({
                    "polarity_x",
                    "polarity_y",
                    "pcp_net_polarity",
                    "pcp_proteins"
                }));
        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
        const auto& cells = model.get_am().cells();
        std::vector<std::size_t> ids{};
        ids.reserve(cells.size());
        std::transform(cells.begin(), cells.end(), std::back_inserter(ids),
                       [](const auto& c) { return c->id(); });
        hdfdataset->add_attribute("coords__id", ids);
    }    
); // end cell position adaptor

/// Datamanager adaptor for edges properties
template <typename ProteinVec>
auto pcp_edges_adaptor = std::make_tuple(

    // name of the task
    "PCP_Proteins",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("PlanarCellPolarity")->open_group("Proteins");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& edges = am.edges();

        std::size_t N =  model.proteins_per_edge();

        std::vector<std::vector<double>> sigma;
        sigma.reserve(edges.size());
        for (const auto& edge : edges) {
            ProteinVec a, b;
            std::vector<double> _s;
            std::tie(a, b) = model.get_polarity_proteins(edge);
            for (std::size_t i = 0; i < N; i++) {
                _s.push_back(a[i]);
            }
            for (std::size_t i = 0; i < N; i++) {
                _s.push_back(b[i]);
            }
            sigma.push_back(_s);
        }

        for (std::size_t i = 0; i < N; i++) {
            dataset->write(edges.begin(), edges.end(),
                [model, i](const auto& edge) {
                    return std::get<0>(model.get_polarity_proteins(edge))(i);
                });
        }
        for (std::size_t i = 0; i < N; i++) {
            dataset->write(edges.begin(), edges.end(),
                [model, i](const auto& edge) {
                    return std::get<1>(model.get_polarity_proteins(edge))(i);
                });
        }
    },
                
    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2 * m.proteins_per_edge(), m.get_am().edges().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        // proteins
        hdfdataset->add_attribute("dim_name__0", "proteins");
        hdfdataset->add_attribute("coords_mode__proteins", "values");
        std::vector<std::string> proteins{};
        for (std::size_t i = 0; i < model.proteins_per_edge(); i++) {
            std::string name = "sigma_alpha__" + std::to_string(i);
            proteins.push_back(name);
        }
        for (std::size_t i = 0; i < model.proteins_per_edge(); i++) {
            std::string name = "sigma_beta__" + std::to_string(i);
            proteins.push_back(name);
        }
        hdfdataset->add_attribute("coords__proteins", proteins);
    
        // ids
        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
        const auto& edges = model.get_am().edges();
        std::vector<std::size_t> ids{};
        ids.reserve(edges.size());
        std::transform(edges.begin(), edges.end(), std::back_inserter(ids),
                       [](const auto& e) { return e->id(); });
        hdfdataset->add_attribute("coords__id", ids);
    }     
); // end edge link adaptor



/// Datamanager adaptor for cell properties
/** \details Properties are
 *      -# x component of polarity vector
 *      -# y component of polarity vector
 */
template <typename SpaceVec>
auto Polarity_adaptor = std::make_tuple(
    // name of the task
    "Polarities",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Polarities");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();

        std::vector<SpaceVec> polarities;
        for (const auto& c : cells) {
            SpaceVec polarity;
            polarity.zeros();
            for (auto [e, flip] : c->custom_links().edges) {
                auto a = e->custom_links().a;
                auto b = e->custom_links().b;
                if (flip) {
                    std::swap(a, b);
                }

                SpaceVec disp = am.displacement(a, b);
                double sigma = am.get_sigma(e, c);
                
                polarity += sigma * disp;
            }
            polarities.push_back(polarity);
        }
        dataset->write(polarities.begin(), polarities.end(),
                       [](auto&& pol) { return pol[0]; });
        dataset->write(polarities.begin(), polarities.end(),
                       [](auto&& pol) { return pol[1]; });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_am().cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "property");
        hdfdataset->add_attribute("coords__property", 
                std::vector<std::string>({
                    "polarity_x",
                    "polarity_y"
                }));
        
        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
        const auto& cells = model.get_am().cells();
        std::vector<std::size_t> ids{};
        ids.reserve(cells.size());
        std::transform(cells.begin(), cells.end(), std::back_inserter(ids),
                       [](const auto& c) { return c->id(); });
        hdfdataset->add_attribute("coords__id", ids);
        
        hdfdataset->add_attribute("AM_time", model.get_am_time());
    }    
); // end cell position adaptor

} // namespace Utopia::Models::PlanarCellPolarity::DataIO

#endif // UTOPIA_MODELS_PlanarCellPolarity_WRITETASKS_HH