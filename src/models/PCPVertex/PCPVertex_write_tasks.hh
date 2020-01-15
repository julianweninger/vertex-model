#ifndef UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH
#define UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"

using namespace Utopia::DataIO;

namespace Utopia::Models::PCPVertex::DataIO{

/// Datamanager adaptor for timepoints
auto time_adaptor = std::make_tuple(

    // name of the task
    "Time",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp;
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_time());
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Time");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "write_time");
    }
); // end time_adaptor

/// Datamanager adaptor for vertex-position
auto vertex_position_adaptor = std::make_tuple(

    // name of the task
    "Vertex_position",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Vertex_position");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& vertices = model.get_vertices();
        dataset->write(vertices.begin(), vertices.end(),
                        [](auto&& vertex) {
                            return static_cast<double>(vertex.lock()->x);
                        });
        dataset->write(vertices.begin(), vertices.end(),
                        [](auto&& vertex) {
                            return static_cast<double>(vertex.lock()->y);
                        });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_vertices().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        // hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("dim_name__0", "idcs");
        // For ids, the dimensions are trivial
        hdfdataset->add_attribute("coords_mode__idcs", "start_and_step");
        hdfdataset->add_attribute("coords__idcs", std::vector<std::size_t>{1, 1});
    }
    
); // end vertex position adaptor

/// Datamanager adaptor for cell position
auto cell_position_adaptor = std::make_tuple(

    // name of the task
    "Cell_position",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Cell_position");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& cells = model.get_cells();
        dataset->write(cells.begin(), cells.end(),
                        [](auto&& cell) {
                            return static_cast<double>(cell.lock()->s->x);
                        });
        dataset->write(cells.begin(), cells.end(),
                        [](auto&& cell) {
                            return static_cast<double>(cell.lock()->s->y);
                        });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        // hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("dim_name__0", "idcs");
        // For ids, the dimensions are trivial
        hdfdataset->add_attribute("coords_mode__idcs", "start_and_step");
        hdfdataset->add_attribute("coords__idcs", std::vector<std::size_t>{1, 1});
    }    
); // end cell position adaptor

/// Datamanager adaptor for position
auto edge_link_adaptor = std::make_tuple(

    // name of the task
    "Edge_link",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Edge_link");
    },

    // writer function
    [](auto& dataset, auto& model) {
        auto vertices = model.get_vertices();
        int running_id = 0;
        for (auto &v : vertices) {
            v.lock()->current_id = running_id++;
        }

        const auto& edges = model.get_edges();
        dataset->write(edges.begin(), edges.end(),
                        [](auto&& edge) {
                            return static_cast<int>(edge.lock()->a->current_id);
                        });
        dataset->write(edges.begin(), edges.end(),
                        [](auto&& edge) {
                            return static_cast<int>(edge.lock()->b->current_id);
                        });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_edges().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        // hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("dim_name__0", "idcs");
        // For ids, the dimensions are trivial
        hdfdataset->add_attribute("coords_mode__idcs", "start_and_step");
        hdfdataset->add_attribute("coords__idcs", std::vector<std::size_t>{1, 1});
    }    
); // end edge link adaptor

} // namespace Utopia::Models::PCPVertex::DataIO


#endif // UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH