#ifndef UTOPIA_MODELS_PCPTOPOLOGY_WRITETASKS_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"

using namespace Utopia::DataIO;

namespace Utopia::Models::PCPVertex::DataIO{

/// Datamanager adaptor for vertex-position
auto cell_neighbourhood_adaptor = std::make_tuple(

    // name of the task
    "Cell_neighbourhood",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_cell_neighbourhood_histogram());
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Cell_neighbourhood", {H5S_UNLIMITED, 9});
    },

    // attribute writer for basegroup
    [](auto& grp, auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");

        hdfdataset->add_attribute("dim_name__1", "bin");
        hdfdataset->add_attribute("coords_mode__bin", "start_and_step");
        hdfdataset->add_attribute("coords__bin", std::vector<std::size_t>{1, 1});
    }
    
); // end cell neighbourhood adaptor

/// Datamanager adaptor for vertex-position
template <bool periodic_bc>
auto cell_area_histogram_adaptor = std::make_tuple(

    // name of the task
    "Cell_area_histogram",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        auto cells = model.get_cells();
        std::array<int, 10> histogram = {0};
        std::array<double, 10> area = {0};
        histogram[0] = cells.size();
        for (auto c_weak : cells) {
            auto c = c_weak.lock();
            int neighbours = c->edges_ordered.size();
            neighbours = std::min(neighbours, 9);
            histogram[neighbours]++;

            auto [Lx, Ly] = model.get_domain_size();
            area[neighbours] += c->template area_abs<periodic_bc>(Lx, Ly);
        }
        for (int i = 1; i < 10; i++) {
            if (histogram[i] == 0) { continue; }
            area[0] += area[i];
            area[i] /= double(histogram[i]);
        }
        area[0] /= double(histogram[0]);
        dataset->write(area);
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Cell_area_histogram", {H5S_UNLIMITED, 10});
    },

    // attribute writer for basegroup
    [](auto& grp, auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");

        hdfdataset->add_attribute("dim_name__1", "num_neighbors");
        hdfdataset->add_attribute("coords_mode__num_neighbors",
                                  "start_and_step");
        hdfdataset->add_attribute("coords__num_neighbors",
                                  std::vector<std::size_t>{0, 1});
    }
    
); // end cell neighbourhood adaptor


} // namespace Utopia::Models::PCPVertex::DataIO


#endif // UTOPIA_MODELS_PCPTOPOLOGY_WRITETASKS_HH