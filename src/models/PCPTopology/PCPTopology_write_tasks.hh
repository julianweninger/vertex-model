#ifndef UTOPIA_MODELS_PCPTOPOLOGY_WRITETASKS_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"

using namespace Utopia::DataIO;

namespace Utopia::Models::PCPVertex::DataIO{

/** Available datatree:
 *      - Energy
 *          - Continuous_time
 *      - Statistics
 *          - Cell_neighbourhood
 *          - Cell_area_histogram
 */

/// Datamanager adaptor for timepoint mapping to PCPVertex model
auto continuous_time_adaptor = std::make_tuple(

    // name of the task
    "Continuous_time",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_continuous_time());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Continuous_time");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end time_adaptor

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
        const auto& cells = model.get_am().cells();

        std::array<int, 7> histogram = {0};
        for (auto c : cells) {
            int neighbours = c->custom_links().edges.size();
            neighbours = std::min(neighbours, 9);
            if (neighbours < 3) {
                throw std::runtime_error("Encountered cell with less than 3 "
                    "neighbours during writing of cell_neighbourhood_adaptor!");
            }
            histogram[neighbours - 3]++;
        }

        dataset->write(histogram);
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Cell_neighbourhood", {H5S_UNLIMITED, 7});
    },

    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");

        hdfdataset->add_attribute("dim_name__1", "num_neighbors");
        hdfdataset->add_attribute("coords_mode__num_neighbors",
                                  "start_and_step");
        hdfdataset->add_attribute("coords__num_neighbors",
                                  std::vector<std::size_t>{3, 1});
    }
    
); // end cell neighbourhood adaptor

/// Datamanager adaptor for vertex-position
auto cell_area_histogram_adaptor = std::make_tuple(

    // name of the task
    "Cell_area_histogram",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();
        std::array<int, 8> histogram = {0};
        std::array<double, 8> area = {0};
        histogram[0] = cells.size();
        for (const auto& c : cells) {
            int neighbours = c->custom_links().edges.size();
            neighbours = std::min(neighbours, 9);
            histogram[neighbours-2]++;

            area[neighbours-2] += am.area_of(c);
        }
        for (int i = 1; i < 8; i++) {
            if (histogram[i] == 0) { continue; }
            area[0] += area[i];
            area[i] /= double(histogram[i]);
        }
        area[0] /= double(histogram[0]);
        dataset->write(area);
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Cell_area_histogram", {H5S_UNLIMITED, 8});
    },

    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");

        hdfdataset->add_attribute("dim_name__1", "num_neighbors");
        hdfdataset->add_attribute("coords_mode__num_neighbors",
                                  "start_and_step");
        hdfdataset->add_attribute("coords__num_neighbors",
                                  std::vector<std::size_t>{2, 1});
    }
    
); // end cell neighbourhood adaptor


} // namespace Utopia::Models::PCPVertex::DataIO


#endif // UTOPIA_MODELS_PCPTOPOLOGY_WRITETASKS_HH