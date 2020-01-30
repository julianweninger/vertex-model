#ifndef UTOPIA_MODELS_PCPTOPOLOGY_WRITETASKS_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"

using namespace Utopia::DataIO;

namespace Utopia::Models::PCPVertex::DataIO{

/// Datamanager adaptor for timepoints
auto time_histogram_adaptor = std::make_tuple(

    // name of the task
    "Time_histogram",

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
        return group->open_dataset("Time_histogram");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) { },

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
    }
); // end time_adaptor

/// Datamanager adaptor for vertex-position
auto cell_neighbourhood_adaptor = std::make_tuple(

    // name of the task
    "Cell_neighbourhood",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp;
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
        hdfdataset->add_attribute("coords__time", "Time_histogram");

        hdfdataset->add_attribute("dim_name__1", "bin");
        hdfdataset->add_attribute("coords_mode__bin", "start_and_step");
        hdfdataset->add_attribute("coords__bin", std::vector<std::size_t>{1, 1});
    }
    
); // end cell neighbourhood adaptor


} // namespace Utopia::Models::PCPVertex::DataIO


#endif // UTOPIA_MODELS_PCPTOPOLOGY_WRITETASKS_HH