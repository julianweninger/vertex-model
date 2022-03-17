#ifndef UTOPIA_MODELS_BENDSPLAY_WRITETASKS_HH
#define UTOPIA_MODELS_BENDSPLAY_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"
#include <utopia/data_io/data_manager/defaults.hh>
#include <complex>
#include <cmath>

using namespace Utopia::DataIO;

namespace Utopia::Models::BendSplay::DataIO{

/// Datamanager adaptor for free energy
auto energy_adaptor = std::make_tuple(

    // name of the task
    "free_energy",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp;
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.free_energy());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("free_energy");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "free_energy_time");
    }
); // end energy_adaptor

/// Datamanager adaptor for timepoints
/** The dataset to which the other energy adaptors link their coordinate time
 */
auto time_energy_adaptor = std::make_tuple(

    // name of the task
    "free_energy_time",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp;
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_time());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("free_energy_time");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
    }
); // end time_energy_adaptor

/// Datamanager adaptor for orientation at every grid point
/** The dataset to which the other energy adaptors link their coordinate time
 */
auto orientations_adaptor = std::make_tuple(

    // name of the task
    "orientations",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp;
    },

    // writer function
    [](auto& dataset, auto& model) {
        const arma::mat& orientations = model.get_orientations();
        
        dataset->write(orientations.begin(), orientations.end(),
                       [](auto&& val) { return val; });
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        const arma::mat& orientations = m.get_orientations();
        return group->open_dataset("orientations",
                                   {H5S_UNLIMITED, orientations.n_elem});
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "orientations_time");
        hdfdataset->add_attribute("dim_name__1", "ids");
    }
); // end time_energy_adaptor

/// Datamanager adaptor for timepoints of grid data
auto orientations_time_adaptor = std::make_tuple(

    // name of the task
    "orientations_time",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp;
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_time());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("orientations_time", {H5S_UNLIMITED});
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
    }
); // end time_energy_adaptor

/// Datamanager adaptor for grid data x-coordinates
auto coords_x_adaptor = std::make_tuple(

    // name of the task
    "coords_x",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp;
    },

    // writer function
    [](auto& dataset, auto& model) {
        const arma::mat coords = model.get_coords_x();
        
        dataset->write(coords.begin(), coords.end(),
                       [](auto&& val) { return val; });
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        const arma::mat coords = m.get_coords_x();
        return group->open_dataset("coords_x", {H5S_UNLIMITED, coords.n_elem});
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "orientations_time");
        hdfdataset->add_attribute("dim_name__1", "ids");
    }
); // end coords_x_adaptor

/// Datamanager adaptor for grid data y-coordinates
auto coords_y_adaptor = std::make_tuple(

    // name of the task
    "coords_y",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp;
    },

    // writer function
    [](auto& dataset, auto& model) {
        const arma::mat coords = model.get_coords_y();
        
        dataset->write(coords.begin(), coords.end(),
                       [](auto&& val) { return val; });
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        const arma::mat coords = m.get_coords_y();
        return group->open_dataset("coords_y", {H5S_UNLIMITED, coords.n_elem});
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "orientations_time");
        hdfdataset->add_attribute("dim_name__1", "ids");
    }
); // end coords_y_adaptor


} // namespace Utopia::Models::BendSplay::DataIO


#endif // UTOPIA_MODELS_BENDSPLAY_WRITETASKS_HH
