#ifndef UTOPIA_MODELS_DIFFERENTIATION_WRITETASKS_HH
#define UTOPIA_MODELS_DIFFERENTIATION_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"

using namespace Utopia::DataIO;

namespace Utopia::Models::Differentiation::DataIO{

/// Datamanager adaptor for timepoints
auto density_time = std::make_tuple(
    "Density_time",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_time()); },
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Time"); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time"); }
); // end time_adaptor

/// Datamanager adaptor for fraction of cells of kind progenitor
auto density_progenitor = std::make_tuple(
    "Density_progenitor",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_densities()[0]); },
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Progenitor"); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_progenitor

/// Datamanager adaptor for fraction of cells of kind hair
auto density_hair = std::make_tuple(
    "Density_hair",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_densities()[1]); },
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Hair"); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_hair

/// Datamanager adaptor for fraction of cells of kind support
auto density_support = std::make_tuple(
    "Density_support",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_densities()[2]); },
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Support"); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_support

/// Datamanager adaptor for fraction of cells that form a rosette
auto density_rosettes = std::make_tuple(
    "Density_rosettes",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_num_rosettes() /
                       static_cast<double>(model.get_cm().cells().size()));
    },
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Rosettes"); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_rosettes

/// Datamanager adaptor for timepoints of the cell manager's data
auto CM_time = std::make_tuple(
    "CM_time",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("CM"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_time()); },
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Time"); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time"); }
); // end CM_time_adaptor

/// Datamanager adaptor for spatial cell types
auto cell_type = std::make_tuple(
    "Cell_type",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("CM"); },
    [](auto& dataset, auto& model) {
        auto cells = model.get_cm().cells();
        dataset->write(cells.begin(), cells.end(),
            [](const auto& cell) {
                return int(cell->state.cell_type); }); },
    [](auto& group, auto& model) -> decltype(auto) {
        return group->open_dataset("Cell_type",
                    {H5S_UNLIMITED, model.get_cm().cells().size()}); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
        hdfdataset->add_attribute("dim_name__1", "cell_id");
        hdfdataset->add_attribute("content", "grid");
        hdfdataset->add_attribute("grid_shape",
                                  model.get_cm().grid()->shape());
        hdfdataset->add_attribute("space_extent",
                                  model.get_cm().grid()->space()->extent);
        hdfdataset->add_attribute("index_order", "F");}
); // end cell_type_adaptor

/// Datamanager adaptor for spatial cell atoh1 levels
auto cell_atoh1 = std::make_tuple(
    "Atoh1",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("CM"); },
    [](auto& dataset, auto& model) {
        auto cells = model.get_cm().cells();
        dataset->write(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.atoh1; }); },
    [](auto& group, auto& model) -> decltype(auto) {
        return group->open_dataset("Atoh1",
                    {H5S_UNLIMITED, model.get_cm().cells().size()}); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
        hdfdataset->add_attribute("dim_name__1", "cell_id");
        hdfdataset->add_attribute("content", "grid");
        hdfdataset->add_attribute("grid_shape", 
                                  model.get_cm().grid()->shape());
        hdfdataset->add_attribute("space_extent",
                                  model.get_cm().grid()->space()->extent);
        hdfdataset->add_attribute("index_order", "F");}
); // end cell_atoh1_adaptor

/// Datamanager adaptor for spatial cell cluster ids
auto cluster_id = std::make_tuple(
    "Cluster_id",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("CM"); },
    [](auto& dataset, auto& model) {
        model.identify_clusters();

        auto cells = model.get_cm().cells();
        dataset->write(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.cluster_id; }); },
    [](auto& group, auto& model) -> decltype(auto) {
        return group->open_dataset("Cluster_id",
                    {H5S_UNLIMITED, model.get_cm().cells().size()}); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
        hdfdataset->add_attribute("dim_name__1", "cell_id");
        hdfdataset->add_attribute("content", "grid");
        hdfdataset->add_attribute("grid_shape", 
                                  model.get_cm().grid()->shape());
        hdfdataset->add_attribute("space_extent",
                                  model.get_cm().grid()->space()->extent);
        hdfdataset->add_attribute("index_order", "F");}
); // end cluster_id_adaptor

} // namespace Utopia::Models::Differentiation::DataIO


#endif // UTOPIA_MODELS_DIFFERENTIATION_WRITETASKS_HH