#ifndef UTOPIA_MODELS_COLLIER_WRITETASKS_HH
#define UTOPIA_MODELS_COLLIER_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"

using namespace Utopia::DataIO;

namespace Utopia::Models::Collier::DataIO{

/// Datamanager adaptor for fraction of cells of kind progenitor
auto density_notch = std::make_tuple(
    "Density_notch",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_protein_densities()[0]); },
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Notch"); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_progenitor

/// Datamanager adaptor for fraction of cells of kind progenitor
auto density_delta = std::make_tuple(
    "Density_delta",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_protein_densities()[1]); },
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Delta"); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_progenitor

/// Datamanager adaptor for fraction of cells of kind progenitor
auto density_nicd = std::make_tuple(
    "Density_nicd",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_protein_densities()[2]); },
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("NICD"); },
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_progenitor

/// Datamanager adaptor for spatial cell atoh1 levels
auto cell_notch = std::make_tuple(
    "Notch",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("CM"); },
    [](auto& dataset, auto& model) {
        auto cells = model.get_cm().cells();
        dataset->write(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.notch; }); },
    [](auto& group, auto& model) -> decltype(auto) {
        return group->open_dataset("Notch",
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
); // end cell_notch_adaptor

/// Datamanager adaptor for spatial cell atoh1 levels
auto cell_delta = std::make_tuple(
    "Delta",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("CM"); },
    [](auto& dataset, auto& model) {
        auto cells = model.get_cm().cells();
        dataset->write(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.delta; }); },
    [](auto& group, auto& model) -> decltype(auto) {
        return group->open_dataset("Delta",
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
); // end cell_delta_adaptor

/// Datamanager adaptor for spatial cell atoh1 levels
auto cell_nicd = std::make_tuple(
    "NICD",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("CM"); },
    [](auto& dataset, auto& model) {
        auto cells = model.get_cm().cells();
        dataset->write(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.nicd; }); },
    [](auto& group, auto& model) -> decltype(auto) {
        return group->open_dataset("NICD",
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
); // end cell_nicd_adaptor

} // namespace Utopia::Models::Collier::DataIO


#endif // UTOPIA_MODELS_COLLIER_WRITETASKS_HH