#ifndef UTOPIA_MODELS_NOTCHDELTA_WRITETASKS_HH
#define UTOPIA_MODELS_NOTCHDELTA_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"

using namespace Utopia::DataIO;

namespace Utopia::Models::NotchDelta::DataIO{

/// Datamanager adaptor for timepoints
auto density_time = std::make_tuple(
    "Density_time",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_time()); },
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Time"); },
    [](auto& grp, auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time"); }
); // end time_adaptor

auto density_progenitor = std::make_tuple(
    "Density_progenitor",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_densities()[0]); },
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Progenitor"); },
    [](auto& grp, auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_progenitor

auto density_hair = std::make_tuple(
    "Density_hair",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_densities()[1]); },
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Hair"); },
    [](auto& grp, auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_hair

auto density_support = std::make_tuple(
    "Density_support",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_densities()[2]); },
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Support"); },
    [](auto& grp, auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_support

auto density_ratio_hair_support = std::make_tuple(
    "Density_ratio_hair_support",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_densities()[1]/model.get_densities()[2]); },
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Ratio_hair_support"); },
    [](auto& grp, auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_ratio

auto number_hair_hair_contacts = std::make_tuple(
    "Number_hair_hair_contacts",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Densities"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_hh_contacts()); },
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Number_hair_hair_contacts"); },
    [](auto& grp, auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time"); }
); // end density_ratio

/// Datamanager adaptor for timepoints
auto CM_time = std::make_tuple(
    "CM_time",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("CM"); },
    [](auto& dataset, auto& model) {
        dataset->write(model.get_time()); },
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Time"); },
    [](auto& grp, auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time"); }
); // end time_adaptor

auto cell_type = std::make_tuple(
    "Cell_type",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("CM"); },
    [](auto& dataset, auto& model) {
        auto cells = model.get_cm()->cells();
        dataset->write(cells.begin(), cells.end(),
            [](const auto& cell) {
                return int(cell->state.cell_type); }); },
    [](auto& group, auto& model) -> decltype(auto) {
        return group->open_dataset("Cell_type",
                                   {H5S_UNLIMITED, model.get_cm()->cells().size()}); },
    [](auto& grp, auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
        hdfdataset->add_attribute("dim_name__1", "cell_id");
        hdfdataset->add_attribute("content", "grid");
        hdfdataset->add_attribute("grid_shape",
                                  model.get_cm()->grid()->shape());
        hdfdataset->add_attribute("space_extent",
                                  model.get_cm()->grid()->space()->extent);
        hdfdataset->add_attribute("index_order", "F");}
); // end density_support

auto cell_atoh1 = std::make_tuple(
    "Atoh1",
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("CM"); },
    [](auto& dataset, auto& model) {
        auto cells = model.get_cm()->cells();
        dataset->write(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->custom_links().env->state.atoh1; }); },
    [](auto& group, auto& model) -> decltype(auto) {
        return group->open_dataset("Atoh1",
                                   {H5S_UNLIMITED, model.get_cm()->cells().size()}); },
    [](auto& grp, auto& m) {},
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
        hdfdataset->add_attribute("dim_name__1", "cell_id");
        hdfdataset->add_attribute("content", "grid");
        hdfdataset->add_attribute("grid_shape", 
                                  model.get_cm()->grid()->shape());
        hdfdataset->add_attribute("space_extent",
                                  model.get_cm()->grid()->space()->extent);
        hdfdataset->add_attribute("index_order", "F");}
); // end density_support


} // namespace Utopia::Models::NotchDelta::DataIO


#endif // UTOPIA_MODELS_NOTCHDELTA_WRITETASKS_HH