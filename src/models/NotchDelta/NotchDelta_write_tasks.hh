#ifndef UTOPIA_MODELS_NOTCHDELTA_WRITETASKS_HH
#define UTOPIA_MODELS_NOTCHDELTA_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"

using namespace Utopia::DataIO;

namespace Utopia::Models::NotchDelta::DataIO{

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

} // namespace Utopia::Models::NotchDelta::DataIO


#endif // UTOPIA_MODELS_NOTCHDELTA_WRITETASKS_HH