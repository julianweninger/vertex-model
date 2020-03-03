#include <iostream>

#include "PCPTopology.hh"

using namespace Utopia::Models::PCPVertex;
using namespace DataIO;
using Utopia::get_as;

/// Factory for model 
template<bool periodic_bc, bool polarity_proteins, typename ParentType>
auto model_factory(ParentType parent) {
    return PCPTopology<periodic_bc, polarity_proteins>("PCPTopology", parent,
        // statistics
        statistics_time_adaptor, cell_neighbourhood_adaptor,
        cell_area_adaptor<periodic_bc>,
        cell_area_histogram_adaptor<periodic_bc>,
        // the position adaptors
        vertex_position_adaptor, cell_position_adaptor<periodic_bc>,
        edge_link_adaptor);
}

int main (int, char** argv) {
    try {
        // Initialize the PseudoParent from config file path
        Utopia::PseudoParent pp(argv[1]);
        auto model_cfg = pp.get_cfg()["PCPTopology"];

        // Initialize the main model instance and directly run it
        if (get_as<bool>("periodic_bc", model_cfg["PCPVertex"]) and 
            get_as<double>("gamma", model_cfg["PCPVertex"]) > 0)
        {
            auto model = model_factory<true, true>(pp);
            model.run();
        }
        else if (get_as<bool>("periodic_bc", model_cfg["PCPVertex"]))
        {
            auto model = model_factory<true, false>(pp);
            model.run();
        }
        else if (get_as<double>("gamma", model_cfg["PCPVertex"]) > 0.)
        {
            auto model = model_factory<false, true>(pp);
            model.run();
        }
        else
        {
            auto model = model_factory<false, false>(pp);
            model.run();
        }    

        // Done
        return 0;
    }
    catch (Utopia::Exception& e) {
        return Utopia::handle_exception(e);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Exception occurred!" << std::endl;
        return 1;
    }
}
