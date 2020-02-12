#include <iostream>

#include "PCPVertex.hh"
#include "PCPVertex_write_tasks.hh"

using namespace Utopia::Models::PCPVertex;
using namespace DataIO;
using Utopia::get_as;

/// Factory for model 
template<bool periodic_bc, bool polarity_proteins, typename ParentType>
auto model_factory(ParentType parent) {
    return PCPVertex<periodic_bc, polarity_proteins>("PCPVertex", parent,
        // the energy adaptors
        time_energy_adaptor, energy_adaptor, linetension_adaptor,
        areaelasticity_adaptor, contractility_adaptor,
        cell_cell_polarity_adaptor, polarity_exclusion_adaptor,
        lagrange_net_polarisation_adaptor, lagrange_const_concentration_adaptor,
        // the position adaptors
        vertex_position_adaptor, cell_position_adaptor<periodic_bc>,
        edge_link_adaptor);
}

int main (int, char** argv) {
    try {
        // Initialize the PseudoParent from config file path
        Utopia::PseudoParent pp(argv[1]);
        auto model_cfg = pp.get_cfg()["PCPVertex"];

        // Initialize the main model instance and directly run it
        if (get_as<bool>("periodic_bc", model_cfg) and 
            get_as<double>("gamma", model_cfg) > 0)
        {
            auto model = model_factory<true, true>(pp);
            model.run();
        }
        else if (get_as<bool>("periodic_bc", model_cfg))
        {
            auto model = model_factory<true, false>(pp);
            model.run();
        }
        else if (get_as<double>("gamma", model_cfg) > 0.)
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
