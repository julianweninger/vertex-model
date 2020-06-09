#include <iostream>

#include "PCPTopology.hh"

using namespace Utopia::Models::PCPVertex;
using namespace DataIO;
using Utopia::get_as;

/// Factory for model 
template<typename ParentType>
auto model_factory(ParentType parent) {
    return PCPTopology("PCPTopology", parent, {}, std::make_tuple(
        // statistics
        statistics_time_adaptor, cell_neighbourhood_adaptor,
        cell_area_adaptor<typename PCPTopology::CellType>,
        cell_area_histogram_adaptor,
        // the position adaptors
        vertices_adaptor,
        cells_adaptor<typename PCPVertex::Space::SpaceVec,
                      typename PCPTopology::CellType>,
        edges_adaptor));
}

int main (int, char** argv) {
    try {
        // Initialize the PseudoParent from config file path
        Utopia::PseudoParent pp(argv[1]);

        auto model = model_factory(pp);
        model.run();

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
