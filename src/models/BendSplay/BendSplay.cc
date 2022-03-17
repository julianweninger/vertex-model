#include <iostream>

#include "BendSplay.hh"
#include "BendSplay_write_tasks.hh"

using namespace Utopia::Models::BendSplay;
using namespace DataIO;


template<typename ParentType>
auto model_factory(ParentType &parent) {
    return BendSplay(
        "BendSplay", parent, {},
        std::make_tuple(
            energy_adaptor,
            time_energy_adaptor,
            orientations_adaptor,
            orientations_time_adaptor,
            coords_x_adaptor,
            coords_y_adaptor
        )
    );
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
