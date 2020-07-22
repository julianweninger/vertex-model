#include <iostream>

#include "Collier.hh"
#include "Collier_write_tasks.hh"
#include "../NotchDelta/Differentiation_write_tasks.hh"

using namespace Utopia::Models::Collier;
using namespace Utopia::Models::Differentiation::DataIO;


int main (int, char** argv) {
    try {
        // Initialize the PseudoParent from a config file path
        Utopia::PseudoParent pp(argv[1]);

        // Initialize the main model instance and directly run it
        Collier("Collier", pp, {},
            std::make_tuple(density_time,  density_progenitor, density_hair,
                density_support, density_rosettes, DataIO::density_notch,
                DataIO::density_delta, DataIO::density_nicd, 
                CM_time, cell_type, cluster_id, DataIO::cell_notch,
                DataIO::cell_delta, DataIO::cell_nicd)
        ).run();

        // Done.
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
