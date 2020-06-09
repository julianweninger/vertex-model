#include <iostream>

#include "NotchDelta.hh"
#include "NotchDelta_write_tasks.hh"

using namespace Utopia::Models::NotchDelta;
using namespace DataIO;


int main (int, char** argv) {
    try {
        // Initialize the PseudoParent from a config file path
        Utopia::PseudoParent pp(argv[1]);

        // Initialize the main model instance and directly run it
        NotchDelta("NotchDelta", pp, {}, std::make_tuple(
            density_time, density_progenitor, 
            density_hair, density_support, density_ratio_hair_support,
            number_hair_hair_contacts,
            CM_time, cell_type, cell_atoh1, cluster_id)
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
