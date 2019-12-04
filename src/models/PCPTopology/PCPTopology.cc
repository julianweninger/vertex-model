#include <iostream>

#include "PCPTopology.hh"

using namespace Utopia::Models::PCPVertex;
using Utopia::get_as;


int main (int, char** argv) {
    try {
        // Initialize the PseudoParent from config file path
        Utopia::PseudoParent pp(argv[1]);
        auto model_cfg = pp.get_cfg()["PCPTopology"];

        // Initialize the main model instance and directly run it
        if (get_as<bool>("periodic_bc", model_cfg["PCPVertex"])) {
            PCPTopology<true>("PCPTopology", pp).run();
        }
        else
        {
            PCPTopology<false>("PCPTopology", pp).run();
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
