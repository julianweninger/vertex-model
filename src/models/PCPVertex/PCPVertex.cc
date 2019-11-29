#include <iostream>

#include "PCPVertex.hh"

using namespace Utopia::Models::PCPVertex;
using Utopia::get_as;


int main (int, char** argv) {
    try {
        // Initialize the PseudoParent from config file path
        Utopia::PseudoParent pp(argv[1]);
        auto model_cfg = pp.get_cfg()["PCPVertex"];

        // Initialize the main model instance and directly run it
        if (get_as<bool>("periodic_bc", model_cfg)) {
            PCPVertex<true>("PCPVertex", pp).run();
        }
        else
        {
            PCPVertex<false>("PCPVertex", pp).run();
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
