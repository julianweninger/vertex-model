#define BOOST_TEST_MODULE MovCotest
#include <boost/test/unit_test.hpp>
#include <utopia/core/types.hh>
#include <utopia/core/agent.hh>
#include "../MovCo.hh"
#include "../MovCo_write_tasks.hh"

using namespace Utopia::Models::MovCo;
using namespace DataIO;
using namespace Utopia;

    // Factory for model 
    template<Mode runmode>
    MovCo<runmode> model_fac() {
        // initialize pp inside to avoid problems with datamanager and memory 
        // access
        PseudoParent pp("./interaction_test_cfg.yml");
        return  MovCo<runmode>("MovCoTEST", pp, payoff_adaptor, 
                            strategy_adaptor, groups_adaptor, resource_adaptor, 
                            position_x_adaptor, position_y_adaptor, 
                            displace_x_adaptor, displace_y_adaptor, 
                            strat_avg_adaptor, time_adaptor, full_groups_adaptor, 
                            full_c_groups_adaptor);
}  


struct MovCoFixture {
    ~MovCoFixture() {model.get_logger()->info("Tear it down");
                     spdlog::drop_all();
                     std::remove("complex_data.h5");}
    // for unit tests the runmode is not important (just changes perform_step)
    MovCo<Mode::basic> model = model_fac<Mode::basic>();
};


BOOST_FIXTURE_TEST_SUITE(all_units, MovCoFixture)

    BOOST_AUTO_TEST_CASE(model_) {
        model.run();
    }

    BOOST_AUTO_TEST_CASE(clean_up)
    {
        // Get the first agent
        auto& agt = model.get_am().agents()[0];
        // Get the cell
        auto& cell = model.get_cm().cells()[0];
        // modify the agent state
        agt->state().payoff = 10.;
        agt->state().received = 42.;
        auto new_agt_state = model._clean_up_agents(agt);
        // check if the state is cleaned up
        BOOST_TEST(new_agt_state.payoff == 0.);
        BOOST_TEST(new_agt_state.received == 0.);
        //modify the cell state
        cell->state().n_agents = 10;
        cell->state().resources = 42.;
        cell->state().n_coop = 15;
        auto new_cell_state = model._clean_up_cells(cell);
        // check wether the cell state is clean now
        BOOST_TEST(new_cell_state.n_agents == 0);
        BOOST_TEST(new_cell_state.resources == 0.);
        BOOST_TEST(new_cell_state.n_coop == 0);
    }

    BOOST_AUTO_TEST_CASE(place_resources_coop)
    {   
        // Get the first agent
        auto agt = model.get_am().agents()[0];
        // Get the cell
        auto cell = model.get_cm().cells()[0];
        // and the agent state 
        auto old_agt_state = agt->state(); 
        // make him a cooperator
        agt->state().strategy = Strategy::coop;
        // let him put a resource unit on the cell
        auto state = model._place_resources(agt);
        BOOST_TEST(state.resources == old_agt_state.resources - 1);
        BOOST_TEST(state.payoff == - 1);
        
        // Get the cell state (only one cell for testing)
        auto cell_state = model.get_cm().cells()[0]->state();
        // check wether the state is correct
        BOOST_TEST(cell_state.n_agents == 1);
        BOOST_TEST(cell_state.resources == 1.9);
    }


    BOOST_AUTO_TEST_CASE(place_resources_defect)
    {
        // Get the first agent
        auto agt = model.get_am().agents()[0];
        //and his state
        auto old_agt_state = agt->state(); 
        // make him a defector
        agt->state().strategy = Strategy::defect;
        // let him put a resource unit on the cell
        auto state = model._place_resources(agt);
        BOOST_TEST(state.resources == old_agt_state.resources);
        BOOST_TEST(state.payoff == 0.);
        
        // Get the cell state (only one cell for testing)
        auto cell_state = model.get_cm().cells()[0]->state();
        // check wether the state is correct
        BOOST_TEST(cell_state.n_agents == 1);
        BOOST_TEST(cell_state.resources == 0.);

    }


    BOOST_AUTO_TEST_CASE(_pay_out_coop_2agt)
    {   
        // Get the first agent
        auto agt = model.get_am().agents()[0];
        // Get the cell
        auto cell = model.get_cm().cells()[0];
        //and his state
        auto old_agt_state = agt->state(); 
        // make him a cooperator
        agt->state().strategy = Strategy::coop;
        // manipulate the cell state 
        cell->state().n_agents = 2;
        cell->state().resources = 2;
        // split the resources among agents (in this case just the one)
        auto state = model._pay_out(agt);
        BOOST_TEST(state.resources == old_agt_state.resources + 1);
        BOOST_TEST(state.payoff == 1.);
        BOOST_TEST(state.received == 1.);
    }


    BOOST_AUTO_TEST_CASE(_pay_out_coop_1agt)
    {
        // Get the first agent
        auto agt = model.get_am().agents()[0];
        // Get the cell
        auto cell = model.get_cm().cells()[0];
        //and his state
        auto old_agt_state = agt->state(); 
        // make him a cooperator
        agt->state().strategy = Strategy::coop;
        // manipulate the cell state 
        cell->state().n_agents = 1;
        // no interactions, no payout
        auto state = model._pay_out(agt);
        BOOST_TEST(state.resources == old_agt_state.resources);
        BOOST_TEST(state.payoff == 0.);
        BOOST_TEST(state.received == 0.);
    }

    BOOST_AUTO_TEST_CASE(_pay_out_defect_1agt)
    {
        // Get the first agent
        auto agt = model.get_am().agents()[0];
        // Get the cell
        auto cell = model.get_cm().cells()[0];
        //and his state
        auto old_agt_state = agt->state(); 
        // make him a cooperator
        agt->state().strategy = Strategy::defect;
        // manipulate the cell state 
        cell->state().n_agents = 1;
        // no interactions, no payout
        auto state = model._pay_out(agt);
        BOOST_TEST(state.resources == old_agt_state.resources);
        BOOST_TEST(state.payoff == 0.);
        BOOST_TEST(state.received == 0.);
    }


    // Test the reproduction
    BOOST_AUTO_TEST_CASE(noreproduction)
    {
        // Get the first agent
        auto agt = model.get_am().agents()[0];
        // Get the cell
        auto cell = model.get_cm().cells()[0];
        // make the agent a cooperator
        agt->state().strategy = Strategy::defect;
        // give him 800 resource units
        agt->state().resources = 800;
        // Get the state
        auto old_agt_state = agt->state();
        // apply repro function
        BOOST_TEST(model.get_am().agents().size() == 1);
        auto state = model._reproduction(agt);
        BOOST_TEST(state.resources == old_agt_state.resources);
        BOOST_TEST(state.payoff == old_agt_state.payoff);
        BOOST_TEST(state.received == old_agt_state.received);
        // make sure there is only one agent
        BOOST_TEST(model.get_am().agents().size() == 1);
    }

    BOOST_AUTO_TEST_CASE(reproduction)
    {
        BOOST_TEST(model.get_cm().cells().size() == 1);
        // Get the first agent
        auto agt = model.get_am().agents()[0];
        // Get the cell
        auto cell = model.get_cm().cells()[0];
        // make the agent a defector
        agt->state().strategy = Strategy::defect;
        // give him 1100 resource units
        agt->state().resources = 1100;
        // Get the state
        auto old_agt_state = agt->state();
        // apply repro function
        auto state = model._reproduction(agt);
        BOOST_TEST(state.resources == old_agt_state.resources / 2);
        BOOST_TEST(state.payoff == old_agt_state.payoff);
        BOOST_TEST(state.received == old_agt_state.received);
        // make sure there are 2 agents now
        BOOST_TEST(model.get_am().agents().size() == 2);
        // get the new agent and test his state
        auto new_agt = model.get_am().agents()[1];
        BOOST_TEST(new_agt->state().resources == old_agt_state.resources / 2);
        BOOST_TEST(new_agt->state().payoff == 0.);
        BOOST_TEST(new_agt->state().received == 0.);
        // and his position
        BOOST_TEST(new_agt->position()[0] == agt->position()[0]);
        BOOST_TEST(new_agt->position()[1] == agt->position()[1]);
    }

    BOOST_AUTO_TEST_CASE(costoflivingandremove)
    {
        // Get the first agent
        auto agt = model.get_am().agents()[0];
        // give him 1100 resource units
        agt->state().resources = 1100;
        // Create another agent
        auto new_agt = model.get_am().add_agent();
        // give him the cost of living as resource units
        new_agt->state().resources = get_as<double>("cost_of_living", model.get_cfg());
        // apply the _cost rule
        apply_rule<false>(model._cost, model.get_am().agents());
        // now check if the resources were set correctly
        BOOST_TEST(agt->state().resources == 1100 - 0.1);
        BOOST_TEST(new_agt->state().resources == 0.);
        // let the newagt be removed
        model._remove_no_res();
        // check if there is only one agent
        BOOST_TEST(model.get_am().agents().size() == 1);
        // Try again for below cost of living resources
        // Create another agent
        new_agt = model.get_am().add_agent();
        // give him the cost of living as resource units
        new_agt->state().resources = get_as<double>("cost_of_living", 
                                                    model.get_cfg()) - 1.;
        // apply the _cost rule
        apply_rule<false>(model._cost, model.get_am().agents());
        // let the newagt be removed
        model._remove_no_res();
        // check if there is only one agent
        BOOST_TEST(model.get_am().agents().size() == 1);
    }

    BOOST_AUTO_TEST_CASE(removeexcessagents)
    {
        // Create another agent
        auto new_agt = model.get_am().add_agent();
        // remove excess agts
        model._remove_exc_agents();
        // check if there is only one agent
        BOOST_TEST(model.get_am().agents().size() == 1);
    }

    BOOST_AUTO_TEST_CASE(removenores)
    {
        // Get the first agent
        auto agt = model.get_am().agents()[0];
        // Create another agent
        auto new_agt = model.get_am().add_agent();
        // Give the first 0 resources and the second below 0 resources
        agt->state().resources = 0.;
        new_agt->state().resources = -42;
        // Remove agents with 0 or below resources
        model._remove_no_res();
        // Check that both agents were removed
        BOOST_TEST(model.get_am().agents().size() == 0);
    }

    BOOST_AUTO_TEST_CASE(get_avg_strat_nan)
    {
        // check what happens with 0 agents
        // Get the agent
        auto agt = model.get_am().agents()[0];
        // remove the agent
        model.get_am().remove_agent(agt);
        // Check if avg strat is nan
        BOOST_TEST(std::isnan(model.get_avg_strat()));
    }

    BOOST_AUTO_TEST_CASE(get_avg_strat)
    {   
        // Get the agent
        auto agt_0 = model.get_am().agents()[0];
        // Create another agent
        auto agt_1 = model.get_am().add_agent();
        // ... and another
        auto agt_2 = model.get_am().add_agent();
        // Set one agent's strategy to coop the other's to defect
        agt_0->state().strategy = Strategy::coop;
        agt_1->state().strategy = Strategy::defect;
        agt_2->state().strategy = Strategy::defect;
        // Check if avg strat is one third
        BOOST_TEST(model.get_avg_strat() == 1. / 3.);
    }

    BOOST_AUTO_TEST_CASE(retrieve)
    {
        // Get the cell 
        auto& cell = model.get_cm().cells()[0];
        // modify the cell state
        cell->state().n_agents = 10;
        cell->state().n_coop = 43;
        // run _retrieve group_data
        model._retrieve_group_data();
        // Check wether the vectors contain the correct data
        BOOST_TEST(model.get_num_agts_cell() == std::vector({10}));
        BOOST_TEST(model.get_num_coop_cell() == std::vector({43}));
        BOOST_TEST(model.get_cell_ids() == std::vector<std::size_t>({0}));
        // clean up cells
        cell->state() = model._clean_up_cells(cell);
        // run _retrieve group_data
        model._retrieve_group_data();
        // Check wether the vectors contain no data
        BOOST_TEST(model.get_num_agts_cell().size() == 0);
        BOOST_TEST(model.get_num_coop_cell().size() == 0);
        BOOST_TEST(model.get_cell_ids().size() == 0);
    }

BOOST_AUTO_TEST_SUITE_END()