"""Tests of the output of the PCPTopology model"""

import numpy as np

import pytest

from utopya.testtools import ModelTest

# Configure the ModelTest class for PCPTopology
mtc = ModelTest("PCPTopology", test_file=__file__)

# Fixtures --------------------------------------------------------------------
# Define fixtures


# Tests -----------------------------------------------------------------------

def test_modes():
    """Test the modes of the model
    
    The model has the following modes:
        - periodic_bc: true / false
        - polarity_proteins: true (active) / false (inactive)
    """
    # Create a Multiverse using the default model configuration
    mv = mtc.create_mv()

    # Run a single simulation
    mv.run_single()

    # Load data using the DataManager and the default load configuration
    mv.dm.load_from_cfg(print_tree=True)

    # Assert that data was loaded, i.e. that data was written
    assert len(mv.dm)


    # Create a Multiverse using the configuration for the original equations
    mv = mtc.create_mv(from_cfg="test_mode_periodic_bc.yml")

    # Run a single simulation
    mv.run_single()

    # Load data using the DataManager and the default load configuration
    mv.dm.load_from_cfg(print_tree=True)

    # Assert that data was loaded, i.e. that data was written
    assert len(mv.dm)  


    # Create a Multiverse using the configuration for the original equations
    mv = mtc.create_mv(from_cfg="test_mode_non_periodic_bc.yml")

    # Run a single simulation
    mv.run_single()

    # Load data using the DataManager and the default load configuration
    mv.dm.load_from_cfg(print_tree=True)

    # Assert that data was loaded, i.e. that data was written
    assert len(mv.dm)  


    # Create a Multiverse using the configuration for the original equations
    mv = mtc.create_mv(from_cfg="test_mode_no_polarity_proteins.yml")

    # Run a single simulation
    mv.run_single()

    # Load data using the DataManager and the default load configuration
    mv.dm.load_from_cfg(print_tree=True)

    # Assert that data was loaded, i.e. that data was written
    assert len(mv.dm)    


    # Create a Multiverse using the default configuration for state_5 dynamics
    mv = mtc.create_mv(from_cfg="test_mode_polarity_proteins.yml")

    # Run a single simulation
    mv.run_single()

    # Load data using the DataManager and the default load configuration
    mv.dm.load_from_cfg(print_tree=True)

    mv.pm.plot_from_cfg(plot_only=["energy_polarity"])

    # Assert that data was loaded, i.e. that data was written
    assert len(mv.dm)

def test_output():
    # Create a Multiverse using the configuration for the original equations
    mv, dm = mtc.create_run_load(from_cfg="test_mode_no_polarity_proteins.yml")

    data = dm['multiverse'][0]['data']['PCPTopology']

    assert 'Energy' in data['PCPVertex']
    assert 'Time' in data['PCPVertex']['Energy']
    assert 'Total' in data['PCPVertex']['Energy']
    assert 'Linetension' in data['PCPVertex']['Energy']
    assert 'Areaelasticity' in data['PCPVertex']['Energy']
    assert 'Contractility' in data['PCPVertex']['Energy']

    # assert that initial state was written
    assert len(data['PCPVertex']['Energy']['Time'].data) == 1801
    # NOTE number of steps iterated is:
    #          1 initial state
    #       600 during prolog (equilibration + 1 cell division)
    #       400 per step (3 steps run)
    #       Total 1801

    assert 'Vertex_position' in data['PCPVertex']
    assert 'Edge_link' in data['PCPVertex']
    assert 'Cell_position' in data['PCPVertex']

    assert len(data['PCPVertex']['Vertex_position']) == 10
    assert len(data['PCPVertex']['Edge_link']) == 10
    assert len(data['PCPVertex']['Cell_position']) == 10
    # NOTE number of writes is:
    #       1 initial state
    #       3 during prolog
    #       2 per step (3 steps run)
    #       Total 10

    assert 'Vertex_position' in data
    assert 'Edge_link' in data
    assert 'Cell_position' in data

    
    # Create a Multiverse using the configuration for the original equations
    mv, dm = mtc.create_run_load(from_cfg="test_mode_polarity_proteins.yml")

    data = dm['multiverse'][0]['data']['PCPTopology']

    assert 'Energy' in data['PCPVertex']
    assert 'Time' in data['PCPVertex']['Energy']
    assert 'Total' in data['PCPVertex']['Energy']
    assert 'Linetension' in data['PCPVertex']['Energy']
    assert 'Areaelasticity' in data['PCPVertex']['Energy']
    assert 'Contractility' in data['PCPVertex']['Energy']
    assert 'Cell_cell_polarity' in data['PCPVertex']['Energy']
    assert 'Polarity_exclusion' in data['PCPVertex']['Energy']
    assert 'Lagrange_net_polarisation' in data['PCPVertex']['Energy']
    assert 'Lagrange_const_concentration' in data['PCPVertex']['Energy']

    # assert that initial state was written
    assert len(data['PCPVertex']['Energy']['Time'].data) == 1801
    # NOTE number of steps iterated is:
    #          1 initial state
    #       600 during prolog (equilibration + 1 cell division)
    #       400 per step (3 steps run)
    #       Total 1801

    assert 'Vertex_position' in data['PCPVertex']
    assert 'Edge_link' in data['PCPVertex']
    assert 'Cell_position' in data['PCPVertex']

    assert len(data['PCPVertex']['Vertex_position']) == 10
    assert len(data['PCPVertex']['Edge_link']) == 10
    assert len(data['PCPVertex']['Cell_position']) == 10
    # NOTE number of writes is:
    #       1 initial state
    #       3 during prolog
    #       2 per step (3 steps run)
    #       Total 10

    assert 'Vertex_position' in data
    assert 'Edge_link' in data
    assert 'Cell_position' in data
