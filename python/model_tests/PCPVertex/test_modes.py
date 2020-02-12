"""Tests of the output of the PCPVertex model"""

import numpy as np

import pytest

from utopya.testtools import ModelTest

# Configure the ModelTest class for PCPVertex
mtc = ModelTest("PCPVertex", test_file=__file__)

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