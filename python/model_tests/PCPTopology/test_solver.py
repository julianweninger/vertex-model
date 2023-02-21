"""Tests of the output of the PCPTopology model"""

import numpy as np

import pytest

from utopya.testtools import ModelTest

# Configure the ModelTest class for PCPTopology
mtc = ModelTest("PCPTopology", test_file=__file__)

# Fixtures --------------------------------------------------------------------
# Define fixtures


# Tests -----------------------------------------------------------------------

def test_solver():
    # Create a Multiverse using the configuration for the original equations
    mv, dm = mtc.create_run_load(from_cfg="test_solver.yml", perform_sweep=True)

    for uni in dm['multiverse']:
        data = dm['multiverse'][uni]['data']['PCPTopology']

        assert 'Energy' in data
        assert 'Time' in data['Energy']
        assert 'Terms' in data['Energy']

        assert 'Energy' in data['PCPVertex']
        assert 'Time' in data['PCPVertex']['Energy']
        assert 'Terms' in data['PCPVertex']['Energy']

        assert 'Vertices' in data
        assert 'Edges' in data
        assert 'Cells' in data
