"""Tests of the output of the PCPVertex model"""

import numpy as np

import pytest

from utopya.testtools import ModelTest

# Configure the ModelTest class for PCPVertex
mtc = ModelTest("PCPVertex", test_file=__file__)

# Fixtures --------------------------------------------------------------------
# Define fixtures


# Tests -----------------------------------------------------------------------
def test_output():
    # Create a Multiverse using the configuration for the original equations
    mv, dm = mtc.create_run_load(from_cfg="test.yml")

    data = dm['multiverse'][0]['data']['PCPVertex']

    assert 'Energy' in data
    assert 'Time' in data['Energy']
    assert 'Total' in data['Energy']
    assert 'Linetension' in data['Energy']
    assert 'Areaelasticity' in data['Energy']
    assert 'Cell_contractility' in data['Energy']
    assert 'Edge_contractility' in data['Energy']

    # assert that initial state was written
    assert len(data['Energy']['Time'].data) == 4

    assert 'Vertices' in data
    assert 'Edges' in data
    assert 'Cells' in data

    assert len(data['Vertices']) == 4
    assert len(data['Edges']) == 4
    assert len(data['Cells']) == 4