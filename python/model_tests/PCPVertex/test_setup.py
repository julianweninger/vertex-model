"""Tests of the output of the PCPVertex model"""

import logging

from utopya.testtools import ModelTest

# Local constants
log = logging.getLogger(__name__)

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


def test_HC_regular_lattice():
    def run_sim(*, HC_structure: str, num_rows: int=20, num_cols: int=10):
        mv, dm = mtc.create_run_load(parameter_space={
            'PCPVertex': {'agent_manager': {'setup_params': { 'hexagonal': {
                'HC_structure': HC_structure,
                'lattice_rows': num_rows,
                'lattice_columns': num_cols}}
            }}})

        return dm['multiverse'][0]['data']['PCPVertex']['Cells'][0]

    def get_cell_data(cells, *, property: str):
        types = cells.sel(property='cell_type')

        prop_data = cells.sel(property=property)
        HC = prop_data.where(types == 1)
        SC = prop_data.where(types == 2)
        prog = prop_data.where(types == 0)

        assert prog.count().data == 0

        return [prop_data, HC, SC]


    log.info("Testing '{}' HC lattice ...", 'ratio_1_to_2')
        
    cells = run_sim(HC_structure='ratio_1_to_2', num_cols=12)
    [data, HC, SC] = get_cell_data(cells, property='num_hair_neighbors')

    assert 2 * HC.count().data == SC.count().data
    assert (HC.dropna(dim='id') == 0).all()
    assert (SC.dropna(dim='id') == 3).all()


    log.info("Testing '{}' HC lattice ...", 'ratio_1_to_3')

    cells = run_sim(HC_structure='ratio_1_to_3')
    [data, HC, SC] = get_cell_data(cells, property='num_hair_neighbors')

    assert 3 * HC.count().data == SC.count().data
    assert (HC.dropna(dim='id') == 0).all()
    assert (SC.dropna(dim='id') == 2).all()


    log.info("Testing '{}' HC lattice ...", 'ratio_1_to_4')

    cells = run_sim(HC_structure='ratio_1_to_4', num_cols=10)
    [data, HC, SC] = get_cell_data(cells, property='num_hair_neighbors')

    assert 4 * HC.count().data == SC.count().data
    assert (HC.dropna(dim='id') == 0).all()
    assert (SC.dropna(dim='id') >= 1).all()
    assert (SC.dropna(dim='id') <= 2).all()


    log.info("Testing '{}' HC lattice ...", 'ratio_1_to_5')

    cells = run_sim(HC_structure='ratio_1_to_5', num_cols=12, num_rows=16)
    [data, HC, SC] = get_cell_data(cells, property='num_hair_neighbors')

    assert 5 * HC.count().data == SC.count().data
    assert (HC.dropna(dim='id') == 0).all()
    assert (SC.dropna(dim='id') >= 1).all()
    assert (SC.dropna(dim='id') <= 2).all()


    log.info("Testing '{}' HC lattice ...", 'ratio_1_to_6')

    cells = run_sim(HC_structure='ratio_1_to_6', num_cols=14, num_rows=28)
    [data, HC, SC] = get_cell_data(cells, property='num_hair_neighbors')

    assert 6 * HC.count().data == SC.count().data
    assert (HC.dropna(dim='id') == 0).all()
    assert (SC.dropna(dim='id') == 1).all()