"""Tests of the output of the PCPVertex model"""

import logging
import math

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

def test_hexagonal_setup():
    def run_sim(*, hexagon_shape: str, num_rows: int=20, num_cols: int=10):
        mv, dm = mtc.create_run_load(parameter_space={
            'PCPVertex': {'agent_manager': {'setup_params': { 'hexagonal': {
                'hexagon_shape': hexagon_shape,
                'lattice_rows': num_rows,
                'lattice_columns': num_cols}}
            }}})

        return dm['multiverse'][0]['data']['PCPVertex']['Cells'][0]

    for shape in ['pointy', 'pointy_top', 'flat', 'flat_top']:
        cells = run_sim(hexagon_shape=shape)
        area = cells.sel(property='area')
        perimeter = cells.sel(property='perimeter')
        q_x = cells.sel(property='q_x')
        q_y = cells.sel(property='q_y')

        assert abs(area.mean().data - 1) < 1.e-3
        assert abs(perimeter.mean().data - 6*(2 / 3**1.5)**0.5) < 1.e-3
        assert abs(q_x.mean().data) < 1.e-3
        assert abs(q_y.mean().data) < 1.e-3

        Lx = cells.attrs['Lx'][0]
        Ly = cells.attrs['Ly'][0]
        assert abs(Lx / Ly - 10. / 20.) < 0.25 * (perimeter.mean().data / 6)

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