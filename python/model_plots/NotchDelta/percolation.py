"""Implements plot functions to analyze the percolation properties of 
the NotchDelta model
"""

import logging
from typing import Union, Tuple, List
import copy

import numpy as np
import xarray as xr
import matplotlib as mpl
import os

import dantro.utils.data_ops as dops
from .data_ops import count_unique
# NOTE this will be available in dantro v13.1 _(WIP)_

from dantro.plot_creators.ext_funcs.generic import facet_grid

from utopya import DataManager, MultiverseGroup
from utopya.plotting import is_plot_func, PlotHelper, MultiversePlotCreator
from utopya.dataprocessing import transform

from utopya.plot_funcs._basic import _errorbar


# Get a logger
log = logging.getLogger(__name__)

@is_plot_func(use_dag=True, required_dag_tags=('cluster_id', 'cells', ),
              supports_animation=True)
def percolation_diagram(*, data: dict, map_coord: str=None, **plot_kwargs):
    """FacetGrid plot of the size of biggest cluster

    This functions manipulates the data to obtain the density of the largest
    cluster and passes the transformed data on to the generic dantro plot 
    function `facet_grid`.
    
    Args:
        data (dict): The data selected by the DAG framework
        map_coord (str, optional): Map this coordinate to the density of hair
            cells. Requires the dag_tag 'cells'.
        **plot_kwargs: Passed on to generic dantro plot function `facet_grid`
    """
    num_cells = data['cluster_id'].count(dim={'x', 'y'})

    clusters = dops.where(data['cluster_id'], ">", 0)
    cluster_sizes = count_unique(data['cluster_id'], dims=['x', 'y'])
    max_cluster_density = cluster_sizes.max(dim='unique') / num_cells

    if map_coord:
        if 'cells' not in data.keys():
            raise ValueError("Percolation diagram requires dag_tag 'cells' when "
                f"mapping coordinates ('map_coord'={map_coord}). Provided "
                f"dag_tags are {data.keys()}")
        
        hair_cells = dops.where(data['cells'], "==", 1)
        num_hair_cells = hair_cells.count(dim={'x', 'y'})
        density = num_hair_cells / num_cells

        max_cluster_density.assign_coords({map_coord: density})
        max_cluster_density.rename({map_coord: 'density'})

    data['data'] = max_cluster_density

    return facet_grid(data=data, **plot_kwargs)
