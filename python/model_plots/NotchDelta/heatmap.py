"""Implements basic generic multiverse plot functions"""

import logging
from typing import Union, Tuple

import numpy as np
import xarray as xr
import matplotlib as mpl
import os

from utopya import DataManager, MultiverseGroup
from utopya.plotting import is_plot_func, PlotHelper, MultiversePlotCreator
from utopya.dataprocessing import transform

# Get a logger
log = logging.getLogger(__name__)

@is_plot_func(creator_type=MultiversePlotCreator)
def heatmap(dm: DataManager, *,
            mv_data: xr.Dataset,
            hlpr: PlotHelper,
            to_plot: str='data',
            transform_data: Tuple[Union[dict, str]]=None,
            transformations_log_level: int=10,
            **plot_kwargs):
    # Get the data and apply transformations, e.g. to reduce dimensionality
    data = mv_data[to_plot]
    
    if transform_data:
        data = transform(data, *transform_data, aux_data=mv_data,
                         log_level=transformations_log_level)

    print(data.coords)
    data.plot()
    # hlpr.ax.imshow(data, **plot_kwargs)
    