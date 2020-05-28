"""PCPTopology-model specific plot function for the state / density"""

import logging
from typing import Tuple

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.patches as mpatches

from utopya import DataManager, UniverseGroup
from utopya.plotting import UniversePlotCreator, is_plot_func, PlotHelper
from utopya.plotting import MultiversePlotCreator
from utopya.plot_funcs.basic_uni import lineplot, lineplots
from utopya.dataprocessing import transform

from ..tools import save_and_close

# -----------------------------------------------------------------------------

def plot_neighbourhood(data, *, hlpr, only_type: str, 
                       helpers_frame_hist: dict=None,
                       helpers_frame_area: dict=None,
                       **plot_kwargs):
    """ Helper function to plot the cell_neighbourhood
    """
    num_neighbors = data.sel(property='num_neighbors')
    area = data.sel(property='area')
    cell_type = data.sel(property='cell_type')

    bins = range(3, 10)

    if only_type == 'hair':
        num_neighbors = num_neighbors[cell_type == 1]
        area = area[cell_type == 1]
    elif only_type == 'support':
        num_neighbors = num_neighbors[cell_type == 2]
        area = area[cell_type == 2]
    elif only_type == 'support_hair':
        num_neighbors = data.sel(property='num_hair_neighbors')
        num_neighbors = num_neighbors[cell_type == 2]
        area = area[cell_type == 2]
        bins = range(0, 7)
    elif only_type == 'support_support':
        num_neighbors = num_neighbors - \
                        data.sel(property='num_hair_neighbors')
        num_neighbors = num_neighbors[cell_type == 2]
        area = area[cell_type == 2]
        bins = range(0, 7)
    elif only_type != 'all':
        raise ValueError("'only_type' unknown, was '{}', but must be "
                    "one of {}"
                    "".format(only_type, ['all', 'hair', 'support',
                                            'support_hair',
                                            'support_support']))

    # the histogram
    hlpr.select_axis(col=0, row=0)
    hlpr.ax.clear()

    hlpr.ax.hist(x=num_neighbors, bins=bins, density=True)

    if helpers_frame_hist:
        for name, args in helpers_frame_hist.items():
            hlpr.invoke_helper(name, **args)

    # the mean area per polygon class
    hlpr.select_axis(col=1, row=0)
    hlpr.ax.clear()

    area_mean = xr.DataArray([area.where(num_neighbors==i)\
                                    .mean().data for i in bins],
                                dims=['num_neighbors'],
                                coords={'num_neighbors': bins})
    area_std = xr.DataArray([area.where(num_neighbors==i)\
                                    .std().data for i in bins],
                                dims=['num_neighbors_std'],
                                coords={'num_neighbors_std': bins})
    
    hlpr.ax.errorbar(x=bins, y=area_mean, yerr=area_std)

    if helpers_frame_area:    
        for name, args in helpers_frame_area.items():
            hlpr.invoke_helper(name, **args)


@is_plot_func(creator_type=UniversePlotCreator, supports_animation=True)
def cell_neighbourhood(dm: DataManager, *, uni: UniverseGroup, hlpr: PlotHelper,
                       datapath = 'PCPTopology/Cells',
                       only_type: str='all',
                       helpers_frame_hist: dict=None,
                       helpers_frame_area: dict=None):
    """Performs a plot of the neighbourhood of the cells and the average area 
        per polygon class
    
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (int): The universe to use
        hlpr (PlotHelper): The PlotHelper
        datapath (str): Path to the data
        only_type (str): If only a single type of cells should be used for 
                         calculation. Can be 'all', 'hair', 'support'
        helpers_frame_hist (dict, optional): Dict passed to helper within every 
                                             frame
        helpers_frame_area (dict, optional): Dict passed to helper within every 
                                             frame
    """

    # Get the group that all datasets are in
    grp = uni['data/'+datapath]

    # Get the shape of the data
    uni_cfg = uni['cfg']

    # Prepare the figure ......................................................
    # Prepare the figure to have as many columns as there are properties
    hlpr.setup_figure(ncols=2)

    def update():
        # grp['cells'] is TimeSeriesGroup -> single dimension time
        for time in grp:
            plot_neighbourhood(grp[time], hlpr=hlpr, only_type=only_type,
                               helpers_frame_hist=helpers_frame_hist,
                               helpers_frame_area=helpers_frame_area)
            
            hlpr.select_axis(col=0, row=0)
            hlpr.invoke_helper('set_title', title="Time {}".format(time))
            
            yield

    hlpr.register_animation_update(update)


@is_plot_func(creator_type=MultiversePlotCreator, use_dag=True)
def cell_neighbourhood_mv(*, data: dict, hlpr: PlotHelper,
                          only_type: str='all',
                          helpers_frame_hist: dict=None,
                          helpers_frame_area: dict=None, **plot_kwargs):
    """A creator-averse plot function using the data transformation
    framework and the plot helper framework.

    Args:
        data: The selected and transformed data, containing specified tags.
        hlpr: The associated plot helper.
        **plot_kwargs: Passed on to matplotlib.pyplot.plot
    """

    # Prepare the figure ......................................................
    # Prepare the figure to have as many columns as there are properties
    hlpr.setup_figure(ncols=2)
    data = data['data'].stack(z=('seed', 'id')).squeeze()

    plot_neighbourhood(data, hlpr=hlpr, only_type=only_type,
                       helpers_frame_hist=helpers_frame_hist,
                       helpers_frame_area=helpers_frame_area)
