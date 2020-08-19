"""PCPTopology-model specific plot function for the state / density"""

import logging
import warnings
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

from ..tools import save_and_close
from ..PCPVertex.state import transitions as transitions_base

# -----------------------------------------------------------------------------

@is_plot_func(creator_type=UniversePlotCreator)
def transitions(dm: DataManager, *, uni: UniverseGroup, hlpr: PlotHelper,
                model_name: str='PCPTopolgy',
                map_to_continuous_time: bool=False,
                map_to_discrete_time: bool=False,
                continuous_time_path: str='PCPTopology/Energy/Continuous_time',
                **plot_kwargs):
    """Performs a plot of the T1 and T2 transitions over time together with 
    the energy
    
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (UniverseGroup): The data for this universe
        hlpr (PlotHelper): The PlotHelper
        model_name (str): The name of the model the data resides in
        path_to_data (str or Tuple[str, str]): The path to the data within the
            model data or the paths to the x and the y data, respectively
        transform_data (dict, optional): Transformations to apply to the data.
            This can be used for dimensionality reduction of the data, but
            also for other operations, e.g. to selecting a slice.
            For available parameters, see
            :py:func:`utopya.dataprocessing.transform`
        transformations_log_level (int, optional): The log level of all the
            transformation operations.
        **plot_kwargs: Passed on to plt.plot
    
    Raises:
        ValueError: On invalid data dimensionality
        ValueError: On mismatch of data shapes
    """
    transitions_base(dm, uni=uni, hlpr=hlpr, model_name=model_name,
                     **plot_kwargs)
    
    continuous_time = uni['data'][continuous_time_path]

    if map_to_continuous_time:
        ax3 = hlpr.ax.twiny()
        ax3.set_xlim(hlpr.ax.get_xlim())
        ax3.set_xticks(continuous_time.time)
        ax3.set_xticklabels(["%.0f" % continuous_time.sel(time=time) for time 
                                in continuous_time.time])

        ax3.set_xlabel("Continuous time")

    if map_to_discrete_time:
        ax3 = hlpr.ax.twiny()
        ax3.set_xlim(hlpr.ax.get_xlim())
        ax3.set_xticks(continuous_time.data)
        ax3.set_xticklabels(["%.0f" % time for time in continuous_time.time])

        ax3.set_xlabel("Time of operations")

def plot_neighbourhood(data, *, hlpr: PlotHelper, only_type: str='all',
                       plot_hist: bool=True,
                       plot_area: bool=True,
                       plot_shape_index: bool=True,
                       helpers_frame_hist: dict=None,
                       helpers_frame_area: dict=None,
                       helpers_frame_shape_index: dict=None,
                       hist_plot_kwargs: dict=None,
                       area_plot_kwargs: dict=None,
                       shape_index_plot_kwargs: dict=None):
    """ Helper function to plot the cell_neighbourhood.

    Plots the properties averaged separately for the polygon classes

    Args:
        data: the data
        hlpr (PlotHelper): The PlotHelper

        only_type (str): If only a single type of cells should be used for 
                         calculation. Can be 'all', 'hair', 'support'

        plot_hist (bool, default: True): Whether to plot the histogram
        plot_area (bool, default: True): Whether to plot the mean area
        plot_shape_index (bool, default: True): Whether to plot the mean shape
                                                index
        helpers_frame_hist (dict, optional): Dict passed to helper within every 
                                             frame
        helpers_frame_area (dict, optional): Dict passed to helper within every 
                                             frame
        helpers_frame_shape_index (dict, optional): Dict passed to helper within every 
                                             frame
        hist_plot_kwargs: passed on to matplotlib.hist (histogram plot)
        area_plot_kwargs: passed on to matplotlib.errorbar (area plot)
        shape_index_plot_kwargs: passed on to matplotlib.errorbar
                                 (shape_index plot)
    """
    num_neighbors = data.sel(property='num_neighbors')
    area = data.sel(property='area')
    shape_index = data.sel(property='shape_index')
    cell_type = data.sel(property='cell_type')

    bins = range(3, 10)

    if only_type == 'hair':
        num_neighbors = num_neighbors[cell_type == 1]
        area = area[cell_type == 1]
        shape_index = shape_index[cell_type == 1]
    elif only_type == 'support':
        num_neighbors = num_neighbors[cell_type == 2]
        area = area[cell_type == 2]
        shape_index = shape_index[cell_type == 2]
    elif only_type == 'support_hair':
        num_neighbors = data.sel(property='num_hair_neighbors')
        num_neighbors = num_neighbors[cell_type == 2]
        area = area[cell_type == 2]
        shape_index = shape_index[cell_type == 2]
        bins = range(0, 7)
    elif only_type == 'support_support':
        num_neighbors = num_neighbors - \
                        data.sel(property='num_hair_neighbors')
        num_neighbors = num_neighbors[cell_type == 2]
        area = area[cell_type == 2]
        shape_index = shape_index[cell_type == 2]
        bins = range(0, 7)
    elif only_type != 'all':
        raise ValueError("'only_type' unknown, was '{}', but must be "
                    "one of {}"
                    "".format(only_type, ['all', 'hair', 'support',
                                            'support_hair',
                                            'support_support']))
    figure_index = 0

    # the histogram
    if plot_hist:
        hlpr.select_axis(col=figure_index, row=0)
        hlpr.ax.clear()
        figure_index += 1

        if (not hist_plot_kwargs):
            hist_plot_kwargs = {}
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning,
                                    message="invalid value encountered in "
                                            "true_divide")
            hlpr.ax.hist(x=num_neighbors, bins=bins, **hist_plot_kwargs)

        if helpers_frame_hist:
            for name, args in helpers_frame_hist.items():
                hlpr.invoke_helper(name, **args)

    # the mean area per polygon class
    if plot_area:
        hlpr.select_axis(col=figure_index, row=0)
        hlpr.ax.clear()
        figure_index += 1

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning,
                                    message="Mean of empty slice")
            warnings.filterwarnings("ignore", category=RuntimeWarning,
                                    message="Degrees of freedom <= 0 for slice.")
            area_mean = xr.DataArray([area.where(num_neighbors==i)\
                                            .mean().data for i in bins],
                                        dims=['num_neighbors'],
                                        coords={'num_neighbors': bins})
            area_std = xr.DataArray([area.where(num_neighbors==i)\
                                            .std().data for i in bins],
                                        dims=['num_neighbors_std'],
                                        coords={'num_neighbors_std': bins})
        
        if (not area_plot_kwargs):
            area_plot_kwargs = {}
        hlpr.ax.errorbar(x=bins, y=area_mean, yerr=area_std, **area_plot_kwargs)

        if helpers_frame_area:    
            for name, args in helpers_frame_area.items():
                hlpr.invoke_helper(name, **args)

    # the mean shape index per polygon class
    if plot_shape_index:
        hlpr.select_axis(col=figure_index, row=0)
        hlpr.ax.clear()
        figure_index += 1

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning,
                                    message="Mean of empty slice")
            warnings.filterwarnings("ignore", category=RuntimeWarning,
                                    message="Degrees of freedom <= 0 for slice.")
            shape_index_mean = xr.DataArray([shape_index.where(num_neighbors==i)\
                                            .mean().data for i in bins],
                                        dims=['num_neighbors'],
                                        coords={'num_neighbors': bins})
            shape_index_std = xr.DataArray([shape_index.where(num_neighbors==i)\
                                            .std().data for i in bins],
                                        dims=['num_neighbors_std'],
                                        coords={'num_neighbors_std': bins})
        
        if (not shape_index_plot_kwargs):
            shape_index_plot_kwargs = {}
        hlpr.ax.errorbar(x=bins, y=shape_index_mean, yerr=shape_index_std,
                        **shape_index_plot_kwargs)

        if helpers_frame_shape_index:    
            for name, args in helpers_frame_shape_index.items():
                hlpr.invoke_helper(name, **args)


@is_plot_func(creator_type=UniversePlotCreator, supports_animation=True)
def cell_neighbourhood(dm: DataManager, *, uni: UniverseGroup, hlpr: PlotHelper,
                       datapath = 'PCPTopology/Cells',
                       only_type: str='all',
                       plot_hist: bool=True,
                       plot_area: bool=True,
                       plot_shape_index: bool=True,
                       helpers_frame_hist: dict=None,
                       helpers_frame_area: dict=None,
                       helpers_frame_shape_index: dict=None,
                       hist_plot_kwargs: dict=None,
                       area_plot_kwargs: dict=None,
                       shape_index_plot_kwargs: dict=None):
    """Performs a plot of the neighbourhood of the cells and the average area 
        per polygon class
    
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (int): The universe to use
        hlpr (PlotHelper): The PlotHelper
        datapath (str): Path to the data
        only_type (str): If only a single type of cells should be used for 
                         calculation. Can be 'all', 'hair', 'support'

        plot_hist (bool, default: True): Whether to plot the histogram
        plot_area (bool, default: True): Whether to plot the mean area
        plot_shape_index (bool, default: True): Whether to plot the mean shape
                                                index
        helpers_frame_hist (dict, optional): Dict passed to helper within every 
                                             frame
        helpers_frame_area (dict, optional): Dict passed to helper within every 
                                             frame
        helpers_frame_shape_index (dict, optional): Dict passed to helper within every 
                                             frame
        hist_plot_kwargs: passed on to matplotlib.hist (histogram plot)
        area_plot_kwargs: passed on to matplotlib.errorbar (area plot)
        shape_index_plot_kwargs: passed on to matplotlib.errorbar (shape_index plot)
    """

    # Get the group that all datasets are in
    grp = uni['data/'+datapath]

    # Get the shape of the data
    uni_cfg = uni['cfg']

    # Prepare the figure ......................................................
    # Prepare the figure to have as many columns as there are properties
    num_subplots = 3
    for prop in [plot_hist, plot_area, plot_shape_index]:
        if not prop:
            num_subplots -= 1
    if num_subplots <= 0:
        raise RuntimeError("Nothing to plot!")
    hlpr.setup_figure(ncols=num_subplots)


    def update():
        # grp['cells'] is TimeSeriesGroup -> single dimension time
        for time in grp:
            plot_neighbourhood(grp[time], hlpr=hlpr, only_type=only_type,
                               plot_hist=plot_hist, plot_area=plot_area,
                               plot_shape_index=plot_shape_index,
                               helpers_frame_hist=helpers_frame_hist,
                               helpers_frame_area=helpers_frame_area,
                               helpers_frame_shape_index=helpers_frame_shape_index,
                               hist_plot_kwargs=hist_plot_kwargs, 
                               area_plot_kwargs=area_plot_kwargs,
                               shape_index_plot_kwargs=shape_index_plot_kwargs)
            
            hlpr.select_axis(col=0, row=0)
            hlpr.invoke_helper('set_title', title="Time {}".format(time))
            
            yield

    hlpr.register_animation_update(update)


@is_plot_func(creator_type=MultiversePlotCreator, use_dag=True)
def cell_neighbourhood_mv(*, data: dict, hlpr: PlotHelper,
                          only_type: str='all',
                          plot_hist: bool=True,
                          plot_area: bool=True,
                          plot_shape_index: bool=True,
                          helpers_frame_hist: dict=None,
                          helpers_frame_area: dict=None,
                          helpers_frame_shape_index: dict=None,
                          hist_plot_kwargs: dict=None,
                          area_plot_kwargs: dict=None,
                          shape_index_plot_kwargs: dict=None):
    """Plot properties of the cell averaged separately for different polygon
    classes.
    Average furthermore over the dimension seed of the data.

    Args:
        only_type (str): If only a single type of cells should be used for 
                         calculation. Can be 'all', 'hair', 'support'

        plot_hist (bool, default: True): Whether to plot the histogram
        plot_area (bool, default: True): Whether to plot the mean area
        plot_shape_index (bool, default: True): Whether to plot the mean shape
                                                index
        helpers_frame_hist (dict, optional): Dict passed to helper within every 
                                             frame
        helpers_frame_area (dict, optional): Dict passed to helper within every 
                                             frame
        helpers_frame_shape_index (dict, optional): Dict passed to helper within every 
                                             frame
        hist_plot_kwargs: passed on to matplotlib.hist (histogram plot)
        area_plot_kwargs: passed on to matplotlib.errorbar (area plot)
        shape_index_plot_kwargs: passed on to matplotlib.errorbar (shape_index plot)
    """

    # Prepare the figure ......................................................
    # Prepare the figure to have as many columns as there are properties
    hlpr.setup_figure(ncols=2)
    data = data['data'].stack(z=('seed', 'id')).squeeze()

    plot_neighbourhood(data, hlpr=hlpr, only_type=only_type,
                       plot_hist=plot_hist, plot_area=plot_area,
                       plot_shape_index=plot_shape_index,
                       helpers_frame_hist=helpers_frame_hist,
                       helpers_frame_area=helpers_frame_area,
                       helpers_frame_shape_index=helpers_frame_shape_index,
                       hist_plot_kwargs=hist_plot_kwargs, 
                       area_plot_kwargs=area_plot_kwargs,
                       shape_index_plot_kwargs=shape_index_plot_kwargs)
