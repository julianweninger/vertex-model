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
    
    continuous_time = uni['data'][continuous_time_path]

    if map_to_continuous_time:
        times = ["%.0f" % continuous_time.sel(time=time) \
                 for time in continuous_time.time]
        ax1 = hlpr.ax.twiny()
        ax1.set_xlim(hlpr.ax.get_xlim())
        ax1.set_xticks(continuous_time.time)
        ax1.set_xticklabels(times, rotation=90)
        ax1.set_xlabel("Continuous time")

    if map_to_discrete_time:
        ax1 = hlpr.ax.twiny()
        ax1.set_xlim(hlpr.ax.get_xlim())
        ax1.set_xticks(continuous_time.data)
        ax1.set_xticklabels(["%.0f" % time for time in continuous_time.time],
                            rotation=90)

        ax1.set_xlabel("Time of operations")

    hlpr.select_axis(0, 0)

    transitions_base(dm, uni=uni, hlpr=hlpr, model_name=model_name,
                     **plot_kwargs)

    hlpr.fig.tight_layout()


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

    bins = range(3, 12)

    # use num_neighbors
    if only_type == 'hair':
        num_neighbors = data.sel(property='num_neighbors')
        num_neighbors = num_neighbors[cell_type == 1]
        area = area[cell_type == 1]
        shape_index = shape_index[cell_type == 1]
    elif only_type == 'support':
        num_neighbors = data.sel(property='num_neighbors')
        num_neighbors = num_neighbors[cell_type == 2]
        area = area[cell_type == 2]
        shape_index = shape_index[cell_type == 2]

    # use num_hair_neighbors
    elif only_type == 'hair_hair':
        num_neighbors = data.sel(property='num_hair_neighbors')
        num_neighbors = num_neighbors[cell_type == 1]
        bins = range(0, 7)
        plot_area = False
        plot_shape_index = False
    elif only_type == 'support_hair':
        num_neighbors = data.sel(property='num_hair_neighbors')
        num_neighbors = num_neighbors[cell_type == 2]
        bins = range(0, 7)
        plot_area = False
        plot_shape_index = False

    # use num_neighbors - num_hair_neighbors
    elif only_type == 'hair_support':
        num_neighbors = num_neighbors - \
                        data.sel(property='num_hair_neighbors')
        num_neighbors = num_neighbors[cell_type == 1]
        bins = range(0, 7)
        plot_area = False
        plot_shape_index = False
    elif only_type == 'support_support':
        num_neighbors = num_neighbors - \
                        data.sel(property='num_hair_neighbors')
        num_neighbors = num_neighbors[cell_type == 2]
        bins = range(0, 7)
        plot_area = False
        plot_shape_index = False
    elif only_type != 'all':
        raise ValueError("'only_type' unknown, was '{}', but must be "
                    "one of {}"
                    "".format(only_type, ['all', 'hair', 'support',
                                          'support_hair',
                                          'hair_support',
                                          'support_support',
                                          'hair_hair']))
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


def _errorbar(*, hlpr: PlotHelper, data: xr.DataArray, std: xr.DataArray,
              min: xr.DataArray=None, max: xr.DataArray=None,
              min_max_kwargs: dict=None,
              fill_between: bool=True, fill_between_kwargs: dict=None,
              **errorbar_kwargs):
    """Given the data and (optionally) the standard deviation data, plots a
    single errorbar line.
    
    Args:
        hlpr (PlotHelper): The helper
        data (xr.DataArray): The data
        std (xr.DataArray): The y-error data
        fill_between (bool, optional): Whether to use plt.fill_between or
            plt.errorbar to plot y-errors
        fill_between_kwargs (dict, optional): Passed on to plt.fill_between
        **errorbar_kwargs: Passed on to plt.errorbar
    
    Raises:
        ValueError: On non-1D data
    """
    # Check dimensionality
    if data.ndim != 1:
        raise ValueError("Requiring 1D data to plot a single errorbar "
                         "line but got {}D data with shape {}:\n{}\n"
                         "Apply dimensionality reducing transformations "
                         "using the `transform_data` argument to arrive "
                         "at plottable data."
                         "".format(data.ndim, data.shape, data))

    elif std is not None and std.ndim != 1:
        raise ValueError("Requiring 1D standard deviation data to plot the "
                         "error markers of a single errorbar line but "
                         "got {}D data with shape {}:\n{}\n"
                         "Apply dimensionality reducing transformations "
                         "using the `transform_std` argument to arrive "
                         "at plottable data."
                         "".format(std.ndim, std.shape, std))
    elif min is not None and min.ndim != 1:
        raise ValueError("Requiring 1D minimum data to plot the "
                         "error markers of a single errorbar line but "
                         "got {}D data with shape {}:\n{}\n"
                         "Apply dimensionality reducing transformations "
                         "using the `transform_std` argument to arrive "
                         "at plottable data."
                         "".format(min.ndim, min.shape, min))
    elif max is not None and max.ndim != 1:
        raise ValueError("Requiring 1D maximum data to plot the "
                         "error markers of a single errorbar line but "
                         "got {}D data with shape {}:\n{}\n"
                         "Apply dimensionality reducing transformations "
                         "using the `transform_std` argument to arrive "
                         "at plottable data."
                         "".format(max.ndim, max.shape, max))

    # Data is ok.
    # Decide on whether yerr is done by errorbar or by fill_between
    yerr = std if not fill_between else None

    # Plot the data against its coordinates, including standard deviation
    ebar = hlpr.ax.errorbar(data.coords[data.dims[0]], data,
                            yerr=yerr, **errorbar_kwargs)

    # Now plot the confidence interval via 
    if fill_between and std is not None:
        # Find out the colour of the error bar line. Get line collection
        lc, _, _ = ebar
        line_color = lc.get_c()
        line_alpha = lc.get_alpha() if lc.get_alpha() else 1.
        line_label = errorbar_kwargs.get('label', None)

        # Prepare kwargs
        fb_kwargs = (copy.deepcopy(fill_between_kwargs) if fill_between_kwargs
                     else {})

        if 'color' not in fb_kwargs:
            fb_kwargs['color'] = line_color
        if 'alpha' not in fb_kwargs:
            fb_kwargs['alpha'] = line_alpha * .2
        if 'label' not in fb_kwargs and line_label:
            fb_kwargs['label'] = line_label + " (std. dev.)"

        # Fill.
        hlpr.ax.fill_between(data.coords[data.dims[0]],
                             y1=(data - std), y2=(data + std),
                             **fb_kwargs)

    for m, id in zip([min, max], ['min', 'max']):
        if m is not None:
            # Find out the colour of the error bar line. Get line collection
            lc, _, _ = ebar
            line_color = lc.get_c()
            line_alpha = lc.get_alpha() if lc.get_alpha() else 1.
            line_label = errorbar_kwargs.get('label', None)

            # Prepare kwargs
            fb_kwargs = (copy.deepcopy(min_max_kwargs) if min_max_kwargs
                        else {})

            if 'color' not in fb_kwargs:
                fb_kwargs['color'] = line_color
            if 'alpha' not in fb_kwargs:
                fb_kwargs['alpha'] = line_alpha * .5
            if line_label:
                fb_kwargs['label'] = line_label + " (" + id + ".)"
            if 'linestyle' not in fb_kwargs and 'ls' not in fb_kwargs:
                fb_kwargs['linestyle'] = '--'

            hlpr.ax.plot(m.coords[data.dims[0]], m, **fb_kwargs)

    # TODO Manually add the legend patch

@is_plot_func(use_dag=True)
def errorbars(*, data: dict, to_plot: dict, hlpr: PlotHelper, property: str,
              average_dim: str=None, cmap: str=None,**errorbar_kwargs):
    """Perform errorbar plots from the selected multiverse data.
    
    This plot, ultimately, requires 1D data, where the remaining dimension is
    plotted on the x-axis. The ``transform_data`` or ``lines_from`` arguments
    can be used to work with higher-dimensional data.

    Creates datasets hair_cells, hair_cells__std, support_cells,
    and support_cells__std in data.
    
    Args:
        data (dict): The data.
        to_plot (dict): A dict of specifications of lines to plot. 
            The keys must be available in data. If key + '__std' is available
            in data, this data is used for errorbars, otherwise simple lineplot
            performed.
            Mapped values can contain these kwargs
                std (str, default: `<name>__std`): name of the std dataset
                plot_std (bool, default: True): Plot the std dataset if
                    available
            The remaining mapped values are passed on to plt.errorbar.
        cmap (str, optional): If given, the lines created from ``to_plot``
            will be colored according to this color map.
        **errorbar_kwargs: Passed on to plt.errorbar
    """

    num_lines = len(data.keys())

    if num_lines > 1:
        hlpr.provide_defaults('set_legend', use_legend=True)

    # Determine whether there will be colours according to a color map
    if cmap is not None:
        cmap = mpl.cm.get_cmap(cmap)
        colors = [cmap(i/max(num_lines-1, 1)) for i in range(num_lines)]
    else:
        colors = [None] * num_lines

    # Iterate over the plot specifications
    for (key, plot_spec), color in zip(to_plot.items(), colors):
        # Prepare additional kwargs
        add_kwargs = dict()

        if color is not None and 'color' not in plot_spec:
            add_kwargs['color'] = color

        if 'label' not in plot_spec:
            add_kwargs['label'] = key

        plot_std = plot_spec.pop('plot_std', True)
        plot_min_max = plot_spec.pop('plot_min_max', False)
        d = data[key]
        prop = d.sel(property=property)
        if plot_std and (property + '__stddev') in d.property.data:
            std = d.sel(property=property+'__stddev')
        else:
            std = None
        if plot_min_max and (property + '__min') in d.property.data:
            min = d.sel(property=property+'__min')
        else:
            min = None
        if plot_min_max and (property + '__max') in d.property.data:
            max = d.sel(property=property+'__max')
        else:
            max = None

        if average_dim:
            if not average_dim in d.coords:
                raise KeyError("No dimension '{}' in data '{}'. "
                               "Available dims: {}"
                               "".format(average_dim, d.name, d.dims))
            prop = prop.sum(dim=average_dim) / prop.count(dim=average_dim)
            if std is not None:
                std = std.sum(dim=average_dim) / std.count(dim=average_dim)
                std = std.squeeze()
            if min is not None:
                min = min.sum(dim=average_dim) / min.count(dim=average_dim)
                min = min.squeeze()
            if max is not None:
                max = max.sum(dim=average_dim) / max.count(dim=average_dim)
                max = max.squeeze()
            prop = prop.squeeze()
   
        _errorbar(hlpr=hlpr, data=prop, std=std, min=min, max=max, **plot_spec,
                  **add_kwargs, **errorbar_kwargs)

@is_plot_func(use_dag=True)
def errorbars_mv(*, data: dict, to_plot: dict, hlpr: PlotHelper,
             cmap: str=None,**errorbar_kwargs):
    """Perform errorbar plots from the selected multiverse data.
    
    This plot, ultimately, requires 1D data, where the remaining dimension is
    plotted on the x-axis. The ``transform_data`` or ``lines_from`` arguments
    can be used to work with higher-dimensional data.

    Creates datasets hair_cells, hair_cells__std, support_cells,
    and support_cells__std in data.
    
    Args:
        data (dict): The data.
        to_plot (dict): A dict of specifications of lines to plot. 
            The keys must be available in data. If key + '__std' is available
            in data, this data is used for errorbars, otherwise simple lineplot
            performed.
            Mapped values can contain these kwargs
                std (str, default: `<name>__std`): name of the std dataset
                plot_std (bool, default: True): Plot the std dataset if
                    available
            The remaining mapped values are passed on to plt.errorbar.
        cmap (str, optional): If given, the lines created from ``to_plot``
            will be colored according to this color map.
        **errorbar_kwargs: Passed on to plt.errorbar
    """

    num_lines = len(data.keys())

    if num_lines > 1:
        hlpr.provide_defaults('set_legend', use_legend=True)

    # Determine whether there will be colours according to a color map
    if cmap is not None:
        cmap = mpl.cm.get_cmap(cmap)
        colors = [cmap(i/max(num_lines-1, 1)) for i in range(num_lines)]
    else:
        colors = [None] * num_lines

    # Iterate over the plot specifications
    for (key, plot_spec), color in zip(to_plot.items(), colors):
        # Prepare additional kwargs
        add_kwargs = dict()

        if color is not None and 'color' not in plot_spec:
            add_kwargs['color'] = color

        if 'label' not in plot_spec:
            add_kwargs['label'] = key

        std=plot_spec.pop('std', None)
        plot_std = plot_spec.pop('plot_std', True)
        if (key + '__std') in data and plot_std and std is None:
            std = data[key + '__std']
   
        _errorbar(hlpr=hlpr, data=data[key], std=std, **plot_spec,
                  **add_kwargs, **errorbar_kwargs)

@is_plot_func(use_dag=True, required_dag_tags=('x', 'y', ))
def scatter_xy_data(*, data: dict, hlpr: PlotHelper,
                    plot_time_color_gradient: bool=False,
                    **plot_kwargs):
    """A scatter-plot of x over y data. 

    Args:
        plot_time_color_gradient (bool, default: false): Whether to use a color
            gradient using the time as 3rd coordinate.
    """
    if plot_time_color_gradient:
        c = data['x'].time
    else:
        c = plot_kwargs.pop("c", None)
        
    hlpr.ax.scatter(data['x'], data['y'], c=c, **plot_kwargs)
    