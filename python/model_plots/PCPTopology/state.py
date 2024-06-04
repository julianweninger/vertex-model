"""PCPTopology-model specific plot function for the state / density"""

from concurrent.futures import thread
import logging
from tokenize import group
import warnings
from typing import List, Tuple

import numpy as np
from sympy import false
import xarray as xr
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import seaborn as sns

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
    
    continuous_time = uni['data'][continuous_time_path].data

    transitions_base(dm, uni=uni, hlpr=hlpr, model_name=model_name,
                     **plot_kwargs)

    if map_to_continuous_time:
        ax = hlpr.ax

        ax1 = ax.twiny()
        ax1.set_xlim(ax.get_xlim())
        
        # use same ticks on both axis
        ticks = ax.get_xticks()
        ax1.set_xticks(ticks)

        # the mapped times
        cticks = []
        for tick in ticks:
            if round(tick) != tick:
                cticks.append("")
            elif tick > continuous_time.time[-1]:
                cticks.append("")
            else:
                cticks.append("%.0f" % continuous_time.sel(time=tick))
        ax1.set_xticklabels(cticks, rotation=90)
        
        ax1.set_xlabel("Continuous time")

    if map_to_discrete_time:
        ax = hlpr.ax

        ax1 = ax.twiny()
        
        # reduce number of ticks
        num_ticks = 15 * (1 -   continuous_time[0]
                              / continuous_time[-1])
        d_len = len(continuous_time)
        slice_step = max(int(d_len/num_ticks), 1)
        
        # the continuous time for every operation, reduced in length
        ctime = continuous_time[::slice_step]
        
        # the tick labels and where to place them
        labels = ["%.0f" % time.time for time in ctime]
        times = [time.data for time in ctime]

        # append the last tick, for correct scaling
        if times[-1] != continuous_time[-1]:
            labels.append("%.0f" % continuous_time[-1].time)
            times.append(continuous_time[-1])

        ax1.set_xlim(ax.get_xlim())
        
        ax1.set_xticks(times)
        ax1.set_xticklabels(labels, rotation=90)

        ax1.set_xlabel("Time of operations")

    hlpr.fig.tight_layout()


def plot_neighbourhood(*, num_neighbors, num_hair_neighbors,
                      area, shape_index, cell_type,
                      hlpr: PlotHelper, only_type: str='all',
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

    bins = range(3, 14)

    # use num_neighbors
    if only_type == 'hair':
        num_neighbors = num_neighbors[cell_type == 1]
        area = area[cell_type == 1]
        shape_index = shape_index[cell_type == 1]
    elif only_type == 'support':
        num_neighbors = num_neighbors[cell_type == 2]
        area = area[cell_type == 2]
        shape_index = shape_index[cell_type == 2]

    # use num_hair_neighbors
    elif only_type == 'hair_hair':
        num_neighbors = num_hair_neighbors[cell_type == 1]
        bins = range(0, 7)
        plot_area = False
        plot_shape_index = False
    elif only_type == 'support_hair':
        num_neighbors = num_hair_neighbors[cell_type == 2]
        bins = range(0, 7)
        plot_area = False
        plot_shape_index = False

    # use num_neighbors - num_hair_neighbors
    elif only_type == 'hair_support':
        num_neighbors = num_neighbors - num_hair_neighbors
        num_neighbors = num_neighbors[cell_type == 1]
        bins = range(0, 7)
        plot_area = False
        plot_shape_index = False
    elif only_type == 'support_support':
        num_neighbors = num_neighbors - num_hair_neighbors
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


@is_plot_func(creator_type=UniversePlotCreator, supports_animation=True,
              use_dag=True,
              required_dag_tags=(
                  'num_neighbors',
                  'area',
                  'shape_index',
                  'cell_type'),
                compute_only_required_dag_tags=False
)
def cell_neighbourhood(*, data: dict, hlpr: PlotHelper,
                       stack_dims: list=None,
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
    num_neighbors = data['num_neighbors']
    if 'num_hair_neighbors' in data:
        num_hair_neighbors = data['num_hair_neighbors']
    else:
        num_hair_neighbors = None
    area = data['area']
    shape_index = data['shape_index']
    cell_type = data['cell_type']
    if stack_dims is not None:
        for dim in stack_dims:
            if not dim in num_neighbors.dims:
                raise ValueError("Dimension {} not in data {}", dim,
                                 num_neighbors.dims)
        stack_dims.append('id')
        num_neighbors = num_neighbors.stack(ids=stack_dims).squeeze()
        if num_hair_neighbors is not None:
            num_hair_neighbors = num_hair_neighbors.stack(ids=stack_dims)\
                                                   .squeeze()
        area = area.stack(ids=stack_dims).squeeze()
        shape_index = shape_index.stack(ids=stack_dims).squeeze()
        cell_type = cell_type.stack(ids=stack_dims).squeeze()

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
        for time in num_neighbors.time:
            plot_neighbourhood(
                num_neighbors=num_neighbors.sel(time=time),
                num_hair_neighbors=num_hair_neighbors.sel(time=time) if \
                    num_hair_neighbors is not None else None,
                area=area.sel(time=time),
                shape_index=shape_index.sel(time=time),
                cell_type=cell_type.sel(time=time),
                hlpr=hlpr, only_type=only_type,
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


@is_plot_func(creator_type=MultiversePlotCreator,
              use_dag=True,
              required_dag_tags=(
                  'hair_cells',
                  'data_experiment'
              ))
def coordination_vs_area(*, data: dict, hlpr: PlotHelper):
    theory = data['hair_cells']

    theory = theory.sort_values(by='area')
    if len(theory['Gamma'].unique()) > 2:
        raise ValueError('Gamma can take a max of 2 values!')

    sns.lineplot(data=theory.loc[theory['Gamma'] == 0.],
                 x='area',
                 y='num_neighbors',
                 color='blue',
                 ci='sd',
                 label='Homogeneous',
    )

    sns.lineplot(data=theory.loc[theory['Gamma'] > 0.],
                 x='area',
                 y='num_neighbors',
                 color='red',
                 ci='sd',
                 label='Heterogeneous',
    )

    experiment = data['data_experiment']
    experiment = experiment.loc[~experiment['is_border_cell']]
    experiment = experiment.loc[experiment['is_HC']]
    experiment = experiment.loc[experiment['stage'] != 'E9']
    experiment = experiment.loc[experiment['stage'] != 'E13']
    experiment = experiment.loc[experiment['stage'] != 'E15']

    stage_ordering = ['E8', 'E10', 'E12', 'E14']


    data = experiment.groupby(
        by=['stage', 'position', 'SI_position', 'PD_position', 'file_id']
    )
    data = data[['HC_normalized_area', 'num_neighbors']].mean()
    data = data.reset_index()

    data['SI_position_long'] = data['SI_position'].where(
        data['SI_position'] != 'S', 'Superior')
    data['SI_position_long'] = data['SI_position_long'].where(
        data['SI_position'] != 'None', 'Mid')
    data['SI_position_long'] = data['SI_position_long'].where(
        data['SI_position'] != 'I', 'Inferior')

    data['Stage'] = data['stage']
    data['S-I position'] = data['SI_position_long']
    data['P-D position'] = data['PD_position']

    sns.scatterplot(data=data,
                    x='HC_normalized_area', y='num_neighbors',
                    hue='Stage', hue_order=stage_ordering,
                    palette='copper',
                    style='P-D position',
                    size='S-I position',
                    size_order=['Superior', 'Mid', 'Inferior'],
                    # markers=True,
                    # ls='',
                    # dashes=False,
                    # ci=None,
    )

    data = experiment.groupby(by='file_id')

    hlpr.ax.errorbar(x=data['normalized_area_cells'].mean(),
                     y=data['num_neighbors'].mean(),
                     xerr=data['normalized_area_cells'].std(),
                     yerr=data['num_neighbors'].std(),
                     ls='', color='black', 
                     linewidth=0.2, alpha=0.8,
                     zorder=0)


@is_plot_func(creator_type=MultiversePlotCreator,
              use_dag=True
             )
def histogram_plot(
    *, data: dict, hlpr: PlotHelper,
    data_tag: str,
    label: str=None,
    data_tag_split: str=None,
    label_split: str=None,
    x: str,
    x_order: List=None,
    x_ticklabels: List=None,
    y: str,
    y_label: str=None,
    y_range: Tuple,
    color: str,
    color_split: str=None,
    plot_means_kwargs: dict=None,
    ):
    """A plot that creates violin plots on discrete data using histogram
    representation.

    Args:
        - data (dict): The data
        - hlpr (PlotHelper): PlotHelper for figure configuration
        - data_tag (str): The key to the dataset in data
        - label (str, optional): The label to data_tag. If none, the key is used
        - data_tag_split (str, optional): Analogous to data_tag. 
            If provided, the histogram is split vertically. 
            Left, corresponds to data with data_tag.
            Right, data with data_tag_split.
        - label_split (str, optional): Analogous to label for data_tag_split
        - x (str): The dimension within data mapped to x-axis
        - x_order (List, optional): An ordering of positions on x-axis
        - x_ticklabels (List, optional): Manually assign labels to x-axis ticks
        - y (str): The dimension within dataset(s) mapped to the y-axis. 
            A frequency count is performed on this dimension.
        - y_label (str, optional): A label for the y-axis
        - y_range (Tuple): A range-expression for datapoints along y used for
            binning. Data needs to be discrete to these points. Forwarded to
            np.arange(*y_range)
        - color (str): The color used for the histogram
        - color_split(str, optional): Color used for split histogram (right)
        - plot_means_kwargs (dict, optional): If given, mean of histogram will
            be indicated and kwargs forwarded to plt.hlines().
    """

    if data_tag is None:
        raise RuntimeError("No data tag provided!")
            

    if data_tag not in data.keys():
        raise RuntimeError("Data tag '{}' not found in data! "
                            "Available tags: {}",
                            data_tag,
                            data.keys())
    if data_tag_split is not None and data_tag_split not in data.keys():
        raise RuntimeError("Data tag '{}' not found in data! "
                            "Available tags: {}",
                            data_tag_split,
                            data.keys())

    x_majorticks = []
    if x_order is not None:
        for _x in x_order:
            x_majorticks.append(_x)
    else:
        for _x in data[data_tag][x].unique():
            x_majorticks.append(_x)

    def x_to_locx(x):
        return 2 * x_majorticks.index(str(x)) + 1

    y_ticks = np.arange(*y_range)
    y_min = y_range[0]
    y_max = y_range[1]

    hlpr.provide_defaults(
        'set_limits',
        x=(0, 2*len(x_majorticks)),
        y=(y_min, y_max-1)
    )
    hlpr.provide_defaults(
        'set_ticks',
        x=dict({
            'major': {
                'locs': list(np.arange(1, 2*len(x_majorticks)+1, 2)),
                'labels': x_ticklabels if x_ticklabels is not None else x_majorticks
            },
            'minor': {
                'locs': list(np.arange(0, 2*len(x_majorticks)+1, 2)),
            },
        })
    )
    

    def __plot_histograms(data, *, 
                         x_groupby: str, groupby_order: List=None,
                         y_data_column: str,
                         plot_left: bool, color: str,
                         plot_means_kwargs: dict=None):
        """Plot histograms

        data: the raw data to create the histogram from
        x_groupby: How to group the data. Every group will be one major x-tick.
        groupby_order: Order of x-ticks
        y_data_column: Column in data to plot -- the y coordinate.
        
        plot_left (bool): If true, plot histogram to the left of the x-tick. 
            Else, plot to the right of x-tick.
        plot_means_kwargs (dict, optional): Forwarded to ax hline with y-position
            at data mean
        """

        norm = data.groupby(by=x_groupby)[y_data_column].count()

        hist = [data.loc[np.round(data[y_data_column]) == i].groupby(by=x_groupby)[y_data_column].count()
                for i in y_ticks]
        hist = pd.concat(hist, axis=1, keys=y_ticks)
        hist = hist.divide(norm, axis=0)

        hist = hist.reset_index()
        hist[x_groupby] = hist[x_groupby].astype("category")
        if groupby_order is not None:
            hist[x_groupby].cat.set_categories(groupby_order)
            hist = hist.sort_values(by=x_groupby)


        sign = 1
        if plot_left:
            sign=-1

        boxes = []
        for grp in hist[x_groupby]:
            d = hist.loc[hist[x_groupby] == grp]
            for N in y_ticks:
                boxes.append(
                    Rectangle(
                        (x_to_locx(grp), N - 0.5),
                        width=sign*min(d[N].values[0], 0.975),
                        height=1,
                        ls='', lw=None
                    )
                )
        
        pc = PatchCollection(boxes, facecolor=color, edgecolor='None')
        hlpr.ax.add_collection(pc)

        if plot_means_kwargs is not None:
            for grp in hist[x_groupby]:
                raw_data = data.loc[data[x_groupby] == grp]
                hist_data = hist.loc[hist[x_groupby] == grp]

                mean = raw_data[y_data_column].mean()
                count = sign*min(hist_data[np.round(mean)].values[0], 0.975)

                hlpr.ax.hlines(
                    y=mean,
                    xmin=x_to_locx(grp),
                    xmax=x_to_locx(grp) + count,
                    **plot_means_kwargs
                )



    __plot_histograms(
        data[data_tag],
        x_groupby=x,
        groupby_order=x_order,
        y_data_column=y,
        plot_left=True,
        color=color,
        plot_means_kwargs=plot_means_kwargs
    )

    if data_tag_split is not None:
        __plot_histograms(
            data[data_tag_split],
            x_groupby=x,
            groupby_order=x_order,
            y_data_column=y,
            plot_left=False,
            color=color_split,
            plot_means_kwargs=plot_means_kwargs
        )
    else:
        __plot_histograms(
            data[data_tag],
            x_groupby=x,
            groupby_order=x_order,
            y_data_column=y,
            plot_left=False,
            color=color,
            plot_means_kwargs=plot_means_kwargs
        )

    if y_label is None:
        y_label = '{} (Frequency)'.format(y)
    hlpr.provide_defaults('set_labels', y=y_label)

    if data_tag_split is not None:
        if label is None:
            label = data_tag
        if label_split is None:
            label_split = data_tag_split
        legend_elements = [
            mpl.patches.Patch(facecolor=color, edgecolor='None', label=label),
            mpl.patches.Patch(facecolor=color_split, edgecolor='None',
                              label=label_split)
        ]
        hlpr.ax.legend(handles=legend_elements, loc='upper left')
    elif label is not None:
        legend_elements = [
            mpl.patches.Patch(facecolor=color, edgecolor='None', label=label),
            mpl.patches.Patch(facecolor=color_split, edgecolor='None',
                              label=label_split)
        ]
        hlpr.ax.legend(handles=legend_elements, loc='upper left')



@is_plot_func(creator_type=MultiversePlotCreator,
              use_dag=True,
              required_dag_tags=(
                  'data_theory',
                  'mapper_theory',
                  'data_experiment'
              ))
def histogram_map_experiment(
        *, data: dict, hlpr: PlotHelper, 
        x: str='stage',
        mapper_experiment: str='normalized_area_cells',
        mapper_theory: str='area',
        **kwargs
    ):
    """Performs a 'histogram_plot' where area coordinates in hair_cell data are
        selected closest to data_experiment HC area.
    """
    
    # theory = pd.concat([data['data_theory'],
    #                     [mapper_theory]],
    #                     axis=1, join='outer')
    theory = data['data_theory']

    experiment = data['data_experiment']

    mapping = pd.DataFrame()
    mapping['experimental'] = experiment.groupby(by=x)[mapper_experiment].mean()

    theory_area = pd.DataFrame()
    theory_area['theoretical'] = data['mapper_theory'][mapper_theory]
        
    for _x in experiment[x]:
        theory_area['difference_to_' + _x] = abs(
              theory_area['theoretical']
            - mapping.loc[_x, 'experimental']
        )
        mapping.loc[_x, 'closest_theoretical'] = theory_area.loc[
            theory_area['difference_to_' + _x].idxmin(),
            'theoretical'
        ]
        
        dt = theory_area['difference_to_' + _x].min()
        mapping.loc[_x, 'difference'] = dt
        time = theory_area['difference_to_' + _x].idxmin()
        mapping.loc[_x, 'theory_time'] = time

        theory.loc[theory['time'] == time, x] = _x
    
    theory = theory.dropna(subset=x)

    print(mapping)

    data['Experiment'] = experiment
    data['Theory'] = theory

    histogram_plot(
        data=data,
        hlpr=hlpr,
        data_tag='Experiment',
        data_tag_split='Theory',
        x=x,
        **kwargs
    )


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
        cmap = plt.cm.get_cmap(cmap)
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
        cnt = d.sel(property="count")
        
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
            prop = (prop * cnt).sum(dim=average_dim) / cnt.sum(dim=average_dim)
            prop = prop.squeeze()
            if std is not None:
                std = (  (std * cnt).sum(dim=average_dim)
                       / cnt.sum(dim=average_dim)).squeeze()
            if min is not None:
                min = min.min(dim=average_dim).squeeze()
            if max is not None:
                max = max.max(dim=average_dim).squeeze()
   
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
        cmap = plt.cm.get_cmap(cmap)
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
    