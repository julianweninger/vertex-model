"""Implements plot functions to analyze the percolation properties of 
the Collier model
"""

import logging
from typing import Union, Tuple, List
import copy

import math
import numpy as np
import xarray as xr
import matplotlib as mpl

from utopya.plotting import is_plot_func, PlotHelper


# Get a logger
log = logging.getLogger(__name__)

@is_plot_func(use_dag=True, required_dag_tags=('delta', 'notch', 'nicd'),
              supports_animation=True)
def protein_concentrations_per_cell(*, data: dict, hlpr: PlotHelper,
                                    cmap: str='seismic', **plot_kwargs):
    """Plots the concentrations of the proteins Notch, Delta, and NICD for every
    cell over time.
    The last value of the delta concentration of every cell is mapped to a color
    in a colormap.

    Args:
        data (dict): The data
        hlpr (PlotHelpr): The plot helper
        cmap (str): The colormap to which to map the delta concentrations
        **plot_kwargs: Passed to matplotlib.plot()

    """
    delta = data['delta']
    delta = delta.stack(id=('x', 'y'))
    notch = data['notch']
    notch = notch.stack(id=('x', 'y'))
    nicd = data['nicd']
    nicd = nicd.stack(id=('x', 'y'))

    norm = mpl.colors.Normalize(vmin=(delta.isel(time=-1).min()),
                                vmax=(delta.isel(time=-1).max()))
    cmap = mpl.cm.get_cmap(cmap)

    for id in delta.id:
        # map the last delta value of this cell to color
        color = cmap(norm((delta.sel(id=id).isel(time=-1))))

        hlpr.select_axis(col=0, row=0)
        delta.sel(id=id).plot(color=color, **plot_kwargs)

        hlpr.select_axis(col=1, row=0)
        notch.sel(id=id).plot(color=color, **plot_kwargs)

        hlpr.select_axis(col=2, row=0)
        nicd.sel(id=id).plot(color=color, **plot_kwargs)


    hlpr.select_axis(col=0, row=0)
    hlpr.provide_defaults('set_title', title=" ")
    hlpr.provide_defaults('set_labels',
                          y=dict(label=delta.name + " concentration"))
    hlpr.select_axis(col=1, row=0)
    hlpr.provide_defaults('set_title', title=" ")
    hlpr.provide_defaults('set_labels',
                          y=dict(label=notch.name + " concentration"))
    hlpr.select_axis(col=2, row=0)
    hlpr.provide_defaults('set_title', title=" ")
    hlpr.provide_defaults('set_labels',
                          y=dict(label=nicd.name + " concentration"))

    hlpr.select_axis(col=0, row=0)
    hlpr.ax.legend([mpl.lines.Line2D([0], [0], color='red', **plot_kwargs),
                    mpl.lines.Line2D([0], [0], color='blue', **plot_kwargs)],
                    ['High delta', 'Low delta'])
