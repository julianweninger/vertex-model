"""SavannaHeterogeneous-model specific plot function for the state / density"""

import logging
from typing import Tuple

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.patches as mpatches

from utopya import DataManager, UniverseGroup
from utopya.plotting import UniversePlotCreator, is_plot_func, PlotHelper
from utopya.plot_funcs.basic_uni import lineplot, lineplots
from utopya.dataprocessing import transform

from ..tools import save_and_close

# -----------------------------------------------------------------------------

@is_plot_func(creator_type=UniversePlotCreator, supports_animation=True)
def cellular_structure(dm: DataManager, *, uni: UniverseGroup, hlpr: PlotHelper,
                       datapath = 'PCPVertex', time: int=0):
    """Performs a plot of the cells, edges and vertices
    
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (int): The universe to use
        hlpr (PlotHelper): The PlotHelper
    """
    # Get the group that all datasets are in
    grp = uni['data/'+datapath]

    # Get the shape of the data
    uni_cfg = uni['cfg']
    # model_cfg = uni_cfg[datapath]

    # Prepare the figure ......................................................
    # Prepare the figure to have as many columns as there are properties
    hlpr.setup_figure()

    for id in grp['vertices']['0'].dim_1:
        v = grp['vertices']['0'].sel(dim_1=id)
        hlpr.ax.scatter(v.data[0], v.data[1])

    steps = []
    for step in grp['vertices']:
        steps.append(int(step))
    steps.sort()

    def update():
        for step in steps:
            hlpr.ax.clear()

            v_data = grp['vertices'][str(step)]
            e_data = grp['edges'][str(step)]
            c_data = grp['cells'][str(step)]

            for v_id in v_data.dim_1:
                v = v_data.sel(dim_1=v_id)
                hlpr.ax.scatter(v.data[0], v.data[1], c='black')

            for e_id in e_data.dim_1:
                e = e_data.sel(dim_1=e_id)
                ax = v_data.sel(dim_1=e[0])[0]
                ay = v_data.sel(dim_1=e[0])[1]
                bx = v_data.sel(dim_1=e[1])[0]
                by = v_data.sel(dim_1=e[1])[1]
                hlpr.ax.arrow(ax, ay, bx-ax, by-ay, 
                         head_width=0., head_length=0., color='black')

            for c_id in c_data.dim_1:
                c = c_data.sel(dim_1=c_id)
                hlpr.ax.scatter(c.data[0], c.data[1], c='red')

            hlpr.invoke_helper('set_title', title="Time {}".format(step))
            # Done with this frame; yield control to the animation framework
            # which will grab the frame...
            # hlpr.invoke_helper('set_limits', x=[0, 3.5], y=[0, 3])
            yield


    hlpr.register_animation_update(update)
