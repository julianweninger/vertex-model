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
def cell_neighbourhood(dm: DataManager, *, uni: UniverseGroup, hlpr: PlotHelper,
                       datapath = 'PCPTopology', time: int=0):
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

    # Prepare the figure ......................................................
    # Prepare the figure to have as many columns as there are properties
    hlpr.setup_figure(ncols=2)

    histogram_data = grp['Cell_neighbourhood']
    area_data = grp['Cell_size']
    area_data = area_data / area_data.sel(bin=0)

    def update():
        for time in histogram_data.time:

            hlpr.select_axis(col=0, row=0)
            hlpr.ax.clear()

            histogram = histogram_data.sel(time=time)
            hlpr.ax.bar(x=histogram.bin, height=histogram)

            hlpr.invoke_helper('set_title', title="Time {}".format(time.data))
            hlpr.invoke_helper('set_labels', x='number of neighbours',
                               y='count')


            hlpr.select_axis(col=1, row=0)
            hlpr.ax.clear()

            area = area_data.sel(time=time)
            hlpr.ax.plot(area.bin[1:], area[1:], '-s')
            hlpr.invoke_helper('set_labels', x='number of neighbours',
                               y='<$A_n$>/<A>')
            yield


    hlpr.register_animation_update(update)
