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
                       datapath: str='PCPVertex', cfgpath: str='PCPVertex',
                       time: int=0):
    """Performs a plot of the cells, edges and vertices
    
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (int): The universe to use
        hlpr (PlotHelper): The PlotHelper
    """
    def scatter_vertex(x, y, ax):
        ax.scatter(x, y, c='black')

    def plot_edge(x0, y0, x1, y1, ax):
        ax.arrow(x0, y0, x1, y1, head_width=0., head_length=0., color='black')


    # Get the group that all datasets are in
    grp = uni['data/'+datapath]
    # Get the shape of the data
    uni_cfg = uni['cfg']
    vertex_cfg = uni_cfg
    for level in cfgpath.split('/'):
        vertex_cfg = vertex_cfg[level]
    hexagon_size = vertex_cfg['hexagon_size']
    num_rows = vertex_cfg['lattice_rows']
    num_columns = vertex_cfg['lattice_columns']
    # model_cfg = uni_cfg[datapath]


    # Prepare the figure ......................................................
    # Prepare the figure to have as many columns as there are properties
    hlpr.setup_figure()

    def update():
        for time in grp['Vertex_position']:
            hlpr.ax.clear()

            v_data = grp['Vertex_position'][time]
            e_data = grp['Edge_link'][time]
            c_data = grp['Cell_position'][time]

            for v_id in v_data.id:
                v = v_data.sel(id=v_id)
                scatter_vertex(v.data[0], v.data[1], hlpr.ax)

            for e_id in e_data.id:
                e = e_data.sel(id=e_id)
                ax = v_data.sel(id=e[0].data)[0]
                ay = v_data.sel(id=e[0].data)[1]
                bx = v_data.sel(id=e[1].data)[0]
                by = v_data.sel(id=e[1].data)[1]

                if (vertex_cfg['periodic_bc']):
                    dx = bx - ax
                    dy = by - ay
                    if (dx >= 0.5 and dy >= 0.5):
                        plot_edge(ax+1., ay+1., dx-1., dy-1., hlpr.ax)
                        plot_edge(ax, ay, dx-1., dy-1., hlpr.ax)
                    elif (dx <= -0.5 and dy <= -0.5):
                        plot_edge(ax-1., ay-1., dx+1., dy+1., hlpr.ax)
                        plot_edge(ax, ay, dx+1., dy+1., hlpr.ax)
                    elif (dx >= 0.5 and dy <= -0.5):
                        plot_edge(ax+1., ay-1., dx-1., dy+1., hlpr.ax)
                        plot_edge(ax, ay, dx-1., dy+1., hlpr.ax)
                    elif (dx <= -0.5 and dy >= 0.5):
                        plot_edge(ax-1., ay+1., dx+1., dy-1., hlpr.ax)
                        plot_edge(ax, ay, dx+1., dy-1., hlpr.ax)

                    elif (dx >= 0.5):
                        plot_edge(ax+1., ay, dx-1., dy, hlpr.ax)
                        plot_edge(ax, ay, dx-1., dy, hlpr.ax)
                    elif (dx <= -0.5):
                        plot_edge(ax-1., ay, dx+1., dy, hlpr.ax)
                        plot_edge(ax, ay, dx+1., dy, hlpr.ax)
                    elif (dy >= 0.5):
                        plot_edge(ax, ay+1., dx, dy-1., hlpr.ax)
                        plot_edge(ax, ay, dx, dy-1., hlpr.ax)
                    elif (dy <= -0.5):
                        plot_edge(ax, ay-1., dx, dy+1., hlpr.ax)
                        plot_edge(ax, ay, dx, dy+1., hlpr.ax)
                    else:
                        plot_edge(ax, ay, dx, dy, hlpr.ax)

                    
                else:
                    plot_edge(ax, ay, bx-ax, by-ay, hlpr.ax)

            for c_id in c_data.id:
                c = c_data.sel(id=c_id)
                hlpr.ax.scatter(c.data[0], c.data[1], c='red')

            hlpr.invoke_helper('set_title', title="Time {}".format(time))

            if (vertex_cfg['periodic_bc']):
                hlpr.invoke_helper('set_limits', x=(-0.1,1.05), y=(-0.1,1.05))
                hlpr.ax.axvline(x=0, ymin=-0.1, ymax = 1.05, c='gray', 
                                linestyle=':')
                hlpr.ax.axvline(x=1., ymin=-0.1, ymax = 1.05, c='gray', 
                                linestyle=':')
                hlpr.ax.axhline(y=0, xmin=-0.1, xmax = 1.05, c='gray', 
                                linestyle=':')
                hlpr.ax.axhline(y=1, xmin=-0.1, xmax = 1.05, c='gray', 
                                linestyle=':')
            else:
                hlpr.invoke_helper('set_limits', x=(-0.1,1.2), y=(-0.1,1.2))
            yield


    hlpr.register_animation_update(update)
