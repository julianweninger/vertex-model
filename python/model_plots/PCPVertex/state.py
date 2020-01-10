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
    def scatter_vertex(x, y, ax):
        ax.scatter(x, y, c='black')

    def plot_edge(x0, y0, x1, y1, ax):
        ax.arrow(x0, y0, x1, y1, head_width=0., head_length=0., color='black')


    # Get the group that all datasets are in
    grp = uni['data/'+datapath]

    # Get the shape of the data
    uni_cfg = uni['cfg']
    vertex_cfg = uni_cfg
    for level in datapath.split('/'):
        vertex_cfg = vertex_cfg[level]
    hexagon_size = vertex_cfg['hexagon_size']
    num_rows = vertex_cfg['lattice_rows']
    num_columns = vertex_cfg['lattice_columns']
  # [4 * sqrt(3) * size_hexagonal, 1.5 * num_columns * size_hexagonal]
    if (vertex_cfg['periodic_bc']):
        periodic_limits = [num_rows * 3**0.5 * hexagon_size, 
                           num_columns * 1.5 * hexagon_size]
    else:
        periodic_limits = None
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
                scatter_vertex(v.data[0], v.data[1], hlpr.ax)

            for e_id in e_data.dim_1:
                e = e_data.sel(dim_1=e_id)
                ax = v_data.sel(dim_1=e[0])[0]
                ay = v_data.sel(dim_1=e[0])[1]
                bx = v_data.sel(dim_1=e[1])[0]
                by = v_data.sel(dim_1=e[1])[1]

                if (periodic_limits):
                    dx = bx - ax
                    dy = by - ay
                    Lx = periodic_limits[0]
                    Ly = periodic_limits[1]
                    if (dx >= Lx/2. and dy >= Ly/2.):
                        plot_edge(ax+Lx, ay+Ly, dx-Lx, dy-Ly, hlpr.ax)
                        plot_edge(ax, ay, dx-Lx, dy-Ly, hlpr.ax)
                    elif (dx <= -Lx/2. and dy <= -Ly/2.):
                        plot_edge(ax-Lx, ay-Ly, dx+Lx, dy+Ly, hlpr.ax)
                        plot_edge(ax, ay, dx+Lx, dy+Ly, hlpr.ax)
                    elif (dx >= Lx/2. and dy <= -Ly/2.):
                        plot_edge(ax+Lx, ay-Ly, dx-Lx, dy+Ly, hlpr.ax)
                        plot_edge(ax, ay, dx-Lx, dy+Ly, hlpr.ax)
                    elif (dx <= -Lx/2. and dy >= Ly/2.):
                        plot_edge(ax-Lx, ay+Ly, dx+Lx, dy-Ly, hlpr.ax)
                        plot_edge(ax, ay, dx+Lx, dy-Ly, hlpr.ax)

                    elif (dx >= Lx/2.):
                        plot_edge(ax+Lx, ay, dx-Lx, dy, hlpr.ax)
                        plot_edge(ax, ay, dx-Lx, dy, hlpr.ax)
                    elif (dx <= -Lx/2.):
                        plot_edge(ax-Lx, ay, dx+Lx, dy, hlpr.ax)
                        plot_edge(ax, ay, dx+Lx, dy, hlpr.ax)
                    elif (dy >= Ly/2.):
                        plot_edge(ax, ay+Ly, dx, dy-Ly, hlpr.ax)
                        plot_edge(ax, ay, dx, dy-Ly, hlpr.ax)
                    elif (dy <= -Ly/2.):
                        plot_edge(ax, ay-Ly, dx, dy+Ly, hlpr.ax)
                        plot_edge(ax, ay, dx, dy+Ly, hlpr.ax)
                    else:
                        plot_edge(ax, ay, dx, dy, hlpr.ax)

                    
                else:
                    plot_edge(ax, ay, bx-ax, by-ay, hlpr.ax)

            for c_id in c_data.dim_1:
                c = c_data.sel(dim_1=c_id)
                hlpr.ax.scatter(c.data[0], c.data[1], c='red')

            hlpr.invoke_helper('set_title', title="Time {}".format(step))
            yield


    hlpr.register_animation_update(update)
