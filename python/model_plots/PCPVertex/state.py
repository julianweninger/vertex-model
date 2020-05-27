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
                       time: int=0, plot_vertices: bool=False):
    """Performs a plot of the cells, edges and vertices
    
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (int): The universe to use
        hlpr (PlotHelper): The PlotHelper
    """
    def scatter_vertex(x, y, ax):
        ax.scatter(x, y, c='black')

    def plot_edge(x0, y0, dx, dy, ax):
        ax.arrow(x0, y0, dx, dy, head_width=0., head_length=0.,
                 color='black')

    def plot_arrow(x0, y0, dx, dy, ax):
        ax.arrow(x0, y0, dx, dy, head_width=0.01, head_length=0.01,
                 color='black')


    def adjustFigAspect(fig,aspect=1):
        '''
        Adjust the subplot parameters so that the figure has the correct
        aspect ratio.
        '''
        xsize,ysize = fig.get_size_inches()
        minsize = min(xsize,ysize)
        xlim = .4*minsize/xsize
        ylim = .4*minsize/ysize
        if aspect < 1:
            xlim *= aspect
        else:
            ylim /= aspect
        fig.subplots_adjust(left=.5-xlim,
                            right=.5+xlim,
                            bottom=.5-ylim,
                            top=.5+ylim)


    # Get the group that all datasets are in
    grp = uni['data/'+datapath]
    # Get the shape of the data
    uni_cfg = uni['cfg']
    vertex_cfg = uni_cfg
    for level in cfgpath.split('/'):
        vertex_cfg = vertex_cfg[level]
    # model_cfg = uni_cfg[datapath]


    # Prepare the figure ......................................................
    # Prepare the figure to have as many columns as there are properties
    hlpr.setup_figure()

    def update():
        for time in grp['Vertices']:
            hlpr.ax.clear()
            hlpr.ax.set_aspect('auto')

            v_data = grp['Vertices'][time]
            e_data = grp['Edges'][time]
            c_data = grp['Cells'][time]
            
            Lx = v_data.attrs["Lx"][0]
            Ly = v_data.attrs["Ly"][0]

            if plot_vertices:
                for v_id in v_data.id:
                    v = v_data.sel(id=v_id)
                    x = v.sel(property="x")
                    y = v.sel(property="y")
                    scatter_vertex(x, y, hlpr.ax)

            for e_id in e_data.id:
                e = e_data.sel(id=e_id)
                vertex_a = e.sel(property="vertex_a")
                vertex_b = e.sel(property="vertex_b")

                ax = v_data.sel(id=vertex_a, property='x')
                ay = v_data.sel(id=vertex_a, property='y')
                bx = v_data.sel(id=vertex_b, property='x')
                by = v_data.sel(id=vertex_b, property='y')

                if (vertex_cfg['space']['periodic']):
                    dx = bx - ax
                    dy = by - ay
                    if (dx >= 0.5 * Lx and dy >= 0.5 * Ly):
                        plot_edge(ax+Lx, ay+Ly, dx-Lx, dy-Ly, hlpr.ax)
                        plot_edge(ax, ay, dx-Lx, dy-Ly, hlpr.ax)
                    elif (dx <= -0.5 * Lx and dy <= -0.5 * Ly):
                        plot_edge(ax-Lx, ay-Ly, dx+Lx, dy+Ly, hlpr.ax)
                        plot_edge(ax, ay, dx+Lx, dy+Ly, hlpr.ax)
                    elif (dx >= 0.5 * Lx and dy <= -0.5 * Ly):
                        plot_edge(ax+Lx, ay-Ly, dx-Lx, dy+Ly, hlpr.ax)
                        plot_edge(ax, ay, dx-Lx, dy+Ly, hlpr.ax)
                    elif (dx <= -0.5 * Lx and dy >= 0.5 * Ly):
                        plot_edge(ax-Lx, ay+Ly, dx+Lx, dy-Ly, hlpr.ax)
                        plot_edge(ax, ay, dx+Lx, dy-Ly, hlpr.ax)

                    elif (dx >= 0.5 * Lx):
                        plot_edge(ax+Lx, ay, dx-Lx, dy, hlpr.ax)
                        plot_edge(ax, ay, dx-Lx, dy, hlpr.ax)
                    elif (dx <= -0.5 * Lx):
                        plot_edge(ax-Lx, ay, dx+Lx, dy, hlpr.ax)
                        plot_edge(ax, ay, dx+Lx, dy, hlpr.ax)
                    elif (dy >= 0.5 * Ly):
                        plot_edge(ax, ay+Ly, dx, dy-Ly, hlpr.ax)
                        plot_edge(ax, ay, dx, dy-Ly, hlpr.ax)
                    elif (dy <= -0.5 * Ly):
                        plot_edge(ax, ay-Ly, dx, dy+Ly, hlpr.ax)
                        plot_edge(ax, ay, dx, dy+Ly, hlpr.ax)
                    else:
                        plot_edge(ax, ay, dx, dy, hlpr.ax)

                    
                else:
                    plot_edge(ax, ay, bx-ax, by-ay, hlpr.ax)

            for c_id in c_data.id:
                c = c_data.sel(id=c_id)
                cell_type = c.sel(property="cell_type")
                x = c.sel(property="x")
                y = c.sel(property="y")
                pol_x = 0 # c.sel(property="polarity_x")
                pol_y = 0 # c.sel(property="polarity_y")
                dx = pol_x / 5.
                dy = pol_y / 5.
                plot_arrow(x - dx/2, y - dy/2., dx, dy, hlpr.ax)
                if (cell_type.data == 1):
                    color = 'red'
                    hlpr.ax.scatter(x, y, c=color, s=15, alpha=0.5)
                elif (cell_type.data == 2):
                    color = 'gray'
                    hlpr.ax.scatter(x, y, c=color, s=15, alpha=0.5)

            hlpr.invoke_helper('set_title', title="Time {}".format(time))

            if (vertex_cfg['space']['periodic']):
                hlpr.invoke_helper('set_limits', x=(-0.1,Lx+.05), y=(-0.1,Ly+.05))
                hlpr.ax.axvline(x=0, ymin=-0.1, ymax = Ly +.05, c='gray', 
                                linestyle=':')
                hlpr.ax.axvline(x=Lx, ymin=-0.1, ymax = Ly +.05, c='gray', 
                                linestyle=':')
                hlpr.ax.axhline(y=0, xmin=-0.1, xmax = Lx+.05, c='gray', 
                                linestyle=':')
                hlpr.ax.axhline(y=Ly, xmin=-0.1, xmax = Lx+.05, c='gray', 
                                linestyle=':')
            else:
                hlpr.invoke_helper('set_limits', x=(-0.1,Lx+.2), y=(-0.1,Ly+.2))
            
            ratio = Ly / Lx
            xleft, xright = hlpr.ax.get_xlim()
            ybottom, ytop = hlpr.ax.get_ylim()
            hlpr.ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
            
            yield


    hlpr.register_animation_update(update)
