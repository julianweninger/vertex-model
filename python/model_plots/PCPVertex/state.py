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
    def scatter_vertex(x, y, ax, Lx, Ly):
        ax.scatter(Lx*x, Ly*y, c='black')

    def plot_edge(x0, y0, dx, dy, ax, Lx, Ly):
        ax.arrow(Lx*x0, Ly*y0, Lx*dx, Ly*dy, head_width=0., head_length=0.,
                 color='black')

    def plot_arrow(x0, y0, dx, dy, ax, Lx, Ly):
        ax.arrow(Lx*x0, Ly*y0, Lx*dx, Ly*dy, head_width=0.01, head_length=0.01,
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
            hlpr.ax.set_aspect('auto')

            v_data = grp['Vertex_position'][time]
            e_data = grp['Edge_link'][time]
            c_data = grp['Cell_position'][time]
            
            Lx = v_data.attrs["Lx"][0]
            Ly = v_data.attrs["Ly"][0]

            for v_id in v_data.id:
                v = v_data.sel(id=v_id)
                x = v.sel(coordinate="x")
                y = v.sel(coordinate="y")
                scatter_vertex(x, y, hlpr.ax, Lx, Ly)

            for e_id in e_data.id:
                e = e_data.sel(id=e_id)
                vertex_a = e.sel(vertex="a")
                vertex_b = e.sel(vertex="b")

                ax = v_data.sel(id=vertex_a, coordinate='x')
                ay = v_data.sel(id=vertex_a, coordinate='y')
                bx = v_data.sel(id=vertex_b, coordinate='x')
                by = v_data.sel(id=vertex_b, coordinate='y')

                if (vertex_cfg['periodic_bc']):
                    dx = bx - ax
                    dy = by - ay
                    if (dx >= 0.5 and dy >= 0.5):
                        plot_edge(ax+1., ay+1., dx-1., dy-1., hlpr.ax, Lx, Ly)
                        plot_edge(ax, ay, dx-1., dy-1., hlpr.ax, Lx, Ly)
                    elif (dx <= -0.5 and dy <= -0.5):
                        plot_edge(ax-1., ay-1., dx+1., dy+1., hlpr.ax, Lx, Ly)
                        plot_edge(ax, ay, dx+1., dy+1., hlpr.ax, Lx, Ly)
                    elif (dx >= 0.5 and dy <= -0.5):
                        plot_edge(ax+1., ay-1., dx-1., dy+1., hlpr.ax, Lx, Ly)
                        plot_edge(ax, ay, dx-1., dy+1., hlpr.ax, Lx, Ly)
                    elif (dx <= -0.5 and dy >= 0.5):
                        plot_edge(ax-1., ay+1., dx+1., dy-1., hlpr.ax, Lx, Ly)
                        plot_edge(ax, ay, dx+1., dy-1., hlpr.ax, Lx, Ly)

                    elif (dx >= 0.5):
                        plot_edge(ax+1., ay, dx-1., dy, hlpr.ax, Lx, Ly)
                        plot_edge(ax, ay, dx-1., dy, hlpr.ax, Lx, Ly)
                    elif (dx <= -0.5):
                        plot_edge(ax-1., ay, dx+1., dy, hlpr.ax, Lx, Ly)
                        plot_edge(ax, ay, dx+1., dy, hlpr.ax, Lx, Ly)
                    elif (dy >= 0.5):
                        plot_edge(ax, ay+1., dx, dy-1., hlpr.ax, Lx, Ly)
                        plot_edge(ax, ay, dx, dy-1., hlpr.ax, Lx, Ly)
                    elif (dy <= -0.5):
                        plot_edge(ax, ay-1., dx, dy+1., hlpr.ax, Lx, Ly)
                        plot_edge(ax, ay, dx, dy+1., hlpr.ax, Lx, Ly)
                    else:
                        plot_edge(ax, ay, dx, dy, hlpr.ax, Lx, Ly)

                    
                else:
                    plot_edge(ax, ay, bx-ax, by-ay, hlpr.ax, Lx, Ly)

            for c_id in c_data.id:
                c = c_data.sel(id=c_id)
                cell_type = c.sel(coordinate="cell_type")
                x = c.sel(coordinate="x")
                y = c.sel(coordinate="y")
                pol_x = c.sel(coordinate="polarity_x")
                pol_y = c.sel(coordinate="polarity_y")
                dx = pol_x / 5.
                dy = pol_y / 5.
                plot_arrow(x - dx/2, y - dy/2.,
                           dx, dy, hlpr.ax, Lx, Ly)
                if (cell_type.data == 1):
                    color = 'red'
                    hlpr.ax.scatter(Lx*x, Ly*y, c=color, s=15, alpha=0.5)
                elif (cell_type.data == 2):
                    color = 'gray'
                    hlpr.ax.scatter(Lx*x, Ly*y, c=color, s=15, alpha=0.5)

            hlpr.invoke_helper('set_title', title="Time {}".format(time))

            if (vertex_cfg['periodic_bc']):
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
