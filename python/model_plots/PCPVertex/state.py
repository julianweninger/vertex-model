"""SavannaHeterogeneous-model specific plot function for the state / density"""

import logging
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

from ..tools import save_and_close

# -----------------------------------------------------------------------------


@is_plot_func(creator_type=UniversePlotCreator)
def transitions(dm: DataManager, *, uni: UniverseGroup, hlpr: PlotHelper,
                model_name: str='PCPVertex', **plot_kwargs):
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
    # Get the data
    energy = uni['data'][model_name]['Energy']['Total']
    linetension = uni['data'][model_name]['Energy']['Linetension']
    contractility = uni['data'][model_name]['Energy']['Contractility']
    areaelasticity = uni['data'][model_name]['Energy']['Areaelasticity']
    num_T1s = uni['data'][model_name]['Statistics/num_T1s']
    num_T1s_attempted = uni['data'][model_name]['Statistics/num_T1s_attempted']
    num_T2s = uni['data'][model_name]['Statistics/num_T2s']

    # Create the line plot of energy
    hlpr.ax.plot(energy.time, energy, label='total')
    hlpr.ax.plot(linetension.time, linetension, label='linetension')
    hlpr.ax.plot(contractility.time, contractility, label='contractility')
    hlpr.ax.plot(areaelasticity.time, areaelasticity, label='areaelasticity')

    hlpr.ax.set_xlabel("Time [steps]")
    hlpr.ax.set_ylabel("Energy [a.u.]")
    hlpr.ax.legend(loc='upper left')
    hlpr.ax.set_xlim(left=energy.time[0], right=energy.time[-1])

    # Plot the T1s
    ax2 = hlpr.ax.twinx()
    ax2.plot(num_T1s_attempted.time, num_T1s_attempted, color='gray', label='#T1 attempted')
    ax2.plot(num_T1s.time, num_T1s, color='black', label='#T1')
    ax2.plot(num_T2s.time, num_T2s, color='seagreen', label='#T2')

    ax2.set_ylabel("Transitions")
    ax2.set_ylim(bottom=0)
    ax2.legend(loc='upper right')


@is_plot_func(creator_type=UniversePlotCreator, supports_animation=True)
def cellular_structure(dm: DataManager, *, uni: UniverseGroup, hlpr: PlotHelper,
                       datapath: str='PCPVertex', cfgpath: str='PCPVertex',
                       cell_marker_size: int=60,
                       plot_vertices: bool=False,
                       property=None, property_marker_size=30,
                       property_cmap='autumn'):
    """Performs a plot of the cells, edges and vertices
    
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (int): The universe to use
        hlpr (PlotHelper): The PlotHelper
        datapath (str): Path to the cell manager's data
        cfgpath (str): Path to the vertex model's configuration
        cell_marker_size (int, default 60): Marker size for the cell-centers.
            Only used for HCs (red) and SCs (gray).
        plot_vertices (bool, default: false): Whether to plot the vertices
        property (str, optional): An additional cell property to plot. Data is 
            cell data where property=property.
        property_marker_size (int, default 30): The marker size for plot of
            the cell's property
        property_cmap (str, default 'autumn'): A colormap to use with the cells
            properties.
    """
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


            ### plot vertices
            if plot_vertices:
                ax.scatter(v_data.sel(property="x"), v_data.sel(property="y"),
                           c="black")


            ### plot edges
            vertex_a = e_data.sel(property="vertex_a")
            vertex_b = e_data.sel(property="vertex_b")

            ax = v_data.sel(id=vertex_a, property='x')
            ay = v_data.sel(id=vertex_a, property='y')
            bx = v_data.sel(id=vertex_b, property='x')
            by = v_data.sel(id=vertex_b, property='y')

            dx = bx - ax
            dy = by - ay

            # map edges crossing periodic boundary
            if (vertex_cfg['space']['periodic']):
                mask = dx >= 0.5 * Lx
                dx += mask * (-Lx)
                mask = dx <= -0.5 * Lx
                dx += mask * Lx
                mask = dy >= 0.5 * Ly
                dy += mask * (-Ly)
                mask = dy <= -0.5 * Ly
                dy += mask * (+Ly)
            
            quiverkwargs = dict(headlength=0., headaxislength=0., headwidth=0.,
                                scale=1, scale_units='xy',
                                color='black')
            hlpr.ax.quiver(ax, ay, dx, dy, **quiverkwargs)

            # need duplicates of periodic edges
            if (vertex_cfg['space']['periodic']):
                vertex_a = e_data.sel(property="vertex_a")
                vertex_b = e_data.sel(property="vertex_b")

                ax = v_data.sel(id=vertex_a, property='x')
                ay = v_data.sel(id=vertex_a, property='y')
                bx = v_data.sel(id=vertex_b, property='x')
                by = v_data.sel(id=vertex_b, property='y')

                # use the inverse arrows
                dx = ax - bx
                dy = ay - by

                # only those edges which cross the boundaries
                # i.e. those which are not whithin the domain
                mask = np.isnan(dx.where(dx <  0.5 * Lx).where(
                                         dx > -0.5 * Lx).where(
                                         dy <  0.5 * Ly).where(
                                         dy > -0.5*Ly))                                         
                bx = bx.where(mask)
                by = by.where(mask)
                dx = dx.where(mask)
                dy = dy.where(mask)
                
                mask = dx >= 0.5 * Lx
                dx += mask * (-Lx)
                mask = dx <= -0.5 * Lx
                dx += mask * Lx
                mask = dy >= 0.5 * Ly
                dy += mask * (-Ly)
                mask = dy <= -0.5 * Ly
                dy += mask * (+Ly)
                
                hlpr.ax.quiver(bx, by, dx, dy, **quiverkwargs)


            ### plot cells
            cell_type = c_data.sel(property="cell_type")
            x = c_data.sel(property="x")
            y = c_data.sel(property="y")
            color = ['red' if d == 1 else 'grey' for d in cell_type ]
            hlpr.ax.scatter(x, y, c=color, s=cell_marker_size,
                            alpha=0.5)

            # # cell polarity plots
            # pol_x = 0 # c_data.sel(property="polarity_x")
            # pol_y = 0 # c_data.sel(property="polarity_y")
            # dx = pol_x / 5.
            # dy = pol_y / 5.
            # hlpr.ax.quiver(x - dx/2, y - dy/2., dx, dy, hlpr.ax)
            
            if property:
                data_copy = c_data

                if property == "hexatic_order":
                    # remove non hair cells
                    # hexatic order of non hair cells is always 0
                    data_copy = data_copy.where(cell_type == 1)

                hlpr.ax.scatter(data_copy.sel(property='x'),
                                data_copy.sel(property='y'),
                                c=data_copy.sel(property=property),
                                s=property_marker_size,
                                cmap=property_cmap,
                                alpha=0.5)


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
            
            hlpr.ax.set_aspect('equal')

            yield


    hlpr.register_animation_update(update)
