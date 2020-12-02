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
from scipy.interpolate import griddata

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
    cell_contractility = uni['data'][model_name]['Energy']['Cell_contractility']
    edge_contractility = uni['data'][model_name]['Energy']['Edge_contractility']
    areaelasticity = uni['data'][model_name]['Energy']['Areaelasticity']
    num_T1s = uni['data'][model_name]['Statistics/num_T1s']
    num_T1s_attempted = uni['data'][model_name]['Statistics/num_T1s_attempted']
    num_T2s = uni['data'][model_name]['Statistics/num_T2s']

    # Create the line plot of energy
    hlpr.ax.plot(energy.time, energy, label='total')
    hlpr.ax.plot(linetension.time, linetension, label='linetension')
    hlpr.ax.plot(cell_contractility.time, cell_contractility,
                 label='cell contractility')
    hlpr.ax.plot(edge_contractility.time, edge_contractility,
                 label='edge contractility')
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
                       property: str=None,
                       property_interpolation_kwargs: dict={},
                       property_interpolation_plot_kwargs: dict={}):
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
        property_interpolation_kwargs (dict, optional): Kwargs passed on to 
            scipy.interpolate.griddata.
        property_interpolation_plot_kwargs (dict, optional): Kwargs passed on to
            imshow on interpolated griddata,
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
            # the ids of the vertices a and b of the edges
            vertex_a = e_data.sel(property="vertex_a")
            vertex_b = e_data.sel(property="vertex_b")

            # the coordinates of vertices a and b in the set of edges
            ax = v_data.sel(id=vertex_a, property='x')
            ay = v_data.sel(id=vertex_a, property='y')
            bx = v_data.sel(id=vertex_b, property='x')
            by = v_data.sel(id=vertex_b, property='y')

            # assign the correct ids to the edges vertices a and b
            ax = ax.assign_coords(id=vertex_a.id)
            ay = ay.assign_coords(id=vertex_a.id)
            bx = bx.assign_coords(id=vertex_b.id)
            by = by.assign_coords(id=vertex_b.id)

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
                # the ids of the vertices a and b of the edges
                vertex_a = e_data.sel(property="vertex_a")
                vertex_b = e_data.sel(property="vertex_b")

                # the coordinates of vertices a and b in the set of edges
                ax_tmp = v_data.sel(id=vertex_a, property='x')
                ay_tmp = v_data.sel(id=vertex_a, property='y')
                bx_tmp = v_data.sel(id=vertex_b, property='x')
                by_tmp = v_data.sel(id=vertex_b, property='y')

                # assign the correct ids to the edges vertices a and b
                ax_tmp = ax.assign_coords(id=vertex_a.id)
                ay_tmp = ay.assign_coords(id=vertex_a.id)
                bx_tmp = bx.assign_coords(id=vertex_b.id)
                by_tmp = by.assign_coords(id=vertex_b.id)

                # use the inverse arrows
                dx_tmp = ax_tmp - bx_tmp
                dy_tmp = ay_tmp - by_tmp

                # only those edges which cross the boundaries
                # i.e. those which are not whithin the domain
                mask = np.isnan(dx_tmp.where(dx_tmp <  0.5 * Lx).where(
                                         dx_tmp > -0.5 * Lx).where(
                                         dy_tmp <  0.5 * Ly).where(
                                         dy_tmp > -0.5*Ly))                                         
                bx_tmp = bx_tmp.where(mask)
                by_tmp = by_tmp.where(mask)
                dx_tmp = dx_tmp.where(mask)
                dy_tmp = dy_tmp.where(mask)
                
                mask = dx_tmp >= 0.5 * Lx
                dx_tmp += mask * (-Lx)
                mask = dx_tmp <= -0.5 * Lx
                dx_tmp += mask * Lx
                mask = dy_tmp >= 0.5 * Ly
                dy_tmp += mask * (-Ly)
                mask = dy_tmp <= -0.5 * Ly
                dy_tmp += mask * (+Ly)
                
                hlpr.ax.quiver(bx_tmp, by_tmp, dx_tmp, dy_tmp, **quiverkwargs)


            ### plot cells
            cell_type = c_data.sel(property="cell_type")
            x = c_data.sel(property="x")
            y = c_data.sel(property="y")
            color = ['red' if d == 1 else 'grey' for d in cell_type]
            hlpr.ax.scatter(x, y, c=color, s=cell_marker_size,
                            alpha=0.5)

            # # cell polarity plots
            # pol_x = 0 # c_data.sel(property="polarity_x")
            # pol_y = 0 # c_data.sel(property="polarity_y")
            # dx = pol_x / 5.
            # dy = pol_y / 5.
            # hlpr.ax.quiver(x - dx/2, y - dy/2., dx, dy, hlpr.ax)
            
            # gather the data for property interpolation
            # for the energies of the cells
            if (   property == 'area_elasticity'
                or property == 'cell_contractility'):
                prop_data = grp['Cell_energies'][time]
                prop_data = prop_data.sel(energy_term=property)
                prop_data = prop_data.assign_coords({'x': x, 'y': y})

                grid_x, grid_y = np.mgrid[0:Lx:1000j, 0:Ly:1000j]
                grid_z1 = griddata((x.data, y.data), prop_data.data,
                                   (grid_x, grid_y))

                hlpr.ax.imshow(grid_z1.T, extent=(0,Lx,0,Ly), origin='lower')
            
            # for the energies of the edges
            elif (   property == 'linetension'
                or property == 'edge_contractility'):
                prop_data = grp['Edge_energies'][time]
                prop_data = prop_data.sel(energy_term=property)
                prop_data = prop_data.assign_coords({'x': (ax + dx / 2.),
                                                     'y': (ay + dy / 2.)})


            # for other cell data
            elif property:
                prop_data = c_data.sel(property=property)
                prop_data = prop_data.assign_coords({'x': x, 'y': y})

                if property == "hexatic_order":
                    # remove non hair cells
                    # hexatic order of non hair cells is always 0
                    prop_data = prop_data.where(cell_type == 1)

            # perform the interpolation
            if property:
                num_data_points = len(prop_data) * 10j
                grid_x, grid_y = np.mgrid[0:Lx:num_data_points,
                                        0:Ly:num_data_points]
                grid_z1 = griddata((prop_data.x, prop_data.y), prop_data,
                                    (grid_x, grid_y),
                                    **property_interpolation_kwargs)

                interpol = hlpr.ax.imshow(grid_z1.T, extent=(0,Lx,0,Ly),
                                          origin='lower',
                                          **property_interpolation_plot_kwargs)
                
                cbar = hlpr.fig.colorbar(interpol, ax=hlpr.ax, extend='both')
                cbar.set_label(label=property)
                cbar.minorticks_on()


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

            # end update here
            yield

            # remove colorbar
            if property:
                cbar.remove()


    hlpr.register_animation_update(update)
