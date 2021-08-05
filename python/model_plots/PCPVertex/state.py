"""SavannaHeterogeneous-model specific plot function for the state / density"""

import logging
from typing import Tuple, Union

from math import floor, ceil
import numpy as np
from numpy.lib.function_base import select
import xarray as xr
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

from utopya import DataManager, UniverseGroup
from utopya.plotting import UniversePlotCreator, is_plot_func, PlotHelper

from ..tools import save_and_close

log = logging.getLogger(__name__)
# -----------------------------------------------------------------------------


@is_plot_func(creator_type=UniversePlotCreator)
def transitions(dm: DataManager, *, uni: UniverseGroup, hlpr: PlotHelper,
                model_name: str='PCPVertex', 
                cumsum_transitions: bool=True, **plot_kwargs):
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
    
    transitions = uni['data'][model_name]['Energy/transitions']
    num_T1s = transitions.sel(property="num_T1s")
    num_T1s_attempted = transitions.sel(property="num_T1s_attempted")
    num_T2s = transitions.sel(property="num_T2s")
    if cumsum_transitions:
        num_T1s = num_T1s.cumsum()
        num_T1s_attempted = num_T1s_attempted.cumsum()
        num_T2s = num_T2s.cumsum()


    # Create the line plot of energy
    ax1 = hlpr.ax
    ax1.plot(energy.time, energy,
             label='total', alpha=0.5)
    ax1.plot(linetension.time, linetension,
             label='linetension', alpha=0.5)
    ax1.plot(cell_contractility.time, cell_contractility,
             label='cell contractility', alpha=0.5)
    ax1.plot(edge_contractility.time, edge_contractility,
             label='edge contractility', alpha=0.5)
    ax1.plot(areaelasticity.time, areaelasticity, 
             label='areaelasticity', alpha=0.5)

    ax1.set_xlabel("Time [steps]")
    ax1.set_ylabel("Energy [a.u.]")
    ax1.legend(loc='upper left')
    ax1.set_xlim(left=energy.time[0], right=energy.time[-1])

    # Plot the T1s
    ax2 = ax1.twinx()
    
    ax2.plot(num_T1s_attempted.time, num_T1s_attempted, color='gray', label='#T1 attempted')
    ax2.plot(num_T1s.time, num_T1s, color='black', label='#T1')

    ax2.set_ylabel("T1 transitions")
    ax2.set_ylim(bottom=0)
    ax2.legend(loc='upper right')

    # Plot the T2s
    ax3 = ax1.twinx()
    
    ax3.plot(num_T2s.time, num_T2s, color='seagreen', label='#T2')
    
    ax3.set_ylabel("T2 transitions")
    ax3.set_ylim(bottom=0)

    ax3.spines['right'].set_position(('outward', 60))
    ax3.yaxis.label.set_color('seagreen')

    hlpr.fig.tight_layout()


@is_plot_func(creator_type=UniversePlotCreator, supports_animation=True)
def cellular_structure(dm: DataManager, *, uni: UniverseGroup, hlpr: PlotHelper,
                       datapath: str='PCPVertex', cfgpath: str='PCPVertex',
                       select_times: list=None,
                       cell_marker_size: int=60,
                       plot_vertices: bool=False,
                       property: str=None,
                       property_hair_cells_only: bool=False,
                       property_support_cells_only: bool=False,
                       property_bulk_cells_only: bool=False,
                       property_interpolation_kwargs: dict={},
                       property_interpolation_plot_kwargs: dict={},
                       edge_property: str=None,
                       quiver_kwargs: dict=None,
                       polarity_HC: bool=False):
    """Performs a plot of the cells, edges and vertices
    
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (int): The universe to use
        hlpr (PlotHelper): The PlotHelper
        datapath (str): Path to the cell manager's data
        cfgpath (str): Path to the vertex model's configuration
        select_times (list, optional): Timepoints to select for plotting.
            If None, all timepoints are plotted
        cell_marker_size (int, default 60): Marker size for the cell-centers.
            Only used for HCs (red) and SCs (gray).
        plot_vertices (bool, default: false): Whether to plot the vertices
        property (str, optional): An additional cell property to plot. Data is 
            cell data where property=property.
        property_hair_cells_only (bool, default: false): Whether to plot only 
            for hair cells
        property_support_cells_only (bool, default: false): Whether to plot 
            only for support cells
        property_bulk_cells_only (bool, default: false): Whether to plot only 
            for non-boundary cells
        property_interpolation_kwargs (dict, optional): Kwargs passed on to 
            scipy.interpolate.griddata.
        property_interpolation_plot_kwargs (dict, optional): Kwargs passed on to
            imshow on interpolated griddata,
        edge_property (std, optional): An additional edge property to plot.
            Data is edge data where property=edge_property or
            Edge_energy where energy_term=edge_property if edge_property
            is one of the edge energies.
        quiver_kwargs (dict, optional): Updates the quiver kwargs used.
        polarity_HC (bool, default: false): Whether to plot the polarity of HC
            as a vector at the HC center of mass.
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
        # the domain extent for non periodic boundaries (centered on (0., 0.))
        domain_size_min_x = 0.
        domain_size_max_x = 0.
        domain_size_min_y = 0.
        domain_size_max_y = 0.

        if (not vertex_cfg['space']['periodic']):
            domain_size_min_x =  1000000.
            domain_size_max_x = -1000000.
            domain_size_min_y =  1000000.
            domain_size_max_y = -1000000.
            for time in grp['Vertices']:
                domain_size_min_x = min(
                    domain_size_min_x,
                    grp['Vertices'][time].sel(property="x").min())
                domain_size_max_x = max(
                    domain_size_max_x,
                    grp['Vertices'][time].sel(property="x").max())
                domain_size_min_y = min(
                    domain_size_min_y,
                    grp['Vertices'][time].sel(property="y").min())
                domain_size_max_y = max(
                    domain_size_max_y,
                    grp['Vertices'][time].sel(property="y").max())

        if select_times is not None:
            times = [str(time) for time in select_times]
        else:
            times = [time for time in grp['Vertices']]

        for time in times:
            hlpr.ax.clear()
            hlpr.ax.set_aspect('auto')

            
            if not time in grp['Vertices']:
                log.warning("Requested time {} is not available in data. "
                            "Available timepoints: {}."
                            "Continuing ..."
                            "".format(time,
                                      [int(time) for time in grp['Vertices']]))
                hlpr.invoke_helper('set_title', title="Time {}".format(time))
                yield
                continue
            
            v_data = grp['Vertices'][time]
            e_data = grp['Edges'][time]
            c_data = grp['Cells'][time]
            
            Lx = v_data.attrs["Lx"][0]
            Ly = v_data.attrs["Ly"][0]
            if (vertex_cfg['space']['periodic']):
                domain_size_max_x = Lx
                domain_size_max_y = Ly


            ### plot vertices
            if plot_vertices:
                hlpr.ax.scatter(v_data.sel(property="x"),
                                v_data.sel(property="y"),
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
                mask = (dx >= 0.5 * Lx)
                dx += mask * (-Lx)
                mask = (dx <= -0.5 * Lx)
                dx += mask * Lx
                mask = (dy >= 0.5 * Ly)
                dy += mask * (-Ly)
                mask = (dy <= -0.5 * Ly)
                dy += mask * (+Ly)
            
            quiver_args = [ax, ay, dx, dy]
            _quiver_kwargs = dict(headlength=0., headaxislength=0.,
                                  headwidth=0., scale=1, scale_units='xy',
                                  color='black')
            if quiver_kwargs:
                _quiver_kwargs.update(quiver_kwargs)
            
            # for the energies of the edges
            if (   edge_property == 'linetension'
                or edge_property == 'edge_contractility'):
                e_prop_data = grp['Edge_energies'][time]
                e_prop_data = e_prop_data.sel(energy_term=edge_property).squeeze()
                e_prop_data = e_prop_data.assign_coords({'x': (ax + dx / 2.),
                                                         'y': (ay + dy / 2.)})

                # append coloring
                quiver_args.append(e_prop_data)

            # for other cell data
            elif edge_property:
                e_prop_data = e_data.sel(property=edge_property)
                e_prop_data = e_prop_data.assign_coords({'x': (ax + dx / 2.),
                                                         'y': (ay + dy / 2.)})

                # append coloring
                quiver_args.append(e_prop_data)

            quiver = hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)

            if (len(quiver_args) == 5):
                cbar = hlpr.fig.colorbar(quiver, ax=hlpr.ax, extend='both')
                cbar.set_label(label=edge_property)
                cbar.minorticks_on()

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
                bx_tmp_prime = bx_tmp
                by_tmp_prime = by_tmp

                # use the inverse arrows
                dx_tmp = ax_tmp - bx_tmp
                dy_tmp = ay_tmp - by_tmp
                dx_tmp_prime = dx_tmp
                dy_tmp_prime = dy_tmp

                # only those edges which cross the boundaries
                # i.e. those which are not whithin the domain
                mask = np.isnan(dx_tmp.where(abs(dx_tmp) <  0.5 * Lx)\
                                      .where(abs(dy_tmp) <  0.5 * Ly))
                bx_tmp = bx_tmp.where(mask)
                by_tmp = by_tmp.where(mask)
                dx_tmp = dx_tmp.where(mask)
                dy_tmp = dy_tmp.where(mask)
                
                mask = (dx_tmp >= 0.5 * Lx)
                dx_tmp += mask * (-Lx)
                mask = (dx_tmp <= -0.5 * Lx)
                dx_tmp += mask * Lx
                mask = (dy_tmp >= 0.5 * Ly)
                dy_tmp += mask * (-Ly)
                mask = (dy_tmp <= -0.5 * Ly)
                dy_tmp += mask * (+Ly)
                
                quiver_args = [bx_tmp, by_tmp, dx_tmp, dy_tmp]
                if edge_property:
                    quiver_args.append(e_prop_data)
                hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)

                # only those edges which cross both boundaries
                # i.e. those which are not whithin the domain
                bx_tmp = bx_tmp_prime
                by_tmp = by_tmp_prime
                dx_tmp = dx_tmp_prime
                dy_tmp = dy_tmp_prime
                # crossing both boundaries
                mask = np.invert(np.isnan(dx_tmp.where(abs(dx_tmp) >  0.5 * Lx)\
                                                .where(abs(dy_tmp) >  0.5 * Ly))
                                )
                mask = (abs(dx_tmp) >  0.5 * Lx)
                mask = (mask.where(abs(dy_tmp) >  0.5 * Ly) == 1)

                bx_tmp = bx_tmp.where(mask)
                by_tmp = by_tmp.where(mask)
                dx_tmp = dx_tmp.where(mask)
                dy_tmp = dy_tmp.where(mask)
                
                mask = (dx_tmp >= 0.5 * Lx)
                dx_tmp += mask * (-Lx)
                mask = (dx_tmp <= -0.5 * Lx)
                dx_tmp += mask * Lx
                mask = (dy_tmp >= 0.5 * Ly)
                dy_tmp += mask * (-Ly)
                mask = (dy_tmp <= -0.5 * Ly)
                dy_tmp += mask * (+Ly)
                
                mask_1 = (bx_tmp < 0.5 * Lx)
                mask_2 = (bx_tmp > 0.5 * Lx)
                mask_3 = (by_tmp < 0.5 * Ly)
                mask_4 = (by_tmp > 0.5 * Ly)

                quiver_args = [bx_tmp + mask_1 * Lx - mask_2 * Lx,
                               by_tmp,
                               dx_tmp,
                               dy_tmp]
                if edge_property:
                    quiver_args.append(e_prop_data)
                hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)
                
                quiver_args = [bx_tmp,
                               by_tmp + mask_3 * Ly - mask_4 * Ly,
                               dx_tmp,
                               dy_tmp]
                if edge_property:
                    quiver_args.append(e_prop_data)
                hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)


            ### plot cells
            cell_type = c_data.sel(property="cell_type")
            x = c_data.sel(property="x")
            y = c_data.sel(property="y")
            color = ['red' if d == 1 else 'grey' for d in cell_type]
            hlpr.ax.scatter(x, y, c=color, s=cell_marker_size,
                            alpha=0.5)
            
            # gather the data for property interpolation
            # for the energies of the cells
            if (   property == 'area_elasticity'
                or property == 'cell_contractility'):
                prop_data = grp['Cell_energies'][time]
                prop_data = prop_data.sel(energy_term=property)
                prop_data = prop_data.assign_coords({'x': x, 'y': y})
            
            # for the energies of the edges
            elif (   property == 'linetension'
                or property == 'edge_contractility'):
                prop_data = grp['Edge_energies'][time]
                prop_data = prop_data.sel(energy_term=property).squeeze()
                prop_data = prop_data.assign_coords({'x': (ax + dx / 2.),
                                                     'y': (ay + dy / 2.)})


            # for other cell data
            elif property:
                prop_data = c_data.sel(property=property)
                prop_data = prop_data.assign_coords({'x': x, 'y': y})

            # perform the interpolation
            if property:
                if abs(prop_data.min().data - prop_data.max().data) < 1.e-12:
                    prop_data.data[0] += 1.e-10
                if property_hair_cells_only:
                    prop_data = prop_data.where(cell_type == 1)
                if property_support_cells_only:
                    prop_data = prop_data.where(cell_type == 2)
                if property_bulk_cells_only:
                    is_boundary = c_data.sel(property="is_boundary").round()
                    prop_data = prop_data.where(is_boundary == 0)
                
                min_x = floor(domain_size_min_x)
                max_x = ceil(domain_size_max_x)
                min_y = floor(domain_size_min_y)
                max_y = ceil(domain_size_max_y)
                grid_x, grid_y = np.mgrid[min_x:max_x:512j, min_y:max_y:512j]
                grid_z1 = griddata((prop_data.x, prop_data.y), prop_data,
                                    (grid_x, grid_y),
                                    **property_interpolation_kwargs)

                interpol = hlpr.ax.imshow(grid_z1.T,
                                          extent=(min_x,
                                                  max_x,
                                                  min_y,
                                                  max_y),
                                          origin='lower',
                                          **property_interpolation_plot_kwargs)
                
                cbar = hlpr.fig.colorbar(interpol, ax=hlpr.ax, extend='both')
                cbar.set_label(label=property)
                cbar.minorticks_on()

            if polarity_HC:
                dx = np.cos(c_data.sel(property="polarity")).where(cell_type == 1)
                dy = np.sin(c_data.sel(property="polarity")).where(cell_type == 1)
                hlpr.ax.quiver(x, y, dx, dy)


            hlpr.invoke_helper('set_title', title="Time {}".format(time))
            hlpr.invoke_helper('set_labels',
                               x=r"$x \ [A_0^{1/2}]$", y=r"$y \ [A_0^{1/2}]$")

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
                hlpr.invoke_helper('set_limits',
                                   x=(domain_size_min_x, domain_size_max_x),
                                   y=(domain_size_min_y, domain_size_max_y))
            
            hlpr.ax.set_aspect('equal')

            # end update here
            yield

    hlpr.register_animation_update(update)
