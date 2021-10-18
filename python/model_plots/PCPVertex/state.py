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
                       property: dict=None,
                       property_path: str=None,
                       property_hair_cells_only: bool=False,
                       property_support_cells_only: bool=False,
                       property_bulk_cells_only: bool=False,
                       property_interpolation_kwargs: dict={},
                       property_interpolation_plot_kwargs: dict={},
                       property_ignore_nan: bool=False,
                       property_resolution: Tuple[int, int]=[512,512],
                       edge_property_path: str=None,
                       edge_property: dict=None,
                       edge_property_split: Tuple[dict, dict]=None,
                       quiver_kwargs: dict=None,
                       vector_property: dict=None,
                       vector_property_path: str=None,
                       vector_property_kwargs: dict=None,
                       vector_property_hair_cells_only: bool=False,
                       vector_property_support_cells_only: bool=False,):
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
        property_grp (str, optional)
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
        property_ignore_nan (bool, default: False): Whether to generate the
            interpolation ignoring nan values, e.g. when selecting only 
            cell type. Nan values will interrupt the interpolation, i.e. 
            the location will remain empty.
        property_resolution (tuple of int, default: [512, 512]): The resolution
            in horizontal and vertical direction of the interpolation meshwork.
        edge_property (std, optional): An additional edge property to plot.
            Data is edge data where property=edge_property or
            Edge_energy where energy_term=edge_property if edge_property
            is one of the edge energies.
        quiver_kwargs (dict, optional): Updates the quiver kwargs used.
        vector_property (dict): Dict with x and y entries for direction of 
            a vector. Otherwise as property. vector_property_kwargs 
            forwarded to mpl.quiver
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
            
            # periodic skrewed boundary condition
            skew_x = e_data.attrs["skew_x"][0]
            skew_y = e_data.attrs["skew_y"][0]


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

            ## map edges crossing periodic boundary
            if (vertex_cfg['space']['periodic']):
                # left boundary
                mask = (dx >= 0.5 * Lx)
                dx += mask * (-Lx)
                dy -= mask * skew_y
                # right boundary
                mask = (dx <= -0.5 * Lx)
                dx += mask * Lx
                dy += mask * skew_y
                # lower boundary
                mask = (dy >= 0.5 * Ly)
                dy += mask * (-Ly)
                dx -= mask * skew_x
                # upper boundary
                mask = (dy <= -0.5 * Ly)
                dy += mask * (+Ly)
                dx += mask * skew_x
            
            quiver_args = [ax, ay, dx, dy]
            _quiver_kwargs = dict(headlength=0., headaxislength=0.,
                                  headwidth=0., scale=1, scale_units='xy',
                                  color='black')
            if quiver_kwargs:
                _quiver_kwargs.update(quiver_kwargs)
            
            # for the energies of the edges
            if edge_property is not None:
                if edge_property_path is not None:
                    if not edge_property_path in grp:
                        raise ValueError("Failed to access property of edges at path "
                            "'data/{}/{}' (relative to '{}'. Available paths: {}"
                            "".format(datapath, edge_property_path, datapath, grp.keys())
                        )
                    if not time in grp[edge_property_path]:
                        raise ValueError('No edge property data available at time {}.'
                                        ''.format(time))
                    e_prop_data = grp[edge_property_path][time]
                else:
                    e_prop_data = e_data

                e_prop_data = e_prop_data.sel(**edge_property).squeeze()
                e_prop_data = e_prop_data.assign_coords(
                    {'x': (ax + dx / 2.).drop('property'),
                     'y': (ay + dy / 2.).drop('property')})   
                if edge_property is not None and edge_property_split is None:
                    # append coloring
                    quiver_args.append(e_prop_data)

                else:
                    __quiver_kwargs = dict(cmap='seismic')
                    __quiver_kwargs.update(_quiver_kwargs)

                    length = (dx**2 + dy**2)**0.5
                    shift_x = 0.03 * -dy / length
                    shift_y = 0.03 *  dx / length
                    hlpr.ax.quiver(ax+shift_x, ay+shift_y, dx, dy, 
                                   e_prop_data.sel(**edge_property_split[0]),
                                   **__quiver_kwargs)

                    quiver = hlpr.ax.quiver(ax-shift_x, ay-shift_y, dx, dy,
                                    e_prop_data.sel(**edge_property_split[1]),
                                    **__quiver_kwargs)

                    cbar = hlpr.fig.colorbar(quiver, ax=hlpr.ax, extend='both')
                    cbar.set_label(label=edge_property)
                    cbar.minorticks_on()

            quiver = hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)

            if (len(quiver_args) == 5):
                cbar = hlpr.fig.colorbar(quiver, ax=hlpr.ax, extend='both')
                cbar.set_label(label=edge_property)
                cbar.minorticks_on()

            ### plot duplicates of periodic edges
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
                ax_tmp_prime = ax_tmp
                ay_tmp_prime = ay_tmp
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
                dy_tmp -= mask * skew_y
                mask = (dx_tmp <= -0.5 * Lx)
                dx_tmp += mask * Lx
                dy_tmp += mask * skew_y
                mask = (dy_tmp >= 0.5 * Ly)
                dy_tmp += mask * (-Ly)
                dx_tmp -= mask * skew_x
                mask = (dy_tmp <= -0.5 * Ly)
                dy_tmp += mask * (+Ly)
                dx_tmp += mask * skew_x
                
                # ax = bx + (ax - bx), the inverse arrow
                quiver_args = [bx_tmp+dx_tmp, by_tmp+dy_tmp, -dx_tmp, -dy_tmp]
                if edge_property is not None and edge_property_split is None:
                    quiver_args.append(e_prop_data)
                elif edge_property is not None:
                    length = (dx_tmp**2 + dy_tmp**2)**0.5
                    shift_x_tmp = 0.03 *  dy_tmp / length
                    shift_y_tmp = 0.03 * -dx_tmp / length
                    hlpr.ax.quiver(bx_tmp + shift_x_tmp + dx_tmp,
                                   by_tmp + shift_y_tmp + dy_tmp,
                                   -dx_tmp, 
                                   -dy_tmp, 
                                   e_prop_data.sel(**edge_property_split[0]),
                                   **__quiver_kwargs)
                    hlpr.ax.quiver(bx_tmp - shift_x_tmp + dx_tmp,
                                   by_tmp - shift_y_tmp + dy_tmp,
                                   -dx_tmp, 
                                   -dy_tmp, 
                                   e_prop_data.sel(**edge_property_split[1]),
                                   **__quiver_kwargs)
                hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)


                ## plot edges that cross 2 boundaries, i.e. corners
                ## first the a->b arrow 
                ax_tmp = ax_tmp_prime
                ay_tmp = ay_tmp_prime
                bx_tmp = bx_tmp_prime
                by_tmp = by_tmp_prime
                dx_tmp = -dx_tmp_prime
                dy_tmp = -dy_tmp_prime

                # crossing both boundaries
                mask = np.invert(np.isnan(dx_tmp.where(abs(dx_tmp) >  0.5 * Lx)\
                                                .where(abs(dy_tmp) >  0.5 * Ly))
                                )
                mask = (abs(dx_tmp) >  0.5 * Lx)
                mask = (mask.where(abs(dy_tmp) >  0.5 * Ly) == 1)

                ax_tmp = ax_tmp.where(mask)
                ay_tmp = ay_tmp.where(mask)
                bx_tmp = bx_tmp.where(mask)
                by_tmp = by_tmp.where(mask)
                dx_tmp = dx_tmp.where(mask)
                dy_tmp = dy_tmp.where(mask)
                
                # left boundary
                mask = (dx_tmp >= 0.5 * Lx)
                dx_tmp += mask * (-Lx)
                dy_tmp -= mask * skew_y
                # left boundary
                mask = (dx_tmp <= -0.5 * Lx)
                dx_tmp += mask * Lx
                dy_tmp += mask * skew_y
                # left boundary
                mask = (dy_tmp >= 0.5 * Ly)
                dy_tmp += mask * (-Ly)
                dx_tmp -= mask * skew_x
                # left boundary
                mask = (dy_tmp <= -0.5 * Ly)
                dy_tmp += mask * (+Ly)
                dx_tmp += mask * skew_x

                ## copy vertices to other side of domain
                # vertices close to left boundary
                ax_tmp += Lx * (ax_tmp_prime < 0.5 * Lx)
                bx_tmp += Lx * (bx_tmp_prime < 0.5 * Lx)
                # vertices close to right boundary
                ax_tmp -= Lx * (ax_tmp_prime > 0.5 * Lx)
                bx_tmp -= Lx * (bx_tmp_prime > 0.5 * Lx)
                # NOTE diagonal copy already has been made
                #      doing horizontal shift, equivalent to vertical shift

                quiver_args = [ax_tmp, ay_tmp, dx_tmp, dy_tmp]
                if edge_property is not None and edge_property_split is None:
                    quiver_args.append(e_prop_data)
                elif edge_property is not None:
                    length = (dx_tmp**2 + dy_tmp**2)**0.5
                    shift_x_tmp = 0.03 *  dy_tmp / length
                    shift_y_tmp = 0.03 * -dx_tmp / length
                    hlpr.ax.quiver(ax_tmp + shift_x_tmp + dx_tmp,
                                   ay_tmp + shift_y_tmp + dy_tmp,
                                   dx_tmp, 
                                   dy_tmp, 
                                   e_prop_data.sel(**edge_property_split[0]),
                                   **__quiver_kwargs)
                    hlpr.ax.quiver(ax_tmp - shift_x_tmp + dx_tmp,
                                   ay_tmp - shift_y_tmp + dy_tmp,
                                   dx_tmp, 
                                   dy_tmp, 
                                   e_prop_data.sel(**edge_property_split[1]),
                                   **__quiver_kwargs)

                hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)
                
                ## the inverse arrow
                quiver_args = [bx_tmp - dx_tmp, by_tmp - dy_tmp, dx_tmp, dy_tmp]
                if edge_property is not None and edge_property_split is None:
                    quiver_args.append(e_prop_data)
                elif edge_property is not None:
                    length = (dx_tmp**2 + dy_tmp**2)**0.5
                    shift_x_tmp = 0.03 *  dy_tmp / length
                    shift_y_tmp = 0.03 * -dx_tmp / length
                    hlpr.ax.quiver(bx_tmp + shift_x_tmp - dx_tmp,
                                   by_tmp + shift_y_tmp - dy_tmp,
                                   -dx_tmp, 
                                   -dy_tmp, 
                                   e_prop_data.sel(**edge_property_split[0]),
                                   **__quiver_kwargs)
                    hlpr.ax.quiver(bx_tmp - shift_x_tmp - dx_tmp,
                                   by_tmp - shift_y_tmp - dy_tmp,
                                   -dx_tmp, 
                                   -dy_tmp, 
                                   e_prop_data.sel(**edge_property_split[1]),
                                   **__quiver_kwargs)
                                   
                hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)


            ### plot cells
            cell_type = c_data.sel(property="cell_type")
            x = c_data.sel(property="x")
            y = c_data.sel(property="y")
            color = ['red' if d == 1 else 'grey' for d in cell_type]
            hlpr.ax.scatter(x, y, c=color, s=cell_marker_size,
                            alpha=0.5)
            
            # plot cell data
            if property is not None:
                if property_path is not None:
                    if not property_path in grp:
                        raise ValueError("Failed to access property data at path "
                            "'data/{}/{}' (relative to '{}'. Available paths: {}"
                            "".format(datapath, property_path, datapath, grp.keys())
                        )
                    if not time in grp[property_path]:
                        raise ValueError('No property data available at time {}.'
                                        ''.format(time))
                    prop_data = grp[property_path][time]
                else:
                    prop_data = c_data
                prop_data = prop_data.sel(**property)
                prop_data = prop_data.assign_coords({'x': x, 'y': y})

                # perform the interpolation
                if abs(prop_data.min().data - prop_data.max().data) < 1.e-12:
                    prop_data.data[0] += 1.e-10
                if property_hair_cells_only:
                    prop_data = prop_data.where(cell_type == 1)
                if property_support_cells_only:
                    prop_data = prop_data.where(cell_type == 2)
                if property_bulk_cells_only:
                    is_boundary = c_data.sel(property="is_boundary").round()
                    prop_data = prop_data.where(is_boundary == 0)
                if property_ignore_nan:
                    prop_data = prop_data.dropna(dim='id')

                min_x = floor(domain_size_min_x)
                max_x = ceil(domain_size_max_x)
                min_y = floor(domain_size_min_y)
                max_y = ceil(domain_size_max_y)
                
                Nx = property_resolution[0]
                Ny = property_resolution[1]
                
                grid_x, grid_y = np.mgrid[min_x:max_x:Nx*1j, min_y:max_y:Ny*1j]
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
            
            if vector_property is not None:
                if vector_property_path is not None:
                    if not vector_property_path in grp:
                        raise ValueError("Failed to access property data at path "
                            "'data/{}/{}' (relative to '{}'. Available paths: {}"
                            "".format(datapath, vector_property_path, datapath, grp.keys())
                        )
                    if not time in grp[vector_property_path]:
                        raise ValueError('No vectoriel property data available at time {}.'
                                        ''.format(time))
                    prop_v_data = grp[vector_property_path][time]
                else:
                    prop_v_data = c_data

                if len(vector_property) == 2:
                    prop_x_data = prop_v_data.sel(vector_property['x'])
                    prop_y_data = prop_v_data.sel(vector_property['y'])
                elif len(vector_property) == 1:
                    prop_x_data = np.cos(prop_v_data.sel(vector_property))
                    prop_y_data = np.sin(prop_v_data.sel(vector_property))
                else:
                    raise ValueError("Expected dict with 2 entries (x and y) "
                                     "or 1 entry (angle), but dict {} has {} "
                                     "entries!".format(vector_property,
                                     len(vector_property)))

                if vector_property_hair_cells_only:
                    prop_x_data = prop_x_data.where(cell_type == 1)
                    prop_y_data = prop_y_data.where(cell_type == 1)
                if vector_property_support_cells_only:
                    prop_x_data = prop_x_data.where(cell_type == 2)
                    prop_y_data = prop_y_data.where(cell_type == 2)

                _v_property_kwargs = dict(angles='xy', scale_units='xy', scale=3.)
                if vector_property_kwargs:
                    _v_property_kwargs.update(vector_property_kwargs)

                hlpr.ax.quiver(x, y, prop_x_data, prop_y_data,
                               **_v_property_kwargs)



            hlpr.invoke_helper('set_title', title="Time {}".format(time))
            hlpr.invoke_helper('set_labels',
                               x=r"$x \ [A_0^{1/2}]$", y=r"$y \ [A_0^{1/2}]$")

            if (vertex_cfg['space']['periodic']):
                hlpr.invoke_helper('set_limits', x=(-.2,Lx+.2), y=(-.2,Ly+.2))
                hlpr.ax.axvline(x=0, c='gray', linestyle=':')
                hlpr.ax.axvline(x=Lx, c='gray', linestyle=':')
                hlpr.ax.axhline(y=0, c='gray', linestyle=':')
                hlpr.ax.axhline(y=Ly, c='gray', linestyle=':')

                if skew_x > 1.e-12:
                    hlpr.ax.axvline(x=skew_x, ymin=1. - 0.5 / Ly, c='gray',
                                    linestyle=':')
                elif skew_x < 1.e-12:
                    hlpr.ax.axvline(x=Lx-skew_x, ymax = 0.5 / Ly, c='gray',
                                    linestyle=':')
                if skew_y > 1.e-12:                    
                    hlpr.ax.axhline(y=skew_y, xmin=1. - 0.5 / Lx, c='gray',
                                    linestyle=':')
                elif skew_y < 1.e-12:                    
                    hlpr.ax.axhline(y=skew_y, xmax = 0.5 / Lx, c='gray', 
                                    linestyle=':')


            else:
                hlpr.invoke_helper('set_limits',
                                   x=(domain_size_min_x, domain_size_max_x),
                                   y=(domain_size_min_y, domain_size_max_y))
            
            hlpr.ax.set_aspect('equal')

            # end update here
            yield

    hlpr.register_animation_update(update)
