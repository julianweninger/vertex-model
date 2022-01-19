"""SavannaHeterogeneous-model specific plot function for the state / density"""

import logging
from typing import Tuple, Union

from math import floor, ceil
import numpy as np
import math as m
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
                       vector_property_is_nematic: bool=False,
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
            curvature = e_data.attrs["curvature"][0]
            max_theta = Lx / 2. * curvature


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

            def displacement(ax, ay, bx, by):
                if not vertex_cfg['space']['periodic']:
                    return bx - ax, by - ay
                
                if abs(curvature) < 1.e-12:
                    dx = bx - ax
                    dy = by - ay

                    # skew
                    dy -= np.round(dx / Lx) * skew_y
                    dx -= np.round(dy / Ly) * skew_x

                    # periodicity
                    dx -= np.round(dx / Lx) * Lx
                    dy -= np.round(dy / Ly) * Ly

                    return dx, dy

                R = 1. / curvature
                _ax = ax - Lx / 2.
                _ay = ay - Ly / 2. + R
                _bx = bx - Lx / 2.
                _by = by - Ly / 2. + R
                
                a_rho = (_ax**2 + _ay**2)**0.5
                a_theta = np.arctan2(_ax, _ay)
                b_rho = (_bx**2 + _by**2)**0.5
                b_theta = np.arctan2(_bx, _by)
                
                # skew
                skew_theta = skew_x * curvature
                b_rho -= np.round((b_theta - a_theta) / 2 / max_theta) * skew_y
                b_theta -= np.round((b_rho - a_rho) / Ly) * skew_theta

                # curvature
                b_theta -= np.round((b_theta-a_theta)/2./max_theta)*2*max_theta
                b_rho -= np.round((b_rho - a_rho) / Ly) * Ly

                _ax += Lx / 2.
                _ay += Ly / 2. - R
                _bx = b_rho * np.sin(b_theta) + Lx / 2.
                _by = b_rho * np.cos(b_theta) + Ly / 2. - R

                return _bx - ax, _by - ay

            dx, dy = displacement(ax, ay, bx, by)

            def quiver_and_colors(x, y, dx, dy, *, colorbar=True):
                quiver_args = [x, y, dx, dy]
                _quiver_kwargs = dict(headlength=0., headaxislength=0.,
                                    headwidth=0., scale=1, scale_units='xy',
                                    color='black')
                if quiver_kwargs:
                    _quiver_kwargs.update(quiver_kwargs)
            
                # for the color of the edges
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
                        {'x': (x + dx / 2.),
                         'y': (y + dy / 2.)})
                    if edge_property is not None and edge_property_split is None:
                        # append coloring
                        quiver_args.append(e_prop_data)

                    else:
                        __quiver_kwargs = dict(cmap='seismic')
                        __quiver_kwargs.update(_quiver_kwargs)


                        if len(edge_property_split) % 2 != 0:
                            raise RuntimeError("Edge split property not mod 2")
                        
                        N = int(len(edge_property_split) / 2)

                        length = (dx**2 + dy**2)**0.5
                        shift_x = 0.03 * -dy / length
                        shift_y = 0.03 *  dx / length


                        for i in range(0, N):
                            _x = x + dx / N * i
                            _y = y + dy / N * i
                            hlpr.ax.quiver(
                                _x+shift_x, _y+shift_y, dx / N, dy / N, 
                                e_prop_data.sel(**edge_property_split[i]),
                                **__quiver_kwargs)

                            quiver = hlpr.ax.quiver(
                                _x-shift_x, _y-shift_y, dx / N, dy / N, 
                                e_prop_data.sel(**edge_property_split[N + i]),
                                **__quiver_kwargs)
                        if colorbar:
                            cbar = hlpr.fig.colorbar(quiver, ax=hlpr.ax,
                                                     extend='both')
                            cbar.set_label(label=edge_property)
                            cbar.minorticks_on()

                quiver = hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)

                if len(quiver_args) == 5 and colorbar:
                    cbar = hlpr.fig.colorbar(quiver, ax=hlpr.ax, extend='both')
                    cbar.set_label(label=edge_property)
                    cbar.minorticks_on()

            quiver_and_colors(ax, ay, dx, dy)

            ### plot duplicates of periodic edges
            if vertex_cfg['space']['periodic']:
                mask = ((abs(bx - ax) > Lx / 2.) | (abs(by - ay) > Ly / 2.))
                # NOTE bitwise or

                _bx = np.where(mask, bx, np.nan)
                _by = np.where(mask, by, np.nan)

                _dx, _dy = displacement(bx, by, ax, ay)

                # use the inverse arrows
                quiver_and_colors((_bx + _dx), (_by + _dy), -_dx, -_dy,
                                  colorbar=False)


                ## plot edges that cross 2 boundaries, i.e. corners
                mask =  ((abs(bx - ax) > Lx / 2.) & (abs(by - ay) > Ly / 2.))
                # NOTE bitwise and
                _ax = ax.where(mask)
                _ay = ay.where(mask)
                _bx = bx.where(mask)
                _by = by.where(mask)

                # rotate corners
                if abs(curvature) < 1.e-12:
                    _ay -= (2. * np.round(_ay / Ly) - 1.) * Ly
                    _ax += (2. * np.round(_ay / Ly) - 1.) * skew_x
                    _by -= (2. * np.round(_by / Ly) - 1.) * Ly
                    _bx += (2. * np.round(_by / Ly) - 1.) * skew_x
                else:
                    R = 1. / curvature
                    _ax = _ax - Lx / 2.
                    _ay = _ay - Ly / 2. + R
                    _bx = _bx - Lx / 2.
                    _by = _by - Ly / 2. + R
                    
                    a_rho = (_ax**2 + _ay**2)**0.5
                    a_theta = np.arctan2(_ay, _ax)
                    b_rho = (_bx**2 + _by**2)**0.5
                    b_theta = np.arctan2(_by, _bx)
                    
                    # skew
                    skew_theta = skew_x * curvature
                    a_theta -= (a_rho - R) / abs(a_rho - R) * skew_theta
                    b_theta -= (b_rho - R) / abs(b_rho - R) * skew_theta

                    # rotate corners
                    a_rho -= (a_rho - R) / abs(a_rho - R) * Ly
                    b_rho -= (b_rho - R) / abs(b_rho - R) * Ly

                    _ax *= 0.
                    _ax += a_rho * np.sin(a_theta) + Lx / 2.
                    _ay *= 0.
                    _ay += a_rho * np.cos(a_theta) + Ly / 2. - R
                    _bx *= 0.
                    _bx += b_rho * np.sin(b_theta) + Lx / 2.
                    _by *= 0.
                    _by += b_rho * np.cos(b_theta) + Ly / 2. - R
                
                
                quiver_and_colors(_ax, _ay, dx, dy,
                                  colorbar=False)
                
                __dx, __dy = displacement(_bx, _by, _ax, _ay)
                quiver_and_colors((_bx + __dx), (_by + __dy), -__dx, -__dy,
                                  colorbar=False)

            ### plot cells
            cell_type = c_data.sel(property="cell_type")
            x = c_data.sel(property="x")
            y = c_data.sel(property="y")
            color = ['red' if d == 1 else 'grey' for d in cell_type]
            hlpr.ax.scatter(x, y, c=color, s=cell_marker_size,
                            alpha=0.5)
            
            # plot cell data
            if property is not None or property_path is not None:
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
                
                # Assign property to every cell
                if not 'id' in prop_data.dims:
                    if 'x' in prop_data.coords and 'y' in prop_data.coords:
                        log.warning("Could not map property data to cells, "
                                    "but found coordinates. Using coordinates "
                                    "provided in property data.")
                    else:
                        raise ValueError("Expected `id` in property_data dims "
                            "(were {}).", prop_data.dims)
                else:
                    prop_data = prop_data.assign_coords({'x': x, 'y': y})
                
                if property is not None:
                    prop_data = prop_data.sel(**property)

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
                grid_z1 = griddata((prop_data.x, prop_data.y),
                                   prop_data.squeeze(),
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

                _v_property_kwargs = dict(angles='xy', scale_units='xy',
                                          scale=3.)
                if vector_property_kwargs:
                    _v_property_kwargs.update(vector_property_kwargs)

                hlpr.ax.quiver(x, y, prop_x_data, prop_y_data,
                               **_v_property_kwargs)
                if vector_property_is_nematic:
                    hlpr.ax.quiver(x, y, -prop_x_data, -prop_y_data,
                                **_v_property_kwargs)




            hlpr.invoke_helper('set_title', title="Time {}".format(time))
            hlpr.invoke_helper('set_labels',
                               x=r"$x \ [A_0^{1/2}]$", y=r"$y \ [A_0^{1/2}]$")

            if (vertex_cfg['space']['periodic']):
                if abs(curvature) > 1.e-12:        
                    R = 1. / curvature
                    Ox = Lx / 2.
                    Oy = Ly / 2. - R
                    
                    if max_theta < m.pi/2. - 0.2:
                        hlpr.invoke_helper('set_limits', 
                            x=((R + Ly/2.) * np.sin(-max_theta)+Lx/2,
                            (R + Ly/2.) * np.sin( max_theta)+Lx/2),
                            y=((R - Ly/2.) * np.cos(-max_theta) - R + Ly/2,
                                Oy + R + Ly/2.))
                    else:
                        hlpr.invoke_helper('set_limits', 
                            x=(-(R + Ly/2.) + Lx/2, (R + Ly/2.) + Lx/2),
                            y=( (R + Ly/2.) * np.cos(-max_theta) - R + Ly/2,
                                Oy + R + Ly/2.))
                    
                    if max_theta < m.pi-0.1:
                        hlpr.ax.quiver(Ox, Oy,
                                    (R + Ly / 2.) * np.sin( max_theta),
                                    (R + Ly / 2.) * np.cos( max_theta),
                                    headlength=0., headaxislength=0.,
                                    headwidth=0., scale=1, scale_units='xy',
                                    color='black')
                        hlpr.ax.quiver(Ox, Oy,
                                    (R + Ly / 2.) * np.sin(-max_theta),
                                    (R + Ly / 2.) * np.cos(-max_theta),
                                    headlength=0., headaxislength=0.,
                                    headwidth=0., scale=1, scale_units='xy',
                                    color='black')
                    circle_inner = plt.Circle((Ox, Oy), R - Ly / 2.,
                                              edgecolor='black', fill=False)
                    circle_outer = plt.Circle((Ox, Oy), R + Ly / 2.,
                                              edgecolor='black', fill=False)

                    hlpr.ax.add_patch(circle_inner)
                    hlpr.ax.add_patch(circle_outer)
                else:                    
                    hlpr.invoke_helper('set_limits', x=(0, Lx), y=(0, Ly))



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
