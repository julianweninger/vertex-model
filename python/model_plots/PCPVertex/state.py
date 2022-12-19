"""PCPVertex-model specific plot function for the state"""

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
from utopya.tools import recursive_update

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
    energy = uni['data'][model_name]['Energy']['Terms']
    
    transitions = uni['data'][model_name]['Energy/Transitions']
    num_T1s = transitions.sel(type="num_T1s")
    num_T1s_attempted = transitions.sel(type="num_T1s_attempted")
    num_T2s = transitions.sel(type="num_T2s")
    if cumsum_transitions:
        num_T1s = num_T1s.cumsum()
        num_T1s_attempted = num_T1s_attempted.cumsum()
        num_T2s = num_T2s.cumsum()


    # Create the line plot of energy
    ax1 = hlpr.ax
    ax1.plot(
        energy.time,
        energy,
        hue='terms',
        alpha=0.5,
    )

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


@is_plot_func(
    creator_type=UniversePlotCreator,
    use_dag=True,
    required_dag_tags=(
        'Vertices',
        'Edges',
        'Cells',
        'periodic',
        'Lx',
        'Ly'
    ),
    compute_only_required_dag_tags=False,
    supports_animation=True
)
def cellular_structure(
    *,
    data: dict,
    hlpr: PlotHelper,
    select_times: list=None,
    cell_center_marker_kwargs: dict=None,
    plot_vertices: bool=False,
    plot_excess_length: float=1.,
    property_interpolation_kwargs: dict=None,
    property_interpolation_plot_kwargs: dict=None,
    property_resolution: Tuple[int, int]=[512,512],
    quiver_kwargs: dict=None,
    vector_property_kwargs: dict=None,
    vector_property_is_nematic: bool=False,
    set_limits: dict=None,
):
    """Performs a plot of the cells, edges and vertices
    
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (int): The universe to use
        hlpr (PlotHelper): The PlotHelper
        datapath (str): Path to the cell manager's data
        cfgpath (str): Path to the vertex model's configuration
        select_times (list, optional): Timepoints to select for plotting.
            If None, all timepoints are plotted
        cell_center_marker_kwargs (dict): Forwarded to ax.scatter for cell
            centers. Can have additional colors argument that will be used to
            define colors of cells
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

    _cell_center_marker_kwargs = dict({
        's': 40,
        'linewidth': 0,
        'alpha': 0.5,
        'colors': ['white', 'r', 'gray']
    })
    if cell_center_marker_kwargs is not None:
        _cell_center_marker_kwargs = recursive_update(
            _cell_center_marker_kwargs,
            cell_center_marker_kwargs
        )
    cell_marker_colors = _cell_center_marker_kwargs.pop('colors')

    _property_interpolation_kwargs=dict({})
    if property_interpolation_kwargs is not None:
        _property_interpolation_kwargs = recursive_update(
            _property_interpolation_kwargs,
            property_interpolation_kwargs
        )
    _property_interpolation_plot_kwargs=dict({})
    if property_interpolation_plot_kwargs is not None:
        _property_interpolation_plot_kwargs = recursive_update(
            _property_interpolation_plot_kwargs,
            property_interpolation_plot_kwargs
        )

    # Prepare the figure ......................................................
    # Prepare the figure to have as many columns as there are properties
    hlpr.setup_figure()

    
    # Dynamically provide some information to the plot helper
    hlpr.provide_defaults('set_title',
                            title="My data at time {}".format(0))
    hlpr.provide_defaults('set_labels', y=dict(label="My data"))
    hlpr.provide_defaults('set_limits', y=[0, 1])

    def update():
        global cbar
        cbar = None
        global edges_cbar
        edges_cbar = None

        # the domain extent for non periodic boundaries (centered on (0., 0.))
        domain_size_min_x = 0.
        domain_size_max_x = 0.
        domain_size_min_y = 0.
        domain_size_max_y = 0.

        if (not data['periodic']):
            domain_size_min_x =  1000000.
            domain_size_max_x = -1000000.
            domain_size_min_y =  1000000.
            domain_size_max_y = -1000000.
            domain_size_min_x = min(
                domain_size_min_x,
                data['Vertices'].sel(property="x").min())
            domain_size_max_x = max(
                domain_size_max_x,
                data['Vertices'].sel(property="x").max())
            domain_size_min_y = min(
                domain_size_min_y,
                data['Vertices'].sel(property="y").min())
            domain_size_max_y = max(
                domain_size_max_y,
                data['Vertices'].sel(property="y").max())

        if select_times is not None:
            times = select_times
        else:
            times = np.unique(data['Vertices'].coords['time'])

        for time in times:
            hlpr.ax.clear()

            if cbar is not None: 
                cbar.remove()

            if edges_cbar is not None: 
                edges_cbar.remove()
            hlpr.ax.set_aspect('auto')

            v_data = data['Vertices'].sel(time=time)
            e_data = data['Edges'].sel(time=time)
            c_data = data['Cells'].sel(time=time)

            Lx = data['Lx'].sel(time=time).data
            Ly = data['Ly'].sel(time=time).data
            if data['periodic']:
                domain_size_max_x = Lx
                domain_size_max_y = Ly
            
            # periodic skewed boundary condition
            if "skew_x" in e_data.attrs:
                skew_x = e_data.attrs["skew_x"][0]
            else:
                skew_x = 0.
            if "skew_y" in e_data.attrs:
                skew_y = e_data.attrs["skew_y"][0]
            else:
                skew_y = 0.
            if "curvature" in e_data.attrs:
                curvature = e_data.attrs["curvature"][0]
            else: 
                curvature = 0.
            max_theta = Lx / 2. * curvature


            ### plot vertices
            if plot_vertices:
                hlpr.ax.scatter(v_data.sel(property="x"),
                                v_data.sel(property="y"),
                                c="black")
                if data['periodic']:
                    ax = v_data.sel(property="x")
                    ay = v_data.sel(property="y")
                    l = plot_excess_length
                    hlpr.ax.scatter(np.where(ax < l, ax, np.nan) + Lx,
                                    ay + skew_y,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ax > Lx - l, ax, np.nan) - Lx,
                                    ay - skew_y,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ay < l, ax, np.nan) + skew_x,
                                    ay + Ly,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ay > Ly - l, ax, np.nan) - skew_x,
                                    ay - Ly,
                                    c="grey")

                    hlpr.ax.scatter(np.where(ax < l, ax, np.nan) + Lx + skew_x,
                                    np.where(ay < l, ay, np.nan) + Ly + skew_y,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ax>Lx-l, ax, np.nan) - Lx + skew_x,
                                    np.where(ay < l, ay, np.nan) + Ly - skew_y,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ax>Lx-l, ax, np.nan) - Lx - skew_x,
                                    np.where(ay>Ly-l, ay, np.nan) - Ly - skew_y,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ax < l, ax, np.nan) + Lx - skew_x,
                                    np.where(ay>Ly-l, ay, np.nan) - Ly + skew_y,
                                    c="grey")



            ### plot edges
            # the ids of the vertices a and b of the edges
            vertex_a = e_data.sel(property="vertex_a").dropna(dim='id').astype(int)
            vertex_b = e_data.sel(property="vertex_b").dropna(dim='id').astype(int)

            # for id in vertex_a.squeeze().data:
            #     if not id in v_data.id:
            #         e_data = e_data.where(e_data!=id, drop=True)
            # for id in vertex_b.squeeze().data:
            #     if not id in v_data.id:
            #         e_data = e_data.where(e_data!=id, drop=True)
            #     # print(v_data.sel(id=vertex_a))
            # e_data = e_data.dropna(dim='id')

            # vertex_a = e_data.sel(property="vertex_a").dropna(dim='id').astype(int)
            # vertex_b = e_data.sel(property="vertex_b").dropna(dim='id').astype(int)

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
                if not data['periodic']:
                    return bx - ax, by - ay
                
                if abs(curvature) < 1.e-12:
                    dx = bx - ax
                    dy = by - ay

                    add_skew_x = - np.round(dy/Ly) * skew_x
                    add_skew_y = - np.round(dx/Lx) * skew_y

                    dx = bx + add_skew_x - ax
                    dy = by + add_skew_y - ay

                    return (dx - np.round(dx/Lx) * Lx,
                            dy - np.round(dy/Ly) * Ly)

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
                add_skew_rho = -np.round((b_theta - a_theta)/2/max_theta)*skew_y
                add_skew_theta = -np.round((b_rho - a_rho) / Ly) * skew_theta

                # curvature
                b_theta += add_skew_theta - np.round((b_theta-a_theta)/2./max_theta)*2*max_theta
                b_rho += add_skew_rho - np.round((b_rho - a_rho) / Ly) * Ly

                _ax += Lx / 2.
                _ay += Ly / 2. - R
                _bx = b_rho * np.sin(b_theta) + Lx / 2.
                _by = b_rho * np.cos(b_theta) + Ly / 2. - R

                return _bx - ax, _by - ay

            dx, dy = displacement(ax, ay, bx, by)

            def quiver_and_colors(x, y, dx, dy, *, colorbar=False):
                global edges_cbar

                if data['periodic']:
                    # only plot within box + excess length
                    x = xr.where(x > -2 * plot_excess_length, x, np.nan)
                    x = xr.where(x < Lx + 2 * plot_excess_length, x, np.nan)
                    y = xr.where(y > -2 * plot_excess_length, y, np.nan)
                    y = xr.where(y < Ly + 2 * plot_excess_length, y, np.nan)

                quiver_args = [x, y, dx, dy]
                _quiver_kwargs = dict(headlength=0., headaxislength=0.,
                                    headwidth=0., scale=1, scale_units='xy',
                                    color='black')
                if quiver_kwargs is not None:
                    _quiver_kwargs.update(quiver_kwargs)
            
                # for the color of the edges
                if 'edge_property' in data:
                    e_prop_data = data['edge_property'].sel(time=time, id=x.id)
                    e_prop_dims = [d for d in e_prop_data.dims if d != 'id']
                    # Assign property to every cell
                    if not 'x' in e_prop_data.coords or not 'y' in e_prop_data.coords:
                        if not 'id' in e_prop_data.dims:
                            raise ValueError("Expected `id` in property_data dims "
                                "(were {}).", e_prop_data.dims)
                    if len(e_prop_dims) == 0:
                        # append coloring
                        quiver_args.append(e_prop_data)

                        quiver = hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)

                        if colorbar:
                            edges_cbar = hlpr.fig.colorbar(
                                quiver, ax=hlpr.ax, extend='both',
                                fraction=0.046, pad=0.04
                            )
                            edges_cbar.set_label(e_prop_data.name)
                            edges_cbar.minorticks_on()

                    elif len(e_prop_dims) == 1:
                        e_prop_dim=e_prop_dims[0]
                        __quiver_kwargs = dict(cmap='seismic')
                        __quiver_kwargs.update(_quiver_kwargs)

                        e_prop_coords=e_prop_data.coords[e_prop_dim].data
                        if len(e_prop_coords) % 2 != 0:
                            raise RuntimeError("Edge split property not mod 2")
                        
                        N = int(len(e_prop_coords) / 2)

                        length = (dx**2 + dy**2)**0.5
                        shift_x = 0.03 * -dy / length
                        shift_y = 0.03 *  dx / length


                        for i in range(0, N):
                            _x = x + dx / N * i
                            _y = y + dy / N * i
                            name_0=e_prop_coords[i]
                            name_1=e_prop_coords[N+i]
                            hlpr.ax.quiver(
                                _x+shift_x, _y+shift_y, dx / N, dy / N, 
                                e_prop_data.sel({e_prop_dim: name_0}),
                                **__quiver_kwargs)

                            quiver = hlpr.ax.quiver(
                                _x-shift_x, _y-shift_y, dx / N, dy / N, 
                                e_prop_data.sel({e_prop_dim: name_1}),
                                **__quiver_kwargs)
                        
                        if colorbar:
                            edges_cbar = hlpr.fig.colorbar(quiver, ax=hlpr.ax,
                                                     extend='both')
                            edges_cbar.set_label(e_prop_data.name)
                            edges_cbar.minorticks_on()
                    else:
                        raise

                else:
                    quiver = hlpr.ax.quiver(*quiver_args, **_quiver_kwargs)


            quiver_and_colors(ax, ay, dx, dy, colorbar=True)

            ### plot duplicates of periodic edges
            if data['periodic'] and curvature < 1.e-12:
                length = max(skew_x, skew_y) + 1 + 1.e-6

                # plot periodic copies of the domain
                for i in range(1, int(np.ceil(length / Lx) + 1)):
                    quiver_and_colors(ax + skew_x, ay + Ly, dx, dy)
                    quiver_and_colors(ax - skew_x, ay - Ly, dx, dy)
                    quiver_and_colors(ax + i * Lx, ay + i * skew_y, dx, dy)
                    quiver_and_colors(ax - i * Lx, ay - i * skew_y, dx, dy)

                    quiver_and_colors(ax + i * Lx + skew_x, 
                                      ay + Ly + i * skew_y, dx, dy)
                    quiver_and_colors(ax - i * Lx + skew_x, 
                                      ay + Ly - i * skew_y, dx, dy)
                    quiver_and_colors(ax - i * Lx - skew_x, 
                                      ay - Ly - i * skew_y, dx, dy)
                    quiver_and_colors(ax + i * Lx - skew_x, 
                                      ay - Ly + i * skew_y, dx, dy)

            ### plot duplicates of periodic edges
            elif data['periodic']:
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
                    _ax += (2. * np.round(_ay / Ly) - 1.) * skew_x
                    _ay -= (2. * np.round(_ay / Ly) - 1.) * Ly
                    _bx += (2. * np.round(_by / Ly) - 1.) * skew_x
                    _by -= (2. * np.round(_by / Ly) - 1.) * Ly
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
            cell_type = c_data.sel(property="cell_type").dropna(dim='id').astype(int)
            x = c_data.sel(property="x").dropna(dim='id')
            y = c_data.sel(property="y").dropna(dim='id')
            hlpr.ax.scatter(
                x, y,
                c=[cell_marker_colors[d.data] for d in cell_type],
                **_cell_center_marker_kwargs
            )
            
            # plot cell data
            if 'property' in data:
                prop_data = data['property'].sel(time=time)

                # Assign property to every cell
                if not 'x' in prop_data.coords or not 'y' in prop_data.coords:
                    if not 'id' in prop_data.dims:
                        raise ValueError("Expected `id` in property_data dims "
                            "(were {}).", prop_data.dims)
                    else:
                        prop_data = prop_data.assign_coords({'x': x, 'y': y})

                # perform the interpolation
                if abs(prop_data.min().data - prop_data.max().data) < 1.e-12:
                    prop_data.data[0] += 1.e-10

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
                                   **_property_interpolation_kwargs)

                interpol = hlpr.ax.imshow(grid_z1.T,
                                          extent=(min_x,
                                                  max_x,
                                                  min_y,
                                                  max_y),
                                          origin='lower',
                                          **_property_interpolation_plot_kwargs)

                cbar = hlpr.fig.colorbar(interpol, ax=hlpr.ax, extend='both')

                cbar.set_label(label=prop_data.name)
                    
                cbar.minorticks_on()
            
            if 'vector_property' in data:
                prop_v_data = data['vector_property'].sel(time=time) #.dropna(dim='id')

                # Assign property to every cell
                if not 'x' in prop_v_data.coords or not 'y' in prop_v_data.coords:
                    if not 'id' in prop_v_data.dims:
                        raise ValueError("Expected `id` in property_data dims "
                            "(were {}).", prop_v_data.dims)
                    else:
                        prop_v_data = prop_v_data.assign_coords({'x': x, 'y': y})

                if len(prop_v_data.dims) == 2:
                    prop_v_dim = [d for d in prop_v_data.dims if d != 'id'][0]
                    prop_x_data = prop_v_data.isel({prop_v_dim: 0})
                    prop_y_data = prop_v_data.isel({prop_v_dim: 1})
                elif len(prop_v_data.dims) == 1:
                    prop_x_data = np.cos(prop_v_data)
                    prop_y_data = np.sin(prop_v_data)
                else:
                    raise ValueError("Expected dict with 2 entries (x and y) "
                                     "or 1 entry (angle), but received {}!"
                                     "".format(prop_v_data))

                _v_property_kwargs = dict(angles='xy', scale_units='xy',
                                          scale=3.)
                if vector_property_kwargs is not None:
                    _v_property_kwargs.update(vector_property_kwargs)

                hlpr.ax.quiver(x, y, prop_x_data, prop_y_data,
                               **_v_property_kwargs)
                if vector_property_is_nematic:
                    hlpr.ax.quiver(x, y, -prop_x_data, -prop_y_data,
                                **_v_property_kwargs)




            hlpr.provide_defaults('set_title', title="Time {}".format(time))
            hlpr.provide_defaults('set_labels',
                               x=r"$x \ [A_0^{1/2}]$", y=r"$y \ [A_0^{1/2}]$")

            if (data['periodic']):
                if abs(curvature) > 1.e-12:        
                    R = 1. / curvature
                    Ox = Lx / 2.
                    Oy = Ly / 2. - R
                    
                    if max_theta < m.pi/2. - 0.2:
                        hlpr.provide_defaults('set_limits', 
                            x=((R + Ly/2.) * np.sin(-max_theta)+Lx/2,
                            (R + Ly/2.) * np.sin( max_theta)+Lx/2),
                            y=((R - Ly/2.) * np.cos(-max_theta) - R + Ly/2,
                                Oy + R + Ly/2.))
                    else:
                        hlpr.provide_defaults('set_limits', 
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
                    hlpr.provide_defaults('set_limits', x=(0, Lx), y=(0, Ly))


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
                hlpr.provide_defaults('set_limits',
                                      x=(domain_size_min_x, domain_size_max_x),
                                      y=(domain_size_min_y, domain_size_max_y))
            
            if set_limits is not None:
                hlpr.provide_defaults('set_limits', **set_limits)
            
            hlpr.ax.set_aspect('equal')

            # end update here
            yield

    hlpr.register_animation_update(update, invoke_helpers_before_grab=True)
