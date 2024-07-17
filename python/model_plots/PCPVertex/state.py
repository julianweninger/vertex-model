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
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from scipy.interpolate import griddata
import skimage

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
    property_interpolation_plot_kwargs: dict=None,
    property_resolution: Tuple[int, int]=[512,512],
    quiver_kwargs: dict=None,
    vector_property_kwargs: dict=None,
    vector_property_is_nematic: bool=False,
    set_limits: dict=None,
    semi_periodic_space: bool=False,
):
    """Performs a plot of the cells, edges and vertices
    
    Args:
        data: dict containing the required data. 
            Required keys: 'Vertices', 'Edges', 'Cells', 'periodic', 'Lx', 'Ly'.
            Optional keys: 'skew_x', 'skew_y'.
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
        # the domain extent for non periodic boundaries (centered on (0., 0.))
        domain_size_min_x = 0.
        domain_size_max_x = 0.
        domain_size_min_y = 0.
        domain_size_max_y = 0.

        Vertices = data['Vertices'].transpose('id', 'time', 'property')
        Edges = data['Edges'].transpose('id', 'time', 'property')
        Cells = data['Cells'].transpose('id', 'time', 'property')

        if (not data['periodic']):
            domain_size_min_x = Vertices.sel(property="x").min().data
            domain_size_max_x = Vertices.sel(property="x").max().data
            domain_size_min_y = Vertices.sel(property="y").min().data
            domain_size_max_y = Vertices.sel(property="y").max().data
        elif semi_periodic_space:
            xs = Vertices.sel(property='x')
            xs = xr.where(xs > data['Lx'] / 2, xs - data['Lx'], xs)
            Vertices.loc[:, :, 'x'] = xs - xs.mean()

            xs = Cells.sel(property='x')
            xs = xr.where(xs > data['Lx'] / 2, xs-data['Lx'], xs)
            Cells.loc[:, :, 'x'] = xs - xs.mean()

            domain_size_min_x = Vertices.sel(property="x").min().data
            domain_size_max_x = Vertices.sel(property="x").max().data
            domain_size_min_y = Vertices.sel(property="y").min().data
            domain_size_max_y = Vertices.sel(property="y").max().data

        if select_times is not None:
            times = select_times
        else:
            times = np.unique(Vertices.coords['time'])

        for time in times:
            hlpr.fig.clear()
            hlpr.ax.set_aspect('auto')

            v_data = Vertices.sel(time=time).dropna(dim='id')
            e_data = Edges.sel(time=time).dropna(dim='id')
            c_data = Cells.sel(time=time).dropna(dim='id')

            Lx = data['Lx'].sel(time=time).data
            Ly = data['Ly'].sel(time=time).data
            if data['periodic'] and not semi_periodic_space:
                domain_size_max_x = Lx
                domain_size_max_y = Ly

            if set_limits is not None:
                __set_limits = set_limits.copy()
                __set_limits['x'] = __set_limits.get(
                    'x', (domain_size_min_x, domain_size_max_x))
                __set_limits['y'] = __set_limits.get(
                    'y', (domain_size_min_y, domain_size_max_y))

            else:
                __set_limits = dict()
                __set_limits['x'] = (domain_size_min_x, domain_size_max_x)
                __set_limits['y'] = (domain_size_min_y, domain_size_max_y)
            
            # periodic skewed boundary condition
            if 'skew_x' in data:
                skew_x = data['skew_x'].sel(time=time).data
            else:
                skew_x = 0.
            if 'skew_y' in data:
                skew_y = data['skew_y'].sel(time=time).data
            else:
                skew_y = 0.


            if not semi_periodic_space:
                L = max(plot_excess_length, max(skew_x, skew_y))
                if skew_x > Lx or skew_y > Ly:
                    log.warn("Skew exceeds the domain size. Not plotting "
                             "periodicity correctly in time {}.".format(time))

            else: 
                L = max(
                    abs(np.nanmin(v_data.sel(property="x").data)),
                    plot_excess_length
                )


            ### plot vertices
            if plot_vertices:
                hlpr.ax.scatter(v_data.sel(property="x"),
                                v_data.sel(property="y"),
                                c="black")
                
                if data['periodic']:
                    ax = v_data.sel(property="x").data
                    ay = v_data.sel(property="y").data

                    hlpr.ax.scatter(np.where(ax < L, ax, np.nan) + Lx,
                                    ay + skew_y,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ax > Lx - L, ax, np.nan) - Lx,
                                    ay - skew_y,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ay < L, ax, np.nan) + skew_x,
                                    ay + Ly,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ay > Ly - L, ax, np.nan) - skew_x,
                                    ay - Ly,
                                    c="grey")

                    hlpr.ax.scatter(np.where(ax < L, ax, np.nan) + Lx + skew_x,
                                    np.where(ay < L, ay, np.nan) + Ly + skew_y,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ax>Lx-L, ax, np.nan) - Lx + skew_x,
                                    np.where(ay < L, ay, np.nan) + Ly - skew_y,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ax>Lx-L, ax, np.nan) - Lx - skew_x,
                                    np.where(ay>Ly-L, ay, np.nan) - Ly - skew_y,
                                    c="grey")
                    hlpr.ax.scatter(np.where(ax < L, ax, np.nan) + Lx - skew_x,
                                    np.where(ay>Ly-L, ay, np.nan) - Ly + skew_y,
                                    c="grey")



            ### plot edges
            # the ids of the vertices a and b of the edges
            vertex_a = e_data.sel(property="vertex_a").dropna(dim='id').astype(int)
            vertex_b = e_data.sel(property="vertex_b").dropna(dim='id').astype(int)

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
                
                dx = bx - ax
                dy = by - ay

                add_skew_x = - np.round(dy/Ly) * skew_x
                add_skew_y = - np.round(dx/Lx) * skew_y

                dx = bx + add_skew_x - ax
                dy = by + add_skew_y - ay

                return (dx - np.round(dx/Lx) * Lx,
                        dy - np.round(dy/Ly) * Ly)


            def quiver_and_colors(x, y, dx, dy, *, colorbar=False, axis=hlpr.ax):
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

                        quiver = axis.quiver(*quiver_args, **_quiver_kwargs)

                        if colorbar:
                            edges_cbar = hlpr.fig.colorbar(
                                quiver, ax=axis, extend='both',
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
                            axis.quiver(
                                _x+shift_x, _y+shift_y, dx / N, dy / N, 
                                e_prop_data.sel({e_prop_dim: name_0}),
                                **__quiver_kwargs)

                            quiver = axis.quiver(
                                _x-shift_x, _y-shift_y, dx / N, dy / N, 
                                e_prop_data.sel({e_prop_dim: name_1}),
                                **__quiver_kwargs)
                        
                        if colorbar:
                            edges_cbar = hlpr.fig.colorbar(quiver, ax=axis,
                                                     extend='both')
                            edges_cbar.set_label(e_prop_data.name)
                            edges_cbar.minorticks_on()
                    else:
                        raise

                else:
                    quiver = axis.quiver(*quiver_args, **_quiver_kwargs)

            def plot_edges(axis=hlpr.ax, colorbar=True):
                dx, dy = displacement(ax, ay, bx, by)
                
                quiver_and_colors(
                    ax, ay, dx, dy, colorbar=colorbar, axis=axis
                )

                ### plot duplicates of periodic edges
                if data['periodic']:
                    _ax = xr.where(ax < L, ax, np.nan) + Lx
                    _dx, _dy = displacement(_ax, ay + skew_y, bx + Lx, by + skew_y)
                    quiver_and_colors(
                        _ax, ay + skew_y, _dx, _dy, colorbar=False, axis=axis
                    )

                    _ax = xr.where(ax > Lx - L, ax, np.nan) - Lx
                    _dx, _dy = displacement(_ax, ay - skew_y, bx - Lx, by - skew_y)
                    quiver_and_colors(
                        _ax, ay - skew_y, _dx, _dy, colorbar=False, axis=axis
                    )

                    _ax = xr.where(ay < L, ax, np.nan) + skew_x
                    _dx, _dy = displacement(_ax, ay + Ly, bx + skew_x, by + Ly)
                    quiver_and_colors(
                        _ax, ay + Ly, _dx, _dy, colorbar=False, axis=axis
                    )

                    _ax = xr.where(ay > Ly - L, ax, np.nan) - skew_x
                    _dx, _dy = displacement(_ax, ay - Ly, bx - skew_x, by - Ly)
                    quiver_and_colors(
                        _ax, ay - Ly, _dx, _dy, colorbar=False, axis=axis
                    )

                    # diagonals
                    _ax = xr.where(ax < L, ax, np.nan) + Lx + skew_x
                    _ay = xr.where(ay < L, ay, np.nan) + Ly + skew_y
                    _dx, _dy = displacement(_ax, _ay, bx+Lx+skew_x, by+Ly+skew_y)
                    quiver_and_colors(
                        _ax, _ay, _dx, _dy, colorbar=False, axis=axis
                    )

                    _ax = xr.where(ax > Lx - L, ax, np.nan) - Lx + skew_x
                    _ay = xr.where(ay < L     , ay, np.nan) + Ly - skew_y
                    _dx, _dy = displacement(_ax, _ay, bx-Lx+skew_x, by+Ly-skew_y)
                    quiver_and_colors(
                        _ax, _ay, _dx, _dy, colorbar=False, axis=axis
                    )

                    _ax = xr.where(ax > Lx - L, ax, np.nan) - Lx - skew_x
                    _ay = xr.where(ay > Ly - L, ay, np.nan) - Ly - skew_y
                    _dx, _dy = displacement(_ax, _ay, bx-Lx-skew_x, by-Ly-skew_y)
                    quiver_and_colors(
                        _ax, _ay, _dx, _dy, colorbar=False, axis=axis
                    )

                    _ax = xr.where(ax < L     , ax, np.nan) + Lx - skew_x
                    _ay = xr.where(ay > Ly - L, ay, np.nan) - Ly + skew_y
                    _dx, _dy = displacement(_ax, _ay, bx+Lx-skew_x, by-Ly+skew_y)
                    quiver_and_colors(
                        _ax, _ay, _dx, _dy, colorbar=False, axis=axis
                    )

                    # Skew exceeds domain size
                    # apply best available periodic copies
                    if skew_x > Lx or skew_y > Ly:
                        mask = ((abs(bx - ax) > Lx / 2.) | (abs(by - ay) > Ly / 2.))

                        _bx = xr.where(mask, bx, np.nan)
                        _by = xr.where(mask, by, np.nan)
                        
                        # use the inverse arrows
                        _dx, _dy = displacement(bx, by, ax, ay)
                        quiver_and_colors((bx + _dx), (by + _dy), -_dx, -_dy,
                                          colorbar=False)

            plot_edges()

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
                x = c_data.sel(property="x")
                y = c_data.sel(property="y")

                # Assign property to every cell
                if not 'x' in prop_data.coords or not 'y' in prop_data.coords:
                    if not 'id' in prop_data.dims:
                        raise ValueError("Expected `id` in property_data dims "
                            "(were {}).", prop_data.dims)
                    else:
                        prop_data = prop_data.assign_coords({'x': x, 'y': y})
                prop_data = prop_data.dropna('id')

                # create a temporary figure that contains all the edges
                tmp_fig, tmp_ax = plt.subplots(
                    1, 1,
                    figsize=hlpr.fig.get_size_inches()
                )

                plot_edges(colorbar=False, axis=tmp_ax)

                tmp_ax.scatter(__set_limits['x'][0], __set_limits['y'][0], c='black', s=10)
                tmp_ax.scatter(__set_limits['x'][1], __set_limits['y'][1], c='black', s=10)

                tmp_ax.set_aspect('equal')
                tmp_ax.set_xlim(__set_limits['x'])
                tmp_ax.set_ylim(__set_limits['y'])

                # Render the plot to a numpy array
                tmp_fig.tight_layout()
                tmp_ax.axis('off')
                tmp_fig.tight_layout(pad=0)

                # To remove the huge white borders
                tmp_ax.margins(0)
                canvas = FigureCanvas(tmp_fig)
                canvas.draw()

                # Convert to a numpy array
                image = np.frombuffer(canvas.tostring_rgb(), dtype='uint8')
                width, height = canvas.get_width_height()
                image = image.reshape((height, width, 3))

                # Convert to grayscale
                # Using the luminosity method: 0.21*R + 0.72*G + 0.07*B
                gray_image = 0.21 * image[:, :, 0] + 0.72 * image[:, :, 1] + 0.07 * image[:, :, 2]

                # Normalize to 0-255 and convert to uint8
                intensity_image = np.uint8(gray_image)

                # crop the image to the y-axis
                sum__ = np.sum(intensity_image, axis=1).max()
                intensity_image = intensity_image[:][np.sum(intensity_image, axis=1) < sum__]
                intensity_image = intensity_image.transpose()
                sum__ = np.sum(intensity_image, axis=1).max()
                intensity_image = intensity_image[:][np.sum(intensity_image, axis=1) < sum__]
                intensity_image = intensity_image.transpose()

                cell_mask = intensity_image < intensity_image.max()
                cell_mask = np.asarray(cell_mask, dtype=int).transpose()
                cell_mask[:,  0] = 1
                cell_mask[:, -1] = 1
                cell_mask[ 0, :] = 1
                cell_mask[-1, :] = 1

                cell_properties = np.zeros_like(cell_mask, dtype=float) + np.nan

                coords_x = prop_data.x - domain_size_min_x
                coords_x /= (domain_size_max_x - domain_size_min_x)
                coords_x *= cell_mask.shape[0]
                coords_x = coords_x.astype(int)

                coords_y = prop_data.y - domain_size_min_y
                coords_y /= (domain_size_max_y - domain_size_min_y)
                coords_y *= cell_mask.shape[1]
                coords_y = cell_mask.shape[1]-1 - coords_y.astype(int)

                for p, coord_x, coord_y in zip(prop_data, coords_x, coords_y):
                    cell_properties = np.where(skimage.segmentation.flood_fill(cell_mask, (coord_x, coord_y), 2) == 2, p.data, cell_properties)
                
                cell_properties = np.where(cell_mask==0, cell_properties, np.nan)

                interpol = hlpr.ax.imshow(
                    cell_properties.transpose(),
                    extent=(domain_size_min_x, domain_size_max_x,
                            domain_size_min_y, domain_size_max_y),
                    **_property_interpolation_plot_kwargs
                )

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

            hlpr.provide_defaults('set_limits', **__set_limits)

            if skew_x > 1.e-12:
                hlpr.ax.axvline(x=skew_x, ymin=1. - 0.5 / Ly, c='gray',
                                linestyle=':', linewidth=0.2)
            elif skew_x < -1.e-12:
                hlpr.ax.axvline(x=Lx-skew_x, ymax = 0.5 / Ly, c='gray',
                                linestyle=':', linewidth=0.2)
            if skew_y > 1.e-12:                    
                hlpr.ax.axhline(y=skew_y, xmin=1. - 0.5 / Lx, c='gray',
                                linestyle=':', linewidth=0.2)
            elif skew_y < -1.e-12:                    
                hlpr.ax.axhline(y=skew_y, xmax = 0.5 / Lx, c='gray', 
                                linestyle=':', linewidth=0.2)
            
            hlpr.ax.set_aspect('equal')

            # end update here
            yield

    hlpr.register_animation_update(update, invoke_helpers_before_grab=True)
