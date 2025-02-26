"""Implements plot functions to analyze the percolation properties of 
the NotchDelta model
"""

import logging
from typing import Union, Tuple, List
import copy

import scipy
import numpy as np
import xarray as xr
import seaborn as sns

import dantro.data_ops.db as dops

from utopya.plotting import is_plot_func, PlotHelper
from dantro.plot.funcs import facet_grid



# Get a logger
log = logging.getLogger(__name__)

def percolation_sigmoidal(xdata, threshold, slope):
    return xdata / (1 + np.exp(-(xdata - threshold) / slope))

def map_percolation_strentgh_to_density(*, cluster_id: xr.DataArray,
        cells: xr.DataArray) -> xr.Dataset:
    """Transforms a dataset of the cluster id per cell and kind per cell to
    the density of the biggest cluster. The map_dim coordinate is replaced by 
    the density at this coordinate.
    """
    
    num_cells = cluster_id.count(dim={'id'})

    clusters = dops.where(cluster_id, ">", 0)
    cluster_sizes = dops.count_unique(clusters, dims=['id'])
    # the density of the biggest cluster 
    perc = cluster_sizes.max(dim='unique') / num_cells
    
    hair_cells = dops.where(cells, "==", 1)
    num_hair_cells = hair_cells.count(dim={'id'})
    density = num_hair_cells / num_cells

    coords = dict()
    for coord in perc.coords:
        coords[coord] = perc[coord]

    data = xr.Dataset({'percolation_strength': perc, 'density': density},
                      coords=coords)
    return data

@is_plot_func(use_dag=True, required_dag_tags=('cluster_id', 'cells', ),
              supports_animation=True)
def percolation_diagram(*, data: dict, hlpr: PlotHelper, map_dim: str=None,
                        **plot_kwargs):
    """FacetGrid plot of the size of biggest cluster

    This functions manipulates the data to obtain the density of the largest
    cluster and passes the transformed data on to the generic dantro plot 
    function `facet_grid`.
    
    Args:
        data (dict): The data selected by the DAG framework
        hlpr (PlotHelper): The plot helper
        map_dim (str, optional): Map this dimension to the density of hair
            cells. Requires the dag_tag 'cells'.
        **plot_kwargs: Passed on to generic dantro plot function `facet_grid`
    """
    if 'x' in data['cells'].dims:
        data['cells'] = data['cells'].stack(id={'x', 'y'})
    if 'x' in data['cluster_id'].dims:
        data['cluster_id'] = data['cluster_id'].stack(id={'x', 'y'})
    
    num_cells = data['cluster_id'].count(dim={'id'})

    clusters = dops.where(data['cluster_id'], ">", 0)
    cluster_sizes = dops.count_unique(clusters, dims=['id'])
    # the density of the biggest cluster 
    perc = cluster_sizes.max(dim='unique') / num_cells

    if map_dim:
        if 'cells' not in data.keys():
            raise ValueError("Percolation diagram requires dag_tag 'cells' when "
                f"mapping coordinates ('map_dim'={map_dim}). Provided "
                f"dag_tags are {data.keys()}")
        
        hair_cells = dops.where(data['cells'], "==", 1)
        num_hair_cells = hair_cells.count(dim={'id'})
        density = num_hair_cells / num_cells

        if len(perc.dims) > 1:
            raise ValueError("Cannot map dimension if more than the mapping "
                "dimension available. Data has dimensions "
                f"{perc.dims}")

        perc = perc.assign_coords({map_dim: density})
        perc = perc.rename({map_dim: 'density'})

    data['data'] = perc

    facet_grid(data=data, hlpr=hlpr, **plot_kwargs)

    [theta, slope], pcov = scipy.optimize.curve_fit(percolation_sigmoidal,
                                                    perc.density,
                                                    perc.data)
    perr = np.sqrt(np.diag(pcov))

    # only plot fit if reasonable
    if theta <= 1 and theta >= 0 and slope >= 0 and slope <= 1:
        x = np.linspace(perc.density.min().data,
                        perc.density.max().data, 100)
        
        fitlabel = r"fit: $\theta$=%5.3f (%5.3f), s=%5.3f (%5.3f)" \
                    "" % (theta, perr[0], slope, perr[1])
        hlpr.ax.plot(x, percolation_sigmoidal(x, threshold=theta, slope=slope),
                     'r-', label=fitlabel)

@is_plot_func(use_dag=True, required_dag_tags=('cluster_id', 'cells', ),
              supports_animation=True)
def percolation_bifurcation(*, data: dict, map_dim: str,
                            p0: Tuple[float]=[0.4, 0.01],
                            min_perc_strength: float=0.,
                            **plot_kwargs):
    """Plot a bifurcation map of the threshold and slope in the 
    bifurcation diagram as fitted to a sigmoidal function
    `percolation_sigmoidal`) along `map_dim`.

    Args:
        data (dict): The data selected by the DAG framework
        map_dim (str): The parameter dimension to map to the density in the
            percolation diagram
        p0 (Tuple[float], default [0.4, 0.01]): Initial guess of params, passed
            on to scipy.optimize.curve_fit
        min_perc_strength (float, optional): Minimum relative percolation
            strength, below which not to fit data. Normalized to density.
        **plot_kwargs: Passed on to generic dantro plot function `facet_grid`
    """
    def _fit_sigmoidal(data: xr.Dataset, stack_dim: str=None):
        """Fits a 1D Dataset to `percolation_sigmoidal`
        """
        dims = [dim for dim in data.dims if len(data[dim]) > 1]
        if len(dims) > 1:
            raise ValueError("Expected 1d data, but dims are "
                                f"{data.dims}")
        
        [theta, slope], pcov = scipy.optimize.curve_fit(
            percolation_sigmoidal,
            data['density'].squeeze(),
            data['percolation_strength'].squeeze(),
            p0=p0)

        # filter to reasobable fits
        if theta > 1 or theta < 0 or slope > 1 or slope < 0:
            theta = np.nan
            slope = np.nan

        # filter sets of data that don't show percolation
        max_perc = (data['percolation_strength'] / data['density']).max().data
        if max_perc < min_perc_strength:
            theta = np.nan
            slope = np.nan

        # have to preserve the stack dim coordinates for unstacking
        if len(data.dims) > 1:            
            return xr.DataArray(data=[[theta, slope]], name='popt',
                                dims=(stack_dim, 'param', ),
                                coords={stack_dim: data[stack_dim],
                                        'param': ['theta', 'slope']})
        else:
            return xr.DataArray(data=[theta, slope], name='popt',
                                dims=('param', ),
                                coords={'param': ['theta', 'slope']})

    if 'x' in data['cells'].dims:
        data['cells'] = data['cells'].stack(id={'x', 'y'})
    if 'x' in data['cluster_id'].dims:
        data['cluster_id'] = data['cluster_id'].stack(id={'x', 'y'})
    
    # The density of the biggest percolation cluster and the density of hair
    # cells
    perc = map_percolation_strentgh_to_density(
        cluster_id=data['cluster_id'], cells=data['cells'])

    # Use split-apply-combine to reduce the data along `map_dim`
    split_dims = [dim for dim in perc.dims if dim != map_dim]
    grp = perc.stack(_stack_pc=split_dims).groupby('_stack_pc')
    popt = grp.map(_fit_sigmoidal, stack_dim='_stack_pc').unstack('_stack_pc')

    data['data'] = popt
    
    facet_grid(data=data, **plot_kwargs)
