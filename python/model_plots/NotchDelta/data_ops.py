"""Implements basic generic multiverse plot functions"""

import logging
from typing import Union, Tuple, List
import copy

import numpy as np
import xarray as xr
import matplotlib as mpl
import os

from dantro.utils.data_ops import count_unique

from utopya import DataManager, MultiverseGroup
from utopya.plotting import is_plot_func, PlotHelper, MultiversePlotCreator
from utopya.dataprocessing import transform

from utopya.plot_funcs._basic import _errorbar


# Get a logger
log = logging.getLogger(__name__)

def my_histogram(data, dims: List[str]=None, **kwargs) -> xr.DataArray:
    """Applies np.histogram to the given data and constructs a xr.DataArray for
    the results.

    The values are left shifted within the bin_edges, i.e. a 0 is appended in
    hist.

    Args:
        data: The data
        dims (list of str, optional): The dimensions along which to apply
            np.unique. The other dimensions will be available after the
            operation.
        **kwargs: Passed on to np.histogram
    """
    def _histogram(data, **kwargs) -> xr.DataArray:
        data = data[~np.isnan(data)]
        hist, bin_edges = np.histogram(data, **kwargs)

        hist = np.append(hist, [0])
        return xr.DataArray(data=hist,
                            name=data.name + " (counts)",
                            dims=('bins',),
                            coords=dict(bins=bin_edges))
    
    if not dims:
        return _histogram(data, **kwargs)
    
    dims = [dim for dim in data.dims if dim not in dims]
    return data.stack(z=dims).groupby('z').map(_histogram, **kwargs).unstack('z')

def replace_dim(data: xr.DataArray, coords: xr.DataArray,
                dim: str, new_dim: str=None) -> xr.DataArray:
    """Assign new coordinates to this object.

    Returns a new object with all the original data in addition to the new
    coordinates.

    If data has single dimension, performs
    data.assign_coords({dim: coords}).rename({dim: new_dim}).

    If data has multiple dimensions, split-apply-combine is used to do this on
    every coordinate except dim; this may cause the data to increase in length
    along this dimension, because filled with NaNs.

    args:
        data (xr.DataArray): The data where the replace the dimension
        coords (xr.DataArray): The coords of the new dimension. Needs to have
            same dimensions as data.
        dim (str): The dim along which to assign the new coords
        new_dim (str, optional): The name of the new dimension. If None the name
        of coords is used.
    """
    def _map_density(data: xr.Dataset, data_variable: str,
                     coords_variable: str,
                     stack_dim: str) -> xr.DataArray:
        """
        """
        # have to preserve the stack dim coordinates for unstacking
        if len(data.dims) > 1:
            return xr.DataArray(data=data[data_variable].data, name=data_variable,
                dims=(stack_dim, coords_variable, ),
                coords={stack_dim: data[stack_dim],
                        coords_variable: data[coords_variable].data})
        else:
            return xr.DataArray(data=data[data_variable].data, name=data_variable,
                dims=(coords_variable, ),
                coords={coords_variable: data[coords_variable].data})

    if not new_dim:
        new_dim = coords.name
    
    dataset = xr.Dataset({data.name: data, new_dim: coords},
                          coords={coord: data[coord] for coord in data.coords})

    # Use split-apply-combine to reduce the data along `map_dim`
    split_dims = [d for d in data.dims if d != dim]
    if len(split_dims) > 0:
        dataset = dataset.stack(_stack_rd=split_dims).groupby('_stack_rd')

        data = dataset.map(_map_density,
                        data_variable=data.name,
                        coords_variable=new_dim,
                        stack_dim='_stack_rd')
        data = data.unstack('_stack_rd')
    else:
        data = data.assign_coords({dim: coords})
        data = data.rename({dim: new_dim})
    
    return data