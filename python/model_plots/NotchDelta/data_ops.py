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

def count_unique(data, dims: List[str]=None) -> xr.DataArray:
    """Applies np.unique to the given data and constructs a xr.DataArray for
    the results.

    NaN values are filtered out.

    NOTE this is a tmp copy from dantro 13.1 _(WIP)_

    Args:
        data: The data
        dims (List[str], optional): The dimensions along which to apply
            np.unique. The other dimensions will be available after the
            operation. If not provided it is applied along all dims.

    """
    def _count_unique(data) -> xr.DataArray:
        unique, counts = np.unique(data, return_counts=True)
        
        # remove np.nan values
        # NOTE np.nan != np.nan, hence np.nan will count 1 for every occurrence,
        #      but duplicate values not allowed in coords.
        counts = counts[~np.isnan(unique)]
        unique = unique[~np.isnan(unique)]

        if isinstance(data, xr.DataArray):
            name = data.name + " (unique counts)"
        else:
            name = "data (unique counts)"

        # Construct a new data array and return
        return xr.DataArray(data=counts,
                            name=name,
                            dims=('unique',),
                            coords=dict(unique=unique))
    
    if not dims:
        return _count_unique(data)

    if not isinstance(data, xr.DataArray):
        raise TypeError("Data needs to be of type xr.DataArray, but was "
                        f"{type(data)}!")
    
    # use split-apply-combine along those dimensions not in dims
    split_dims = [dim for dim in data.dims if dim not in dims]
    
    if len(split_dims) == 0:
        return _count_unique(data)

    data = data.stack(_stack_cu=split_dims).groupby('_stack_cu')
    return data.map(_count_unique).unstack('_stack_cu')

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

