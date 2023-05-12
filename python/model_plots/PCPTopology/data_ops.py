"""Generic data operations"""

import logging
import warnings
from typing import Tuple

import numpy as np
import xarray as xr
import pandas as pd
import math

log = logging.getLogger(__name__)

# -----------------------------------------------------------------------------

def stack(d, *a, **k):
    """Forward to xarray.DataArray.stack"""
    return d.stack(*a, **k)
    
def unstack(d, *a, **k):
    """Forward to xarray.DataArray.unstack"""
    return d.unstack(*a, **k)

def rolling(d: xr.DataArray, *a, **k):
    """Forward to xarray.DataArray.rolling"""
    return d.rolling(*a, **k)

def reindex(d: xr.DataArray, *a, **k):
    """Forward to xarray.DataArray.reindex"""
    return d.reindex(*a, **k)

def dropna(d: xr.DataArray, *a, **k):
    """Forward to xarray.DataArray.dropna"""
    return d.dropna(*a, **k)

def drop_duplicates(d: xr.DataArray, *a, **k):
    """Forward to xarray.DataArray.drop_duplicates"""
    return d.drop_duplicates(*a, **k)

def sortby(d: xr.DataArray, *a, **k):
    """Forward to xarray.DataArray.sortby"""
    return d.sortby(*a, **k)

def apply(d, f, *a, **k):
    return d.apply(f, *a, **k)

def to_dataframe(d: xr.DataArray):
    """Converts a xarray.DataArray to pandas dataframe
    and transfers all dimensions
    """
    _d = d.to_dataframe()
    
    rename = {}
    for i, c in enumerate(_d.columns):
        _c = c.split("_")
        if _c[0] == "level":
            rename[c] = d.dims[int(_c[1])]
        if c in _d.index.names:
            _d = _d.drop(columns=c)
    
    _d = _d.rename(rename, axis=1)

    return _d.reset_index()


def classify_edges(edges):
    """Classify edges according to the cell type on either side.

    Classification:
        - 0: same type, except hair-hair
        - 1: different type
        - 2: hair-hair
        - (-1): boundary

    """
    type_alpha = edges.sel(property='type_alpha')
    type_beta = edges.sel(property='type_beta')

    # Same type edges are type 0
    type = xr.zeros_like(type_alpha, dtype='int64')
    type = type.where(type_alpha == type_beta, 1)

    # map H-H edges to 2
    type += 2 * (type_alpha == 1)

    # map boundary edges to -1
    type = type.where(type_alpha >= 0, -1)
    type = type.where(type_beta >= 0, -1)

    return type


def displacement(TimeSeries, dim, coords, *, shift=1):
    """Calculate the displacement in periodic space space of a TimeSeries.
    Calculates the difference of data with itself:
    displ = data - data.shift({'time': shift})

    And corrects for periodic boundary conditions. Maximum displacement between
    two displacements must be < 1/2 domain size.


    Args:
        TimeSeries: The Time series containing the data and time information
        dim: The dimension in TimeSeries from which to select the data
        coords: The coordinates in dim to select; 
                i.e. data = TimeSeries.sel({dim: coords})
    """
    space = [[TimeSeries[t].attrs['Lx'][0], TimeSeries[t].attrs['Ly'][0]]\
             for t in TimeSeries]
    time = [float(t) for t in TimeSeries]

    data = TimeSeries.sel({dim: coords})
    
    space = xr.DataArray(space, dims=['time', dim],
                         coords={dim: coords, 'time': time})
    space = space.reindex({'time': data.time})

    # use shifted data before as reference
    displ = data - data.shift({'time': shift})

    # correct periodicity
    displ = xr.where(displ > -space / 2., displ, displ + space)
    displ = xr.where(displ <  space / 2., displ, displ - space)

    # cumulate displacement over time
    displ = (displ - displ.mean('id')).cumsum('time')

    return displ


def map_directional_time(d: xr.DataArray):
    if 'direction' not in d.coords:
        raise RuntimeError("'direction' must be in dims of data {}", d.coords)
    if 'direction' not in d.dims:
        d = d.expand_dims(dim='direction')
    d = d.rename(time="_time")
    d = d.stack(time=("_time", "direction"))

    
    _time = d._time
    direction = d.direction

    time = (  d._time * (d.direction == 'forward') 
            - 1. * d._time * (d.direction == 'backward')
            - 1.e-5 * (d._time == 0) * (d.direction == 'backward'))
            
    d = d.assign_coords(time=time.data)
    d = d.assign_coords(_time=("time", _time.data))
    d = d.assign_coords(direction=("time", direction.data))

    return d.sortby("time")

def groupby_bins(HC, *, bins: int, **kwargs):
    """Forward to xarray.DataArray.groupby_bins using the beginning of the bins
    as coordinate """
    if kwargs.pop("labels", None) is not None:
        log.warning("Overwriting labels entry!")

    kwargs["labels"] = np.linspace(HC.x.min(), HC.x.max(), bins, endpoint=False)
    return HC.groupby_bins("x", bins=bins, **kwargs)

def spatial_binning(x: xr.DataArray,
                    data: xr.DataArray,
                    co_data: xr.DataArray=None,
                    *, bins: int):
    """Group spatial system into bins by their x-axis
    
    x: the coordinates along x-axis
    data: The data to be grouped
    co_data (optional): A second set of data
    bins (int): The number of bins
    """
    _bins = np.linspace(0.5 / bins, 1. - 0.5 / bins, bins)

    x = to_dataframe(x)

    data = to_dataframe(data)
    data['x'] = x['x']
    data = data.dropna(subset='x')

    data['x_bins'] = data['x'].apply(
        lambda x: _bins[min(int(np.floor(x * _bins.size)), _bins.size-1)]
    )

    if co_data is not None:
        __co_data = to_dataframe(co_data)
        data[co_data.name] = __co_data[co_data.name]
    

    result = None
    for time in data.time.unique():
        frame = data.loc[data['time'] == time]

        __result = frame._get_numeric_data().groupby(by='x_bins').mean()
        __result = __result.reset_index()
        __result['time'] = time

        if result is None:
            result = __result
        else:
            result = pd.concat([result, __result])
    
    result = result.drop(columns='id')
    if 'seed' in result.columns:
        result = result.drop(columns='seed')

    if co_data is not None:
        return result.sort_values(co_data.name).reset_index()

    else:
        return result.reset_index()
    

def map_stage(data: xr.DataArray, area: xr.DataArray, *,
              exp_area: xr.DataArray=None, **kwargs_interp):
    if not exp_area:
        exp_area = xr.DataArray([0.58, 0.88, 0.99, 1.39, 3.09, 3.25, 3.47],
                                dims=("stage"),
                                coords={'stage': [8, 9, 10, 12, 14, 15, 17]})
    
    # drop duplicates using np.unique 
    # area = area.isel(time=np.unique(area.data, return_index=True)[1])

    time_to_area = xr.DataArray(area.time, dims="area",
                                coords={"area": area.data})
    time_to_stage = time_to_area.interp(area=exp_area, method="nearest",
                                        **kwargs_interp)
    
    return data.sel(time=time_to_stage.dropna(dim="stage"))

def align_polarity(data: xr.DataArray):
    def flip(d: xr.DataArray):
        mean = d.mean(dim='id')
        if (mean.data < 0):
            d = (d + math.pi) % (2 *math.pi)
        else:
            d = (d + 2 * math.pi) % (2 * math.pi)
        return d
    data = data.unstack(dim='ids')
    data = data.groupby('seed').apply(flip)
    return data


def concat_exp_coordination(data: xr.DataArray, *, coordins: xr.DataArray=None):
    """Define a dataset with the experimental data from the basilar pappillar.
    Contains coordination numbers at different stages.

    Credits to Anubhav Prakash (2020, NCBS).
    """
    if not coordins:
        coordins = xr.DataArray([
                [5, 5, 4, 4, 4, 5, 4, 4, 4, 4, 5, 4, 4, 5, 5, 4, 4, 3, 4, 4, 4, 4, 4, 5, 4, 3, 4, 5, 4, 5, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 5, 4, 5, 4, 4, 4, 4, 3, 3, 4, 5, 4, 4, 5, 4, 4, 5, 4, 5, 3, 5, 4, 5, 4, 3, np.nan, np.nan, np.nan, np.nan],
                [5, 5, 5, 5, 5, 5, 4, 5, 6, 5, 5, 6, 5, 5, 5, 5, 4, 5, 5, 6, 5, 4, 5, 5, 5, 5, 4, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 5, 5, 4, 6, 4, 6, 5, 5, 5, 5, 5, 4, 5, 5, 5, 6, 5, 5, 5, 5, 4, 5, 5, 6, 6, 5, 5, 5, 5, 6, 5, 4, 5, 6, 5, 5, 5, 5],
                [6, 6, 6, 5, 6, 7, 6, 6, 6, 7, 7, 6, 7, 6, 6, 6, 7, 6, 5, 7, 6, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 6, 6, 7, 6, 5, 5, 6, 5, 6, 5, 6, 4, 6, 6, 6, 5, 6, 6, 6, 5, 5, 6, 6, 5, 6, 5, 6, 6, 6, 6, 6, 6, 5, 6, 6, 6, 5, 6, 6, 6, 6, 6, 6, 5, np.nan, np.nan, np.nan, np.nan],
                [9, 10, 9, 9, 9, 9, 9, 9, 10, 8, 9, 9, 8, 9, 9, 9, 10, 9, 9, 10, 9, 9, 9, 9, 10, 9, 10, 9, 9, 10, 9, 11, 9, 10, 9, 10, 9, 10, 9, 9, 8, 9, 9, 9, 10, 9, 9, 9, 9, 8, 9, 8, 9, 9, 8, 8, 9, 8, 9, 9, 9, 9, 9, 9, 9, 9, 8, 10, 9, 9, 9, 9, 9, 9, 9, 9, np.nan, np.nan, np.nan, np.nan]
            ],
            name="num_neighbors",
            dims=["stage", "id"],
            coords={'stage': [8, 10, 12, 15], "id": range(80)}
        )
        coordins = coordins.expand_dims("type").assign_coords(type=["experiment"])
        
        if "ids" in data.dims:
            coordins = coordins.expand_dims("seed").assign_coords(seed=[0])
            coordins = coordins.stack(ids=["id", "seed"])

    data = data.expand_dims("type").assign_coords(type=["simulation"])
    data = to_dataframe(data)
    coordins = to_dataframe(coordins)

    return pd.concat([data, coordins])

def define_exp_coordin_area():
    """Define a dataset with the experimental data from the basilar pappillar.
    Contains, averages of Coordination number and area at different stages 
    E8-17 and different positions on BP.

    Credits to Anubhav Prakash (2020, NCBS).
    """
    num_neighbors = xr.DataArray([
            4.184, np.nan, 4.795, 6.250, np.nan, 8.395, np.nan,
            4.250, np.nan, 5.013, 5.947, np.nan, np.nan, np.nan,
            4.250, np.nan, 5.013, 5.947, np.nan, 9.066, np.nan,
            4.074, np.nan, 5.649, 6.948, np.nan, np.nan, np.nan,
            4.074, np.nan, 5.649, 6.948, np.nan, np.nan, np.nan,
            4.341, np.nan, 5.169, np.nan, np.nan, np.nan, np.nan,
            4.341, np.nan, 5.169, np.nan, np.nan, np.nan, np.nan,
            4.303, np.nan, 4.886, np.nan, np.nan, np.nan, np.nan
        ],
        name="num_neighbors",
        dims=["area"],
        coords={"area": [
            0.773, 0.796, 1.252, 1.536, 2.800, 3.018, 3.376,
            0.444, 0.800, 0.959, 1.202, 2.749, 3.020, 3.487,
            0.576, 0.876, 0.992, 1.386, 3.089, 3.249, 3.474,
            0.431, 0.762, 1.197, 1.072, np.nan, np.nan, 2.780,
            0.539, 0.700, 1.159, 1.607, np.nan, np.nan, np.nan,
            0.509, 0.864, 0.861, np.nan, np.nan, np.nan, np.nan,
            0.512, 0.848, 1.224, np.nan, np.nan, np.nan, np.nan,
            0.491, 0.736, 1.053, np.nan, np.nan, np.nan, np.nan
        ]}
    )

    num_neighbors = num_neighbors.assign_coords(
        stage=("area", [
            8, 9, 10, 12, 14, 15, 17,
            8, 9, 10, 12, 14, 15, 17,
            8, 9, 10, 12, 14, 15, 17,
            8, 9, 10, 12, 14, 15, 17,
            8, 9, 10, 12, 14, 15, 17,
            8, 9, 10, 12, 14, 15, 17,
            8, 9, 10, 12, 14, 15, 17,
            8, 9, 10, 12, 14, 15, 17
        ]),
        position=("area", [
            "0", "0", "0", "0", "0", "0", "0",
            "25 S", "25 S", "25 S", "25 S", "25 S", "25 S", "25 S", 
            "25 I", "25 I", "25 I", "25 I", "25 I", "25 I", "25 I", 
            "50 S", "50 S", "50 S", "50 S", "50 S", "50 S", "50 S", 
            "50 I", "50 I", "50 I", "50 I", "50 I", "50 I", "50 I", 
            "75 S", "75 S", "75 S", "75 S", "75 S", "75 S", "75 S", 
            "75 I", "75 I", "75 I", "75 I", "75 I", "75 I", "75 I", 
            "100", "100", "100", "100", "100", "100", "100"
        ]))

    return num_neighbors

def read_csv(*a, **k):
    return pd.read_csv(*a, **k)

def stack_xls_sheets(path: str, *, label: str, label_sheets: str,
                     label_columns: str, 
                     reduce_mean: bool=False,
                     convert_sheet_labels_to_float: bool=False, **kwargs):
    """Read a excel (xls) file and stack sheets. Returns a pandas dataframe.
    """
    if path is None:
        log.warning("No path provided to stack_xls_sheets! Returning empty "
                    "dataframe.")
        return pd.DataFrame({label_sheets: [],
                             'id': [],
                             label_columns: [],
                             label: []})

    data = pd.read_excel(path, sheet_name=None, **kwargs)

    def stack_sheet(d):
        if reduce_mean:
            d = d.mean(axis=0)
            d = d.reset_index()
            d = d.rename(columns={"index": label_columns, 0: label})

        else:
            d = d.stack()
            d = d.reset_index()
            d = d.rename(columns={"level_0": "id", "level_1": label_columns,
                                  0: label})
            d = d.drop(columns=["id"])
        
        if convert_sheet_labels_to_float:
            d[label_columns] = pd.to_numeric(d[label_columns], downcast="float")

        return d


    data = pd.concat([stack_sheet(data[k]) for k in data.keys()],
                     keys=data.keys())
    data = data.reset_index()
    data = data.rename(columns={"level_0": label_sheets, "level_1": "id"})

    if reduce_mean:
        data = data.drop(columns=["id"])

    return data

def dataframe_where(d, *, v, k):
    """Select data where d[key] == value"""
    return d[d[k] == v]

def concat_exp_and_data(exp, d):
    """Concatenate two dataframes. Applies only to particular arrangement of
    data.
    """
    exp = exp.drop(columns=['type', 'id'])
    stage_coords = exp['stage']
    exp = xr.DataArray(exp['num_neighbors'], dims=['id'])
    exp = exp.assign_coords(stage=('id', stage_coords)).sortby('stage')

    d = d.drop(columns=['time', '_time', 'direction', 'area', 'seed', 'id'])
    stage_coords = d['stage']
    d = xr.DataArray(d['num_neighbors'], dims=['id'])
    d = d.assign_coords(stage=('id', stage_coords)).sortby('stage')

    d = d.expand_dims("type").assign_coords(type=["simulation"])
    exp = exp.expand_dims("type").assign_coords(type=["experiment"])

    return xr.concat([d, exp], dim="type")