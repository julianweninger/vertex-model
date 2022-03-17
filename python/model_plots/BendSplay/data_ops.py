"""Generic data operations"""

import logging
import warnings

import numpy as np
import pandas as pd
import xarray as xr
import math

log = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
def to_dataframe(d: xr.DataArray):
    """Converts a xarray.DataArray to pandas dataframe
    and transfers all dimensions
    """
    _d = d.to_dataframe().reset_index()
    
    rename = {}
    for i, c in enumerate(_d.columns):
        _c = c.split("_")
        if _c[0] == "level":
            rename[c] = d.dims[int(_c[1])]
    
    _d = _d.rename(rename, axis=1)

    return _d