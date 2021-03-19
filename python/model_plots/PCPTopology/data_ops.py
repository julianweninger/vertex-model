"""Generic data operations"""

import logging
import warnings
from typing import Tuple

import numpy as np
import xarray as xr

# -----------------------------------------------------------------------------

def stack(d, *a, **k):
    """Forward to xarray.DataArray.stack"""
    return d.stack(*a, **k)