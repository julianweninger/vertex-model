"""BendSplay-model specific plot function for the state"""

import logging
from typing import Tuple, Union

import math as m
import numpy as np
import xarray as xr
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt


from utopya import DataManager, UniverseGroup
from utopya.plotting import UniversePlotCreator, is_plot_func, PlotHelper


from ..tools import save_and_close

log = logging.getLogger(__name__)
# -----------------------------------------------------------------------------

@is_plot_func(creator_type=UniversePlotCreator, supports_animation=True)
def vector_field(dm: DataManager, *,uni: UniverseGroup, hlpr: PlotHelper,
                 datapath: str='BendSplay',
                 quiver_kwargs: dict=None, N_rho: int=10, N_phi: int=30):

    # Get the group that all datasets are in
    grp = uni['data/'+datapath]


    # Prepare the figure ......................................................
    # Prepare the figure to have as many columns as there are properties
    hlpr.setup_figure()

    def update():
        data = grp['angles'].data
        coords_x = grp['coords_x'].data
        coords_y = grp['coords_y'].data

        scale = 10
        _quiver_kwargs = dict(scale=scale, scale_units='xy', color='black')
        if quiver_kwargs:
            _quiver_kwargs.update(quiver_kwargs)


        for time in data.time:
            hlpr.ax.clear()

            x = coords_x.sel(time=time)
            y = coords_y.sel(time=time)

            angles = data.sel(time=time)
            phi = np.arctan2(y, x)

            qx = np.cos(angles) * np.cos(phi) - np.sin(angles) * np.sin(phi)
            qy = np.cos(angles) * np.sin(phi) + np.sin(angles) * np.cos(phi)


            hlpr.ax.quiver(x - qx/2./scale, y - qy/2./scale, qx, qy,
                           **_quiver_kwargs)

            
            hlpr.ax.set_aspect('equal')

            hlpr.ax.set_title("Time {}".format(time.data))

            # end update here
            yield

    hlpr.register_animation_update(update)
