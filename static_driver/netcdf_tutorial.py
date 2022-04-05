#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 16:54:44 2022

@author: rdevinen
"""
# https://www.youtube.com/watch?v=VH-PCQ991fw
# https://www.youtube.com/watch?v=mnoGmS2XtDg
import datetime
import math
import pytz

# from netCDF4 import Dataset
import numpy as np

import netCDF4 as nc
import xarray as xr
# from warnings import filterwarnings
# filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')


#%%


import xarray as xr
import matplotlib as plt


#%%
odr = xr.open_dataset(
    "/home/rdevinen/Documents/GitHub/heatGUIde_palm/static_driver/files/testing_static", engine="netcdf4")

ndr = xr.open_dataset(
    "/home/rdevinen/Documents/GitHub/heatGUIde_palm/static_driver/files/practice_static", engine="netcdf4")

#%%
odr["surface_fraction"].sel(nsurface_fraction=slice(0, 2)).plot(col="nsurface_fraction", col_wrap=2, robust=True);











