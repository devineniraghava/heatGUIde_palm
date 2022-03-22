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

#%% Write netCDF data

file = "test.nc"
# Create file
ds = nc.Dataset(file, "w", format="NETCDF4")
# Create Dimensions
time = ds.createDimension("time", None)

lat = ds.createDimension("lat", 10)
lon = ds.createDimension("lon", 10)
# Create Variables 
times = ds.createVariable("time", "f4", ("time",))
lats = ds.createVariable("lat", "f4", ("lat",))
lons = ds.createVariable("lon","f4",("lon", ))
value = ds.createVariable("value", "f4", ("time", "lat", "lon",))
# Create metadata for variables
value.units = "Unknown"

# Populate Dimensions
lats[:] = np.arange(40,50,1)
lons[:] = np.arange(-110, -100, 1)

# Populate Variables
value[0,:,:] = np.random.uniform(0, 100, size=(10,10))

xval = np.linspace(0.5, 5, 10)
yval = np.linspace(0.5, 5, 10)

value[1,:,:] = np.array(xval.reshape(-1,1) + yval)

# Close netCDF file
ds.close()


#%% Read netCDF data

fn = "example_static_file.nc"

ds = nc.Dataset(fn)


















#%% Xarray

ds = xr.open_dataset("example_static_file.nc", engine="netcdf4")

















#%%














