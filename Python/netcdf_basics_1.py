#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 14:51:06 2022

@author: rdevinen
"""

import netCDF4 as nc

path = '/home/rdevinen/palm/current_version/JOBS/test_urban_2/INPUT/test_urban_dynamic'

# path = '/home/rdevinen/Documents/PALM_Seminar/NetCDF/Python/example_static_file.nc'

# path = '/home/rdevinen/Desktop/compare/test_urban_dynamic'

ds = nc.Dataset(path)

print(ds.dimensions)

# a = ds.variables['u'][:]

# b = a[20]






#%%


