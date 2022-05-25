#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 16:21:47 2022

@author: rdevinen

Note:
This is the PART 1 simulation of CAMPUS OFFENBURG..
The steps of the static driver are:
1) Create RIZ in green area
2) Input Pavement
3) Input Trees
    
This script is mostly a practice to get correct masked data for pavement and vegitation
The idea is to implement the correct veg and pavement to the netcdf files, wherer building is masked 
The sum of the vegitation and pavement should be equal to total, maksing the building.
This is implemented in the this version 5

"""

import datetime
from PIL import Image
from numpy import asarray
import pytz

from netCDF4 import Dataset
import numpy as np
import xarray as xr
# import hvplot.xarray
import numpy.ma as ma

from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

def prRed(skk): print("\033[31;1;m {}\033[00m" .format(skk)) 
def prYellow(skk): print("\033[33;1;m {}\033[00m" .format(skk)) 

import os
try:
    os.remove("files/riz_4_static")
    print("file deleted")
except FileNotFoundError:
    print("file not found")


#%% netCDF main attributes
#%%% Optional attribbutes
path = '/home/rdevinen/Documents/GitHub/heatGUIde_palm/static_driver/files/riz_5_static'
path1 = '/home/rdevinen/Documents/GitHub/heatGUIde_palm/static_driver/files/part1.png'
nc_file = Dataset(path, 'w', format='NETCDF4')
nc_file.title = 'HS Offenburg PART 1'
nc_file.author = 'devineni and haag'
nc_file.institution = 'Institut für nachhaltige Energiesysteme,' \
      'Hochschule Offenburg'
nc_file.comment = 'Generic crossing example'
nc_file.creation_date = \
    pytz.utc.localize(datetime.datetime.utcnow()).strftime('%y-%m-%d %H:%M:%S %z')[:-2]
nc_file.history = ''
nc_file.keywords = 'example, PALM-4U'
nc_file.license = ''
nc_file.palm_version = ''
nc_file.references = ''
nc_file.source = 'PALM trunk'
nc_file.version = '1'

#%%% Mandatory global attributes
# ---------------------------
nc_file.Conventions = 'CF-1.7'
nc_file.origin_lat = 48.459638  # (overwrite initialization_parameters) Approx RIZ building
nc_file.origin_lon = 7.941967  # Used to initialize Coriolis parameter
nc_file.origin_time = '2022-02-22 10:00:00 +00'
nc_file.origin_x = 421779.645 # Source: https://coordinates-converter.com/en/decimal/48.459638,7.941967?karte=OpenStreetMap&zoom=19
nc_file.origin_y = 5367929.731
nc_file.origin_z = 0.0
nc_file.rotation_angle = 0.0

#%% Coordinates
nx = 259
ny = 145
nz = 50
dx = 1
dy = 1
dz = 1

prRed("netCDF main attributes done!")

#%% Create dimensions

#%%% X dimension
nc_file.createDimension('x', nx+1)
x = nc_file.createVariable('x', 'f4', ('x',))
x.long_name = 'distance to origin in x-direction'
x.units = 'm'
x.axis = 'X'
x[:] = np.arange(0, (nx+1)*dx, dx) + 0.5 * dx

#%%% Y dimension
nc_file.createDimension('y', ny+1)
y = nc_file.createVariable('y', 'f4', ('y',))
y.long_name = 'distance to origin in y-direction'
y.units = 'm'
y.axis = 'Y'
y[:] = np.arange(0, (ny+1)*dy, dy) + 0.5 * dy

#%%% Z dimension
# NOTE if your simulation uses a stretched vertical grid, you
# need to modify the z coordinates, e.g. z_array = (...)

z_array = np.append(0, np.arange(dz/2, (nz)*dz, dz))
nc_file.createDimension('z', nz+1)
z = nc_file.createVariable('z', 'f4', ('z',))
z.long_name = 'height above origin'
z.units = 'm'
z.axis = 'Z'
z.positive = 'up'
z[:] = z_array

#%%% nbuilding_pars dimension
nbuilding_pars_array = np.arange(150)
nc_file.createDimension('nbuilding_pars', len(nbuilding_pars_array))
nbuilding_pars = nc_file.createVariable('nbuilding_pars', 'i4', ('nbuilding_pars',))
nbuilding_pars[:] = nbuilding_pars_array


#%%% zlad dimension (height above ground)
zlad_array = z[:6]
nc_file.createDimension('zlad', len(zlad_array))
zlad = nc_file.createVariable('zlad', 'f4', ('zlad',))
zlad.long_name = 'height above ground'
zlad.units = 'm'
zlad.axis = 'Z'
zlad.positive = 'up'
zlad[:] = zlad_array

#%%% zsoil dimension 
dz_soil = np.array((0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86))
zsoil_fullLayers = np.zeros_like(dz_soil)
zsoil_fullLayers = np.around(
    [np.sum(dz_soil[:zs]) for zs in np.arange(1, len(dz_soil)+1)], 2)
zsoil_array = zsoil_fullLayers - dz_soil/2.

#############################################################################

nc_file.createDimension('zsoil', len(zsoil_array))
zsoil = nc_file.createVariable('zsoil', 'f4', ('zsoil',))
zsoil.long_name = 'depth in the soil'
zsoil.units = 'm'
zsoil.axis = 'Z'
zsoil.positive = 'down'
zsoil[:] = zsoil_array
#%%% nsurface_fraction dimension
nc_file.createDimension('nsurface_fraction', 3)
nsurface_fraction = nc_file.createVariable(
    'nsurface_fraction', 'i4', ('nsurface_fraction',))
nsurface_fraction[:] = np.arange(3)

prRed("Create dimensions done!")

#%% Create Variables
#Variables Topography set up

#%%% create building_id variable
nc_building_id = nc_file.createVariable(
    'building_id', 'i4', ('y', 'x'), fill_value=-9999)
nc_building_id.long_name = "building id number"
nc_building_id.units = "1"
nc_building_id[:, :] = nc_building_id._FillValue


#%%% create buildings_2d variable
nc_buildings_2d = nc_file.createVariable(
    'buildings_2d', 'f4', ('y', 'x'), fill_value=-9999.0)
nc_buildings_2d.long_name = "building height"
nc_buildings_2d.units = "m"
nc_buildings_2d.lod = np.int32(1)
nc_buildings_2d[:, :] = nc_buildings_2d._FillValue

#%%% create buildings_3d variable
nc_buildings_3d = nc_file.createVariable(
    'buildings_3d', 'i1', ('z', 'y', 'x'), fill_value=-127)
nc_buildings_3d.long_name = "building flag"
nc_buildings_3d.units = "1"
nc_buildings_3d.lod = np.int32(2)
nc_buildings_3d[:, :, :] = nc_buildings_3d._FillValue

#%%% create zt variable (terrain height)
nc_zt = nc_file.createVariable(
    'zt', 'f4', ('y', 'x'), fill_value=-9999.0)
nc_zt.long_name = 'terrain height'
nc_zt.units = 'm'
nc_zt[:, :] = nc_zt._FillValue



#%%% create building_type variable


nc_building_type = nc_file.createVariable(
            'building_type', 'i1', ('y', 'x'), fill_value=-127)
nc_building_type.long_name = "building type classification"
nc_building_type.units = "1"
nc_building_type[:, :] = nc_building_type._FillValue


#%%% create pavement_type variable

nc_pavement_type = nc_file.createVariable(
    'pavement_type', 'i1', ('y', 'x'), fill_value=-127)
nc_pavement_type.long_name = "pavement type classification"
nc_pavement_type.units = "1"
nc_pavement_type[:, :] = nc_pavement_type._FillValue


#%%% create soil_type variable

nc_soil_type = nc_file.createVariable(
    'soil_type', 'i1', ('y', 'x'), fill_value=-127)
nc_soil_type.long_name = "soil type classification"
nc_soil_type.units = "1"
nc_soil_type.lod = np.int32(1)
nc_soil_type[:, :] = nc_soil_type._FillValue

#%%% create street_type variable

nc_street_type = nc_file.createVariable(
    'street_type', 'i1', ('y', 'x'),
    fill_value=-127)
nc_street_type.long_name = "street type classification"
nc_street_type.units = "1"
nc_street_type[:, :] = nc_street_type._FillValue


#%%% create surface_fraction variable

nc_surface_fraction = nc_file.createVariable(
    'surface_fraction', 'f4', ('nsurface_fraction', 'y', 'x'), fill_value=-9999.0)
nc_surface_fraction.long_name = "surface fraction"
nc_surface_fraction.units = "1"
nc_surface_fraction[0, :, :] = nc_surface_fraction._FillValue  # vegetation fraction
nc_surface_fraction[1, :, :] = nc_surface_fraction._FillValue  # pavement fraction
nc_surface_fraction[2, :, :] = nc_surface_fraction._FillValue  # water fraction


#%%% create vegetation_type variable

nc_vegetation_type = nc_file.createVariable(
    'vegetation_type', 'i1', ('y', 'x'), fill_value=-127)
nc_vegetation_type.long_name = "vegetation type classification"
nc_vegetation_type.units = "1"
nc_vegetation_type[:, :] = nc_vegetation_type._FillValue



#%%% create water_type variable

nc_water_type = nc_file.createVariable(
    'water_type', 'i1', ('y', 'x'), fill_value=-127)
nc_water_type.long_name = "water type classification"
nc_water_type.units = "1"
nc_water_type[:, :] = nc_water_type._FillValue

#%%% create leaf area density variable

nc_lad = nc_file.createVariable(
    'lad', 'f4', ('zlad', 'y', 'x'), fill_value=-9999.0)
nc_lad.long_name = "leaf area density"
nc_lad.units = "m2 m-3"
nc_lad[:, :, :] = nc_lad._FillValue

#%%% create building_pars variable
nc_building_pars = nc_file.createVariable(
    'building_pars', 'f4', ('nbuilding_pars', 'y', 'x'), fill_value=-9999.0)
nc_building_pars.long_name = "building parameters"
nc_building_pars.units = ""
nc_building_pars.res_orig = 1.0
nc_building_pars.lod = "E_UTM N_UTM lon lat"
nc_building_pars[:, :, :] = nc_building_pars._FillValue

prRed("Create Variables done!")

#%% Populate Variables

# new code
img_n = Image.open(path1) 
imgg_n = img_n.convert(mode='RGB')
numpydata_n = asarray(imgg_n)

r_arr, g_arr, b_arr = numpydata_n[:,:,0], numpydata_n[:,:,1], numpydata_n[:,:,2]


# =============================================================================
# # Just to check for understanding  
# # shape
# print(numpydata.shape)
#   
# # Below is the way of creating Pillow
# # image from our numpyarray
# 
# pilImage = Image.fromarray(numpydata)
# print(type(pilImage))
# 
# # Let us check  image details
# # The array is like a mirror of the original image
# # you can check that from the following command
# print(pilImage.mode)
# print(pilImage.size)
# 
# =============================================================================


#%%% populate nc_building_id 
#                y      x
nc_building_id[50:94, 55:85] = 1 
nc_building_id[50:94, 55:65] = 2 
nc_building_id[50:94, 75:85] = 3 
# nc_building_id[-5:, -5:] = 4


#%%% populate nc_buildings_2d 

nc_buildings_2d[:, :] = np.where(
    nc_building_id[:, :] == 1,
    nc_building_id[:, :] * 15.0,
    nc_buildings_2d[:, :])

nc_buildings_2d[:, :] = np.where(
    nc_building_id[:, :] == 2,
    (nc_building_id[:, :] * 19.0)/2,
    nc_buildings_2d[:, :])

nc_buildings_2d[:, :] = np.where(
    nc_building_id[:, :] == 3,
    (nc_building_id[:, :] * 19.0)/3,
    nc_buildings_2d[:, :])


#%%% populate nc_buildings_3d 

for k in range(nz+1):
    nc_buildings_3d[k, :, :] = np.where(nc_buildings_2d[:, :] > z[k], 1, 0)

nc_building_pars[0, :, :] = np.where(nc_buildings_2d[:, :] > 0, 0.4, -9999.0) #Wall fraction above the ground floor level
nc_building_pars[1, :, :] = np.where(nc_buildings_2d[:, :] > 0, 0.6, -9999.0) #Window fraction above the ground floor level
    
nc_zt[:, :] = 4.0


#%%% populate nc_building_type 

# Set surface types
# -----------------
nc_building_type[:, :] = np.where(
nc_building_id[:, :] > nc_building_id._FillValue,
nc_building_id[:, :] * 1, 
nc_building_type[:, :])

#%%% populate nc_pavement_type 

pavement_arr = np.where(g_arr > 200,1, 0)
pavement_arr = np.ma.masked_where(pavement_arr == 1, pavement_arr)


pavement_arr_flip = np.flip(pavement_arr,0)

nc_pavement_type[:,:] = pavement_arr_flip

#%%% populate nc_vegetation_type 

vegetation_arr = np.where(g_arr > 200,1, 0)
vegetation_arr = np.ma.masked_where(vegetation_arr == 0, vegetation_arr)


vegetation_arr_flip = np.flip(vegetation_arr,0)
    
nc_vegetation_type[:,:] = vegetation_arr_flip


#%%% populate nc_water_type 

nc_water_type[17:29, 47:] = 2
nc_water_type[50:61, 47:] = 2
#%%% populate nc_soil_type 
# soil type is dependent on pavement
nc_soil_type[:, :] = np.where(
    nc_vegetation_type[:, :] > nc_vegetation_type._FillValue,
    2,
    nc_soil_type[:, :])
nc_soil_type[:, :] = np.where(
    nc_pavement_type[:, :] > nc_pavement_type._FillValue,
    2,
    nc_soil_type[:, :])
#%%% populate nc_street_type 

nc_street_type[:, :] = np.where(nc_pavement_type[:, :] == 1, 11, nc_street_type[:, :])
nc_street_type[:, :] = np.where(nc_pavement_type[:, :] == 2, 13, nc_street_type[:, :])
#%%% populate nc_surface_fraction 

nc_surface_fraction[0, :, :] = np.where(
    nc_building_id[:, :] > nc_building_id._FillValue,
    nc_surface_fraction[0, :, :],
    0)
nc_surface_fraction[2, :, :] = nc_surface_fraction[1, :, :] = nc_surface_fraction[0, :, :]
nc_surface_fraction[0, :, :] = np.where(
    nc_vegetation_type[:, :] > nc_vegetation_type._FillValue,
    1,
    nc_surface_fraction[0, :, :])
nc_surface_fraction[1, :, :] = np.where(
    nc_pavement_type[:, :] > nc_pavement_type._FillValue,
    1,
    nc_surface_fraction[1, :, :])
nc_surface_fraction[2, :, :] = np.where(
    nc_water_type[:, :] > nc_water_type._FillValue,
    1,
    nc_surface_fraction[2, :, :])


#%%% populate leaf area density 

# Set vegetation
        # --------------
lad_profile = [0.0, 0.01070122, 0.1070122, 0.3130108, 0.3879193, 0.1712195]
for k in range(len(zlad)):
    nc_lad[k, 38:43, 51:57] = lad_profile[k]
    # nc_lad[k, 11:13, 15:17] = lad_profile[k]
    

a=nc_buildings_3d[:,:, :]    
nc_file.close()
prYellow("File Should be saved in: {}".format(path))



#%% Plotting

ds = xr.open_dataset(path, engine="netcdf4")
ds["pavement_type"].plot()
print(ds.info())

ds.close()




















#%%




