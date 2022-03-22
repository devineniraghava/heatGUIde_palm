#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 11:32:36 2022

@author: rdevinen
"""
import datetime
import math
import pytz

from netCDF4 import Dataset
import numpy as np

from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')


#%% netCDF main attributes

#%%% Optional attribbutes
nc_file = Dataset('/home/rdevinen/palm/current_version/JOBS/practice/INPUT/practice_static', 'w', format='NETCDF4')
nc_file.title = 'Example PALM static driver'
nc_file.author = 'devineni'
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
nc_file.origin_lat = 48.459638  # (overwrite initialization_parameters) 
nc_file.origin_lon = 7.941967  # Used to initialize Coriolis parameter
nc_file.origin_time = '2022-02-22 10:00:00 +00'
nc_file.origin_x = 421779.645 # Source: https://coordinates-converter.com/en/decimal/48.459638,7.941967?karte=OpenStreetMap&zoom=19
nc_file.origin_y = 5367929.731
nc_file.origin_z = 0.0
nc_file.rotation_angle = 0.0

#%% Coordinates
nx = 63
ny = 71
nz = 40
dx = 1
dy = 1
dz = 1


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

#%% Populate Variables

#%%% populate nc_building_id 
#                y      x
nc_building_id[17:61, 17:47] = 1
nc_building_id[17:61, 17:27] = 2
nc_building_id[17:61, 37:47] = 3
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

nc_zt[:, :] = 4.0


#%%% populate nc_building_type 

# Set surface types
# -----------------
nc_building_type[:, :] = np.where(
nc_building_id[:, :] > nc_building_id._FillValue,
nc_building_id[:, :] * 1, 
nc_building_type[:, :])

#%%% populate nc_pavement_type 

nc_pavement_type[17:61, 9:17] = 1
nc_pavement_type[8:17, 9:] = 2


#%%% populate nc_vegetation_type 

nc_vegetation_type[:, :9] = 3
nc_vegetation_type[61:, 9:] = 1
nc_vegetation_type[:8, 9:] = 1
nc_vegetation_type[29:50, 47:] = 3


# nc_vegetation_type[5:15, 15:] = 3
# nc_vegetation_type[7:9, 15:17] = 1
# nc_vegetation_type[11:13, 15:17] = 1
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
    nc_lad[k, 7:9, 15:17] = lad_profile[k]
    nc_lad[k, 11:13, 15:17] = lad_profile[k]
    

a=nc_buildings_3d[:,:, :]    
nc_file.close()












#%%






