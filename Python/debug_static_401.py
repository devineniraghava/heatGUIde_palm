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

#%%create a netCDF file

nc_file = Dataset('testing_static', 'w', format='NETCDF4')
nc_file.title = 'Example PALM static driver'
nc_file.author = 'PALM user'
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

# Mandatory global attributes
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
nx = 39
ny = 39
nz = 40
dx = 2
dy = 2
dz = 2


nc_file.createDimension('x', nx+1)
x = nc_file.createVariable('x', 'f4', ('x',))
x.long_name = 'distance to origin in x-direction'
x.units = 'm'
x.axis = 'X'
x[:] = np.arange(0, (nx+1)*dx, dx) + 0.5 * dx

nc_file.createDimension('y', ny+1)
y = nc_file.createVariable('y', 'f4', ('y',))
y.long_name = 'distance to origin in y-direction'
y.units = 'm'
y.axis = 'Y'
y[:] = np.arange(0, (ny+1)*dy, dy) + 0.5 * dy


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


zlad_array = z[:6]
nc_file.createDimension('zlad', len(zlad_array))
zlad = nc_file.createVariable('zlad', 'f4', ('zlad',))
zlad.long_name = 'height above ground'
zlad.units = 'm'
zlad.axis = 'Z'
zlad.positive = 'up'
zlad[:] = zlad_array

nc_file.createDimension('nsurface_fraction', 3)
nsurface_fraction = nc_file.createVariable(
    'nsurface_fraction', 'i4', ('nsurface_fraction',))
nsurface_fraction[:] = np.arange(3)


#%% Variables Topography set up


nc_building_id = nc_file.createVariable(
    'building_id', 'i4', ('y', 'x'), fill_value=-9999)
nc_building_id.long_name = "building id number"
nc_building_id.units = "1"
nc_building_id[:, :] = nc_building_id._FillValue

nc_buildings_2d = nc_file.createVariable(
    'buildings_2d', 'f4', ('y', 'x'), fill_value=-9999.0)
nc_buildings_2d.long_name = "building height"
nc_buildings_2d.units = "m"
nc_buildings_2d.lod = np.int32(1)
nc_buildings_2d[:, :] = nc_buildings_2d._FillValue

nc_buildings_3d = nc_file.createVariable(
    'buildings_3d', 'i1', ('z', 'y', 'x'), fill_value=-127)
nc_buildings_3d.long_name = "building flag"
nc_buildings_3d.units = "1"
nc_buildings_3d.lod = np.int32(2)
nc_buildings_3d[:, :, :] = nc_buildings_3d._FillValue

nc_zt = nc_file.createVariable(
    'zt', 'f4', ('y', 'x'), fill_value=-9999.0)
nc_zt.long_name = 'terrain height'
nc_zt.units = 'm'
nc_zt[:, :] = nc_zt._FillValue



############


nc_building_type = nc_file.createVariable(
            'building_type', 'i1', ('y', 'x'), fill_value=-127)
nc_building_type.long_name = "building type classification"
nc_building_type.units = "1"
nc_building_type[:, :] = nc_building_type._FillValue

nc_pavement_type = nc_file.createVariable(
    'pavement_type', 'i1', ('y', 'x'), fill_value=-127)
nc_pavement_type.long_name = "pavement type classification"
nc_pavement_type.units = "1"
nc_pavement_type[:, :] = nc_pavement_type._FillValue

nc_soil_type = nc_file.createVariable(
    'soil_type', 'i1', ('y', 'x'), fill_value=-127)
nc_soil_type.long_name = "soil type classification"
nc_soil_type.units = "1"
nc_soil_type.lod = np.int32(1)
nc_soil_type[:, :] = nc_soil_type._FillValue

####################

nc_street_type = nc_file.createVariable(
    'street_type', 'i1', ('y', 'x'),
    fill_value=-127)
nc_street_type.long_name = "street type classification"
nc_street_type.units = "1"
nc_street_type[:, :] = nc_street_type._FillValue

nc_surface_fraction = nc_file.createVariable(
    'surface_fraction', 'f4', ('nsurface_fraction', 'y', 'x'), fill_value=-9999.0)
nc_surface_fraction.long_name = "surface fraction"
nc_surface_fraction.units = "1"
nc_surface_fraction[0, :, :] = nc_surface_fraction._FillValue  # vegetation fraction
nc_surface_fraction[1, :, :] = nc_surface_fraction._FillValue  # pavement fraction
nc_surface_fraction[2, :, :] = nc_surface_fraction._FillValue  # water fraction

nc_vegetation_type = nc_file.createVariable(
    'vegetation_type', 'i1', ('y', 'x'), fill_value=-127)
nc_vegetation_type.long_name = "vegetation type classification"
nc_vegetation_type.units = "1"
nc_vegetation_type[:, :] = nc_vegetation_type._FillValue

nc_water_type = nc_file.createVariable(
    'water_type', 'i1', ('y', 'x'), fill_value=-127)
nc_water_type.long_name = "water type classification"
nc_water_type.units = "1"
nc_water_type[:, :] = nc_water_type._FillValue

#####################

nc_lad = nc_file.createVariable(
    'lad', 'f4', ('zlad', 'y', 'x'), fill_value=-9999.0)
nc_lad.long_name = "leaf area density"
nc_lad.units = "m2 m-3"
nc_lad[:, :, :] = nc_lad._FillValue

##################

nc_building_id[:5, :5] = 1
nc_building_id[-5:, :5] = 2
nc_building_id[:5, -5:] = 3
nc_building_id[-5:, -5:] = 4

nc_buildings_2d[:, :] = np.where(
    nc_building_id[:, :] > nc_building_id._FillValue,
    nc_building_id[:, :] * 10.0,
    nc_buildings_2d[:, :])

for k in range(nz+1):
    nc_buildings_3d[k, :, :] = np.where(nc_buildings_2d[:, :] > z[k], 1, 0)

nc_zt[:, :] = 4.0

# Set surface types
# -----------------
nc_building_type[:, :] = np.where(
nc_building_id[:, :] > nc_building_id._FillValue,
nc_building_id[:, :] * 1,
nc_building_type[:, :])


nc_pavement_type[9:11, :5] = 1
nc_pavement_type[:, 7:13] = 2

nc_vegetation_type[:, 5:7] = 3
nc_vegetation_type[:, 13:15] = 3
nc_vegetation_type[5:15, 15:] = 3
nc_vegetation_type[7:9, 15:17] = 1
nc_vegetation_type[11:13, 15:17] = 1

nc_water_type[5:9, :5] = 2
nc_water_type[11:15, :5] = 2

nc_soil_type[:, :] = np.where(
    nc_vegetation_type[:, :] > nc_vegetation_type._FillValue,
    2,
    nc_soil_type[:, :])
nc_soil_type[:, :] = np.where(
    nc_pavement_type[:, :] > nc_pavement_type._FillValue,
    2,
    nc_soil_type[:, :])

nc_street_type[:, :] = np.where(nc_pavement_type[:, :] == 1, 11, nc_street_type[:, :])
nc_street_type[:, :] = np.where(nc_pavement_type[:, :] == 2, 13, nc_street_type[:, :])

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

# Set vegetation
        # --------------
lad_profile = [0.0, 0.01070122, 0.1070122, 0.3130108, 0.3879193, 0.1712195]
for k in range(len(zlad)):
    nc_lad[k, 7:9, 15:17] = lad_profile[k]
    nc_lad[k, 11:13, 15:17] = lad_profile[k]

nc_file.close()












#%%






