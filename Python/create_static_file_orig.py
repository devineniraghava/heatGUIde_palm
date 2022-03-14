#------------------------------------------------------------------------------#
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2020 Leibniz Universitaet Hannover
#------------------------------------------------------------------------------#

import math
import netCDF4
from   netCDF4 import Dataset
import numpy as np
import datetime


class StaticDriver:
    """ This is an example script to generate static drivers for PALM. 
    
    You can use it as a starting point for creating your setup specific
    driver.  
    """
    
    
    def __init__(self):
        """ Open the static driver as NetCDF4 file. Here, you have to give the 
        full path to the static driver that shall be created. Existing file 
        with same name is deleted. 
        """
        print('Opening file...')
        self.nc_file = Dataset("/home/rdevinen/Documents/PALM_Seminar/NetCDF/Python/basic_static_driver_1.nc", 'w', format='NETCDF4')


    def write_global_attributes(self):
        """ Write global attributes to static driver. """
        print("Writing global attributes...")
        
        # mandatory global attributes
        self.nc_file.origin_lon = 55.0 # used to initialize coriolis parameter
        self.nc_file.origin_lat = 0.0 # (overwrite initialization_parameters)
        self.nc_file.origin_time = '2018-07-08 00:00:00 +00'
        self.nc_file.origin_x = 308124
        self.nc_file.origin_y = 6098908
        self.nc_file.origin_z = 0.0
        self.nc_file.rotation_angle = 0.0
        
        # optional global attributes
        self.nc_file.author = 'Matthias Suehring'
        self.nc_file.comment = 'Miscellaneous information about the data ' \
                               'or methods to produce it.'
        self.nc_file.creation_date = str(datetime.datetime.now())
        self.nc_file.institution = 'Pecanode'
        self.nc_file.history = ''
        self.nc_file.palm_revision = ''
        self.nc_file.title = 'Static driver for performance test'


    def define_dimensions(self):
        """ Set dimensions on which variables are defined. """
        print("Writing dimensions...")
        
        # specify general grid parameters
        # these values must equal those set in the initialization_parameters
        self.nx = 199
        self.ny = 199
        self.nz = 20
        self.dx = 10
        self.dy = 10
        self.dz = 5
        
        # create soil grid (only relevant if land surface module is used)
        dz_soil = np.array((0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86))
        zsoil_fullLayers = np.zeros_like(dz_soil)
        zsoil_fullLayers = np.around([np.sum(dz_soil[:zs]) for zs in np.arange(1,len(dz_soil)+1)],2)
        zsoil_array = zsoil_fullLayers - dz_soil/2.
        
        # create plant canopy grid (only relevant of plant canopy module is used)
        zlad_array = np.arange(0,8*self.dz,self.dz) #np.arange(0,10*self.dz,self.dz)
        zlad_array[1:len(zlad_array)] = zlad_array[1:len(zlad_array)] - 0.5*self.dz
        
        # mandatory dimensions
        self.nc_file.createDimension('x' ,self.nx+1)
        self.x = self.nc_file.createVariable('x', 'f8', ('x',))
        self.x.units = 'm'
        self.x.standard_name = 'x coordinate of cell centers'
        self.x[:] = np.arange(0,(self.nx+1)*self.dx,self.dx)
        
        self.nc_file.createDimension('xu', self.nx)
        self.xu = self.nc_file.createVariable('xu', 'f8', ('xu',))
        self.xu.units = 'm'
        self.xu.standard_name = 'x coordinate of cell edges'
        self.xu[:] = np.arange(self.dx/2.,self.nx*self.dx,self.dx)
        
        self.nc_file.createDimension('y', self.ny+1)
        self.y = self.nc_file.createVariable('y', 'f8', ('y',))
        self.y.units = 'm'
        self.y.standard_name = 'y coordinate of cell centers'
        self.y[:] = np.arange(0,(self.ny+1)*self.dy,self.dy)
        
        self.nc_file.createDimension('yv', self.ny)
        self.yv = self.nc_file.createVariable('yv', 'f8', ('yv',))
        self.yv.units = 'm'
        self.yv.standard_name = 'y coordinate of cell edges'
        self.yv[:] = np.arange(self.dy/2.,self.ny*self.dy,self.dy)
        
        # if your simulation uses a stretched vertical grid, you need to 
        # modify the z and zw coordinates, 
        # e.g. z_array = (...) and zw_array = (...)
        z_array = np.append(0, np.arange(self.dz/2,(self.nz+1)*self.dz,self.dz))
        self.nc_file.createDimension('z',self.nz+2)
        self.z = self.nc_file.createVariable('z', 'f8', ('z',))
        self.z.units = 'm'
        self.z.standard_name = 'z coordinate of cell centers'
        self.z[:] = z_array

        zw_array = np.arange(0,(self.nz+2)*self.dz,self.dz)
        self.nc_file.createDimension('zw', self.nz+2)
        self.zw = self.nc_file.createVariable('zw', 'f8', ('zw',))
        self.zw.units = 'm'
        self.zw.standard_name = 'z coordinate of cell edges'
        self.zw[:] = zw_array
        
        # optional dimensions, uncomment if needed
        self.nc_file.createDimension('zsoil', len(zsoil_array))
        self.zsoil = self.nc_file.createVariable('zsoil', 'f8', ('zsoil',))
        self.zsoil.positive = 'down'
        self.zsoil.units = 'm'
        self.zsoil.standard_name = 'depth below land'
        self.zsoil[:] = zsoil_array
        
        self.nc_file.createDimension('zlad', len(zlad_array))
        self.zlad = self.nc_file.createVariable('zlad', 'f8', ('zlad',))
        self.zlad.units = 'm'
        self.zlad.standard_name = 'z coordinate of resolved plant canopy'
        self.zlad[:] = zlad_array
        

        self.nc_file.createDimension('nsurface_fraction', 3)
        self.nsurface_fraction = self.nc_file.createVariable(
            'nsurface_fraction', 'i1', ('nsurface_fraction',))
        self.nsurface_fraction[:] = np.arange(0,3)
        

    def add_variables(self):

        print("Writing variables...")
        fill_val_real = -9999.9
        fill_val_int8 = -127
        
        # topography variables
        nc_buildings_2d = self.nc_file.createVariable(
             'buildings_2d', 'f4', ('y','x'), fill_value=fill_val_real)
        nc_buildings_2d.lod = 1
        nc_buildings_2d[:,:] = fill_val_real
        
        nc_building_id = self.nc_file.createVariable(
            'building_id', 'i1', ('y','x'), fill_value=fill_val_int8)
        nc_building_id[:,:] = fill_val_int8

        nc_building_type = self.nc_file.createVariable(
            'building_type', 'i1', ('y','x'), fill_value=fill_val_int8)
        nc_building_type[:,:] = fill_val_int8


        
        ## terrain height
        nc_zt = self.nc_file.createVariable(    \
            'zt', 'f4', ('y','x'), fill_value=-fill_val_real)
        nc_zt.long_name = 'orography'
        nc_zt.units = 'm'

        ## resolved plant canopy
        nc_lad = self.nc_file.createVariable(
            'lad', 'f4', ('zlad','y','x'), fill_value=fill_val_real)
        nc_lad.long_name = 'leaf area density'
        nc_lad[:,:,:] = fill_val_real

        # surface types
        nc_pavement_type = self.nc_file.createVariable(
            'pavement_type', 'i1', ('y','x'), fill_value=fill_val_int8)
        nc_pavement_type[:,:] = fill_val_int8

        nc_vegetation_type = self.nc_file.createVariable( \
            'vegetation_type', 'i1', ('y','x'), fill_value=fill_val_int8)
        nc_vegetation_type[:,:] = fill_val_int8

        nc_water_type = self.nc_file.createVariable(
            'water_type', 'i1',('y','x'), fill_value=fill_val_int8)
        nc_water_type[:,:] = fill_val_int8

        nc_soil_type = self.nc_file.createVariable(
            'soil_type', 'i1', ('y','x'), fill_value=fill_val_int8)
        nc_soil_type.lod = 1 # if set to 1, adjust dimensions of soil_type

        nc_surface_fraction = self.nc_file.createVariable(
            'surface_fraction', 'f', ('nsurface_fraction','y','x'), fill_value=fill_val_real)

        # Create terrain
        nc_zt[:,:] = 0.0

        #
#       Create canopy profile
        lad_array = np.zeros(len(self.zlad ))
        lad_array = [ 0.1, 0.2, 0.4, 0.4, 0.4, 0.3, 0.2, 0.1 ]

        # Set index ranges for buildings and courtyards
        build_inds_start_x = [ 5,  25, 45, 65, 85, 105, 125, 145, 165, 185 ]
        build_inds_end_x   = [ 15, 35, 55, 75, 95, 115, 135, 155, 175, 195 ]

        court_inds_start_x = [ 8,  28, 48, 68, 88, 108, 128, 148, 168, 188 ]
        court_inds_end_x   = [ 12, 32, 52, 72, 92, 112, 132, 152, 172, 192 ]

        build_inds_start_y = build_inds_start_x
        build_inds_end_y   = build_inds_end_x

        court_inds_start_y = court_inds_start_x
        court_inds_end_y   = court_inds_end_x

        # Set index ranges for streets and tree lines
        pave_inds_start_x = [ 0, 18, 38, 58, 78, 98,  118, 138, 158, 178, 198 ]
        pave_inds_end_x   = [ 2, 22, 42, 62, 82, 102, 122, 142, 162, 182, 199 ]
        pave_inds_start_y = pave_inds_start_x
        pave_inds_end_y = pave_inds_end_x


        for i in range(0,self.nx+1):
           for j in range(0,self.ny+1):
              # Pre-set vegetation and soil type, as well as fractios. Note, later on these
              # are partly re-set.
              nc_vegetation_type[j,i] = 3
              nc_soil_type[j,i] = 1

              nc_surface_fraction[0,j,i] = 1.0
              nc_surface_fraction[1,j,i] = 0.0
              nc_surface_fraction[2,j,i] = 0.0

              # Set topography
              nc_zt[j,i] = 25.0 * ( np.sin(2.0 * np.pi * (i-0.5)*self.dx / (self.dx * self.nx+1))**2 \
                                  + np.cos(2.0 * np.pi * (j-0.5)*self.dy / (self.dy * self.ny+1))**2 )

              # Set building height
              for ii in range(0,len(build_inds_start_x)):
                 for jj in range(0,len(build_inds_start_y)):

                    if ( i >= build_inds_start_x[ii] and i <= build_inds_end_x[ii] ) and \
                       ( j >= build_inds_start_y[jj] and j <= build_inds_end_y[jj] ) and \
                       not ( ( i >= court_inds_start_x[ii] and i <= court_inds_end_x[ii] ) and \
                             ( j >= court_inds_start_y[jj] and j <= court_inds_end_y[jj] ) ) :

                       nc_buildings_2d[j,i] = 50.0
                       nc_building_type[j,i] = 2

                       # For each (ii,jj)-pair create an unique ID
                       nc_building_id[j,i] = ii + len(build_inds_start_x) * jj

                       nc_vegetation_type[j,i] = fill_val_int8
                       nc_soil_type[j,i] = fill_val_int8
                       nc_surface_fraction[0:2,j,i] = fill_val_real

                    # Set LAD in courtyard area
                    if ( i >= court_inds_start_x[ii] and i <= court_inds_end_x[ii] ) and \
                       ( j >= court_inds_start_y[jj] and j <= court_inds_end_y[jj] ):
                       nc_lad[:,j,i] = lad_array[:]

              # Set pavement and tree lines
              for ii in range(0,len(pave_inds_start_x)):
                 for jj in range(0,len(pave_inds_start_y)):
                    if ( i >= pave_inds_start_x[ii] and i <= pave_inds_end_y[ii] ) or \
                       ( j >= pave_inds_start_y[jj] and j <= pave_inds_end_y[jj] ):
                       nc_lad[:,j,i] = lad_array[:]

                       nc_pavement_type[j,i] = 2
                       nc_vegetation_type[j,i] = fill_val_int8
                       nc_surface_fraction[0,j,i] = 0.0
                       nc_surface_fraction[1,j,i] = 1.0
                       nc_surface_fraction[2,j,i] = 0.0

    def finalize(self):
        """ Close file """
        print("Closing file...")
        
        self.nc_file.close()


if __name__ == '__main__':
    driver = StaticDriver()
    driver.write_global_attributes()
    driver.define_dimensions()
    driver.add_variables()
    driver.finalize()
    
