#!/usr/bin/env python
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import numpy as np
from osgeo import osr
from osgeo import gdal
import time as t

# Define KM_PER_DEGREE
KM_PER_DEGREE = 111.32

# GOES-16 Spatial Reference System
sourcePrj = osr.SpatialReference()
sourcePrj.ImportFromProj4('+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027489803406172 +lat_0=0.0 +lon_0=-89.5 +sweep=x +no_defs')

# Lat/lon WSG84 Spatial Reference System
targetPrj = osr.SpatialReference()
targetPrj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

#def exportImage(image,path):
#    driver = gdal.GetDriverByName('netCDF')
#    return driver.CreateCopy(path,image,0)

def getGeoT(extent, nlines, ncols):
    # Compute resolution based on data dimension
    resx = (extent[2] - extent[0]) / ncols
    resy = (extent[3] - extent[1]) / nlines
    return [extent[0], resx, 0, extent[3] , 0, -resy]

def getScaleOffset (resolution): #(path):
    #values are based on GOES-R products documentation for ABI Full disk GOES-EAST
    #if these values are in your files of interest then load them in from files instead
    offset_y = [0.151865, 0.151858, 0.151844, 0.151816, 0.151900]
    offset_x = [-0.151865, -0.151858, -0.151844, -0.151816, -0.151900]
    scale_factor_y = [-0.000014, -0.000028, -0.000056, -0.0000112, -0.000280]
    scale_factor_x = [0.000014, 0.000028, 0.000056, 0.0000112, 0.000280]
    if resolution ==0.5 : 
        i=0
    elif resolution == 1 :
        i=1
    elif resolution == 2 :
        i=2
    elif resolution == 4 :
        i=3
    elif resolution == 10 :
        i=4
    else: 
        print('Offset and Scale factor are not known for this resolution.')
    offset = offset_x[i]
    scale = scale_factor_x[i]
    
    #file_ob = fs.open(path)
    #nc = xr.open_dataset(file_ob)
    #scale = nc.variables['CMI'].scale_factor
    #offset = nc.variables['CMI'].add_offset
    #nc = Dataset(path, mode='r')
    #nc.close()
    return scale, offset
    
def remap(path, extent, resolution, x1, y1, x2, y2): 
    
    # GOES-16 Extent (satellite projection) [llx, lly, urx, ury]
    GOES16_EXTENT = [x1, y1, x2, y2]
    
    # Setup NetCDF driver
    gdal.SetConfigOption('GDAL_NETCDF_BOTTOMUP', 'NO')
        
    # Read scale/offset from file
    scale, offset = getScaleOffset(resolution) #(path)
         
    try:  
        connectionInfo = 'NETCDF:\"' + path + '\":CMI'    
        # Open NetCDF file (GOES-16 data)  
        raw = gdal.Open(connectionInfo)
    except:
        connectionInfo = 'HDF5:\"' + path + '\"://CMI'    
        # Open NetCDF file (GOES-16 data)  
        raw = gdal.Open(connectionInfo)    
                
    # Setup projection and geo-transformation
    raw.SetProjection(sourcePrj.ExportToWkt())
    #raw.SetGeoTransform(getGeoT(GOES16_EXTENT, raw.RasterYSize, raw.RasterXSize))
    raw.SetGeoTransform(getGeoT(GOES16_EXTENT, raw.RasterYSize, raw.RasterXSize))  
  
    #print (KM_PER_DEGREE)
    # Compute grid dimension
    sizex = int(((extent[2] - extent[0]) * KM_PER_DEGREE) / resolution)
    sizey = int(((extent[3] - extent[1]) * KM_PER_DEGREE) / resolution)
    
    # Get memory driver
    memDriver = gdal.GetDriverByName('MEM')
   
    # Create grid
    grid = memDriver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)
    
    # Setup projection and geo-transformation
    grid.SetProjection(targetPrj.ExportToWkt())
    grid.SetGeoTransform(getGeoT(extent, grid.RasterYSize, grid.RasterXSize))

    # Perform the projection/resampling 

    print ('Remapping', path)
        
    start = t.time()
    
    gdal.ReprojectImage(raw, grid, sourcePrj.ExportToWkt(), targetPrj.ExportToWkt(), gdal.GRA_NearestNeighbour, options=['NUM_THREADS=ALL_CPUS']) 
    
    print ('- finished! Time:', t.time() - start, 'seconds')
    
    # Close file
    raw = None
        
    # Read grid data
    array = grid.ReadAsArray()
    
    # Mask fill values (i.e. invalid values)
    np.ma.masked_where(array, array == -1, False)
    
    # Apply scale and offset
    array = array * scale + offset
    
    grid.GetRasterBand(1).SetNoDataValue(-1)
    grid.GetRasterBand(1).WriteArray(array)

    return grid
