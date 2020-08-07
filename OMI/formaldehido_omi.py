#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:25:53 2020

@author: noelia
"""

import numpy as np
import h5py
import time
import calendar
from glob import glob

#open several data
#inputt='/home/noelia/Downloads/'

INPUT = '/data/noelia/imagen_data/omi/OMHCHO_ori/2017/'
output=INPUT+'figures/2G/'

files = sorted(glob(INPUT+'data/2G/'+'*.he5'))
#### camino para obtener datos del NO2  #########
###################################################
#path = 'HDFEOS/SWATHS/OMI Total Column Amount HCHO/'
# LAT_NAME = path + 'Geolocation Fields/Latitude'
# LON_NAME = path + 'Geolocation Fields/Longitude'
# TIME = path + 'Geolocation Fields/Time'
path = 'HDFEOS/GRIDS/OMI Total Column Amount HCHO/Data Fields/'
DATAFIELD_NAME = path + 'ColumnAmountHCHO'
LAT_NAME = path + 'Latitude'
LON_NAME = path + 'Longitude'
TIME = path + 'Time'

####################################################
### leer todas las imagenes ###############
for FILE_NAME in files:
    FILE_NAME=FILE_NAME.strip()
    with h5py.File(FILE_NAME, mode='r') as f:
#        print(f[path+'/Data Fields'].keys())
        aa = []
        dset = f[DATAFIELD_NAME]
        data =dset[:].astype(np.float64)
        # Read lat/lon data.
        lat = f[LAT_NAME][:]
        lon = f[LON_NAME][:]
        times = f[TIME][:]

        # Get attributes needed for the plot.
        # String attributes actually come in as the bytes type and should
        # be decoded to UTF-8 (python3).
        title = dset.attrs['Title'].decode()
        units = dset.attrs['Units'].decode()
        
        missing_value = dset.attrs['MissingValue']
        _FillValue = dset.attrs['_FillValue']
        add_offset = dset.attrs['Offset']
        scale_factor = dset.attrs['ScaleFactor']
        UniqueField = dset.attrs['UniqueFieldDefinition']
        
        # Handle fill value.
        lat[lat == _FillValue] = np.nan
        lon[lon == _FillValue] = np.nan
            
        times[times == -1.00000000e+30]= np.nan
        data[data == _FillValue] = np.nan
        
        data = data * scale_factor + add_offset
        data[data < 0] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)
        # Subset data at nCandidate = 0
        lat = lat[0,:,:]
        lon = lon[0,:,:]
        data = data[0,:,:]*(10**-16)
        times = times[0,:,:]
        temp = []
        temp_min=time.localtime(np.nanmin(times)+calendar.timegm(time.strptime('Dec 31, 1992 @ 23:59:59 UTC','%b %d, %Y @ %H:%M:%S UTC')))
        temp.append(time.strftime('%d-%m-%Y %H:%M', temp_min))
        temp_max=time.localtime(np.nanmax(times)+calendar.timegm(time.strptime('Dec 31, 1992 @ 23:59:59 UTC','%b %d, %Y @ %H:%M:%S UTC')))
        temp.append(time.strftime('%d-%m-%Y %H:%M', temp_max))

        import basemap_hcho as bm_hcho
        unidad = str('$10^{16} moleculas/cm^{2}$')
        bm_hcho.plot_map(output, data, lat, lon, temp, title, unidad, FILE_NAME)



