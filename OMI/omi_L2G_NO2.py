#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 13:44:09 2019

@author: noelia
"""

##################################################################
##   autor: rnoeliab                                             #
#    fecha: March - 2019                                         #
#    funciona: codigo para leer varias imagenes OMI con formato  #
#    HDFEOS5 con nivel 2LG para el MO2  (2014)                   #
##################################################################

import numpy as np
import h5py
import time
import calendar
from glob import glob

#open several data
inputt='/data/noelia/imagen_data/omi/OMNO2G_ori/'
output=inputt+'figures/'

files = sorted(glob(inputt+'data/'+'*.he5'))
#### camino para obtener datos del NO2  #########
###################################################
path = 'HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/'
DATAFIELD_NAME = path + 'ColumnAmountNO2Trop/'
CLOUD_FRACTION = path + 'CloudFraction'
LAT_NAME = path + 'Latitude'
LON_NAME = path + 'Longitude'
QUALIY = path + 'MeasurementQualityFlags'
TIME = path + 'Time'
ORBIT_NUMBER = path + 'OrbitNumber'
####################################################
### leer todas las imagenes ###############
for FILE_NAME in files:
    FILE_NAME=FILE_NAME.strip()
    with h5py.File(FILE_NAME, mode='r') as f:
        aa = []
        dset = f[DATAFIELD_NAME]
        data =dset[:].astype(np.float64)
        # Read lat/lon data.
        lat = f[LAT_NAME][:]
        lon = f[LON_NAME][:]
        times = f[TIME][:]
        on = f[ORBIT_NUMBER][:]     
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
            
        times[times == _FillValue] = np.nan
        data[data == _FillValue] = np.nan
        
        data = data * scale_factor + add_offset
        data[data < 0] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)
        # Subset data at nCandidate = 0
        data = (data[0,:,:])*(10**-16)
        lon = lon[0,:,:]
        lat = lat[0,:,:]
        tiempo = times[0,:,:]
        temp = []
        temp_min=time.localtime(np.nanmin(tiempo)+calendar.timegm(time.strptime('Dec 31, 1992 @ 23:59:59 UTC','%b %d, %Y @ %H:%M:%S UTC')))
        temp.append(time.strftime('%d-%m-%Y %H:%M', temp_min))
        temp_max=time.localtime(np.nanmax(tiempo)+calendar.timegm(time.strptime('Dec 31, 1992 @ 23:59:59 UTC','%b %d, %Y @ %H:%M:%S UTC')))
        temp.append(time.strftime('%d-%m-%Y %H:%M', temp_max))

        import basemap_no2 as bm_no2
        unidad = str('$10^{16} moleculas/cm^{2}$')
        bm_no2.plot_map(output, data, lat, lon, temp, title, unidad, FILE_NAME)





        
