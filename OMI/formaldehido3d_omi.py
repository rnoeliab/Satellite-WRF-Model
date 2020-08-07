#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 18:59:46 2020

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
output=INPUT+'figures/3d/'
files = sorted(glob(INPUT+'data/3d/'+'*.nc.nc4'))
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
#        print(f.keys())
        aa = []
        dset = f['key_science_data_column_amount']
        data =dset[:].astype(np.float64)
        # Read lat/lon data.
        latitude = f['latitude'][:]
        longitude = f['longitude'][:]
        
        lon, lat = np.meshgrid(longitude,latitude)
        # Handle fill value.
        lat[lat == -1.00000000e+30] = np.nan
        lon[lon == -1.00000000e+30] = np.nan
        data[data == -1.00000000e+30] = np.nan
        
        data[data < 0] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)
        # Subset data at nCandidate = 0
        data = data*(10**-16)
    
        import basemap_hcho as bm_hcho
        title = str('HCHO column total')
        unidad = str('$10^{16} moleculas/cm^{2}$')
        temp = []
        temp.append(FILE_NAME[-38:-29])
        bm_hcho.plot_map(output, data, lat, lon, temp, title, unidad, FILE_NAME)

