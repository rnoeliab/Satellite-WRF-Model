#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 10:23:06 2020

@author: noelia
"""
########
# conda create -n gdal_test python=3.5
# activate gdal_test
# conda install gdal
#

import os
import gdal
import numpy as np
from PIL import Image
import netCDF4 as nc4
from datetime import datetime
from glob import glob

def find_nearest_x(longitude, point_x):
    longitude = np.asarray(longitude)
    idx = (np.abs(longitude - point_x)).argmin()
    return idx


def find_nearest_y(latitude, point_y):
    latitude = np.asarray(latitude)
    dist = (np.abs(latitude - point_y))
    idy = np.where(dist ==dist.min())[0][1]
    return idy


def create_nc(path,file_names,lat,lon,data,metadata):
    f = nc4.Dataset(OUT_PATH+file_names+".nc","w",format='NETCDF4')
    f.description = "Data extraida del MOD04 y MY04"
    database = f.createGroup('AOD_3K')
    database.createDimension('lon',lon.shape[1])
    database.createDimension('lat',lat.shape[0])
    database.createDimension('time', None)
    
    ## Building variables
    longitude = database.createVariable('Longitude','f4',('lon',))
    latitude = database.createVariable('Latitude','f4',('lat',))
    aod = database.createVariable('Optical_Depth_Land_And_Ocean','f4',('lat','lon'))
    time = database.createVariable('Time','i4','time')
  
    longitude[:] = lon[0,:] #The "[:]" at the end of the variable instance is necessary
    latitude[:] = lat[:,0]
    aod[:,:] = data
    
    date_time_str = str(metadata['RANGEENDINGDATE'] + " " + metadata['RANGEBEGINNINGTIME'])
    date_time_obj = datetime.strptime(date_time_str, '%Y-%m-%d %H:%M:%S.%f').toordinal()
    time[0] = date_time_obj
    f.history = "Created" + " " + str(date_time_str)
    
    longitude.units = "degrees"
    latitude.units = "degrees"
    time.units = "day since Jan 01, 0001"
    aod.units = " "
    
    f.close()

#open several data
INPUT_PATH = '/data/noelia/imagen_data/modis/DATA/SP/3K/2017/'
OUT_PATH = '/data/noelia/imagen_data/modis/DATA/SP/3K_nc/2017/'

mylist=[]
listdirr = []

files = sorted(glob(INPUT_PATH+"*.hdf"))  ## archivos wrfout
for i in range(len(files)):
    if i == (len(files)-1):
        break;      
    else:
        if files[i][-44:-27] == files[i+1][-44:-27]:
            listdirr.append(files[i])          
            if i == (len(files)-2):
                listdirr.append(files[len(files)-1])
                mylist.append(listdirr)               
        else:
            listdirr.append(files[i])
            mylist.append(listdirr)
            listdirr = []

################### AREA STUDY
x0 = -49.0
x1 = -44.0
y0 = -21.0
y1 = -26.0
################    

#var_name = 'AOD_550_Dark_Target_Deep_Blue_Combined'   ## (1)
#var_name = 'Optical_Depth_Land_And_Ocean'              ## (2)   
#var_name = 'Deep_Blue_Aerosol_Optical_Depth_550_Land'  ## (3)

var_name = 'mod04:Optical_Depth_Land_And_Ocean'
for n in range(len(mylist)):
#    print n    
#    count=0
    for FILE_NAME in mylist[n]:
        print(FILE_NAME)  
        os.system('gdalwarp -of GTIFF -tps -t_srs EPSG:4326 HDF4_EOS:EOS_SWATH:"{0}":{1} teste.hdf'.format(FILE_NAME, var_name))    
        OUT = os.getcwd()
        print(OUT)
        ifile = OUT+os.sep+'teste.hdf' 
        gc = gdal.Open(ifile)
        metadata = gc.GetMetadata()
        width = gc.RasterXSize
        height = gc.RasterYSize
        if width <= 1000 and height <=1000:
            print(width, height)
            gt = gc.GetGeoTransform()
            bands = gc.RasterCount
            projInfo = gc.GetProjection()
            valid_range = [-100, 5000]
            fv = -9999
            scale_factor=0.00100000004749745
            minx = gt[0]
            maxx = gt[0] + width*gt[1] + height*gt[2]
            miny = gt[3] + width*gt[4] + height*gt[5]
            maxy = gt[3]
            print ("la latitude esta en el rango:" + str(miny) + ":" + str(maxy))
            print ("la longitude esta en el rango:"+ str(minx) + ":" + str(maxx))
        
            lat = np.linspace(miny, maxy, height)[::-1]
            lon = np.linspace(minx, maxx, width)
            
            x,y = np.meshgrid(lon,lat)
            im1 = Image.open(OUT+os.sep+'teste.hdf')
            im1 = np.array(im1)
            data = np.float64(im1)            
            data[data == fv] = np.nan
            data[data <= 0.0] = np.nan
            data = np.ma.masked_array(data, np.isnan(data))
            data_f = data * scale_factor
            os.remove(OUT+os.sep+'teste.hdf')
            
            if ((y0 > miny and y0 < maxy) or (y1 > miny and y1 < maxy)):
                print("dentro del rango latitud")
                if ((x0 > minx and x0 < maxx) or (x1 > minx and x1 < maxx)):
                    print("dentro del rango longitud")
                    posx0 = find_nearest_x(x,x0)
                    posx1 = find_nearest_x(x,x1)
                    posy0 = find_nearest_y(y,y0)
                    posy1 = find_nearest_y(y,y1)
                    
                    new_x0 = x[0,posx0]
                    new_x1 = x[0,posx1]
                    new_y0 = y[posy0,0]
                    new_y1 = y[posy1,0]
                    
                    y = y[posy0:posy1+1,posx0:posx1+1]
                    x = x[posy0:posy1+1,posx0:posx1+1]
                    
                    data_final = data_f[posy0:posy1+1,posx0:posx1+1]
                    # data_final[data_fil <= 0.000001] = np.nan
                    # mask = np.ma.getmask(data_final)
                    time = FILE_NAME[18:20]+':'+ FILE_NAME[20:22]
                    a,b = data_final.shape
#                    print data_final.shape
#                    print data_final.max(), data_final.min()  
                    if a<=1 or b<=1:
                        print("Area muy pequena")
                        continue 
                    if  np.nansum(data_final) != 0.0:
                        print("tem dados")
                        file_names = FILE_NAME[-44:-4]
                        create_nc(OUT_PATH,file_names,y,x,data_final,metadata)
                    else:
                        print("nao tem dados")                  
                else:
                    print("fuera del rango longitud")
            else:
                 print("fuera del rango latitud")   
            
        else:
            print("imagen demasiado grande")
            os.remove(OUT+os.sep+'teste.hdf')

                    
        
        
    

