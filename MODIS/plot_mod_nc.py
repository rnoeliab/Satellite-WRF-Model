#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 22:14:43 2020

@author: noelia
"""
import netCDF4 as nc4
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import date
from glob import glob
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator

#open several data
INPUT_PATH = '/data/noelia/imagen_data/modis/DATA/SP/3K_nc/2017/'
OUT_PATH = '/data/noelia/imagen_data/modis/results/plot_nc_modis_rmsp/2017/'
mylist=[]
listdirr = []

DATAFIELD_NAME = 'Optical_Depth_Land_And_Ocean' 

files = sorted(glob(INPUT_PATH+"*.nc"))  ## archivos wrfout

for i in range(len(files)):
    if i == (len(files)-1):
        break;      
    else:
        if files[i][-43:-26] == files[i+1][-43:-26]:
            listdirr.append(files[i])          
            if i == (len(files)-2):
                listdirr.append(files[len(files)-1])
                mylist.append(listdirr)               
        else:
            listdirr.append(files[i])
            mylist.append(listdirr)
            listdirr = []

levels = MaxNLocator(nbins=18).tick_values(0,0.6)
cmap = cm.get_cmap("Spectral_r",lut=25)
cmap.set_bad("w")
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
for n in range(len(mylist)):
    fig = plt.figure(figsize=(15, 8))
    ax=plt.axes()
    for FILENAME in mylist[n]:
        print(FILENAME)
        fh = nc4.Dataset(FILENAME, mode = "r")
        df = fh.groups["AOD_3K"]
        # print df.variables.keys()
        lons = df.variables['Longitude'][:]
        lats = df.variables['Latitude'][:]
        aod = df.variables['Optical_Depth_Land_And_Ocean'][:]
        aod_units = df.variables['Optical_Depth_Land_And_Ocean'].units
        lon, lat = np.meshgrid(lons, lats)
        print(lat.max(),lat.min())
        print(lon.max(), lon.min())
        time = df.variables['Time'][:]
        dt = date.fromordinal(time[0]).strftime('%d-%m-%Y')

        m = Basemap(projection='cyl', resolution='h',
                    llcrnrlat=-26.0, urcrnrlat=-21.0,
                    llcrnrlon=-49.0, urcrnrlon=-44.0)
        m.drawcoastlines(linewidth=1.5)
        m.drawstates(linewidth=1.5)    
        m.drawparallels(np.arange(-90., 120., 1), labels=[1, 0, 0, 0],fontsize=15)
        m.drawmeridians(np.arange(-180., 181., 1), labels=[0, 0, 0, 1],fontsize=15) 
        m.readshapefile('/data/noelia/shapefile/RM_Sao_Paulo/transformed','sp')
        x, y = m(lon, lat)
        trend=m.pcolormesh(x, y, aod, cmap=cmap, norm = norm)
    cbar = m.colorbar(trend, location='right', pad="5%", ticks=levels)
    # label colorboar
    cbar.set_label('None', fontsize=16)
    ax.set_title("Sao Paulo Metropolitan Region" + '\n' +
         "AOD  from satellite "+'\n'+ 
         "for "+ str(dt)+ " "+ str(FILENAME[-25:-23])+":"+str(FILENAME[-23:-21])  
         + " (Local Time)", fontsize=20)
    fig.savefig(OUT_PATH+str(FILENAME[-43:-21])+'.png')
    plt.show()

#
#fig = plt.figure(figsize=(10,8))
#
#for n in range(len(mylist)):
#    for FILENAME in mylist[n]:
#        print(FILENAME)
#        fh = nc4.Dataset(FILENAME, mode = "r")
#        df = fh.groups["AOD_3K"]
#        # print df.variables.keys()
#        lons = df.variables['Longitude'][:]
#        lats = df.variables['Latitude'][:]
#        aod = df.variables['AOD'][:]
#        aod_units = df.variables['AOD'].units
#        lon, lat = np.meshgrid(lons, lats)
#        time = df.variables['Time'][:]
#        dt = date.fromordinal(time[0])
#        print lat.min(), lon.min()
#        print lat.max(), lon.max()
#        
#        
#
#levels = MaxNLocator(nbins=25).tick_values(0.0, 0.5)
#cmap = cm.get_cmap("jet",lut=25)
#cmap.set_bad("w")
#norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
#
#m = Basemap(projection='cyl',area_thresh=10000, resolution='l' ,llcrnrlat=-26.0, urcrnrlat=-21.0,
#    llcrnrlon=-49.0, urcrnrlon=-44.0) 
#m.drawcoastlines(linewidth=1.5)
#m.drawstates(linewidth=1.5)
#m.drawparallels(np.arange(-90., 120., 1), labels=[1, 0, 0, 0],fontsize=13)
#m.drawmeridians(np.arange(-180., 181., 1), labels=[0, 0, 0, 1],fontsize=13) 
#m.readshapefile('/media/noelia/TOSHIBA EXT/doctorado/usp/shapefile/RM_Sao_Paulo/transformed','sp')
#x, y = m(lon, lat)
#trend=m.pcolormesh(x, y, aod, cmap=cmap, norm = norm)
#cbar = m.colorbar(trend, location='right', pad="8%", ticks=levels)
#
## label colorboar
#cbar.set_label('AOD', fontsize=16)
#
## title the plot
#plotTitle = FILENAME[-43:-35]
#plt.title(plotTitle +"\n" + 'Time: '+str(dt) + " " + str(FILENAME[-25:-23]) + "-" + str(FILENAME[-23:-21]))
#fig.savefig(OUT_PATH+FILENAME[-43:-21]+'_rmsp.tif')
## Show the plot window.
#plt.show()

