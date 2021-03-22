#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 17:16:10 2020

@author: noelia
"""

from __future__ import print_function
from glob import glob
from netCDF4 import Dataset
from wrf import to_np, getvar, ALL_TIMES
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
from datetime import datetime
import pandas as pd
from metpy.calc import wind_components
from metpy.units import units



wrf_file = sorted(glob(".../wrfout_d01_*"))  ## archivos wrfout
output = "output path to save created files"

ncfiles = [Dataset(x) for x in wrf_file]  ## read the wrfout files
#print ncfiles

print("winds")
wind = getvar(ncfiles, 'uvmet10_wspd_wdir', timeidx=ALL_TIMES, method='cat')
ws = wind.sel(wspd_wdir='wspd')
wd = wind.sel(wspd_wdir='wdir')
###############################################################################
##################### change from m/s to knots ###############################
ws10 = (ws.values*units('m/s')).to('knots')
wd10 = wd.values*units.degree
##################### calculate the U and V components##########################
u10,v10 = wind_components(ws10, wd10)

tau1 = getvar(ncfiles, "TAUAER1", timeidx=ALL_TIMES, method='cat')   ### AOD in 300nm
tau2 = getvar(ncfiles, "TAUAER2", timeidx=ALL_TIMES, method='cat')   ### AOD in 400nm
tau3 = getvar(ncfiles, "TAUAER3", timeidx=ALL_TIMES, method='cat')   ### AOD in 600nm
tau4 = getvar(ncfiles, "TAUAER4", timeidx=ALL_TIMES, method='cat')   ### AOD in 1000nm 

lat = ncfiles[0].variables['XLAT'][:]
lon = ncfiles[0].variables['XLONG'][:]

angstrom = np.zeros((tau1.shape[0],tau1.shape[1],tau1.shape[2],tau1.shape[3]))
aod_550 = np.zeros((tau1.shape[0],tau1.shape[1],tau1.shape[2],tau1.shape[3]))
###############################################################################
#################### Calculate AOD in 550nm #####################
start=datetime.now()
for t in range(tau1.shape[0]):
    print("t=",t)
    for l in range(tau1.shape[1]):
        angstrom[t,l,:,:] = np.log(to_np(tau1[t,l,:,:])/to_np(tau4[t,l,:,:]))/(np.log((1000./300.)))
        aod_550[t,l,:,:] = to_np(tau2[t,l,:,:])*np.power((550./400.),-1*angstrom[t,l,:,:])
print(datetime.now()-start)

aod_550[np.isnan(aod_550)] = 0.0             
###############################################################################
print("Integrando en la columna vertical")
aod_550_col = np.zeros((tau1.shape[0],tau1.shape[2],tau1.shape[3]))
for t in range(tau1.shape[0]):
    for l in range(tau1.shape[1]):
        aod_550_col[t,:,:] = aod_550_col[t,:,:] + aod_550[t,l,:,:]
aod_550_col[aod_550_col <= 0.0] = np.nan

date = pd.DataFrame({'date': tau1.Time.values})
date['local'] = date['date'].dt.tz_localize('UTC').dt.tz_convert('America/Sao_Paulo')

levels = MaxNLocator(nbins=18).tick_values(0,np.nanmax(aod_550_col))
cmap = cm.get_cmap("Spectral_r",lut=25)
cmap.set_bad("w")
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

print("Plotting")
for t in range(tau1.Time.shape[0]):
    print(t)   
    
    fig = plt.figure(figsize=(15, 8))
    ax=plt.axes()
    m = Basemap(projection='cyl', resolution='h' ,llcrnrlat=lat.min(), urcrnrlat=lat.max(),
                llcrnrlon=lon.min(), urcrnrlon=lon.max())                
    m.drawcoastlines(linewidth=1.5)
    m.drawstates(linewidth=1.5)
    m.drawparallels(np.arange(-90., 120., 1), labels=[1, 0, 0, 0],fontsize=15)
    m.drawmeridians(np.arange(-180., 181., 1), labels=[0, 0, 0, 1],fontsize=15)
    m.readshapefile('/data/noelia/shapefile/RM_Sao_Paulo/transformed','sp')
    x,y = m(to_np(lon)[0], to_np(lat)[0])
    yy = np.arange(0,y.shape[0],8)
    xx = np.arange(0, x.shape[1],8)
    u = u10[t,:,:] 
    v = v10[t,:,:]
    points = np.meshgrid(yy, xx)
    trend=m.pcolormesh(x, y, aod_550_col[t,:,:], cmap=cmap, norm = norm)
    m.quiver(x[points],y[points],u[points],v[points], latlon=True, scale=30, scale_units='inches')        
    cbar = m.colorbar(trend, location='right', pad="5%", ticks=levels)            
    cbar.set_label("None", fontsize=16)
    ax.set_title("Sao Paulo Metropolitan Region" + '\n' +
          "AOD integrated in the troposphere "+'\n'+
         " for "+ str(date['local'].dt.strftime('%d-%m-%Y %H:%M').values[t])
         + " (Local Time)", fontsize=20)
    fig.savefig(output+'AOD_col_'+str(wrf_file[t][-30:-6])+'.png')
    plt.show()
#    plt.clf()

####################### convert png to gif ################################
#filenames = ["formaldehido_"+str(i)+".png" for i in range(385)]
#images = " ".join([output + filename for filename in filenames])
#path = "/data/noelia/modelo/plot_gif/" 
#os.system('convert -alpha deactivate -verbose -delay 50 -loop 0 -density 90 {} {}formaldehido.gif'.format(images, path))


