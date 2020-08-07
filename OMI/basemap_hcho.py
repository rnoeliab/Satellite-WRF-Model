#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:39:04 2020

@author: noelia
"""

from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import numpy as np

def plot_map(OUTPUT,gas,lat,lon,time,name,units,filename):
    # gas :  tipo do poluente
    # lat: latitude (matriz)
    # lon: longitude (matriz)
    # name: nome do poluente 
    # units: unidade do poluente
    lev = MaxNLocator(nbins=18).tick_values(0.0,7.0)
    cmap = cm.get_cmap("Spectral_r",lut=25)
    cmap.set_bad("w")
    norm = BoundaryNorm(lev, ncolors=cmap.N, clip=True)
    fig = plt.figure(figsize=(15, 8))
    ax=plt.axes()
    m = Basemap(projection='cyl', resolution='h', llcrnrlat=-26.0, 
                urcrnrlat=-21.0,llcrnrlon=-49.0, urcrnrlon=-44.0)
    # m = Basemap(projection='cyl', resolution='l',llcrnrlat=np.nanmin(lat), urcrnrlat = np.nanmax(lat),
    #             llcrnrlon=np.nanmin(lon), urcrnrlon = np.nanmax(lon))
    m.drawcoastlines(linewidth=1.5)
    m.drawstates(linewidth=1.5)
    m.drawparallels(np.arange(-90., 120., 15), labels=[1, 0, 0, 0],fontsize=13)
    m.drawmeridians(np.arange(-180., 181., 25), labels=[0, 0, 0, 1],fontsize=13)
    m.readshapefile('/data/noelia/shapefile/RM_Sao_Paulo/transformed','sp')
#    x, y = m(lon, lat)
    trend=m.contourf(lon, lat, gas, levels = lev, cmap=cmap, norm=norm)
    cbar = m.colorbar(trend, location='right', pad="5%", ticks=lev, ax=ax)        
    cbar.set_label(str(units), fontsize=16)
    ax.set_title("Sao Paulo Metropolitan Region" + '\n' +
         str(name)+ '\n'+ " for "+ str(time[0])+ ' (Local Time)' , fontsize=20)
    fig.savefig(OUTPUT+str(name)+'_col_'+str(filename[-58:-7])+'.png')
    plt.show()
