#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:30:50 2021

@author: noelia
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator

levels = MaxNLocator(nbins=10).tick_values(0,0.6)
cmap = cm.get_cmap("Spectral_r",lut=25)
cmap.set_bad("w")
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
def plot_map_month(output,dat_mod,dat_wrf,lat,lon,satelite,p,n):
    fig = plt.figure(figsize=(22, 11))
    rect = fig.patch
    rect.set_facecolor('lightgoldenrodyellow')
    ax0 = fig.add_subplot(111, frame_on=False)
    ax0.set_xticks([])
    ax0.set_yticks([])
    ax = fig.add_subplot(121)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(3.0)
    m = Basemap(projection='cyl', resolution='h', llcrnrlat=-25.2, urcrnrlat=-21.5,
                llcrnrlon=-49.0, urcrnrlon=-44.3)
    m.drawcoastlines(linewidth=1.5)
    m.drawstates(linewidth=1.5)    
    m.drawparallels(np.arange(-90., 120., 1), labels=[1, 0, 0, 0],fontsize=20)
    m.drawmeridians(np.arange(-180., 181., 1), labels=[0, 0, 0, 1],fontsize=20) 
    m.readshapefile('../shapefile/RM_Sao_Paulo/transformed','sp')
    trend=m.pcolormesh(lon, lat, dat_wrf, cmap=cmap, norm = norm)
    ax.set_title("Sao Paulo Metropolitan Region" + '\n' +
                  "AOD  from WRF model "+'\n'+ 
                  "for June    n =" + str(n), fontsize=20)

    ax = fig.add_subplot(122)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(3.0)
    m = Basemap(projection='cyl', resolution='h', llcrnrlat=-25.2, urcrnrlat=-21.5,
                llcrnrlon=-49.0, urcrnrlon=-44.3)
    m.drawcoastlines(linewidth=1.5)
    m.drawstates(linewidth=1.5)    
    m.drawparallels(np.arange(-90., 120., 1), labels=[1, 0, 0, 0],fontsize=20)
    m.drawmeridians(np.arange(-180., 181., 1), labels=[0, 0, 0, 1],fontsize=20) 
    m.readshapefile('../shapefile/RM_Sao_Paulo/transformed','sp')
    trend2 = m.pcolormesh(lon,lat, dat_mod, cmap=cmap, norm = norm)            
    cbar = m.colorbar(trend, location='right', pad="5%", ticks=levels)
    cbar.set_label('None', fontsize=19)
    cbar.ax.tick_params(labelsize=18) 
    ax.set_title("Sao Paulo Metropolitan Region" + '\n' +
                  "AOD  from Regridded "+str(satelite)+" satellite "+'\n'+ 
                  "for June    n =" + str(n), fontsize=20)
    fig.savefig(output+"plot_periodo/p2/june"+"_"+str(satelite)+"_"+str(p)+'.png',bbox_inches='tight')
    plt.show()

def plot_map_day(output,dat_mod,dat_wrf,lat,lon,time,satelite):
    fig = plt.figure(figsize=(22, 11))
    rect = fig.patch
    rect.set_facecolor('lightgoldenrodyellow')
    ax0 = fig.add_subplot(111, frame_on=False)
    ax0.set_xticks([])
    ax0.set_yticks([])
    ax = fig.add_subplot(121)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(3.0)
    m = Basemap(projection='cyl', resolution='h', llcrnrlat=-25.2, urcrnrlat=-21.5,
                llcrnrlon=-49.0, urcrnrlon=-44.3)
    m.drawcoastlines(linewidth=1.5)
    m.drawstates(linewidth=1.5)    
    m.drawparallels(np.arange(-90., 120., 1), labels=[1, 0, 0, 0],fontsize=20)
    m.drawmeridians(np.arange(-180., 181., 1), labels=[0, 0, 0, 1],fontsize=20) 
    m.readshapefile('../shapefile/RM_Sao_Paulo/transformed','sp')
    trend=m.pcolormesh(lon, lat, dat_wrf, cmap=cmap, norm = norm)
    ax.set_title("Sao Paulo Metropolitan Region" + '\n' +
                  "AOD  from WRF model "+'\n'+ 
                  "for "+ str(time), fontsize=20)

    ax = fig.add_subplot(122)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(3.0)
    m = Basemap(projection='cyl', resolution='h', llcrnrlat=-25.2, urcrnrlat=-21.5,
                llcrnrlon=-49.0, urcrnrlon=-44.3)
    m.drawcoastlines(linewidth=1.5)
    m.drawstates(linewidth=1.5)    
    m.drawparallels(np.arange(-90., 120., 1), labels=[1, 0, 0, 0],fontsize=20)
    m.drawmeridians(np.arange(-180., 181., 1), labels=[0, 0, 0, 1],fontsize=20) 
    m.readshapefile('../shapefile/RM_Sao_Paulo/transformed','sp')
    trend2 = m.pcolormesh(lon,lat, dat_mod, cmap=cmap, norm = norm)            
    cbar = m.colorbar(trend, location='right', pad="5%", ticks=levels)
    cbar.set_label('None', fontsize=19)
    cbar.ax.tick_params(labelsize=18) 
    ax.set_title("Sao Paulo Metropolitan Region" + '\n' +
                  "AOD  from Regridded "+str(satelite)+" satellite "+'\n'+ 
                  "for "+ str(time), fontsize=20)
    fig.savefig(output+"plot_day/p2/"+str(time.replace("-","_"))+"_"+str(satelite)+'.png',bbox_inches='tight')
    plt.show()

