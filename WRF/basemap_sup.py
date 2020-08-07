#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 18:00:40 2020

@author: noelia
"""

from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import numpy as np
from wrf import to_np
import pandas as pd

def plot_map(OUTPUT,wrf_file,gas,u10,v10,lat,lon,tiempo,name,unidad):
    # gas :  tipo do poluente
    # l: nivel de la atmosfera
    # lat: latitude (matriz)
    # lon: longitude (matriz)
    # tiempo : para cada instante de la imagen (gas.Time)
    # name: nome do poluente (gas.description)
    # unidad : unidade do poluente (gas.units)
    date = pd.DataFrame({'date': tiempo.values})
    date['local'] = date['date'].dt.tz_localize('UTC').dt.tz_convert('America/Sao_Paulo')    
    if str(len(gas.shape)) == '3':
        print('El gas tiene 3 dimensiones')
        for t in range(gas.shape[0]):
            levels = MaxNLocator(nbins=18).tick_values(gas[t].min(),gas[t].max())
            cmap = cm.get_cmap("Spectral_r",lut=25)
            cmap.set_bad("w")
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
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
            trend=m.pcolormesh(x, y, gas[t,:,:], cmap=cmap, norm = norm)
            m.quiver(x[points],y[points],u[points],v[points], latlon=True, scale=30, scale_units='inches')        
            cbar = m.colorbar(trend, location='right', pad="5%", ticks=levels)
            cbar.set_label(str(unidad), fontsize=16)
            ax.set_title("Sao Paulo Metropolitan Region" + '\n' +
                 str(name)+ " on the surface  "+'\n'+ 
                 "for "+ str(date['local'].dt.strftime('%d-%m-%Y %H:%M').values[t])
                             + " (Local Time)", fontsize=20)
            fig.savefig(OUTPUT+str(name)+"_sup_"+str(wrf_file[t][-30:-6])+'.png')
            plt.show()
#            plt.clf()
    else:
        print('ERRORR!!!!')
        print('Este script solo plotea en la superficie')
####################### convertir png para gif ################################
#filenames = ["formaldehido_"+str(i)+".png" for i in range(385)]
#images = " ".join([output + filename for filename in filenames])
#path = "/data/noelia/modelo/plot_gif/" 
#os.system('convert -alpha deactivate -verbose -delay 50 -loop 0 -density 90 {} {}formaldehido.gif'.format(images, path))
