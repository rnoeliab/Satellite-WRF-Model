#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 14:28:00 2021

@author: noelia
"""

import glob
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
import xarray as xr
import pandas as pd
import basemap_modis as basemod

dire_mod = "../DATA/SP/regridded_2017/"
dire_wrf = '../DATA/'

output = "../results/plot_wrf_modis/"
########################### READ WRF IMAGES ###############################
ds_disk = xr.open_dataset(dire_wrf+"june_aod.avg.column.550_p1.nc")
aod_wrf = ds_disk["aod_055"].values
lon = ds_disk.XLONG.values
lat = ds_disk.XLAT.values
times = pd.DataFrame({'date': ds_disk.Time.values})
times['date_local'] = times['date'].dt.tz_localize('UTC').dt.tz_convert('America/Sao_Paulo')
times['date_qualar'] = times['date_local'].dt.strftime('%Y-%m-%d %H:%M:%S')
time_wrf = []
for ll in str(times['date_qualar'].tolist()).split(","):
    time_wrf.append(ll[2:15])

serie_time = pd.DataFrame(pd.date_range(start='2017-06-01', end='2017-06-16'), columns=["Date"])  
serie_time["Date"] = serie_time["Date"].dt.strftime('%Y-%m-%d') 
################### eliminando los dias repetidos ###############################
mj3 = []
mj2 = []   ### almacenando solo los datos aod de los dias no repetidos
for i in range(len(time_wrf)):
    if i == (len(time_wrf)-1):
        break;
    else:
        if time_wrf[i] == time_wrf[i+1]:
            pass;
        else:
            mj3.append(time_wrf[i]) 
            mj2.append(aod_wrf[i])
    if i == (len(time_wrf)-2):
        if time_wrf[i] == time_wrf[i+1]:
            mj3.append(time_wrf[len(time_wrf)-1])
            mj2.append(aod_wrf[len(time_wrf)-1])
        else:
            pass;
########################### READ MODIS IMAGES #############################
mylist=[]
listdirr = []
files1 = sorted(glob.glob(dire_mod+"*_regrid.nc"))  ## archivos wrfout
for i in range(len(files1)):
    if i == (len(files1)-1):
        break;      
    else:
        if files1[i][-50:-33] == files1[i+1][-50:-33]:  #### para la misma hora
            listdirr.append(files1[i])          
            if i == (len(files1)-2):
                listdirr.append(files1[len(files1)-1])
                mylist.append(listdirr)               
        else:
            listdirr.append(files1[i])
            mylist.append(listdirr)
            listdirr = []

##### Realizando un mosaico para los datos MODIS para la misma hora
aod_mod_total = np.zeros((len(mylist),aod_wrf.shape[1],aod_wrf.shape[2]))   
time_mod_total = []
for n in range(len(mylist)):
    if len(mylist[n]) != 1: #### cuando hay mas de 1 dato para la misma hora
        m=[]
        for FILE_NAME in mylist[n]:
            df = xr.open_dataarray(FILE_NAME)
            times_mod = df.attrs["time"] 
            aod_mod = df.values
            m.append(aod_mod)
        aod_mod_total[n] = np.nansum(np.dstack(m),2)  ### mosaico 
        time_mod_total.append(times_mod)
    else:
        for FILE_NAME in mylist[n]:
            df = xr.open_dataarray(FILE_NAME)
            times_mod = df.attrs["time"] 
            aod_mod_total[n] = df.values
            time_mod_total.append(times_mod)
            
######## comparando tiempo del modelo y de las imagenes satelitales ###########
aod_mod_total[aod_mod_total==0]=np.nan
list_aod_mod = []; list_aod_wrf = []; time_wrf_total = []
for i in range(len(time_mod_total)):
    for j in range(len(mj3)):
        if str(time_mod_total[i][0:13]) == str(mj3[j]):  ##### seleccion de tiempo
            time_wrf_total.append(mj3[j])
            for a in range(lon.shape[0]):
                for b in range(lat.shape[1]):
                    if str(aod_mod_total[i,a,b]) == str(np.nan):
                        mj2[j][a,b] = np.nan
            list_aod_mod.append(aod_mod_total[i])
            list_aod_wrf.append(mj2[j])
                        
########################### TERRA Y AQUA ######################################        
dat_wrf_t_day = [];dat_wrf_a_day = [];dat_mod_t_day = [];dat_mod_a_day = []
for n in range(len(serie_time)):
    day_wrf_t = [];day_wrf_a = [];day_mod_t = [];day_mod_a = [] 
    for o in range(len(time_mod_total)):
        ######################## media por dia ################################
        if str(serie_time["Date"][n]) == str(time_mod_total[o][0:10]):
            if int(time_mod_total[o][11:13]) < 12:   #### Terra
                day_mod_t.append(list_aod_mod[o])
                day_wrf_t.append(list_aod_wrf[o])
            else:
                day_mod_a.append(list_aod_mod[o])   #### Aqua
                day_wrf_a.append(list_aod_wrf[o])
    if len(day_mod_t) != 0:
        dat_mod_t_day.append(np.nanmean(day_mod_t,axis=0))
        dat_wrf_t_day.append(np.nanmean(day_wrf_t,axis=0))
    else:
        dat_mod_t_day.append(np.nan)
        dat_wrf_t_day.append(np.nan)
    if len(day_mod_a) != 0:
        dat_mod_a_day.append(np.nanmean(day_mod_a,axis=0))
        dat_wrf_a_day.append(np.nanmean(day_wrf_a,axis=0))
    else:
        dat_mod_a_day.append(np.nan)
        dat_wrf_a_day.append(np.nan)
    ################# day ###############################
    ############### TERRA ######################################
    if type(dat_mod_t_day[n]) == np.ndarray:
        basemod.plot_map_day(output,dat_mod_t_day[n],dat_wrf_t_day[n],lat,lon,serie_time["Date"][n],"Terra")
    ############### AQUA ######################################
    if type(dat_mod_a_day[n]) == np.ndarray:    
        basemod.plot_map_day(output,dat_mod_a_day[n],dat_wrf_a_day[n],lat,lon,serie_time["Date"][n],"Aqua")


################ seleccionando por area ################################
area = list_aod_mod[0].shape[0]*list_aod_mod[0].shape[1]
porcentaje = [80,70,50,30,20,10] #### porcentajes de numeros no nan
for i in range(len(porcentaje)):
    por = porcentaje[i]*(area)/100  ### cantidad de nan aceptables
    dat_mt = []; dat_wt = []; dat_ma = []; dat_wa = []
    dat_wrf_t_day = [];dat_wrf_a_day = [];dat_mod_t_day = [];dat_mod_a_day = []
    for n in range(len(serie_time)):
        day_wrf_t = [];day_wrf_a = [];day_mod_t = [];day_mod_a = [] 
        for o in range(len(time_mod_total)):
            ######################## media por dia ################################
            if str(serie_time["Date"][n]) == str(time_mod_total[o][0:10]):
                if int(time_mod_total[o][11:13]) < 12:   #### Terra
                    day_mod_t.append(list_aod_mod[o])
                    day_wrf_t.append(list_aod_wrf[o])
                else:
                    day_mod_a.append(list_aod_mod[o])   #### Aqua
                    day_wrf_a.append(list_aod_wrf[o])
        if len(day_mod_t) != 0:
            dat_mod_t_day.append(np.nanmean(day_mod_t,axis=0))
            dat_wrf_t_day.append(np.nanmean(day_wrf_t,axis=0))
        else:
            dat_mod_t_day.append(np.nan)
            dat_wrf_t_day.append(np.nan)
        if len(day_mod_a) != 0:
            dat_mod_a_day.append(np.nanmean(day_mod_a,axis=0))
            dat_wrf_a_day.append(np.nanmean(day_wrf_a,axis=0))
        else:
            dat_mod_a_day.append(np.nan)
            dat_wrf_a_day.append(np.nan)

        if type(dat_mod_t_day[n]) == np.ndarray:
            if (area - np.isnan(dat_mod_t_day[n]).sum() >= int(por)):
                dat_mt.append(dat_mod_t_day[n])
            if (area - np.isnan(dat_wrf_t_day[n]).sum() >= int(por)):
                dat_wt.append(dat_wrf_t_day[n])
        if type(dat_mod_a_day[n]) == np.ndarray:    
            if (area - np.isnan(dat_mod_a_day[n]).sum() >= int(por)):
                dat_ma.append(dat_mod_a_day[n])
            if (area - np.isnan(dat_wrf_a_day[n]).sum() >= int(por)):
                dat_wa.append(dat_wrf_a_day[n])

    if (len(dat_mt) != 0 and len(dat_wt) != 0):      
        ############### TERRA ######################################
        dat_mod_t = np.nanmean(dat_mt,axis=0)
        dat_wrf_t = np.nanmean(dat_wt,axis=0)
        basemod.plot_map_month(output,dat_mod_t,dat_wrf_t,lat,lon,"Terra",porcentaje[i],len(dat_mt))
    if (len(dat_ma) != 0 and len(dat_wa) != 0):      
        ############### AQUA ######################################
        dat_mod_a = np.nanmean(dat_ma,axis=0)
        dat_wrf_a = np.nanmean(dat_wa,axis=0)
        basemod.plot_map_month(output,dat_mod_a,dat_wrf_a,lat,lon,"Aqua",porcentaje[i],len(dat_ma))


######################## media por periodo ################################
for n in range(len(dat_mod_t_day)):    
    if type(dat_mod_t_day[n]) != np.ndarray:
        del dat_mod_t_day[n]
        del dat_wrf_t_day[n]
for n in range(len(dat_mod_a_day)):  
    if type(dat_mod_a_day[n]) != np.ndarray:
        del dat_mod_a_day[n]
        del dat_wrf_a_day[n]
dat_mod_t = np.nanmean(dat_mod_t_day,axis=0)
dat_wrf_t = np.nanmean(dat_wrf_t_day,axis=0)
basemod.plot_map_month(output,dat_mod_t,dat_wrf_t,lat,lon,"Terra","100",len(dat_mod_t_day))

dat_mod_a = np.nanmean(dat_mod_a_day,axis=0)
dat_wrf_a = np.nanmean(dat_wrf_a_day,axis=0)
basemod.plot_map_month(output,dat_mod_a,dat_wrf_a,lat,lon,"Aqua","100",len(dat_mod_a_day))


diff_t = dat_mod_t - dat_wrf_t
diff_a = dat_mod_a - dat_wrf_a

levels = MaxNLocator(nbins=10).tick_values(0,0.5)
cmap = cm.get_cmap("Blues",lut=25)
cmap.set_bad("w")
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig = plt.figure(figsize=(11, 11))
rect = fig.patch
rect.set_facecolor('lightgoldenrodyellow')
ax0 = fig.add_subplot(111, frame_on=False)
ax0.set_xticks([])
ax0.set_yticks([])
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.0)
m = Basemap(projection='cyl', resolution='h', llcrnrlat=-25.2, urcrnrlat=-21.5,
            llcrnrlon=-49.0, urcrnrlon=-44.3)
m.drawcoastlines(linewidth=1.5)
m.drawstates(linewidth=1.5)    
m.drawparallels(np.arange(-90., 120., 1), labels=[1, 0, 0, 0],fontsize=18)
m.drawmeridians(np.arange(-180., 181., 1), labels=[0, 0, 0, 1],fontsize=18) 
trend=m.pcolormesh(lon, lat, diff_t, cmap=cmap, norm = norm)
cbar = m.colorbar(trend, location='right', pad="5%", ticks=levels)
cbar.set_label('None', fontsize=19)
ax.set_title("Difference between Terra Satellite and Model Data" +"\n" +
             "over SPMR for June", fontsize=20)

fig = plt.figure(figsize=(11, 11))
rect = fig.patch
rect.set_facecolor('lightgoldenrodyellow')
ax0 = fig.add_subplot(111, frame_on=False)
ax0.set_xticks([])
ax0.set_yticks([])
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.0)
m = Basemap(projection='cyl', resolution='h', llcrnrlat=-25.2, urcrnrlat=-21.5,
            llcrnrlon=-49.0, urcrnrlon=-44.3)
m.drawcoastlines(linewidth=1.5)
m.drawstates(linewidth=1.5)    
m.drawparallels(np.arange(-90., 120., 1), labels=[1, 0, 0, 0],fontsize=18)
m.drawmeridians(np.arange(-180., 181., 1), labels=[0, 0, 0, 1],fontsize=18) 
trend=m.pcolormesh(lon, lat, diff_a, cmap=cmap, norm = norm)
cbar = m.colorbar(trend, location='right', pad="5%", ticks=levels)
cbar.set_label('None', fontsize=19)
ax.set_title("Difference between Aqua Satellite and Model Data" +"\n" +
             "over SPMR for June", fontsize=20)
plt.show()


