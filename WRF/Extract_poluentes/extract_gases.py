#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 12:46:50 2021

@author: noelia
"""

from __future__ import print_function
from glob import glob
from netCDF4 import Dataset
import xarray as xr
from wrf import to_np, getvar, ALL_TIMES
import numpy as np
from datetime import datetime
import pandas as pd

f = open('../omi/OMNO2G_ori/data/2017_list_omi.txt')
satellite_txt = f.readlines()

wrf_file = sorted(glob("../wrfout_d02_*"))  ## archivos wrfout

output = "../modelo/DATA/"

#### find the intersection between satellite and WRF times 
times = [datetime.strptime(i[-19:-9], '%Y-%m-%d') for i in wrf_file] 
listdir = []
for n in range(len(times)):
    name = '_2017m'+str(times[n])[5:7]+str(times[n])[8:10]
    mod_wrf = sorted(list(filter(lambda a: str(name) in a,satellite_txt)))
    listdir += mod_wrf
with open('../list_mod_wrf.txt', 'w') as f:
     for item in listdir:
         f.write(item)   
###############################################################################    
############# Read WRF data, only the variable of interest  ###########
ncfiles = [Dataset(x) for x in wrf_file]  ## leer todos los wrfout
#print ncfiles

tau1 = getvar(ncfiles, "TAUAER1", timeidx=ALL_TIMES, method='cat')   ### AOD in 300nm
tau2 = getvar(ncfiles, "TAUAER2", timeidx=ALL_TIMES, method='cat')   ### AOD in 400nm
tau3 = getvar(ncfiles, "TAUAER3", timeidx=ALL_TIMES, method='cat')   ### AOD in 600nm
tau4 = getvar(ncfiles, "TAUAER4", timeidx=ALL_TIMES, method='cat')   ### AOD in 1000nm 

angstrom = np.zeros((tau1.shape[0],tau1.shape[1],tau1.shape[2],tau1.shape[3]))
aod_550 = np.zeros((tau1.shape[0],tau1.shape[1],tau1.shape[2],tau1.shape[3]))
###############################################################################
######################### Calculate AOD in 550nm ##############################
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
aod_550_col[aod_550_col < 0.0] = np.nan

date_wrf = pd.DataFrame({'date': tau1.Time.values})   ##### Times in UTC

aod_col = np.zeros((len(listdir),aod_550_col.shape[1],aod_550_col.shape[2]))

ddt = []; dtt = []; dttt = []
#### find the intersection between satellite and WRF times (all in julian days)
count = 0
for n in range(len(djulian)):
    hours_wrf = date_wrf['date'].dt.strftime('%H:%M')[n*24:(n*24)+24]  ### 00hrs - 24hrs
    name = '_3K.A2017'+str(djulian[n])+'.'
    mod_wrf = sorted(list(filter(lambda a: str(name) in a,satellite_txt)))
    ################# Selecting only the hours of interest ####################
    for j in mod_wrf:
        hours_mod = j[18:20] + ':00'
        h_mod_wrf = hours_wrf[j[18:20] + ':00' == hours_wrf].index.values[0]       
        h_mod_wrf_1 = hours_wrf[str(int(j[18:20])+1) + ':00' == hours_wrf].index.values[0]
        print(h_mod_wrf,h_mod_wrf_1)
        ddt.append(tau4.coords['XTIME'].values[h_mod_wrf])
        dtt.append(tau4.coords['Time'].values[h_mod_wrf])
        dttt.append(tau4.coords['datetime'].values[h_mod_wrf])
        aod_col[count,:,:] = (aod_550_col[h_mod_wrf,:,:] + aod_550_col[h_mod_wrf_1,:,:])/2
        count = count + 1
###############################################################################
########### save the information in other netcdf with xarray ##################
da = xr.DataArray(
    data=aod_col,
    dims = ('Time', 'south_north', 'west_east'),
    coords = dict(
        XLONG = (['south_north', 'west_east'], tau4.coords['XLONG'].values),
        XLAT = (['south_north', 'west_east'], tau4.coords['XLAT'].values),
        XTIME = (['Time'],np.array(ddt)),
        Time = ('Time', dtt),
        datetime = ('Time', dttt)
        ),
    attrs = tau4.attrs)
da.name = "AOD_550nm"
da.attrs['description'] = '550nm optical thickness'
del da.attrs['projection']
da.to_netcdf(output+"aod.avg.column.550.nc")

