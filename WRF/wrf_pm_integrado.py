#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 22:55:13 2020

@author: noelia
"""

###############################################################################
# Calculando la integral en la columna atmosferica para algunos poluentes    ##
#   se necesita correr el script basemap_integrado.py para plotear las       ##
#  imagenes: PM25, PM10. Para la superficie esta en wrf_sup ##
###############################################################################

from __future__ import print_function
from glob import glob
from netCDF4 import Dataset
from wrf import to_np, getvar, ALL_TIMES
import numpy as np
from metpy.calc import wind_components
from metpy.units import units

wrf_file = sorted(glob("/data/noelia/modelo/wrf_chem_10/wrfout_d01_*"))  ## archivos wrfout
output = "/data/noelia/modelo/plot_wrf_satelite/poluentes_integrado/"

print("leyendo los datos")
ncfiles = [Dataset(x) for x in wrf_file]  ## leer todos los wrfout
#print ncfiles
print("vientos")
wind = getvar(ncfiles, 'uvmet10_wspd_wdir', timeidx=ALL_TIMES, method='cat')
ws = wind.sel(wspd_wdir='wspd')
wd = wind.sel(wspd_wdir='wdir')
###############################################################################
##################### cambiar de m/s para knots ###############################
ws10 = (ws.values*units('m/s')).to('knots')
wd10 = wd.values*units.degree
##################### calcular los componentes U y V ##########################
u10,v10 = wind_components(ws10, wd10)

######################### obtener los gases y particulas ######################
print("leyendo PM25")
pm25 = getvar(ncfiles, "PM2_5_DRY", timeidx=ALL_TIMES, method='cat')   ### pm25 [ug/m3]
pm10 = getvar(ncfiles, "PM10", timeidx=ALL_TIMES, method='cat')   ### pm10 [ug/m3]

lat = ncfiles[0].variables['XLAT'][:]
lon = ncfiles[0].variables['XLONG'][:]

###################### Integrando en la columna ###############################
# [ppm](z) = (R / M_i)* T(z) * [ug/m3](z)
# columna = Sumatoria (diferencia)[ppm](z)*((diferencia)[P](z)/g)
pm25_u_col = np.zeros((pm25.shape[0],pm25.shape[2],pm25.shape[3]))
pm10_u_col = np.zeros((pm10.shape[0],pm10.shape[2],pm10.shape[3]))

print("Integrando PM25 en la columna vertical")
for t in range(pm25.shape[0]):
    for l in range(pm25.shape[1]):
        pm25_u_col[t,:,:] = pm25_u_col[t,:,:] + to_np(pm25[t,l,:,:])
  
print("Integrando PM10 en la columna vertical")
for t in range(pm10.shape[0]):
    for l in range(pm10.shape[1]):
        pm10_u_col[t,:,:] = pm10_u_col[t,:,:] + to_np(pm10[t,l,:,:])


###############################################################################
import basemap_integrado as bm_inter

print("Ploteando PM25")
unidad = str("$\mu g/m^{3}$")
descripcion = "PM25"
tiempo = pm25.Time
bm_inter.plot_map(output,wrf_file,pm25_u_col,u10,v10,lat,lon,tiempo,descripcion,unidad)

print("Ploteando PM10")
unidad = str("$\mu g/m^{3}$")
descripcion = "PM10"
tiempo = pm10.Time
bm_inter.plot_map(output,wrf_file,pm10_u_col,u10,v10,lat,lon,tiempo,descripcion,unidad)
