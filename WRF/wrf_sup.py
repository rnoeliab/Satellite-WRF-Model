#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 11:44:59 2020

@author: noelia
"""

from __future__ import print_function
from glob import glob
from netCDF4 import Dataset
from wrf import to_np, getvar, ALL_TIMES
from metpy.calc import wind_components
from metpy.units import units
 
print('leyendo el directorio')
INPUT = '/data/noelia/modelo/DATA/wrfout_'
name_out = str(input("Put the name file [mechanism]_[regional/local]: "))

wrf_file = sorted(glob(INPUT+name_out+"/wrfout_d01_*"))  ## archivos wrfout
output = "/data/noelia/modelo/figures/plot_wrf_satelite/pol_sup_"+str(name_out)+"/"

ncfiles = [Dataset(x) for x in wrf_file]  ## leer todos los wrfout
#print ncfiles

print('leyendo los poluentes')
######################### obtener los gases y particulas ######################
tc = getvar(ncfiles, 'T2', timeidx=ALL_TIMES, method='cat')
psfc = getvar(ncfiles, 'PSFC', timeidx=ALL_TIMES, method='cat')

print('vientos')
wind = getvar(ncfiles, 'uvmet10_wspd_wdir', timeidx=ALL_TIMES, method='cat')
ws = wind.sel(wspd_wdir='wspd')
wd = wind.sel(wspd_wdir='wdir')

print('PM')
pm25 = getvar(ncfiles, "PM2_5_DRY", timeidx=ALL_TIMES, method='cat')   ### pm25 [ug/m3]
pm25_s = pm25.isel(bottom_top=0)

pm10 = getvar(ncfiles, "PM10", timeidx=ALL_TIMES, method='cat')   ### pm10 [ug/m3]
pm10_s = pm10.isel(bottom_top=0)

print('gases')
no = getvar(ncfiles, "no", timeidx=ALL_TIMES, method='cat')   ### NO [ppmv]
no_s = no.isel(bottom_top=0)

no2 = getvar(ncfiles, "no2", timeidx=ALL_TIMES, method='cat')   ### NO2 [ppmv]
no2_s = no2.isel(bottom_top=0)

o3 = getvar(ncfiles, "o3", timeidx=ALL_TIMES, method='cat')   ### O3 [ppmv]
o3_s = o3.isel(bottom_top=0)

co = getvar(ncfiles, "co", timeidx=ALL_TIMES, method='cat')   ### CO [ppmv]
co_s = co.isel(bottom_top=0)

print('Formaldeido')
hcho = getvar(ncfiles, "hcho", timeidx=ALL_TIMES, method='cat')  ### HCHO [ppmv]
hcho_s =hcho.isel(bottom_top=0)

print('BC')
bc1 = getvar(ncfiles, "BC1", timeidx=ALL_TIMES, method='cat')   ### [ug/kg-dry air] 
bc2 = getvar(ncfiles, "BC2", timeidx=ALL_TIMES, method='cat')   ### [ug/kg-dry air] 
bc1_s = bc1.isel(bottom_top=0)
bc2_s = bc2.isel(bottom_top=0)

bc = to_np(bc1) + to_np(bc2)
bc_s = to_np(bc1_s) + to_np(bc2_s)

print('OC')
oc1 = getvar(ncfiles, "OC1", timeidx=ALL_TIMES, method='cat')   ### [ug/kg-dry air]
oc2 = getvar(ncfiles, "OC2", timeidx=ALL_TIMES, method='cat')   ### [ug/kg-dry air]
oc1_s = oc1.isel(bottom_top=0)
oc2_s = oc2.isel(bottom_top=0)

ec = to_np(oc1) + to_np(oc2)
ec_s = to_np(oc1_s) + to_np(oc2_s)

lat = ncfiles[0].variables['XLAT'][:]
lon = ncfiles[0].variables['XLONG'][:]

print('convirtiendo las unidades')
###############################################################################
################## convirtiendo [ug/kg-dry air] para [ug/m3] ##################
#
# (P*Mi)/R*T = rho_i (densidad del aire) 
# sustancia [ug/m3] = rho_i * sustancia [ug/kg-dry air]
#  m_ar = 28.97*(10**-3)    ### [kg/mol]
mol = 6.023*(10**23)
R = 8.3143
bc_u = bc_s*to_np(psfc)*(28.97*(10**-3))/(R*to_np(tc))
ec_u = ec_s*to_np(psfc)*(28.97*(10**-3))/(R*to_np(tc))
########################## Changing ppm to ug/m3 ##############################
# [ug/m3] = [ppm] * P * M_i / (R * T) 
# R = 8.3143 J/K mol
# P in Pa
# T in K
# WRF-Chem gas units in ppmv
o3_u = to_np(o3_s)*to_np(psfc)*(16*3)/(R*to_np(tc))
no_u = to_np(no_s)*to_np(psfc)*(14+16)/(R*to_np(tc))
no2_u = to_np(no_s)*to_np(psfc)*(14+2*16)/(R*to_np(tc))
hcho_u = to_np(hcho_s)*to_np(psfc)*(30.03)/(R*to_np(tc)) 

o3_moleculas = to_np(o3_s)*(to_np(psfc)*(10**-6)/(R*to_np(tc)))*mol**(10**-6)
no_moleculas = to_np(no_s)*(to_np(psfc)*(10**-6)/(R*to_np(tc)))*mol**(10**-6)
no2_moleculas = (to_np(no2_s)*(to_np(psfc)*(10**-6)/(R*to_np(tc)))*mol)*(10**-6)
co_moleculas = (to_np(co_s)*(to_np(psfc)*(10**-6)/(R*to_np(tc)))*mol)*(10**-6)
hcho_moleculas = (to_np(hcho_s)*(to_np(psfc)*(10**-6)/(R*to_np(tc)))*mol)*(10**-6)

###############################################################################
##################### cambiar de m/s para knots ###############################
ws10 = (ws.values*units('m/s')).to('knots')
wd10 = wd.values*units.degree
##################### calcular los componentes U y V ##########################
u10,v10 = wind_components(ws10, wd10)

################################## ploteando ##################################
print("ploteando")
import basemap_sup as bm_sup

unidad = str("$\mu g/m^{3}$")
descripcion = "PM25"
tiempo = pm25.Time
bm_sup.plot_map(output,wrf_file,to_np(pm25_s),u10,v10,lat,lon,tiempo,descripcion,unidad)

unidad = str("$\mu g/m^{3}$")
descripcion = "PM10"
tiempo = pm10.Time
bm_sup.plot_map(output,wrf_file,pm10_s,u10,v10,lat,lon,tiempo,descripcion,unidad)

unidad = str("$\mu g/m^{3}$")
descripcion = "NO"
tiempo = no.Time
bm_sup.plot_map(output,wrf_file,no_u,u10,v10,lat,lon,tiempo,descripcion,unidad)

no2_mole = no2_moleculas*(10**-12)
unidad = str("$x 10^{12} moleculas/cm^{3}$")
descripcion = "NO2"
tiempo = no2.Time
bm_sup.plot_map(output,wrf_file,no2_mole,u10,v10,lat,lon,tiempo,descripcion,unidad)

co_mole = co_moleculas*(10**-14)
unidad = str("$x 10^{14} moleculas/cm^{3}$")
descripcion = "CO"
tiempo = co.Time
bm_sup.plot_map(output,wrf_file,co_mole,u10,v10,lat,lon,tiempo,descripcion,unidad)

unidad = str("$ppmv$")
descripcion = "CO"
tiempo = co.Time
bm_sup.plot_map(output,wrf_file,to_np(co_s),u10,v10,lat,lon,tiempo,descripcion,unidad)

unidad = str("$\mu g/m^{3}$")
descripcion = "O3"
tiempo = o3.Time
bm_sup.plot_map(output,wrf_file,o3_u,u10,v10,lat,lon,tiempo,descripcion,unidad)

unidad = str("$\mu g/m^{3}$")
descripcion = "BC"
tiempo = bc1.Time
bm_sup.plot_map(output,wrf_file,bc_u,u10,v10,lat,lon,tiempo,descripcion,unidad)

unidad = str("$\mu g/m^{3}$")
descripcion = "EC"
tiempo = oc1.Time
bm_sup.plot_map(output,wrf_file,ec_u,u10,v10,lat,lon,tiempo,descripcion,unidad)

hcho_mole = hcho_moleculas*(10**-12)
unidad = str("$x 10^{12} moleculas/cm^{3}$")
descripcion = "HCHO"
tiempo = hcho.Time
bm_sup.plot_map(output,wrf_file,hcho_mole,u10,v10,lat,lon,tiempo,descripcion,unidad)
