#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 22:58:42 2020

@author: noelia
"""
###############################################################################
# Calculando la integral en la columna atmosferica para algunos poluentes    ##
#   se necesita correr el script basemap_integrado.py para plotear las       ##
#  imagenes: PM25, PM10, CO, O3, NO, NO2. Para la superficie esta en wrf_sup ##
###############################################################################

from __future__ import print_function
from glob import glob
from netCDF4 import Dataset
from wrf import to_np, getvar, ALL_TIMES
import numpy as np
from metpy.calc import wind_components
from metpy.units import units

wrf_file = sorted(glob("/data/noelia/modelo/wrf_chem_9/wrfout_d01_*"))  ## archivos wrfout
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
print("leyendo BC")
bc1 = getvar(ncfiles, "BC1", timeidx=ALL_TIMES, method='cat')   ### [ug/kg-dry air] 
bc2 = getvar(ncfiles, "BC2", timeidx=ALL_TIMES, method='cat')   ### [ug/kg-dry air] 
bc = to_np(bc1) + to_np(bc2)
print("leyendo OC")
oc1 = getvar(ncfiles, "OC1", timeidx=ALL_TIMES, method='cat')   ### [ug/kg-dry air]
oc2 = getvar(ncfiles, "OC2", timeidx=ALL_TIMES, method='cat')   ### [ug/kg-dry air]
ec = to_np(oc1) + to_np(oc2)

lat = ncfiles[0].variables['XLAT'][:]
lon = ncfiles[0].variables['XLONG'][:]
###############################################################################
#################### calcular presion total ###################################
p = getvar(ncfiles, "P", timeidx=ALL_TIMES, method='cat') #perturbado [Pa]
pb = getvar(ncfiles, "PB", timeidx=ALL_TIMES, method='cat')  # estado base [Pa]
pt = to_np(p)+to_np(pb)    # presion total [Pa]
########################## calcular la temperatura ############################
# theta = tk*(p_0/pt)^kappa
# kappa = R/c_p = 0.2854
# tk = theta*(pt/p_0)^kappa
# p_0 = 1000  o 100000 Pa
the = getvar(ncfiles, "T", timeidx=ALL_TIMES, method='cat') # temperatura potencial perturbado [K]
the0 = 300  ##### temperatura potencial estado base [K]
theta = to_np(the) + the0  #####
p_0 = 100000 ## [Pa]
kappa = 0.2854
tk = theta*((pt/p_0)**kappa)   #### [K]
###############################################################################
################## convirtiendo [ug/kg-dry air] para [ug/m3] ##################
#
# (P*Mi)/R*T = rho_i (densidad del aire) 
# sustancia [ug/m3] = rho_i * sustancia [ug/kg-dry air]
#  
##################### calculo de la densidad del aire seco ####################
qvapor = getvar(ncfiles, "QVAPOR", timeidx=ALL_TIMES, method='cat')  # vapor de agua [kg/kg]
r = 286.9   #### [Jkg-1K-1]
pv = pt*to_np(qvapor) ### presion de vapor 
p_dry= pt - pv ### presion del aire seco
rho_dry=p_dry/(tk*r)     ### densidad de aire seco

ec_u = ec*rho_dry*(10**9)  # [ug/m3]
bc_u = bc*rho_dry*(10**9)  # [ug/m3]

###################### Integrando en la columna ###############################
# [ppm](z) = (R / M_i)* T(z) * [ug/m3](z)
# columna = Sumatoria (diferencia)[ppm](z)*((diferencia)[P](z)/g)
bc_u_col = np.zeros((bc1.shape[0],bc1.shape[2],bc1.shape[3]))
ec_u_col = np.zeros((oc1.shape[0],oc1.shape[2],oc1.shape[3]))

print("Integrando BC en la columna vertical")
for t in range(bc1.shape[0]):
    for l in range(bc1.shape[1]):
        bc_u_col[t,:,:] = bc_u_col[t,:,:] + bc[t,l,:,:]
bc_u_col[bc_u_col <= 0.0] = np.nan

print("Integrando EC en la columna vertical")
for t in range(oc1.shape[0]):
    for l in range(oc1.shape[1]):
        ec_u_col[t,:,:] = ec_u_col[t,:,:] + ec[t,l,:,:]
ec_u_col[ec_u_col <= 0.0] = np.nan

###############################################################################
import basemap_integrado as bm_inter

print("Ploteando BC")
unidad = str("$\mu g/m^{3}$")
descripcion = "BC"
tiempo = bc1.Time
bm_inter.plot_map(output,wrf_file,bc_u_col,u10,v10,lat,lon,tiempo,descripcion,unidad)

print("Ploteando EC")
unidad = str("$\mu g/m^{3}$")
descripcion = "EC"
tiempo = oc1.Time
bm_inter.plot_map(output,wrf_file,ec_u_col,u10,v10,lat,lon,tiempo,descripcion,unidad)


