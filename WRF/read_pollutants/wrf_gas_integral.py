#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 22:58:40 2020

@author: noelia
"""
###############################################################################
# Calculando la integral en la columna atmosferica para algunos poluentes    ##
#   se necesita correr el script basemap_integrado.py para plotear las       ##
#  imagenes: CO, O3, NO, NO2, HCHO. Para la superficie esta en wrf_sup ##
###############################################################################

from __future__ import print_function
from glob import glob
from netCDF4 import Dataset
from wrf import to_np, getvar, ALL_TIMES
import numpy as np
from metpy.calc import wind_components
from metpy.units import units

print('leyendo el directorio')
INPUT = '/data/noelia/modelo/DATA/wrfout_'
name_out = str(input("Put the name file [mechanism]_[regional/local]: "))

wrf_file = sorted(glob(INPUT+name_out+"/wrfout_d01_*"))  ## archivos wrfout
output = "/data/noelia/modelo/figures/plot_wrf_satelite/pol_integ_"+str(name_out)+"/1_dominio/"

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
print("leyendo gases")
co = getvar(ncfiles, "co", timeidx=ALL_TIMES, method='cat')   ###  CO [ppmv]
no = getvar(ncfiles, "no", timeidx=ALL_TIMES, method='cat')   ###  NO [ppmv]
no2 = getvar(ncfiles, "no2", timeidx=ALL_TIMES, method='cat') ### NO2 [ppmv]
o3 = getvar(ncfiles, "o3", timeidx=ALL_TIMES, method='cat')   ###  O3 [ppmv]
hcho = getvar(ncfiles, "hcho", timeidx=ALL_TIMES, method='cat') ### HCHO [ppmv]

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
######################## convertir [ppmv] to [ug/m3] or [moleculas/m3] ##########################
# PV = nnRT
# 1 mol = 6.023x10^23 moleculas
# [ug/m3] = [ppm] * P * M_i / (R * T) 
# [moleculas/m3] = [ppm] * (P * (10**-6) / (R * T)) * mol
#R = 8.3143 J/K mol
# P in Pa
# T in K
# WRF-Chem gas units in ppmv to moleculas/m3
mol = 6.023*(10**23)
R = 8.314472
print("convirtiendo de ppmv para moleculas")
o3_u = to_np(o3)*(pt*(16 * 3)/(R*tk))
no_u = to_np(no)*(pt*(14 + 16)/(R*tk))
no2_u = to_np(no2)*(pt*(14 + 2*16)/(R*tk))
co_u = to_np(co)*(pt*(12+ 16)/(R*tk))
hcho_u = to_np(hcho)*(pt*(30.03)/(R*tk))  

o3_moleculas = to_np(o3)*(pt*(10**-6)/(R*tk))*mol
no_moleculas = to_np(no)*(pt*(10**-6)/(R*tk))*mol
no2_moleculas = to_np(no2)*(pt*(10**-6)/(R*tk))*mol
co_moleculas = to_np(co)*(pt*(10**-6)/(R*tk))*mol
hcho_moleculas = to_np(hcho)*(pt*(10**-6)/(R*tk))*mol

###################### Integrando en la columna ###############################
# [ppm](z) = (R / M_i)* T(z) * [ug/m3](z)
# columna = Sumatoria (diferencia)[ppm](z)*((diferencia)[P](z)/g)
no2_mole_m2 = np.zeros((no2.Time.shape[0],no2.shape[2],no2.shape[3]))
o3_mole_m2 = np.zeros((o3.Time.shape[0],o3.shape[2],o3.shape[3]))
co_mole_m2 = np.zeros((co.Time.shape[0],co.shape[2],co.shape[3]))
hcho_mole_m2 = np.zeros((hcho.Time.shape[0],hcho.shape[2],hcho.shape[3]))

no2_cm2 = np.zeros((no2.Time.shape[0],no2.shape[2],no2.shape[3]))
o3_du = np.zeros((o3.Time.shape[0],o3.shape[2],o3.shape[3]))
co_cm2 = np.zeros((co.Time.shape[0],co.shape[2],co.shape[3]))
hcho_cm2 = np.zeros((hcho.Time.shape[0],hcho.shape[2],hcho.shape[3]))

print("Integrando NO2 en la columna vertical")
del_ppm = np.zeros((no2.shape[2],no2.shape[3]))
del_press = np.zeros((no2.shape[2],no2.shape[3]))
for t in range(no2.Time.shape[0]):
    for l in range(no2.shape[1]-1):
        del_ppm[:,:] = (R/(46*(10**-3)*9.8))*(((tk[t,l+1,:,:]*no2_u[t,l+1,:,:]*(10**-6))/pt[t,l+1,:,:]) 
                                         - ((tk[t,l,:,:]*no2_u[t,l,:,:]*(10**-6))/pt[t,l,:,:]))
        del_press[:,:] = pt[t,l+1,:,:]-pt[t,l,:,:]    ### [Pa]
        total = abs(del_ppm[:,:]*del_press[:,:]*mol/46)
        no2_mole_m2[t,:,:] = no2_mole_m2[t,:,:] + total     #### moleculas/m2
    no2_cm2[t,:,:] = no2_mole_m2[t,:,:]*(10**-4)     #### moleculas/cm2 

print("Integrando O3 en la columna vertical")
del_ppm = np.zeros((o3.shape[2],o3.shape[3]))
del_press = np.zeros((no.shape[2],no.shape[3]))
for t in range(o3.Time.shape[0]):
    for l in range(o3.shape[1]-1):
        del_ppm[:,:] = (R/(48*(10**-3)*9.8))*(((tk[t,l+1,:,:]*o3_u[t,l+1,:,:]*(10**-6))/pt[t,l+1,:,:]) 
                                         - ((tk[t,l,:,:]*o3_u[t,l,:,:]*(10**-6))/pt[t,l,:,:]))
        del_press[:,:] = pt[t,l+1,:,:]-pt[t,l,:,:]    ### [Pa]
        total = abs(del_ppm[:,:]*del_press[:,:]*mol/48)
        o3_mole_m2[t,:,:] = o3_mole_m2[t,:,:] + total     #### moleculas/m2
    o3_du[t,:,:] = o3_mole_m2[t,:,:]/(2.69*(10**20))     #### du 

print("Integrando CO en la columna vertical")
del_ppm = np.zeros((co.shape[2],co.shape[3]))
del_press = np.zeros((co.shape[2],co.shape[3]))
for t in range(co.Time.shape[0]):
    for l in range(co.shape[1]-1):
        del_ppm[:,:] = (R/(28*(10**-3)*9.8))*(((tk[t,l+1,:,:]*co_u[t,l+1,:,:]*(10**-6))/pt[t,l+1,:,:]) 
                                         - ((tk[t,l,:,:]*co_u[t,l,:,:]*(10**-6))/pt[t,l,:,:]))
        del_press[:,:] = pt[t,l+1,:,:]-pt[t,l,:,:]    ### [Pa]
        total = abs(del_ppm[:,:]*del_press[:,:]*mol/28)
        co_mole_m2[t,:,:] = co_mole_m2[t,:,:] + total     #### moleculas/m2
    co_cm2[t,:,:] = co_mole_m2[t,:,:]*(10**-4)     #### moleculas/cm2 

print("Integrando HCHO en la columna vertical")
del_ppm = np.zeros((hcho.shape[2],hcho.shape[3]))
del_press = np.zeros((hcho.shape[2],hcho.shape[3]))
for t in range(hcho.Time.shape[0]):
    for l in range(hcho.shape[1]-1):
        del_ppm[:,:] = (R/(30.03*(10**-3)*9.8))*(((tk[t,l+1,:,:]*hcho_u[t,l+1,:,:]*(10**-6))/pt[t,l+1,:,:]) 
                                         - ((tk[t,l,:,:]*hcho_u[t,l,:,:]*(10**-6))/pt[t,l,:,:]))
        del_press[:,:] = pt[t,l+1,:,:]-pt[t,l,:,:]    ### [Pa]
        total = abs(del_ppm[:,:]*del_press[:,:]*mol/30.03)
        hcho_mole_m2[t,:,:] = hcho_mole_m2[t,:,:] + total     #### moleculas/m2
    hcho_cm2[t,:,:] = hcho_mole_m2[t,:,:]*(10**-4)     #### moleculas/cm2 

###############################################################################
import basemap_integrado as bm_inter

print("Ploteando NO2")
no2_cm2_1 = no2_cm2*(10**-16)
unidad = str("$x 10^{16} moleculas/cm^{2}$")
descripcion = str("NO2")
tiempo = no2.Time
bm_inter.plot_map(output,wrf_file,no2_cm2_1,u10,v10,lat,lon,tiempo,descripcion,unidad)

print("Ploteando O3")
unidad = str("$DU$")
descripcion = str("O3")
tiempo = o3.Time
bm_inter.plot_map(output,wrf_file,o3_du,u10,v10,lat,lon,tiempo,descripcion,unidad)

print("Ploteando CO")
co_cm2_1 = co_cm2*(10**-16)
unidad = str("$x 10^{16} moleculas/cm^{2}$")
descripcion = "CO"
tiempo = co.Time
bm_inter.plot_map(output,wrf_file,co_cm2_1,u10,v10,lat,lon,tiempo,descripcion,unidad)

print("Ploteando HCHO")
hcho_cm2_1 = hcho_cm2*(10**-16)
unidad = str("$x 10^{16} moleculas/cm^{2}$")
descripcion = "HCHO"
tiempo = hcho.Time
bm_inter.plot_map(output,wrf_file,hcho_cm2_1,u10,v10,lat,lon,tiempo,descripcion,unidad)
