# Comparison between WRF-Chem and MODIS sensor
* Before extracting the AOD data from the WRF-Chem model, we let's first to create a "txt" where the  AOD data name from the modis sensor in "hdf" format are found. 
```python
import os
satellite_path = '../DATA/SP/3K/2017/'  ## path where the hdf files are found
satellite_txt = sorted(os.listdir(satellite_path))  
with open('../DATA/SP/3K/2017_list_modis.txt', 'w') as f:
     for item in satellite_txt:
         f.write("%s\n" % item)
```
## Extract the AOD data from the WRF-Chem model 
* The next step would be to identify the start and end date of our comparison, in this case, both the "hdf" files and the outputs of the WRF-Chem model must match the start and end date. 
* Now we are going to explain script [1.Extract_aod_wrf_to_nc.py](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/MODIS/AOD/WRF-Chem_Modis/1.Extract_aod_wrf_to_nc.py) a bit:
* Install the libraries we need. For that, we first need to have created a project through anaconda (see [Installing_anaconda](https://github.com/rnoeliab/Installing_anaconda)):
```
conda create --name py37 python=3.7 matplotlib basemap gdal 
conda activate py37
```
* As we need to install some libraries, we are going to write the following command in the terminal:
```
 conda install -c anaconda netcdf4 
 conda install -c anaconda xarray 
 conda install -c conda-forge wrf-python 
 conda install -c anaconda pandas 
 conda install -c unidata metpy 
```
* Now, we can run this part of script:
```python
from __future__ import print_function
from glob import glob
from netCDF4 import Dataset
import xarray as xr
from wrf import to_np, getvar, ALL_TIMES
import numpy as np
from datetime import datetime
import pandas as pd
from metpy.calc import wind_components
from metpy.units import units
```
* Then, we let's go to read the created "txt" file and read the outputs from WRF-Chem model:
```python
f = open('../DATA/SP/3K/2017_list_modis.txt')
satellite_txt = f.readlines()

wrf_file = sorted(glob("../wrfout_d02_*"))  # archivos wrfout
output = "../wrf_to_nc/"
```
* Find the intersection between satellite and WRF times (all in julian days):
```python
djulian = [datetime.strptime(i[-19:-9], '%Y-%m-%d').timetuple().tm_yday for i in wrf_file]  # JD = tm_yday
listdir = []
print(djulian[0], djulian[len(djulian)-1])

for n in range(len(djulian)):
    name = '_3K.A2017'+str(djulian[n])+'.'
    mod_wrf = sorted(list(filter(lambda a: str(name) in a, satellite_txt)))
    listdir += mod_wrf
with open('../wrf_to_nc/list_mod_wrf.txt', 'w') as f:
    for item in listdir:
        f.write(item)
```
* Read WRF data, only the variable of interest
```python
ncfiles = [Dataset(x) for x in wrf_file]  # read all wrfout
# print ncfiles

print("winds")
wind = getvar(ncfiles, 'uvmet10_wspd_wdir', timeidx=ALL_TIMES, method='cat')
ws = wind.sel(wspd_wdir='wspd')
wd = wind.sel(wspd_wdir='wdir')
##################### change m/s to knots ###############################
ws10 = (ws.values*units('m/s')).to('knots')
wd10 = wd.values*units.degree
##################### calculate the U and V components ##################
u10, v10 = wind_components(ws10, wd10)

############################### AOD #####################################
tau1 = getvar(ncfiles, "TAUAER1", timeidx=ALL_TIMES, method='cat')  # AOD in 300nm
tau2 = getvar(ncfiles, "TAUAER2", timeidx=ALL_TIMES, method='cat')  # AOD in 400nm
tau3 = getvar(ncfiles, "TAUAER3", timeidx=ALL_TIMES, method='cat')  # AOD in 600nm
tau4 = getvar(ncfiles, "TAUAER4", timeidx=ALL_TIMES, method='cat')  # AOD in 1000nm
```
* Calculate AOD in 550nm :
```python
angstrom = np.zeros((tau1.shape[0], tau1.shape[1], tau1.shape[2], tau1.shape[3]))
aod_550 = np.zeros((tau1.shape[0], tau1.shape[1], tau1.shape[2], tau1.shape[3]))

start = datetime.now()
for t in range(tau1.shape[0]):
    print("t=", t)
    for l in range(tau1.shape[1]):
        angstrom[t, l, :, :] = np.log(to_np(tau1[t, l, :, :])/to_np(tau4[t, l, :, :]))/(np.log((1000./300.)))
        aod_550[t, l, :, :] = to_np(tau2[t, l, :, :])*np.power((550./400.), -1*angstrom[t, l, :, :])
print(datetime.now()-start)

aod_550[np.isnan(aod_550)] = 0.0
```
* Integrating in the vertical column
```python
aod_550_col = np.zeros((tau1.shape[0], tau1.shape[2], tau1.shape[3]))
for t in range(tau1.shape[0]):
    for l in range(tau1.shape[1]):
        aod_550_col[t, :, :] = aod_550_col[t, :, :] + aod_550[t, l, :, :]
aod_550_col[aod_550_col < 0.0] = np.nan

date_wrf = pd.DataFrame({'date': tau1.Time.values})  # Times in UTC

aod_col = np.zeros((len(listdir), aod_550_col.shape[1], aod_550_col.shape[2]))
u = np.zeros((len(listdir), aod_550_col.shape[1], aod_550_col.shape[2]))
v = np.zeros((len(listdir), aod_550_col.shape[1], aod_550_col.shape[2]))
```
* Find the intersection between satellite and WRF times (all in julian days)
```python
ddt = []
dtt = []
dttt = []

count = 0
for n in range(len(djulian)):
    hours_wrf = date_wrf['date'].dt.strftime('%H:%M')[n*24:(n*24)+24]  # 00hrs - 24hrs
    name = '_3K.A2017'+str(djulian[n])+'.'
    mod_wrf = sorted(list(filter(lambda a: str(name) in a, satellite_txt)))
    ################# Selecting only the hours of interest ####################
    for j in mod_wrf:
        hours_mod = j[18:20] + ':00'
        h_mod_wrf = hours_wrf[j[18:20] + ':00' == hours_wrf].index.values[0]
        h_mod_wrf_1 = hours_wrf[str(int(j[18:20])+1) + ':00' == hours_wrf].index.values[0]
        print(h_mod_wrf, h_mod_wrf_1)
        ddt.append(tau4.coords['XTIME'].values[h_mod_wrf])
        dtt.append(tau4.coords['Time'].values[h_mod_wrf])
        dttt.append(tau4.coords['datetime'].values[h_mod_wrf])
        aod_col[count, :, :] = (aod_550_col[h_mod_wrf, :, :] + aod_550_col[h_mod_wrf_1, :, :])/2
        v[count, :, :] = (v10[h_mod_wrf, :, :] + v10[h_mod_wrf_1, :, :])/2
        u[count, :, :] = (u10[h_mod_wrf, :, :] + u10[h_mod_wrf_1, :, :])/2
        count = count + 1
```
* Save the information in other netcdf with xarray 
```python
da = xr.Dataset(
    data_vars = dict(
    aod_055=(['Time', 'south_north', 'west_east'],aod_col),
    uu=(['Time', 'south_north', 'west_east'],u),
    vv=(['Time', 'south_north', 'west_east'],v),
    ),
    coords=dict(
        XLONG=(['south_north', 'west_east'], tau4.coords['XLONG'].values),
        XLAT=(['south_north', 'west_east'], tau4.coords['XLAT'].values),
        XTIME=(['Time'], np.array(ddt)),
        Time=('Time', dtt),
        datetime=('Time', dttt)
    ),
    attrs=tau4.attrs)
da.attrs['description'] = '550nm optical thickness'
del da.attrs['projection']

da.to_netcdf(output+"june_aod.avg.column.550_p1.nc")
```

## Regridding the AOD data from the MODIS sensor to the same resolution as the WRF-Chem model 
* After extracting the variable of interest in a "netCDF" format, we are going to read it and regridding it to the same resolution as the MODIS sensor [2.Regridding_modis_wrf_aod.py](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/MODIS/AOD/WRF-Chem_Modis/2.Regridding_modis_wrf_aod.py).
* Therefore, you need to install a "xesmf" library: 
```
 conda install -c conda-forge xesmf 
```
* Later, you can run these libraries on the spyder platform: 
```python
import os
import gdal
import time
import calendar
import numpy as np
import xesmf as xe
import xarray as xr
from PIL import Image
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Polygon
```
* Now, we are going to read the directories where the data of the WRF-Chem model and the MODIS sensor are located. 
```python
dire_mod = '../DATA/SP/3K/2017/'
dire_wrf = '../DATA/'
input_mod = open(dire_wrf+'list_mod_wrf.txt').readlines()
output = "..s/DATA/SP/regridded_2017/"
```
* Then, we are going to read the data from the WRF-Chem model and organize it properly:
```python
ds_disk = xr.open_dataarray(dire_wrf+"june_aod.avg.column.550_p2.nc")
ds = xr.DataArray(
    data=ds_disk.variable.values[0],
    dims = ("y","x"),
    coords = dict(
        lon = (["y","x"],ds_disk.XLONG.values),
        lat = (["y","x"],ds_disk.XLAT.values),),)  
```
* And read the data from MODIS sensor. The "var_name = mod04:Optical_Depth_Land_And_Ocean" means the name of the variable of interest.
* From command "os.system ('gdalwarp -of GTIFF -tps -t_srs EPSG: 4326 HDF4_EOS: EOS_SWATH:" {0} ": {1} teste.hdf'.format (FILE_NAME, var_name))" to "im1 = Image.open(OUT+os.sep+'teste.hdf')" is used to reproject the MODIS sensor dice from sinusoidal to Cartesian. 
```python
var_name = 'mod04:Optical_Depth_Land_And_Ocean'
for enu,n in enumerate(input_mod):
    FILE_NAME = dire_mod+n.rstrip()
    print(FILE_NAME)
    os.system('gdalwarp -of GTIFF -tps -t_srs EPSG:4326 HDF4_EOS:EOS_SWATH:"{0}":{1} teste.hdf'.format(FILE_NAME, var_name))    
    OUT = os.getcwd()
    ifile = OUT+os.sep+'teste.hdf' 
    gc = gdal.Open(ifile)
    metadata = gc.GetMetadata()
    width = gc.RasterXSize
    height = gc.RasterYSize 
    
    if width <= 1000 and height <=1000:
        print(width, height)
        gt = gc.GetGeoTransform()
        bands = gc.RasterCount
        projInfo = gc.GetProjection()
        valid_range = [-100, 5000]
        fv = -9999
        scale_factor=0.00100000004749745
        minx = gt[0]
        maxx = gt[0] + width*gt[1] + height*gt[2]
        miny = gt[3] + width*gt[4] + height*gt[5]
        maxy = gt[3]     
#        print ("la latitude esta en el rango:" + str(miny) + ":" + str(maxy))
#        print ("la longitude esta en el rango:"+ str(minx) + ":" + str(maxx))
    
        lat = np.linspace(miny, maxy, height)[::-1]
        lon = np.linspace(minx, maxx, width)
        x,y = np.meshgrid(lon,lat)

        im1 = Image.open(OUT+os.sep+'teste.hdf')
        im1 = np.array(im1)
        data = np.float64(im1)            
```
* After reprojecting the MODIS sensor data, we will select only the data for our area of interest: 
```python
################### AREA STUDY
x0 = -49.0
x1 = -44.0
y0 = -21.0
y1 = -26.0
################    
def find_nearest_x(longitude, point_x):
    longitude = np.asarray(longitude)
    idx = (np.abs(longitude - point_x)).argmin()
    return idx

def find_nearest_y(latitude, point_y):
    latitude = np.asarray(latitude)
    dist = (np.abs(latitude - point_y))
    idy = np.where(dist ==dist.min())[0][1]
    return idy

if ((y0 > miny and y0 < maxy) or (y1 > miny and y1 < maxy)):
  if ((x0 > minx and x0 < maxx) or (x1 > minx and x1 < maxx)):
      posx0 = find_nearest_x(x,x0)
      posx1 = find_nearest_x(x,x1)
      posy0 = find_nearest_y(y,y0)
      posy1 = find_nearest_y(y,y1)

      y = y[posy0:posy1+1,posx0:posx1+1]
      x = x[posy0:posy1+1,posx0:posx1+1]

      data_final = data_f[posy0:posy1+1,posx0:posx1+1]
      a,b = data_final.shape
      print(a,b)
```
* The next step would be to convert UTC time to local time with the command "time.localtime ()":
```python
date = metadata['RANGEBEGINNINGDATE']
times = metadata['RANGEBEGINNINGTIME'][0:8]
time_tuple = time.strptime(date+' '+times, "%Y-%m-%d %H:%M:%S")
t = calendar.timegm(time_tuple)
aa = time.localtime(t)
time_local = time.strftime('%Y-%m-%d %H:%M:%S', aa)
```
* Finally, the necessary data is saved: Data, latitude and longitude.
```python
da = xr.DataArray(
    data=data_final,
    dims = ("y","x"),
    coords = dict(
        lon = (["y","x"],x),
        lat = (["y","x"],y),),)  
```
* And regridding the Data and save it in a "netCDF" format:
```python
regridder = xe.Regridder(da, ds, 'nearest_s2d')
regridder.clean_weight_file()
print(regridder)

data_out = regridder(da.values) 

df = xr.DataArray(
    data=data_out,
    dims = ('y', 'x'),
    coords = dict(
        lon = (['y', 'x'], ds_disk.XLONG.values),
        lat = (['y', 'x'], ds_disk.XLAT.values),
        ),
    attrs = dict(
        description = '550nm optical thickness',
        units = 'None',
        time = str(time_local),
        time_units = 'Local Time'
        ),)
df.to_netcdf(output+str(file_names)+"_regrid.nc")
```
## Comparing the model WRF-Chem with the MODIS sensor
* The last step would be to compare the regridded data from the WRF-Chem model with the data from the MODIS sensor.
* To run the [3.TandA_wrf_compared.py](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/MODIS/AOD/WRF-Chem_Modis/3.TandA_wrf_compared.py) script, it is necessary to run script [basemap_modis.py](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/MODIS/AOD/WRF-Chem_Modis/basemap_modis.py) first.
```python
dire_mod = "../DATA/SP/regridded_2017/"
dire_wrf = '../DATA/'
output = "../results/plot_wrf_modis/"

########################### READ WRF IMAGES ###############################
ds_disk = xr.open_dataset(dire_wrf+"june_aod.avg.column.550_p2.nc")
aod_wrf = ds_disk["aod_055"].values
.
.
.
########################### READ MODIS IMAGES #############################
files1 = sorted(glob.glob(dire_mod+"*_regrid.nc"))  ## archivos wrfout
for i in range(len(files1)):
    if i == (len(files1)-1):
        break;      
    else:
        if files1[i][-50:-33] == files1[i+1][-50:-33]:  #### to the same time
            listdirr.append(files1[i])          
            if i == (len(files1)-2):
                listdirr.append(files1[len(files1)-1])
                mylist.append(listdirr)               
        else:
            listdirr.append(files1[i])
            mylist.append(listdirr)
            listdirr = []
.
.
.                        
########################### TERRA Y AQUA ################################      
######################## mean by periodo ################################
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
```
![Alt text](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/MODIS/AOD/figures/2017_06_04_Terra.png)

