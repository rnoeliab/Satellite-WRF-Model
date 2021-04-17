## MODIS

### 1. Read a MODIS image:
* To read only a MODIS image you can use the [read_image.py](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/MODIS/AOD/read_modis.py) script. This script is a easy script.

#### How to use? 

* First you need to have Anaconda 3 and python 3.6 or more installed on your computer! For that go to [Installing Anaconda](https://github.com/rnoeliab/Installing_anaconda)

* Also, you need to install some libraries:
```
conda install -c anaconda basemap
conda install -c conda-forge pyhdf
```
* Then, you can run the first lines of the script :
```python
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt
from pyhdf.SD import SD, SDC
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm

```
* Use `INPUT_PATH` and `FILE_NAME` to add the path and archive of the MODIS image. The `SD` command is to read the MODIS image, to print the SDS datasets use `print hdf.datasets()`. In this SDS datasets, we can find the datafield name `Optical_Depth_Land_And_Ocean` if you are analyzed with 3km resolution spatial but if you are analyzing with a 10km resolution spatial, you can use  `Optical_Depth_Land_And_Ocean` or/and 'Deep_Blue_Aerosol_Optical_Depth_550_Land' algorithm.

``` python
# read a MODIS data
hdf = SD(INPUT_PATH+FILE_NAME, SDC.READ)

# Read dataset.
DATAFIELD_NAME='Optical_Depth_Land_And_Ocean'
sds_obj = hdf.select(DATAFIELD_NAME)
```
* To get the coordinates: Lat: `min_lat=latitude.min()` and `max_lat=latitude.max()` Lon: `min_lon=longitude.min()` and `max_lon=longitude.max()`
``` python
#get latitude and longitude
lat=hdf.select('Latitude')
latitude=lat[:]
lon=hdf.select('Longitude')
longitude=lon[:]
```
* To get the scale factor, the valid range for AOD SDS and create a matrix:
``` python
attributes=sds_obj.attributes()
scale_factor=attributes['scale_factor']

range=sds_obj.getrange()
min_range=min(range)
max_range=max(range)

dim2d=latitude.shape
data=np.zeros((dim2d[0],dim2d[1]))
```
* To get SDS data:
``` python
data= sds_obj.get()

attrs = sds_obj.attributes(full=1)
fillvalue=attrs['_FillValue']
fv = fillvalue[0]

data=data.astype(float)
data[data == fv] = np.nan
data = np.ma.masked_array(data, np.isnan(data))
data_final = data * scale_factor
```
* To plot the data:
``` python
fig = plt.figure(figsize=(8, 8))

###### creating colorbar
levels = MaxNLocator(nbins=15).tick_values(data_final.min(), data_final.max())
cmap = cm.get_cmap("jet",lut=10)
cmap.set_bad("w")
norm = BoundaryNorm(levels, ncolors=cmap1.N, clip=True)

m = Basemap(projection='cyl', area_thresh=10000,resolution='l' ,llcrnrlat=min_lat, urcrnrlat=max_lat,
            llcrnrlon=min_lon, urcrnrlon=max_lon)
m.drawcoastlines(linewidth=1.5)
m.drawstates(linewidth=1.5)
m.drawparallels(np.arange(-90., 120., 5), labels=[1, 0, 0, 0],fontsize=15)
m.drawmeridians(np.arange(-180., 181., 5), labels=[0, 0, 0, 1],fontsize=15)

x, y = m(longitude, latitude)
t = m.pcolormesh(x, y, data_final, cmap=cmap, norm=norm)

cbar = m.colorbar(t, location='bottom', pad="15%")

cbar.set_label('AOD', fontsize=12)

long_name = attributes['long_name']
time = FILE_NAME[18:20]+':'+FILE_NAME[20:22]
plt.title(long_name +'\n'+ 'Time: '+time)
fig = plt.gcf()
wsize = 300
hsize = 200
#fig.savefig(output+'image_modis.png',dpi=500)
plt.show()
```
* This script show the following image:

![Alt text](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/MODIS/AOD/figures/2017_06_04_Terra.png)

