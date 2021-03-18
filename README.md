# Satellite-WRF-Model
This project is create to work with satellite and model data. The programming language is python.

## MODIS
### 1. Read a MODIS image:
* To read only a MODIS image you can use the [read_image.py](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/MODIS/AOD/read_modis.py) script. This script is a easy script.
#### How to use?  
* First you need install some libraries:
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

``` python
#get latitude and longitude
lat=hdf.select('Latitude')
latitude=lat[:]
min_lat=latitude.min()
max_lat=latitude.max()
lon=hdf.select('Longitude')
longitude=lon[:]
min_lon=longitude.min()
max_lon=longitude.max()

#get scale factor for AOD SDS
attributes=sds_obj.attributes()
scale_factor=attributes['scale_factor']

#get valid range for AOD SDS
range=sds_obj.getrange()
min_range=min(range)
max_range=max(range)

#create a matrix
dim2d=latitude.shape
data=np.zeros((dim2d[0],dim2d[1]))

#get SDS data
data= sds_obj.get()

#get data within valid range
valid_data=data.ravel()
valid_data=[x for x in valid_data if x>=min_range]
valid_data=[x for x in valid_data if x<=max_range]
valid_data=np.asarray(valid_data)

#scale the valid data
valid_data=valid_data*scale_factor

#find the average
average=sum(valid_data)/len(valid_data)

#find the standard deviation
stdev=np.std(valid_data)
attrs = sds_obj.attributes(full=1)
fillvalue=attrs['_FillValue']

# fillvalue[0] is the attribute value (-9999)
fv = fillvalue[0]

#turn fillvalues to NaN
data=data.astype(float)
data[data == fv] = np.nan
data = np.ma.masked_array(data, np.isnan(data))
data_final = data * scale_factor

######## ploting the image 
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

# create colorbar
cbar = m.colorbar(t, location='bottom', pad="15%")

# label colorboar
cbar.set_label('AOD', fontsize=12)

# title the plot
long_name = attributes['long_name']
time = FILE_NAME[18:20]+':'+FILE_NAME[20:22]
plt.title(long_name +'\n'+ 'Time: '+time)
fig = plt.gcf()
wsize = 300
hsize = 200
#fig.savefig(output+'image_modis.png',dpi=500)
# Show the plot window.
plt.show()
```
* sss
![Alt text](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/MODIS/AOD/figures/Figure%202021-03-18%20124501.png)

## OMI
