import os
import gdal
import numpy as np
from PIL import Image
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Polygon
import xesmf as xe
import time
import xarray as xr
import calendar

dire_mod = '../DATA/SP/3K/2017/'
dire_wrf = '../DATA/'
input_mod = open(dire_wrf+'june_list_mod_wrf.txt').readlines()
output = "..s/DATA/SP/regridded_2017/"

ds_disk = xr.open_dataarray(dire_wrf+"june_aod.avg.column.550_p2.nc")
ds = xr.DataArray(
    data=ds_disk.variable.values[0],
    dims = ("y","x"),
    coords = dict(
        lon = (["y","x"],ds_disk.XLONG.values),
        lat = (["y","x"],ds_disk.XLAT.values),),)  

x00,x01 = ds_disk.XLONG.values[0,0],ds_disk.XLONG.values[0,-1]
x10,x11 = ds_disk.XLONG.values[-1,0],ds_disk.XLONG.values[-1,-1]
y00,y01 = ds_disk.XLAT.values[-1,0],ds_disk.XLAT.values[-1,-1]
y10,y11 = ds_disk.XLAT.values[0,0],ds_disk.XLAT.values[0,-1]

def find_nearest_x(longitude, point_x):
    longitude = np.asarray(longitude)
    idx = (np.abs(longitude - point_x)).argmin()
    return idx


def find_nearest_y(latitude, point_y):
    latitude = np.asarray(latitude)
    dist = (np.abs(latitude - point_y))
    idy = np.where(dist ==dist.min())[0][1]
    return idy

################### AREA STUDY
x0 = -49.0
x1 = -44.0
y0 = -21.0
y1 = -26.0
################    
            
levels = MaxNLocator(nbins=18).tick_values(0,1.2)
cmap = cm.get_cmap("Spectral_r",lut=25)
cmap.set_bad("w")
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

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
        data[data == fv] = np.nan
        data[data <= 0.0] = np.nan
#        data = np.ma.masked_array(data, np.isnan(data))
        data_f = data * scale_factor
        os.remove(OUT+os.sep+'teste.hdf')

        if ((y0 > miny and y0 < maxy) or (y1 > miny and y1 < maxy)):
            if ((x0 > minx and x0 < maxx) or (x1 > minx and x1 < maxx)):
                posx0 = find_nearest_x(x,x0)
                posx1 = find_nearest_x(x,x1)
                posy0 = find_nearest_y(y,y0)
                posy1 = find_nearest_y(y,y1)
                
                y = y[posy0:posy1+1,posx0:posx1+1]
                x = x[posy0:posy1+1,posx0:posx1+1]
                
                data_final = data_f[posy0:posy1+1,posx0:posx1+1]
#                data_final[data_f <= 0.000001] = np.nan
                # mask = np.ma.getmask(data_final)
                a,b = data_final.shape
                print(a,b)
#                    print data_final.shape
#                    print data_final.max(), data_final.min()  
                if a<=1 or b<=1:
                    print("Area muy pequena")
                    continue 
                if  np.nansum(data_final) != 0.0:
                    file_names = FILE_NAME[-44:-4]
                    print(file_names)
                    date = metadata['RANGEBEGINNINGDATE']
                    times = metadata['RANGEBEGINNINGTIME'][0:8]
                    time_tuple = time.strptime(date+' '+times, "%Y-%m-%d %H:%M:%S")
                    t = calendar.timegm(time_tuple)
                    aa = time.localtime(t)
                    time_local = time.strftime('%Y-%m-%d %H:%M:%S', aa)
    
                    da = xr.DataArray(
                        data=data_final,
                        dims = ("y","x"),
                        coords = dict(
                            lon = (["y","x"],x),
                            lat = (["y","x"],y),),)                
                      
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

                    # fig = plt.figure(figsize=(22, 11))
                    # rect = fig.patch
                    # rect.set_facecolor('lightgoldenrodyellow')
                    # ax0 = fig.add_subplot(111, frame_on=False)
                    # ax0.set_xticks([])
                    # ax0.set_yticks([])
                    # ax = fig.add_subplot(121)
                    # for axis in ['top','bottom','left','right']:
                    #     ax.spines[axis].set_linewidth(3.0)
                    # m = Basemap(projection='cyl', resolution='h', llcrnrlat=-26.0, urcrnrlat=-21.0,
                    #             llcrnrlon=-49.0, urcrnrlon=-44.0)
                    # m.drawcoastlines(linewidth=1.5)
                    # m.drawstates(linewidth=1.5)    
                    # m.drawparallels(np.arange(-90., 120., 2), labels=[1, 0, 0, 0],fontsize=18)
                    # m.drawmeridians(np.arange(-180., 181., 2), labels=[0, 0, 0, 1],fontsize=18) 
                    # trend=m.pcolormesh(x, y, data_final, cmap=cmap, norm = norm)
                    # poly = Polygon([m(x00,y00),m(x01,y01),m(x11,y11),m(x10,y10)],facecolor='none',edgecolor='red',linewidth=5)
                    # plt.gca().add_patch(poly)       	           
                    # ax.set_title("Sao Paulo Metropolitan Region" + '\n' +
                    #               "AOD  from satellite "+'\n'+ 
                    #               "for "+ str(time_local)+ " (Local Time)", fontsize=20)
                    
                    # ax = fig.add_subplot(122)
                    # for axis in ['top','bottom','left','right']:
                    #     ax.spines[axis].set_linewidth(3.0)
                    # m = Basemap(projection='cyl', resolution='h', llcrnrlat=-25.2, urcrnrlat=-21.5,
                    #             llcrnrlon=-49.0, urcrnrlon=-44.3)
                    # m.drawcoastlines(linewidth=1.5)
                    # m.drawstates(linewidth=1.5)    
                    # m.drawparallels(np.arange(-90., 120., 1), labels=[1, 0, 0, 0],fontsize=18)
                    # m.drawmeridians(np.arange(-180., 181., 1), labels=[0, 0, 0, 1],fontsize=18) 
                    # x2,y2 = ds_disk.XLONG.values,ds_disk.XLAT.values
                    # trend2 = m.pcolormesh(x2,y2, data_out, cmap=cmap, norm = norm)
                    # poly = Polygon([m(x00,y00),m(x01,y01),m(x11,y11),m(x10,y10)],facecolor='none',edgecolor='red',linewidth=5)
                    # plt.gca().add_patch(poly)       	           
                
                    # cbar = m.colorbar(trend, location='right', pad="5%", ticks=levels)
                    # cbar.set_label('None', fontsize=19)
                    
                    # ax.set_title("Sao Paulo Metropolitan Region" + '\n' +
                    #               "AOD  from Regridded satellite "+'\n'+ 
                    #               "for "+ str(time_local)+ " (Local Time)", fontsize=20)
                    
                    # plt.show()

