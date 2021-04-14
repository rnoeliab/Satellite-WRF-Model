# Comparison between WRF-Chem and MODIS sensor
* Before extracting the AOD data from the WRF-Chem model, we let's first to create a "txt" where the  AOD data name from the modis sensor in "hdf" format are found. 
```python
import os
satellite_path = '../DATA/SP/3K/2017/'
satellite_txt = sorted(os.listdir(satellite_path))
with open('../DATA/SP/3K/2017_list_modis.txt', 'w') as f:
     for item in satellite_txt:
         f.write("%s\n" % item)
```
## Extract the AOD data from the WRF-Chem model 

## Regridding the AOD data from the MODIS sensor to the same resolution as the WRF-Chem model 
