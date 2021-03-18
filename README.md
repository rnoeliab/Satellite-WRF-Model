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
* Then you can run the first lines of the script :
```python
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt
from pyhdf.SD import SD, SDC
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm

```
* fddd
![Alt text](https://github.com/rnoeliab/Satellite-WRF-Model/blob/master/MODIS/AOD/figures/Figure%202021-03-18%20124501.png)

## OMI
