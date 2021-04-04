## OMI

### 1. Read a OMI image:

* To read only a OMI image you can use the read_image.py script. This script is a easy script.

How to use?
* First you need to have Anaconda 3 and python 3.6 or more installed on your computer! For that go to Installing Anaconda
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
  
