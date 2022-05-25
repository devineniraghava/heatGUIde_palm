# -*- coding: utf-8 -*-
"""
Created on Tue May 24 09:56:54 2022

@author: rdevinen
"""

import datetime

import pytz

from netCDF4 import Dataset
import numpy as np
import xarray as xr
# import hvplot.xarray
import numpy.ma as ma
from PIL import Image
from numpy import asarray
import matplotlib.pyplot as plt

#%%


img = Image.open('C:/Users/rdevinen/Downloads/part1.png')

# pixel = img.load()
#%%
# for row in range(img.size[0]):
#     print(row)
#     for column in range(img.size[1]):
#         print(column)
        

from PIL import Image
from numpy import asarray
  
  
# load the image and convert into 
# numpy array
img = Image.open('C:/Users/rdevinen/Downloads/part1.png')
imgg = img.convert(mode='L')
numpydata = asarray(imgg)

print(type(numpydata))
  
#  shape
print(numpydata.shape)
  
# Below is the way of creating Pillow 
# image from our numpyarray
pilImage = Image.fromarray(numpydata)
print(type(pilImage))
  
# Let us check  image details
print(pilImage.mode)
print(pilImage.size)













#%%