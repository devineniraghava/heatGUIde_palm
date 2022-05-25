# -*- coding: utf-8 -*-
"""
Created on Mon May 23 13:28:14 2022

@author: rdevinen
"""

import numpy as np
from PIL import Image

#%%

an_image = Image.open("C:/Users/rdevinen/Downloads/cat1.jpg")


image_sequence = an_image.getdata()
image_array = np.array(image_sequence)
numpydata = np.asarray(image_sequence)

print(image_sequence.size)
print(image_sequence.mode)

#%%
from PIL import Image
from numpy import asarray
  
  
# load the image and convert into 
# numpy array
img = Image.open('C:/Users/rdevinen/Downloads/cat1.jpg')
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

arr1 = np.where(numpydata > 100, 1, 0)