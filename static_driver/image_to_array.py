# -*- coding: utf-8 -*-
"""
Created on Mon May 23 13:28:14 2022

@author: rdevinen
"""

import numpy as np
from PIL import Image

from numpy import asarray
  


#%%

  
# load the image and convert into 
# numpy array
img = Image.open('/home/rdevinen/Documents/GitHub/heatGUIde_palm/static_driver/files/part1.png')
imgg = img.convert(mode='RGB') # convert image to array in RGB 
numpydata = np.asarray(imgg) # the array is an 3D array with 3values in Z axis in the order or R G B



  
#shape
print(numpydata.shape)
r_arr, g_arr, b_arr = numpydata[:,:,0], numpydata[:,:,1], numpydata[:,:,2]
 

#%%

import numpy as np 
import matplotlib.pyplot as plt

plt.imshow(g_arr)
plt.legend()
# plt.show()
#%%

import matplotlib.pyplot as plt
import numpy as np

# Fixing random state for reproducibility
np.random.seed(19680801)

plt.subplot(211)
plt.imshow(np.random.random((100, 100)))
plt.subplot(212)
plt.imshow(np.random.random((100, 100)))

plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = plt.axes([0.85, 0.1, 0.075, 0.8])
# plt.colorbar(cax=cax)

plt.show()



#%%