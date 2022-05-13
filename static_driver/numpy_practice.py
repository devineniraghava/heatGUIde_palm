#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:59:59 2022

@author: rdevinen
"""

import numpy as np

import datetime
import pytz

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')
#%%


x = np.array(range(26))

y = -(1/1.83) * pow(x,1/0.6) + 117

y = np.around(y).astype(int)
# =============================================================================
# 
# plt.plot(x,y)
# 
# 
# # Add a title
# plt.title('My first Plot with Python')
# 
# # Add X and y Label
# plt.xlabel('x axis')
# plt.ylabel('y axis')
# 
# # Add a grid
# plt.grid(alpha=.4,linestyle='--')
# 
# # Add a Legend
# plt.legend()
# 
# # Show the plot
# plt.show()
# 
# =============================================================================


a = np.zeros((200,200))

for i, j in zip(x,y):
    a[i,:j+1] = 1




















#%%