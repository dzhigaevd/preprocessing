# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 15:28:20 2020

@author: Admin
"""

#!/usr/bin/env python
import h5py 
import numpy as np 
import math
from matplotlib import pyplot as plt

pre_path     = "C:/WORK_DIRECTORY/dzhigd/CurrentProject/FIB_STO_BCDI/APS/data_calibration/Sample3/Sample3__105/Sample3C_ES.bin"
data_path    = "C:/WORK_DIRECTORY/dzhigd/CurrentProject/NanoSponges/Simulation/model_source/10nm.bin"
        
data = np.fromfile(data_path,dtype='double')

#data = np.multiply(amp, np.exp(1j * ph))

# Test conversion
data = np.reshape(data, (259,259,259))\
data = np.transpose(data, (0,2,1))

plt.figure(1)
plt.imshow(np.abs(data[48,:,:]))
plt.show()

plt.figure(2)
plt.imshow(np.log10(data[25,:,:]))
plt.show()

np.savez_compressed('NS_10nm_data_51_256_256.npz', amp= data)
