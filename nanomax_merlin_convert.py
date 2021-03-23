#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 13:08:32 2020

@author: dzhigd
"""
import os
import hdf5plugin
import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as scio
from tifffile import imsave

def read_masterh5_nanomax(path):
    h5file = h5py.File(path,'r')
    command = str(h5file['entry']['description'][()])[3:-2] # Reading only useful symbols
    motor_positions = {
            # Detector positions
            "delta": h5file['entry']['snapshot']['delta'][()],
            "gamma": h5file['entry']['snapshot']['gamma'][()],           
            "gonphi": h5file['entry']['snapshot']['gonphi'][()],
            "gontheta": h5file['entry']['snapshot']['gontheta'][()],
            "radius": h5file['entry']['snapshot']['radius'][()],
            "energy": h5file['entry']['snapshot']['energy'][()]
            }
    scan_keys = list(h5file['entry']['measurement'].keys())
    
    # Check which motor was scanned
    if 'gonphi' in scan_keys:        
        scan_info = {'gonphi': h5file['entry']['measurement']['gonphi'][()]}  
    elif 'gontheta' in scan_keys:
        scan_info = {'gontheta': h5file['entry']['measurement']['gontheta'][()]}  
    return command, motor_positions, scan_info

def read_data_merlin(data_path,roi):
    h5file = h5py.File(data_path, 'r')
    data = h5file['entry']['measurement']['Merlin']['data'][:,roi[0]:roi[1],roi[2]:roi[3]]       
    return data

def save_bin_data_nanomax(save_path,data):
    fid = open(save_path,'w+b')
    fid.write(data)
    fid.close()
    
def read_data_xspress3(data_path,roi):
    h5file = h5py.File(data_path, 'r')
    entry_keys = list(h5file['/'])
    ii = 0
    data = np.zeros((91, len(entry_keys)))
    for key in entry_keys:
        data_temp = h5file[key]['measurement']['xspress3']['data'][:,3,roi[0]:roi[1]]
        data[:,ii] = np.sum(data_temp,1)       
        ii = ii+1
    return data

#inputs
year = "2020"                                                                  #The year for the experiemnt
beamtimeID ="2020101408"                                                       #The beamtimeID
sample_name = r"0002_sample_P246_AB"                                           #The name for the p10 newfile
scan = np.linspace(109,128,20)                                                 #The scan numbers, can be the list
xrf_roi = [870,970]
xrd_roi = [0,230,250,400]

for ii in range(0,len(scan)):
    print(('--Processing scan %d--')%(int(scan[ii])))
    #determing the paths"
    path=r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/raw/%s"%(beamtimeID,sample_name)
    
#    h5_path=os.path.join(path, "%06d.h5"%(int(scan[ii])))
    data_xrd_path=os.path.join(path, "scan_%06d_merlin.hdf5")%(int(scan[ii]))
    data_xrf_path=os.path.join(path, "scan_%06d_xspress3.hdf5")%(int(scan[ii]))
    
    # Load the data ###########################################################
    data_xrd =  read_data_merlin(data_xrd_path,xrd_roi)
    data_xrf =  read_data_xspress3(data_xrf_path,xrf_roi)
    ###########################################################################
    
    # Save the xrd data #######################################################
#    save_path = os.path.join(path,"converted", "scan_%06d_merlin_%d_%d_%d.bin")%(int(scan[ii]),(data_xrd.shape[0]),(data_xrd.shape[1]),(data_xrd.shape[2]))      #The path to save the results
#    save_bin_data_nanomax(save_path, data_xrd)
#    print(('Saved scan %d to %s')%(int(scan[ii]),save_path))
    
    save_path = os.path.join(path,"converted", ("scan_%06d_merlin.mat")%(int(scan[ii])))      #The path to save the results
    data_dic = {"data":data_xrd}
    scio.savemat(save_path, data_dic)    
    print(('Saved scan %d to %s')%(int(scan[ii]),save_path))
    ###########################################################################

    # Save the xrf data #######################################################
#    save_path = os.path.join(path,"converted","scan_%06d_xspress3_%d_%d_%d.bin")%(int(scan[ii]),(data.shape[0]),(data.shape[1]),(data.shape[2]))      #The path to save the results
#    save_bin_data_nanomax(save_path, data_xrf)
#    print(('Saved scan %d to %s')%(int(scan[ii]),save_path))
    
#    save_path = os.path.join(path,"converted", ("scan_%06d_xspress3.mat")%(int(scan[ii])))      #The path to save the results
#    data_dic = {"data":data_xrf}
#    scio.savemat(save_path, data_dic)    
#    print(('Saved scan %d to %s')%(int(scan[ii]),save_path))
    
    data_xrf = np.float32(data_xrf/np.max(data_xrf))
    
    imsave(os.path.join(path,"converted", ("scan_%06d_xspress3.tif")%(int(scan[ii]))), data_xrf)
    ###########################################################################

    # NW correct the drift in 2D ##############################################
    