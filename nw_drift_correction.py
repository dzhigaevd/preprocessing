#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 13:08:32 2020

@author: dmitry dzhigaev and tomas stankevic
"""

import os
import hdf5plugin
import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as scio
from tifffile import imsave
import tifffile as tf
from scipy import ndimage, interpolate
from skimage.transform import resize
import matplotlib.mlab as ml
from scipy.optimize import minimize

def get_q_coordinates(xrd_roi,data_size,motor_positions,beamline): 
    # Constants. They are needed for correct labeling of axes
    h = 4.1357e-15                                                            # Plank's constant
    c = 2.99792458e8                                                          # Speed of light in vacuum
    
    wavelength = h*c/motor_positions['energy']
    
    k = 2*np.pi/wavelength; # wave vector
    
    dq = k*2*np.arctan(beamline['detector_pitch']/(2*nanomax['radius']))      # q-space pitch at the detector plane
    
    hd,vd = np.meshgrid(np.linspace(-data_size[1]/2,(data_size[1]/2-1),data_size[1]),np.linspace(-data_size[0]/2,(data_size[0]/2-1),data_size[0]))
    
    hd = (hd+(data_size[1]/2-beamline['direct_beam'][1]))*beamline['detector_pitch'];
    vd = (vd+(data_size[0]/2-beamline['direct_beam'][0]))*beamline['detector_pitch'];
    zd = np.ones(size(vd))*beamline['radius']
    
    hd = hd[xrd_roi(2):xrd_roi[3],xrd_roi[0]:xrd_roi[1],:];
    vd = vd[scan.roi(2):scan.roi(2)+scan.roi(4),scan.roi(1):scan.roi(1)+scan.roi(3),:];
    zd = zd[scan.roi(2):scan.roi(2)+scan.roi(4),scan.roi(1):scan.roi(1)+scan.roi(3),:];

    d = [np.concatenate(hd),np.concatenate(vd),np.concatenate(zd)]
    d = np.transpose(d)

    r = np.squeeze(np.sqrt(np.sum(d**2,0)));

    hq = k*(d[1,:]/r);
    vq = k*(d[2,:]/r);
    zq = k*(1-d[3,:]/r);

    q = [hq,vq,zq];
    
    return q

# Map correction methods
def generate_ref_image(fluo_image, ref_image, scan_pixel_size, reference_pixel_size, tl_offset):
    """
    modifies reference image to match pixel size and roi of the scan area
    fluo_path - path to fluorescence image
    reference_path - path to reference image
    scan_pixel_size - [#, #] pixel size in y and x
    reference_pixel_size - [#,#] pixel size of original reference image
    tl_offset - [y, x] offset of the top left corner the scan area with respect to the reference image
    
    returns:
    referene image, matching to the scan area
    mask, zero where there should not be any objects accordint to the design
    """

    scan_size = fluo_image.shape # measured image size in pixels
    padsize = 1000
    beam_footprint = 0.08*2/scan_pixel_size[1]
      
    ref_image[ref_image>0]=255 
    ref_image[ref_image==0]=1
    ref_image_pad = np.pad(ref_image[:,:,0], pad_width=padsize,mode='constant')
    
    #scan range (um)
    scan_range_top_left = ((np.array(tl_offset))/orig_pixel_size).astype(np.int32)+np.array([padsize,padsize])
    scan_range_bottom_right = ((np.array(tl_offset) + np.array(scan_size)*np.array(scan_pixel_size))/np.array(orig_pixel_size)).astype(np.int32)+np.array([padsize,padsize])

    ref_image_cropped = ref_image_pad[scan_range_top_left[0]:scan_range_bottom_right[0],scan_range_top_left[1]:scan_range_bottom_right[1]]
    mask = ref_image_cropped.copy()
    mask[mask<1] = 0
    mask[mask>0] = 255
    mask = resize(mask, scan_size)
    
    ref_image_cropped = resize(ref_image_cropped, scan_size)
    
    smoothed = ndimage.gaussian_filter(ref_image_cropped, beam_footprint)
    return smoothed, mask

def normalize(data):
    """
    normalize data for cross correlation
    """
    return (data - np.mean(data, axis=0)) / (np.std(data, axis=0))

def transform(y, params):
    """
    polynomial transformation of x coordinate into x' and interpolation of the input data y
    params - polynomial coefficients [scaling, scaling gradient]
    returns distorted data y
    """
    x = np.linspace(0,len(y),len(y))
    # generate new x values - xprime by adding a polynomial
    if len(params)==1:
        xprime = x +  params[0]*x 
    if len(params)==2:
        xprime = x +  params[0]*x+ params[1]*x**2
    if len(params)==3:
        xprime = x +  params[0]*x+ params[1]*x**2 + params[2]*x**3
        
    if len(params)==4:
        xprime = x +  params[0]*x+ params[1]*x**2 + params[2]*x**3 + params[3]*x**4
        
    return np.interp(xprime,x,y), xprime

def cross_corr(y1, y2):
    """Calculates the cross correlation and lags without normalization.

    The definition of the discrete cross-correlation is in:
    https://www.mathworks.com/help/matlab/ref/xcorr.html

    Args:
    y1, y2: Should have the same length.

    Returns:
    max_corr: Maximum correlation without normalization.
    lag: The lag in terms of the index.
    """
    if len(y1) != len(y2):
        raise ValueError('The lengths of the inputs should be the same.')

    y1_auto_corr = np.dot(y1, y1) / len(y1)
    y2_auto_corr = np.dot(y2, y2) / len(y1)
    corr = np.correlate(y1, y2, mode='same')
    # The unbiased sample size is N - lag.
    unbiased_sample_size = np.correlate(
    np.ones(len(y1)), np.ones(len(y1)), mode='same')
    corr = corr / unbiased_sample_size / np.sqrt(y1_auto_corr * y2_auto_corr)
    shift = len(y1) // 2

    max_corr = np.max(corr)
    argmax_corr = np.argmax(corr)#[shift-5:shift+5])
    return max_corr, argmax_corr - shift

def errorf(params, y1, y2):
    """
    caucluates error function for minimization. 
    params - distortion polynomial coefficients
    y1 - data
    y2 - reference
    minimization parameter is negative max of cross-correlation
    """
    y,_ = transform(y1, params)
    maxcorr, lag = cross_corr(y, y2)
    return -maxcorr

def find_distortion(y1, y2, x0 = [0,0]):
    """
    y1 - data
    y2 - reference
    finds optimal distortion parameters based on minimization of cross-correlation function
    returns optimial distortion parameters
    """
    x0 = [0,0]
    res = minimize(errorf, x0, method='Nelder-Mead', tol=1e-8, args=(y1,y2))
    return res.x

def image_to_grid(image,positions):
    X = positions[1,:]
    Y = positions[0,:]
    xi = np.linspace(np.min(X), np.max(X), image.shape[0])
    yi = np.linspace(np.min(Y), np.max(Y), image.shape[1])
    xi, yi = np.meshgrid(xi, yi)
    zi = interpolate.griddata((X.transpose().ravel(),Y.transpose().ravel()),    
                                   image.transpose().ravel(),          
                                   (xi,yi),    
                                   method='linear').transpose()
    zi[zi==np.nan] = 0
    return zi

def grid_to_image(image,positions):
    X = positions[1,:]
    Xr = np.reshape(X,image.transpose().shape)
    xi = np.linspace(np.min(X), np.max(X), image.shape[0])
    image_int = X_shift.copy()
    i=0
    
#    plt.figure()
    for line in X_shift.transpose():
        x = Xr[i,:]
        line_int = np.interp(x,xi,line)
        image_int[:,i] = line_int
        i=i+1
        #image_int[image_int==np.nan] = 0
    return image_int
    

def align_image(fluo_image, reference_image, positions, scan_pixel_size, scaleY=False):
    fluo_image = image_to_grid(fluo_image,positions)
    fluo_aligned = fluo_image.copy() #placeholder array
    X_shift = fluo_image.copy()
    Y_shift = fluo_image.copy()
    # find y shift
    fluo_profile_x = normalize(np.diff(np.nansum(fluo_image,axis=0)))
    fluo_ref_profile_x = normalize(np.diff(np.nansum(reference_image,axis=0)))
    corr,Ylag = cross_corr(fluo_profile_x,fluo_ref_profile_x) # determine shift by xcorrelation
    fluo_image = np.roll(fluo_image,-Ylag, axis=1)
    i=0
    for line_data in fluo_image.transpose():
        line_data[np.isnan(line_data)] = 0
        line_ref = normalize(reference_image[:,i]) #normalize reference
        x = np.linspace(0,len(line_data),len(line_data)) # original x coordinates
        if np.nanmean(line_data)/np.nanstd(line_data) < 2: # check if signal to noise ratio is good (to avoid excessive random shifts of empty lines)
            line_data_norm = normalize(line_data) # normalize data
            params = find_distortion(line_data_norm,line_ref,[0,0,0,0]) # find distortion parameters
            line_data_aligned, xprime = transform(line_data, params) # transform data based on found parameters. does not include shift
            line_data_aligned_norm = normalize(line_data_aligned) # normalize transformed data
            corr,lag = cross_corr(line_data_aligned_norm,line_ref) # determine shift by xcorrelation
        else:
            xprime = x
            lag = 0
        xprime = xprime + lag # add shift to xprime coordinates
        line_data_aligned = np.interp((xprime), x, line_data) # interpolate data from x to xprime
        fluo_aligned[:,i] = line_data_aligned
        X_shift[:,i] = (xprime-x)*scan_pixel_size[0]
        i = i+1
    
    Y_shift = Y_shift*0 + Ylag*scan_pixel_size[1]
    return fluo_aligned, X_shift, Y_shift

# End of correction methods #

def read_data_meta(path):
    h5file = h5py.File(path,'r')
    command = str(h5file['entry']['description'][()])[3:-2] # Reading only useful symbols
    motor_positions = {
            # Detector positions
            "delta":    h5file['entry']['snapshot']['delta'][()],
            "gamma":    h5file['entry']['snapshot']['gamma'][()],           
            "gonphi":   h5file['entry']['snapshot']['gonphi'][()],
            "gontheta": h5file['entry']['snapshot']['gontheta'][()],
            "radius":   h5file['entry']['snapshot']['radius'][()],
            "energy":   h5file['entry']['snapshot']['energy'][()]
            }
    
    scan_position_x = h5file['entry']['measurement']['pseudo']['x'][()]
#    scan_position_y = h5file['entry']['measurement']['pseudo']['y'][()]
    scan_position_z = h5file['entry']['measurement']['pseudo']['z'][()]
    
    return command, motor_positions, scan_position_x, scan_position_z

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
    nPoints = len(h5file[entry_keys[0]]['measurement']['xspress3']['data'][()])
#    nLines  = np.shape(entry_keys)
    data = np.zeros((nPoints, len(entry_keys)))
    for key in entry_keys:
        data_temp = h5file[key]['measurement']['xspress3']['data'][:,3,roi[0]:roi[1]]
        data[:,ii] = np.sum(data_temp,1)       
        ii = ii+1
    return data

# INPUTS ######################################################################
year        = "2020"                                                           #The year for the experiemnt
beamtimeID  ="2020101408"                                                      #The beamtimeID
sample_name = r"0002_sample_P246_AB"                                           #The name for the p10 newfile
#scan_number = np.linspace(109,128,20)                                          #The scan numbers, can be the list
#scan_number = np.linspace(137,156,20)                                          #The scan numbers, can be the list
scan_number = [146]
rocking_motor = "gontheta" # "gonphi"

xrd_analysis     = True
xrf_analysis     = True
q_space_analysis = False
save_data_xrd    = False
save_data_xrf    = False
align_scans      = True
align_data_xrd   = True
save_data_xrd_interpolate = False
save_data_xrf_interpolate = False

xrf_roi     = [870,970]
#xrd_roi     = [100,200,275,375]
xrd_roi     = [0,500,0,500]
scan_pixel_size = [0.08,0.08] # um from the motor positions in the data
###############################################################################

#determing the paths"
path_root       = r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/raw/%s"%(beamtimeID,sample_name)
processing_path = r"/home/dzhigd/work/projects/Qdevs_2020_NanoMAX/data/scan_%d_%d"%(scan_number[0],scan_number[-1])

try:
    os.mkdir(processing_path)
    print("--Processing folder created--")
except:
    print("--Processing folder already exists--")
    
# DATA IMPORT #################################################################
orig_pixel_size = [0.01,0.01] # um pixel size of synthetic image
rocking_angle = np.zeros((len(scan_number),1))

##### Start here
offset_top_left = [-2.7,0]

for ii in range(0,len(scan_number)):
    print(('--Processing scan %d--')%(int(scan_number[ii])))
    #determing the paths"
    
    data_meta_path = os.path.join(path_root, "%06d.h5"%(int(scan_number[ii])))
    data_xrd_path  = os.path.join(path_root, "scan_%06d_merlin.hdf5")%(int(scan_number[ii]))
    data_xrf_path  = os.path.join(path_root, "scan_%06d_xspress3.hdf5")%(int(scan_number[ii]))
    
    # Load the data ###########################################################
    [command, motor_positions, scan_positions_x, scan_positions_z] = read_data_meta(data_meta_path) 
    scan_positions_xz = np.array([scan_positions_x,scan_positions_z])
    
    if xrd_analysis == True:
        data_xrd  = read_data_merlin(data_xrd_path,xrd_roi)
    if xrf_analysis == True:
        data_xrf  = read_data_xspress3(data_xrf_path,xrf_roi)

    rocking_angle[ii] = motor_positions[rocking_motor]
    ###########################################################################
    
    # Load the data ###########################################################
#    figure(num=13)
#    imshow(np.log10(np.sum(data_xrd,0)))
    
    # Save the original xrd data ############################################## 
    if save_data_xrd == True:
        save_path = os.path.join(path_root,"converted", "scan_%06d_merlin_%d_%d_%d.npz")%(int(scan_number[ii]),(data_xrd.shape[0]),(data_xrd.shape[1]),(data_xrd.shape[2]))      #The path to save the results
    #    save_bin_data_nanomax(save_path, data_xrd)
        np.savez_compressed(save_path,data_xrd)       
        save_path = os.path.join(path_root,"converted", "scan_%06d_scan_positionsxz.npz")%(int(scan_number[ii]))      #The path to save the results
        np.savez_compressed(save_path,scan_positions_xz)
        print(('Saved scan %d to %s')%(int(scan_number[ii]),save_path))
    
    if save_data_xrf == True:        
        data_xrf = np.float32(data_xrf/np.max(data_xrf))   
        imsave(os.path.join(processing_path, ("scan_%06d_xspress3.tif")%(int(scan_number[ii]))), data_xrf)
        print(('Saved scan %d to %s')%(int(scan_number[ii]),processing_path))
    
    print("--Import Done!--")
    
    if align_scans == True:
        # make reference image   
        ref_path = r"/home/dzhigd/work/projects/Qdevs_2020_NanoMAX/data/reference_%d_%d.tif"%(scan_number[0],scan_number[-1])
        
        with tf.TiffFile(ref_path) as tif:
            ref_image = tif.asarray()
        # Do for each angle:
        if ii == 0:
            reference_image,mask = generate_ref_image(data_xrf,ref_image,scan_pixel_size,orig_pixel_size,offset_top_left)
        
        # find misalignments
        fluo_aligned, X_shift, Y_shift = align_image(data_xrf, reference_image, scan_positions_xz, scan_pixel_size)
        
        if align_data_xrd == True:
            # interpolate the misalignments from regular grid onto the motor positions so they can be just added
            X_shift = grid_to_image(X_shift, scan_positions_xz)
            #Y_shift = grid_to_image(Y_shift, positions)
            
            # add misalignments to positions
            positions_new = scan_positions_xz + np.array([Y_shift.transpose().ravel(),X_shift.transpose().ravel()])
            
            # find interpolant
            F = interpolate.LinearNDInterpolator(scan_positions_xz.transpose(), data_xrd, fill_value=np.nan, rescale=True)
            
            # interpolate, this is the output of aligned diffraction data
            xrd_interp = F(positions_new.transpose())
            
            # sum data to plot total xrd map
            xrd_int = np.sum(np.sum(xrd_interp,axis=2),axis=1)
            xrd_int.shape
            xrd_image = np.reshape(xrd_int,data_xrf.transpose().shape).transpose()
            
        #    plt.figure(num=1)
        #    plt.imshow(reference_image)
        
#            plt.figure(num=1,figsize=(15,15))
#            plt.imshow(fluo_aligned)
#            plt.grid()
#            
#            plt.figure(num=3,figsize=(15,15))
#            plt.imshow(xrd_image)
#            plt.grid()
        # Q-space #################################################################
        
        # Save the interpolated xrd data ########################################## 
        if save_data_xrd_interpolate== True and align_data_xrd == True:
            # Save NPZ compressed numpy array
#            save_path = os.path.join(processing_path, "scan_%06d_merlin.npz")%(int(scan_number[ii]))      #The path to save the results
#            np.savez_compressed(save_path,xrd_interp)
            
            # Save binary uncompressed array
    #        save_bin_data_nanomax(processing_path, xrd_interp)
            
            # Save mat compressed array for matlab                 
            save_path = os.path.join(processing_path, ("scan_%06d_merlin.mat")%(int(scan_number[ii])))    #The path to save the results
            data_dic = {"data":xrd_interp[4040:-1,:,:], "command": command, "motor_positions":motor_positions, "scan_positions_x":scan_positions_x, "scan_positions_z":scan_positions_z }
            scio.savemat(save_path, data_dic) 
        if save_data_xrf_interpolate == True:            
            fluo_aligned = np.float32(fluo_aligned/np.max(fluo_aligned))   
            imsave(os.path.join(processing_path, ("scan_%06d_xspress3_aligned.tif")%(int(scan_number[ii]))), fluo_aligned)
            print(('Saved scan %d to %s')%(int(scan_number[ii]),processing_path))
        
    print(('--Processing scan %d done!--')%(int(scan_number[ii])))
    
print("--Overall Processing Done!--")