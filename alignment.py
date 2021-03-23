import phidl.geometry as pg
from phidl.device_layout import _parse_layer
from matplotlib import pyplot as plt
from PIL import Image, ImageFilter
from matplotlib.pyplot import imshow
import aggdraw
import numpy as np
import yaml
from scipy import ndimage
import tifffile as tf
import json
import cv2
from scipy import ndimage

# fluorescence image

datapath = r"C:\Users\tostanke\Dropbox\Nanomax2020\scan109_128\scan_000113_xspress3.tif"

# reference image
with tf.TiffFile(r"C:\DDzhigaev\Projects\ACTIVE\Qdev_Microsoft\Experiments\NanoMAX_14102020\vertical_array10.tif") as tif:
    image = tif.asarray()
    
fluo_image = tf.imread(datapath)
plt.figure(figsize=(5,5))
imshow(fluo_image)

orig_pixel_size = [0.01,0.01] # um pixel size of synthetic image
image_size = fluo_image.shape # measured image size in pixels
pixelsize = [0.083, 0.15] # x/y pixel size of scan um
tl_offset = [-2.5,0] # top left offset in um. adjust so that reference looks like scan

padsize = 1000

##############################################

#scan range (um)
sr_tl = ((np.array(tl_offset))/orig_pixel_size).astype(np.int32)+np.array([padsize,padsize])
sr_br = ((np.array(tl_offset) + np.array(image_size)*np.array(pixelsize))/np.array(orig_pixel_size)).astype(np.int32)+np.array([padsize,padsize])

image = cv2.rotate(image,cv2.ROTATE_90_COUNTERCLOCKWISE)
image = cv2.copyMakeBorder( image, 1000, 1000, 1000, 1000, cv2.BORDER_CONSTANT,0)

image_cropped = image[sr_tl[0]:sr_br[0],sr_tl[1]:sr_br[1]]
    
ret,binary = cv2.threshold(image_cropped,127,255,cv2.THRESH_BINARY)

kernel = np.ones((5,5),np.uint8)
binary = cv2.dilate(binary,kernel,iterations = 5)

smoothed = ndimage.gaussian_filter(binary, 20)
resized = cv2.resize(smoothed,(image_size[1],image_size[0]))

plt.figure(figsize=(5,5))
imshow(resized)

reference = resized[:,:,0]

data = fluo_image

from scipy.optimize import minimize
# normalization before ACF
def normalize(data):
    return (data - np.mean(data, axis=0)) / (np.std(data, axis=0))

def transform(y, params):
    x = np.linspace(0,len(y),len(y))
    # generate new x values - xprime by adding a polynomial
    xprime = x +  params[0]*x + params[1]*(x-np.mean(x))**2 #+ params[2]*x**3
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
    # caucluates error function for minimization. minimization parameter is negative max of cross-correlation
    y,_ = transform(y1, params)
    maxcorr, lag = cross_corr(y, y2)
    return -maxcorr

def find_distortion(y1, y2):
    #finds optimal distortion parameters based on minimization of cross-correlation function
    x0 = [0,0]
    res = minimize(errorf, x0, method='Nelder-Mead', tol=1e-6, args=(y1,y2))
    return res.x

# run through all the lines
aligned_data = data_upscale.copy() #placeholder array
plt.figure()
i=0
for line_data in data.transpose():
    line_data_norm = normalize(line_data) # normalize data
    line_ref = normalize(reference_upscale[:,i]) #normalize reference
    params = find_distortion(line_data_norm,line_ref) # find distortion parameters
    line_data_aligned, xprime = transform(line_data, params) # transform data based on found parameters. does not include shift
    line_data_aligned_norm = normalize(line_data_aligned) # normalize transformed data
    corr,lag = cross_corr(line_data_aligned_norm,line_ref) # determine shift by xcorrelation
    x = np.linspace(0,len(line_data),len(line_data)) # original x coordinates
    xprime = xprime + lag # add shift to xprime coordinates
    line_data_aligned = np.interp((xprime), x, line_data) # interpolate data from x to xprime
    plt.plot(line_data_aligned)
    aligned_data[:,i] = line_data_aligned
    i = i+1

plt.figure()
plt.figure(figsize=(20,20))
imshow(data_upscale)
plt.figure(figsize=(20,20))
imshow(aligned_data)