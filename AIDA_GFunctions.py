################################################################################
#
#   File:       AIDA_GFunctions.py
#
#   Summary:    Set of "G"eneric functions used in the Adaptive Image 
#               Deconvolution Algorithm (AIDA) package (AIDAP) and useful in
#               interactive sessions
#
#   Authors:    Erik F.Y. Hom (Sedat Lab, UCSF) with assistance from
#               Sebastian Haase (Sedat Lab, UCSF) and
#               Timothy K. Lee (Undergraduate Intern (UCB), Sedat Lab)
##
#   Other:      See 'AIDA_version.py' for date and version details
#               See 'LICENSE.txt' for use, license, and reference details
#
################################################################################

# from Priithon.all import U
# Used from priithon:
#  - N (numpy)
#  - U (useful) [ nb , findMax, max2d, fitGaussian, fitPoly, histogram, 
#                 generalhistogram, saveImg, saveImg_seq, findMin ]
#  - F (fftfuncs) : noiseArr
#  - Mrc
import numpy as N
from Priithon_heritage import useful as U
from Priithon_heritage import fftfuncs as F
from Priithon_heritage import Mrc

#from Priithon import fftw
from Priithon_heritage import fftw

import os, string, time, types
# import pyfits  Clement: removed
import astropy.io.fits as iofits

from PIL import Image, ImageOps     # import PIL components to also output JPEG files of results
import scipy.misc
from astropy.io.fits.hdu.image import PrimaryHDU


#######  CLASS: N_random_array  #######
####    [checked]
class N_random_array:
    def normal(mean, std, shape):
        return N.random.normal(mean,std,shape)
    def poisson(mean, shape=1):
        return N.random.poisson(mean,shape)
####


#######  FUNCTION: PrintTime  #######
####    [checked]
def PrintTime(timeinseconds, denominator=1.): # used in 'AIDA_Functions.py'
    """
    Prints 'timeinseconds' in secs, mins, or hrs, whichever is the most
    appropriate denomination

    if denominator != 1: divide timeinseconds by it
    check first if denominator is zero
    """
    if denominator == 0:
        print "<cannot-divide-by-zero>"
        return
    timeinseconds = timeinseconds / float(denominator)
    

    if timeinseconds <= 60:

        print '%.4g' %timeinseconds, 'secs',
    elif 60 < timeinseconds <= 3600:

        print '%.4g' %(timeinseconds/60.), 'mins',
    else:

        print '%.4g' %(timeinseconds/3600.), 'hrs', 
####


#######  FUNCTION: CalculateImageData  #######
####    [checked]
def CalculateImageData(image, background=0., sigma_det=None, wiener=None,
        dtype=N.float64): # used in 'AIDA_Functions.py'
    """
    Returns 'background' subtracted 'image' along with 'inv_w' 
    (inverse of weights), 'sigma_det' (Gaussian std of detector noise),
    and an estimate for the 'wiener' coefficient if requested (returns 'None')
    for this variable otherwise. 'dtype' sets the return dtype for 'inv_w'.

    Returns:
    (image, inv_w, sigma_det, wiener)
    """

    ## substract off baseline (arbitrary correction at the moment)
    image -= background
    (inv_w, sigma_det) = ImageNoiseWeights(image, baseline=0., 
            sigma_det=sigma_det, dtype=dtype) #@

    if wiener == 'wiener':

        noise = N_random_array.normal(mean=0., std=sigma_det, 
                shape=image.shape)
        ## use object ~ image
        wiener = WienerWeight(noise=noise, object=image) #@
    else:

        wiener = None

    return (image, inv_w, sigma_det, wiener)
####


#######  FUNCTION: ImageNoiseWeights  #######
####    [checked]
def ImageNoiseWeights(image, baseline=0., sigma_det=None, dtype=N.float64):
        
    """
    Computes the weights, w(r), for use in calculating the data fidelity
    term, Jn.  Note that w(r) should be squared for the Jn calculation

    'baseline' = 0. is the default, which assumes there are intensity values
    less than this from which to calculate the detection noise statistics

    Returns a tuple of (1/weights, noise, and sigma_det)
    where 'noise' is to be used in calculating the Weiner parameter and
    'sigma_det' is the standard deviation for the Gaussian readout noise

    Returns:
    ((1./w).astype(dtype), N.sqrt(sigma_det2))
    """

    if sigma_det:
    
        sigma_det2 = sigma_det*sigma_det
    else:

        neg_pixels = N.extract(N.less_equal(image, baseline), image) - \
                baseline

        if len(neg_pixels) < 2:

            message = "\n'image' does not have sufficient number of " + \
                      "negative pixels!\nCheck that 'baseline' is not " + \
                      "lower than the minimum image pixel value"
            raise ValueError, message
            
        mean_neg_pixels = neg_pixels.mean()
        ## sigma_det2 is the (constant) gaussian detection error/pixel
        ## from formula sigma_det**2 = 0.5*N.pi*(noise.mean())**2
        sigma_det2 = 0.5*N.pi*(mean_neg_pixels*mean_neg_pixels)

### EHom (20130625): this is a mistake in the JOSA A paper and is not correct.
### sigma_det2 simply must be positive given the quantal nature of detection, and so
### cannot < 0 and Eq. (17) in Hom et al. is still well defined.  However, it is not well defined
### for Eq. (32).  lambda_object will be negative or ill defined for sqrt(2*pi*sigma_det) <= 1
### Thus, sigma_det must be > 1/(2*pi) and sigma_det2 must be > (1/2*pi)**2

### Old code:
    if sigma_det2 < (N.pi / 2.):
   
        message = "\n\nWarning: sigma_det is " + str(N.sqrt(sigma_det2)) + ", which is below the " + \
                  "theoretical limit of sqrt(pi/2)!\nResetting to theoretical limit...\n\n"
        print message
       
        sigma_det2 = N.pi / 2.
### New code:
#     sigma_det2_mintheory = (2.*N.pi)**(-2)   # this value is 0.0253302959106,
#                                              # i.e., sigma_det_mintheory = 0.159154943092
#                                              # Equivalence of Eq. (30) in Hom et al. works for 
#                                              # lambda_object ~< 10, which means
#                                              # (sqrt(2*pi*sigma_det) - 1)**-1 ~< 10
#                                              # or sigma_det ~> 1.01/(2*pi)=0.160746492523
#                                              # or sigma_det2 ~> 0.0258394348584
#                                              # this is greater (thankfully, than sigma_det2_mintheory
#                                              # by ~0.0005; use 0.001 as epsilon below to set
#                                              # minimal practical sigma_det2 value
#                                              
#     if sigma_det2 <= sigma_det2_mintheory:
#     
#         message = "\n\nWarning: sigma_det is " + str(N.sqrt(sigma_det2)) + ", which is below the " + \
#                   "theoretical limit of (1/2*pi)+epsilon!\nResetting to theoretical limit...\n\n"
#         print message
#         
#         sigma_det2 = sigma_det2_mintheory + 0.001   # value should be greater than minimum or lambda_object will explode
        

    sigma_photon2 = N.maximum(image-baseline, 0.)
    w = sigma_det2 + sigma_photon2

    if w.min() == 0:

        w_min = (N.extract(sigma_photon2 > 0, sigma_photon2)).min()
        w = N.where(w==0, w_min, w)

    return ((1./w).astype(dtype), N.sqrt(sigma_det2))
####


#######  FUNCTION: UpdateWeights  #######
####    [checked]
def UpdateWeights(noiseless_image, sigma_det):  # not used...
        
    """
    Updates weights, w(r), for use in calculating the data fidelity
    term, Jn. 
    
    Returns: (1./w)
    """

    sigma_photon2 = N.maximum(noiseless_image, 0.)
    w = sigma_det*sigma_det + sigma_photon2

    if w.min() == 0:

        w_min = (N.extract(sigma_photon2 > 0, sigma_photon2)).min()
        w = N.where(w==0, w_min, w)

    return (1./w)
####


#######  FUNCTION: WienerWeight  #######
####    [checked]
def WienerWeight(dark_image, object):   
    """
    Returns the Wiener filter weight parameter (as N.float64) for use
    in the function 'WienerFilter':
    
    weight = <||NOISE||^2> / <||OBJECT||^2>, which is
        ~the inverse of the signal-to-noise ratio   
    
    'dark_image' is the detection noise
    'OBJECT' is the realFT of the object (~image)

     If 'dark_image' is a number, this function will assume this is the
        standard deviation of the detector noise

    Returns:
    (absNOISE*absNOISE).mean()/(absOBJECT*absOBJECT).mean() 
    """

    shape = object.shape[:-1] + (object.shape[-1]/2 + 1,)

    if type(dark_image) in (float,int):
    
       dark_image = F.noiseArr(shape=shape, stddev=dark_image, mean=0.0, dtype=N.float32)
    
    elif isinstance(dark_image,N.ndarray):
    
        dark_image = (dark_image - dark_image.mean()).astype(N.float32)
    else:
    
        message = "\n'dark_image' must be a number or an array!"
        raise ValueError, message

    OBJECT = N.empty(shape=shape, dtype=N.complex64)
    NOISE = N.empty(shape=shape, dtype=N.complex64)
    fftw.rfft(a=object.astype(N.float32), af=OBJECT, inplace=False)
    fftw.rfft(a=dark_image, af=NOISE, inplace=False)
    absNOISE = N.abs(NOISE)
    absOBJECT = N.abs(OBJECT)


    ## returns the ratio of the noise:image power spectral densities
    return (absNOISE*absNOISE).mean()/(absOBJECT*absOBJECT).mean()
####


#######  FUNCTION: ProcessPSF  #######
####    
def ProcessPSF(PSF, cropshape=None, center=None, exclude_radius=5, 
        clean=(1,1,1,1), background_percent=0., nsigmas=2.,
        threshold_percent=1e-3, fill=0., dtype=N.float64):
    """
    Returns a cleaned, centered, resized PSF for use in AIDA given a raw PSF
    """

    if cropshape is None:
    
        cropshape = PSF.shape
    elif len(cropshape) < PSF.ndim:
    
        message = "\n'dimension' cannot be larger than length of 'cropshape'!"
        raise RuntimeError, message

    ## locate the center of the PSF image
    if center is None:

        (new_center, pixel_center) = LocatePSFcenter(PSF=PSF, 
                xyCentroidSize=None) #@
    elif center in ('origin', 'array_center'):
    
        new_center = center
    else:
    
        new_center = tuple(center)
        center = None       # reset 'center' flag to None for below shifting

    ## resize PSF if _bigger_ than 'cropshape'
    (tempPSF, new_center, pixel_center) = ConditionalResizePSF(PSF=PSF, 
            cropshape=cropshape, center=new_center, condition='1', fill=fill) #@

    ## clean-up the PSF image
    if N.sum(clean):
    
         tempPSF[:] = CleanPSF(PSF=tempPSF, clean=clean, center=new_center, 
                exclude_radius=exclude_radius, delta=0.3, 
                background_percent=background_percent, nsigmas=nsigmas, 
                threshold_percent=1e-3, threshold_value=None, fill=fill,
                dtype=dtype) #@
    else:

        tempPSF[:] = tempPSF.astype(dtype)

    if center is None:      # either inputted center or deduced center

        ## shift PSF center to origin with subresolution centering
        ## with order=5 spline interpolation

    
    #Clement: shift workaround
        tempPSF = U.shift2D_workaround(tempPSF.copy(),-N.array(new_center))
         
#         U.nd.shift(tempPSF.copy(), shift=-N.array(new_center), output=tempPSF, 
#                 order=5, mode="wrap")
#         

        ##  N.B.  Multiplying cOTF by PhaseShift function leads to "cross-hair"
        ##  FT boundary artifacts; therefore, resorted to U.nd.shift operation 
        ##  with spline interpolation instead

        ## re-threshold PSF after centering to origin to prevent non-negative
        ## values
        if N.sum(clean[-2:]):
    
             tempPSF[:] = CleanPSF(PSF=tempPSF, clean=(0,0)+clean[-2:], 
                    center=None, exclude_radius=exclude_radius, delta=0.3, 
                    background_percent=background_percent, nsigmas=nsigmas, 
                    threshold_percent=1e-3, threshold_value=None, fill=fill,
                    dtype=dtype) #@
    elif center == 'array_center':
    
        tempPSF = ArrayCenter2Origin(tempPSF)


    ## resize PSF after cleaning if _smaller_ than 'shape'
    (tempPSF, new_center, pixel_center) = ConditionalResizePSF(PSF=tempPSF, 
            cropshape=cropshape, center=new_center, condition='2', fill=0) #@

    return tempPSF
####


#######  FUNCTION: LocatePSFcenter  #######
####    [checked]
def LocatePSFcenter(PSF, xyCentroidSize=None):
    """
    Locates the true center and approximate pixel center of a PSF using
    centroid fitting
    
    return (center, pixel_center)
    """

    center = FindPSFcentroid(PSF, xyCentroidSize=None) #@
    pixel_center = tuple(N.around(center).astype(N.int32))

    return (center, pixel_center)
####


#######  FUNCTION: FindPSFcentroid  #######
####    [[checked]]
def FindPSFcentroid(PSF, xyCentroidSize=None):
    """
    Returns the center coordinates of a 2D _or_ 3D 'PSF' array, assuming a
    somewhat Gaussian PSF whose center can be determined by centroid fitting.
    Assumes that there are no bad pixels - so it is recommendend that the PSF
    supplied be cleaned or median-filtered first.
    
    'xyCentroidSize' is the side-length of the square area around the xy-peak 
    'xyCentroidSize' should be an odd number
    'xyCentroidSize' should be (at least) 7 times the sigma of a "gaussian 
                     shaped" peak
                     
    Returns: (zc, yc, xc)
    """
    
    if xyCentroidSize is None:
    
        xyCentroidSize = PSF.shape[-1] / 4
    
    if IsEven(xyCentroidSize): #@
    
        xyCentroidSize -= 1
    
    if PSF.ndim == 2:
    
        (junk, junk, y0, x0) = U.findMax(PSF)

        percentage = xyCentroidSize**2. / N.product(PSF.shape[-2:]) * 100
        #U.DEBUG_HERE()
        #print "#DEBUG:", percentage
        #PSF.tofile("seb_PSF")
        #import cPickle
        #cPickle.dump(PSF, file('seb_PSF', 'w'), 2)
        #levelAtRim = U.topPercentile(N.around(PSF).astype(N.int32), 
        #        percentage)
        #print "#DEBUG:      ", levelAtRim
        # workaround for broken Priithon's U.topPercentile
        levelAtRim = U__topPercentile(PSF, percentage)

        yxC = tuple(U.findMax(PSF)[-2:])
        
        if yxC == (0,0):
        
            return N.asarray((0,0))
        else:
        
            return CentroidOverSquare(array=PSF, center=yxC, 
                    size=xyCentroidSize, background_threshold=levelAtRim) #@
                

    else:       ## then 3D psf
    
        zprofile = U.max2d(PSF)     ## max val for each x-y plane

        (PSFmax, junk, junk, z0) = U.findMax(zprofile)      
        zprofile_baseline = U__topPercentile(zprofile, 99)             ## value of 99%
        FWHM = len(N.where(zprofile > (PSFmax/2))[0])/2    ## FWHM/2

        #fit = U.fitGaussian1D((z0, FWHM, PSFmax - zprofile_baseline), 
        #        zprofile - zprofile_baseline)
        fit = U.fitGaussian(zprofile - zprofile_baseline, (z0, FWHM, PSFmax - zprofile_baseline))
        zc = fit[0][0]      ## PSF center final z

        ## sandwich above & below: find (2D) yx center in both planes
        za = int(zc)    ## bottom
        zb = za + 1     ## top
        weight_b = zc - za
        weight_a = 1. - weight_b

        (junk, junk, y0, x0) = U.findMax(PSF[za])

        percentage = xyCentroidSize**2. / N.product(PSF.shape[-2:]) * 100
        levelAtRim = U__topPercentile(PSF[za:za+2],
                                      percentage)

        yxC = U.findMax(PSF[za])[-2:]
        
        if yxC == (0,0):
        
            return (0, yc, xc)      ## assumes zmax is at zero
        else:

            (ya, xa) = CentroidOverSquare(array=PSF[za], center=yxC, 
                    size=xyCentroidSize, background_threshold=levelAtRim) #@
            yxC = U.findMax(PSF[zb])[-2:]
            (yb, xb) = CentroidOverSquare(array=PSF[zb], center=yxC, 
                    size=xyCentroidSize, background_threshold=levelAtRim) #@

            xc = weight_a*xa + weight_b*xb
            yc = weight_a*ya + weight_b*yb

            return (zc, yc, xc)
####


#######  FUNCTION: IsEven  #######
####
def IsEven(number):

    return (number % 2 == 0)
####


#######  FUNCTION: IsOdd  #######
####
def IsOdd(number):

    return (number % 2 != 0)
####


#######  FUNCTION: CentroidOverSquare  #######
####
def CentroidOverSquare(array, center=None, size=None, background_threshold=0):
    """
    Returns a tuple of coordinates for the centroid of an 'array' 

    Returns:
    tuple(center + relative_center - (size - 1)/2)
    """
    
    if center is None:
    
        center = N.array(array.shape)/2
    else:

        center = N.asarray(center)

    if size is None:
    
        size = center.copy()
    else:
    
        size = N.asarray(size)
                
    if center.ndim == 0:

        center.shape = 1
    if len(center) < array.ndim:

        #seb center.resize(array.ndim)
        center = N.resize(center, array.ndim)
    
    if size.ndim == 0:

        size.shape = 1
        
    if len(size) < array.ndim:
        #seb size.resize(array.ndim)
        size = N.resize(size, array.ndim)
    for s in size:
        
        if IsEven(s): #@
    
            s -= 1

    # CHECK should we raise if lower or upper outside array !?
    lower = map(int, center-size/2.+1)
    upper = map(int, center+size/2.+1)
    
    for i in range(array.ndim):

        array = array[lower[i]:upper[i]]
        array = N.transpose(array, [array.ndim-1] + range(array.ndim-1))

    #array = U.thrsh(array, min=background_threshold, force0Base=1)
    array = N.where(  array > background_threshold, array-background_threshold, 0   )
    relative_center = U.nd.center_of_mass(array)

    return tuple(center + relative_center - (size - 1)/2)
####    


#######  FUNCTION: ConditionalResizePSF  #######
####    [checked]
def ConditionalResizePSF(PSF, cropshape, center, condition, fill=0.):
    """
    Returns tuple of resized PSF along with 'center' and 'pixel_center'
    depending on 'cropshape' and whether 'condition' case is satisfied
    by the input variables:
    
    'condition' = 1:        PSF has larger shape than cropshape
    'condition' = 2:        PSF has smaller shape than cropshape

    Returns:
    (PSF, center, pixel_center)
    """
    
    if type(center) not in (types.IntType, types.FloatType) and \
           not isinstance(center,  N.number) : #seb
    
            pixel_center = center
    else:
    
        pixel_center = tuple(N.around(center).astype(N.int32))

    ## first condition, PSF is bigger than cropshape
    ## second condition, PSF is smaller than cropshape
    if (condition == '1' and N.product(PSF.shape) > N.product(cropshape)) or \
            (condition == '2' and N.product(PSF.shape) < N.product(cropshape)):

        tempPSF = ResizePSF(PSF=PSF, new_shape=cropshape, 
                    pixel_center=pixel_center, fill=fill) #@

        if type(center) is not types.StringType:    # i.e. 'origin' or
                                                    # 'array_center'
            ## Re-locate the center of the PSF image
            (center, pixel_center) = LocatePSFcenter(PSF=tempPSF, 
                    xyCentroidSize=None) #@

        return (tempPSF, center, pixel_center)

    else:
    
        return(PSF, center, pixel_center)
####


#######  FUNCTION: ResizePSF  #######
####    [checked]
def ResizePSF(PSF, new_shape, pixel_center=None, fill=0.):   
    """
    Takes PSF data in an array and resizes it to 'new_shape'.  Then:
    
    (1) PSF is centered about the maximum value (we assume all cosmic ray pixels
        have been removed!)
    (2) PSF image/array is truncated if it is larger than 'new_shape'
    (3) Values beyond the PSF array for PSFs with a shape smaller than 
        'new_shape' are filled in (so that PSF.shape = 'new_shape') using a 
        mean value calculated from edge values of the PSF image
    (4) A new PSF array is returned that matches 'new_shape'
    """

    PSF.shape = tuple(N.compress(N.array(PSF.shape) > 1, PSF.shape))
    oldshape = PSF.shape

    if len(oldshape) != len(new_shape):
    
        message = "\n'new_shape' dimension does now match PSF input!"
        raise ValueError, message   

    array_center = tuple(N.floor(N.array(oldshape)/2.))
    tempPSF = N.empty(shape=oldshape, dtype=PSF.dtype)
    loweroffset = N.abs(N.floor((N.array(oldshape) - \
            N.array(new_shape))/2.).astype(N.int))
    upperoffset = N.abs(N.floor((N.array(new_shape) - \
            N.array(PSF.shape))/2.).astype(N.int))

    if pixel_center is None:        # estimate center (in pixels) based on 
                                    # centroid approach 
        if len(new_shape) == 2:

            pixel_center = (ycenter, xcenter) = tuple(N.around(
                    FindPSFcentroid(PSF, xyCentroidSize=None)).astype(N.int32))
        else:
        
            pixel_center = (zcenter, ycenter, xcenter) = tuple(N.around(
                    FindPSFcentroid(PSF, xyCentroidSize=None)).astype(N.int32))
    
    elif type(pixel_center) is not types.StringType:    # i.e. 'origin' or
                                                        # 'array_center'
        if len(pixel_center) != len(new_shape):
    
            message = "\n'pixel_center' coordinates do not match " + \
                    "'new_shape' dimensions!"
            raise ValueError, message

    if pixel_center in ('array_center', array_center): 
    
        tempPSF = PSF
        pixel_center = 'array_center'
    elif pixel_center in ('origin', (0,)*PSF.ndim):
    
        tempPSF = Origin2ArrayCenter(PSF)
        pixel_center = 'origin' 
    else:
        # shift to array_center based on estimated pixel_center position
        shift = (N.array(array_center)-N.array(pixel_center))
        #Clement: shift workaround
        tempPSF = U.shift2D_workaround(PSF,shift)
#         U.nd.shift(PSF, shift=shift, output=tempPSF, order=3, mode="wrap")
    
    ## crop if oldshape is > newshape
    if oldshape[-1] > new_shape[-1]:    ## x-dimension

        tempPSF = tempPSF[...,loweroffset[-1]:(oldshape[-1]-upperoffset[-1])]

    if oldshape[-2] > new_shape[-2]:    ## y-dimension

        tempPSF = tempPSF[...,loweroffset[-2]:(oldshape[-2]-upperoffset[-2]),:]

    if len(new_shape)==3 and (PSF.shape[-3] > new_shape[-3]):   ## z-dimension

        tempPSF = tempPSF[loweroffset[-3]:(oldshape[-3]-upperoffset[-3]), ...]

    ## expand if oldshape < newshape
    if tempPSF.shape != new_shape:

        if len(new_shape) == 2:

            edge_length = N.floor(N.sqrt(tempPSF.shape)).astype(N.int)
            y_edge_range = [edge_length[-2], tempPSF.shape[-2]-edge_length[-2]]
            x_edge_range = [edge_length[-1], tempPSF.shape[-1]-edge_length[-1]]
            loweroffset = N.floor((N.array(new_shape) - \
                    N.array(tempPSF.shape))/2.).astype(N.int)
            upperoffset = N.abs(N.floor((N.array(tempPSF.shape) - \
                    N.array(new_shape))/2.).astype(N.int))
            PSFout = N.zeros(new_shape, dtype=PSF.dtype) + fill
            PSFout[loweroffset[-2]:(new_shape[-2]-upperoffset[-2]),
                    loweroffset[-1]:(new_shape[-1]-upperoffset[-1])] = tempPSF
        elif len(new_shape) ==3:

            edge_length = N.floor(N.sqrt(tempPSF.shape)).astype(N.int)
            z_edge_range = [edge_length[-3], tempPSF.shape[-3]-edge_length[-3]]
            y_edge_range = [edge_length[-2], tempPSF.shape[-2]-edge_length[-2]]
            x_edge_range = [edge_length[-1], tempPSF.shape[-1]-edge_length[-1]]
            loweroffset = N.floor((N.array(new_shape) - \
                    N.array(tempPSF.shape))/2.).astype(N.int)
            upperoffset = N.abs(N.floor((N.array(tempPSF.shape) - \
                    N.array(new_shape))/2.).astype(N.int))    
            PSFout = N.zeros(new_shape, dtype=PSF.dtype) + fill
            PSFout[loweroffset[-3]:(new_shape[-3]-upperoffset[-3]),
                    loweroffset[-2]:(new_shape[-2]-upperoffset[-2]), \
                    loweroffset[-1]:(new_shape[-1]-upperoffset[-1])] = tempPSF

        tempPSF = PSFout

    ## now shift back according to pixel_center input value:
    ## if 'array_center', leave alone
    ## if 'origin', shift back to origin
    ## if None or something else explicit, shift to array_center
    if pixel_center == 'origin':
    
        return ArrayCenter2Origin(tempPSF)
    else:
    
        return tempPSF
####


#######  FUNCTION: CleanPSF  #######
####    [checked]
def CleanPSF(PSF, clean=(1,1,1,1), center=None, exclude_radius=5, delta=0.3,
        background_percent=0., nsigmas=0., threshold_percent=1e-3,
        threshold_value=None, fill=0., dtype=N.float64):   
    """
    Returns a cleaned-up version of a normalized PSF image for use in AIDA.
    Note that this function assumes that these images have been pre-processed
    already by dark current subtraction and flat-fielding.
        
    'clean' =   (1,0,0,0)       - remove bad pixels
                (0,1,0,0)       - subtract background 
                                    ('clean_threshold_percent'/'clean_sigma')
                (0,0,1,0)       - threshold to zero
                (0,0,0,1)       - normalize by sum

    Returns: PSFout.astype(dtype)
    """

    dimension = PSF.ndim
    shape = N.array(PSF.shape)
        
    if len(PSF.shape) == 3 and PSF.shape[-3] == 1:
    
        PSF.shape = PSF.shape[-2:]
        shape = shape[-2:]

    elif len(PSF.shape) not in (2, 3):
    
        message = "\n'PSF' shape must be 2D or 3D!"
        raise RuntimeError, message 

    if len(clean) != 4:
    
        message = "\n'clean_level' must be be a list or tuple of 0's and " + \
                "1's of length 4 (e.g., (1,1,1,1))!"
        raise ValueError, message
    else:
    
        for i in range(4):
        
            if clean[i] not in (0,1):
            
                message = "\n'clean_level' must be be a list or tuple of " + \
                        "0's and 1's\nentry '" + str(i) + \
                        "' has a value of " + str(clean[i]) + "!"
                raise ValueError, message

    PSFout = PSF.copy()

    ## remove bad pixels
    if clean[-4]:
        
        PSFout[:] = RemoveBadPSFpixels(PSF=PSFout, center=center,
                exclude_radius=exclude_radius, delta=delta) #@
        
    ## subtract background
    if clean[-3]:

        PSFout[:] = SubtractPSFbackground(PSF=PSFout, 
                background_percent=background_percent, nsigmas=nsigmas) #@
        
        if PSFout.min() < 0:
            print "#DEBUG: PSF negative after background subtraction"
            print "#DEBUG: will do thesholding now -- yes or no ?", 'yes' if clean[-2] else 'no'

                
    ## threshold low PSF values to zero
    if clean[-2]:
        print "#DEBUG: calling   ThresholdPSF_inplace"
        ThresholdPSF_inplace(PSF=PSFout, 
                threshold_percent=threshold_percent, threshold_value=None,
                fill=0.) #@

        if PSFout.min() < 0:
            print "#DEBUG: PSF negative after ThresholdPSF_inplace"

    ## normalize PSF by sum         
    if clean[-1]:
    
        PSFout[:] = NormalizePSF(PSF=PSFout) #@
        
    return PSFout.astype(dtype)
####


#######  FUNCTION: RemoveBadPSFpixels  #######
####    [[checked]]
def RemoveBadPSFpixels(PSF, center=None, exclude_radius=5, delta=0.3):
    """
    Removes bad pixels (cosmic rays, dead/hot pixels) from a 2D _or_ 3D 'PSF'
    by comparison with a median filtered version.  Intensity values that are
    'delta' percent higher than the median filtered PSF will be replaced with
    a median filtered value.  However, points within 'exclude_radius' from the
    'center' of the PSF will not be touched.  If 'center' is 'none', 'center'
    will be determined using the 'FindPSFcentroid' function.
    
    Returns: PSFout
    """

    if center not in (None, 'origin', 'array_center') and \
            PSF.ndim != len(center):
    
        message = "\n'center' coordinates do not match dimension of PSF!"
        raise RuntimeError, message

    if PSF.ndim == 2:
    
        footprint = [[1]*3]*3
    else:
    
        footprint = [[[1]*3]*3]
        
    PSF_mf = U.nd.median_filter(PSF, footprint=footprint, mode="reflect")

    if center is None:
    
        new_center = FindPSFcentroid(PSF=PSF_mf, xyCentroidSize=51) #@
    elif center == 'origin':
    
        new_center = (0,)*PSF.ndim
    elif center == 'array_center':
    
        new_center = tuple(N.array(PSF.shape)/2)
    else:
    
        new_center = tuple(center)

    threshold = (1. + delta)*PSF_mf

    outliers = N.asarray(N.where(PSF > threshold))
    
    if outliers.shape[-1] > 0:
    
        dist2 = N.sum((N.transpose(outliers) - new_center)**2, -1)
        exclude_radius2 = exclude_radius**2
        PSFout = PSF.copy()

        for i in range(len(dist2)):
    
            if dist2[i] > exclude_radius2:
        
                p = tuple(outliers[:,i])
                PSFout[p] = PSF_mf[p]

        return PSFout
    else:
    
        return PSF
####


#######  FUNCTION: SubtractBackground  #######
####    [checked]
def SubtractPSFbackground(PSF, background_percent=0., nsigmas=0.):
    """
    Returns a background subtracted PSF.  Level subtracted is either
    'background_percent' of the PSF maximum _or_ 'nsigmas' of Gaussian 
    background noise (calculated from negative pixels) subtracted.
    'nsigmas' subtraction takes precedence over 'background_percent'

    Returns:
    (PSF - background)
    """

    if background_percent != 0 and nsigmas != 0:
        message = "\nonly one of 'background_percent' and 'nsigmas' can be non-zero"
        raise ValueError, message

    if background_percent == 0 and nsigmas == 0:
    
        return PSF

    elif nsigmas > 0:
        if PSF.min() >= 0:

            message = "WARNING: 'nsigma' > 0 for PSF cleaning but PSF does not have any " +\
                      "negative values to clean\nPlease CHECK that your PSF obeys AIDA's " +\
                      "noise model (with negative pixels!)\n" +\
                      "Proceeding without subtracting background from PSF"
 #                     "Going to proceed with " +\
 #                     "subtracting off 'nsigma'*(theoretical minimal 'sigma_det'=sqrt(pi/2)) " +\
 #                     "from your PSF in cleaning step"
 #           print message
 #           print U.mmms(PSF)
 #           print nsigmas*N.sqrt(N.pi/2)
 #           print U.mmms(PSF - nsigmas*N.sqrt(N.pi/2))
 #           return (PSF - nsigmas*N.sqrt(N.pi/2))
            return (PSF)
        else:
            
            neg_pixels = N.extract(N.less_equal(PSF, 0), PSF)
            ## take the mean, ignoring the extreme most points
            mean_neg_pixels = (N.sort(neg_pixels)[1:-1]).mean()
            sigma = N.sqrt(0.5*N.pi*(mean_neg_pixels*mean_neg_pixels))

            if PSF.max() > nsigmas*sigma:

                return (PSF - nsigmas*sigma)
            else:

                print "WARNING: 'nsigmas'*sigma is larger than max(PSF)!"
                print "Proceeding w/o subtracting off 'nsigmas'*sigma"

    elif background_percent > 0:
        
        background = PSF.max()*background_percent*0.01
        
        return (PSF - background)
    else:

        message = "\n'background_percent' and 'nsigmas' must be > 0!"
        raise ValueError, message
####


#######  FUNCTION: ThresholdPSF_inplace  #######
####    [checked]
def ThresholdPSF_inplace(PSF, threshold_percent, threshold_value=None, fill=0.):
        ## used in 'AIDA_Functions.py'
    """
    Returns a PSF that is thresholded with 'threshold_percent'*PSF.max() _or_
    with 'threshold_value'.  Values below this threshold are replaced with
    'fillvalue'

    now changes PSF in-place
    old version:   Returns:
    old version:   N.where(PSF < threshold_value, fill, PSF)
    """
    
    if threshold_value is None:
    
        threshold_value = threshold_percent*PSF.max()
    
    if threshold_value < 0:
        print "#DEBUG: warning: threshold_percent*PSF.max() resulted in negative threshold", threshold_value
        print "#DEBUG: warning: threshold -- PSF.max() =", PSF.max()
        threshold_value = 0

    PSF[ PSF < threshold_value ] = fill
#    print U.mmms(PSF)
####


#######  FUNCTION: NormalizePSF  #######
####    [checked]
def NormalizePSF(PSF):
    """
    Normalizes 'PSF' by the sum
    
    return (PSF / normalization)
    """
    normalization = float(N.sum(PSF.flat))
        
    if normalization <= 0.:
        
        ## adjustment in cases where large noisy background regions
        ## cause sum of PSF to be negative
        
        normalization = float(
                N.sum((PSF - (normalization/len(PSF.flat)) - 1.).flat))
        
    return (PSF / normalization)
####


'''unused seb 20091118
#######  FUNCTION: PSF2cOTF  #######
####    [checked]
def PSF2OTF(PSF, OTF_threshold_percent=1e-5, fill=0., fullFT=False):
    """
    Returns an OTF of an origin-centered PSF, after thresholding

    Returns: cOTF
    """
    
    if PSF.ndim not in (2,3):

        message = "\nPSF must be 2d or 3d!"
        raise RuntimeError, message

    if PSF.dtype.type in (N.int32, N.int64, N.float32):
    
        dtype = N.complex64
    elif PSF.dtype.type == N.float64:
    
        dtype = N.complex128
    
    else:   
        message = "\n'dtype' of 'PSF' is not correct!"
        raise TypeError, message

    if fullFT:
    
        OTF = N.empty(shape=PSF.shape, dtype=dtype)
    else:

        OTF = N.empty(shape=(PSF.shape[:-1]+(PSF.shape[-1]/2 + 1,)), 
                dtype=dtype)
    
        fftw.rfft(a=PSF, af=OTF, inplace=False)
        OTF.flat[0] = 1.        ## ensure DC value is 1.

        ## threshold OTF values below certain threshold percent
        return CleanRealOTF(realOTF=OTF, cutoff_radius=None, 
                OTF_threshold_percent=OTF_threshold_percent, fill=fill) #@
####
'''#unused seb 20091118

#######  FUNCTION: CleanRealOTF  #######
####
def CleanRealOTF(realOTF, cutoff_radius=None, OTF_threshold_percent=1e-5, 
        fill=0.): # used in 'AIDA_Functions.py'
    """
    Given a 'realOTF', this sets values below 'OTF_threshold_percent' and
    pixels outside 'cutoff_radius' from the origin to 'fill' values

    Returns: cOTF
    """
    
    if cutoff_radius is None:

        cutoff_radius = realOTF.shape[-1] - 1

    ## clean up outside a circle (2D) or cylinder (3D)
    ## assumes realOTF is a single 2D array, a stack of 2D arrays, or
    ## a single 3D array
    OTFout = N.where(RadialArray(realOTF.shape, lambda r:r, origin=0., 
            wrap=(1,0)) > cutoff_radius, fill, realOTF) #@

    return N.where(N.abs(OTFout) < OTF_threshold_percent*OTFout.flat[0], fill,
            OTFout)
####


#######  FUNCTION: RadiallyAverage2D  #######
####    [checked]
def FitOTFzero(OTF, npoints=20): # used in 'AIDA_Functions.py'
    """
    Returns a polynomial fit value for the zero-value of 'OTF'
    """

    if OTF.ndim == 3 and OTF.shape[0] == 1:
    
        OTF.shape = OTF.shape[-2:]
    
    if OTF.ndim == 2:
    
        radialMTF = RadiallyAverage2D(N.abs(OTF), FT=True, origin=None, 
                wrap=(False,), subtract_background=False, fill=0.)[0]
    else:

        radialMTF = RadiallyAverage3D(N.abs(OTF), FT=True, origin=None, 
                wrap=(False,), subtract_background=False, fill=0.)[0]
        radialMTF = RadiallyAverage2D(radialMTF, FT=True, origin=None, 
                wrap=(False,), subtract_background=True, fill=0.)[0]


    ## cubic fit using the first 19 points
    
    intercept = U.fitPoly(radialMTF[1:npoints], (1,-1,-1,-1))[0][0] 
    
    if intercept < radialMTF[1]: # zero-pt fit intercept is below next point
                                 # try other fits, choose closest larger value
        test = N.empty(shape=(3,), dtype=N.float)
        test[0] = U.fitPoly(radialMTF[1:npoints], (1,-1,-1,-1,-1))[0][0]
                # quartic
        test[1] = U.fitPoly(radialMTF[1:npoints], (1,-1,-1))[0][0] # quadratic
        test[2] = U.fitPoly(radialMTF[1:npoints], (1,-1))[0][0] # line
    
        try:
        
            intercept = N.extract(test > radialMTF[1], test).min()
        except:
        
            intercept = radialMTF[1]*1.01

    return intercept
####


#######  FUNCTION: RadiallyAverage2D  #######
####    [checked]
def RadiallyAverage2D(array, FT=True, origin=None, wrap=(False,), 
        subtract_background=True, fill=0.): # used in 'AIDA_Functions.py'  
    """
    Computes 2D radial average of an array up to radius 'array.shape[-2]/2+1'
    Beyond that, values are set to 'fill'.  'subtract_background' will subtract
    an average 'array-background' value computed from regions beyond this
    radius from the entire radial average, before filling with value 'fill'
    
    Returns:
    (radial_mD[...,:array.shape[-2]/2+1], radial_nD, array_background)
    """

    if array.ndim not in (2,3):

        message = "\nDimension of 'array' is neither 2D or 3D!"
        raise ValueError, message

    if array.ndim == 2:
    
        array.shape = (1,) + array.shape
    
    if array.dtype.type in (N.int32, N.int64, N.float32, N.float64):
    
        dtype = N.float32
    else:
    
        dtype = N.complex64
    
    if FT:
        
        if array.shape[-1] == array.shape[-2]/2 + 1:    ## real FT
        
            origin = 0.
            wrap = (True,)*(array.ndim-1) + (False,)

        else:   ## full FT
        
            origin = 0.
            wrap = (True,)*array.ndim
    elif origin is None:
    
        origin = N.array(array.shape)/2
        wrap = wrap*array.ndim
    
    ## N.B.  use a 2D distance_array for _averaging_ even for 3D data (~PSF)
    ## since z-direction averaging should be handled differently and averaged 
    ## across the xy-plane.  This distance_array should have the same shape as
    ## the array input
    distance_array = RadialArray(array.shape[-2:], lambda r:r, origin=origin,
            wrap=wrap[-2:]) #@

    array_background = 0.

    ## use values beyond a circle of radius = length to estimate the
    ## background noise statistics of MTF
    for slice in array:

        array_background += (N.extract(distance_array > array.shape[-2]/2.+1,
                slice)).mean()

    array_background /= array.shape[-3]

    length_z = array.shape[-3]
    length = N.array(array.shape[-2:]).max()
    radial_mD = N.zeros(shape=((length_z,) + (length,)), dtype=dtype)
    radial_nD = N.zeros(shape=array.shape, dtype=dtype)

    for z in range(length_z):
        
        radial_mD[z] = AverageRadially(array[z], origin=origin, wrap=wrap[-2:])
                #@
        radial_nD[z] = RadialArray(array.shape[-2:], 
                lambda r:radial_mD[z][N.round(r).astype(N.int)], origin=origin,
                wrap=wrap[-2:]) #@
                
    ## set all values outside of length cutoff or less than array_background
    ## to floor value
    if subtract_background:

        radial_nD = N.where(N.abs(radial_nD) > 0, 
                radial_nD - array_background, fill)
        radial_mD = N.where(N.abs(radial_mD) > 0, 
                radial_mD - array_background, fill)

    if length_z == 1:
        
        #seb radial_mD.resize(length)
        #seb radial_nD.nd.resize(array.shape[-2:])
        radial_mD = N.resize(radial_mD, length)
        radial_nD  = N.resize(radial_nD, array.shape[-2:])
            
    return (radial_mD[...,:array.shape[-2]/2+1], radial_nD, array_background)
####


#######  FUNCTION: RadiallyAverage3D  #######
####    [checked]
def RadiallyAverage3D(array, FT=True, origin=None, wrap=(False,), 
        subtract_background=True, fill=0.): # used in 'AIDA_Functions.py' 
    """
    Returns 2D and 3D radial average arrays along with estimate of array
    background calculated from entries outside of 1/2 length of array in the y
    dimension.  First calculates radial average about xy, then averages across
    z dimension separately.
    
    If 'FT' is set true, average across z dimension is computed with appropriate
    index matching for a Fourier Transform array
    
    Returns:
    (radial_2D, radial_3D, array_background)
    """

    if array.ndim != 3 or array.shape[-3] == 1:

        message = "\nDimension of 'array' must be 3D!"
        raise ValueError, message

    length_z = array.shape[-3]/2

    (radial_2D, radial_3D, array_background) = RadiallyAverage2D(array, FT=FT, 
            origin=origin, wrap=wrap[-2:], 
            subtract_background=subtract_background, fill=fill) #@

    if FT:
    
        ## average 2 halves over z appropriately for Fourier Transform array
        radial_2D[1:length_z] = radial_2D[-length_z+1:][::-1] = \
                0.5*(radial_2D[1:length_z] + radial_2D[-length_z+1:][::-1])
        radial_3D[1:length_z] = radial_3D[-length_z+1:][::-1] = \
                0.5*(radial_3D[1:length_z] + radial_3D[-length_z+1:][::-1])
    else:
    
        ## average simply as even or odd shaped array in z
        radial_2D[:length_z] = radial_2D[-length_z:][::-1] = \
                0.5*(radial_2D[:length_z] + radial_2D[-length_z:][::-1])
        radial_3D[:length_z] = radial_3D[-length_z:][::-1] = \
                0.5*(radial_3D[:length_z] + radial_3D[-length_z:][::-1])
    
    return (radial_2D, radial_3D, array_background)
####


#######  FUNCTION: AverageRadially  #######
####    [checked]
def AverageRadially(array, origin=0., wrap=(True,)):
    """
    Returns a n-1 D radial average of array
    
    Returns:
    (histogram / frequency)
    """
    
    distance_array = RadialArray(array.shape, lambda r:r, origin=origin, 
            wrap=wrap[-2:])
    
    nbins = N.array(array.shape).max()
    histogram = N.empty(shape=nbins, dtype=array.dtype)
    ## uses Sebastian Haase's general histogram functions
    frequency = U.histogram(distance_array, nBins=nbins, amin=0, amax=nbins)
    ## to prevent division by zero
    eps = 1./N.product(array.shape)
    frequency = N.where(frequency==0, eps, frequency)
        
    histogram = U.generalhistogram(distance_array, array, nBins=nbins, amin=0,
            amax=nbins)

    return histogram / frequency
####
    

#######  FUNCTION: RadialArray  #######
####    [checked]
def RadialArray(shape, lambda_function, origin=0., wrap=(True,)):
    """
    Generates and returns a radially symmetric function sampled in 
    volume(image) of shape
    if 'origin' is 'None' the origin defaults to the center
    func is a 1D function with 2 paramaters: r0,r
    dtype is Float32
    """

    if len(wrap) == 1:
    
        wrap = N.ones(shape=len(shape)) * wrap[0]
    elif len(wrap) > len(shape):
    
        message = "\n'wrap' needs to be a list of True or False values " + \
                "that match length of 'shape'"
        raise ValueError, message
    elif len(shape) == 3 and shape[-3] == 1 and len(wrap) == 2:
    
        wrap = (1,) + wrap
    else:
    
        wrap = (0,) + wrap
        
    if origin is None:

        origin = N.array(shape) / 2.
    else:

        try:

            origin = N.ones(shape=len(shape)) * float(origin)
        except:

            pass

    if len(shape) != len(origin):

        message = "\n'shape' and 'origin' not same dimension"
        raise ValueError, message

    if len(shape) == 1:

        x0 = origin
        nx = shape[-1]

        return N.fromfunction(
                lambda x: lambda_function(
                WrapIt(N.abs(x-x0), nx, wrap[-1])), shape, dtype=N.float32) #@
    elif len(shape) == 2:

        (y0, x0) = origin
        (ny, nx) = shape
        
        return N.fromfunction(
                lambda y,x: lambda_function(N.sqrt((
                WrapIt((x-x0), nx, wrap[-1]))**2 + \
                (WrapIt((y-y0), ny, wrap[-2]))**2) ), shape, dtype=N.float32) #@
    elif len(shape) == 3:

        (z0, y0, x0) = origin
        (nz, ny, nx) = shape

        return N.fromfunction(
                lambda z,y,x: lambda_function(N.sqrt(
                (WrapIt((x-x0), nx, wrap[-1]))**2 + \
                (WrapIt((y-y0), ny, wrap[-2]))**2 + \
                (WrapIt((z-z0), nz, wrap[-3]))**2) ), shape, dtype=N.float32) #@
    else:

        message = "\nSorry - only works for arrays with dimension <= 3!"
        raise RuntimeError, message
####


#######  FUNCTION: WrapIt  #######
####    
def WrapIt(q, nq, wrap):
    """
    Used for lambda type functions to wrap array, 'q', of size 'nq'.  If 
    'wrap' = 1, then wrap, otherwise, don't wrap.  Works on a per axis basis
    """
    
    if wrap:

        return N.where(q > nq/2, q-nq, q)
    else:

        return q
####


#######  FUNCTION:  Output2File  #######
####    [checked]
def Output2File(data_array, filebase, format, hdr=None, shape=None):
        # used by 'AIDA_Functions.py'
    """
    Writes out 'data_array' as a '.fits' or '.mrc' file or '.tif'files

    For format='m', writes out 'data_array' in .mrc format, including a header
        containing array 'shape' information
    For format='f', writes out 'data_array' in .fits format, with shape
        information stored in the [1].data  slot (default)
    For format='t', writes out 'data_array' in .tif format
	
	Also outputs PNG versions of the results, for easy viewing.
    """
    
    # below is old
    #if shape is None:
    #
    #    shape = data_array.shape
    
    ### EHom (20130625): adding line to shape data_array according to shape input parameter
    ### Should have been here before
    if (shape != None):
        data_array.shape = shape
        
    import matplotlib.pyplot as plt
    #plt.figure()
    #plt.imshow(data_array)
    #plt.title(data_array[0,0])
    #plt.show()
    
    if format == 'm':

        Mrc.save(data_array, filebase + '.mrc', ifExists="overwrite")
        
        # below is old way - Mrc.bindArr no longer exists in Priithon
        #rs = ''
        #
        #for i in shape:
        #   
        #    rs += '%d ' %i
        #
        #dtype = data_array.dtype
        #
        #temp = Mrc.bindArr(filebase + '.mrc', data_array.astype(N.float32))
        ## can only write out as single precision
        #fileheader = temp.Mrc.hdrArray[0]
        #fileheader.setfield('NumTitles',1)
        #fileheader.field('title')[0] = 'Shape: ' + rs
        #temp.Mrc.close()
        ## STILL NEED TO PROVIDE A WAY OF SETTING HEADER INFO FROM INPUT
        
    elif format == 'f':

        if os.path.exists(filebase + '.fits') == 1:

            os.remove(filebase + '.fits')

        # Clement: using astropy.io.fits now
        
        fits_file = iofits.HDUList()
        datahdu = PrimaryHDU()
        datahdu.data = data_array
        
        
        iofits.append(filebase + '.fits',data_array,header=hdr)
        
    elif format == 't':
        if os.path.exists(filebase + '.tiff') == 1:

            os.remove(filebase + '.tiff')
        
        img = scipy.misc.toimage(data_array)
        img.save(filebase + '.tiff')
        
    elif format == 't2':
        if os.path.exists(filebase + '.tif') == 1:

            os.remove(filebase + '.tif')
        
        img = scipy.misc.toimage(data_array)
        img.save(filebase + '.tif')
        
# Clement: Old version using pyfits (deprecated)
#         fits_file = pyfits.HDUList()
#         datahdu = pyfits.PrimaryHDU()
#         datahdu.data = data_array
#         
#         ## STILL NEED TO PROVIDE A WAY OF SETTING HEADER INFO FROM INPUT
#         #if type(hdr) is not types.NoneType:
#         #
#         #    datahdu.header = hdr
#         #    
#         #   print hdr
#         
#         # Provide header info from the original fits file.
#         
#         
#         fits_file.append(datahdu)
#         fits_file.writeto(filebase + '.fits')
        
#     else:   # format must be .tiff
#     
#         #!!!! TENTATIVE !!!!
#         # make sure orientation of TIFF file matches convention
#         if len(data_array.shape) == 2:
#         
#             U.saveImg(data_array[...,::-1,...], filebase + ".tiff")
#         elif len(data_array.shape) == 3:
#         
#             U.saveImg_seq(data_array[...,::-1,...], filebase + ".tiff")
#         else:
#         
#             message = "\n'data_array' shape is not 2 or 3!  Cannot write " + \
#                     "out TIFF file!"
#             raise ValueError, message

    ### EHom (20130616): also output results (if 2D) as an 8-bit JPEG files using PIL
    ### In the division of 255, I hack the addition of a small value to avoid 
    ### a divide by zero in a true_divide call
    if len(data_array.shape) == 2:

        min = data_array.min()
        max = data_array.max()
        #print data_array.min()
        #print data_array.max()
        #print data_array.mean()
        rescaled = N.where(data_array > min, data_array-min, 0.)
        if ((max - min) == 0):
            message = "\nMax Min problem in outputting array!  Cannot write JPEG file\n"
            print message
        else:
            rescaled *= (255.0 / (max - min))
            # Clement: we don't need to save the jpeg
            # im = ImageOps.flip(Image.fromarray(rescaled.astype(N.uint8)))
            # rescale and flip vertically to properly register image with FITS output
            # im.save(filebase + '.jpeg')
####


#######  FUNCTION: LoadFile  #######
####    [checked]
def LoadFile(filename, dtype=N.float64): # used by 'AIDA_Functions.py'

    if os.path.isfile(filename) != 1:

        message = "\nFile: " + "'" + str(filename) + "' does not exist!"
        raise RuntimeError, message

    ext = string.split(filename,sep='.')[-1]

    if ext.lower() == 'fits':

        ## Load as simple FITS file from Primary HDU
        # Clement: New version usin astropy.io.fits
        fits = iofits.open(filename)
        fits[0].verify(option='silentfix')
        data = fits[0].data.astype(dtype)
        hdr = fits[0].header.copy()
        
        import matplotlib.pyplot as plt
        #plt.figure()
        #plt.imshow(data)
        #plt.title(data[0,0])
        #plt.show()
        
    
    elif ext.lower() == 'tiff' or ext.lower() == 'tif':
        
        im = Image.open(filename)
        
        data_tmp = N.asarray(im).astype(dtype).copy() + 0.000001
        data = data_tmp.astype(N.float64)
        import matplotlib.pyplot as plt
        #plt.figure()
        #plt.imshow(data)
        #plt.title(data[0,0])
        #plt.show()
        
        hdr = []
        
        
# Clement: Old version using pyfits
#         fits = pyfits.open(filename)
#         fits[0].verify(option='silentfix')  ## fix data & headers if necessary
#         data = fits[0].data.astype(dtype)
#         hdr = fits[0].header.copy()
    elif ext.lower() == 'mrc':

        data = Mrc.bindFile(filename).astype(dtype)
        hdr = []
    
    ## FIX!! - NEED TO FIGURE OUT HOW TO ACCESS AND SAVE HEADER INFO
    
    else:   ## default, if extension is not ".fits" or ".mrc"
    
        try:
        
            data = Mrc.bindFile(filename).astype(dtype)
            hdr = []
        except:
        
            message = "\nFile is not a .fits or .mrc file!  Cannot open!"
            raise RuntimeError, message
    
    ## compress data if one of the axes is of dimension size 1
    shape = N.array(data.shape)
    data.shape = tuple(N.compress(shape > 1, shape))

    return (data, hdr)
####


#######  FUNCTION: DetermineFileFormat  #######
####    [checked]
def DetermineFileFormat(directory, filenamebase): # used in 'AIDA_Functions.py'
    """
    Changed by Clement 2016-08-16 : handling of TIFF files added
    Returns the file format of a file(s) within 'directory' identified by 
    'filenamebase'.
    
    Assumes the file(s) is either in .fits ('f') or .mrc ('m') or >tif ('t') format and will
    select the format that characterizes the majority of files that are located;
    if the same number of .fits and .mrc files are identified, the format will
    default to 'f'
    """

    fits_check = len(PickFiles(directory, prefix=filenamebase, 
            suffix='.fits')) #@
    mrc_check = len(PickFiles(directory, prefix=filenamebase, suffix='.mrc')) 
    
    tif_check = len(PickFiles(directory, prefix=filenamebase, suffix='.tif')) 
    tiff_check = len(PickFiles(directory, prefix=filenamebase, suffix='.tiff')) 

    if fits_check >= mrc_check and fits_check >= tif_check and fits_check >= tiff_check:
    
        format = 'f'
    elif tiff_check >= mrc_check and tiff_check > tif_check:
        format = 't' 
    elif tif_check >= mrc_check:
        format = 't2' 
    elif mrc_check >= 0:

        format = 'm'
    else:

        message = "\nNo appropriate file(s) found with filename base:\n\t" + \
                "'" + str(filenamebase) + "'\n" + \
                "(should be a .mrc or .fits file)"
        raise RuntimeError, message
     

    return format
####


#######  FUNCTION: PickFiles  #######
####    [checked]
def PickFiles(directory, prefix='psf', suffix='.fits'):
    # used in 'AIDA_Functions.py'
    """
    Returns list of files in directory with specified 'prefix' and 'suffix'
    """
    

    allfiles = os.listdir(directory)
    allfiles.sort()
    filelist = []
    prefix_length = len(prefix)
    suffix_length = len(suffix)

    for file in allfiles:

        if file[-suffix_length:] == suffix:

            if file[:prefix_length] == prefix:
            
                filelist.append(file)

    return filelist
####


#######  FUNCTION: ResizeImage  #######
####    [checked]
def ResizeImage(image, dimension=2, zresize=False, dtype=N.float64):
    # used in 'AIDA_Functions.py'
    """
    Resizes (crops or pads) images so that dimensions are set to the closest
    length of 2^n
    
    Returns: new_image
    """
    
    size = N.array((4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
    shape = tuple(N.compress(N.array(image.shape) > 1, N.array(image.shape)))

    if len(shape) < dimension:
    
        message = "\n'dimension' is larger than the shape if 'image'!"
        raise ValueError, message

    if dimension > 3 or dimension < 2:
    
        message = "\n'image' must be either 2D or 3D!'"
        raise RuntimeError, message
    
    image.shape = shape
    new_shape = N.array(shape)
    new_shape[-2] = size[ U.findMin(N.abs(shape[-2] - size))[-1] ]
    new_shape[-1] = size[ U.findMin(N.abs(shape[-1] - size))[-1] ]
        
    if dimension == 3 and zresize:
        
        new_shape[-3] = size[ U.findMin(N.where((shape[-3] - size) >= 0,
                (shape[-3] - size), 9999999))[-1] + 1 ]     ## round up 

    offset = (N.array(shape) - N.array(new_shape))/2.
    new_image = N.empty(shape=tuple(N.maximum(shape, new_shape)), dtype=dtype)
    lower = N.empty(shape=(dimension,), dtype=N.int)
    upper = N.empty(shape=(dimension,), dtype=N.int)
    
    for i in range(1,dimension+1):
    
        lower[-i] = int(N.floor(N.abs(offset[-i])))
        upper[-i] = lower[-i] + N.minimum(shape[-i], new_shape[-i])
        
        if offset[-i] > 0:      ## crop
        
            image = N.take(image, indices=range(lower[-i], upper[-i]),
                    axis=-i)    
    
    if dimension == 2:
    
        if len(shape) == 3:

            new_image[:,lower[-2]:upper[-2], lower[-1]:upper[-1]] = image       
        else:
        
            new_image[lower[-2]:upper[-2], lower[-1]:upper[-1]] = image
    else:
    
        new_image[lower[-3]:upper[-3], lower[-2]:upper[-2], 
                lower[-1]:upper[-1]] = image    
    
    for i in range(1,dimension+1):
        
        if offset[-i] > 0:
        
            new_image = N.take(new_image, indices=range(lower[-i],
                    lower[-i] + new_shape[-i]), axis=-i)
        else:
        
            new_image = N.take(new_image, indices=range(0, \
                    new_shape[-i]), axis=-i)

    return new_image
####


#######  FUNCTION: SelectKernel  #######
####    [checked]
def SelectKernel(dimension, operator, shift, rho=1., zeta=3., dtype=N.float64):
    # used in 'AIDA_Functions.py'
    """
    Returns the appropriate 'kernel' for use in calculating the gradient
    of an object image via convolution (see CalcGradient function).
    
    Kernels available here correspond to: (1) a simple directed finite-
    difference approach advocated by Mugnier et al. (2001) "Myopic 
    deconvolution from wavefront sensing" JOSA A 18:862-872; (2) a symmetric
    finite-difference approach as described by Press et al. (1992) Numerical 
    Recipes (2nd Ed) pp.186-188; or (3) an isotropic differencing scheme
    suggested by Frei and Chen (1977) "Fast boundary detection: a generalization
    and a new algorithm" IEEE Trans. Computers C-26:988
    
    'dimension' specifies the dimension of the data for which a kernel should be
            returned (currently '2' or '3')
    'operator' specifies the type of gradient kernel to select.  Options are:
            'pixel' (directed finite difference)
            'symmetric' (symmetric finite difference)
            'FC' (Frei-Chen isotropic finite difference scheme)
    'shift' is the finite-difference orientation.
            Default is '-1' meaning: object(x) - object(x-1)  ("backwards")
            Alternatively for '1':   object(x) - object(x+1)  ("forwards")
    """

    if dimension == 2:

        if operator == 'pixel':

            if shift == -1:

                kernel = N.array(( ((0.,0.,0.), (0.,rho,0.), (0.,-rho,0.)), 
                        ((0.,0.,0.), (0.,rho,-rho), (0.,0.,0.)) ), dtype=dtype)
                        # (y,x)
            elif shift == +1:

                kernel = N.array(( ((0.,-rho,0.), (0.,rho,0.), (0.,0.,0.)), 
                        ((0.,0.,0.), (-rho,rho,0.), (0.,0.,0.)) ), dtype=dtype)
                        # (y,x)
        elif operator == 'symmetric':

            a = rho*0.5
            
            if shift == -1:

                kernel = N.array(( ((0.,a,0.), (0.,0.,0.), (0.,-a,0.)), 
                        ((0.,0.,0.), (a,0.,-a), (0.,0.,0.)) ), dtype=dtype)
            elif shift == +1:

                kernel = N.array(( ((0.,-a,0.), (0.,0.,0.), (0.,a,0.)), 
                        ((0.,0.,0.), (-a,0.,a), (0.,0.,0.)) ), dtype=dtype)
        elif operator == 'FC':

            a = 1./(2. + (N.sqrt(2.)*rho) )
            a_sqrt2 = a*N.sqrt(2.)*rho

            if shift == -1:

                kernel = N.array(( ((a,a_sqrt2,a), (0.,0.,0.),
                        (-a,-a_sqrt2,-a)), ((a,0.,-a), (a_sqrt2,0.,-a_sqrt2),
                        (a,0.,-a)) ), dtype=dtype)
            elif shift == +1:

                kernel = N.array(( ((-a,-a_sqrt2,-a), (0.,0.,0.),
                        (a,a_sqrt2,a)), ((-a,0.,a), (-a_sqrt2,0.,a_sqrt2),
                        (-a,0.,a)) ), dtype=dtype)
        elif operator == 'laplacian':

            if shift == 1:

                a = -1./(4.*rho)
                kernel = N.array(( (0.,a,0.), (a,1.,a), (0.,a,0.) ),
                        dtype=dtype)
            elif shift == 2:

                a = -1./(8.*rho)
                kernel = N.array(( (a,a,a), (a,1.,a), (a,a,a) ), dtype=dtype)
            elif shift == 3:

                a = -1./(4 + (8*rho) )
                b = 2.*rho*a
                kernel = N.array(( (a,b,a), (b,1.,b), (a,b,a) ), dtype=dtype)
            else:

                message = "\n'laplacian_operator' must be 0, 1, 2, or 3!"
                raise ValueError, message
                
    elif dimension == 3:

        if operator == 'pixel':

            if shift == -1:

                kernel = N.array(( (((0.,0.,0.),(0.,0.,0.),(0.,0.,0.)),  ## k_z
                                ((0.,0.,0.),(0.,zeta*rho,0.), (0.,0.,0.)),
                                ((0.,0.,0.),(0.,-zeta*rho,0.),(0.,0.,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_y
                                ((0.,0.,0.), (0.,rho,0.), (0.,-rho,0.)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_x
                                ((0.,0.,0.), (0.,rho,-rho), (0.,0.,0.)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))) ),
                                dtype=dtype)
            elif shift == +1:

                kernel = N.array(( (((0.,0.,0.),(0.,-zeta*rho,0.),(0.,0.,0.)),
                                                                        ## k_z
                                ((0.,0.,0.),(0.,zeta*rho,0.), (0.,0.,0.)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_y
                                ((0.,-rho,0.), (0.,rho,0.), (0.,0.,0.)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_x
                                ((0.,0.,0.), (-rho,rho,0.), (0.,0.,0.)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))) ), 
                                dtype=dtype)
        elif operator == 'symmetric':

            a = rho*0.5
            
            if shift == -1:
                                                                        ## k_z
                kernel = N.array(( (((0.,0.,0.),(0.,zeta*a,0.),(0.,0.,0.)),    
                                ((0.,0.,0.),(0.,0.,0.), (0.,0.,0.)),
                                ((0.,0.,0.),(0.,-zeta*a,0.),(0.,0.,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_y
                                ((0.,a,0.), (0.,0.,0.), (0.,-a,0.)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_x
                                ((0.,0.,0.), (a,0.,-a), (0.,0.,0.)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))) ), 
                                dtype=dtype)
            elif shift == +1:
                                                                        ## k_z
                kernel = N.array(( (((0.,0.,0.),(0.,-zeta*a,0.),(0.,0.,0.)), 
                                ((0.,0.,0.),(0.,0.,0.), (0.,0.,0.)),
                                ((0.,0.,0.),(0.,zeta*a,0.),(0.,0.,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_y
                                ((0.,-a,0.), (0.,0.,0.), (0.,a,0.)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_x
                                ((0.,0.,0.), (-a,0.,a), (0.,0.,0.)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))) ), 
                                dtype=dtype)
        elif operator == 'FC':

            a = 1./(2. + (N.sqrt(2.)*rho) )
            a_sqrt2 = a*N.sqrt(2.)*rho

            if shift == -1:

                kernel = N.array(( (((0.,zeta*a,0.),(0.,zeta*a_sqrt2,0.),
                                            (0.,zeta*a,0)),             ## k_z
                                ((0.,0.,0.),(0.,0.,0.), (0.,0.,0.)),
                                ((0.,-zeta*a,0.),(0.,-zeta*a_sqrt2,0.),
                                            (0.,-zeta*a,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_y
                                ((a,a_sqrt2,a), (0.,0.,0.), (-a,-a_sqrt2,-a)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_x
                                ((a,0.,-a), (a_sqrt2,0.,-a_sqrt2), (a,0.,-a)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0.))) ),
                                dtype=dtype)
            elif shift == +1:

                kernel = N.array(( (((0.,-zeta*a,0.),(0.,-zeta*a_sqrt2,0.),
                                            (0.,-zeta*a,0)),            ## k_z
                                ((0.,0.,0.),(0.,0.,0.), (0.,0.,0.)),
                                ((0.,zeta*a,0.),(0.,zeta*a_sqrt2,0.),
                                            (0.,zeta*a,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_y
                                ((-a,-a_sqrt2,-a), (0.,0.,0.), (a,a_sqrt2,a)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0))),
                              (((0.,0.,0.),(0.,0.,0.),(0.,0.,0)),       ## k_x
                                ((-a,0.,a), (-a_sqrt2,0.,a_sqrt2), (-a,0.,a)),
                                ((0.,0.,0.),(0.,0.,0.),(0.,0.,0.))) ),
                                dtype=dtype)
        elif operator == 'laplacian':

            if shift == 1:
                a = -1./((4.+zeta*2.)*rho)
                kernel = N.array(( ((0.,0.,0.), (0.,zeta*a,0.), (0.,0.,0.)),
                            ((0.,a,0.), (a,1.,a), (0.,a,0.)),
                            ((0.,0.,0.,), (0.,zeta*a,0.), (0.,0.,0.)) ), 
                            dtype=dtype)
            elif shift == 2:

                a = -1./((22.+zeta*2.)*rho)
                kernel = N.array(( ((a,a,a), (a,zeta*a,a), (a,a,a)), 
                            ((a,a,a), (a,1.,a), (a,a,a)),
                            ((a,a,a), (a,zeta*a,a), (a,a,a)) ), dtype=dtype)
            elif shift == 3:

                ## original normalization: a = -1./32.
                a = -1./(16.*(1+zeta)*rho)
                b = 2.*a
                kernel = N.array(( ((a,a,a), (a,zeta*b,a), (a,a,a)), 
                            ((zeta*a,zeta*b,zeta*a), (zeta*b,1.,zeta*b), 
                                    (zeta*a,zeta*b,zeta*a)),
                            ((a,a,a), (a,zeta*b,a), (a,a,a)) ), dtype=dtype)
            else:

                message = "\n'laplacian_operator' must be 0, 1, 2, or 3"
                raise ValueError, message

    return kernel
####


#######  FUNCTION: WienerFilter  #######
####    [checked]
def WienerFilter(IMAGE, OTF, weight): # used in 'AIDA_Functions.py'    
    """
    Returns a Wiener filtered image according to:
    
    object = invFT(conjugate(OTF)*IMAGE / (||OTF||^2 + weight))
    
    'weight' is the scalar Wiener Filter parameter quantity
        Note that this represents a "squared" quantity:
            weight = (sigma_noise/sigma_image)**2
    """

    fullshape = IMAGE.shape[:-1] + ((IMAGE.shape[-1]-1)*2,)

    if IMAGE.dtype == N.complex128:
    
        filtered = N.empty(shape=fullshape, dtype=N.float64)
    elif IMAGE.dtype == N.complex64:
    
        filtered = N.empty(shape=fullshape, dtype=N.float32)
    else:
    
        message = "\n'IMAGE' must be complex dtype!"
        raise ValueError, message

    if IMAGE.dtype != OTF.dtype:
    
        OTF = OTF.astype(IMAGE.dtype)

    conjOTF = N.conjugate(OTF)
    fftw.irfft(af=( (conjOTF*IMAGE / (conjOTF*OTF + weight)) / \
            N.product(fullshape)), a=filtered, inplace=False, copy=False)

    return filtered
####


'''unused seb 20091118
#######  FUNCTION: WienerFilterReal  #######
####    [checked]
def WienerFilterReal(image, PSF, weight):    
    """
    Returns a Wiener filtered image according to:
    
    object = invFT(conjugate(OTF)*IMAGE / (||OTF||^2 + weight))
    
    'weight' is the scalar Wiener Filter parameter quantity
        Note that this represents a "squared" quantity:
            weight = (sigma_noise/sigma_image)**2
    
    This is the version that takes a real image and psf as input
    """

    halfshape = image.shape[:-1] + (image.shape[-1]/2 + 1,)
    fullshape = image.shape
        
    if image.dtype.type == N.float64:
    
        IMAGE = N.empty(shape=halfshape, dtype=N.complex128)
    elif image.dtype.type == N.float32:
    
        IMAGE = N.empty(shape=halfshape, dtype=N.complex64)
    elif image.dtype.type in (N.int32, N.uint32, N.int64, N.uint64):
    
        image = image.astype(N.float32)
        IMAGE = N.empty(shape=halfshape, dtype=N.complex64)
    else:
    
        message = "\n'image' must be an array of numbers!"
        raise ValueError, message

    filtered = N.empty(shape=image.shape, dtype=image.dtype)

    if image.dtype != PSF.dtype:
    
        PSF = PSF.astype(image.dtype)
    
    fftw.rfft(a=image, af=IMAGE, inplace=False)
    OTF = PSF2OTF(PSF, OTF_threshold_percent=1e-5, fill=0., fullFT=False)
    conjOTF = N.conjugate(OTF)
    fftw.irfft(af=( (conjOTF*IMAGE / (conjOTF*OTF + weight)) / \
            N.product(fullshape)), a=filtered, inplace=False, copy=False)

    return filtered
####
'''#unused seb 20091118

######  FUNCTION: SetEdgeValues  #######
####    [checked]
def SetEdgeValues(data, ndim, value=0):
    """
    Returns 'data' with edges along dimensions up to 'ndim' set to 'value'
    
    'ndim' is assumed to be <= data.ndim (no safeguard in place)
    """

    axes = range(1, data.ndim) + [0]    

    for dimension in range(ndim):       ## N.B. first transpose starts
                                        ## positions at fastest index
        data.transpose(axes)
        data[0]=data[-1]=value

    for dimension in range(data.ndim - ndim):    ## loop to transpose to the
                                                 ## original configuration
        data.transpose(axes)
####


######  FUNCTION: Num2Str  #######
####    [checked]
def Num2Str(num, digits=5): # used in 'AIDA_Functions.py'

    form = '%.' + str(int(digits)) + 'g'

    if type(num) in (types.IntType, types.FloatType) or isinstance(num, 
            N.number):
    
        basestr = form %num
        split = string.split(basestr, '.')

        if len(split) == 2:
    
            strg = split[0] + 'p' + split[1]
        elif len(split) == 1:

            strg = split[0]
        else:

            message = "\nError in number to string conversion!"
            raise RuntimeError, message

    elif type(num) in (types.ListType, types.TupleType):
    
        strg = '{'
        
        ## hard limit
        if len(num) > 10:

            for i in (0, len(num)-1):
            
                basestr = form %num[i]
                split = string.split(basestr, '.')

                if len(split) == 2:
    
                    strg += split[0] + 'p' + split[1]
                elif len(split) == 1:

                    strg += split[0]
                else:   

                    message = "\nError in number to string conversion!"
                    raise RuntimeError, message

                if i < len(num) - 1:
            
                    strg += "--"
                else:
            
                    strg += "}"
        else:
        
            for i in range(len(num)):

                basestr = form %num[i]
                split = string.split(basestr, '.')

                if len(split) == 2:
    
                    strg += split[0] + 'p' + split[1]
                elif len(split) == 1:

                    strg += split[0]
                else:   

                    message = "\nError in number to string conversion!"
                    raise RuntimeError, message
    
                if i < len(num) - 1:
            
                    strg += ","
                else:
            
                    strg += "}"
    else:
    
        message = "\n'num' entry is not a number or a list!"
        raise ValueError, message
                
    return strg
####

'''unused seb 20091118
######  FUNCTION: ComputeVarianceMean  #######
####
def ComputeMeanVariance(x, data_dim=None, dtype=N.float32, conjugate=False):
        ## used in 'AIDA_Functions.py'
    """
    Returns variance and mean of a sequence 'x' along the slowest axis 
    using the Knuth-Welford algorithm (similar to West's algorithm)
    
    Uses 'data_dim' to make sure sequence has more than 1 slice
    """

    x = N.asarray(x)
    
    if x.ndim == 0:     # if single integer
    
        return (x, 0)
    elif x.ndim == data_dim:
    
        return (x, N.sqrt(x))      # return variance = poisson noise
    else:
    
        n = 0
        mean = N.zeros(shape=x.shape[-x.ndim+1:], dtype=dtype)
        delta = N.zeros(shape=x.shape[-x.ndim+1:], dtype=dtype)
        S = N.zeros(shape=x.shape[-x.ndim+1:], dtype=dtype)
    
        if conjugate:       ## use for complex 'x' to generate power spectral
                            ## density instead of the straight variance     
            for p in x:
    
                n += 1
                delta = p - mean
                mean += delta / n
                S += delta * N.conjugate(delta)
        else:
        
            for p in x:
    
                n += 1
                delta = p - mean
                mean += delta / n
                S += delta * delta

        return (mean, S/(n-1))  # return mean and variance
####


#######  FUNCTION: ComputeVarianceMean  #######
####
def ComputeMeanVariance3(x):
    """
    Returns variance and mean of a sequence using the pairwise algorithm
    of Chan et al. (1979) "Updating formulae and pairwise algorithm for
    computing sample variance" Stanford Computer Science Report STAN-CS-79-773
    Doesn't seem to work as well as Knuth-Welford, even though it is suppose to!
    Wonder why...
    """
    
    n = len(x)
    
    if x.ndim == 1:
    
        shape = (int(N.round(N.log(n)/N.log(2)) + 1),)
    else:

        shape = (int(N.round(N.log(n)/N.log(2)) + 1),) + x.shape[-x.ndim+1:]

    terms = N.zeros(shape=shape[0], dtype=N.int)
    sum_pair = N.zeros(shape=shape, dtype=N.float64)
    S_pair = N.zeros(shape=shape, dtype=N.float64)
    
    terms[0] = 0
    top = 1
    
    for i in range(0,n,2):
    
        ## compute sum and sum of sqaures for the next 2 data points
        ## in sequence; put these quantities on top of stack
        sum_pair[top] = x[i] + x[i+1]
        S_pair[top] = 0.5*(x[i+1] - x[i])**2 
        terms[top] = 2
        
        while terms[top] == terms[top-1]:
            ## top 2 elements on stack contain quantities computed from the
            ## same number of data points; combine them
            top -= 1
            terms[top] *= 2
            S_pair[top] += S_pair[top+1] + \
                    ((sum_pair[top] - sum_pair[top+1])**2) / terms[top]
            sum_pair[top] += sum_pair[top+1]
            
        top += 1
    
    top -= 1
    
    if IsOdd(n):
    
        # n is odd; put last point on stack
        top += 1
        terms[top] = 1
        sum_pair[top] = x[n]
        S_pair[top] = 0.
    
    T = terms[top]
    Nsum = sum_pair[top]
    NS = S_pair[top]
    
    if top >= 2:

        ## n is not a power of 2; the stack contains more than one element;
        ## combine them
        for j in range(2, top):
        
            i = top + 2 - j
            NS += S_pair[i] + (T*(terms[i]*Nsum/T - sum_pair[i])**2) / \
                    (terms[i]*(terms[i] + T))
            Nsum += sum_pair[i]
            T += terms[i]
    
    return (Nsum/T, NS/(T-1))   # return mean and variance
####


######  FUNCTION: v  #######
####    [checked]
def v():

    Y.view(LoadFile(FN())[0])
####


#######  FUNCTION: l  #######
####    [checked]
def l(dtype=None):   

    if dtype is None:
    
        return LoadFile(FN(), N.float64)[0]
    else:
    
        return LoadFile(FN(), dtype)[0]

####


######  FUNCTION: ArrayCenter2Origin  #######
####    [checked]
def ArrayCenter2Origin(array2D):

    shape = array2D.shape
    
    if len(array2D.shape) != 2 or shape[0] == 1:
    
        message = "\nSorry - works only for 2D arrays at the moment!"
        raise RuntimeError, message
    
    (y0, x0) = N.array(shape)/2
    newarray = N.empty(shape=shape, dtype=array2D.dtype)
    
    newarray[0:y0, 0:x0] = array2D[y0:, x0:]
    newarray[y0:, 0:x0] = array2D[0:y0, x0:]
    newarray[0:y0, x0:] = array2D[y0:, 0:x0]
    newarray[y0:, x0:] = array2D[0:y0, 0:x0]
    
    return newarray
####
'''#unused seb 20091118

######  FUNCTION: Origin2ArrayCenter  #######
####    [checked]
def Origin2ArrayCenter(array2D):

    shape = array2D.shape
    
    if len(array2D.shape) != 2 or shape[0] == 1:
    
        message = "\nSorry - works only for 2D arrays at the moment!"
        raise RuntimeError, message
    
    (y0, x0) = N.array(shape)/2
    newarray = N.empty(shape=shape, dtype=array2D.dtype)
    
    newarray[y0:, x0:] = array2D[0:y0, 0:x0]
    newarray[0:y0, x0:] = array2D[y0:, 0:x0]
    newarray[y0:, 0:x0] = array2D[0:y0, x0:]
    newarray[0:y0, 0:x0] = array2D[y0:, x0:]
    
    return newarray
####


######  FUNCTION: ArrayCenter2Origin  #######
####    [checked]
def ArrayCenter2Origin(array2D):

    shape = array2D.shape
    
    if len(array2D.shape) != 2 or shape[0] == 1:
    
        message = "\nSorry - works only for 2D arrays at the moment!"
        raise RuntimeError, message
    
    (y0, x0) = N.array(shape)/2
    newarray = N.empty(shape=shape, dtype=array2D.dtype)
    
    newarray[0:y0, 0:x0] = array2D[y0:, x0:]
    newarray[y0:, 0:x0] = array2D[0:y0, x0:]
    newarray[0:y0, x0:] = array2D[y0:, 0:x0]
    newarray[y0:, x0:] = array2D[0:y0, 0:x0]
    
    return newarray
####


######  FUNCTION: realFT2fullFT  #######
####
def realFT2fullFT(rFT):     ## used in 'AIDA_Functions.py'
    """
    Converts a real Fourier Transform array to a full Fourier Transform
    """

    for i in range(len(rFT.shape[:-1])):
    
        if IsOdd(rFT.shape[i]): #@
        
            message = "\nrFT array must be of an even shaped array!"
            raise RuntimeError, message

    temp = N.empty(shape=(rFT.shape[:-1] + (rFT.shape[-1]-2,)), 
            dtype=rFT.dtype)

    if rFT.ndim == 2:

        temp[0,:] = rFT[0,1:-1][::-1].real - 1j*rFT[0,1:-1][::-1].imag
        temp[1:,:] = rFT[1:,1:-1][::-1,::-1].real - \
                1j*rFT[1:,1:-1][::-1,::-1].imag
    elif rFT.ndim == 3:

        temp[0,0,:] = rFT[0,0,1:-1][::-1].real - 1j*rFT[0,0,1:-1][::-1].imag
        temp[0,1:,:] = rFT[0,1:,1:-1][::-1,::-1].real - \
                1j*rFT[0,1:,1:-1][::-1,::-1].imag
        temp[1:,0,:] = rFT[1:,0,1:-1][::-1,::-1].real - \
                1j*rFT[1:,0,1:-1][::-1,::-1].imag
        temp[1:,1:,:] = rFT[1:,1:,1:-1][::-1,::-1,::-1].real - \
                1j*rFT[1:,1:,1:-1][::-1,::-1,::-1].imag
    else: 

        message = "\nSorry, rFT must be ndim 2 or 3!"
        raise RuntimeError, message
        
    return N.concatenate((rFT[...,:rFT.shape[-1]], temp), axis=-1)
####


'''unused seb 20091118
######  FUNCTION: fullFT2realFT  #######
####    [checked]
def fullFT2realFT(fFT):
    
    for i in range(len(fFT.shape)):
    
        if IsOdd(fFT.shape[i]): #@
        
            message = "\nfFT array must be of an even shaped array!"
            raise RuntimeError, message
    
    return fFT[...,:fFT.shape[-1]/2+1]
####


######  FUNCTION: GenImageDataSet  #######
####    [checked]
def GenMonoImageDataSet(filebase, true_object, true_PSF, directory='.',
        Nimages=1, snr=[10], snr_type='Ivar', dimension=2, format='f',
        output_dtype= N.int16):
    
    if os.path.isdir(directory):

        directory = os.path.abspath(directory) + '/'
    else:
    
        os. makedirs(os.path.abspath(directory))
        directory = os.path.abspath(directory) + '/'
        
    ## BELOW IS TEMP
    #if true_object.shape != true_PSF.shape:
    #            
    #    max, pixel_center = LocatePSFcenter(true_PSF, xyCentroidSize=None)
    #    true_PSF = ResizePSF(true_PSF, new_shape=true_object.shape, 
    #            pixel_center=None, fill=0.)
    #            
    #max, pixel_center = LocatePSFcenter(true_PSF, xyCentroidSize=None)
    #U.nd.shift(input=true_PSF, shift=pixel_center, output=true_PSF, order=3, 
    #        mode="wrap")

    shape = true_object.shape
    realFFTshape = shape[:-1] + (shape[-1]/2+1,)

    scale = 1./N.product(shape)

    ## normalize PSF
    #true_PSF /= float(N.sum(true_PSF.flat))

    OBJECT = N.empty(shape=realFFTshape, dtype=N.complex128)
    noiseless_image = N.empty(shape=shape, dtype=N.float64)
    fftw.rfft(a=true_object.astype(N.float64), af=OBJECT, inplace=False)
    OTF = PSF2OTF(true_PSF, OTF_threshold_percent=1e-5, fill=0., fullFT=False)
    fftw.irfft(af=(OBJECT*OTF*scale), a=noiseless_image, inplace=False,
            copy=False)

    Output2File(format=format, data_array=noiseless_image,
            filebase=directory + 'noiseless_' + filebase, hdr=None, shape=None)
            
    
    (Imin, Imax, Iave, Istd) = U.mmms(noiseless_image)
    std_list = []
    
    for i in range(len(snr)):

        for j in range(Nimages):
    
            # add poisson noise first
            image = N.where(noiseless_image<=0, 0, N_random_array.poisson(
                    mean=noiseless_image)).astype(output_dtype)
    
            # then add gaussian noise based on snr input
            if snr_type == 'Imax':
        
                # snr = Imax/sqrt(<w>)
                if snr[i] > Imax/N.sqrt(Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
        
                flag = 'Imax'
                std = N.sqrt((Imax/snr[i])**2 - Iave)  # commonly used in astro
                image += N_random_array.normal(mean=0., std=std,
                        shape=image.shape)
            elif snr_type == 'Ivar':

                # snr = (Istd**2)/<w>
                if snr[i] > (Istd**2)/Iave:
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
    
                flag = 'Ivar'   
                std = N.sqrt((Istd**2)/snr[i] - Iave)
                                                    # ~Bertero & Boccacci p. 53
                image += N_random_array.normal(mean=0., std=std,
                        shape=image.shape)
            elif snr_type == 'Bert':
        
                # snr = 10*log10((Istd**2)/<w>)     
                if snr[i] > 10*N.log10((Istd**2)/Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
        
                flag = 'Bert'
                std = N.sqrt((Istd**2)/(10**(snr[i]/10)) - Iave)
                                                    # Bertero & Boccacci p. 53
                image += N_random_array.normal(mean=0., std=std,
                        shape=image.shape)
            elif snr_type == 'Iave':

                # snr = Iave/sqrt(<w>)
                if snr[i] > N.sqrt(Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 

                flag = 'Iave'
                std = N.sqrt((Iave/snr[i])**2 - Iave)
                                                    # ~CCD Eq of Howell p.54
                image += N_random_array.normal(mean=0., std=std,
                        shape=image.shape)
            elif snr_type == 'Astr':

                # snr = Imax/<sigma_det>
                if snr[i] == 0:
            
                    message = "\n'snr' value of " + str(snr[i]) + " must " + \
                            "be > 0!"
                    raise RuntimeError, message 

                flag = 'Astr'
                std = (Imax/snr)                    # Practical Astro Imaging
                image += N_random_array.normal(mean=0., std=std,
                        shape=image.shape)
            
            if Nimages > 1:
            
                strg = string.split('%.2f' %(j/100.), '.')[1]
                filename = directory +  filebase + '_' + strg + '_' + flag + \
                        'snr=' + Num2Str(snr[i], digits=4) + '_Gstd=' + \
                        Num2Str(std, digits=4)
            else:
            
                filename = directory +  filebase + '_' + flag + 'snr=' + \
                        Num2Str(snr[i], digits=4) + '_Gstd=' + \
                        Num2Str(std, digits=4)              
            
            Output2File(format=format, data_array=image, filebase=filename,
                    hdr=None, shape=None)   
                    
            print filename, ': Gaussian std =', std
            std_list.append(std)
            
    print std_list
####                    


######  FUNCTION: Gen2DGaussPSFs  #######
####    [checked]
def Gen2DGaussPSFs(filebase='psf2D', directory='.', fwhm=4, peakVal=10000,
        Npsfs=20, snr=[5], snr_type='Ivar', length=1024, format='f',
        output_dtype=N.int16):

    if os.path.isdir(directory):

        directory = os.path.abspath(directory) + '/'
    else:
    
        message = "\n'directory' provided is not a valid directory path!"
        raise ValueError, message
        
    sigma = fwhm / (2*N.sqrt(2*N.log(2)))
    noiseless_PSF = F.gaussianArr(shape=(length,length), sigma=sigma, 
            integralScale=None, peakVal=peakVal, orig=None) 
    Output2File(format=format, data_array=noiseless_PSF,
            filebase=directory + 'noiseless_' + filebase + '_fwhm=' + \
            str(fwhm) + '_sh=' + str(length), hdr=None, shape=None)   

    print directory

    (Imin, Imax, Iave, Istd) = U.mmms(noiseless_PSF)

    for i in range(len(snr)):

        for j in range(Npsfs):

            # add poisson noise first   
            PSF = N.where(noiseless_PSF<=0, 0, N_random_array.poisson( 
                    mean=noiseless_PSF)).astype(output_dtype)
    
            # then add gaussian noise based on snr input
            if snr_type == 'Imax':
        
                # snr = Imax/sqrt(<w>)
                if snr[i] > Imax/N.sqrt(Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
        
                flag = 'Imax'
                std = N.sqrt((Imax/snr[i])**2 - Iave)  # commonly used in astro
                PSF += N_random_array.normal(mean=0., std=std, \
                        shape=PSF.shape)
            elif snr_type == 'Ivar':

                # snr = (Istd**2)/<w>
                if snr[i] > (Istd**2)/Iave:
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
    
                flag = 'Ivar'   
                std = N.sqrt((Istd**2)/snr[i] - Iave)
                                                    # ~Bertero & Boccacci p. 53
                PSF += N_random_array.normal(mean=0., std=std, \
                        shape=PSF.shape)
            elif snr_type == 'Bert':

                # snr = 10*log10((Istd**2)/<w>)     
                if snr[i] > 10*N.log10((Istd**2)/Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
        
                flag = 'Bert'
                std = N.sqrt((Istd**2)/(10**(snr[i]/10)) - Iave)
                                                    # Bertero & Boccacci p. 53
                PSF += N_random_array.normal(mean=0., std=std, \
                        shape=PSF.shape)
            elif snr_type == 'Iave':

                # snr = Iave/sqrt(<w>)
                if snr[i] > Iave:
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 

                flag = 'Iave'
                std = N.sqrt((Iave/snr[i])**2 - Iave)
                                                    # ~CCD Eq of Howell p.54
            elif snr_type == 'Astr':

                # snr = Imax/<sigma_det>
                if snr[i] == 0:
            
                    message = "\n'snr' value of " + str(snr[i]) + " must " + \
                            "be > 0!"
                    raise RuntimeError, message 

                flag = 'Astr'
                std = (Imax/snr)                    # Practical Astro Imaging

            if Npsfs > 1:
            
                strg = string.split('%.2f' %(j/100.), '.')[1]
                filename = directory + filebase + '_' + strg + '_fwhm=' + \
                        Num2Str(fwhm, digits=2) + '_' + flag + 'snr=' + \
                        Num2Str(snr[i], digits=4) + '_sh=' + str(length)
            else:

                filename = directory + filebase + '_fwhm=' + Num2Str(fwhm, \
                        digits=2) + '_' + flag + 'snr=' + Num2Str(snr[i], \
                        digits=4) + '_sh=' + str(length)
            
            print filename

            Output2File(format=format, data_array=PSF, filebase=filename,
                    hdr=None, shape=None)   

    print 'noiseless_PSF stats:', (Imin, Imax, Iave, Istd)
    
    return noiseless_PSF
####                    


######  FUNCTION: Gen2DGaussPSFs  #######
####    [checked]
def Gen2DAberGaussPSFs(filebase='psf2D', directory='.', fwhm=4, sig_fwhm=1,
        peakVal=10000, Npsfs=20, snr=[5], snr_type='Ivar', length=1024,
        format='f', output_dtype=N.int16):

    if os.path.isdir(directory):

        directory = os.path.abspath(directory) + '/'
    else:
    
        message = "\n'directory' provided is not a valid directory path!"
        raise ValueError, message

    rfwhm = N_random_array.normal(mean=fwhm, std=sig_fwhm)
    sigma = rfwhm / (2*N.sqrt(2*N.log(2)))
    rand_noiseless_PSF = F.gaussianArr(shape=(length,length), sigma=sigma,
            integralScale=None, peakVal=peakVal, orig=None) 
    Output2File(format=format, data_array=rand_noiseless_PSF,
            filebase=directory + 'rand_aber_noiseless_' + filebase + \
            '_fwhm=' + Num2Str(rfwhm, digits=3) + '_sh=' + str(length), 
            hdr=None, shape=None)   

    print directory

    (Imin, Imax, Iave, Istd) = U.mmms(rand_noiseless_PSF)

    for i in range(len(snr)):

        for j in range(Npsfs):

            rfwhm = N_random_array.normal(mean=fwhm, std=sig_fwhm)
            sigma = rfwhm / (2*N.sqrt(2*N.log(2)))
            noiseless_PSF = F.gaussianArr(shape=(length,length), sigma=sigma,
                    integralScale=None, peakVal=peakVal, orig=None) 
            (Imin, Imax, Iave, Istd) = U.mmms(noiseless_PSF)

            # add poisson noise first   
            PSF = N.where(noiseless_PSF<=0, 0, N_random_array.poisson(
                    mean=noiseless_PSF)).astype(output_dtype)
    
            # then add gaussian noise based on snr input
            if snr_type == 'Imax':
        
                # snr = Imax/sqrt(<w>)
                if snr[i] > Imax/N.sqrt(Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
        
                flag = 'Imax'
                std = N.sqrt((Imax/snr[i])**2 - Iave)  # commonly used in astro
                PSF += N_random_array.normal(mean=0., std=std,
                        shape=PSF.shape)
            elif snr_type == 'Ivar':

                # snr = (Istd**2)/<w>
                if snr[i] > (Istd**2)/Iave:
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
    
                flag = 'Ivar'   
                std = N.sqrt((Istd**2)/snr[i] - Iave)
                                                    # ~Bertero & Boccacci p. 53
                PSF += N_random_array.normal(mean=0., std=std,
                        shape=PSF.shape)
            elif snr_type == 'Bert':

                # snr = 10*log10((Istd**2)/<w>)     
                if snr[i] > 10*N.log10((Istd**2)/Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
        
                flag = 'Bert'
                std = N.sqrt((Istd**2)/(10**(snr[i]/10)) - Iave)
                                                    # Bertero & Boccacci p. 53
                PSF += N_random_array.normal(mean=0., std=std,
                        shape=PSF.shape)
            elif snr_type == 'Iave':

                # snr = Iave/sqrt(<w>)
                if snr[i] > Iave:
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 

                flag = 'Iave'
                std = N.sqrt((Iave/snr[i])**2 - Iave)
                                                    # ~CCD Eq of Howell p.54
            elif snr_type == 'Astr':

                # snr = Imax/<sigma_det>
                if snr[i] == 0:
            
                    message = "\n'snr' value of " + str(snr[i]) + " must " + \
                            "be > 0!"
                    raise RuntimeError, message 

                flag = 'Astr'
                std = (Imax/snr)                    # Practical Astro Imaging

            if Npsfs > 1:
            
                strg = string.split('%.2f' %(j/100.), '.')[1]
                filename = directory + filebase + '_' + strg + '_fwhm=' + \
                        Num2Str(fwhm, digits=2) + '_' + flag + 'snr=' + \
                        Num2Str(snr[i], digits=4) + '_sh=' + str(length)
            else:

                filename = directory + filebase + '_fwhm=' + Num2Str(fwhm,
                        digits=2) + '_' + flag + 'snr=' + Num2Str(snr[i],
                        digits=4) + '_sh=' + str(length)
            
            print filename

            Output2File(format=format, data_array=PSF, filebase=filename,
                    hdr=None, shape=None)   

    print 'rand_noiseless_PSF stats:', (Imin, Imax, Iave, Istd)
    
    return rand_noiseless_PSF
####                    



######  FUNCTION: GenPupil2PSF  #######
####    [checked]
def GenPupil2PSF(shape=(256,256), PSF_radius=4, zernike_mode=[0],
        zernike_coef=[1], kappa=3):

    if shape[-2] != shape[-1]:
        
        raise ValueError, "shape must be square in xy!"
        
    for i in (-2, -1):
    
        if N.remainder(shape[i], 2) != 0:
        
            raise RuntimeError, "shape in xy must be even!"
        
    if len(shape) == 3 and shape[-3] < 2:
    
        shape = shape[-2:]

    if len(zernike_mode) != len(zernike_coef):
    
        raise ValueError, "list lengths of 'zernike_mode' and " + \
                "'zernike_coef' do not match!" 

    coPSF = N.empty(shape=shape, dtype=N.complex128)
    PSF = N.zeros(shape=shape, dtype=N.float64)
    scale = 1./N.product(shape[:2])
    (ny, nx) = N.array(shape[-2:])/2

    if len(shape) == 2:
    
        pupil_radius = shape[-1]/float(2.*PSF_radius)
        pupil_aperture = F.zzernikeArr(shape=shape, no=0, crop=1,
                radius=pupil_radius, orig=None).astype(N.float64)
        phase = N.zeros(shape=shape, dtype=N.float64)
        
        for i in range(len(zernike_mode)):
            
            phase += zernike_coef[i] * F.zzernikeArr(shape=shape, 
                    no=zernike_mode[i], crop=1, radius=pupil_radius, 
                    orig=None) * 2*N.pi
    
        fftw.ifft(af=(scale*pupil_aperture*N.exp(1j*phase)), a=coPSF, 
                inplace=False)
        
        ## shift to array origin
        PSF[0:ny, 0:nx] = N.abs(coPSF[ny:, nx:])**2
        PSF[0:ny, nx:] = N.abs(coPSF[ny:, 0:nx])**2
        PSF[ny:, 0:nx] = N.abs(coPSF[0:ny, nx:])**2
        PSF[ny:, nx:] = N.abs(coPSF[0:ny, 0:nx])**2

    elif len(shape) == 3:

        (nZ, nY, nX) = shape
        (nz, ny, nx) = N.array(shape)/2
        inv_nZ = 1./nZ
        inv_nY = 1./nY
        inv_nX = 1./nX
        pupil_radius = shape[-1]/float(2.*PSF_radius)
        pupil_aperture = F.zzernikeArr(shape=(nY,nX), no=0, crop=1,
                radius=pupil_radius, orig=None).astype(N.float64)
        phase = N.zeros(shape=shape[-2:], dtype=N.float64)
        Y.view(pupil_aperture)
        
        for i in range(len(zernike_mode)):
            
            phase += zernike_coef[i] * F.zzernikeArr(shape=shape[-2:],
                    no=zernike_mode[i], crop=1, radius=pupil_radius,
                    orig=None) * 2*N.pi
        
        pupil_function = pupil_aperture*N.exp(1j*phase)
        coPSF = N.empty(shape=shape, dtype=N.complex128)
        rad2 = kappa*kappa + 2
        zextent = (shape[-3])/2
        
        def kz(ky, kx):

            return pupil_radius * N.where(ky <= nY/2, N.where(kx <= nX/2, \
                    (N.sqrt(rad2 - \
                                (kx*inv_nX)**2 - (ky*inv_nY)**2)), \
                                                            # ky, kx
                    (N.sqrt(rad2 - \
                                (kx*inv_nX - 1)**2 - (ky*inv_nY)**2)) ), \
                                                            # ky, kx-nX
                    N.where(kx <= nX/2,
                            (N.sqrt(rad2 - \
                            (kx*inv_nX)**2 - (ky*inv_nY - 1)**2)), \
                                                            # ky-nY, kx
                            (N.sqrt(rad2 - \
                            (kx*inv_nX - 1)**2 - (ky*inv_nY - 1)**2)) ) )
                                                            # ky-nY, kx-nX

        kz_zero_centered = 2j*N.pi*N.fromfunction(kz, (nY, nX))
        kz_array_centered = N.empty(shape=(nY, nX), dtype=N.complex128)
        ## assumes array in xy is even!
        kz_array_centered[0:ny, 0:nx] = kz_zero_centered[ny:, nx:]
        kz_array_centered[0:ny, nx:] = kz_zero_centered[ny:, 0:nx]
        kz_array_centered[ny:, 0:nx] = kz_zero_centered[0:ny, nx:]
        kz_array_centered[ny:, nx:] = kz_zero_centered[0:ny, 0:nx]

        zrange = N.arange(shape[-3]) - zextent

        # SHOULD REPLACE BELOW WITH ARRAY MULTIPLICATION APPROACH TO SAVE
        # TIME ON FOR LOOP?  NEED TO PROFILE AND SEE WHERE THE MAJOR COST
        # IS IN COMPUTING COPSF
        for z in zrange:
        
            zfactor = N.exp(kz_array_centered*float(z))            
            fftw.ifft(af=(scale*pupil_function*zfactor), a=coPSF[z], \
                    inplace=False)

        if nz == 0:
        
            PSF[:, 0:ny, 0:nx] = N.abs(coPSF[:, ny:, nx:])**2
            PSF[:, 0:ny, nx:] = N.abs(coPSF[:, ny:, 0:nx])**2
            PSF[:, ny:, 0:nx] = N.abs(coPSF[:, 0:ny, nx:])**2
            PSF[:, ny:, nx:] = N.abs(coPSF[:, 0:ny, 0:nx])**2
            PSF[:, 0:ny, 0:nx] = N.abs(coPSF[:, ny:, nx:])**2
            PSF[:, 0:ny, nx:] = N.abs(coPSF[:, ny:, 0:nx])**2
            PSF[:, ny:, 0:nx] = N.abs(coPSF[:, 0:ny, nx:])**2
            PSF[:, ny:, nx:] = N.abs(coPSF[:, 0:ny, 0:nx])**2
        else:
        
            PSF[0:nz, 0:ny, 0:nx] = N.abs(coPSF[nz:, ny:, nx:])**2
            PSF[0:nz, 0:ny, nx:] = N.abs(coPSF[nz:, ny:, 0:nx])**2
            PSF[0:nz, ny:, 0:nx] = N.abs(coPSF[nz:, 0:ny, nx:])**2
            PSF[0:nz, ny:, nx:] = N.abs(coPSF[nz:, 0:ny, 0:nx])**2
            PSF[nz:, 0:ny, 0:nx] = N.abs(coPSF[0:nz, ny:, nx:])**2
            PSF[nz:, 0:ny, nx:] = N.abs(coPSF[0:nz, ny:, 0:nx])**2
            PSF[nz:, ny:, 0:nx] = N.abs(coPSF[0:nz, 0:ny, nx:])**2
            PSF[nz:, ny:, nx:] = N.abs(coPSF[0:nz, 0:ny, 0:nx])**2

    else:
    
        raise ValueError, "shape must be 2D or 3D!"

    return PSF/N.sum(PSF.flat)
####


######  FUNCTION: GenAverZernPSFs  #######
####    [checked]
def GenAberZernPSFs(directory=None, shape=(256,256), Npsfs=20, PSF_radius=4,
        aber_zern_modes=[1,2,3], zern_coef_mean=[0,0,0],
        zern_coef_sigma=[0.2,0.2,0.2], Imax=1000, snr=[50], 
        snr_type='Ivar', format='f', kappa=3, output_dtype=N.int16):

    if directory:
    
        if os.path.isdir(directory):

            directory = os.path.abspath(directory) + '/'
        else:
    
            message = "\n'directory' provided is not a valid directory path!"
            raise ValueError, message

    total_zern_modes = len(aber_zern_modes)
    
    if total_zern_modes != len(zern_coef_mean):
    
        message = "\n'zern_coef_mean' does not match length of 'aber_zern_modes'!"
        raise ValueError, message

    elif total_zern_modes != len(zern_coef_sigma):
    
        message = "\n'zern_coef_sigma' does not match length of 'aber_zern_modes'!"
        raise ValueError, message

    noiseless_PSF = N.empty(shape=(Npsfs,)+shape, dtype=N.float32)
    PSF = N.empty(shape=(Npsfs,)+shape, dtype=N.float32)
    
    for i in range(Npsfs):
    
#       noiseless_PSF[i] = GenPupil2PSF(shape=shape, \
#               PSF_radius=PSF_radius, \
#               zernike_mode=tuple(N.arange(0,total_zern_modes+1)), \
#               zernike_coef=tuple(N_random_array.normal(0, \
#               std=zern_coef_sigma, shape=[total_zern_modes + 1])) )

        coefficients = tuple(N_random_array.normal(mean=zern_coef_mean,
                std=zern_coef_sigma, shape=(total_zern_modes,)))
        print i, coefficients
        noiseless_PSF[i] = GenPupil2PSF(shape=shape,
                PSF_radius=PSF_radius, 
                zernike_mode=tuple(aber_zern_modes),
                zernike_coef=coefficients, kappa=kappa )
        noiseless_PSF[i] *= Imax/noiseless_PSF[i].max()
                
#   chosen_PSF = noiseless_PSF[0]

#   Output2File(format=format, data_array=chosen_PSF, filebase=directory + \
#           'noiseless_rand_aberPSF_Z' + str(total_zern_modes) + '_coefsig=' + \
#           Num2Str(zern_coef_sigma, digits=3) + '_rad=' + \
#           Num2Str(PSF_radius, digits=3) + '_sh=' + str(shape[-1]), hdr=None,
#           shape=None)   
#   Output2File(format=format, data_array=chosen_PSF, filebase=directory + \
#           'chosen_noiseless_rand_aberPSF_Z' + Num2Str(aber_zern_modes) + 
#           '_coefave=' + Num2Str(zern_coef_mean) + '_coefsig=' + \
#           Num2Str(zern_coef_sigma, digits=3) + '_rad=' + \
#           Num2Str(PSF_radius, digits=3) + '_sh=' + str(shape[-1]), hdr=None,
#           shape=None)   

    if directory:
    
        Output2File(format=format, data_array=noiseless_PSF, filebase=directory + \
                'noiseless_rand_aberPSF_Z' + Num2Str(aber_zern_modes) + 
                '_coefave=' + Num2Str(zern_coef_mean) + '_coefsig=' + \
                Num2Str(zern_coef_sigma, digits=3) + '_rad=' + \
                Num2Str(PSF_radius, digits=3) + '_sh=' + str(shape[-1]), hdr=None,
                shape=None)   

    for i in range(len(snr)):

        for j in range(Npsfs):

            #Y.view(noiseless_PSF[j])

            (Imin, Imax, Iave, Istd) = U.mmms(noiseless_PSF[j])
            #print (Imin, Imax, Iave, Istd)
            
            #median = N.sort(noiseless_PSF[j].flat)[len(noiseless_PSF[j].flat)/2]
            #test1 = N.compress(noiseless_PSF[j] > median, noiseless_PSF[j])
            #test2 = N.compress(noiseless_PSF[j] > Iave, noiseless_PSF[j])
            #print U.mmms(test1)
            #print U.mmms(test2)

            # add poisson noise first   
            PSF[j] = N_random_array.poisson(mean=noiseless_PSF[j]).astype(output_dtype)
    
            # then add gaussian noise based on snr input
            if snr_type == 'Imax':
        
                # snr = Imax/sqrt(<w>)
                if snr[i] > Imax/N.sqrt(Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
        
                flag = 'Imax'
                std = N.sqrt((Imax/snr[i])**2 - Iave)  # commonly used in astro
                PSF[j] += N_random_array.normal(mean=0., std=std,
                        shape=PSF[j].shape)
            elif snr_type == 'Ivar':

                # snr = (Istd**2)/<w>
                if snr[i] > (Istd**2)/Iave:
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
    
                flag = 'Ivar'   
                std = N.sqrt((Istd**2)/snr[i] - Iave)
                                                    # ~Bertero & Boccacci p. 53
                PSF[j] += N_random_array.normal(mean=0., std=std, \
                        shape=PSF[j].shape)
            elif snr_type == 'Bert':

                # snr = 10*log10((Istd**2)/<w>)     
                if snr[i] > 10*N.log10((Istd**2)/Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
        
                flag = 'Bert'
                std = N.sqrt((Istd**2)/(10**(snr[i]/10)) - Iave)
                                                    # Bertero & Boccacci p. 53
                PSF[j] += N_random_array.normal(mean=0., std=std,
                        shape=PSF[j].shape)
            elif snr_type == 'Iave':

                # snr = Iave/sqrt(<w>)
                if snr[i] > N.sqrt(Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 

                flag = 'Iave'
                std = N.sqrt((Iave/snr[i])**2 - Iave)
                                                    # ~CCD Eq of Howell p.54
                PSF[j] += N_random_array.normal(mean=0., std=std,
                        shape=PSF[j].shape)
            elif snr_type == 'Astr':

                # snr = Imax/<sigma_det>
                if snr[i] == 0:
            
                    message = "\n'snr' value of " + str(snr[i]) + " must " + \
                            "be > 0!"
                    raise RuntimeError, message 

                flag = 'Astr'
                std = (Imax/snr)                    # Practical Astro Imaging
                PSF[j] += N_random_array.normal(mean=0., std=std,
                        shape=PSF[j].shape)
                        
        if directory:
        
            if Npsfs > 1:
                
                #strg = string.split('%.2f' %(j/100.), '.')[1]
                strg = 'n' + str(Npsfs)
#               filename = directory + 'rand_aberPSF_Z' + \
#                       str(total_zern_modes) + '_coefsig=' + \
#                       Num2Str(zern_coef_sigma, digits=3) + '_' + strg + \
#                       '_rad=' + Num2Str(PSF_radius, digits=3) + '_' + flag + \
#                       'snr=' + Num2Str(snr[i], digits=4) + '_sh=' + \
#                       str(shape[-1])
                filename = directory + 'rand_aberPSF_Z' + \
                        Num2Str(aber_zern_modes) + \
                        '_coefave=' + Num2Str(zern_coef_mean) +'_coefsig=' + \
                        Num2Str(zern_coef_sigma, digits=3) + '_' + strg + \
                        '_rad=' + Num2Str(PSF_radius, digits=3) + '_' + flag + \
                        'snr=' + Num2Str(snr[i], digits=4) + '_sh=' + \
                        str(shape[-1])
            else:

#               filename = directory + 'rand_aberPSF_Z' + \
#                       str(total_zern_modes) + '_coefsig=' + \
#                       Num2Str(zern_coef_sigma, digits=3) + '_rad=' + \
#                       Num2Str(PSF_radius, digits=3) + '_' + flag + 'snr=' + \
#                       Num2Str(snr[i], digits=4) + '_sh=' + str(shape[-1])
                filename = directory + 'rand_aberPSF_Z' + \
                        Num2Str(aber_zern_modes) + \
                        '_coefave=' + Num2Str(zern_coef_mean) +'_coefsig=' + \
                        Num2Str(zern_coef_sigma, digits=3) + '_rad=' + \
                        Num2Str(PSF_radius, digits=3) + '_' + flag + 'snr=' + \
                        Num2Str(snr[i], digits=4) + '_sh=' + str(shape[-1])
            
            print filename

            Output2File(format=format, data_array=PSF, filebase=filename, hdr=None,
                    shape=None)   
    
#   return (chosen_PSF, noiseless_PSF, PSF)
    return (noiseless_PSF, PSF)
####                    


def PlotLogRadialFourierProfile(array, new=1, normalize_flag=1, c="-"):

    if array.ndim == 3:
    
        zstacks = array.shape[-3]
    else:
    
        array.shape = (1,) + array.shape
        zstacks = 1

    if new == 1:

        Y.plotFigure()
        
    for i in range(zstacks):

        if normalize_flag:

            r1d = RadiallyAverage2D(array=N.log10(N.abs( 
                    F.fft2d(array[i]/N.sum(array[i].flat)))), FT=True, 
                    origin=None, wrap=(False,), subtract_background=False,
                    fill=0.)
        else:
        
            r1d = RadiallyAverage2D(array=N.log10(N.abs( 
                    F.fft2d(array[i]/N.sum(array[i].flat)))), FT=True, 
                    origin=None, wrap=(False,), subtract_background=False,
                    fill=0.)[0] #@

        if i == 0:
        
            Y.ploty(r1d, c=c)
            Y.plothold()
        else:
        
            Y.ploty(r1d, c=c)
        
    #return r1d
            

def PlotRadialFourierProfile(array, normalize_flag=1, c="-"):

    if normalize_flag:
        r1d = RadiallyAverageAboutOrigin(array=(N.abs( 
                F.fft2d(array/N.sum(array.flat)))), dimensions_to_average=2, 
                floor=0.0, fullFT_flag=1)[0]
    else:
        r1d = RadiallyAverageAboutOrigin(array=(N.abs( 
                F.fft2d(array))), dimensions_to_average=2, 
                floor=0.0, fullFT_flag=1)[0]
    Y.ploty(r1d, c=c)
    return r1d

### TEMP FUNCTION - needs generalization
def PlotRadialFourierProfileDiff(ref_array, dirprefix, var_list_order=[], 
        filebase='Rpsf', normalize_flag=1, c='-'):

    import glob
    
    names = os.listdir(os.getcwd())
    alldir = []
    
    for k in names:
    
        if k[:len(dirprefix)] == dirprefix:
        
            alldir.append(k)

    commonprefix = os.path.commonprefix(alldir)
    dirlist = []
    test = 0
    Y.plotFigure()

    print commonprefix
    
    if normalize_flag:

        true_r1d = RadiallyAverage2D(array=(N.abs( 
                F.fft2d(ref_array/N.sum(ref_array.flat)))), FT=True, 
                origin=None, wrap=(False,), subtract_background=False,
                fill=0.)[0]

    else:

        true_r1d = RadiallyAverage2D(array=(N.abs(F.fft2d(ref_array))), 
                FT=True, origin=None, wrap=(False,), subtract_background=False,
                fill=0.)[0]

    for i in var_list_order:
    
        dirname = os.getcwd() + '/' + commonprefix + str(i) + 'p*'
        dirlist.append(glob.glob(dirname)[0])
        
    for j in dirlist:
    
        allfiles = os.listdir(j)
        file = glob.glob(j + '/' + filebase + '*')[0]
        data = LoadFile(file)[0]
        
        array = data
        #array = N.abs(reference_array - data)
        
        if normalize_flag:

            r1d = RadiallyAverage2D(array=(N.abs( 
                    F.fft2d(array/N.sum(array.flat)))), FT=True, 
                    origin=None, wrap=(False,), subtract_background=False,
                    fill=0.)[0]                 
            diff_r1d = N.abs(r1d-true_r1d)

        else:

            r1d = RadiallyAverage2D(array=(N.abs(F.fft2d(array))), FT=True, 
                    origin=None, wrap=(False,), subtract_background=False,
                    fill=0.)[0]
            diff_r1d = N.abs(r1d-true_r1d)

        if test == 0:
        
            Y.ploty(diff_r1d, c=c)
            Y.plothold()
            test = 1
        else:
            Y.ploty(diff_r1d, c=c)


### TEMP FUNCTION - needs generalization
def PlotRadialFourierProfileDiff2(ref_array, array, normalize_flag=1, c='-'):
    
    if normalize_flag:

        true_r1d = RadiallyAverage2D(array=(N.abs( 
                F.fft2d(ref_array/N.sum(ref_array.flat)))), FT=True, 
                origin=None, wrap=(False,), subtract_background=False,
                fill=0.)[0]
        r1d = RadiallyAverage2D(array=(N.abs( 
                F.fft2d(array/N.sum(array.flat)))), FT=True, 
                origin=None, wrap=(False,), subtract_background=False,
                fill=0.)[0]                 
        diff_r1d = N.abs(r1d-true_r1d)

    else:

        true_r1d = RadiallyAverage2D(array=(N.abs(F.fft2d(ref_array))), 
                FT=True, origin=None, wrap=(False,), subtract_background=False,
                fill=0.)[0]
        r1d = RadiallyAverage2D(array=(N.abs(F.fft2d(array))), FT=True, 
                origin=None, wrap=(False,), subtract_background=False,
                fill=0.)[0]
        diff_r1d = N.abs(r1d-true_r1d)

    Y.ploty(diff_r1d, c=c)
#####


######  FUNCTION: GenAberZernPSFs  #######
####    [checked]
def Gen2DharmonicOTFNoisyPSFs(originPSF, stdOTF, directory=None, nn=20, Imax=1000, 
        snr=[50], snr_type='Ivar', format='f', output_dtype=N.int16):

    if directory:
    
        if os.path.isdir(directory):

            directory = os.path.abspath(directory) + '/'
        else:
    
            message = "\n'directory' provided is not a valid directory path!"
            raise ValueError, message

    shape = originPSF.shape
    #real_shape = shape[:-1] + (shape[-1]/2+1,)
    PSFs = N.empty(shape=(nn,)+shape, dtype=N.float32)
    cleanPSFs = N.empty(shape=(nn,)+shape, dtype=N.float32)
    #OTFs = N.empty(shape=(nn,)+real_shape, dtype=N.complex64)
    OTFs = N.empty(shape=(nn,)+shape, dtype=N.complex64)

    nPSF = originPSF/N.sum(originPSF.flat)
    scaling = Imax/nPSF.max()
    PSF = nPSF*scaling
    
    OTF = F.fft2d(nPSF)
    phase = U.phase(OTF)
    
#   if directory:
#   
#       Output2File(format=format, data_array=PSF, filebase=directory + \
#               'noiseless_PSF' + '_sh=' + str(shape[-1]), hdr=None,
#               shape=None)   

    for i in range(len(snr)):

        for j in range(nn):

            real = U.nd.median_filter(input=N_random_array.normal(mean=0., 
                    std=stdOTF*N.cos(phase), shape=OTF.shape), size=3, 
                    footprint=None, output=None, mode="reflect", cval=0.0, 
                    origin=0)
            imag = U.nd.median_filter(input=N_random_array.normal(mean=0., 
                    std=stdOTF*N.sin(phase), shape=OTF.shape), size=3, 
                    footprint=None, output=None, mode="reflect", cval=0.0, 
                    origin=0)
        
#           OTFs[j] = U.nd.median_filter(input=(OTF.real + real), size=7, 
#                   footprint=None, output=None, mode="wrap", cval=0.0,
#                   origin=0)) + \
#                   1j * U.nd.median_filter(input=(OTF.imag + imag), size=7, 
#                   footprint=None, output=None, mode="wrap", cval=0.0,
#                   origin=0)) 

            OTFs[j] = U.nd.gaussian_filter(input=(OTF.real + real), sigma=2, 
                    mode="wrap") + \
                    1j * U.nd.gaussian_filter(input=(OTF.imag + imag), sigma=2, 
                    mode="wrap") 

            PSFs[j] = F.ifft2d(OTFs[j])
            cleanPSFs[j] = N.where(PSFs[j]>0, PSFs[j], 0)
            (Imin, Imax, Iave, Istd) = U.mmms(scaling*cleanPSFs[j])

            # add poisson noise first
            PSFs[j] = N_random_array.poisson(
                    mean=scaling*cleanPSFs[j]).astype(output_dtype)
    
            # then add gaussian noise based on snr input
            if snr_type == 'Imax':
        
                # snr = Imax/sqrt(<w>)
                if snr[i] > Imax/N.sqrt(Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
        
                flag = 'Imax'
                std = N.sqrt((Imax/snr[i])**2 - Iave)  # commonly used in astro
                PSFs[j] += N_random_array.normal(mean=0., std=std,
                        shape=PSFs[j].shape)
            elif snr_type == 'Ivar':

                # snr = (Istd**2)/<w>
                if snr[i] > (Istd**2)/Iave:
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
    
                flag = 'Ivar'   
                std = N.sqrt((Istd**2)/snr[i] - Iave)
                                                    # ~Bertero & Boccacci p. 53
                PSFs[j] += N_random_array.normal(mean=0., std=std,
                        shape=PSFs[j].shape)
            elif snr_type == 'Bert':

                # snr = 10*log10((Istd**2)/<w>)     
                if snr[i] > 10*N.log10((Istd**2)/Iave):
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 
        
                flag = 'Bert'
                std = N.sqrt((Istd**2)/(10**(snr[i]/10)) - Iave)
                                                    # Bertero & Boccacci p. 53
                PSFs[j] += N_random_array.normal(mean=0., std=std,
                        shape=PSFs[j].shape)
            elif snr_type == 'Iave':

                # snr = Iave/sqrt(<w>)
                if snr[i] > Iave:
            
                    message = "\n'snr' value of " + str(snr[i]) + " is too " + \
                            "high leads to an imaginary standard deviation " + \
                            "for Gaussian noise!"
                    raise RuntimeError, message 

                flag = 'Iave'
                std = N.sqrt((Iave/snr[i])**2 - Iave)
                                                    # ~CCD Eq of Howell p.54
                PSFs[j] += N_random_array.normal(mean=0., std=std,
                        shape=PSFs[j].shape)
            elif snr_type == 'Astr':

                # snr = Imax/<sigma_det>
                if snr[i] == 0:
            
                    message = "\n'snr' value of " + str(snr[i]) + " must " + \
                            "be > 0!"
                    raise RuntimeError, message 

                flag = 'Astr'
                std = (Imax/snr)                    # Practical Astro Imaging
                PSFs[j] += N_random_array.normal(mean=0., std=std,
                        shape=PSFs[j].shape)
                        
        if directory:

            for i in range(len(snr)):

                for j in range(nn):
                
                    if nn > 1:
                
                        if len(str(j)) == 1:
                        
                            strg = '_n0' + str(j) +'_'
                        else:
                        
                            strg = '_n' + str(j) + '_'
                            
                        filename = directory + 'PSF' + strg + flag + 'snr=' + \
                                Num2Str(snr[i], digits=4) + '_sh=' + \
                                str(shape[-1])
                    else:

                        filename = directory + 'PSF' + flag + 'snr=' + \
                                Num2Str(snr[i], digits=4) + '_sh=' + \
                                str(shape[-1])
            
                    print filename

                    Output2File(format=format, data_array=ArrayCenter2Origin(
                            PSFs[j]), filebase=filename, hdr=None, shape=None)   
    
#   return (chosen_PSF, noiseless_PSF, PSF)
    return (PSF, PSFs, cleanPSFs, OTFs)
####                    



def sim_v(shape=(128,128), no=[1,2,3], radius=30):

    sum=0.
    
    for i in no:
    
        sum += F.zzernikeArr(shape=shape, no=i, crop=1, radius=radius)
    
    sum=N.abs(sum)
    apod = 1 - F.zzernikeArr(shape=shape, no=3, crop=0, radius=radius-2)
    apod /= apod.max()
    apod = N.where(apod > 0, apod, 0)
    
    out = ArrayCenter2Origin(apod*sum)
    Y.view(out)
    return out
####

'''#unused seb 20091118

def U__topPercentile(a, percentile):
    """
    workaround for broken Priithon's U.topPercentile():
    always scale a to max 65535 - ignore values < 0
    """
    #_lastErrSettings = N.seterr(under="ignore")
    maa = (1<<16) - 1

    lower = a.min()
    upper = a.max()
    factor = float(maa)/(float(upper)-lower)
    th = U.topPercentile((U.asFloat32(a)-lower)*factor, percentile)/factor + lower
    #N.seterr(**_lastErrSettings)
    return th
