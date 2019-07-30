################################################################################
#
#   File:       useful.py
#
#   Summary:    Contains useful shortcut functions, from Priithon
#
#   Authors:    Clement Chalumeau (SETI, interning for Franck Marchis) 
#
#        
################################################################################

# Clement: not needed anymore
# from Priithon_heritage_bin import seb as S

try:
    from scipy import ndimage as nd
except ImportError:
    pass

import numpy as N


def max(arr):
    arr = N.asarray(arr)
    return arr.max()

def shift2D_workaround(A, shiftTuple):
    """
    A : array to be shifted
    
    """
    return nd.geometric_transform(A,lambda coords:((coords[0]-shiftTuple[0])%A.shape[0],(coords[1]-shiftTuple[1])%A.shape[1]),mode='wrap')


def _getGoodifiedArray(arr):
    """
    return "well behaved" version of a numpy array
    1) convert lists or tuple to numpy-array
    2) make copy of numpy arrays if non-contigous or non-native

    (used in conjunction with SWIGed functions)
    """
    try:
        if arr.dtype.isnative:
            arr = N.ascontiguousarray(arr)
        else:
            arr = N.ascontiguousarray(arr, arr.dtype.newbyteorder('='))
    except AttributeError:
            arr = N.ascontiguousarray(arr)

    if arr.dtype == N.bool:  # no SWIGed function for bool, use uint8
        arr = arr.view(N.uint8)

    return arr

def mm(arr):
    """
    returns min,max of arr
    """

    arr = N.asarray(arr)
    return (N.minimum.reduce(arr.flat), N.maximum.reduce(arr.flat))


def mmm(arr):
    """
    returns min,max,mean of arr
    """
    arr = _getGoodifiedArray(arr)
    #TODO: make nice for memmap
    # Clement: Changed here
    #m = S.mean(arr)
    m = N.mean(arr)
    return (N.minimum.reduce(arr.flat), N.maximum.reduce(arr.flat), m)

# Clement: Old version
# def mmms(arr):
#     """
#     returns min,max,mean,stddev of arr
#     """
#     arr = _getGoodifiedArray(arr)
#     #TODO: make nice for memmap
#     mi,ma,me,st = S.mmms( arr )
#     return (mi,ma,me,st)

def mmms(A):
    min = A.max()
    max = A.min()
    mean = A.mean()
    stdd = N.std(A)
    return (min,max,mean,stdd)

# Clement: Old version
# def findMax(arr):
#     """returns value and position of maximum in arr
#     assumes 3D array, so returned is a 4-tuple: [val,z,y,x]
#     for 2D or 1D z,y would be respectively 0
#     """
#     arr = _getGoodifiedArray(arr)
#     return S.findMax( arr )

def findMax(A):
    # Return (val,z,y,x) tuple
    imax = N.unravel_index(A.argmax(), A.shape)
    if len(imax)==1:
        return[A.max(),0,0,imax[0]]
    elif len(imax)==2:
        return[A.max(),0,imax[0],imax[1]]
    else:
        return[A.max(),imax[0],imax[1],imax[2]]
        
# Clement: Old version
# def findMin(arr):
#     """returns value and position of minimum in arr
#     assumes 3D array, so returned is a 4-tuple: [val,z,y,x]
#     for 2D or 1D z,y would be respectively 0
#     """
#     arr = _getGoodifiedArray(arr)
#     return S.findMin( arr )

def findMin(A):
    # Return (val,z,y,x) tuple
    imin = N.unravel_index(A.argmin(), A.shape)
    if len(imin)==1:
        return[A.min(),0,0,imin[0]]
    elif len(imax)==2:
        return[A.min(),0,imin[0],imin[1]]
    else:
        return[A.min(),imin[0],imin[1],imin[2]]

def max2d(arr, outtype=None):
    """returns an array of shape arr.shape[:-2] and dtype outtype
    if outtype=None it uses arr.dtype
    """
    if outtype is None:
        outtype = arr.dtype

    b = N.empty(shape=arr.shape[:-2], dtype=outtype)
    bb = b.view()
    bb.shape = (-1,)
    aarr = arr.view()
    aarr.shape = (-1,) + arr.shape[-2:]

    for i in range( bb.shape[0] ):
        bb[i] = aarr[i].max()

    return b

def fitAny(f, parmTuple0, data):
    '''
    data should be list of (x,y)  tuples
    TODO: or (x,y,deltaY)
    (instead of 'list' you can of course have an array w/
    shape=(n,2) or shape=(n,3), n beeing the number of data points

    if data.ndim == 1 or data.shape = (n,1) it fits assuming x=0,1,2,3,...n-1

    f is your 'model' function that takes two arguments:
    a tuple of parameters and x
    
    The function returns a list containing the optimal parameter values
    and the chi-squared value describing the quality of the fit.
    '''
    from scipy.optimize import leastsq

    data = N.asarray(data, dtype=N.float64)

    if len(data.shape) == 1:
        data = N.transpose(N.array([N.arange(len(data)), data]))
    elif data.shape[1] == 1:
        data = N.transpose(N.array([N.arange(len(data)), data][0]))

    x,y = data.T
    def func(p):
        return f(p, x)-y
    
    x0 = parmTuple0
    return leastsq(func, x0)#, args=(), Dfun=None,
                   #full_output=0, col_deriv=0,
                   #ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None)
 

def yGaussian(parms=(10,100), t=0):
    '''
    t can be a scalar or a vector
    returns y value(s) of a 1D-gaussian model

    parms can be tuple of ltength 2,3 or 4, with
    2: tuple is [sigma, peakVal]
    3: tuple is [x0, sigma, peakVal]
    4: tuple is [y0, x0, sigma, peakVal]

    x0 is center of gaussian (default 0)
    y0 is baseline offset gaussian (default 0)
    sigma is sigma (stddev) of gaussian
    peakval is  "center height" above baseline
    '''
    import fftfuncs as F

    if len(parms) == 4:
        y0,x0 = parms[:2]
    elif len(parms) == 3:
        y0,x0 = 0.0, parms[0]
    else:
        y0,x0 = 0.0, 0.0
    sigma, peakVal = parms[-2:]

    return y0+F.gaussian(t-x0, dim=1, sigma=sigma, peakVal=peakVal)
    
    
def fitGaussian(data, p=(0,10,100)):
    '''
    see yGaussian.
    p: initial guess
    data: vector of data points to be fit
    '''
    return fitAny(yGaussian, p, data)
    

def yPoly(parms=(1,1,0), t=0):
    '''
    t can be a scalar or a vector
    returns y value(s) of a polygon model
    parms:
      baseline, first-order coeff, 2nd, ...
    '''
    r = 0.0
    for i in range(len(parms)):
        r = r + parms[i]*N.power(t, i)
    return r
    
def fitPoly(data, p=(1,1,1)):
    """
    see yPoly

    data should be list of y or (x,y)- or (x,y,deltaY)-tuples
    (instead of 'list' you can of course have an array w/
    shape=(n,2) or shape=(n,3), n beeing the number of data points

    uses polynomial 'model' ( U.yPoly )
    
    The function returns a list containing the optimal parameter values
    and the chi-squared value describing the quality of the fit.
    """


    return fitAny(yPoly, p, data)


# Clement Old version
# def histogram(a, nBins=None, amin=None,amax=None, histArr=None, norm=False, returnTuple=False):
#     """
#     creates/returns  array with nBins int32 entries
#        fills it with histogram of 'a'
#     if amin and/or amax is None it calculates the min/max of a and uses that
#     if nBins is None:
#         nBins = int(amax-amin+1)
#         if narr is of float dtype  Bins < 100:
#             nBins = 100
#     if histArr is given it is used to fill in the histogram values
#         then nBins must be None and histArr of dtype N.int32
# 
#     if norm:
#        divide bins (=histArr) by sum of bins and convert to float64
#     if returnTuple:
#         return (histArr, nBins, amin, amax)
#     """
#     a = N.asarray(a)
#     
#     if amin is None and amax is None:
#         amin = a.min()  
#         amax = a.max()
#     elif amin is None:
#         amin = a.min()
#     elif amax is None:
#         amax = a.max()
# 
#     if histArr is not None:
#         if nBins is not None:
#             raise "only one of histArr and nBins can be given"
#         if histArr.dtype != N.int32:
#             raise "histArr must of dtype N.int32"
#         if not histArr.flags.carray or  not histArr.dtype.isnative:
#             raise RuntimeError, 'histArr must be a "native c(ordered)-array"'
#         nBins = len(histArr)
#     else:
#         if nBins is None:
#             nBins = int(amax-amin+1)
#             if N.issubdtype(float, a.dtype) and nBins < 100:
#                 nBins = 100
# 
#         histArr = N.empty( shape=(nBins,), dtype=N.int32 )
# 
#     a = _getGoodifiedArray(a)
# 
#     # NOTE: S.histogram *ignores* all values outside range (it does not count amax !!)
#     #       it only count amin<= val < amax
#     
#     amaxTweaked = amin+nBins*(amax-amin)/(nBins-1.)
#     # CHECK numpy - why type(a.min())=numpy.float32 not SWIG compatible to float!
#     S.histogram(a, float(amin),float(amaxTweaked), histArr)
# 
#     if norm:
#         histArrNormed = N.empty( shape=(nBins,), dtype=N.float64 )
#         histArrNormed[:] = histArr
#         histArrNormed /= histArr.sum()
#         histArr = histArrNormed
# 
#     if returnTuple:
#         return (histArr, nBins, amin, amax)
#     else:
#         return histArr

def histogram(a, nBins=None, amin=None,amax=None, histArr=None, norm=False, returnTuple=False):
    """
    creates/returns  array with nBins int32 entries
       fills it with histogram of 'a'
    if amin and/or amax is None it calculates the min/max of a and uses that
    if nBins is None:
        nBins = int(amax-amin+1)
        if narr is of float dtype  Bins < 100:
            nBins = 100
    if histArr is given it is used to fill in the histogram values
        then nBins must be None and histArr of dtype N.int32

    if norm:
       divide bins (=histArr) by sum of bins and convert to float64
    if returnTuple:
        return (histArr, nBins, amin, amax)
    """
    a = N.asarray(a)
    
    if amin is None and amax is None:
        amin = a.min()  
        amax = a.max()
    elif amin is None:
        amin = a.min()
    elif amax is None:
        amax = a.max()

    if histArr is not None:
        if nBins is not None:
            raise "only one of histArr and nBins can be given"
        if histArr.dtype != N.int32:
            raise "histArr must of dtype N.int32"
        if not histArr.flags.carray or  not histArr.dtype.isnative:
            raise RuntimeError, 'histArr must be a "native c(ordered)-array"'
        nBins = len(histArr)
    else:
        if nBins is None:
            nBins = int(amax-amin+1)
            if N.issubdtype(float, a.dtype) and nBins < 100:
                nBins = 100

        histArr = N.empty( shape=(nBins,), dtype=N.int32 )

    a = _getGoodifiedArray(a)

    # NOTE: S.histogram *ignores* all values outside range (it does not count amax !!)
    #       it only count amin<= val < amax
    
    amaxTweaked = amin+nBins*(amax-amin)/(nBins-1.)
    # CHECK numpy - why type(a.min())=numpy.float32 not SWIG compatible to float!
    # Old S.histogram(a, float(amin),float(amaxTweaked), histArr)
    histArray_temp, bins = N.histogram(a,bins = nBins, range = (float(amin),float(amaxTweaked)))
    histArr[:] = histArray_temp
    if norm:
        histArrNormed = N.empty( shape=(nBins,), dtype=N.float64 )
        histArrNormed[:] = histArr
        histArrNormed /= histArr.sum()
        histArr = histArrNormed

    if returnTuple:
        return (histArr, nBins, amin, amax)
    else:
        return histArr

# Clement: Old version
# def generalhistogram(a, weightImg, nBins=None, amin=None,amax=None):
#     """
#     creates/returns ("histogram") array with nBins entries of same dtype as weightImg
#     while for a standard histogram one adds up 1s in bins for
#           each time you encouter a certain value in a
#     generalhistogram  adds the pixel value found in weightImg 
#           each time it encouters a certain value in a (for that pixel)
#     
#     if amin and/or amax is None it calculates the min/max of a and uses that
#     if nBins is None:
#         nBins = int(amax-amin+1)
#         if a is of float dtype   and nBins < 100:
#              nBins = 100
#     """
#     if amin is None and amax is None:
#         amin = a.min()
#         amax = a.max()
#     elif amin is None:
#         amin = a.min()
#     elif amax is None:
#         amax = a.max()
# 
#     if nBins is None:
#         nBins = int(amax-amin+1)
#         if N.issubdtype(float, a.dtype) and nBins < 100:
#             nBins = 100
#     b = N.empty( shape=(nBins,), dtype=weightImg.dtype )
# 
#     a = _getGoodifiedArray(a)
#     weightImg = _getGoodifiedArray(weightImg)
# 
#     # NOTE: S.histogram *ignores* all values outside range (it does not count amax !!)
#     #       it only count amin<= val < amax
#     
#     amaxTweaked = amin+nBins*(amax-amin)/(nBins-1)
#     # CHECK numpy - why type(a.min())=numpy.float32 not SWIG compatible to float!
#     S.generalhist(a, weightImg, float(amin),float(amaxTweaked), b)
# 
#     return b

def generalhistogram(a, weightImg, nBins=None, amin=None,amax=None):
    
   if amin is None and amax is None:
       amin = a.min()
       amax = a.max()
   elif amin is None:
       amin = a.min()
   elif amax is None:
       amax = a.max()
 
   if nBins is None:
       nBins = int(amax-amin+1)
       if N.issubdtype(float, a.dtype) and nBins < 100:
           nBins = 100
   # NOTE: S.histogram *ignores* all values outside range (it does not count amax !!)
   #       it only count amin<= val < amax
     
   amaxTweaked = amin+nBins*(amax-amin)/(nBins-1)
   # CHECK numpy - why type(a.min())=numpy.float32 not SWIG compatible to float!
   hist, bin = N.histogram(a,range =(float(amin),float(amaxTweaked)),weights=weightImg, bins=nBins)
   return hist  
   


def array2image(a, rgbOrder="rgba"):
    """Convert numpy array to image
       a must be of ndim 2 and dtype UInt8,Float32 or UInt16
       if a.ndim ==3:
          a.dtype must be uint8
          the first axis is interpreted as RGB color axis -- for fewer "sections" in a, remaining are assumed to be zero
          rgbOrder: order in which axes are mapped to RGB(A) channels
    """
    import Image

    # translate string to "inverse map" aAxis->rgbaIndex
    # e.g. "br" -> [1, 2, 0, 3]
    rgbOrder = rgbOrder.lower()
    rgbOrder = [rgbOrder.find(col) for col in "rgba"]
    fillAx=max(rgbOrder)+1
    for i,ax in enumerate(rgbOrder):
        if ax<0:
            rgbOrder[i] = fillAx
            fillAx+=1

    if a.ndim == 3:
        if   a.shape[0] == 1:
            assert a.dtype==N.uint8
            a22 = N.transpose(a,(1,2,0)) # .copy()
            import fftfuncs as F
            a22 = N.append(a22,F.zeroArr(a22.dtype,a22.shape[:2]+(2,)), -1)
            a22 = a22[:,:,rgbOrder[:3]]
            ii = Image.fromstring("RGB", (a.shape[-1],a.shape[-2]), a22.tostring())
            return ii
        elif   a.shape[0] == 2:
            assert a.dtype==N.uint8
            a22 = N.transpose(a,(1,2,0)) # .copy()
            import fftfuncs as F
            a22 = N.append(a22,F.zeroArr(a22.dtype,a22.shape[:2]+(1,)), -1)
            a22 = a22[:,:,rgbOrder[:3]]
            ii = Image.fromstring("RGB", (a.shape[-1],a.shape[-2]), a22.tostring())
            return ii
        elif a.shape[0] == 3:
            assert a.dtype==N.uint8
            a22 = N.transpose(a,(1,2,0)) # .copy()
            a22 = a22[:,:,rgbOrder[:3]]
            ii = Image.fromstring("RGB", (a.shape[-1],a.shape[-2]), a22.tostring())
            return ii
        elif a.shape[0] == 4:
            assert a.dtype==N.uint8
            a22 = N.transpose(a,(1,2,0)) # .copy()
            a22 = a22[:,:,rgbOrder[:4]]
            ii = Image.fromstring("RGBA", (a.shape[-1],a.shape[-2]), a22.tostring())
            return ii
        else:
            raise ValueError, "only 2d greyscale or 3d (RGB[A]) supported"
    # else:  (see return above)
    if a.ndim != 2: 
        raise ValueError, "only 2d greyscale or 3d (RGB[A]) supported"

    if a.dtype.type == N.uint8:
        mode = "L"
    elif a.dtype.type == N.float32:
        mode = "F"
    elif a.dtype.type in ( N.int16, N.uint16 ):
        mode = "I;16"
    else:
        raise ValueError, "unsupported array datatype"
    return Image.fromstring(mode, (a.shape[1], a.shape[0]), a.tostring())
    #20040929 todo: try this:   return Image.frombuffer(mode, (a.shape[1], a.shape[0]), a._data)




def saveTiffMultipage(arr, fn, rescaleTo8bit=False, rgbOrder="rgba", **params):
    """
    extension to PIL save TIFF
    **params is directly forwarded to PIL save function
    """
    if arr.ndim == 4:
        if arr.shape[1] not in (1,2,3,4):
            raise ValueError, "can save 4d arrays (color) only with second dim of len 1..4 (RG[B[A]])"
    elif arr.ndim != 3:
        raise ValueError, "can only save 3d (grey) or 4d (color) arrays"

    fp = open(fn, 'w+b')

    ifd_offsets=[]

    if rescaleTo8bit:
        mi,ma = float(arr.min()), float(arr.max())
        ra = ma-mi        

    params["_debug_multipage"] = True
    for z in range(arr.shape[0]):
        if rescaleTo8bit:
            a=(arr[z]-mi)*255./ra
            ii = array2image(a.astype(N.uint8), rgbOrder=rgbOrder)
        else:
            ii = array2image(arr[z], rgbOrder=rgbOrder)

        fp.seek(0,2) # go to end of file
        if z==0:
            # ref. PIL  TiffImagePlugin
            # PIL always starts the first IFD at offset 8
            ifdOffset = 8
        else:
            ifdOffset = fp.tell()

        ii.save(fp, format="TIFF", **params)
        
        if z>0: # correct "next" entry of previous ifd -- connect !
            ifdo = ifd_offsets[-1]
            fp.seek(ifdo)
            ifdLength = ii._debug_multipage.i16(fp.read(2))
            fp.seek(ifdLength*12,1) # go to "next" field near end of ifd
            fp.write(ii._debug_multipage.o32( ifdOffset ))

        ifd_offsets.append(ifdOffset)
    fp.close()




def saveImg(arr, fn, forceMultipage=False, rgbOrder="rgba"):
    """
    Saves data array as image file (format from    extension !! .tif,.jpg,...)
    tries to use number format of 'arr'
    also supports multipage TIFF:
        3D arrays: grey (if more than 4 z-secs or forceMultipage==True)
        4D arrays: color (second dim must be of len 2..4 (RG[B[A]])"

    for multi-color images:
         rgbOrder: order in which axes are mapped to RGB(A) channels
      
    !!be careful about up-down orientation !!
    """

    arr = N.asarray(arr)
    if (arr.ndim == 3 and (len(arr)>4 or forceMultipage)) or \
            arr.ndim == 4:
        return saveTiffMultipage(arr, fn, rescaleTo8bit=False, rgbOrder=rgbOrder)

    im = array2image(arr, rgbOrder)
    im.save(fn)


def _saveSeq_getFixedFN(fn, n):
    """
    check that fn contains a '%02d'-like part'
    autofix if necessary (add enough digits to fit n filenames)
    """
    try:
        __s = fn % 1 # test if fn contains '%d'
    except TypeError:
        import os
        fnf = os.path.splitext(fn)
        fns = '_%0' + '%d'%(int(N.log10(n-1))+1) +'d'
        fn = fnf[0] + fns + fnf[1]
    return fn
    
def saveImg_seq(arr, fn, rgbOrder="rgba"):
    """
    Saves 3D data array as 8-bit gray image file sequence (format from  extension !! .tif,.jpg,...)
    filename should contain a "template" like %02d - use '%%' otherwise inplace of single '%'
    template gets replaced with 00,01,02,03,...

    for multi-color images:
         rgbOrder: order in which axes are mapped to RGB(A) channels
      
    !!be careful about up-down orientation !!
    """
    arr = N.asarray(arr)
    #if arr.ndim != 3:
    #    raise "can only save 3d arrays"
    if not (arr.ndim == 3 or (arr.ndim == 4 and arr.shape[1] in (1,2,3,4))):
        raise ValueError, "can only save 3d arrays or 4d with second dim of len 1..4 (RG[B[A]])"

    fn = _saveSeq_getFixedFN(fn, len(arr))

    for i in range(arr.shape[0]):
        saveImg(arr[i], fn % i, rgbOrder=rgbOrder)
    
# Clement: old version
# def topPercentile(img, percentile=1):
#     """find Intens. for highest percentile
# 
#         slow!! ****** might only work for uint16 arr *********
#     """
#     a = N.empty( shape=(1<<16), dtype=N.int32 ) # bins
#     (mi,ma,mean,stddev) = S.histogram2(img, 0, (1<<16), a)
#     nPix = N.prod( img.shape )
# 
# 
#     a = _getGoodifiedArray(a)
# 
#     tp = S.toppercentile(a, nPix, int(ma), percentile)
#     return tp

def topPercentile(img, percentile=1):
    return N.percentile(img,100-percentile)


def asFloat32(a):
    """returns N.asarray(a, N.Float32)"""
    return N.asarray(a, N.float32)



