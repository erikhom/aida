################################################################################
#
#   File:       fftfuncs.py
#
#   Summary:    Contains fft functions, to replace the use of Priithon
#
#   Authors:    Clement Chalumeau (SETI, interning for Franck Marchis) 
#
#        
################################################################################

import numpy as N

def zeroArr(dtype, *shape):
    """allocates and returns array of given dtype and shape"""
    if type(shape[0]) == tuple:
        shape = shape[0]
    return N.zeros(dtype=dtype, shape=shape)

def zeroArrD(*shape):
    """allocates and returns 'double prec. float' array of given shape"""
    return zeroArr(N.float64, *shape)

defshape = (256,256)  #use this default shape to make testing easy


def noiseArr(shape=defshape, stddev=1., mean=0.0, dtype=N.float32):
    return N.random.normal(mean, stddev,shape).astype(dtype)

def rfft(a, minFdtype=N.float32):
    '''
    calculate nd fourier transform
    performs real- fft, i.e. the return array has shape with last dim halfed+1

    `a` should be a real array,
      otherwise it gets converted to
      minFdtype
    '''
    if a.dtype.type not in (N.float32, N.float64):
        a = N.asarray(a, minFdtype)
    
    import fftw
    return fftw.rfft(a)


def gaussian(r, dim=1, sigma=1., integralScale=None, peakVal=None):
    """returns n-dim Gaussian centered at 0
    if integralScale is not None
         result.sum() == integralScale
    if peakVal       is not None
         result max is peakVal
    if both are None
         results defaults to integralScale==1
    """
    #sigma = _normalize_ndparm(sigma)
    # sigma = float(sigma)
    if integralScale is None and peakVal is None:
        integralScale = 1.
    if integralScale is not None:
        # *N.prod(sigma))
        f = 1./ ( N.sqrt(2.*N.pi) * sigma)**dim   * integralScale
    elif peakVal is not None:
        f = peakVal

    f2 = 2. * sigma**2
    return f*N.exp(-r**2/f2)

