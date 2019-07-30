### EHom (20130611): deprecate old workaround using numpy.fft
### Use pyFFTW module instead: http://hgomersall.github.io/pyFFTW/sphinx/tutorial.html
### Using the quick and easy numpy.fft "drop-in" approach for now
### Could probably speed up by using pyFFTW syntax directly
### See Henry Gomerall's post about this wrapper code: http://hgomersall.wordpress.com/2012/02/01/the-joys-of-cython-numpy-and-a-nice-fftw-api/
### If using a Mac, see  post about pyFFTW and FFTW on a Mac: http://dawes.wordpress.com/2012/03/31/python-fourier-mac/
#

import numpy as _N

def numpy_rfft(a,af=None, inplace=0):
    if af is None:
        return _N.fft.rfftn(a)
    else:
        af[:] = _N.fft.rfftn(a)
def numpy_irfft(af, a=None, inplace=0, copy=1):
    #if a is None:
    #    return _N.fft.irfftn(af)
    #else:
        a[:] = _N.fft.irfftn(af)
        a *= _N.product(a.shape) # to fake FFTW's behavior of not normalizing
def numpy_fft(a,af=None, inplace=0):
    if af is None:
        return _N.fft.fftn(a)
    else:
        af[:] = _N.fft.fftn(a)
def numpy_ifft(af, a=None, inplace=0):
    #if a is None:
    #    return _N.fft.ifftn(af)
    #else:
        a[:] = _N.fft.ifftn(af)
        a *= _N.product(a.shape) # to fake FFTW's behavior of not normalizing

### NEED TO ACTIVATE AFTER CONFERRING WITH SEB
### This requires that FFTW is installed on the computer
### One can do this with (check for latest FFTW source here: ftp://ftp.fftw.org/pub/fftw/):
#        wget http://www.fftw.org/fftw-3.3.3.tar.gz
#         tar -zxvf fftw-3.3.3.tar.gz
#        cd fftw-3.3.3/
#        ./configure
#        sudo make install

import AIDA_Settings as Set

try:
    import pyfftw
    Set.fft_type = 'pyFFTW'
    def pyfftw_rfft(a,af=None, inplace=0):
        if af is None:
             return pyfftw.interfaces.numpy_fft.rfftn(a)
        else:
             af[:] = pyfftw.interfaces.numpy_fft.rfftn(a)
    def pyfftw_irfft(af, a=None, inplace=0, copy=1):
            a[:] = pyfftw.interfaces.numpy_fft.irfftn(af)
            a *= _N.product(a.shape) # to fake FFTW's behavior of not normalizing
    def pyfftw_fft(a,af=None, inplace=0):
        if af is None:
           return pyfftw.interfaces.numpy_fft.fftn(a)
        else:
            af[:] = pyfftw.interfaces.numpy_fft.fftn(a)
    def pyfftw_ifft(af, a=None, inplace=0):
            a[:] = pyfftw.interfaces.numpy_fft.ifftn(af)
            a *= _N.product(a.shape) # to fake FFTW's behavior of not normalizing
    
    # Turn on the cache for optimum performance
    pyfftw.interfaces.cache.enable()
    
    # Keep cache around for 1 second
    pyfftw.interfaces.cache.set_keepalive_time(1)
    
except ImportError:
    print(" ** pyFFTW cannot be loaded/found **")
    Set.fft_type = 'numpy_fft'

