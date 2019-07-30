################################################################################
#
#   File:       AIDA_CostGradFunctions.py
#
#   Summary:    Contains Cost Function and functions to compute the necessary
#               gradients for conjugate gradient minimization in the
#               Adaptive Image Deconvolution Algorithm (AIDA) package (AiDAP)
#
#   Authors:    Erik F.Y. Hom (Sedat Lab, UCSF) with assistance from
#               Timothy K. Lee (Summer Intern, Sedat Lab)
#
#   Other:      See 'AIDA_version.py' for date and version details
#               See 'LICENSE.txt' for use, license, and reference details
#
#   Modifications : Clement Chalumeau (SETI intern) 2016/04/06
#                Remove all priithon occurences 
################################################################################
# calculate all cost functions + gradients  equations 15-16 and 22-23 of Hom et al 2007

#from Priithon.all import *
#Used from Priithon.all:
# - N (numpy)
# - U (useful): only nd (ndimages from scipy)
import numpy as N
import Priithon_heritage.useful as U

#from Priithon import fftw
from Priithon_heritage import fftw

# Clement
from scipy.signal import fftconvolve

import AIDA_Settings as Set
import AIDA_GFunctions as AGF
### line below is no longer needed, since we have bundled the CCG optimizer
### with Priithon; see that for CCG calling details
#import AIDA_CCG as CCG         # this is a SWIG python interface that calls
                                # Dennis Goodman's (LLNL) CCG C++ code

#<< vvv BELOW IS FOR LATER INTERGRATION WITH AIDA SI CODE vvv >>#
if Set.decon_type in ('si', 'siclassical'):

    from AIDA_SI_Functions import Calculate_h_J_and_gradJ_SI
    from AIDA_SI_Functions import Calculate_o_J_and_gradJ_SI
#<< ^^^^^^^ >>#


#######  FUNCTION: _CostFunction  #######
####    [checked]
def CostFunction(xo, grad):
    """
    Use this function as input to the CCG code of D.M. Goodman @ LLNL

    Original void function needed by the CCG code:
        "funct(nc, xo, fn, gn, iqflag)"

    Original Input:
    'nc' = Nd = number of variables (unchanged), this is used in
        Calculate_h_J_and_gradJ
    'xo' = 1D vector of variables = "object" OR "PSF" (input values)

    Rewrite below...
    Additional Set parameter inputs:
    'shape' = shape of data - to convert from 1D xo to nD data
    'costfunction_type' = specification of the costfunction to use for
            different deconvolution types (classical, myopic, nobjects,
            npsfs)
    'OTF' = Fourier Transform of PSF estimate
    'grad' = gradient array (= gn; see below)
    'OBJECT' = Fourier Transform of object estimate
    'IMAGE' = Fourier Transform of image
    'axesFT' = axes in which to perform FFTs

    (Original) Returns:
    'fn' = function value at xo = "J"
    'gn' = 'grad' = 1D vector of the gradient at xn
    'iqflag' = inquiry flag
    """

    grad[:] = 0.
    Set.J[:] = 0.

    ###  Minimize PSF = xo
    if Set.costfunction_type == 0:

        #<< testing and debugging block >>#
        #try:
        #
        #   xo /= N.sum(xo)        ## normalize and scale PSF
        #except:
        #
        #   pass

        grad.shape = xo.shape = Set.shape
        fftw.rfft(a=xo, af=Set.OTF, inplace=False)
        fn = Calculate_h_J_and_gradJ(grad) #@

    ###  Minimize Object = 'xo'
    elif Set.costfunction_type == 1:

        #<< testing and debugging line >>#
        #Set.J[0] = Set.J[1] = 0.
        
        grad.shape = xo.shape = Set.shape
        fftw.rfft(a=xo, af=Set.OBJECT, inplace=False)
        fn = Calculate_o_J_and_gradJ(xo, grad) #@

    ###  Minimize PSF = 'xo' in the Presence of Nobjects
    elif Set.costfunction_type == 2:

        #<< testing and debugging block >>#
        #Set.J[0] = Set.J[2] = 0.
        #try:
        #
        #   xo /= N.sum(xo)        ## normalize and scale PSF
        #except:
        #
        #   pass

        grad.shape = xo.shape = Set.shape
        fn = 0.
        fftw.rfft(a=xo, af=Set.OTF, inplace=False)

        if Set.memory_usage_level > 0:

            for i in range(len(Set.image_list)):

                Set.image = Set.image_array[i]
                Set.inv_w = Set.inv_w_array[i]
                Set.OBJECT = Set.OBJECT_array[i]
                ## cost function is calculated as sum over each PSF estimate
                fn += Calculate_h_J_and_gradJ(grad) #@
        else:

            for i in range(len(Set.image_list)):

                Set.image = Set.image_array[i]
                Set.inv_w = Set.inv_w_array[i]
                fftw.rfft(a=Set.object_array[i], af=Set.OBJECT, inplace=False)
                ## cost function is calculated as sum over each PSF estimate
                fn += Calculate_h_J_and_gradJ(grad) #@

    ###  Minimize Object = 'xo' in the Presence of Npsfs
    elif Set.costfunction_type == 3:

        #<< testing and debugging line >>#
        #Set.J[0] = Set.J[1] = 0.
        
        grad.shape = xo.shape = Set.shape
        fn = 0.
        fftw.rfft(a=xo, af=Set.OBJECT, inplace=False)

        if Set.memory_usage_level > 0:

            for i in range(len(Set.image_list)):

                Set.image = Set.image_array[i]
                Set.inv_w = Set.inv_w_array[i]
                Set.inv_theta = Set.inv_theta_center_array[i]
                Set.lambda_object = Set.lambda_object_Narray[i]
                Set.mu = Set.mu_Narray[i]
                Set.OTF = Set.OTF_array[i]
                ## cost function is calculated as sum over each object estimate
                fn += Calculate_o_J_and_gradJ(xo, grad) #@
        else:

            for i in range(len(Set.image_list)):

                Set.image = Set.image_array[i]
                Set.inv_w = Set.inv_w_array[i]
                Set.inv_theta = Set.inv_theta_center_array[i]
                Set.lambda_object = Set.lambda_object_Narray[i]
                Set.mu = Set.mu_Narray[i]
                fftw.rfft(a=Set.PSF_array[i], af=Set.OTF,
                        inplace=False)
                ## cost function is calculated as sum over each object estimate
                fn += Calculate_o_J_and_gradJ(xo, grad) #@

#<< vvv BELOW IS TEMPORARY UNTIL WE INTEGRATE IT WITH SI CODE vvv >>#
    ###  Minimize Object = 'xo' for Structured Illumination
    elif Set.costfunction_type == 4:

        grad.shape = xo.shape = Set.shape
        #fftw.rfft(a=xo,af=Set.OBJECT, inplace=False)
        fn = Calculate_o_J_and_gradJ_SI(xo, grad) #@

    ###  Minimize PSF = 'xo' for Structured Illumination
    elif Set.costfunction_type == 5:

        grad.shape = xo.shape = Set.shape_small
        fftw.rfft(a=xo, af=Set.OTF_array[Set.PSForder], inplace=False)
        fn = Calculate_h_J_and_gradJ_SI(grad) #@
#<< ^^^^^^ >>#

    else:

        print "Error!!  Check `type` input flag!"

        fn = 0

    iqflag = 0 # dummy flag in this case

    return (fn, iqflag)
####



########  FUNCTION: Calculate_h_J_and_gradJ  #######
####    [checked]
def Calculate_h_J_and_gradJ(grad):
    """
    Returns the cost function, J and modifies its derivative with respect to
    the PSF, h, in the case where the _PSF_ is to be minimized

    Note that not all cost_function terms are calculated; only those that
    vary with respect to what is being differentiated/minimized are computed
    """

    ###  Data Fidelity Term
    if Set.terms[0] == 1:

        # cost function: Jn
        fftw.irfft(af=((Set.OBJECT*Set.OTF)*Set.inv_Nd), a=Set.temp,
                inplace=False, copy=False)
        # # Clement: test
        # obj = N.zeros((256,256))
        # fftw.irfft(af = Set.OBJECT , a = obj, inplace =False)
        # Set.temp = fftconvolve(Set.PSF,obj, mode = 'same')
        Set.temp -= Set.image
        Set.J[0] += 0.5*N.sum( ((Set.temp * Set.temp) * Set.inv_w).flat )

        ## gradient: Grad_h_Jn
        fftw.rfft(a=(Set.temp*(Set.inv_w)), af=Set.temp_realFFT, inplace=False)
        fftw.irfft(af=(N.conjugate(Set.OBJECT)*Set.temp_realFFT*Set.inv_Nd),
                a=Set.temp, inplace=False, copy=False)
        grad += Set.temp


    ###  OTF Harmonic Constraint Term
    if Set.terms[2] == 1:

        ## gradient: Grad_o_Jn
        ## using _corrected_ derivative expression:
        ## deriv = (1/2)*(Nd+1)*invFT[(OTF-mOTF)/v]

        fftw.irfft(af=(Set.lambda_OTF * Set.OTF_deriv_scale * \
                (Set.OTF - Set.mOTF) * Set.inv_v), a=Set.temp, inplace=False,
                copy=False)

        ## N.B.: below is if scaling is already included in FFT
        ##      this is not the case for FFTW, so Set.Nd ends up
        ##      canceling out above
        grad += Set.temp

        ## cost function: Jh
        Set.temp_realFFT_real = N.abs(Set.OTF - Set.mOTF)
        Set.temp_realFFT_real *= Set.temp_realFFT_real * Set.inv_v
        Set.J[2] +=  0.5 * Set.lambda_OTF * \
                (2. * N.sum(Set.temp_realFFT_real.flat) - \
                N.sum(Set.temp_realFFT_real[...,0].flat) - \
                N.sum(Set.temp_realFFT_real[...,-1].flat))

    ### PSF Harmonic Constraint Term (can be changed to some other constraint
    ### in the future, e.g., autocorrelation ala Pixon/Puetter approach)
    if Set.terms[3] == 1:
        # print '\n*****\n*****\n*****\n*****\n*****\n*****\n Calcul de la fonction pSF\n'

    ## THIS TERM APPEARS TO BE A STRONGER CONSTRAINT THEN THE OTF CONSTRAINT
    ## OTF CONSTRAINT SEEMS TO BE PRETTY TOLERANT OF DIFFERENCES IN THE OTF
    ## (PERHAPS TOO TOLERANT)
        Set.temp = Set.lambda_PSF * (Set.PSF - Set.mPSF) * Set.inv_u
        grad += Set.temp
        # Clement: error here !!!
        # Set.J[3] += 0.5 * N.sum((grad * (Set.PSF - Set.mPSF)).flat)
        Set.J[3] += 0.5 * N.sum((Set.temp * (Set.PSF - Set.mPSF)).flat)

    ## N.B.  PSF is being estimated, therefore ignore Jo calculation
    return N.sum(Set.J)
####


#######  FUNCTION: Calculate_o_J_and_gradJ  #######
####    [checked]
def Calculate_o_J_and_gradJ(object, grad):
    """
    Returns the cost function of J and modified its derivative with respect
    to the object, o

    Note that not all cost_function terms are calculated; only those that
    vary with respect to what is being differentiated/minimized are computed
    """

    Set.Delta1[:] = 0.
    Set.Laplacian[:] = 0.

    ###  Data Fidelity Term
    if Set.terms[0] == 1:

        ## cost function: Jn
        fftw.irfft(af=((Set.OBJECT*Set.OTF)*Set.inv_Nd), a=Set.temp,
                inplace=False, copy=False)
        # Clement: test
        # Set.temp = fftconvolve(Set.PSF,object, mode = 'same')
        Set.temp -= Set.image
        Set.J[0] += 0.5*N.sum( ((Set.temp * Set.temp) * Set.inv_w).flat )


        ## gradient: Grad_o_Jn
        fftw.rfft(a=(Set.temp*(Set.inv_w)), af=Set.temp_realFFT, inplace=False)
        fftw.irfft(af=(N.conjugate(Set.OTF)*Set.temp_realFFT*Set.inv_Nd),
                a=Set.temp, inplace=False, copy=False)
        grad += Set.temp


    ### Edge-preserving object prior
    if Set.terms[1] == 1:

        if Set.laplacian_operator == 0:

            if Set.norm == 'sqrt':      ## MISTRAL object gradient expression

                ## norm = sqrt[ (grad_x)^2 + (grad_y)^2 + ...
                for ax in range(Set.dimension):

                    U.nd.convolve(object, Set.kernel_grad[ax],
                            output=Set.object_gradient, mode='wrap')

                   #<< testing and debugging line >>#
                   #U.nd.convolve(Set.test_true_object, Set.kernel_grad[ax],
                   #         output=Set.object_gradient, mode='wrap')

                    ## temp Delta1:
                    Set.Delta1 += Set.object_gradient * Set.object_gradient
                    U.nd.convolve(Set.object_gradient, Set.kernel_lapl[ax],
                            output=Set.temp, mode='wrap')
                    Set.Laplacian += Set.temp
                ## final Delta1:
                Set.Delta1 = 1. + N.sqrt(Set.Delta1)*Set.inv_theta
            elif Set.norm == 'abs':         # approximation of gradient

                ## norm ~ (1/sqrt[axes]) * [ abs(grad_x) + abs(grad_y) + ... ]
                for ax in range(Set.dimension):

                    U.nd.convolve(object, Set.kernel_grad[ax],
                            output=Set.object_gradient, mode='wrap')
                    ## temp Delta1:
                    Set.Delta1 += N.abs(Set.object_gradient)
                    U.nd.convolve(Set.object_gradient, Set.kernel_lapl[ax],
                            output=Set.temp, mode='wrap')
                    Set.Laplacian += Set.temp
                ## final Delta1:
                Set.Delta1 = 1. + ( (1/N.sqrt(Set.dimension)) * \
                                    Set.Delta1 * Set.inv_theta )
        else:

            if Set.norm == 'sqrt':

                for ax in range(Set.dimension):

                    U.nd.convolve(object, Set.kernel_grad[ax],
                                output=Set.object_gradient, mode='wrap')

                   #<< testing and debugging line >>#
                   #U.nd.convolve(Set.test_true_object, Set.kernel_grad[ax],
                   #         output=Set.object_gradient, mode='wrap')

                    ## temp Delta1:
                    Set.Delta1 += Set.object_gradient * Set.object_gradient
                ## final Delta1:
                Set.Delta1 = 1. + (N.sqrt(Set.Delta1) * Set.inv_theta)
            elif Set.norm == 'abs':

                for ax in range(Set.dimension):

                    U.nd.convolve(object, Set.kernel_grad[ax],
                            output=Set.object_gradient, mode='wrap')
                    ## temp Delta1:
                    Set.Delta1 += N.abs(Set.object_gradient)
                ## final Delta1
                Set.Delta1 = 1. + ( (1/N.sqrt(Set.dimension)) * \
                                    Set.Delta1 * Set.inv_theta )

            U.nd.convolve(object, Set.kernel_lapl, output=Set.Laplacian,
                    mode='wrap')

        if Set.zero_edges_flag:

            ndim = Set.Delta1.ndim
            AGF.SetEdgeValues(Set.Delta1, (ndim-Set.ndim_lowering), value=1)
                    #@ AGF
            AGF.SetEdgeValues(Set.Laplacian, (ndim-Set.ndim_lowering), value=0)
                    #@ AGF

        ## gradient: Grad_o_Jo
        Set.temp = Set.Laplacian / Set.Delta1
        grad += Set.mu*Set.temp

        ## cost function: Jo
        Set.J[1] += N.sum( (Set.lambda_object*((Set.Delta1-1.) - \
                N.log(Set.Delta1))).flat )

    # N.B.  object is being estimated, therefore ignore Jh calculation
    return N.sum(Set.J)
####



