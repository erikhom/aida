################################################################################
#
#   File:       AIDA_Default_Settings.py
#
#   Summary:    This file contains the COMPLETE (and commented) set of default 
#               settings for AIDA.  These defaults are overwritten by the values
#               in 'AIDA_Settings.py' which can be overwritten by values in a 
#               settings file supplied by the user and specified by a '-S' 
#               cammand-line option
#
#				~/anaconda/bin/python ~/AIDA1.4/AIDA.py -S AIDA_Settings.py
#
#   Author:     Erik F.Y. Hom (Sedat Lab, UCSF)
#
#   Other:      See 'AIDA_version.py' for date and version details
#               See 'LICENSE.txt' for use, license, and reference details
#
################################################################################
##############  DO NOT MODIFY UNLESS YOU KNOW WHAT YOU'RE DOING!!  #############
################################################################################


#########################################
##11111  GENERAL AIDA PARAMETERS  11111##

#>>> \(PREFIX OF) IMAGE FILE TO DECONVOLVE (INCLUDE PATH)
image_filenames = ''         # (default = '')


#>>> \(PREFIX OF) PSF FILE TO USE FOR DECONVOLUTION (INCLUDE PATH)
PSF_filenames = ''           # (default = '')


#>>> \KEY TO LABEL DECONVOLUTION FILES
dataset_label = None            # (default = None)


#>>> \FILE_ARRAY OF TUPLES FOR DECONVOLUTION
dataset_list = None         # format should be a list of tuples:
                            # [('image_file', 'PSF_file', {'key'}), ...]
                            # where 'key' can be an optional 3rd tuple entry 

#>>> AMOUNT OF INFO TO OUTPUT DURING THE RUN
#info_level = 0     # no printing
#info_level = 1     # minimal printing
info_level = 2      # moderate printing (default)
#info_level = 3     # full printing

#>>> AMOUNT OF MEMORY USAGE THAT SHOULD BE USED
#memory_usage_level = 0         # use the smallest amount of memory
                                # possible (slowest)
#memory_usage_level = 1         # IMAGE(s) are kept in memory
#memory_usage_level = 2         # image stats arrays are kept in memory
                                    # as well
memory_usage_level = 3          # object and/or PSF output results are 
                                # kept in memory; particularly relevant
                                # in 'Nobjects' and 'Npsfs' decon
                                # (default)

#>>> COST FUNCTION TERMS TO USE
terms = (1,1,1,0)       # (default = (1,1,1,0))
                        # (1,0,0,0): data fidelity term
                        # (0,1,0,0): object regularization term
                        # (0,0,1,0): OTF harmonic constraint term
                        # (0,0,0,1): autocorrelation term (not 
                        # implemented (yet and may never be!?))

#>>> WHERE TO STORE THE RESULTS
results_directory = '.'     # (default = '.' 
                            # (current working directory))

#>>> WHETHER TO STORE RESULTS IN TIMESTAMPED RESULTS_DIRECTORY
timestamp_directory_flag = True     # (default = True)

#>>> TYPE OF DECONVOLUTION TO RUN
#decon_type = 'classical'   # object varies and PSF is fixed
decon_type = 'myopic'       # myopic decon: object and PSF 
                            # vary (default)
#decon_type = 'nobjects'    # multiframe myopic decon of 
                            # N different objects with one
                            # common PSF
#decon_type = 'npsfs'       # multiframe myopic decon of
                            # one common object and N 
                            # different PSFs

#>>> INITIAL STARTING PSF FOR OPTIMIZATION		### EHom (20130611): added this to allow initial PSF to be inputted PSF value
initial_psf = 'mean'		# uses mPSF: average of PSFs in PSF list
# initial_psf = 'initial'	# uses input PSFs as intial PSF guesses

#>>> NUMERICAL PRECISION TO USE
precision = 'double'        # 64-bit (default)
#precision = 'single'       # 32-bit

#>>> \PATIAL DIMENSION OF IMAGE DATA
dimension = 2               # either 2D or 3D (default = 2)

#>>> \OUTPUT FILE FORMAT FOR DECONVOLUTION RESULTS
output_format = None        # same as image data format
                            # options: 'f' for FITS; 'm' for MRC; 't' for TIFF
                            # (default = None)

#>>> WHETHER TO OUTPUT CLEAN PSFS, MEAN PSF, AND V FILE
output_intermediate_files_level = 0     # (default; don't output)
#output_intermediate_files_level = 1    # (output mPSF and v)
#output_intermediate_files_level = 2    # (output all, including cleaned PSFs)



#####################################################
##22222  IMAGE & NOISE PROCESSING PARAMETERS  22222##

#>>> WHETHER TO RESIZE THE IMAGES IN XY-DIMENSIONS TO BE 2^n
resizeimage_flag = True         # default = True


#>>> WHETHER TO RESIZE THE IMAGES IN Z-DIMENSION FOR 3D DATA TO BE 2^n
zresize_flag = False            # default = False

#>>> INITIAL OBJECT GUESS
initial_object_guess = 'zero'       # (default)
#initial_object_guess = 'wiener'    # wiener filter guess
#initial_object_guess = 'image'     # image as initial guess


#>>> FACTOR TO MULTIPLY WIENER PARAMETER WHEN INITIAL GUESS IS 'WIENER'
wiener_multiply_factor = 10.        # (default = 10.)

#>>> BACKGROUND SUBTRACTION VALUE; IDEALLY, THIS HAS ALREADY BEEN DONE
background = 0.     # images should be "cleaned"; should have the average
                        # of detector noise subtracted so that there are 
                        # negative pixels from which to calculate the Gaussian
                        # statistics of detector noise (default = 0.)

#>>> DARK IMAGE FILE; ONE OPTION USED TO COMPUTE SIGMA_DETECTOR VALUE
dark_image = None       # (default = None (i.e., not to be used))

#>>> SIGMA_DETECTOR VALUE
sigma_det = None        # value for the standard deviation of the Gaussian
                        # noise component
                        # (default = None (i.e., not to be used))



###########################################
##33333  PSF PROCESSING PARAMETERS  33333##

#>>> FOOTPRINT TO USE FOR CENTERING FIT
exclude_radius = 5          # radius from center of PSF to exclude from
                            # clean-up procedure (default = 5)

#>>> HOW PSFS SHOULD BE CLEANED BY AIDA
cleaning = (1,1,1,1)        # controls how PSF images provided are to be
                            # cleaned to generate 'theoretical' PSFs
                            #       (0,0,0,1) - normalize by sum
                            #       (0,0,1,0) - threshold to zero
                            #       (0,1,0,0) - subtract 'clean_sigma'/ 
                            #                   'clean_threshold_percent'
                            #       (1,0,0,0) - remove bad pixels
                            # (default = (1,1,1,1))

#>>> CENTERING STATUS OF SAMPLED PSFS
PSF_centering = 'unknown'       # 'unknown' (default) will compute center
                                #           based on a centroid approximation
                                # 'origin' will assume PSFs are centered at
                                #           (0,0, (0))
                                # 'array_center' will assume PSFs have centers
                                #           at the center of the image array

#>>> PERCENT OF PSFMAX FOR WHICH VALUES OF THE PSF THAT ARE EQUAL OR BELOW
#    SHOULD BE SET TO 'SET.FILL'
PSF_subtract_background_percent = 0.    # used in the cleaning process for PSFs
                                        # (default = 0.)

#>>> NUMBER OF BACKGROUND NOISE STANDARD DEVIATIONS TO SUBTRACT OFF OF A PSF
#    IMAGE
nsigmas = 2.        # calculated from negative pixels of the PSF
                    # after thresholding with 'clean_threshold_percent'
                    # used in the cleaning process for PSFs
                    # (default = 2)

#>>> LOWER BOUND PERCENT OF MAX BELOW WHICH PSF VALUES WILL BE SET TO 'SET.FILL'
#>>> USED IN CLEANING PROCESS
PSF_threshold_percent = 1e-7    # default = 1e-7
                                # presume 1/100 of 16-bit (2**16) dynamic range
                                # based on CCD camera sensitivity (1e-5)

#>>> LOWER BOUND PERCENT OF MAX BELOW WHICH OTF MAGNITUDE VALUES WILL BE SET
#>>> TO 'SET.FILL'; USED IN CLEANING PROCESS
OTF_threshold_percent = 1e-7    # default = 1e-7

#>>> VALUE USED TO FILL IN FLOOR VALUES OF PSF AND OTF
fill = 0.           # default = 0.

#>>> LOWER BOUND VALUE FOR THE PSF SAMPLE VARIANCE; FOR REGIONS WHERE U = 0
#>>> EXPRESSED AS A PER DATA-DIMENSION QUANTITY
u_floor_per_dim = 1e-6          # (default = 1e-6)

#>>> LOWER BOUND VALUE FOR THE POWER SPECTRAL DENSITY (VARIANCE) OF THE 
#>>> SAMPLED OTFs; ESPECIALLY FOR REGIONS WHERE V COULD = 0
#>>> EXPRESSED AS A PER DATA-DIMENSION QUANTITY
v_floor_per_dim = 1e-6          # (default = 1e-6)

#>>> LOWER BOUND OF WHICH V MAGNITUDE VALUES WILL BE SET TO 'SET.FILL';
v_threshold = 0.01          # default = 0.01

#>>> WHETHER TO USE A RADIALLY AVERAGED POWER SPECTRAL DENSITY FOR THE 
#    OTF HARMONIC CONSTRAINT
dimensions_to_radially_average_v = 0        # is/are the dimensions in which to 
                                            # radially average the OTF variance
                                            # this effectively makes the 
                                            # constraint a function of just
                                            # _frequency_ = |k| (and not each 
                                            # k component of the k-vector)
                                            # Values can be 0-3
                                            # (default = 0 - no averaging)

#>>> WHETHER TO OUTPUT PSFs AS ORIGIN CENTERED OR IMAGE CENTERED
origin_centered_PSF_flag = False        # if 'True', all PSFs written to file
                                        # will be centered about the origin
                                        # if 'False', the PSFs (and weights 'v')
                                        # will be centered about 
                                        # ([N/2,] N/2, N/2), the pixel center
                                        # of the image array
                                        # (default = False)



#############################################################################
##44444  DERIVATIVE & LAPLACIAN OPERATOR PARAMETERS FOR OBJECT PRIOR  44444##

#>>> DERIVATIVE OPERATOR 
#derivative_operator = 'pixel'          # nearest neighbor asymmetric finite
                                        # difference
#derivative_operator = 'symmetric'      # using symmetric finite difference
derivative_operator = 'FC'              # using Frei-Chen approximation
                                        #  (default)
                                        # For further details see
                                        # 'SelectKernel' function in
                                        # 'AIDA_Functions.py'

#>>> LAPLACIAN OPERATOR
laplacian_operator = 3          # should be 0-3
                                # 0: same as derivative operator (uses
                                #    object_gradient array)
                                # 1-3: use laplacian kernel operators 1-3
                                # (uses object array directly)
                                # (default = 3)
                                # For further details see 'SelectKernel'
                                # function in 'AIDA_Functions.py'

#>>> RHO RESOLUTION COMPENSATION FACTOR
rho = 1.            # ratio to compensate/weight gradient/laplacian calculation
                    # in the xy direction relative to a 3-pixel sized PSF
                    # Change only for extended PSFs for which default setting
                    # of 1.0 leads to over-regularized results

#>>> ZETA RATIO OF Z-RESOLUTION LOSS:XY-RESOLUTION
zeta = 3.           # ratio to compensate/weight gradient/laplacian calculation
                    # in the z direction relative to the xy-directions for
                    # 3D image deconvolution

#>>> GRADIENT NORM FORMULA
norm = 'sqrt'       # standard norm definition (default)
#norm = 'abs'       # approximation to the norm as sum of abs gradient

#>>> WHETHER EDGES SHOULD BE SET TO ZERO AFTER TAKING THE GRADIENT
zero_edges_flag = False     # (default = False)

#>>> USED IN CONJUNCTION WITH 'ZERO_EDGES' - SEE AIDA_FUNCTIONS.PY FOR DETAILS
ndim_lowering = 0       # (default = 0)



########################################################################
##55555  OPTIMIZATION AND PARTIAL CONJUGATE GRADIENT PARAMETERS  55555##

#>>> TOLERANCE FOR CONSECUTIVE PARTIAL CG ESTIMATES OF THE OBJECT
object_PCG_tolerance = 0.1              # default = 0.1
                                        # true tolerance
                                        # = object_PCG_tolerance 
                                        #   * 1/N_total_pixels
                                        #   # mean(image)

#>>> TOLERANCE FOR CONSECUTIVE PARTIAL CG ESTIMATES OF THE PSF
PSF_PCG_tolerance = 0.1                 # default = 0.1
                                        # true tolerance
                                        # = PSF_PCG_tolerance 
                                        #   * 1/N_total_pixels
                                        #   # mean(PSF)

#>>> PARTIAL CONJUGATE GRADIENT RESTART PROTOCOLS
PCG_iter_array = (1,3,5,7,9,11,13,15)   # default = (1,3,5,7,9,11,13,15)
object_PCG_iter_array = None            # default = None;
                                        # will be set to PCG_iter_array
PSF_PCG_iter_array = None               # default = None;
                                        # will be set to PCG_iter_array

#>>> TOTAL NUMBER OF PARTIAL CONJUGATE GRADIENT BLOCKS FOR DECON=CLASSICAL
max_classical_PCG_iter = (24,)          # (default = (24,))

#>>> TOTAL NUMBER OF PARTIAL CONJUGATE BLOCKS FOR OBJECT/PSF ESTIMATE
max_total_PCG_blocks = len(PCG_iter_array)
                                        # default = len(PCG_iter_array)

#>>> CAP ON NUMBER OF SEQUENTIAL 'STOP' SIGNALS ENCOUNTERED BEFORE BREAKING OUT 
#    OF A PARTIAL CONJUGATE GRADIENT BLOCK
max_sequential_PCG_stops = 2            # default = 2

#>>> CAP ON THE NUMBER OF TIMES OBJECT ESTIMATE CAN INCREASE CONSECUTIVELY
#    IN A PARTIAL CONJUGATE GRADIENT BLOCK
max_uphill_object_PCG_steps = 2         # default = 2

#>>> CAP ON THE NUMBER OF TIMES PSF ESTIMATE CAN INCREASE CONSECUTIVELY
#    IN A PARTIAL CONJUGATE GRADIENT BLOCK
max_uphill_PSF_PCG_steps = 2            # default = 2

#>>> TOLERANCE RATIO OF CONTINUING BEYOND CURRENT ESTIMATE AS A FUNCTION OF 
#    THE DIFFERENCE BETWEEN PREVIOUS TWO STEPS
rising_tol_ratio = 1.03                 # default = 1.03 (3%)

#>>> CAP ON THE NUMBER OF SEQUENTIAL PCG STOP SIGNALS ENCOUNTERED IN MINIMIZING
#   A PARTICULAR VARIABLE - BEYOND THIS, VARIABLE WILL NO LONGER BE OPTIMIZED
max_optimization_stops = 3              # default = 3

#>>> CAP ON THE NUMBER OF SEQUENTIAL 'RISING ESTIMATE' PCG SIGNALS ENCOUNTERED
#   IN MINIMIZING A PARTICULAR VARIABLE - BEYOND THIS, VARIABLE WILL NO LONGER
#   BE OPTIMIZED
max_rising_stops = 3                    # default = 3

#>>> FOR MULTIFRAME DECON, FRACTION OF FRAMES THAT MUST EACH SATISFY AIDA 
#    CONVERGENCE TOLERANCE CRITERIA BEFORE STOPPING
Nframes_fractional_convergence = 0.99   # default = 0.99



###########################################################
##66666  CONSTRAINED CONJUGATE GRADIENT PARAMETERS  66666##

#N.B. Details can be found in 'ccg.cpp' and 'cgg.h'

#>>> OBJECT CONSTRAINED CONJUGATE GRADIENT TOLERANCE
object_CCG_tolerance = 1e-7     # (default = 1e-7)

#>>> PSF CONSTRAINED CONJUGATE GRADIENT TOLERANCE
PSF_CCG_tolerance = 1e-7        # (default = 1e-7)

#>>> CAP ON OBJECT CONSTRAINED CONJUGATE GRADIENT ITERATIONS
max_object_CCG_iter = 24        # (default = 24)

#>>> CAP ON PSF CONSTRAINED CONJUGATE GRADIENT ITERATIONS
max_PSF_CCG_iter = 24       # (default = 24)

#>>> MINIMUM VALUE ALLOWED FOR VARIABLE TO BE OPTIMIZED
xmin_value = 0.     # (default = 0.)

#>>> FLAG INDICATING HOW TO HANDLE OPTIMIZATION OF VARIABLE 
ivec_value = 0      # see 'ccg.h' for more details
                    # (default = 0) ("bounded_no constraint" specification)



##################################################
##77777  HYPERPARAMETER RELATED PARAMETERS 77777##

#>>> \LAMBDA_OBJECT INPUT VALUE
lambda_object_input = None      # 'lambda_object' by default is calculated
                                # using a formula of E. Hom
                                # (default = None)

#>>> \THETA INPUT VALUE
theta_input = None          # 'theta' by default is calculated
                            # using a formula of E. Hom
                            # (default = None)

#>>> LAMBDA_OTF VALUE
lambda_OTF_input = None     # default = None; 
                            # will be set to 1./(number of pixels)

#>>> LAMBDA_OTF VALUE
lambda_PSF_input = None     # default = None; will be set to 1.

#>>> (EXPONENTIAL) GRID SEARCH: CENTER VALUE FOR LAMBDA_OBJECT GRID SEARCH
orig_lambda_object_center = 0.01    # default center value for lambda_object
                                    # used for grid search

#>>> (EXPONENTIAL) GRID SEARCH: MULTIPLICATION FACTOR FOR LAMBDA_OBJECT
lambda_object_multiply_factor = 2.      # (default = 2.)

#>>> (EXPONENTIAL) GRID SEARCH: EXPONENT FOR MULTIPLICATION FACTOR FOR 
#    LAMBDA_OBJECT
lambda_object_above_below_exp = 0.      # (default = 0.)
                                        # (i.e., NO grid search)

# N.B: theta grid searching is only implemented for classical and myopic 
# decon types
#>>> (EXPONENTIAL) GRID SEARCH: MULTIPLICATION FACTOR FOR THETA
theta_multiply_factor = 2.      # (default = 2.)

#>>>(EXPONENTIAL) GRID SEARCH: EXPONENT FOR MULTIPLICATION FACTOR FOR 
### THETA
theta_above_below_exp = 0.      # (default = 0.)
                                # (i.e., NO grid search)

# N.B. grid search over lambda_OTF is not implemented yet
#>>> (EXPONENTIAL) GRID SEARCH: MULTIPLICATION FACTOR FOR THETA
lambda_OTF_multiply_factor = 2.     # (default = 2.)

#>>>(EXPONENTIAL) GRID SEARCH: EXPONENT FOR MULTIPLICATION FACTOR FOR 
### THETA
lambda_OTF_above_below_exp = 0.     # (default = 0.)
                                # (i.e., NO grid search)

#>>> WHETHER TO RECENTER GRID SEARCH ON LAMBDA_OBJECT ESTIMATE
grid_search_lambda_object_estimate_flag = False     # if a grid search is
                                                    # specified, uses
                                                    # 'lambda_object'
                                                    # estimate as the center
                                                    # of the grid search
                                                    # (default = False)

#>>> SCALING FACTOR TO MULTIPLY LAMBDA_OBJECT ESTIMATE
lambda_object_scaling = 1.0     # in cases where estimation scheme 
                                # routinely over- or under-regularizes,
                                # lambda_object can be scaled by this
                                # factor (default = 1.0)

#>>> SCALING FACTOR TO MULTIPLY THETA ESTIMATE
theta_scaling = 1.0     # in cases where estimation scheme may over- or 
                        # under-regularize, theta can be scaled by this
                        # factor (default = 1.0)

#>>> SCALING FACTOR TO MULTIPLY LAMBDA_OTF ESTIMATE
lambda_OTF_scaling = 1.0        # in cases where estimation scheme may over- or 
                                # under-regularize/constrain the OTF, 
                                # lambda_OTF can be scaled by this factor
                                # (default = 1.0)

#>>> SCALING FACTOR TO MULTIPLY LAMBDA_PSF ESTIMATE
lambda_PSF_scaling = 1.0        # in cases where estimation scheme may over- or 
                                # under-regularize/constrain the PSF, 
                                # lambda_PSF can be scaled by this factor
                                # (default = 1.0)



################################################
##88888  MISCELLANEOUS OTHER PARAMETERS  88888##

#>>> USER SPECIFIED NOTES
notes = ''

#>>> WHETHER TO RUN IN DEBUG MODE (NOT FULLY IMPLEMENTED YET)
debug = False           # default = False
debug_numpy_error_settings = {'all': 'raise', 'under':'ignore'}	# EHom added 'under':'ignore' (20130605)

#>>> USE BANDLIMITED CONSTRAINT ONLY INSTEAD OF HARMONIC OTF CONSTRAINT
band_limited_constraint = False             # default = False

#>>> USE TRUE_OBJECT FOR DEBUGGING/TESTING (NOT FULLY IMPLEMENTED YET)
true_object_file = None                         # default = None

#>>> USE TRUE_OBJECT FOR DEBUGGING/TESTING (NOT FULLY IMPLEMENTED YET)
true_PSF_file = None                            # default = None

imageHdr = None


