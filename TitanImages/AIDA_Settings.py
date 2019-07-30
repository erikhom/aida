################################################################################
#
#   File:       AIDA_Settings.py
#
#   Summary:    This file contains user-defined settings to run AIDA
#               It can be copied elsewhere and called using the '-S' 
#               command line option, e.g.:
#
#                      ~/Priithon_25_mac/priithon_script ~/AIDA_1.2.2i/AIDA.py -S ./AIDA_Settings.py
#
#               or this file could be modified directly (but must reside with
#               the rest of the AIDA code.  Settings here will be used to
#               overwrite the default settings in 'AIDA_Default_Settings.py'
#               Any user specified settings file called with the '-S' option
#               will have priority and used to overwrite settings.  Settings
#               specified via command-line options (see 'AIDA.py' file)
#               will have the highest overwrite priority
#
#   Author:     Erik F.Y. Hom (Sedat Lab, UCSF)
#
#   Other:      See 'AIDA_version.py' for date and version details
#               See 'LICENSE.txt' for use, license, and reference details
#
################################################################################


from AIDA_Default_Settings import *


################################################################################
###############   UNCOMMENT / MODIFY SETTINGS BELOW AS NEEDED   ################
################################################################################

notes = ''

path = './'
results_directory = path + 'Results'

#image_filenames = ['j8fb01rqq_rec.fits']
image_filenames = ['titanhe_153_IF_scaled']

#PSF_filenames = ['psf00_psf.fits']
PSF_filenames = ['psf_titanhe']

dataset_label = ['myopic_Python2.5']
#dataset_label = ['dec_class']

decon_type = 'myopic'       # ['classical', 'myopic', 'nobjects', 'npsfs']
#decon_type = 'classical'

#background = 0.
#sigma_det = 5
#dark_image = None

### Alternative entry format, though not thoroughly tested:
#dataset_list = None        # format should be a list of tuples:
#                           # [('image_file', 'PSF_file', {'key'}), ...]
#                           # where 'key' can be an optional 3rd tuple entry 

#timestamp_directory_flag = True    # default = True
#dimension = 2                      # either 2D or 3D (default = 2)
#zeta = 3.          # ratio to compensate/weight gradient/laplacian calculation
                    # in the z direction relative to the xy-directions for
                    # 3D image deconvolution


################################################################################


#lambda_object_scaling = 1.0        # default = 1.0
#theta_scaling = 1.0                # default = 1.0
#lambda_OTF_scaling = 1.0           # default = 1.0
#lambda_PSF_scaling = 1.0           # default = 1.0


################################################################################
# Other options
################################################################################
### 11111 ###
#info_level = 2                         # default = 2  [0, 1]
#memory_usage_level = 3                 # default = 3  [0, 1, 2]
#terms = (1,1,1,0)                      # default = (1,1,1,0)
#precision = 'double'                   # default = 'double'  ['single']
#output_format = None                   # default = None  ['f', 'm']
#output_intermediate_files_level = 0    # default = 0

### 22222 ###
#resizeimage_flag = True                # default = True
#zresize_flag = False                   # default = False
#initial_object_guess = 'zero'          # default = 'zero'  ['image', 'wiener']
#...

### 33333 ###
#exclude_radius = 5                     # default = 5
#cleaning = (1,1,1,1)                   # default = (1,1,1,1)
#PSF_centering = 'unknown'              # default = 'unknown'  
                                        #           ['origin', 'array_center']
#PSF_subtract_background_percent = 0.   # default = 0.
#nsigmas = 2.                           # default = 2
#PSF_threshold_percent = 1e-7           # default = 1e-7
#OTF_threshold_percent = 1e-7           # default = 1e-7
#fill = 0.                              # default = 0.
#u_floor_per_dim = 1e-6                 # default = 1e-6
#v_floor_per_dim = 1e-6                 # default = 1e-6
#v_threshold = 0.01                     # default = 0.01
#dimensions_to_radially_average_v = 0   # default = 0  [2, 3]
#origin_centered_PSF_flag = False       # default = False

### 44444 ###
#derivative_operator = 'FC'             # default = 'FC  ['pixel', 'symmetric']
#laplacian_operator = 3                 # default = 3  [0, 1, 2]
#rho = 1.                               # default = 1.
#...

### 55555 ###
#object_PCG_tolerance = 0.1             # default = 0.1
#PSF_PCG_tolerance = 0.1                # default = 0.1
#PCG_iter_array = (1,3,5,7,9,11,13,15)  # default = (1,3,5,7,9,11,13,15)
#object_PCG_iter_array = None           # default = None; set to PCG_iter_array
#PSF_PCG_iter_array = None              # default = None; set to PCG_iter_array

#max_classical_PCG_iter = (24,)                 # default = (24,)
#max_total_PCG_blocks = len(Set.PCG_iter_array) # default = len(PCG_iter_array)
#max_sequential_PCG_stops = 2                   # default = 2
#max_uphill_object_PCG_steps = 2                # default = 2
#max_uphill_PSF_PCG_steps = 2                   # default = 2
#rising_tol_ratio = 1.03                        # default = 1.03
#max_optimization_stops = 3                     # default = 3
#max_rising_stops = 3                           # default = 3
#Nframes_fractional_convergence = 0.99          # default = 0.99

### 66666 ###
#object_CCG_tolerance = 1e-7        # default = 1e-7
#PSF_CCG_tolerance = 1e-7           # default = 1e-7
#max_object_CCG_iter = 24           # default = 24
#max_PSF_CCG_iter = 24              # default = 24
#xmin_value = 0.                    # default = 0.
#...

### 77777 ###
#lambda_object_input = None     # default = None
#theta_input = None             # default = None
#lambda_OTF_input = None        # default = None; set to 1/N_total_pixels
#lambda_PSF_input = None        # default = None; set to 1.

#orig_lambda_object_center = 0.01           # default = 0.01
#lambda_object_multiply_factor = 2.         # default = 2.
#lambda_object_above_below_exp = 0.         # default = 0.; i.e. no grid search
#theta_multiply_factor = 2.                 # default = 2.
#theta_above_below_exp = 0.                 # default = 0.
#lambda_OTF_multiply_factor = 2.            # default = 2.
#lambda_OTF_above_below_exp = 0.            # default = 0.; i.e. no grid search
#grid_search_lambda_object_estimate_flag = False    # default = False

### 88888 ###
#band_limited_constraint = False            # default = False
#true_object_file = None                    # default = None    
#true_PSF_file = None                       # default = None

### 99999 ###
#debug = False             # default = False
