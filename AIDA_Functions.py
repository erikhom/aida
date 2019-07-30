################################################################################
#
#   File:       AIDA_Functions.py
#
#   Summary:    Set of functions used in the Adaptive Image Deconvolution 
#               Algorithm (AIDA) package (AIDAP)
#
#   Authors:    Erik F.Y. Hom (Sedat Lab, UCSF) with assistance from
#               Timothy K. Lee (Undergraduate Intern (UCB), Sedat Lab)
#
#   Other:      See 'AIDA_version.py' for date and version details
#               See 'LICENSE.txt' for use, license, and reference details
#
#   Modifications : Clement Chalumeau (SETI intern) 2016/04/06
#                Remove all priithon occurences 
################################################################################


# from Priithon.all import *
#Used from Priithon.all:
#  - N (numpy)
#  - U (useful): nd, max and mmms
#  - F (fftfuncs): rfft and ndArrD
import numpy as N
from Priithon_heritage import useful as U
from Priithon_heritage import fftfuncs as F
from PIL import Image #Used to open TIFF files

#from Priithon import fftw
try:
    from Priithon_heritage import fftw
except:
    pass
#from Priithon import ccg

# from Priithon_heritage import ccg

import os, string, time, types, gc
#import pyfits #Clement: removed
import astropy.io.fits as iofits
import AIDA_Version as Ver
import AIDA_Settings as Set
import AIDA_GFunctions as AGF
import AIDA_CostGradFunctions


#######  FUNCTION: AIDA_Setup  #######
####    [checked]
def AIDA_Setup():
    """
    Used in 'AIDA.py' to process, set-up, and check input variables supplied
    by the user
    
    Loads user specified variables from 'AIDA_Settings.py' and/or the command
    line and  checks their validity.  Also processes specified image and PSF 
    files in preparation for deconvolution, writing out mPSF and v files
    """ 
    
    SetupResultsDirectoryOutputLog()

    print Ver.header
    print "\t\t(code location: " + str(os.path.dirname(AGF.__file__)) + ")"
    print
    print "###"
    
    if not Set.debug:

        print '|\tAIDA_' + Set.timestamp + '.log'
    
    if Set.dataset_label != '':
    
        print '|\tDataset label: "' + str(Set.dataset_label) + '"'

    if Set.notes != '':
    
        print '|\tNotes: "' + str(Set.notes) + '"'
    
    print "###"

    ### EHom (20130611): edited this to use pyFFTW instead of numpy's default FFT routines
    ### http://hgomersall.github.io/pyFFTW/sphinx/tutorial.html
    ### This is a necessary workaround since Seb's SWIGed FFTW code is broken with recent FFTW upgrades
    # Test SWIG-based FFTW
    a_ = F.zeroArrD(8)
    a_[0]=1
    # Use workaround anyway
    # if any(F.rfft(a_).real != [ 1.,  1.,  1.,  1.,  1.]):
    if True:
        print 
        print "SWIGed FFTW appears broken - use workaround module"
        import fftw_workaround
        #import sys,
        #sys.modules['Priithon.fftw'] = fftw_workaround
### EHom (20130612): TEMP  conditional until we work out pyFFTW integration

        if (Set.fft_type == 'pyFFTW'):
            fftw.fft = fftw_workaround.pyfftw_fft
            fftw.ifft = fftw_workaround.pyfftw_ifft
            fftw.rfft = fftw_workaround.pyfftw_rfft
            fftw.irfft = fftw_workaround.pyfftw_irfft
# ###
        elif(Set.fft_type == 'numpy_fft'):
            fftw.fft = fftw_workaround.numpy_fft
            fftw.ifft = fftw_workaround.numpy_ifft
            fftw.rfft = fftw_workaround.numpy_rfft
            fftw.irfft = fftw_workaround.numpy_irfft

        if any(F.rfft(a_).real != [ 1.,  1.,  1.,  1.,  1.]):
            raise RuntimeError, "FFT still broken."

        else:
            print " ** using " + Set.fft_type + " routines **"
            print

    print
    print "AIDA Start Time:", 
    print time.asctime(time.localtime(Set.start_setuptime))
    print
    
    
    print "  * processing decon specifications..."
    ProcessDeconSpecs() #@
    ProcessCostFunctionTerms() #@
    ProcessSetImageFilename() #@
    ProcessSetPSFFilename() #@
    
    print "  * setting up temporary working arrays..."
    SetupShapeParamsTempArrays() #@

    print "  * processing PSF data..."
    ProcessPSFs() #@
    
    print "\n  * converting to OTFs and processing data..."
    CalculatePSFOTFstats()
    
    print "  * processing background and noise settings..."
    ProcessImageBackgroundNoiseSettings() #@
    
    print "  * setting up optimization parameters..."
    SetupOptimizationParams() #@
    SetupObjectDerivativeLaplacian() #@
    SetupHyperparameterSearch() #@
    
    gc.collect()
    
    print "  * ...done with setup "
    print
####


#######  FUNCTION: SetupResultsDirectoryOutputLog()  #######
####    [checked]
def SetupResultsDirectoryOutputLog():   

    Set.start_setuptime = time.time()
    Set.start_localtime = time.localtime()
    Set.timehms = time.strftime('%H%M%S', Set.start_localtime)
    Set.timestamp = time.strftime('%y%m%d_', Set.start_localtime) + Set.timehms + '-' + Set.decon_type
    ### EHom (20130605): added decon tupe as a suffix to timestamp; this helps with directory identification

    if Set.timestamp_directory_flag:        ## use timestamped directories

        if Set.dataset_label:
    
            Set.results_directory = \
                    os.path.abspath(Set.results_directory) + '/' + \
                    Set.dataset_label + '_' + Set.timestamp + '/'
        else:
    
            Set.results_directory = \
                    os.path.abspath(Set.results_directory) + '/' + \
                    Set.timestamp + '/'
    else:                                   ## do not timestamp directories
        
        if Set.dataset_label:
    
            Set.results_directory = \
                    os.path.abspath(Set.results_directory) + '/' + \
                    Set.dataset_label + '/'     
        else:
        
            Set.results_directory = \
                    os.path.abspath(Set.results_directory) + '/'

    if os.path.isdir(Set.results_directory) == False:
    
        os.makedirs(Set.results_directory, mode=0777)

    ## open output log buffer
    Set.output = Set.results_directory + 'AIDA_' + Set.timestamp + '.log'

    if Set.debug == False:

        # write prints only into log file
        import sys
        sys.stdout = open(Set.output, 'w', 0)
    else:
        import sys
        # from: http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python
        # from: http://mail.python.org/pipermail/python-list/2007-May/442737.html
        class Tee(object):
            def __init__(self, name, mode, bufsize=0):
                self.file = open(name, mode, bufsize)
                self.stdout = sys.stdout
                sys.stdout = self
            def __del__(self):
                sys.stdout = self.stdout
                self.file.close()
            def write(self, data):
                self.file.write(data)
                self.stdout.write(data)
        Set._tee = Tee(Set.output, 'w', 0)
        '''
        if hasattr(sys.stdout, "fileno"): # not the case inside PyShell
            #print >>sys.stderr, "DUP2"

            # write prints both into log file and straight on the screen
            # open our log file
            so = se = open(Set.output, 'w', 0)

            # re-open stdout without buffering
            sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

            # redirect stdout and stderr to the log file opened above
            os.dup2(so.fileno(), sys.stdout.fileno())
            os.dup2(se.fileno(), sys.stderr.fileno())
        #else:
        #    print >>sys.stderr, "NO DUP2"
        '''
####


####### FUNCTION: ProcessDeconSpecs  #######
####    [checked]
def ProcessDeconSpecs():   

    length = N.maximum(len(Set.object_PCG_iter_array), 
            len(Set.PSF_PCG_iter_array))

    if len(Set.object_PCG_iter_array) < length:
    
        diff = length - len(Set.object_PCG_iter_array)
        Set.object_PCG_iter_array = tuple(Set.object_PCG_iter_array) + diff*(0,)
    
    if len(Set.PSF_PCG_iter_array) < length:
    
        diff = length - len(Set.PSF_PCG_iter_array)
        Set.PSF_PCG_iter_array = tuple(Set.PSF_PCG_iter_array) + diff*(0,)

    if Set.decon_type == 'classical':

        Set.object_PCG_iter_array = Set.max_classical_PCG_iter
        Set.PSF_PCG_iter_array = length*(0,)
        Set.max_total_PCG_blocks = 1
        ####


#######  FUNCTION: ProcessCostFunctionTerms #######
####    [checked]
def ProcessCostFunctionTerms():   

    if Set.terms[1] == 1 and Set.terms[0] == 0:

        Set.terms[1] = 0
        print "\tN.B. Redundant terms. Terms changed to '" + str(Set.terms) + "'"
        
    elif N.sum(Set.terms[2:]) > 0 and Set.decon_type == 'classical':
    
        message = "\tN.B. 'terms'='" + str(Set.terms) + \
                "' are invalid for 'decon_type'='classical'!" + \
                "\n\tResetting 'terms' to (1,1,0,0)" 
        print message
        
        Set.terms = (1,1,0,0)
    elif Set.terms in ((0,0,1,0), (0,0,0,1)):
    
        if Set.decon_type == 'npsfs':

            print '\tN.B. Npsfs support for given terms not yet implemented.'
            print "\t'decon_type' changed to 'myopic'"

        elif Set.decon_type != 'myopic':

            print '\tN.B. Invalid decon_type for given terms!'
            print "\t'decon_type' changed to 'myopic'"

        Set.decon_type = 'myopic'

        temp = tuple(N.floor(N.array(Set.max_classical_PCG_iter)))
        Set.PSF_PCG_iter_array= 10*temp
        Set.object_PCG_iter_array = len(Set.PSF_PCG_iter_array)*(0,)
        Set.max_total_PCG_blocks = 2
####


#######  FUNCTION: ProcessSetImageFilename  #######
####    
def ProcessSetImageFilename():   

    if type(Set.PSF_filenames) in (types.ListType, types.TupleType):
    
        ## for Nframes decon only
        if Set.decon_type not in ('nobjects', 'npsfs'):
        
            message = "\n'image_filenames' cannot be a 'list of a list' " + \
                    "for mono-frame deconvolution!\nPlease remove any " + \
                    "extraneous brackets/parentheses for this entry"
            raise ValueError, message

        imagefilelist = []
        Set.image_directory = []
        Set.imagefilebase = []
        
        for i in range(len(Set.image_filenames)):
            
            Set.image_directory.append(os.path.abspath(
                    os.path.dirname(Set.image_filenames[i])) + '/')
            Set.imagefilebase.append(os.path.basename(Set.image_filenames[i]))
            
            if i == 0:
                ## assuming all files are the same format
                Set.image_format = AGF.DetermineFileFormat(
                        directory=Set.image_directory[i], 
                        filenamebase=Set.imagefilebase[i]) #@ AGF

                if Set.output_format is None:
                
                    Set.output_format = Set.image_format



            if Set.image_format == 'f':

                temp_list = AGF.PickFiles(directory=Set.image_directory[i],
                        prefix=Set.imagefilebase[i], suffix='.fits') #@ AGF 
            elif Set.image_format == 'm':

                temp_list = AGF.PickFiles(directory=Set.image_directory[i],
                        prefix=Set.imagefilebase[i], suffix='.mrc') #@ AGF
            
            elif Set.image_format == 't':

                temp_list = AGF.PickFiles(directory=Set.image_directory[i],
                        prefix=Set.imagefilebase[i], suffix='.tiff') #@ AGF
                
            elif Set.image_format == 't2':

                temp_list = AGF.PickFiles(directory=Set.image_directory[i],
                        prefix=Set.imagefilebase[i], suffix='.tif') #@ AGF
            
            if len(temp_list) == 0:

                message = "\nNo image files were found that match:\n\t" + \
                    "imagefilebase = " + "'" + str(Set.imagefilebase[i]) + "'!"
                raise RuntimeError, message

            elif len(temp_list) == 1:
            
                ## test to make sure files have only a single image
                if Set.image_format == 'f':
    
                    temp = N.array(iofits.open(Set.image_directory[i] + \
                            temp_list[0])[0].data)
                elif Set.image_format == 'm':

                    temp = Mrc.bindFile('%s%s' %(Set.image_directory[i], 
                            temp_list[0]))
                    
                elif Set.image_format == 't':
                    temp = np.asarray(Image.open('%s%s' %(Set.image_directory[i], 
                            temp_list[0]))).astype(N.float64).copy()
                    
                elif Set.image_format == 't2':
                    temp = np.asarray(Image.open('%s%s' %(Set.image_directory[i], 
                            temp_list[0]))).astype(N.float64).copy()

                Set.shape = N.squeeze(temp).shape
                        
                if len(Set.shape) > Set.dimension:

                    if len(Set.image_filenames) > 1:
                    
                        message = "\nAIDA does not handle the " + \
                                "deconvolution of multiple files\n " + \
                                "containing multiple image frames this " + \
                                "time\nPlease run AIDA on each " + \
                                "multiple-framed file at a time"
                        raise RuntimeError, message

                else:
                                
                    imagefilelist.append(temp_list[0])

            elif len(temp_list) > 1:    ## use input common file prefix
            
                if len(Set.image_filenames) > 1:
            
                    message = "\nMultiple image files were detected " + \
                            "using: " + str(Set.image_filenames[i]) + \
                            " !\nOnly single file handles should be " + \
                            "provided \nin a 'list of list specification " + \
                            "of 'image_filename'\nN.B. Nframes decon " + \
                            "cannot be run with multiple multi-image files!"
                    raise ValueError, message

                else:
                
                    imagefilelist = temp_list[:]    ## assign file list
                    ## reset image_directory and imagefilebase
                    Set.image_directory = len(temp_list) * \
                            [Set.image_directory[i]]
                    Set.imagefilebase = imagefilelist[:]

        ## set-up data arrays
        ## assumes set of images in single decon run are ~the same 
        ## shape so that shape parameters, etc. may be set using the LAST
        ## file in imagefilelist
        
        if Set.image_format == 'f':
    
            temp = N.array(iofits.open(Set.image_directory[-1] + \
                    imagefilelist[-1])[0].data)
        elif Set.image_format == 'm':

            temp = Mrc.bindFile('%s%s' %(Set.image_directory[-1], 
                    imagefilelist[-1]))\
                    
        elif Set.image_format == 't':
            temp = np.asarray(Image.open('%s%s' %(Set.image_directory[-1], 
                    imagefilelist[-1]))).astype(N.float64).copy()
            
        elif Set.image_format == 't2':
            temp = np.asarray(Image.open('%s%s' %(Set.image_directory[-1], 
                    imagefilelist[-1]))).astype(N.float64).copy()

        if Set.resizeimage_flag:
            Set.temp_image = AGF.ResizeImage(image=temp, 
                    dimension=Set.dimension, zresize=Set.zresize_flag, 
                    dtype=N.float32) #@ AGF
        else:
        
            Set.temp_image = temp.copy()

        Set.shape = N.squeeze(Set.temp_image).shape
        Set.images_in_file = 'single'
    else:       ## 'Set.image_filenames' is a single string filename
                ## for mono-frame decon
        
            # KBU: case when image_finlenames is string or list *** add [0]
            # NO ! only when image_filenames is a list! (Clement Chalumeau)
        Set.image_directory = \
                os.path.abspath(os.path.dirname(Set.image_filenames)) + '/'
        Set.imagefilebase = os.path.basename(Set.image_filenames)


        Set.image_format = AGF.DetermineFileFormat(
                directory=Set.image_directory, filenamebase=Set.imagefilebase)
                #@ AGF

        if Set.output_format is None:

            Set.output_format = Set.image_format

        if Set.image_format == 'f':
            imagefilelist = AGF.PickFiles(directory=Set.image_directory,
                    prefix=Set.imagefilebase, suffix='.fits') #@ AGF 
        elif Set.image_format == 'm':
            imagefilelist = AGF.PickFiles(directory=Set.image_directory,
                    prefix=Set.imagefilebase, suffix='.mrc') #@ AGF '
        elif Set.image_format == 't':
            imagefilelist = AGF.PickFiles(directory=Set.image_directory,
                    prefix=Set.imagefilebase, suffix='.tiff') #@ AGF '
        elif Set.image_format == 't2':
            imagefilelist = AGF.PickFiles(directory=Set.image_directory,
                    prefix=Set.imagefilebase, suffix='.tif') #@ AGF '
            
            
        if len(imagefilelist) == 0:
            message = "\nNo image files were found that match:\n\t" + \
                "'Set.imagefilebase'=" + "'" + str(Set.imagefilebase) + "'!"
            raise RuntimeError, message

        ## set-up data arrays
        ## assumes set of images in one decon run are ~the same shape
        ## so that shape parameters, etc. may be set using the LAST file
        ## in imagefilelist
        if Set.image_format == 'f':
    
            temp = N.array(iofits.open(Set.image_directory + \
                    imagefilelist[-1])[0].data)
        elif Set.image_format == 'm':

            temp = Mrc.bindFile('%s%s' %(Set.image_directory, 
                    imagefilelist[-1]))
        
        else :#TIFF or TIF
            temp = N.asarray(Image.open('%s%s' %(Set.image_directory, 
                    imagefilelist[-1]))).astype(N.float64).copy()
        
        if Set.resizeimage_flag:

            Set.temp_image = AGF.ResizeImage(image=temp, 
                    dimension=Set.dimension, zresize=Set.zresize_flag, 
                    dtype=N.float32) #@ AGF
        else:
        
            Set.temp_image = temp.copy()

        Set.shape = N.squeeze(Set.temp_image).shape
        Set.images_in_file = 'single'   ## if necessary, this will be reset below

    ## create Set.image_list
    if Set.dimension == 3:
        
        if len(Set.shape) == 4: ## limited to multiple 3D images in a single file

            if len(imagefilelist) == 1:

                Set.images_in_file = 'multiple'
                filename = imagefilelist[0]
                Set.image_list = []

                for i in range(Set.shape[-4]):
                    Set.image_list.append(filename + '.%g'%i)

            elif len(imagefilelist) > 1:

                message = "\nAIDA does not handle the deconvolution of " + \
                        "multiple files\ncontaining multiple image frames " + \
                        "this time\nPlease run AIDA on each " + \
                        "multiple-framed file at a time"
                raise RuntimeError, message
                
            Set.shape = Set.shape[-3:]
        elif len(Set.shape) == 3:

            Set.image_list = imagefilelist
        else:

            Set.image_list = imagefilelist
            Set.dimension = 2
            print '\tN.B. Your so-called 3D data is really 2D!'
            print '\tSetting dimension to 2 instead'
    else:

        if len(Set.shape) > 2:

            if len(imagefilelist) == 1:

                Set.images_in_file = 'multiple'
                filename = imagefilelist[0]
                Set.image_list = []

                for i in range(Set.shape[-3]):

                    Set.image_list.append(filename + '.%g'%i)

            elif len(imagefilelist) > 1:

                message = "\nAIDA does not handle the deconvolution of " + \
                        "multiple files\ncontaining multiple image frames " + \
                        "this time\nPlease run AIDA on each multiple-framed" + \
                        " file at a time"
                raise RuntimeError, message

            Set.shape = Set.shape[-2:]
        else:

            Set.image_list = imagefilelist

    if Set.images_in_file == 'single' and len(imagefilelist) == 1:

        if Set.decon_type == ('nobjects', 'npsfs'):

            print "\tOnly one image file with a single image was provided"
            print "\tOne can't run", Set.decontype, "deconvolution on this"
            print "\tDefaulting to myopic deconvolution"
            print "\tPlease read the instructions on how to run AIDA!"

            Set.decon_type = 'myopic'

    ## check to make sure grid search is not paired with multiple images
    ## in a file
    if Set.images_in_file == 'multiple' and \
            (Set.lambda_object_above_below_exp != 0 or \
            Set.theta_above_below_exp != 0):
            
        message = "\nSorry, you cannot run a grid search on files that " + \
                "contain multiple images!"
        raise ValueError, message
        
    ###TEMP!!!
    if Set.true_object_file is not None:
    
        Set.true_object = AGF.LoadFile(Set.true_object_file)[0] #@ AGF
        print "\tTEMP!!!  TRUEOBJECT_FILE: ", Set.true_object_file

#   if Set.test_true_object_file is not None:
#   
#       Set.test_true_object = AGF.LoadFile(Set.test_true_object_file)[0] #@ AGF
#       print "\tTEMP!!!  TRUEOBJECT_FILE: ", Set.true_object_file
####


#######  FUNCTION: ProcessSetPSFFilename  #######
####
def ProcessSetPSFFilename():   
    
    if type(Set.PSF_filenames) in (types.ListType, types.TupleType):
        
        ## for Nframes decon only
        if Set.decon_type not in ('nobjects', 'npsfs'):

            message = "\n'PSF_filenames' cannot be a 'list of a list' " + \
                    "for mono-frame deconvolution!\nPlease remove any " + \
                    "extraneous brackets/parentheses for this entry"
            raise ValueError, message

        PSFfilelist = []
        Set.PSF_directory = []
        Set.PSFfilebase = []
        
        for i in range(len(Set.PSF_filenames)):

            Set.PSF_directory[i] = os.path.abspath( \
                    os.path.dirname(Set.PSF_filenames[i])) + '/'
            Set.PSFfilebase[i] = os.path.basename(Set.PSF_filenames[i])

            if i == 0:
                ## assuming all files are the same format
                Set.PSF_format = \
                        AGF.DetermineFileFormat(directory=Set.PSF_directory[i],
                        filenamebase=Set.PSFfilebase[i]) #@ AGF

            if Set.image_format == 'f':

                temp_list = AGF.PickFiles(directory=Set.PSF_directory[i],
                        prefix=Set.PSFfilebase[i], suffix='.fits') #@ AGF
            elif Set.image_format == 'm':

                temp_list = AGF.PickFiles(directory=Set.PSF_directory[i],
                        prefix=Set.PSFfilebase[i], suffix='.mrc') #@ AGF
    
            elif Set.image_format == 't':
                temp_list = AGF.PickFiles(directory=Set.PSF_directory[i],
                        prefix=Set.PSFfilebase[i], suffix='.tiff') #@ AGF
            elif Set.image_format == 't2':
                temp_list = AGF.PickFiles(directory=Set.PSF_directory[i],
                        prefix=Set.PSFfilebase[i], suffix='.tif') #@ AGF
    
    
            if len(temp_list) == 0:

                message = "\nNo PSF files were found that match:\n\t" + \
                    "PSFfilebase = " + "'" + str(Set.PSFfilebase[i]) + "'!"
                raise RuntimeError, message

            elif len(temp_list) > 1:
            
                message = "\nMultiple PSF files were detected using: " + \
                        str(Set.PSF_filenames[i]) + " !\nOnly single " + \
                        "file handles should be provided \nin a 'list of " + \
                        "list specification of 'PSF_filenames'!\nN.B. " + \
                        "Nframes decon cannot be run with multiple " + \
                        "multi-PSF files\nMake sure you do not have " + \
                        "extraneous brackets/parentheses"
                raise ValueError, message
            else:
            
                PSFfilelist.append(temp_list)
    else:       ## PSF_filenames is a string

        Set.PSF_directory = \
                os.path.abspath(os.path.dirname(Set.PSF_filenames)) + '/'
        Set.PSFfilebase = os.path.basename(Set.PSF_filenames)
        Set.PSF_format = AGF.DetermineFileFormat(directory=Set.PSF_directory,
                filenamebase=Set.PSFfilebase) #@ AGF
    
        if Set.PSF_format == 'f':

            PSFfilelist = AGF.PickFiles(directory=Set.PSF_directory,
                    prefix=Set.PSFfilebase, suffix='.fits') #@ AGF
        elif Set.image_format == 'm':

            PSFfilelist = AGF.PickFiles(directory=Set.PSF_directory,
                    prefix=Set.PSFfilebase, suffix='.mrc') #@ AGF
    
        else :
            PSFfilelist = AGF.PickFiles(directory=Set.PSF_directory,
                    prefix=Set.PSFfilebase, suffix='.tif') #@ AGF
    
        if len(PSFfilelist) == 0:

            message = "\nNo PSF files were found that match:\n\t" + \
                "'PSFfilebase' = " + "'" + str(Set.PSFfilebase) + "'!"
            raise RuntimeError, message
####


#######  FUNCTION: SetupShapeParamsTempArrays  #######
####    [checked]
def SetupShapeParamsTempArrays():   

    Set.realFFTshape = Set.shape[:-1] + (Set.shape[-1]/2 + 1,)
    Set.Nd = N.product(Set.shape)      ## Nd is nc in CCG
    Set.inv_Nd = 1./Set.Nd              ## used to scale inverse FFTW

    ## Temp working buffer arrays to be used throughout AIDA
    Set.temp = N.empty(shape=Set.shape, dtype=Set.dtype)
    Set.temp_single = N.empty(shape=Set.shape, dtype=N.float32)
    Set.temp_realFFT = N.empty(shape=Set.realFFTshape, dtype=Set.complex_type)
    Set.temp_realFFT_single = N.empty(shape=Set.realFFTshape, 
            dtype=N.complex64)
    Set.temp_realFFT_real = N.empty(shape=Set.realFFTshape, dtype=Set.dtype)
    Set.temp_realFFT_real_single = N.empty(shape=Set.realFFTshape, 
            dtype=N.float32)
####


#######  FUNCTION: ProcessPSFs  #######
####
def ProcessPSFs():   

    Set.start_PSF_processing_time = time.time()
    Set.PSF_processing_details = "PSF files processed"

    if N.sum(Set.cleaning):        ## check cleaning tuple specs
    
        Set.PSF_processing_details += " (and cleaned):"
    else:

        Set.PSF_processing_details += ":"

    ### Determine PSF file list
    if Set.PSF_format == 'f':
    
        PSFfilelist = AGF.PickFiles(Set.PSF_directory, Set.PSFfilebase, 
                '.fits') #@ AGF
        
    
    elif Set.PSF_format == 'm':

        PSFfilelist = AGF.PickFiles(Set.PSF_directory, Set.PSFfilebase, '.mrc')
                #@ AGF
    elif Set.PSF_format == 't' :
        PSFfilelist = AGF.PickFiles(Set.PSF_directory, Set.PSFfilebase, '.tiff')
    elif Set.PSF_format == 't2' :
        PSFfilelist = AGF.PickFiles(Set.PSF_directory, Set.PSFfilebase, '.tif')
        

    if len(PSFfilelist) == 0:

        message = "\nNo image files were found for:\n\t" + \
                "PSFfilebase = '" + str(Set.PSFfilebase) + "'!\n" + \
                "Please check that filebase and file extension are valid"
        raise RuntimeError, message
    
    ### Process PSFs ###
    Set.Nsample_PSFs = (len(PSFfilelist))   

    if len(PSFfilelist) == 1:   ## assume PSFs are all in one file

        (Set.PSFs, Set.hdr_PSF) = AGF.LoadFile(Set.PSF_directory + \
                PSFfilelist[0], dtype=N.float32) #@ AGF

        if Set.PSFs.ndim > Set.dimension:   ## multiple PSFs in one file
        
            Set.Nsample_PSFs = 1./Set.PSFs.shape[0]

            print "\t",
            
            for i in range(Set.PSFs.shape[0]):

                print "[%s]" %(i+1),
                
                Set.PSF_processing_details += "\n\t\t" + PSFfilelist[0] + \
                        " : #" + str(i)
                
                ## replace PSF with centered, cleaned, resized PSF
                Set.PSFs[i] = AGF.ProcessPSF(Set.PSFs[i], cropshape=Set.shape, 
                        center=Set.center, exclude_radius=Set.exclude_radius, 
                        clean=Set.cleaning, background_percent= \
                        Set.PSF_subtract_background_percent, 
                        nsigmas=Set.nsigmas,
                        threshold_percent=Set.PSF_threshold_percent, 
                        fill=Set.fill, dtype=N.float32) #@ AGF
            
        else:       ## only one PSF file
            
            Set.Nsample_PSFs = 1.
            Set.PSF_processing_details += "\n\t\t" + PSFfilelist[0]

            ## replace PSF with centered, cleaned, resized PSF
            Set.PSFs = AGF.ProcessPSF(Set.PSFs, cropshape=Set.shape, 
                    center=Set.center, exclude_radius=Set.exclude_radius, 
                    clean=Set.cleaning, background_percent= \
                    Set.PSF_subtract_background_percent, nsigmas=Set.nsigmas,
                    threshold_percent=Set.PSF_threshold_percent, 
                    fill=Set.fill, dtype=N.float32) #@ AGF

        if N.sum(Set.cleaning) and Set.output_intermediate_files_level > 1:

            if Set.origin_centered_PSF_flag:

                shift = None
            else:
            
                shift = N.array((0,)*(Set.PSFs.ndim-Set.dimension) + \
                        Set.shape)/2

            OutputPSFs(filebase=PSFfilelist[0], PSF=Set.PSFs, shift=shift)


    else:       ## PSFs are in different files
        #print "We're running here"
        Set.PSFs = N.empty(shape=((len(PSFfilelist),) + Set.shape), 
                dtype=N.float32)

        print "\t",
        
        for i in range(len(PSFfilelist)):
        
            print "[%s]" %(i+1),
            
            (temp, Set.hdr_PSF) = AGF.LoadFile(Set.PSF_directory + \
                    PSFfilelist[i], dtype=N.float32) #@ AGF

            PSF_shape = Set.PSFs.shape[-Set.dimension:]
            file_dim = len(Set.PSFs.shape)

            Set.PSF_processing_details += "\n\t\t" + PSFfilelist[i] + " : #" + \
                    str(i) #Edited by TBM May 13th 2013

            ## replace PSF with centered, cleaned, resized PSF
            Set.PSFs[i] = AGF.ProcessPSF(temp, cropshape=Set.shape, 
                    center=Set.center, exclude_radius=Set.exclude_radius, 
                    clean=Set.cleaning, background_percent= \
                    Set.PSF_subtract_background_percent, nsigmas=Set.nsigmas,
                    threshold_percent=Set.PSF_threshold_percent, 
                    fill=Set.fill, dtype=N.float32) #@ AGF

            if N.sum(Set.cleaning) and Set.output_intermediate_files_level > 1:

                if i == 0:      ## set shift only on first pass
                
                    if Set.origin_centered_PSF_flag:

                        shift = None
                    else:
            
                        shift = N.array(Set.shape)/2   ## PSF is Set.shape

                OutputPSFs(filebase=PSFfilelist[i], PSF=Set.PSFs[i], 
                        shift=shift)
####


#######  FUNCTION: Output_PSFs #######
####
def OutputPSFs(filebase, PSF, shift=None):

    cleanPSF_filename = GenerateCleanPSFfilename(filebase=filebase) #@
    Set.temp_hdr = PrepareOutputHeader(input_hdr=Set.hdr_PSF)

    if shift is not None:
        ## shift to array center and output
        AGF.Output2File(format=Set.output_format, data_array=(U.nd.shift(PSF, 
                shift=shift, order=3, mode="wrap")), filebase=cleanPSF_filename, 
                hdr=Set.temp_hdr, shape=None) #@ AGF
    else:

        AGF.Output2File(format=Set.output_format, data_array=PSF,
                filebase=cleanPSF_filename, hdr=Set.temp_hdr, shape=None) #@ AGF  
####


#######  FUNCTION: GenerateCleanPSFfilename  #######
####    
def GenerateCleanPSFfilename(filebase):

    if Set.origin_centered_PSF_flag:

        if Set.dataset_label:

            return '%s%s_cleanPSF_o_%s_%s' %( Set.results_directory, 
                    Set.dataset_label, Set.timehms, string.split(filebase, 
                    sep='.')[0])
        else:

            return '%scleanPSF_o_%s_%s' %( Set.results_directory, Set.timehms,
                    string.split(filebase, sep='.')[0])     
    else:

        if Set.dataset_label:

            return '%s%s_cleanPSF_%s_%s' %( Set.results_directory, 
                    Set.dataset_label, Set.timehms, string.split(filebase, 
                    sep='.')[0])
        else:

            return '%scleanPSF_%s_%s' %( Set.results_directory, Set.timehms,
                    string.split(filebase, sep='.')[0])
####


######  FUNCTION: CalculatePSFOTFstats  #######
####
def CalculatePSFOTFstats():
    """
    Returns (mPSF, inv_u, mOTF, inv_v)
    Based on ComputeMeanVariance single pass function
    """
    
    ## set-up output arrays
    Set.mPSF = N.zeros(shape=Set.shape, dtype=Set.dtype)
    Set.inv_u = N.zeros(shape=Set.shape, dtype=Set.dtype)
    Set.mOTF = N.zeros(shape=Set.realFFTshape, dtype=Set.complex_type)
    Set.inv_v = N.zeros(shape=Set.realFFTshape, dtype=Set.dtype)
    temp_realFFT = N.empty(shape=Set.realFFTshape, dtype=Set.complex_type)

    ## set-up u_floor and v_floor values
    Set.u_floor = (Set.u_floor_per_dim)**Set.dimension
    Set.v_floor = (Set.v_floor_per_dim)**Set.dimension
    
    ## compute statistics
    if Set.Nsample_PSFs == 1:
    
        Set.mPSF = Set.PSFs.copy()	# use .copy() instead of [:] = 
        
        if Set.initial_psf == 'mean':
	        del Set.PSFs		### EHom (20130611): commented out to allow use as initial PSFs guess
	        
        Set.inv_u = N.sqrt(Set.mPSF)   ## really u - temp assignment
        fftw.rfft(a=Set.mPSF, af=Set.mOTF, inplace=False)
        if Set.mOTF.flat[:100].real.max() > 1e100:
            message = "\n'Set.mOTF' has exploded - broken FFTW !?"
            raise RuntimeError, message             
        ## clean otf
        Set.mOTF = AGF.CleanRealOTF(realOTF=Set.mOTF, cutoff_radius=None, 
                OTF_threshold_percent=Set.OTF_threshold_percent, fill=Set.fill)
                #@ AGF
        Set.inv_v = N.sqrt(N.abs(Set.mOTF))   ## really v - temp assignment
        # EHom (20130712) TEMP 
        # print 'stats inv_v: '
        # print U.mmms(Set.inv_v)
        n = 2.
    else:
    
        n = 0.
 
        for psf in Set.PSFs:  # using the Knuth/Welford algorithm
    
            n += 1
            Set.temp = psf - Set.mPSF
            fftw.rfft(a=Set.temp, af=Set.temp_realFFT, inplace=False)
            ## clean otf difference
            Set.temp_realFFT = AGF.CleanRealOTF(realOTF=Set.temp_realFFT, 
                    cutoff_radius=None, OTF_threshold_percent= \
                    Set.OTF_threshold_percent, fill=Set.fill) #@ AGF
            Set.mPSF += Set.temp / n
            Set.mOTF += Set.temp_realFFT / n
            
            temp = (psf - Set.mPSF)     # compute new difference with new mean
            fftw.rfft(a=temp, af=temp_realFFT, inplace=False)
            temp_realFFT = AGF.CleanRealOTF(realOTF=temp_realFFT, 
                    cutoff_radius=None, OTF_threshold_percent= \
                    Set.OTF_threshold_percent, fill=Set.fill) #@ AGF
            Set.inv_u += Set.temp * temp    # temp inv_u - need to invert
            Set.inv_v += (Set.temp_realFFT * N.conjugate(temp_realFFT)).real

            # EHom (20130712) TEMP 
            # print 'stats inv_v: '
            # print U.mmms(Set.inv_v)

        if Set.initial_psf == 'mean':
            del Set.PSFs		### EHom (20130611): commented out to allow use as initial PSFs guess

        del temp, temp_realFFT

        ## make sure mOTF is properly normalized
        # Clement: problem with AGF.FitOTZero
#         Set.mOTF.flat[0] = AGF.FitOTFzero(OTF=Set.mOTF, npoints=20)
        Set.mOTF.flat[0] = AGF.FitOTFzero(OTF=Set.mOTF, npoints=20)
        Set.mOTF /= Set.mOTF.flat[0] # scale mOTF to start at 1.0

        ## zero out values of mPSF and u (temp here as inv_u) below threshold
        AGF.ThresholdPSF_inplace(PSF=Set.mPSF, threshold_percent= \
                Set.PSF_threshold_percent, threshold_value=None, fill=0.) #@ AGF
        AGF.ThresholdPSF_inplace(PSF=Set.inv_u, threshold_percent= \
                Set.PSF_threshold_percent, threshold_value=None, fill=0.) #@ AGF
                
#       ## make sure mPSF is normalized properly first
#       Set.mPSF /= N.sum(Set.mPSF.flat)       
#       fftw.rfft(a=Set.mPSF, af=Set.mOTF, inplace=False)

#       ## make sure mOTF is normalized properly second and regenerate PSF
#       Set.mOTF.flat[0] = AGF.FitOTFzero(OTF=Set.mOTF, npoints=20)
#       Set.mOTF /= Set.mOTF.flat[0] # scale mOTF to start at 1.0
#       fftw.irfft(af=Set.mOTF*Set.inv_Nd, a=Set.mPSF, inplace=False)       

#       ## zero out values of mPSF and u (temp here as inv_u) below threshold
#       Set.mPSF = AGF.ThresholdPSF(PSF=Set.mPSF, threshold_percent= \
#               Set.PSF_threshold_percent, threshold_value=None, fill=0.) #@ AGF
#       Set.inv_u = AGF.ThresholdPSF(PSF=Set.inv_u, threshold_percent= \
#               Set.PSF_threshold_percent, threshold_value=None, fill=0.) #@ AGF
                
#       ## normalized again and regenerate PSF
#       fftw.rfft(a=Set.mPSF, af=Set.mOTF, inplace=False)
#       Set.mOTF.flat[0] = AGF.FitOTFzero(OTF=Set.mOTF, npoints=20)     
##      Set.mOTF -= N.abs(Set.mOTF).min()
#       Set.mOTF /= Set.mOTF.flat[0] # scale mOTF to start at 1.0
#       fftw.irfft(af=Set.mOTF*Set.inv_Nd, a=Set.mPSF, inplace=False)       

        ## radially average v (temp here as inv_v)
        if Set.dimensions_to_radially_average_v == 2:

            Set.inv_v = AGF.RadiallyAverage2D(array=Set.inv_v, FT=True, 
                    origin=None, wrap=(False,), subtract_background=False, 
                    fill=0.)[1].astype(N.float64) #@ AGF
        elif Set.dimensions_to_radially_average_v == 3:

            Set.inv_v = AGF.RadiallyAverage3D(array=Set.inv_v, FT=True, 
                    origin=None, wrap=(False,False,True), subtract_background=False, 
                    fill=0.)[1].astype(N.float64)  #@ AGF  

    ## clean-up v values below set threshold
    if Set.band_limited_constraint:

        # allow for using a simple band-limited constraint
        #average = Set.inv_v.mean()
        Set.inv_v = N.where(Set.inv_v > Set.v_threshold, 1.0, 0.)
    else:
    
        Set.inv_v = N.where(Set.inv_v > Set.v_threshold, Set.inv_v, Set.fill)
    
    ## output mPSF, u, and v
    ## NOTE: Set.inv_u and Set.inv_v are 'u' and 'v' at this point!
    Output_mPSF_u_v() #@

    ## invert for true inv_u and inv_v; set inverse to zero if u or v is smaller
    ## than u_floor or v_floor, respectively
    #seb N.Error.pushMode(dividebyzero="ignore")
    _lastErrSettings = N.seterr(divide="ignore")

    Set.inv_u = N.where(Set.inv_u > Set.u_floor, (n-1)/Set.inv_u, 
            1./Set.u_floor) 

    if Set.band_limited_constraint:

        # allow for using a simple band-limited constraint
        Set.inv_v = N.where(Set.inv_v > 0., 1., 1./Set.v_floor)

    else:
    
        Set.inv_v = N.where(Set.inv_v > Set.v_floor, (n-1)/Set.inv_v, 
                1./Set.v_floor) 

    #seb N.Error.popMode()      # turn warning back on
    N.seterr(**_lastErrSettings)
    ## below used to scale temporary derivative variable in OTF constraint
    ## calculation.  See 'AIDA_CostGradFunctions.py'
    Set.OTF_deriv_scale = (0.5 + 0.5*Set.inv_Nd)

    if Set.true_PSF_file is not None:
    
        Set.true_PSF = AGF.LoadFile(Set.true_PSF_file)[0].astype(Set.dtype) #@ AGF
#         U.nd.shift(Set.true_PSF.copy(), shift=tuple(
#                 N.array(Set.true_PSF.shape)/2), output=Set.true_PSF, order=3,
#                 mode="wrap")
        #Clement: shift workaround
        Set.true_PSF = U.shift2D_workaround(Set.true_PSF.copy(),tuple(N.array(Set.true_PSF.shape)/2))
                        
    Set.PSF_processing_time = (time.time() - Set.start_PSF_processing_time)
####


#######  FUNCTION: Output_mPSF_u_v  #######
####
def Output_mPSF_u_v():
    """
    Writes mPSF, u, and v to file
    """

    ## Note that at the point this function is called, Set.inv_u and Set.inv_v
    ## are really u and v, respectively; output at this point to prevent
    ## division later
    u = Set.inv_u
    v = Set.inv_v
    
    if Set.output_format == 'f':
            
        Set.temp_hdr = None     # change me!
    else:
            
        Set.temp_hdr = None

    if Set.output_intermediate_files_level:
    
        if Set.origin_centered_PSF_flag:

            if Set.dataset_label:
        
                mPSF_filename = '%s%s_mPSF_o_%s' %(Set.results_directory,
                        Set.dataset_label, Set.timehms)
                u_filename = '%s%s_u_o_%s' %(Set.results_directory,
                        Set.dataset_label, Set.timehms)
                v_filename = '%s%s_v_o_%s' %(Set.results_directory,
                        Set.dataset_label, Set.timehms)
            else:
        
                mPSF_filename = '%smPSF_o_%s' %(Set.results_directory,
                        Set.timehms)
                u_filename = '%su_o_%s' %(Set.results_directory, Set.timehms)
                v_filename = '%sv_o_%s' %(Set.results_directory, Set.timehms)

            AGF.Output2File(format=Set.output_format, data_array=Set.mPSF,
                    filebase=mPSF_filename, hdr=Set.temp_hdr, shape=None) #@ AGF 

            #seb N.Error.pushMode(dividebyzero="ignore")
            _lastErrSettings = N.seterr(divide="ignore")
            AGF.Output2File(format=Set.output_format, data_array= \
                    N.where(Set.inv_u > Set.u_floor, 
                    Set.inv_u.astype(N.float64), 0.),
                    filebase=u_filename, hdr=Set.temp_hdr, shape=None) #@ AGF 
            AGF.Output2File(format=Set.output_format, data_array= \
                    AGF.realFT2fullFT(N.where(Set.inv_v > Set.v_floor, 
                    Set.inv_v.astype(N.float64), 0.)),
                    filebase=v_filename, hdr=Set.temp_hdr, shape=None) #@ AGF
            #seb N.Error.popMode()
            N.seterr(**_lastErrSettings)
        else:

            if Set.dataset_label:
        
                mPSF_filename = '%s%s_mPSF_%s' %(Set.results_directory,
                        Set.dataset_label, Set.timehms)
                u_filename = '%s%s_u_%s' %(Set.results_directory,
                        Set.dataset_label, Set.timehms)
                v_filename = '%s%s_v_%s' %(Set.results_directory,
                        Set.dataset_label, Set.timehms)
            else:
        
                mPSF_filename = '%smPSF_%s' %(Set.results_directory,
                        Set.timehms)
                u_filename = '%su_%s' %(Set.results_directory, Set.timehms)
                v_filename = '%sv_%s' %(Set.results_directory, Set.timehms)

            array_center = N.around(N.array(Set.shape)/2)
            #Clement: shift workaround
            Set.temp = U.shift2D_workaround(Set.mPSF.astype(N.float64),array_center)
#             U.nd.shift(Set.mPSF.astype(N.float64), shift=array_center, 
#                     output=Set.temp, order=3, mode="wrap")
            AGF.Output2File(format=Set.output_format, data_array= \
                    Set.temp, filebase=mPSF_filename, hdr=Set.temp_hdr, 
                    shape=None) #@ AGF  
            
            #seb N.Error.pushMode(dividebyzero="ignore")
            _lastErrSettings = N.seterr(divide="ignore")
            #Clement: shift workaround
            Set.temp = U.shift2D_workaround(Set.inv_u.astype(N.float64),array_center)
#             U.nd.shift(Set.inv_u.astype(N.float64), shift=array_center, 
#                     output=Set.temp, order=3, mode="wrap")
            AGF.Output2File(format=Set.output_format, data_array= \
                    N.where(Set.temp > Set.u_floor, Set.temp, 0.), 
                    filebase=u_filename, hdr=Set.temp_hdr, shape=None) #@ AGF  

            shift = array_center #; shift[-1] = 0
            #Clement: shift workaround
            Set.temp = U.shift2D_workaround(AGF.realFT2fullFT(Set.inv_v.astype(N.float64)),shift)
#             U.nd.shift(AGF.realFT2fullFT(Set.inv_v.astype(N.float64)), 
#                     shift=shift, output=Set.temp, order=3, mode="wrap") #@ AGF
            AGF.Output2File(format=Set.output_format, data_array= \
                    N.where(Set.temp > Set.v_floor, Set.temp, 0.), 
                    filebase=v_filename, hdr=Set.temp_hdr, shape=None) #@ AGF
            #seb N.Error.popMode()
            N.seterr(**_lastErrSettings)
####


#######  FUNCTION: ProcessImageBackgroundNoiseSettings  #######
####
def ProcessImageBackgroundNoiseSettings():   

    ## Set.background
    if type(Set.background) in (types.ListType, types.TupleType):
    
        ## check whether list/tuple is appropriate
        if Set.decon_type not in ('npsfs', 'nobjects') and \
                Set.images_in_file != 'multiple':
        
            message = "\n'background' array ENTRY ('-a' option) supplied:\n" + \
                    str(Set.background) + "\nis a list/tuple and can not " + \
                    "be used for mono-frame deconvolution unless multiple " + \
                    "images are in a single file!"
            raise ValueError, message
            
        elif len(Set.background) != len(Set.image_list):

            message = "\n'background' array ENTRY ('-a' option) supplied:\n" + \
                    str(Set.background) + "\n" + "does not match the " + \
                    "length of 'Set.image_list': " + \
                    str(len(Set.image_list)) + "!"
            raise ValueError, message
            
        elif Set.images_in_file == 'multiple':
        
            Set.background_array = len(Set.image_list)*[Set.background]
        else:
        
            Set.background_array = Set.background

    ## Set.dark_image
    ## ...for Nframes decon or Nimages/file
    if type(Set.dark_image) in (types.ListType, types.TupleType):

        ## check whether list/tuple is appropriate
        if Set.decon_type not in ('npsfs', 'nobjects') and \
                Set.images_in_file != 'multiple':
        
            message = "\n'dark_image' array ENTRY ('-a' option) supplied:\n" + \
                    str(Set.dark_image) + "\nis a list/tuple and can not " + \
                    "be used for mono-frame deconvolution unless multiple " + \
                    "images are in a single file!"
            raise ValueError, message

        elif len(Set.dark_image) != len(Set.image_list):
    
            message = "\n'dark_image' array ENTRY ('-a' option) supplied:\n" + \
                    str(Set.dark_image) + "\n" + "does not match the " + \
                    "length of 'Set.image_list': " + \
                    str(len(Set.image_list)) + "!"
            raise ValueError, message

        elif Set.sigma_det is None:
            ## check each list/tuple entry and see that it exists and that the
            ## list length matches that of 'Set.image_list'
            Set.sigma_det_array = []
        
            for i in range(len(Set.dark_image)):

                if Set.dark_image[i] is None:
            
                    Set.sigma_det_array.append(Set.sigma_det)
                elif os.path.isfile(Set.dark_image[i]):

                    dark = AGF.LoadFile(filename=Set.dark_image, 
                            dtype=N.float32)[0] #@ AGF
                             
                    Set.sigma_det_array.append(U.mmms(dark.flat)[3])
                else:
            
                    message = "\n'dark_image' array file ('-a' option) " + \
                            "specified:\n" + str(Set.dark_image[i]) + \
                            " is not valid or does not exist!"
                    raise ValueError, message
                    
            Set.sigma_det = Set.sigma_det_array[0]

    ## next two conditionals for mono-frame or single image/file
    ## sigma_det is the same for each image in the dataset
    elif Set.dark_image is not None:
    
        # compute Set.sigma_det for the dataset from a single dark_image file
        dark = AGF.LoadFile(filename=Set.dark_image, dtype=N.float32)[0] #@ AGF  
        Set.sigma_det = U.mmms(dark.flat)[3]

    ## Set.sigma_det
    ## for Nframes decon or Nimages/file
    if type(Set.sigma_det) in (types.ListType, types.TupleType):
    
        ## check whether list/tuple is appropriate
        if Set.decon_type not in ('npsfs', 'nobjects') and \
                Set.images_in_file != 'multiple':
        
            message = "\n'sigma_det' array ENTRY ('-a' option) supplied:\n" + \
                    str(Set.sigma_det) + "\nis a list/tuple and can not " + \
                    "be used for mono-frame deconvolution unless multiple " + \
                    "images are in a single file!"
            raise ValueError, message
            
        elif len(Set.sigma_det) != len(Set.image_list):

            message = "\n'sigma_det' array ENTRY ('-a' option) supplied:\n" + \
                    str(Set.sigma_det) + "\n" + "does not match the " + \
                    "length of 'Set.image_list': " + \
                    str(len(Set.image_list)) + "!"
            raise ValueError, message
            
        else:
        
            Set.sigma_det_array = Set.sigma_det[:]
            Set.sigma_det = Set.sigma_det_array[0]
    else:
        # sigma_det is the same for each image in the dataset
        Set.sigma_det_array = len(Set.image_list)*[Set.sigma_det]
####


####### FUNCTION: SetupOptimizationParams  #######
####    [checked]
def SetupOptimizationParams():   

    ## AIDA algorithm set-up parameters (outer loops)
    Set.xmax_PSF = 10.*U.max(Set.mPSF)
    Set.xmax_object = 10.*U.max(Set.temp_image)
    Set.object_PCG_tolerance *= (1./Set.Nd) * (Set.temp_image).mean()
    Set.PSF_PCG_tolerance *= (1./Set.Nd) * (Set.mPSF).mean()
    del Set.temp_image
    
    ## this is created only for mono-frame decon; otherwise, these variables
    ## will serve as a reference to an array element in multi-frame decon
    if Set.decon_type not in ('npsfs', 'nobjects'):

        Set.object = N.empty(shape=Set.shape, dtype=Set.dtype)
        Set.OBJECT = N.empty(shape=Set.realFFTshape, dtype=Set.complex_type)
        Set.PSF = Set.mPSF.copy().astype(Set.dtype)
        Set.OTF = N.empty(shape=Set.realFFTshape, dtype=Set.complex_type)
#       del Set.mPSF        # use for terms[3]
    
    Set.object_gradient = N.empty(shape=Set.shape,dtype=Set.dtype)
    Set.Delta1 = N.empty(shape=Set.shape, dtype=Set.dtype)
    Set.Laplacian = N.empty(shape=Set.shape, dtype=Set.dtype)
    Set.J = N.zeros(shape=len(Set.terms),dtype=Set.dtype)
    Set.obj_diff = []
    
    ## Other useful minimization arrays
    Set.old_estimate = N.empty(shape=Set.shape, dtype=Set.dtype)
    Set.old_estimate_difference = N.empty(shape=Set.shape, dtype=Set.dtype)
    
    if Set.decon_type in ('classical', 'myopic'):

        if Set.images_in_file == 'multiple' and Set.memory_usage_level > 0: 

            StoreImageData()   
        else:

            Set.memory_usage_level = 0

            print
            print '     Note: No gain in computation time for higher ' + \
                    'memory usage levels\n           for mono-frame decon ' + \
                    'with images in separate files'
            print '           Memory usage level defaulted to zero\n'

    elif Set.decon_type in ('npsfs', 'nobjects', 'si') and \
            Set.memory_usage_level > 0:

        StoreImageData() #@

    ## CCG set-up parameters (inner loops)
    Set.xmin = Set.xmin_value   ## minimum possible solution values (scalar!)
    # array of constraint flags
    Set.ivec = Set.ivec_value*N.ones(shape=Set.Nd, dtype=N.uint8)
####


#######  FUNCTION: SetupObjectDerivativeLaplacian  #######
####    [checked]
def SetupObjectDerivativeLaplacian():   

    Set.kernel_grad = AGF.SelectKernel(dimension=Set.dimension,
            operator=Set.derivative_operator, shift=-1, rho=Set.rho,
            zeta=Set.zeta, dtype=Set.dtype) #@ AGF

    if Set.laplacian_operator == 0:

        Set.kernel_lapl = AGF.SelectKernel(dimension=Set.dimension,
                operator=Set.derivative_operator, shift=+1, rho=Set.rho,
                zeta=Set.zeta, dtype=Set.dtype) #@ AGF
                
    else:

        Set.kernel_lapl = AGF.SelectKernel(dimension=Set.dimension,
                operator='laplacian', shift=Set.laplacian_operator, rho=Set.rho,
                zeta=Set.zeta, dtype=Set.dtype) #@ AGF
####


#######  FUNCTION: SetupHyperparameterSearch  #######
####    [checked]
def SetupHyperparameterSearch():   

    ### create hyperparameter search grid ###
    Set.lambda_object_exp_grid = \
            N.arange(-Set.lambda_object_above_below_exp,
                     Set.lambda_object_above_below_exp+1)
    Set.theta_exp_grid = N.arange(-Set.theta_above_below_exp,
                                  Set.theta_above_below_exp+1)

    if not Set.grid_search_lambda_object_estimate_flag:

        ## lambda_object_center set here only if not centering on
        ## lambda_object estimate
        Set.lambda_object_center = Set.orig_lambda_object_center

            
    ## script set-up steps
    Set.cum_CG_time = N.zeros(shape=3, dtype=N.float32)
    Set.cum_prerun_time = 0.
    Set.cum_decon_time = 0.
    Set.cum_deconvolutions = 0
    Set.extra_time = 0.
    Set.cum_CG_itns = N.zeros(shape=3, dtype=N.int)
    Set.cum_CostFunction_itns = N.zeros(shape=3, dtype=N.int)

    ## set-up hyperparameter estimation flags ###
    if Set.lambda_object_input:

        Set.lambda_object = Set.lambda_object_input
    else:
    
        Set.lambda_object = Set.orig_lambda_object_center
        
    if Set.lambda_OTF_input:
    
        Set.lambda_OTF = Set.lambda_OTF_input
    else:
    
        Set.lambda_OTF = Set.lambda_OTF_scaling*Set.inv_Nd
        ## note that Set.lambda_OTF_scaling and Set.lambda_OTF input are
        ## assumed to be mutually exclusive settings
        ###!! need to build in grid search over lambda_OTF !!###

    if Set.lambda_PSF_input:
    
        Set.lambda_PSF = Set.lambda_PSF_input
    else:
    
        Set.lambda_PSF = Set.lambda_PSF_scaling * 1.0   ## default = 1
        ## note that Set.lambda_PSF_scaling and Set.lambda_PSF input are
        ## assumed to be mutually exclusive settings
####


####  ASK SID RE BELOW...WORK IN PROGRESS
def ProcessOTF_EM():
    ''' Calculates inv_v for the spring damping term in the decon equation
    '''
    #Num=N.abs(Set.OTF_param-Set.OTF_meanparam)**2
    Den=Set.OTF_param**2-Set.OTF_meanparam**2
    Set.inv_v=Den
####    


#######  FUNCTION: PrepareOutputHeader  #######
####
def PrepareOutputHeader(input_hdr):

    if input_hdr is None or input_hdr == []:
    
        fitsobj = iofits.HDUList()
        fitsobj.append(iofits.PrimaryHDU())
        hdr = fitsobj[0].header
    else:
    
        hdr = input_hdr.copy()

    s000 = " AIDA Processed: " + time.asctime(Set.start_localtime)
    s001 = " timestamp = " + Set.timestamp
    hdr.update('AIDA_000', s000)
    hdr.update('AIDA_001', s001)
    
    return hdr
####


#######  FUNCTION: StoreImageData  #######
####    [checked]
def StoreImageData():   
    """
    Creates arrays for IMAGE and associated image weights for memory 
    usage levels greater than zero.  The data stored in these arrays are 
    called by GetImageData()
    """

    ## create arrays for image data storage
    Set.image_array = N.empty(shape=(len(Set.image_list),) + Set.shape,
            dtype=Set.dtype)
    Set.inv_w_array = N.empty(shape=(len(Set.image_list),) + Set.shape,
            dtype=N.float32)
    Set.inv_theta_center_array = N.empty(shape=(len(Set.image_list),) + \
            Set.shape, dtype=Set.dtype)



    ## create wiener filtered initial guesses for image array
    if Set.initial_object_guess == 'wiener':

        Set.wiener_array = N.empty(shape=len(Set.image_list), dtype=N.float32)

    ## create OBJECT and OTF storage arrays if enough memory is specified
    if Set.memory_usage_level > 1:

        if Set.decon_type == 'nobjects':
        
            Set.OBJECT_array = N.empty(shape=(len(Set.image_list),) + \
                    Set.realFFTshape, dtype=Set.complex_type)
            Set.OTF_array = N.empty(shape=(1,) + Set.realFFTshape, 
                    dtype=Set.complex_type)
        elif Set.decon_type == 'npsfs':
        
            Set.OBJECT_array = N.empty(shape=(1,) + \
                    Set.realFFTshape, dtype=Set.complex_type)
            Set.OTF_array = N.empty(shape=(len(Set.image_list),) + \
                    Set.realFFTshape, dtype=Set.complex_type)

    ## calculate image data and populate created arrays
    print '-'*30
    print "  ", Set.images_in_file, " images per file have been provided..."
    if Set.images_in_file == 'multiple':

        filecheck = string.split(Set.image_list[0], sep='.')
        imagefilename = '%s.%s'%(filecheck[0], filecheck[1])
        
        if Set.decon_type in ('npsfs', 'nobjects'):

            ## assuming Nframes decon limited to 1 multiple image file per round    
            imagefile = AGF.LoadFile(Set.image_directory[0] + imagefilename,
                    dtype=Set.dtype)[0] #@ AGF 
        else:

            imagefile = AGF.LoadFile(Set.image_directory + imagefilename,
                    dtype=Set.dtype)[0] #@ AGF

        if Set.resizeimage_flag:

            image = AGF.ResizeImage(image=imagefile, dimension=Set.dimension,
                    zresize=Set.zresize_flag, dtype=Set.dtype) #@ AGF
        else:
        
            image = imagefile.copy()

        for i in range(image.shape[0]):

            (Set.image_array[i], Set.inv_w_array[i], Set.sigma_det_array[i],
                    Set.wiener) = AGF.CalculateImageData(image[i],
                    background=Set.background, sigma_det=Set.sigma_det,
                    wiener=Set.initial_object_guess, dtype=Set.dtype) #@ AGF

            if Set.initial_object_guess == 'wiener':

                Set.wiener_array[i] = Set.wiener
                
            # Set theta
            Set.inv_theta_center_array[i] = \
                    1. / (CalculateTheta(w=(1./Set.inv_w_array[i]),
                    sigma_det=Set.sigma_det_array[i]) ) #@
    else:       # images in single files

        for i in range(len(Set.image_list)):

            #KBU test
            #print Set.image_directory[i]
            #print Set.image_list[i]
            
            # KBU: case when image_directory is string or list
            if type(Set.image_directory)==type(str()):
                imagefile = AGF.LoadFile(Set.image_directory + Set.image_list[i],
                dtype=Set.dtype)[0] #@ AGF
            else:
                imagefile = AGF.LoadFile(Set.image_directory[i] + Set.image_list[i],
                dtype=Set.dtype)[0] #@ AGF
                    
            if Set.resizeimage_flag:

                image = AGF.ResizeImage(image=imagefile, 
                        dimension=Set.dimension, zresize=Set.zresize_flag,
                        dtype=Set.dtype) #@ AGF
            else:
            
                image = imagefile.copy()

            (Set.image_array[i], Set.inv_w_array[i], Set.sigma_det_array[i],
                    Set.wiener) = AGF.CalculateImageData(image,
                    background=Set.background, sigma_det=Set.sigma_det_array[i],
                    wiener=Set.initial_object_guess, dtype=Set.dtype) #@ AGF
                    
                
            if Set.initial_object_guess == 'wiener':
        
                Set.wiener_array[i] = Set.wiener    

            # Set theta
            Set.inv_theta_center_array[i] = \
                    1. / (CalculateTheta(w=(1./Set.inv_w_array[i]),
                    sigma_det=Set.sigma_det_array[i]) ) #@
####


####### FUNCTION: CalculateTheta  #######
####    [checked]
def CalculateTheta(w, sigma_det):   
    """
    Returns 'theta' given data fidelity weights 'w' and 'sigma_det'
    
    Uses:
        'Set.theta_input'
        'Set.theta_definition'
        'Set.dimension'
    """

    if Set.theta_input:
    
        theta = Set.theta_input
    else:

        theta = Set.theta_scaling * N.sqrt( w / sigma_det )
    
    return theta
####


#######  FUNCTION: CCG:init_global_array_pointers  #######
####    [this is the SWIG generated py code linking to the CCG code]


#######  FUNCTION: Start_Print  #######
####    [checked]
def Start_Print():   
    """ 
    Prints user supplied info specified in AIDA.py, consistent with
    the specified 'infolevel' flag
    
    Used in 'AIDArun.py'
    """
    
    if Set.info_level >= 1:

        print '<<< Key User Input Variables >>>'
        print 'Results directory =', Set.results_directory
        print 'Deconvolution type =', Set.decon_type
        
        if type(Set.image_directory) in (types.ListType, types.TupleType) and \
                len(Set.image_directory) > 1:
        
            print 'Image file base =',
            print Set.image_directory[0] + Set.imagefilebase[0]
            
            for i in range(1, len(Set.image_directory)):
            
                print '\t\t ', Set.image_directory[i] + Set.imagefilebase[i]
        else:
        
            print 'Image file base =', Set.image_directory + Set.imagefilebase

        if type(Set.PSF_directory) in (types.ListType, types.TupleType) and \
                len(Set.PSF_directory) > 1:

            print 'PSF file base =',
            print Set.PSF_directory[0] + Set.PSFfilebase[0]
            print 'Initial PSF guess ='
            print Set.initial_psf
            
            for i in range(1, len(Set.PSF_directory)):
            
                print '\t\t', Set.PSF_directory[i] + Set.PSFfilebase[i]
                
        else:
        
            print 'PSF file base =', Set.PSF_directory + Set.PSFfilebase

        print 'Spatial dimension of data =', Set.dimension
        
        if Set.dark_image:
        
            print 'Dark image file =', Set.dark_image

        elif Set.sigma_det:
        
            print 'Sigma detector value =', Set.sigma_det
            
        else:
        
            print 'Background to subtract =', Set.background
            
        if Set.lambda_object_input:
        
            print 'lambda_object input =', Set.lambda_object_input
            
        if Set.theta_input:
        
            print 'theta input =', Set.theta_input

        if Set.lambda_OTF_input:

            print 'lambda_OTF input =', Set.lambda_OTF_input

        if Set.lambda_PSF_input:

            print 'lambda_PSF input =', Set.lambda_PSF_input

        print 'Results output format =', Set.output_format
        print 'Computational precision =', Set.precision
        print '-   -   -   -   -   -   -   -   -'
        print 'Initial object guess =', Set.initial_object_guess,
        
        if Set.initial_object_guess == 'wiener':
        
            print '\tWiener multiply factor =', Set.wiener_multiply_factor
            
        else:
        
            print
        
        print 'Object background subtraction =', Set.background         

        if Set.sigma_det_array:

            print 'sigma_det =', Set.sigma_det_array, ' (input)'
        
        elif Set.sigma_det:

            print 'sigma_det = %g (input)' %Set.sigma_det
            
        elif Set.dark_image:

            print 'sigma_det = %g (from dark image)' %Set.sigma_det

#       if Set.object_PCG_tolerance == 'default':
#       
#           Set.object_PCG_tolerance = (1./Set.Nd)
#       elif type(Set.object_PCG_tolerance) is types.StringType:
#       
#           message = "\n'object_tolerance' input must be 'default' " + \
#                   "or a number!"
#           raise ValueError, message
            
        print 'Input object tolerance =', Set.original_object_PCG_tolerance

#       if Set.PSF_PCG_tolerance == 'default':
#       
#           Set.PSF_PCG_tolerance = \
#                   Set.object_PCG_tolerance * Set.object_PCG_tolerance
#       elif type(Set.PSF_PCG_tolerance) is types.StringType:
#       
#           message = "\n'PSF_tolerance' input must be 'default' " + \
#                   "or a number!"
#           raise ValueError, message
        
        print 'Input PSF tolerance =', Set.original_PSF_PCG_tolerance
        print 'object PCG iteration array =', Set.object_PCG_iter_array
        print 'PSF PCG iteration array =', Set.PSF_PCG_iter_array
        print 'Maximum uphill object moves allowed before breaking =', 
        print Set.max_uphill_object_PCG_steps
        print 'Maximum uphill PSF moves allowed before breaking =', 
        print Set.max_uphill_PSF_PCG_steps
        print 'Rising tolerance ratio between estimate steps =',
        print Set.rising_tol_ratio
        print 'CCG object minimization tolerance =', Set.object_CCG_tolerance
        print 'CCG PSF minimization tolerance =', Set.PSF_CCG_tolerance
        print 'CCG object iterations per PCG block =', Set.max_object_CCG_iter
        print 'CCG PSF iterations per PCG block =', Set.max_PSF_CCG_iter
        print 'Number of sequential AIDA stop alerts before breaking =',
        print Set.max_sequential_PCG_stops
        print 'Maximum PCG stop signals before skipping optimization =',
        print Set.max_optimization_stops 
        print 'Maximum uphill signals before skipping optimization =',
        print Set.max_rising_stops
        
        if len(Set.lambda_object_exp_grid) > 1:

            print 'lambda_object exponent array = ', Set.lambda_object_exp_grid
            print 'lambda_object multipier =',
            print Set.lambda_object_multiply_factor

        if len(Set.theta_exp_grid) > 1:

            print 'theta exponent array = ', Set.theta_exp_grid
            print 'theta multiplier =', Set.theta_multiply_factor

        if Set.lambda_object_input is None:
                
            print "'lambda_object' will be calculated as: ",
            print "1./ < sqrt(2 pi w(r))/theta - 1 >"
                
    if Set.info_level >= 2:

        print '---'
        print 'Infolevel =', Set.info_level
        print 'Decon terms =', Set.terms
        print 'Object derivative and Laplacian operator parameters:'
        print '\tgradient operator = ', Set.derivative_operator,
        print '\tlaplacian operator = ', Set.laplacian_operator,
        
        if Set.dimension == 3:
        
            print '\tzeta = ', Set.zeta
        else:
        
            print ''
            
        print '\tgradient norm type =', Set.norm,
        print '\tzero edges:', Set.zero_edges_flag
        print 'Variable lower bound for CCG =', Set.xmin_value
        print 'PSF processing parameters:'
        print '\texclude_radius =', Set.exclude_radius
        print '\tcleaning flags =', Set.cleaning
        print '\tPSF centering =', Set.PSF_centering
        print '\tsubtract PSF background percent =',
        print Set.PSF_subtract_background_percent
        print '\tclean nsigmas to subtract =', Set.nsigmas
        print '\tPSF percent threshold set to fill value =', 
        print Set.PSF_threshold_percent
        print '\tOTF percent threshold set to fill value =', 
        print Set.OTF_threshold_percent
        print '\tfill value =', Set.fill
        print '\tfloor for PSF constraint weight =', Set.u_floor
        print '\tfloor for OTF constraint weight =', Set.v_floor
        print '\tradially averaged OTF constraint weight (dimensions):',
        print Set.dimensions_to_radially_average_v
        print Set.PSF_processing_details
        print 'lambda_OTF Value =', Set.lambda_OTF
        print 'lambda_PSF Value =', Set.lambda_PSF

    if Set.info_level > 1:

        print '<<< -   -   -   -   -   -   - >>>'

    Set.start_time = time.time()
    Set.setup_time = (Set.start_time - Set.start_setuptime)
    
    print "\nSet-up time =", 
    AGF.PrintTime(Set.setup_time) #@ AGF
    print "\n\n======="
    print "======="
    print
    print "Deconvolution Start Time:",
    print time.asctime(time.localtime(Set.start_time))
    print
####


#######  FUNCTION: DeconSingleFrame  #######
####    [checked]
def DeconSingleFrame():   
    """
    Runs AIDA in 'classical' or 'myopic' mode
    
    'decon_type' specifies 'classical' or 'myopic' decon
    'infolevel' controls the amount of information printed
    'memusagelevel' controls the number of variables stored in memory
        
    DeconSingleFrame works on image files sequentially.
    Procedure is as follows:
    (1) Load image and prepare Image data/stats
    (2) Prepare arrays for results storage
    (3) Run PCGSingleframe and generate
        deconvolved object and PSF
    (4) Output object and PSf result to disk.   
    """

    ### For loop over files
    ## Prepare arrays to memory map multiple results to a single file
    ResultsArrayPrep()   

    if Set.terms[0] == 1:

        Set.obj_diff_array = N.empty(shape=( (len(Set.image_list),) + \
                (len(Set.lambda_object_exp_grid),) + \
                (len(Set.theta_exp_grid),) + \
                (N.maximum(Set.max_total_PCG_blocks, 
                N.maximum(N.sum(Set.object_PCG_iter_array),
                N.sum(Set.PSF_PCG_iter_array)) ),) ), dtype=Set.dtype)

    for i in range(len(Set.image_list)):

        if i > 0 and Set.images_in_file == 'single':
        
            ResultsArrayPrep() #@       ## re-initialize arrays

        if Set.memory_usage_level > 0:



            Set.image = Set.image_array[i]
            Set.inv_w = Set.inv_w_array[i]
            Set.sigma_det = Set.sigma_det_array[i]
            filename = Set.image_list[i]
            Set.inv_theta = Set.inv_theta_center_array[i]
            Set.inv_theta_center = Set.inv_theta_center_array[i]
            
            if Set.info_level >= 1:

                print 'vvvvvvv'
                print '[[ ', Set.image_list[i], ' ]]', '\tshape =', Set.shape

            if Set.info_level >= 2:

                if Set.initial_object_guess == 'wiener':

                    print 'wiener = %.6g' %Set.wiener, 

                print 'sigma_det = %.6g' %Set.sigma_det
                print 'w (mmms): %.6g (min), %.6g (max), %.6g (mean), ' + \
                        '%.6g (std)' %U.mmms(1./Set.inv_w)
                print '-------'
        else:
        
            ProcessImageData(index=i) #@

        ### Hyperparameter for loop search (over object hyperparameters)
        hypersearch_total_CG_itns = N.zeros(shape=3, dtype=N.int)
        hypersearch_total_CostFunction_itns = N.zeros(shape=3, dtype=N.int)
        start_hypersearch_time = time.time()
        results_index = [0,0]
        
        ### For loop over hyperparameter values (lambda_object and theta)
        for j in range(len(Set.lambda_object_exp_grid)):

            Set.lambda_object_exp = Set.lambda_object_exp_grid[j]

            if not Set.grid_search_lambda_object_estimate_flag:

                Set.lambda_object = (Set.lambda_object_center) * \
                        (Set.lambda_object_multiply_factor**( \
                        Set.lambda_object_exp))
                ## N.B. centering on estimated lambda_object is handled
                ## in PCG_Single_Frame directly
                
            if len(Set.lambda_object_exp_grid) > 1:
            
                print "... lambda_object grid_search exponent =", 
                print Set.lambda_object_exp, "of", 
                print Set.lambda_object_exp_grid, "..."
            
            for k in range(len(Set.theta_exp_grid)):

                Set.inv_theta = (Set.inv_theta_center) * \
                        (Set.theta_multiply_factor**(\
                        Set.theta_exp_grid[k]))
                #Set.mu = (Set.lambda_object*Set.inv_theta**2)
                ## handle Set.mu setting in PCG_SingleFrame below
                
                if len(Set.theta_exp_grid) > 1:
            
                    print "... theta grid_search exponent =", 
                    print Set.theta_object_exp, "of", 
                    print Set.theta_object_exp_grid, "..."

                PCG_SingleFrame(results_index) #@ ,<- first iteration
                
                if Set.images_in_file == 'single':      ## grid search permitted
                
                    if Set.terms[0] == 1:
                        
                        Set.object_results[j,k] = Set.object
                        
                    if Set.decon_type == 'myopic':
                
                        Set.PSF_results[j,k] = Set.PSF

                hypersearch_total_CG_itns += Set.decon_total_CG_itns
                hypersearch_total_CostFunction_itns += \
                        Set.decon_total_CostFunction_itns

                if Set.terms[0] == 1:
                    
                    Set.obj_diff_array[i,j,k,:len(Set.obj_diff)] = \
                            N.array(Set.obj_diff, dtype=N.float32)
                results_index[1] += 1

            results_index[1] = 0
            results_index[0]  += 1

        Set.cum_CG_itns += hypersearch_total_CG_itns
        Set.cum_CostFunction_itns += hypersearch_total_CostFunction_itns
        
        ### Arrange results for output to disk
        if Set.images_in_file == 'single':
        
            # output each image in separate file to match input
            ProcessOutput(file=i) #@
        else:       # 'multiple'
            
            if Set.terms[0] == 1:
                    
                Set.object_results[i] = Set.object

            if Set.decon_type == 'myopic':
            
                Set.PSF_results[i] = Set.PSF
        
        File_Decon_Print(hypersearch_total_CG_itns,
                hypersearch_total_CostFunction_itns, start_hypersearch_time) #@

    if Set.images_in_file == 'multiple':

        # write out multiple images in one file
        ProcessOutput(file='') #@
####


#######  FUNCTION: OutputFilePrep  #######
####    [checked]
def ResultsArrayPrep():   
    """
    Prepares storage arrays for object and PSF results that are later Set. 
    for file output.

    Variables in the 'Set' module modified by this function:


    'Set' module variables used as input:


    """
    
    if Set.decon_type in ('classical', 'myopic'):

        if Set.images_in_file == 'single':
        
            Set.results_shape = tuple(N.array( \
                    (len(Set.lambda_object_exp_grid),) + \
                    (len(Set.theta_exp_grid),) + Set.shape))
        else:   # multiple images in a single file - NO grid search allowed
        
            Set.results_shape = tuple((len(Set.image_list),) + Set.shape)

        if Set.terms[0] == 1:
            # save data in single precision
            Set.object_results = N.empty(shape=Set.results_shape,
                    dtype=N.float32)
        else:

            Set.object_results = None
        if Set.decon_type == 'myopic':

            Set.PSF_results = N.empty(shape=Set.results_shape, 
                    dtype=N.float32)
        else:

            Set.PSF_results = None
    elif Set.decon_type == 'nobjects':      # no grid search allowed
                                            # converted to N.float32 on output
        Set.results_shape = N.array((len(Set.image_list),) + Set.shape)
        Set.object_results = N.empty(shape=Set.results_shape, dtype=N.float64)
        Set.PSF_results = N.empty(shape=(1,)+Set.shape, dtype=N.float64)
    elif Set.decon_type == 'npsfs':         # no grid search allowed
                                            # converted to N.float32 on output
        Set.results_shape = N.array((len(Set.image_list),) + Set.shape)
        Set.PSF_results = N.empty(shape=Set.results_shape, dtype=N.float64)
        Set.object_results = N.empty(shape=(1,)+Set.shape, dtype=N.float64)
####


#######  FUNCTION: ProcessImageData  #######
####    [checked]
def ProcessImageData(index=0):   
    """
    This function provides Image data necessary to calculate the Cost 
    Function

    'IMAGE' is the Fourier Transform of the image
    'inv_theta' is a regularization factor used in the object prior term
    'sigma_det' is a measure of Gaussian noise present in the image
    'inv_w' is a regularization factor used in the data fidelity term
    'wiener' (optional) can be used to generate a wiener-filtered initial 
            guess for the object

    ProcessImageData either calculates these variables or alternatively
    references existing arrays containing these variables for higher
    memory usage levels.
    """

    ### Calculate image data and populate created arrays
    if Set.images_in_file == 'multiple':
    
        filecheck = string.split(Set.image_list[0], sep='.')
        imagefilename = '%s.%s'%(filecheck[0], filecheck[1])

        # Clement
        # imagefile = AGF.LoadFile(Set.image_directory + imagefilename,
        #         dtype=Set.dtype)[0] #@ AGF
        result = AGF.LoadFile(Set.image_directory + imagefilename,
                dtype=Set.dtype)#@ AGF
        imagefile = result[0]
        Set.imageHdr = result[1]


        if Set.resizeimage_flag:

            image = AGF.ResizeImage(image=imagefile[index], 
                    dimension=Set.dimension, zresize=Set.zresize_flag, 
                    dtype=Set.dtype) #@ AGF
        else:
        
            image = imagefile[index].copy()

        (Set.image, Set.inv_w, Set.sigma_det, Set.wiener) = \
                AGF.CalculateImageData(image, background=Set.background,
                sigma_det=Set.sigma_det, wiener=Set.initial_object_guess,
                dtype=Set.dtype) #@ AGF
            
        ### Set theta       
        Set.inv_theta_center = 1. / (CalculateTheta( w=(1./Set.inv_w), \
                sigma_det=Set.sigma_det) ) #@   
    else:       # images in single files
        # Clement
        # imagefile = AGF.LoadFile(Set.image_directory + Set.image_list[index],
        #         dtype=Set.dtype)[0] #@ AGF
        result = AGF.LoadFile(Set.image_directory + Set.image_list[index],
                dtype=Set.dtype) #@ AGF
        imagefile = result[0]
        Set.imageHdr = result[1]

        if Set.resizeimage_flag:

            image = AGF.ResizeImage(image=imagefile, dimension=Set.dimension,
                    zresize=Set.zresize_flag, dtype=Set.dtype) #@ AGF 
        else:
        
            image = imagefile.copy()
    
        (Set.image, Set.inv_w, Set.sigma_det, Set.wiener) = \
                AGF.CalculateImageData(image, background=Set.background, 
                sigma_det=Set.sigma_det, wiener=Set.initial_object_guess,
                dtype=Set.dtype) #@ AGF

        ### Set theta
        Set.inv_theta_center = 1. / (CalculateTheta( w=(1./Set.inv_w), \
                sigma_det=Set.sigma_det) ) #@ 

    if Set.info_level >= 1:

        print 'vvvvvvv'
        print '[[ ', Set.image_list[index], ' ]]', '\tshape =', Set.shape

    if Set.info_level >= 2:

        if Set.initial_object_guess == 'wiener':

            print 'wiener = %.6g' %Set.wiener, 

        print 'sigma_det = %.6g' %Set.sigma_det
        print 'w (mmms): %.6g (min), %.6g (max), %.6g (mean), %.6g (std)' \
                %U.mmms(1./Set.inv_w)
        print '-------'
####


#######  FUNCTION: PCG_SingleFrame  #######
####    [checked]
def PCG_SingleFrame(results_index):   
    """
    The core routine of AIDA in 'myopic' mode
    
    [ORDER OF AIDA Minimization]

        INITIAL PASS:
            (1) Fix PSF = invFT(mOTF) ==>> PSF(i=0)
            (2) Set object(j=0) = 0 (option: Weiner filtered image)
    
        SUBSEQUENT PASSES:
            (3) Using PSF(i), solve for object(j+1) using CCG
            (4) If |object(j+1) - object(j)| <= object_tolerance,
                    check_object=1
            (5) Using object(j+1), solve for PSF(i+1) using CCG
            (6) If |PSF(i+1) - PSF(i)| <= PSF_tolerance, check_PSF=1
            (7) If (check_object == 'stop' && check_PSF == 'stop'), STOP
            (8) Else, go to (1) with i==>>i+1 and j==>>j+1
                
    For Set.terms[0] == 0 a PSF is fit to a sampled PSF distribution
    Currently used for testing purposes
        
    'decon_type' specifes 'nobjects' or 'npsfs' decon
    'infolevel' controls the level of information printed
    'memusagelevel' specifies the number of variables
        stored to memory
    'object_results' and 'PSF_results' stores the deconvolved
        object and PSF for output
    'results_index' indexes the results files based on the
        iteration in the hyperparameter search
    """

    ### Reset/set some AID variables
    if Set.initial_object_guess == 'wiener': 
    
        IMAGE = N.empty(shape=Set.realFFTshape,dtype=Set.complex_type)
        fftw.rfft(a=Set.image,af=IMAGE,inplace=None)
        Set.object = AGF.WienerFilter(IMAGE=IMAGE, OTF=Set.mOTF,
                weight=Set.wiener*Set.wiener_multiply_factor) #@ AGF
    elif Set.initial_object_guess == 'image':

        Set.object = Set.image.copy()	# use .copy() instead of [:] = 
    elif Set.initial_object_guess == 'zero':

        Set.object[:] = 0.
    else:
    
        message = "\n'initial_object' input must be 'zero', 'wiener', " + \
                "or 'image'!"
        raise ValueError, message


### MAYBE INSERT CONDITIONAL FOR FIRST PASS
###vvv TEMP vvv###
    if Set.true_PSF_file is not None:

        print "TEST!!! SET PSF TO TRUEPSF"
        # option 1
        #Set.PSF = Set.true_PSF.copy()      # initial PSF with true_PSF
        # option 2
        #Set.PSF = Set.mPSF.copy()          # initialize PSF with mean PSF
        # option 3
        Set.PSF = N.zeros(shape=Set.shape, dtype=Set.dtype)
        # re-set PSF initial guess to zero

        Set.mPSF = Set.true_PSF.copy()      # re-set mean PSF to true_PSF; use .copy() instead of [:] = 
        fftw.rfft(a=Set.true_PSF, af=Set.mOTF,inplace=None)
    else:

        ### create/initialize 'Set.PSF'
        Set.PSF = Set.mPSF.copy()
        ## or...re-set PSF initial guess to zero
        #Set.PSF = N.zeros(shape=Set.shape, dtype=Set.dtype)
    
    ## Initialize 'Set.OTF' for object estimation using 'Set.mOTF'
    Set.OTF = Set.mOTF.copy()	# use .copy() instead of [:] = 

    
##TEMP
#   U.nd.shift(F.gaussianArr(shape=Set.mPSF.shape, sigma=4,
#           integralScale=None, peakVal=None, orig=None), 
#           shift=tuple(N.array(Set.mPSF.shape)/2), output=Set.PSF[:], 
#           order=3, mode="wrap")
#   print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#   print "USING GAUSSIAN INITIAL GUESS!!!)"
#
#   tempname = '%sinitPSF_%s' %(Set.results_directory, Set.timehms)
#   Output2File(format=Set.output_format, data_array=Set.PSF,
#               filebase=tempname, hdr=None, shape=None)   

    Set.obj_diff = []

    ### Reset/set some CG counters
    optimization_round = 1
    check_object = check_PSF = 'go'
    object_stops = PSF_stops = 0
    Set.decon_total_CG_itns = N.zeros(shape=3, dtype=N.int)
    Set.decon_total_CostFunction_itns = N.zeros(shape=3, dtype=N.int)
    Set.extra_time = 0.

    if Set.terms[0] == 0:       ## stop check_object if no data fidelity term

        check_object = 'stop'

    if Set.decon_type == 'classical':

        check_PSF = 'stop'

    ### Estimate 'lambda_object'; if grid search, set center
    if Set.lambda_object_input is None:

        (Set.lambda_object, Set.inv_theta_center, Set.mu) = \
                CalculateLambdaObject(Set.object, Set.PSF) #@

        print "Using pixel-based object regularization"

        if Set.grid_search_lambda_object_estimate_flag:
        
            Set.lambda_object_center = Set.lambda_object
            Set.lambda_object = (Set.lambda_object_center) * \
                    (Set.lambda_object_multiply_factor**( \
                    Set.lambda_object_exp))
                    
    ## handle setting Set.mu here (instead of in DeconSingleFrame)
    ## to accommodate grid centering on lambda_object estimate
    Set.mu = (Set.lambda_object*Set.inv_theta**2)

    ### Print out data for current image
    PCG_Print(index=results_index)

    ### AID Outer Loop ###
    start_decon_time = time.time()
    seq_object_stops = seq_PSF_stops = 0
    seq_object_risings = seq_PSF_risings = 0
    test_object_stops = test_PSF_stops = 0

    while (optimization_round <= Set.max_total_PCG_blocks):
    #this is eta in Fig 1B

        if Set.terms[0] == 1 and \
                seq_object_stops < Set.max_optimization_stops and \
                seq_object_risings < Set.max_rising_stops and \
                check_object == 'go':

            ### Object Minimization Loop ###
            (object_stops, object_risings, fn) = \
                    Minimize_Object(optimization_round) #@
                    
            Minimization_Print(fn) #@

            if object_stops > 0:
            
                if object_stops == Set.max_sequential_PCG_stops:

                    seq_object_stops += 1
                else:
                
                    test_object_stops = object_stops
            elif object_risings >= Set.max_uphill_object_PCG_steps:
            
                seq_object_risings += 1
                
                if test_object_stops > 0:
                
                    seq_object_stops += 1
                    check_object = 'stop'

            if seq_object_stops == Set.max_optimization_stops or \
                    seq_object_risings == Set.max_rising_stops:
            
                check_object = 'stop'

        ## stop while loop if BOTH object and PSF estimates converge to
        ## the specified tolerance or not improving (max stops or rising)
        if check_object == 'stop' and check_PSF == 'stop':

            Convergence_Print(check_object, seq_object_stops, 
                    seq_object_risings, object_stops, check_PSF, seq_PSF_stops,
                    seq_PSF_risings, PSF_stops) #@

            break
            
        if Set.decon_type == 'myopic' and \
                seq_PSF_stops < Set.max_optimization_stops and \
                seq_PSF_risings < Set.max_rising_stops and check_PSF == 'go':

            ### PSF Minimization Loop ###
            if Set.true_object_file is not None:
            
                Set.object = Set.true_object.copy()	# use .copy() instead of [:] = 
                print "TEST!!! SET OBJECT TO TRUEOBJECT"
    
            (PSF_stops, PSF_risings, fn) = Minimize_PSF(optimization_round) #@

            Minimization_Print(fn) #@

            if  PSF_stops > 0:
            
                if PSF_stops == Set.max_sequential_PCG_stops:

                    seq_PSF_stops += 1
                else:
                
                    test_PSF_stops = PSF_stops 
            elif PSF_risings == Set.max_uphill_PSF_PCG_steps:
            
                seq_PSF_risings += 1
                
                if test_PSF_stops > 0:
                
                    seq_PSF_stops += 1
                    check_PSF = 'stop'

            if seq_PSF_stops == Set.max_optimization_stops or \
                    seq_PSF_risings == Set.max_rising_stops:
            
                check_PSF = 'stop'

        ### Global Convergence Check ###
        if check_object == 'stop' and check_PSF == 'stop':

            Convergence_Print(check_object, seq_object_stops, 
                    seq_object_risings, object_stops, check_PSF, seq_PSF_stops,
                    seq_PSF_risings, PSF_stops) #@

            break

        optimization_round += 1

    ### Print out Decon info
    PCG_Decon_Print(results_index, fn, optimization_round, start_decon_time) #@
####


#######  FUNCTION: PCG_Print  #######
####    [checked]
def PCG_Print(index):   
    """
    Prints partial conjugate gradient results consistent with user-specified
    'infolevel' flag
    
    Used in 'AIDA_Functions.doPCG_Classical' and 'AIDA_Functions.doPCG_Npairs'
    """

    if Set.info_level >= 1:

        if Set.decon_type not in ('npairs', 'nobjects') and \
                len(Set.lambda_object_exp_grid) > 1 and \
                len(Set.theta_exp_grid) > 1:

            print '\t<<< results index = ', index, '>>>'

        print 'terms:', Set.terms

        try:
        
            print 'lambda_object: %.6g (min), %.6g (max), %.6g (mean), %.6g (std)' \
                    %U.mmms(Set.lambda_object),
        except:
        
            print 'lambda_object: %.6g' %Set.lambda_object,

        if Set.lambda_object_input is None:
        
            print ' [adaptive]'
        else:
        
            print ''
        
        try:
        
            print 'theta: %.6g (min), %.6g (max), %.6g (mean), %.6g (std)' \
                    %U.mmms(1./Set.inv_theta)
        except:
        
            print 'theta: %.6g' %(1./Set.inv_theta)
                    
        try:
        
            print 'mu: %.6g (min), %.6g (max), %.6g (mean), %.6g (std)' \
                    %tuple(U.mmms(Set.mu))
        except:
        
            print 'mu: %.6g' %(Set.mu)

        if Set.lambda_OTF_input is None:
            
            print 'lambda_OTF: %.6g  [adaptive]' %Set.lambda_OTF
        else:

            print 'lambda_OTF: %.6g' %Set.lambda_OTF
        
        if Set.lambda_PSF_input is None:
            
            print 'lambda_PSF: %.6g  [adaptive]' %Set.lambda_PSF
        else:
        
            print 'lambda_PSF: %.6g' %Set.lambda_PSF

        if Set.wiener:
        
            print 'wiener: %.6g' %Set.wiener

        print '-------'
####


#######  FUNCTION: CalculateLambdaObject  #######
####    [checked]
def CalculateLambdaObject(object, PSF):
    """
    Computes lambda_object using formula derived by E. Hom:
    
        lambda_object = 1./ < (sqrt(2*pi*w)/theta) - 1 >
        
    Works on a single image:PSF pair
    """

    lambda_object = Set.lambda_object_scaling / \
            (N.sqrt(2*N.pi*Set.sigma_det) - 1)    # note: scalar!

##!!! TEMP 060411
##  w = 1./Set.inv_w
##  theta = 1./Set.inv_theta
##  den = F.irfft(F.rfft(PSF)*F.rfft(theta))
##  lambda_object = Set.lambda_object_scaling / \
##          (N.sqrt(2*N.pi*w)/den - 1)
##!!!

    if lambda_object < 0.:
    
        message = "\n'lambda_object' is negative!"
        raise RuntimeError, message

    inv_theta = Set.inv_theta
    mu = lambda_object * inv_theta * inv_theta

#   print "Using pixel-based object regularization"
    
#   try:
#       
#       print "calculated theta value = (%.6g, %.6g, %.6g, %.6g)" \
#               %U.mmms(1./inv_theta)           
#   except:
#
#       print "calculated theta value = %.6g" %(1./inv_theta)
#
#   try:
#
#       print "calculated lambda_object value = (%.6g, %.6g, %.6g, %.6g)" \
#               %U.mmms(lambda_object)
#               
    try:
    
        if (lambda_object).min() < 0.:
            
            print "Warning! Some 'lambda_object' values are negative!"
            print "Resetting these to the theoretical lower limit!"
                
            temp = 1./(N.sqrt(N.pi*Set.sigma_det*N.sqrt(2*N.pi)) - 1)
            lambda_object = N.where(lambda_object < 0, temp, lambda_object)

    except:

        pass
#       print "calculated lambda_object value = %.6g" %lambda_object
        
    print '---'

    return (lambda_object, inv_theta, mu)
####


#######  FUNCTION: CCG:doCCG()  #######
####    [this is the SWIG generated py code linking to the CCG code]
        
####



#######  FUNCTION: Minimize_Object  #######
####    [checked]
def Minimize_Object(optimization_round):
    """
    Prepares variables and setings for input into CCG code for minimization 
    Checks returned CCG data for convergence.
    Used in PCGSingleFrame and PCGMultiFrame
    
    'check_PSF' flags when PSF difference is beneath tolerance
    'PSF_stops' is the number of consecutive PSF differences
            beneath tolerance
    """

    if Set.decon_type in ('classical', 'myopic', 'nobjects'):

        Set.costfunction_type = 1   ## minimize object  
    elif Set.decon_type == 'npsfs':

        Set.costfunction_type = 3   ## minimize object in presence of npsfs
    elif Set.decon_type in ('si','siclassical'):
        
        Set.costfunction_type = 4   ## SI related; NEED INTEGRATION with SI code
    
    ### Initialization of CCG Variables  ###
    Set.old_estimate = Set.object.copy()	# use .copy() instead of [:] = 
    Set.old_estimate_difference[:] = 0.
    itn = 0;  ifn = 0;  fn = 0.;  fmin = 0.
    Set.ivec[:] = 0;  df0 = 0.  ## df0, scalar float; input and output
    Nclsrch = 0;  istop = 100
    rising_test_count = 0
    old_test = 100000.  ## to determine if solution steps go "uphill"
    object_stops = 0
    
    ###  Minimization Loop Over Object  ###
    for i in range(Set.object_PCG_iter_array[optimization_round-1]):
	#this is PIo
        startCGtime = time.time()
        
        if i == 0:

            print '[obj] PCG iter:', i+1, 'of', 
            print Set.object_PCG_iter_array[optimization_round-1],
            print ' in optimization round ', 
            print optimization_round, ' out of ',
            print Set.max_total_PCG_blocks, '(max)'
        else:

            print '[obj] PCG iter:', i+1, ' of ', 
            print Set.object_PCG_iter_array[optimization_round-1]
        
         # Pour test on fait avec scipy opt.minimize
        import scipy.optimize as opt
        def function_cost(X):
            temp_grad = X.copy()
            (res,dummy) = AIDA_CostGradFunctions.CostFunction(X.copy(),temp_grad)
            return(res,temp_grad.flatten())
        myObject = Set.object.copy()
        # pour test, ne sert a rien en fait
        (res,grad_res) = function_cost(myObject)
        # Tentative d'optimisation
        
        resOptim2 = opt.minimize(function_cost,myObject.flatten(),jac = True,method='L-BFGS-B',
                                bounds=N.asarray([(Set.xmin,Set.xmax_object)]*len(myObject.flatten())),
                                options={'disp': False,'maxiter': Set.max_object_CCG_iter/2
                                         })
        
        itn = resOptim2.nit
        ifn = resOptim2.nfev
        istop = 0
        fn = resOptim2.fun
        df0 = 0
        Nclsrch = 0
        Set.object[:]= resOptim2.x.reshape(Set.object.shape)
        
        
# #seb        (itn, ifn, istop, fn, df0, Nclsrch) = CCG.doCCG(Set.object, Set.xmin,
# #seb                Set.xmax_object, Set.ivec, Set.max_object_CCG_iter, fmin, df0,
# #seb                Set.object_CCG_tolerance, Nclsrch) #@ CCG           
#         (istop, itn, ifn, Nclsrch, fn, df0) = \
#                 ccg.getsol( Set.object, Set.xmin,
#                 Set.xmax_object, Set.ivec, Set.max_object_CCG_iter, fmin, df0,
#                 Set.object_CCG_tolerance,
#                 AIDA_CostGradFunctions.CostFunction)
        
                            
       
                

        
        if Set.info_level >= 1:
        
            Set.cum_CG_time[0::2] += time.time() - startCGtime
            Set.decon_total_CG_itns[0::2] += itn
            Set.decon_total_CostFunction_itns[0::2] += ifn

            if Set.info_level >= 2:

                if i == 0 :

                    try:

                        print '<lambda_object>: %.7g' %Set.lambda_object.mean(),
                    except:
                    
                        print 'lambda_object: %.7g' %Set.lambda_object,

                    try:
                    
                        print '   <mu>: %.7g' %Set.mu.mean(),
                    except:
                    
                        print '   mu: %.7g' %Set.mu,

                    try:
                    
                        print '   <theta>: %.7g' %(1./Set.inv_theta).mean()
                    except:
                    
                        print '   theta: %.7g' %(1./Set.inv_theta)

                print '\tCG itns', itn, '   istop', istop, '   ifn', ifn, 
                print '   df0 %.6g' %df0, '   Nclsrch', Nclsrch,
                print '   CGtime %.6g' %(time.time() - startCGtime)

        Nclsrch = 0

        if optimization_round == 3:

            Set.start = False

        ###  Check Global Solution Convergence of Object Estimate  ###
        test = (N.abs(N.abs(Set.object - Set.old_estimate) - \
                Set.old_estimate_difference)).mean()
        Set.obj_diff.append(test)

        if test <= Set.object_PCG_tolerance:

            old_test = 0.
            object_stops += 1

            if Set.info_level >= 2:

                print '\t    obj diff: ', test, '\tcheck_object = stop'

            if object_stops == Set.max_sequential_PCG_stops:

                if Set.info_level >= 2:

                    print '\t\t*** max object_stops reached ***'

                break
        else:

            object_stops = 0
                    
            if itn < (Set.max_object_CCG_iter+1) and \
                    i > Set.object_PCG_iter_array[optimization_round-1]+1 and \
                    Set.decon_type != 'classical':

                if Set.info_level >= 2:

                    print '\t    obj diff: ', test, '\tcheck_object = go'
                    print '\t    >> max specified object iterations reached <<'

                break

            if i != 0 and test >= Set.rising_tol_ratio * old_test:

                rising_test_count += 1

                if rising_test_count == Set.max_uphill_object_PCG_steps:

                    if Set.info_level >= 2:

                        print '\t    obj diff: ', test, '\tcheck_object = go'
                        print '\t    >> max rising test count encountered <<'

                    break   ## break out if test > old_test occurs
                            ## more than max_uphill_steps, consecutively
                else:
                    
                    if Set.info_level >= 2:

                        print '\t    obj diff: ', test, '\tcheck_object = go'
            else:

                Set.old_estimate_difference = N.abs(Set.object - Set.old_estimate).copy()
                                              # use .copy() instead of [:] = 
                old_test = test

                if Set.info_level >= 2:

                    print '\t    obj diff: ', test, '\tcheck_object = go'

        ###  Swap Old Object With Current Object Estimate 'xo'  ###
        Set.old_estimate = Set.object.copy()	# use .copy() instead of [:] = 

    return (object_stops, rising_test_count, fn)
####



#######  FUNCTION: Minimization_Print  #######
####    [checked]
def Minimization_Print(fn=0, Nelements=1.):
    """
    Prints results about the object minimization step in AIDA
    
    Used in 'AIDA_Functions.???', 'AIDA_Functions.???d',...
    """

    if Set.info_level > 1 :

        startextratime = time.time()
        temp = N.array((fn,) + tuple(Set.J), dtype=Set.dtype)
        temp /= Nelements
        Set.J_final = temp.copy()

        print "DEBUG seb: TEMP is nan:", N.any(N.isnan(temp))
        print 'J = %.6g \tJn: %.6g \tJo: %.6g \tJhk: %.6g \tJhr: %.6g' \
                %(temp[0], temp[1], temp[2], temp[3], temp[4])

        if Set.terms[0] == 1:       ## if Jn exists, print out relative
                                    ## contributions of other terms
            temp /= temp[1]

            print ' \t\tJn: %.6g \t\tJo: %.6g \tJhk: %.6g \tJhr: %.6g' \
                    %(temp[1], temp[2], temp[3], temp[4])

        print '-------'

        Set.extra_time += time.time() - startextratime
####


#######  FUNCTION: Convergence_Print  #######
####    [checked]
def Convergence_Print(check_object, seq_object_stops, seq_object_risings,
        object_stops, check_PSF, seq_PSF_stops, seq_PSF_risings, PSF_stops):
    """
    Prints convergence status of AIDA consistent with 'infolevel' flag
    """
    
    if Set.decon_type == 'classical':

        if Set.info_level >= 2:

            print '\t\t>>> obj::', check_object, '\tseq_stops:', 
            print seq_object_stops, '  seq_rising:', seq_object_risings, 
            print '  stops:', object_stops, '<<<'
    elif Set.terms[0] == 0:

        if Set.info_level >= 2:

            print '\t\t>>> PSF::', check_PSF, '\tseq_stops:', seq_PSF_stops, 
            print '  seq_rising:', seq_PSF_risings, 
            print '  stops:', PSF_stops, '<<<'
    else:

        if Set.info_level >= 2:

            print '\t\t>>> obj::', check_object, '\tseq_stops:', 
            print seq_object_stops, '  seq_rising:', seq_object_risings, 
            print '  stops:', object_stops, '<<<'
            print '\t\t>>> PSF::', check_PSF, '\tseq_stops:', seq_PSF_stops, 
            print '  seq_rising:', seq_PSF_risings, '  stops:', PSF_stops, '<<<'
    
    if Set.info_level >= 2:
    
        print "-------"
####


#######  FUNCTION: Convergence_Print  #######
####    [checked]
def Multiframe_Convergence_Print(check_Nobjects, Nobject_stops, check_Npsfs, 
        Npsf_stops):   
    """
    Prints convergence status of AIDA consistent with 'infolevel' flag
    """
    
    if Set.info_level > 1:

        print '\t\t>>> obj::', check_Nobjects, '\tstops:', Nobject_stops, '<<<'
        print '\t\t>>> PSF::', check_Npsfs, '\tstops:', Npsf_stops, '<<<'
        print "-------"
####


#######  FUNCTION: Minimize_PSF  #######
####    [checked]
def Minimize_PSF(optimization_round):
    """
    Prepares variables and setings for input into CCG code for 
    minimization.  Checks returned CCG data for convergence
    Used in PCGSingleFrame and PCGMultiFrame
    
    'check_PSF' flags when PSF difference is beneath tolerance
    'PSF_stops' is the number of consecutive PSF differences
            beneath tolerance
    """

    if Set.decon_type in ('classical', 'myopic', 'npsfs'):

        Set.costfunction_type = 0  ## minimize PSF
    elif Set.decon_type == 'nobjects':

        Set.costfunction_type = 2  ## minimize PSF in the presence of nobjects
    elif Set.decon_type == 'si':
        
        Set.costfunction_type = 5   ## NEEDS INTEGRATION with SI code

#   if optimization_round > 1: # for optimization round 1, this is taken
#                              # care of by the initialization in function
#                              # CalculatePSTOTFstats for PSF=mPSF
#                              # otherwise, normalization is assumed to be
#                              # okay for user inputted initial guess (e.g., 0)
#       ## first normalize OTF and regenerate PSF (do per optimization round)
#       fftw.rfft(a=Set.PSF, af=Set.OTF, inplace=False)
#       
#       print "HEY!", Set.OTF.flat[0]
#       
#       Set.OTF.flat[0] = AGF.FitOTFzero(OTF=Set.OTF, npoints=20)
#       Set.OTF -= N.abs(Set.OTF).min()
#       Set.OTF /= Set.OTF.flat[0] # scale OTF to start at 1.0
#       fftw.irfft(af=Set.OTF*Set.inv_Nd, a=Set.PSF, inplace=False)     
#       print "YO!", Set.OTF.flat[0]

    normalization = float(N.sum(Set.PSF.flat))

    if normalization > 0:
    
        Set.PSF /= normalization # renormalize everytime Minimize_PSF is called

    ###  Initialization of CCG Variables  ###
    Set.old_estimate = Set.PSF.copy()	# use .copy() instead of [:] = 
    Set.old_estimate_difference[:] = 0.
    itn = 0;  ifn = 0;  fn = 0.;  fmin = 0.
    Set.ivec[:] = 0;  df0 = 0.  # df0, scalar float; input and output
    Nclsrch = 0; istop = 100
    rising_test_count = 0
    old_test = 100000.      ## to determine if solution steps go "uphill"
    PSF_stops = 0

    ###  Minimization Loop Over PSF  ###
    for i in range(Set.PSF_PCG_iter_array[optimization_round-1]):
# PIh iterations here see Fig1B
        startCGtime = time.time()

        if i == 0:
        
            print '[PSF] PCG iter:', i+1, 'of', 
            print Set.PSF_PCG_iter_array[optimization_round-1],
            print ' in optimization round ', 
            print optimization_round, ' out of ', 
            print Set.max_total_PCG_blocks, '(max)'
        else:

            print '[PSF] PCG iter:', i+1, ' of ', 
            print Set.PSF_PCG_iter_array[optimization_round-1]

#seb        (itn, ifn, istop, fn, df0, Nclsrch) = CCG.doCCG(Set.PSF, Set.xmin,
#seb                Set.xmax_PSF, Set.ivec, Set.max_PSF_CCG_iter, fmin, df0,
#seb                Set.PSF_CCG_tolerance, Nclsrch) #@ CCG
#         (istop, itn, ifn, Nclsrch, fn, df0) = \
#                 ccg.getsol( Set.PSF, Set.xmin,
#                 Set.xmax_PSF, Set.ivec, Set.max_PSF_CCG_iter, fmin, df0,
#                 Set.PSF_CCG_tolerance,
#                 AIDA_CostGradFunctions.CostFunction)
                
        import scipy.optimize as opt

        
        def function_cost(X):
            temp_grad = X.copy()
            (res,dummy) = AIDA_CostGradFunctions.CostFunction(X.copy(),temp_grad)
            return(res,temp_grad.flatten())
        myPSF = Set.PSF.copy()
        
        if not hasattr(Set, 'evaluationTime'):
            Set.evaluationTime = 0
            Set.nb_Evaluations = 0
        import time as tt
        start_t = tt.time()
        res_osef = function_cost(myPSF.copy())
        end_t = tt.time()
        Set.evaluationTime = (Set.nb_Evaluations * Set.evaluationTime + 1000*(end_t - start_t))/(1+Set.nb_Evaluations)
        Set.nb_Evaluations += 1
        

        print "Duration of last evaluation: ",1000*(end_t - start_t),"ms"
        print "Mean: ",Set.evaluationTime,"ms"
        print "nb_evals: ", Set.nb_Evaluations
        
        resOptim2 = opt.minimize(function_cost,myPSF.flatten(),jac = True,method='L-BFGS-B',
                                bounds=N.asarray([(Set.xmin,Set.xmax_PSF)]*len(myPSF.flatten())),
                                options={'disp': False,'maxiter': Set.max_PSF_CCG_iter/2
                                         })
        
        itn = resOptim2.nit
        ifn = resOptim2.nfev
        istop = 0
        fn = resOptim2.fun
        df0 = 0
        Nclsrch = 0
        Set.PSF[:]= resOptim2.x.reshape(Set.PSF.shape)

        if Set.info_level >= 1:

            Set.cum_CG_time[1:] += time.time() - startCGtime
            Set.decon_total_CG_itns[1:] += itn
            Set.decon_total_CostFunction_itns[1:] += ifn

            if Set.info_level >= 2:

                if i == 0 :
                    try:
                        print '<lambda_OTF>: %.7g' %Set.lambda_OTF.mean(),
                    except:
                    
                        print 'lambda_OTF: %.7g' %Set.lambda_OTF,

                    try:
                        print '   <lambda_PSF>: %.7g' %Set.lambda_PSF.mean()
                    except:
                    
                        print '   lambda_PSF: %.7g' %Set.lambda_PSF

                print '\tCG itns', itn, '   istop', istop, '   ifn',
                print ifn, '   df0 %.6g' %df0, '   Nclsrch', Nclsrch,
                print '   CGtime %.6g' %(time.time() - startCGtime)

        Nclsrch = 0

        ###  Check Global Solution Convergence of PSF Estimate  ###
        test = (N.abs(N.abs(Set.PSF - Set.old_estimate) - \
                Set.old_estimate_difference)).mean()

        if test <= Set.PSF_PCG_tolerance:

            old_test = 0.
            PSF_stops += 1

            if Set.info_level >= 2:

                print '\t    PSF diff: ', test, '\tcheck_PSF = stop'

            if PSF_stops == Set.max_sequential_PCG_stops:

                if Set.info_level >= 2:

                    print '\t\t*** max PSF_stops reached ***'

                break
        else:

            PSF_stops = 0
            
            if itn < (Set.max_PSF_CCG_iter+1) and \
                    i > Set.PSF_PCG_iter_array[optimization_round-1]+1 and \
                    Set.decon_type != 'classical':

                if Set.info_level >= 2:

                    print '\t    PSF diff: ', test, '\tcheck_PSF = go'
                    print '\t    >> max specified PSF iterations reached <<'

                break

            if i != 0 and test >= Set.rising_tol_ratio * old_test:

                rising_test_count += 1

                if rising_test_count == Set.max_uphill_PSF_PCG_steps:

                    if Set.info_level >= 2:
                        
                        print '\t    PSF diff: ', test, '\tcheck_PSF = go'
                        print '\t    >> max rising test count encountered <<'
                    
                    break   ## break out if test > old_test occurs
                            ## more than max_uphill_steps, consecutively
                else:

                    if Set.info_level >= 2:

                        print '\t    PSF diff: ', test, '\tcheck_PSF = go'
            else:

                Set.old_estimate_difference = N.abs(Set.PSF - Set.old_estimate).copy()
                                              # use .copy() instead of [:] = 
                old_test = test
                
                if Set.info_level >= 2:

                    print '\t    PSF diff: ', test, '\tcheck_PSF = go'

        ###  Swap Old PSF With Current PSF Estimate 'xo'  ###
        Set.old_estimate = Set.PSF.copy()	# use .copy() (deep copy) instead of [:] (view)
        fmin=0.; df0=0.
    
    return (PSF_stops, rising_test_count, fn)
####


#######  FUNCTION: PCG_Decon_Print  #######
####    [checked]
def PCG_Decon_Print(results_index, fn, optimization_round, start_decon_time):
    """
    Prints details about the AIDA deconvolution for a given dataset image
    
    Used in 'AIDA_Functions.???', 'AIDA_Functions.???', etc.
    """

    if Set.info_level >= 1:
    
        Set.cum_deconvolutions +=1
        decon_time = (time.time() - start_decon_time) - Set.extra_time
        Set.cum_decon_time += decon_time

    if Set.decon_type == 'nobjects':

        print
        print '======='
        print 'Image Decon List'

        for i in range(len(Set.image_list)):

            J = Set.Nframes_J[i].copy()

            if Set.info_level > 1:

                startextratime = time.time()

                print '\t%g. %s'%(i, Set.image_list[i]), '\tshape =', \
                        Set.shape
                print '\t  J = %.5g' %(N.sum(J))
                print '\t\tJn: %.5g\tJo: %.5g\tJhk: %.5g\tJhr: %.5g' \
                        %(J[0], J[1], J[2], J[3])
                        
                J /= J[0]

                print '\t\tJn: %.4g \t\t\tJo: %.4g \tJhk: %.4g \tJhr: %.4g' \
                        %(J[0], J[1], J[2], J[3])
                Set.extra_time += time.time() - startextratime

        print '\t---------------'

		### EHom (20130605): added "axis=0" below to sum correctly; noticed by 
		### SHaase in Dec. 2010 and again by TMilnes recently (see below)
        ave_J = N.sum(Set.Nframes_J, axis=0)/len(Set.image_list)   

        print '\tJ_avg = %.5g' %(N.sum(ave_J))
        print '\t\tJn_avg: %.5g\tJo_avg: %.5g\tJhk_avg: %.5g\tJhr_avg: %.5g' \
                %(ave_J[0], ave_J[1], ave_J[2], ave_J[3])

        ave_J /= ave_J[0]

        print '\t\tJn_avg: %.4g' %(ave_J[0])
        print '\t\tJo_avg: %.4g \tJhk_avg: %.4g \t\tJhr_avg: %.4g' \
                %(ave_J[1], ave_J[2], ave_J[3])

        print '======='
    elif Set.decon_type == 'npsfs':

        print
        print '======='
        print 'Image Decon List'

        for i in range(len(Set.image_list)):
        
            J = Set.Nframes_J[i].copy()

            if Set.info_level > 1:
                startextratime = time.time()

                print '\t%g. %s'%(i, Set.image_list[i]), '\tshape =', Set.shape
                print '\t  J = %.5g' %(N.sum(J))
                print '\t\tJn: %.5g\tJo: %.5g\tJhk: %.5g\tJhr: %.5g' \
                        %(J[0], J[1], J[2], J[3])

                J /= J[0]

                print '\t\tJn: %.4g \tJo: %.4g \tJhk: %.4g \tJhr: %.4g' \
                        %(J[0], J[1], J[2], J[3])

                Set.extra_time += time.time() - startextratime

        print '\t----------------'

		# Edited by TBM--needed to sum along 0th axis...
		### EHom (20130605): yes--see above for Nobjects too
        ave_J = N.sum(Set.Nframes_J, axis=0)/len(Set.image_list)  

        print '\t  J_avg = %.5g' %(N.sum(ave_J))
        print '\t\tJn_avg: %.5g\tJo_avg: %.5g\tJhk_avg: %.5g\tJhr_avg: %.5g' \
                %(ave_J[0], ave_J[1], ave_J[2], ave_J[3])

        ave_J /= ave_J[0]

        print '\t\tJn_avg: %.4g \tJo_avg: %.4g \tJhk_avg: %.4g \tJhr_avg: %.4g'\
                %(ave_J[0], ave_J[1], ave_J[2], ave_J[3])

        print '======='
    else:

        if Set.info_level > 1:

            startextratime = time.time()

            print 'J = %.5g\t\tJn: %.5g\tJo: %.5g\tJhk: %.5g \tJhr: %.5g' \
                    %(Set.J_final[0], Set.J_final[1], Set.J_final[2],
                    Set.J_final[3], Set.J_final[4])

            if Set.terms[0] == 1:               ## if no Jn term, print out
                Set.J_final /= Set.J_final[1]   ## J terms relative to Jn

                print '\t\t\tJn: %.4g \t\tJo: %.4g \tJhk: %.4g \tJhr: %.4g' \
                        %(Set.J_final[1], Set.J_final[2], Set.J_final[3],
                        Set.J_final[4])

            print '^^^^^^^'

            Set.extra_time += time.time() - startextratime
    
    if Set.info_level >= 1:

        try:

            print '\n<lambda_object>: %.6g' %Set.lambda_object.mean(),
        except:
        
            print '\nlambda_object: %.6g' %Set.lambda_object,
        
#       try:
#
#           print '   mu: %.6g' %Set.mu.mean(),
#       except:
#       
#           print '   mu: %.6g' %Set.mu,
#   
#       try:        
#       
#           print '   theta: %.6g' %(1./Set.inv_theta).mean(),
#       except:
#
#           print '   theta: %.6g' %(1./Set.inv_theta),
#
        try:

            print '   <lambda_OTF>: %.6g' %Set.lambda_OTF.mean(),
        except:
        
            print '   lambda_OTF: %.6g' %Set.lambda_OTF,

        try:

            print '   <lambda_PSF>: %.6g' %Set.lambda_PSF.mean()
        except:
        
            print '   lambda_PSF: %.6g' %Set.lambda_PSF
            
        print
        print '| total number of while-loop decon iterations: ', 
        print optimization_round - 1
        print '| total number of CG calls/iterations: ',
        print Set.decon_total_CG_itns, '(obj PSF tot)'
        print '| total number of Cost Function calculations: ',
        print Set.decon_total_CostFunction_itns, '(obj PSF tot)'
        print '| total time of deconvolution: ', 
        AGF.PrintTime(decon_time) #@ AGF
        print
        
        if len(Set.lambda_object_exp_grid) > 1 and len(Set.theta_exp_grid) > 1:

            print str(results_index)*7

        print
####


#######  FUNCTION: ProcessOutput  #######
####    [checked]
def ProcessOutput(object_results=None, PSF_results=None, file=''):
    """
    Takes AIDA decon object and PSF result arrays and prepares for output to
    file.

    'object_results' = deconvolved object, can be a single object or an
            array of objects)
    'PSF_results' = deconvolved PSF, can be a single PSF or an array of PSFs)
    'file' = indexing reference Set. to indicate the corresponding image file
            for which object and PSF results are being outputted.
            If '', multiple results are outputted into a single file 
            
    N.B. 'decon_type'=classical and 'file'=None outputs the mean PSF to file, 
    the only PSF Set. in the classical decon scheme.  This is done only once 
    per decon.
    
    [Variables in the 'Set' module modified by this function:
        Set.output_results_shape
        Set.images_in_file
        Set.output_array_shape
        Set.object_results_array
        Set.PSF_results_array
        
    'Set' module variables used as input:
        Set.decon_type
        Set.terms
        Set.images_in_file
        Set.output_results_shape
        Set.output_format
        Set.results_directory
        Set.image_list
        Set.output_array_shape
    ]
    """
    
    outputshape = tuple(N.compress(N.array(Set.results_shape) > 1,
            Set.results_shape))
    (objresultsfile, PSFresultsfile) = GenOutputFilenames(file) #@


    # Clement: output HDR
    if Set.imageHdr is not None:
        outputHdr = PrepareOutputHeader(Set.imageHdr)
    else:
        outputHdr = PrepareOutputHeader()


    ###  Single Frame Decon Output  ####
    ## outputs one file each loop over files    
    if (Set.decon_type in ('classical', 'myopic')) and \
            Set.images_in_file == 'single':

        ## output object
        if Set.terms[0] == 1:


            Set.object_results.shape = outputshape
            AGF.Output2File(format=Set.output_format, 
                    data_array=Set.object_results, filebase=objresultsfile, 
                    hdr=outputHdr, shape=outputshape) #@ AGF


        ## output PSF
        if Set.decon_type == 'myopic':
            
            if not Set.origin_centered_PSF_flag:

                array_center = N.around(N.array(Set.shape)/2)

                ## loop over possible grid search; see '_ResultsArrayPrep()'
                for j in range(Set.results_shape[0]):
                
                    for k in range(Set.results_shape[1]):
                        
                        #Clement: shift workaround
                        Set.PSF_results[j,k] = U.shift2D_workaround(Set.PSF_results[j,k].copy(),array_center)   

#                         U.nd.shift(Set.PSF_results[j,k].copy(),
#                                 shift=array_center, output=Set.PSF_results[j,k],
#                                 order=3, mode="wrap")



            AGF.Output2File(format=Set.output_format, 
                    data_array=Set.PSF_results, filebase=PSFresultsfile, 
                    hdr=None, shape=outputshape) #@ AGF
            ### FIX hdr!! above

    ### Creates one file with multiple images consistent with input type
    elif (Set.decon_type in ('classical', 'myopic')) and \
            Set.images_in_file == 'multiple':   

        ## output object
        if Set.terms[0] == 1:


            Set.object_results.shape = outputshape

            AGF.Output2File(format=Set.output_format, 
                    data_array=Set.object_results, filebase=objresultsfile, 
                    hdr=outputHdr, shape=outputshape) #@ AGF
            ### FIX hdr!! above

        ## output psf
        if not Set.origin_centered_PSF_flag:
        
            array_center = N.array((0,) + tuple(N.around(N.array( \
                    Set.shape)/2)[-Set.dimension:]) )
            #Clement: shift workaround
            Set.PSF_results = U.shift2D_workaround(Set.PSF_results.copy(),array_center)
#             U.nd.shift(Set.PSF_results.copy(), shift=array_center,
#                     output=Set.PSF_results, order=3, mode="wrap")

        AGF.Output2File(format=Set.output_format, data_array=Set.PSF_results,
                filebase=PSFresultsfile, hdr=None, shape=outputshape) #@ AGF 
        ### FIX hdr!! above

    ###  nobjects Decon Output  ###
    elif Set.decon_type == 'nobjects':

        ## output objects as N.float32
        if Set.images_in_file == 'single':

            ## loop over Set.image_list
            for i in range(outputshape[0]):

                # Clement : fix HDR!!!! only one for all the frames :(  
                AGF.Output2File(format=Set.output_format, 
                        data_array=Set.object_results[i].astype(N.float32),
                        filebase=objresultsfile[i], hdr=outputHdr, 
                        shape=outputshape) #@ AGF
                ### FIX hdr!! above
        elif Set.images_in_file == 'multiple':


            AGF.Output2File(format=Set.output_format, 
                    data_array=Set.object_results.astype(N.float32),
                    filebase=objresultsfile, hdr=outputHdr, shape=outputshape) #@ AGF
            ### FIX hdr!! above

        ## output single psf as N.float32
        if not Set.origin_centered_PSF_flag:

            array_center = N.array((0,) + tuple(N.around(N.array(
                    Set.shape)/2)[-Set.dimension:]) )
            #Clement: shift workaround
            Set.PSF_results = U.shift2D_workaround(Set.PSF_results.copy(),array_center)
#             U.nd.shift(Set.PSF_results.copy(), shift=array_center,
#                     output=Set.PSF_results, order=3, mode="wrap")
                    
        Set.PSF_results.shape = Set.shape
        AGF.Output2File(format=Set.output_format, 
                data_array=Set.PSF_results.astype(N.float32),
                filebase=PSFresultsfile, hdr=None, shape=outputshape) #@ AGF 
            ## FIX hdr!! above

    ### npsfs Decon Output ###
    elif Set.decon_type == 'npsfs':


        ## output single object as N.float32
        Set.object_results.shape = Set.shape
        AGF.Output2File(format=Set.output_format, 
                data_array=Set.object_results.astype(N.float32),
                filebase=objresultsfile, hdr=outputHdr, shape=None) #@ AGF
#                filebase=objresultsfile, hdr=None, shape=outputshape) #@ AGF
        ## FIX hdr!! above

        ## output psfs as N.float32
        if Set.images_in_file == 'single':

            if Set.origin_centered_PSF_flag:
                
                for i in range(outputshape[0]):



            
                    AGF.Output2File(format=Set.output_format, 
                            data_array=Set.PSF_results[i].astype(N.float32),
                            filebase=PSFresultsfile[i], hdr=None, 
                            shape=outputshape) #@ AGF
                    ### FIX hdr!! above
            else:

                for i in range(outputshape[0]):
            
                    array_center = N.around(N.array( \
                            Set.shape)/2)[-Set.dimension:]
                    #Clement: shift workaround
                    Set.PSF_results[i] = U.shift2D_workaround(Set.PSF_results[i].copy(),array_center)       
#                     U.nd.shift(Set.PSF_results[i].copy(), shift=array_center,
#                             output=Set.PSF_results[i], order=3, mode="wrap")
                    AGF.Output2File(format=Set.output_format, 
                            data_array=Set.PSF_results[i].astype(N.float32),
                            filebase=PSFresultsfile[i], hdr=None, 
                            shape=None) #@ AGF
#                            shape=outputshape) #@ AGF
                    ## FIX hdr!! above
        elif Set.images_in_file == 'multiple':

            if not Set.origin_centered_PSF_flag:
            
                array_center = N.array( (0,) + tuple(N.around(N.array( \
                        Set.shape)/2)[-Set.dimension:]) )
                #Clement: shift workaround
                Set.PSF_results = U.shift2D_workaround(Set.PSF_results.copy(),array_center)
#                 U.nd.shift(Set.PSF_results.copy(), shift=array_center,
#                         output=Set.PSF_results, order=3, mode="wrap")

            AGF.Output2File(format=Set.output_format, 
                    data_array=Set.PSF_results.astype(N.float32), 
                    filebase=PSFresultsfile, hdr=None, shape=outputshape) #@ AGF
            ## FIX hdr!! above
####


####  FUNCTION: GenOutputFilename  ####
####    [checked]
def GenOutputFilenames(file):
    """
    Generates output filename for the function '_ProcessOutput'
    """

    if len(Set.lambda_object_exp_grid) <= 1 and len(Set.theta_exp_grid) <= 1:

        if Set.lambda_object_input:

            lambdastr = AGF.Num2Str(Set.lambda_object) #@ AGF
            
            if Set.lambda_OTF_input:
            
                lambdastr += '_O' + AGF.Num2Str(Set.lambda_OTF) #@ AGF
                
            if Set.lambda_PSF_input:
            
                lambdastr += '_P' + AGF.Num2Str(Set.lambda_PSF) #@ AGF
        else:

            try:
    
                lambdastr = 'a' + AGF.Num2Str(Set.lambda_object.mean()) #@ AGF
            except:
            
                lambdastr = 'a' + AGF.Num2Str(Set.lambda_object) #@ AGF

            if Set.lambda_OTF_input:
            
                lambdastr += '_O' + AGF.Num2Str(Set.lambda_OTF) #@ AGF

            if Set.lambda_PSF_input:
            
                lambdastr += '_P' + AGF.Num2Str(Set.lambda_PSF) #@ AGF

    elif len(Set.lambda_object_exp_grid) > 1 and \
            Set.grid_search_lambda_object_estimate_flag and \
            Set.lambda_object_input is None:

            lambdastr = 'aG' + AGF.Num2Str(Set.lambda_object_center) #@ AGF

            if Set.lambda_OTF_input:
            
                lambdastr += '_O' + AGF.Num2Str(Set.lambda_OTF) #@ AGF

            if Set.lambda_PSF_input:
            
                lambdastr += '_P' + AGF.Num2Str(Set.lambda_PSF) #@ AGF
    else:

        try:

            lambdastr = 'grid' + AGF.Num2Str(Set.lambda_object_center.mean()) 
                    #@ AGF
        except:
        
            lambdastr = 'grid' + AGF.Num2Str(Set.lambda_object_center) #@ AGF

        if Set.lambda_OTF_input:
            
            lambdastr += '_O' + AGF.Num2Str(Set.lambda_OTF) #@ AGF

        if Set.lambda_PSF_input:
            
            lambdastr += '_P' + AGF.Num2Str(Set.lambda_PSF) #@ AGF

    ###  Mono Frame Decon Output  ###
    if (Set.decon_type in ('classical', 'myopic')):

        if file == '' and Set.images_in_file == 'single':

            ## object
            if Set.dataset_label:
            
                objresultsfile = '%sRobj_%s_%s_%s_L%s' %(
                        Set.results_directory, Set.dataset_label, Set.timehms,
                        string.split(Set.image_list, sep='.')[0], lambdastr)
            else:
            
                objresultsfile = '%sRobj_%s_%s_L%s' %(Set.results_directory,
                        Set.timehms, string.split(Set.image_list, sep='.')[0], 
                        lambdastr)

            ## PSF
            if Set.decon_type != 'classical':
        
                if Set.origin_centered_PSF_flag:

                    if Set.dataset_label:
                
                        PSFresultsfile = '%sRpsf_o_%s_%s_%s' %(
                                Set.results_directory, Set.dataset_label, 
                                Set.timehms, string.split(Set.image_list, 
                                sep='.')[0])
                    else:

                        PSFresultsfile = '%sRpsf_o_%s_%s' %(
                                Set.results_directory, Set.timehms,
                                string.split(Set.image_list, sep='.')[0])
                else:

                    if Set.dataset_label:
            
                        PSFresultsfile = '%sRpsf_%s_%s_%s' %(
                                Set.results_directory, Set.dataset_label, 
                                Set.timehms, string.split(Set.image_list, 
                                sep='.')[0])
                    else:

                        PSFresultsfile = '%sRpsf_%s_%s' %( 
                                Set.results_directory, Set.timehms,
                                string.split(Set.image_list, sep='.')[0])
            else:
            
                PSFresultsfile = None
                
        elif file == '' and Set.images_in_file == 'multiple':
        
            ## object
            if Set.dataset_label:
            
                objresultsfile = '%sRobj_%s_%s_%s_L%s' %(
                        Set.results_directory, Set.dataset_label, Set.timehms,
                        string.split(Set.image_list[0], sep='.')[0], lambdastr)
            else:
            
                objresultsfile = '%sRobj_%s_%s_L%s' %(Set.results_directory,
                        Set.timehms, string.split(Set.image_list[0], 
                        sep='.')[0], lambdastr)

            ## PSF
            if Set.decon_type != 'classical':
        
                if Set.origin_centered_PSF_flag:

                    if Set.dataset_label:
                
                        PSFresultsfile = '%sRpsf_o_%s_%s_%s' %(
                                Set.results_directory, Set.dataset_label, 
                                Set.timehms, string.split(Set.image_list[0], 
                                sep='.')[0])
                    else:

                        PSFresultsfile = '%sRpsf_o_%s_%s' %(
                                Set.results_directory, Set.timehms,
                                string.split(Set.image_list[0], sep='.')[0])
                else:

                    if Set.dataset_label:
            
                        PSFresultsfile = '%sRpsf_%s_%s_%s' %( \
                                Set.results_directory, Set.dataset_label, 
                                Set.timehms, string.split(Set.image_list[0],
                                sep='.')[0])
                    else:

                        PSFresultsfile = '%sRpsf_%s_%s' %( \
                                Set.results_directory, Set.timehms,
                                string.split(Set.image_list[0], sep='.')[0])
        else:
        
            ## object
            if Set.dataset_label:
            
                objresultsfile = '%sRobj_%s_%s_%s_L%s' %(
                        Set.results_directory, Set.dataset_label, Set.timehms,
                        string.split(Set.image_list[file], sep='.')[0], 
                        lambdastr)
            else:
            
                objresultsfile = '%sRobj_%s_%s_L%s' %(Set.results_directory,
                        Set.timehms, string.split(Set.image_list[file],
                        sep='.')[0], lambdastr)

            ## PSF
            if Set.decon_type != 'classical':
        
                if Set.origin_centered_PSF_flag:

                    if Set.dataset_label:
                
                        PSFresultsfile = '%sRpsf_o_%s_%s_%s' %(
                                Set.results_directory, Set.dataset_label,
                                Set.timehms, string.split(Set.image_list[file],
                                sep='.')[0])
                    else:

                        PSFresultsfile = '%sRpsf_o_%s_%s' %(
                                Set.results_directory, Set.timehms,
                                string.split(Set.image_list[file], sep='.')[0])
                else:

                    if Set.dataset_label:
            
                        PSFresultsfile = '%sRpsf_%s_%s_%s' %(
                                Set.results_directory, Set.dataset_label, 
                                Set.timehms, string.split(Set.image_list[file], 
                                sep='.')[0])
                    else:

                        PSFresultsfile = '%sRpsf_%s_%s' %(
                                Set.results_directory, Set.timehms,
                                string.split(Set.image_list[file], sep='.')[0])
            else:
            
                PSFresultsfile = None

    ###  nobjects Decon Output  ###
    elif Set.decon_type == 'nobjects':

        # object
        if Set.images_in_file == 'single':

            objresultsfile = []
        
            for i in range(len(Set.image_list)):

                if Set.dataset_label:

                    objresultsfile.append('%sRobj_%s_%s_%s_L%s_%s' %(
                            Set.results_directory, Set.dataset_label, 
                            Set.timehms, string.split(Set.image_list[i], 
                            sep='.')[0], AGF.Num2Str(
                            Set.lambda_object_Narray[i]), i)) #@ AGF
                else:
                
                    objresultsfile.append('%sRobj_%s_%s_L%s_%s' %(
                            Set.results_directory, Set.timehms,
                            string.split(Set.image_list[i], sep='.')[0],
                            AGF.Num2Str(Set.lambda_object_Narray[i]), i)) #@ AGF

        elif Set.images_in_file == 'multiple':

            if Set.dataset_label:
            
                objresultsfile = '%sRobj_%s_%s_%s_%s' %(Set.results_directory,
                        Set.dataset_label, Set.timehms, Set.imagefilebase,
                        'Nobj')
            else:
            
                objresultsfile = '%sRobj_%s_%s_%s' %(Set.results_directory,
                        Set.timehms, Set.imagefilebase, 'Nobj')

        # PSF
        if type(Set.imagefilebase) in (types.ListType, types.TupleType) and \
                len(Set.imagefilebase) > 1:
        
            temp = string.split(Set.imagefilebase[0], sep='.')[0] + '---' + \
                    string.split(Set.imagefilebase[-1], sep='.')[0]
        else:
        
            temp = Set.imagefilebase

        if Set.origin_centered_PSF_flag:
        
            if Set.dataset_label:
            
                PSFresultsfile = '%sRpsf_o_%s_%s_%s_%s' %(
                        Set.results_directory, Set.dataset_label, Set.timehms,
                        temp, 'Nobj')
            else:
            
                PSFresultsfile = '%sRpsf_o_%s_%s_%s' %(Set.results_directory,
                        Set.timehms, temp, 'Nobj')
        else:

            if Set.dataset_label:
            
                PSFresultsfile = '%sRpsf_%s_%s_%s_%s' %(Set.results_directory,
                        Set.dataset_label, Set.timehms, temp, 'Nobj')
            else:
            
                PSFresultsfile = '%sRpsf_%s_%s_%s' %(Set.results_directory,
                        Set.timehms, temp, 'Nobj')

    ### npsfs Decon Output ###
    elif Set.decon_type == 'npsfs':

        # object
        if type(Set.imagefilebase) in (types.ListType, types.TupleType) and \
                len(Set.imagefilebase) > 1:
        
            temp = string.split(Set.imagefilebase[0], sep='.')[0] + '---' + \
                    string.split(Set.imagefilebase[-1], sep='.')[0]
        else:
        
            temp = Set.imagefilebase

        if Set.dataset_label:

            objresultsfile = '%sRobj_%s_%s_%s_%s' %(Set.results_directory,
                    Set.dataset_label, Set.timehms, temp, 'Npsf')
        else:
        
            objresultsfile = '%sRobj_%s_%s_%s' %(Set.results_directory,
                    Set.timehms, temp, 'Npsf')

        # PSF
        if Set.images_in_file == 'single':

            PSFresultsfile = []
        
            if Set.origin_centered_PSF_flag:
                
                for i in range(len(Set.image_list)):

                    if Set.dataset_label:
        
                        PSFresultsfile.append('%sRpsf_o_%s_%s_%s_%s' %(
                                Set.results_directory, Set.dataset_label,
                                Set.timehms, string.split(Set.image_list[i], 
                                sep='.')[0], i))
                    else:

                        PSFresultsfile.append('%sRpsf_o_%s_%s_%s' %(
                                Set.results_directory, Set.timehms,
                                string.split(Set.image_list[i], sep='.')[0], i))
            else:

                for i in range(len(Set.image_list)):
            
                    if Set.dataset_label:
        
                        PSFresultsfile.append('%sRpsf_%s_%s_%s_%s' %(
                                Set.results_directory, Set.dataset_label,
                                Set.timehms, string.split(Set.image_list[i],
                                sep='.')[0], i))
                    else:

                        PSFresultsfile.append('%sRpsf_%s_%s_%s' %(
                                Set.results_directory, Set.timehms,
                                string.split(Set.image_list[i], sep='.')[0], i))

        elif Set.images_in_file == 'multiple':

            if Set.origin_centered_PSF_flag:
            
                if Set.dataset_label:
                
                    PSFresultsfile = '%sRpsf_o_%s_%s_%s_%s' %(
                            Set.results_directory, Set.dataset_label,
                            Set.timehms, Set.imagefilebase, 'Npsf')
                else:
                
                    PSFresultsfile = '%sRpsf_o_%s_%s_%s' %(
                            Set.results_directory, Set.timehms,
                            Set.imagefilebase, 'Npsf')          
            else:
            
                if Set.dataset_label:
            
                    PSFresultsfile = '%sRpsf_%s_%s_%s_%s' %(
                            Set.results_directory, Set.dataset_label,
                            Set.timehms, Set.imagefilebase, 'Npsf')
                else:
                
                    PSFresultsfile = '%sRpsf_%s_%s_%s' %(
                            Set.results_directory, Set.timehms,
                            Set.imagefilebase, 'Npsf')

    return (objresultsfile, PSFresultsfile)
####

    
#######  FUNCTION: File_Decon_Print  #######
####    [checked]
def File_Decon_Print(hypersearch_total_CG_itns, 
        hypersearch_total_CostFunction_itns, start_hypersearch_time):
    """
    Prints details about the AIDA deconvolution for a given dataset
    
    Used in 'AIDA_Functions.DeconClassical', 'AIDA_Functions.DeconNpairs',
    'AIDA_Functions.DeconNobjects', and 'AIDA_Functions.DeconNpsfs'
    """
    
    if Set.decon_type == 'classical':

        if Set.info_level >= 1 and len(Set.lambda_object_exp_grid) > 1 and \
                len(Set.theta_exp_grid) > 1:

            print '*** *** *** *** *** *** ***'
            print '|| TOTAL CG ITERATIONS FOR HYPERPARAMETER SEARCH:',
            print hypersearch_total_CG_itns[0], '(object only)'
            print '|| TOTAL COST FUNCTION ITERATIONS FOR HYPERPARAMETER',
            print 'SEARCH:',
            print hypersearch_total_CostFunction_itns[0], '(object only)'
            print '|| TOTAL TIME FOR HYPERPARAMETER SEARCH: ',
            print time.time() - start_hypersearch_time, 'secs'
            print '*** *** *** *** *** *** ***'
            print           
    elif Set.decon_type == 'myopic' and Set.terms[0] == 0:

        if Set.info_level >= 1 and len(Set.lambda_object_exp_grid) > 1 and \
                len(Set.theta_exp_grid) > 1:

            print '*** *** *** *** *** *** ***'
            print '|| TOTAL CG ITERATIONS FOR HYPERPARAMETER SEARCH:',
            print hypersearch_total_CG_itns[1], '(PSF only)'
            print '|| TOTAL COST FUNCTION ITERATIONS FOR HYPERPARAMETER',
            print 'SEARCH:',
            print hypersearch_total_CostFunction_itns[1], '(PSF only)'
            print '|| TOTAL TIME FOR HYPERPARAMETER SEARCH: ',
            print time.time() - start_hypersearch_time, 'secs'
            print '*** *** *** *** *** *** ***'
            print
    else:

        if Set.info_level >= 1 and len(Set.lambda_object_exp_grid) > 1 and \
                len(Set.theta_exp_grid) > 1:

            print '*** *** *** *** *** *** ***'
            print '|| TOTAL CG ITERATIONS FOR HYPERPARAMETER SEARCH:',
            print hypersearch_total_CG_itns, '(object, PSF, total)'
            print '|| TOTAL COST FUNCTION ITERATIONS FOR HYPERPARAMETER',
            print 'SEARCH:',
            print hypersearch_total_CostFunction_itns, '(object, PSF, total)'
            print '|| TOTAL TIME FOR HYPERPARAMETER SEARCH: ',
            print time.time() - start_hypersearch_time, 'secs'
            print '*** *** *** *** *** *** ***'
            print
####


#######  FUNCTION: DeconMultiFrame  #######
####    [checked]
def DeconMultiFrame():
    """
    Runs AIDA in 'nobjects' or 'npsfs' mode.
    
    'decon_type' specifies nobjects or npsfs decon
    'infolevel' controls the amount of information printed
    'memusagelevel' controls the number of variables stored to memory   

    This function prepares arrays for results storage, receives decon results
    from PCGMultiframe and outputs those results to file
    
    If specified by the user, a hyperparameter search over lambda_object and 
    and theta is executed here.  Theta search is centered around an analytical 
    expresion for theta (using the Hom formula) and expanded using 
    inv_theta_multiplier
    """

    ## Grid search may only work for npsfs (single object).  Note that output
    ## files for grid (npsfs) will need to be altered.  If grid search over
    ## Nobj, problem size increases multiplicatively...
    ## AS OF NOW, ASSUMES NO GRID SEARCH POSSIBLE FOR NFRAMES DECON
    ## (see ResultsArrayPrep)

    ###  [[ Hyperparameter For-Loop Search ]]  ###
    ## Set/Reset Hyperparameter search counters/timers
    hypersearch_total_CG_itns = N.zeros(shape=3, dtype=N.int)
    hypersearch_total_CostFunction_itns = N.zeros(shape=3, dtype=N.int)
    start_hypersearch_time = time.time()

    ###  Initialize Hyperparameter Arrays for Multiple Objects
    Set.lambda_object_Narray = N.empty(shape=(len(Set.image_list),),
            dtype=Set.dtype)
    Set.mu_Narray = N.empty(shape=Set.inv_theta_center_array.shape, 
            dtype=Set.dtype)

    if Set.lambda_object_input:

        Set.lambda_object = Set.lambda_object_center        ## no grid search

    ### Prepare arrays for output to file
    ResultsArrayPrep() #@ 
    PCG_MultiFrame() #@
    hypersearch_total_CG_itns += Set.decon_total_CG_itns
    hypersearch_total_CostFunction_itns += Set.decon_total_CostFunction_itns
    Set.cum_CG_itns += hypersearch_total_CG_itns
    Set.cum_CostFunction_itns += hypersearch_total_CostFunction_itns
    
    ### Arrange results for output to file
    ProcessOutput() #@
    File_Decon_Print(hypersearch_total_CG_itns,
            hypersearch_total_CostFunction_itns, start_hypersearch_time) #@  
####


#######  FUNCTION: PCG_MultiFrame  #######
####    [checked]
def PCG_MultiFrame():
    """
    The core routine of AIDA in 'nobjects' or 'npsfs' mode
    
        [ORDER OF AIDA Minimization]

            INITIAL PASS:
                (1) Fix PSF = invFT(mOTF) ==>> PSF(i=0)
                (2) Set object(j=0) = 0 (option: Weiner filtered image)
    
            SUBSEQUENT PASSES:
                (3) Using PSF(i), solve for object(j+1) using CCG
                (4) If |object(j+1) - object(j)| <= object_tolerance,
                        check_object=1
                (5) Using object(j+1), solve for PSF(i+1) using CCG
                (6) If |PSF(i+1) - PSF(i)| <= PSF_tolerance, check_PSF=1
                (7) If (check_object == 'stop' && check_PSF == 'stop'), STOP
                (8) Else, go to (1) with i==>>i+1 and j==>>j+1
    
        'decon_type' specifes 'nobjects' or 'npsfs' decon
        'infolevel' controls the level of information printed
        'memusagelevel' specifies the number of variables
            stored to memory
        'object_results' and 'PSF_results' stores the deconvolved
            object and PSF for output
        'results_index' indexes the results files based on the
            iteration in the hyperparameter search
    """
    ### N.B.  requires enough memory for Set.object_results and Set.PSF_results 
    
    results_index = 0

    ###  Reset/Set Some CG Counters/Timers
    optimization_round = 1
    check_Nframes = ['go','go']
    Nframes_array_stops = [0,0]
    Set.Nframes_J = N.zeros(shape=(len(Set.image_list),) + (len(Set.terms),), 
            dtype=Set.dtype)
    Set.decon_total_CG_itns = N.zeros(shape=3, dtype=N.int)
    Set.decon_total_CostFunction_itns = N.zeros(shape=3, dtype=N.int)

    ### Set-up arrays and initialize
    Set.PSF_array = Set.PSF_results # synonyms (see ResultsArrayPrep)
    Set.object_array = Set.object_results # synonyms (see ResultsArrayPrep)

### Correction by EHom (20130605)
### 	Commented below TEMP code out.  TEMP code was certainly for testing purposes
### 	but left behind, has been replacing all mPSF estimates in Npsfs/Nobjects decon with a 
###		gaussian mPSF instead of one based on the inputted PSF data and calculated in the
###		function: CalculatePSFOTFstats() !
###
#	###TEMP 091706
#    Set.mPSF = F.gaussianArr(shape=(256, 256), sigma=2, integralScale=None, 
#            peakVal=None, orig=orig(0,0), wrap=1).astype(N.float64)	# added replaced orig=None with orig=orig(0,0) due to bad shifting errors
#	# Edited above line by TBM...if this is the solution it's because of the hardcoded variable change to match the size of my test images...this size should instead be set to reflect whatever the PSF or Image size is...
#    #Set.mPSF = F.gaussianArr(shape=(512, 512), sigma=2, integralScale=None, 
#    #        peakVal=None, orig=None, wrap=1).astype(N.float64)
#    fftw.rfft(a=Set.mPSF, af=Set.mOTF, inplace=False)
#	###
    if Set.initial_psf == 'mean':
        Set.PSF_array = Set.mPSF.copy() 	# set for each PSF_array slice - make deep copy

    elif Set.initial_psf == 'initial':
        Set.PSF_array = Set.PSFs.copy() 	# set for each PSF_array slice - make deep copy
        del Set.PSFs					# now can delete this since loaded into Set.PSF_array

	### EHom (20130611): don't think Set.OTF needs to be set here for MultiframeDecon--will use Set.OTF_array
	# Set.OTF is used for decons other than Npsfs, where Set.OTF_array is calculated and used instead, based on Set.PSF_array
	#Set.OTF = Set.mOTF.copy() 	#N.empty(shape=Set.realFFTshape,dtype=Set.complex_type)
    Set.OBJECT = N.empty(shape=Set.realFFTshape,dtype=Set.complex_type)

##TEMP
#           if i == 0:
#           
#               junk = F.gaussianArr(shape=Set.mPSF.shape, sigma=5+3*i,
#                       integralScale=None, peakVal=None, orig=None)
#               U.nd.shift(junk, shift=tuple(N.array(Set.mPSF.shape)/2), 
#                       output=Set.PSF_array[i], order=3, mode="wrap", cval)
#               print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#               print "USING GAUSSIAN INITIAL GUESS!!!)"
#
#           tempname = '%sinitPSF_%s%s' %(Set.results_directory, Set.timehms, i)
#           Output2File(format=Set.output_format, data_array=Set.PSF_array[i],
#                   filebase=tempname, hdr=None, shape=None)   

##TEMP
#       U.nd.shift(F.gaussianArr(shape=Set.mPSF.shape, sigma=5,
#               integralScale=None, peakVal=None, orig=None),
#               shift=tuple(N.array(Set.mPSF.shape)/2), 
#               output=Set.PSF_array[:], order=3, mode="wrap")
#       print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#       print "USING GAUSSIAN INITIAL GUESS!!!)"
#
#       tempname = '%sinitPSF_%s' %(Set.results_directory, Set.timehms)
#       Output2File(format=Set.output_format, data_array=Set.PSF_array[i],
#               filebase=tempname, hdr=None, shape=None)   


    ### Setup Variables Specific For npsfs or nobjects ###
    if Set.decon_type == 'npsfs':
    
        check_Nframes_array = (N.ones(shape=(1,1), dtype=N.int),) + \
                (N.ones(shape=len(Set.image_list), dtype=N.int),)
        PSF_min_type = 'npsfs'
        obj_min_type = 'Object'
        Nobj = 1
        Npsf = len(Set.image_list)

        ### Estimate hyperparameters from set of images
        for i in range(len(Set.image_list)):

            if Set.initial_object_guess in ('zero', 'wiener', 'image'):

                if Set.memory_usage_level > 0:
    
                    Set.image = Set.image_array[i]
                    Set.inv_w = Set.inv_w_array[i]
                    Set.sigma_det = Set.sigma_det_array[i]
                    filename = Set.image_list[i]
                    Set.inv_theta = Set.inv_theta_center = \
                            Set.inv_theta_center_array[i]
                else:
        
                    ProcessImageData(index=i)   

            ###  Reset/Set Some AIDA Variables
            if Set.initial_object_guess == 'wiener':
            
                IMAGE = N.empty(shape=Set.realFFTshape,dtype=Set.complex_type)
                fftw.rfft(a=Set.image,af=IMAGE,inplace=None)
                Set.object_array[0] = AGF.WienerFilter(IMAGE=IMAGE, 
                        OTF=Set.mOTF, weight=Set.wiener_array[i] * \
                        Set.wiener_multiply_factor) #@ AGF
            elif Set.initial_object_guess == 'image':
    
                Set.object_array[0] = Set.image
            elif Set.initial_object_guess == 'zero':
                
                Set.object_array[0] = 0.

            if Set.lambda_object_input is None:

                (Set.lambda_object_Narray[i], Set.inv_theta_center_array[i],
                        Set.mu_Narray[i]) = \
                        CalculateLambdaObject(Set.object_array[0],
                        Set.PSF_array[i]) #@

                if i == 0:
                
                    print "Using pixel-based object regularization"
            else:
            
                Set.lambda_object_Narray[i] = Set.lambda_object
                Set.mu_Narray[i] = (Set.lambda_object*Set.inv_theta**2)
    elif Set.decon_type == 'nobjects':

        check_Nframes_array = (N.ones(shape=len(Set.image_list),
                dtype=N.int),) + (N.ones(shape=(1,1),dtype=N.int),)
        PSF_min_type = 'PSF'
        obj_min_type = 'nobjects'
        Nobj = len(Set.image_list)
        Npsf = 1

        for i in range(len(Set.image_list)):

            if Set.initial_object_guess in ('zero', 'wiener', 'image'):

                if Set.memory_usage_level > 0:
        
                    Set.image = Set.image_array[i]
                    Set.inv_w = Set.inv_w_array[i]
                    Set.sigma_det = Set.sigma_det_array[i]
                    filename = Set.image_list[i]
                    Set.inv_theta = Set.inv_theta_center = \
                            Set.inv_theta_center_array[i]
                else:
        
                    ProcessImageData(index=i) #@
                
            ###  Reset/Set Some AIDA Variables
            if Set.initial_object_guess == 'wiener':
            
                IMAGE = N.empty(shape=Set.realFFTshape,dtype=Set.complex_type)
                fftw.rfft(a=Set.image,af=IMAGE,inplace=None)
                Set.object_array[i] = AGF.WienerFilter(IMAGE=IMAGE, 
                        OTF=Set.mOTF, weight=Set.wiener_array[i] * \
                        Set.wiener_multiply_factor) #@ AGF
            elif Set.initial_object_guess == 'image':
    
                Set.object_array[i] = Set.image
            elif Set.initial_object_guess == 'zero':
            
                Set.object_array[i]  = 0.
                
            if Set.lambda_object_input is None:

                (Set.lambda_object_Narray[i],
                        Set.inv_theta_center_array[i], Set.mu_Narray[i]) = \
                        CalculateLambdaObject(Set.object_array[i],
                        Set.PSF_array[0]) #@
                
                if i == 0:
                
                    print "Using pixel-based object regularization"
            else:
            
                Set.lambda_object_Narray[i] = Set.lambda_object
                Set.mu_Narray[i] = (Set.lambda_object*Set.inv_theta**2)

    ## list images and image data
    Imagelist_Print() #@

    print '---'
    
    ###  AIDA Outer Loop
    start_decon_time = time.time()
    
    check_object = Nobj*['go']
    check_PSF = Npsf*['go']
    object_stops = N.zeros(shape=(Nobj,))
    PSF_stops = N.zeros(shape=(Npsf,))
    seq_object_stops = N.zeros(shape=(Nobj,))
    seq_PSF_stops = N.zeros(shape=(Npsf,))
    object_risings = N.zeros(shape=(Nobj,))
    PSF_risings = N.zeros(shape=(Npsf,))
    seq_object_risings = N.zeros(shape=(Nobj,))
    seq_PSF_risings = N.zeros(shape=(Npsf,))
    test_object_stops = N.zeros(shape=(Nobj,))
    test_PSF_stops = N.zeros(shape=(Npsf,))

    while (optimization_round <= Set.max_total_PCG_blocks):
# this is the eta loop
        if check_Nframes[0] == 'go':

            ### Object Minimization Loop ###
            print
            print '[[[ %s Minimization ]]]' %(obj_min_type)
            print

            if Set.memory_usage_level > 1: ## lower memory usage not implemented yet

                for i in range(Npsf):
    
                    ## FFT here to save doing repeatedly in CostGradFunction
                    fftw.rfft(a=Set.PSF_array[i], af=Set.OTF_array[i], 
                            inplace=False)

                if Set.decon_type == 'nobjects': ## only 1 PSF permitted now    
                        
                    Set.OTF = Set.OTF_array[0]
            else: ## Set.OTF_array not stored

                if Set.decon_type == 'nobjects':
                
                    fftw.rfft(a=Set.PSF_array[0], af=Set.OTF, inplace=False)

            ## PSF fixed, loop over N objects
            for n in range(Nobj):       ## generally for nobjects decon
                    
                if seq_object_stops[n] < Set.max_optimization_stops and \
                        seq_object_risings[n] < Set.max_rising_stops and \
                        check_object[n] == 'go':
        
                    Set.image = Set.image_array[n]
                    Set.inv_w = Set.inv_w_array[n].astype(Set.dtype)
                    Set.inv_theta = Set.inv_theta_center_array[n]
                    Set.lamba_object = Set.lambda_object_Narray[n]
                    Set.mu = Set.mu_Narray[n]
                    Set.object = Set.object_array[n]  ## synonym
    
                    if Set.decon_type == 'nobjects':

                        print '\t(usingimage:  %g. %s)' %(n, Set.image_list[n])

                    (object_stops[n], object_risings[n], fn) = \
                            Minimize_Object(optimization_round) #@
                                                
                    if object_stops[n] > 0:
                    
                        if object_stops[n] == Set.max_sequential_PCG_stops:
                        
                            seq_object_stops[n] += 1
                        else:
                        
                            test_object_stops[n] = object_stops[n]
                    elif object_risings[n] >= Set.max_uphill_object_PCG_steps:
                    
                        seq_object_risings[n] += 1
                        
                        if test_object_stops[n] > 0:
                        
                            seq_object_stops[n] += 1
                            check_object[n] = 'stop'
                            
                    if seq_object_stops[n] == Set.max_optimization_stops or \
                            seq_object_risings[n] == Set.max_rising_stops:
                            
                            check_object[n] = 'stop'
                    
                    if check_object[n] == 'stop':
        
                        check_Nframes_array[0][n] = 0
                    else:
        
                        check_Nframes_array[0][n] = 1
        
                    if Nobj > 1:    ## n is index over Nobj
        
                        Set.Nframes_J[n,0] = Set.J[0]   ## data fidelity
                        Set.Nframes_J[n,1] = Set.J[1]   ## object regularization
                    else:       ## n is index over npsfs
                    
                        ## Nobj = 1, but npsfs > 1
                        Set.Nframes_J[:,0] = Set.J[0]
                        Set.Nframes_J[:,1] = Set.J[1]

                    Minimization_Print(fn, Npsf)   
                    
        ###  Global Convergence Check  ###
        ## stop while loop if BOTH object and PSF estimates converge to the
        ## specified tolerance or do not improve (max stops or rising)
        for i in range(2):

            if (check_Nframes_array[i]).mean() <= \
                    (1. - Set.Nframes_fractional_convergence):

                check_Nframes[i] = 'stop'
                Nframes_array_stops[i] += 1
            else:

                check_Nframes[i] = 'go'
                Nframes_array_stops[i] = 0
            
        ## stop while loop if BOTH object and PSF estimates converge to
        ## the specified tolerance
        if check_Nframes[0] == 'stop' and check_Nframes[1] == 'stop' and \
                Nframes_array_stops[0] >= Set.max_sequential_PCG_stops and \
                Nframes_array_stops[1] >= Set.max_sequential_PCG_stops:

            Multiframe_Convergence_Print(check_Nobjects=check_Nframes[0],
                    Nobject_stops=Nframes_array_stops[0], 
                    check_Npsfs=check_Nframes[1],
                    Npsf_stops=Nframes_array_stops[1]) #@

            if Set.images_in_file == 'multiple':
            
                print
                print str(results_index)*7

            break

        if check_Nframes[1] == 'go':

            ###  PSF Minimization Loop  ###         
            print
            print '[[[ %s Minimization ]]]' %(PSF_min_type)
            print



### Correction by EHom (20130611)
### 	Commented below TEMP code out.  TEMP code was certainly for testing purposes
### 	but left behind, has been replacing all PSF estimates in Npsfs/Nobjects decon with a 
###		gaussian PSF instead of one based on the inputted PSF data and calculated in the
###		function: CalculatePSFOTFstats() !
###
## TEMP!!! 060517
#             if optimization_round == 1:
#             
#                 for i in range(Npsf):
#     
#                     Set.PSF_array[i] = F.gaussianArr(shape=Set.shape, sigma=1,
#                             integralScale=1.0, peakVal=None, orig=(0,0), wrap=1)
#                     ## FFT here to save doing repeatedly in CostGradFunction
#                     fftw.rfft(a=Set.PSF_array[i], af=Set.OTF_array[i], 
#                             inplace=False)
### TEMP

            if Set.memory_usage_level > 1: ## lower memory usage not implemented yet

                for i in range(Nobj):

                    ## FFT here to save doing repeatedly in CostGradFunction                
                    fftw.rfft(a=Set.object_array[i], af=Set.OBJECT_array[i], 
                            inplace=False)

                if Set.decon_type == 'npsfs': ## only 1 object permitted now
                
                    Set.OBJECT = Set.OBJECT_array[0]
            else: ## Set.OBJECT_array not stored

                if Set.decon_type == 'npsfs':
                
                    fftw.rfft(a=Set.object_array[0], af=Set.OBJECT,
                            inplace=False)
    
            ## object fixed, loop over N PSFs
            for n in range(Npsf):       # generally for npsfs decon
    
                if seq_PSF_stops[n] < Set.max_optimization_stops and \
                        seq_PSF_risings[n] < Set.max_rising_stops and \
                        check_PSF[n] == 'go':
    
                    Set.image = Set.image_array[n]
                    Set.inv_w = Set.inv_w_array[n].astype(Set.dtype)
                    Set.PSF = Set.PSF_array[n]  ## synonym - will change Set.PSF_array if Set.PSF changes

                    if Set.decon_type == 'npsfs':
        
                        print '\t(using image:  %g. %s)' %(n, Set.image_list[n])
                                
                    (PSF_stops[n], PSF_risings[n], fn) = \
                            Minimize_PSF(optimization_round)  

#                   ## generate OTF_array for next round of object minimization
#                   if Set.memory_usage_level > 0:
#                   
#                       fftw.rfft(af=Set.OTF_array[n]*Set.inv_Nd, a=Set.PSF, 
#                               inplace=False)
                    
                    if PSF_stops[n] > 0:
                    
                        if PSF_stops[n] == Set.max_sequential_PCG_stops:
                        
                            seq_PSF_stops[n] += 1
                        else:
                        
                            test_PSF_stops[n] = PSF_stops[n]
                    elif PSF_risings[n] == Set.max_uphill_PSF_PCG_steps:
                    
                        seq_PSF_risings[n] += 1
                        
                        if test_PSF_stops[n] > 0:
                        
                            seq_PSF_stops[n] += 1
                            check_PSF[n] = 'stop'
                            
                    if seq_PSF_stops[n] == Set.max_optimization_stops or \
                            seq_PSF_risings[n] == Set.max_rising_stops:
                        
                        check_PSF[n] = 'stop'
        
                    if check_PSF[n] == 'stop':
        
                        check_Nframes_array[1][n] = 0
                    else:
        
                        check_Nframes_array[1][n] = 1
        
                    if Npsf > 1:    ## n is index over npsfs
        
                        Set.Nframes_J[n,0] = Set.J[0]   ## data fidelity
                        Set.Nframes_J[n,2] = Set.J[2]   ## OTF constraint
                        Set.Nframes_J[n,3] = Set.J[3]   ## PSF constraint
                    else:       ## n is index over Nobj
                    
                        ## Npsf = 1, but Nobj > 1
                        Set.Nframes_J[:,0] = Set.J[0]
                        Set.Nframes_J[:,2] = Set.J[2]
                        Set.Nframes_J[:,3] = Set.J[3]
                        
                    Minimization_Print(fn, Nobj)   
    
        ### Global convergence check
        if check_Nframes[0] == 'stop' and check_Nframes[1] == 'stop' and \
                Nframes_array_stops[0] >= Set.max_sequential_PCG_stops and \
                Nframes_array_stops[1] >= Set.max_sequential_PCG_stops:

            Multiframe_Convergence_Print(check_Nobjects=check_Nframes[0],
                    Nobject_stops=Nframes_array_stops[0], 
                    check_Npsfs=check_Nframes[1],
                    Npsf_stops=Nframes_array_stops[1]) #@
    
            print
            print str(results_index)*7
            break
    
        optimization_round += 1     ## increase counter for each loop

    ### Store Minimization Results
    Set.PSF_results = Set.PSF_array         ## synonym
    Set.object_results = Set.object_array   ## synonym
    
    PCG_Decon_Print(results_index, fn, optimization_round, start_decon_time) #@
            
####


#######  FUNCTION: Imagelist_print  #######
####    [checked]
def Imagelist_Print():
    """
    Prints the list of images to be used in the deconvolution, along with
    associated AIDA parameters determined from the properties of each image
    Used ONLY in PCG_Multiframe
    """

    if Set.info_level >=1:
    
        startextratime = time.time()

        print
        print '======='

        if Set.decon_type == 'npsfs':

            print 'Image list for npsfs minimization'

        elif Set.decon_type =='nobjects':

            print 'Image list for nobjects minimization'
        
        print 'terms:', Set.terms
        print '======='

        for i in range(len(Set.image_list)):

            print '\t', i, '. ', Set.image_list[i], '\tshape =', Set.shape

            if hasattr(Set.lambda_object_Narray[i], 'shape') == 1:
    
                print '\t\tlambda_object: (%.6g, %.6g, %.6g, %.6g)' \
                        %U.mmms(Set.lambda_object_Narray[i])
            else:
                    
                print '\t\tlambda_object: ', Set.lambda_object_Narray[i]
            
            if hasattr(Set.inv_theta_center_array[i], 'shape') == 1:
    
                print '\t\ttheta: (%.6g, %.6g, %.6g, %.6g)' \
                        %U.mmms(1./Set.inv_theta_center_array[i])

            else:

                print '\t\ttheta: ', (1./Set.inv_theta_center_array[i])
                    
            if hasattr(Set.mu_Narray[i], 'shape') == 1:
    
                print '\t\tmu: (%.6g, %.6g, %.6g, %.6g)' \
                        %U.mmms(Set.mu_Narray[i])
            else:

                print '\t\tmu: ', (Set.mu_Narray[i])
            
            if hasattr(Set.lambda_OTF, 'shape') == 1:
                    
                print '\t\tlambda_OTF: (%.6g, %.6g, %.6g, %.6g)' \
                        %U.mmms(Set.lambda_OTF)
            else:
                    
                print '\t\tlambda_OTF: ', Set.lambda_OTF

            if hasattr(Set.lambda_PSF, 'shape') == 1:
                    
                print '\t\tlambda_PSF: (%.6g, %.6g, %.6g, %.6g)' \
                        %U.mmms(Set.lambda_PSF)
            else:
                    
                print '\t\tlambda_PSF: ', Set.lambda_PSF

            print '\t\tw: (%.6g, %.6g, %.6g, %.6g)' \
                    %U.mmms(1./Set.inv_w_array[i])
            print '\t\tsigma det: ', Set.sigma_det_array[i]

            if Set.initial_object_guess == 'wiener':

                print '\t\twiener =', Set.wiener

        print '-------'

        Set.extra_time += time.time() - startextratime
####


#######  FUNCTION: Final_Print  #######
####    [checked]
def Final_Print():
    """
    Prints summary AIDA deconvolution statistics
    
    Used in 'AIDA.py'
    """
    
    if Set.info_level >= 1:
    
        print
        print '======='
        print '======='
    
        if Set.decon_type == 'classical':

            print 'Total CG iterations =', Set.cum_CG_itns[0]
            print 'Total Cost Function calculations =', 
            print Set.cum_CostFunction_itns[0]
            print 'Average CG iteration time: (object only)',
            AGF.PrintTime(Set.cum_CG_time[0], Set.cum_CG_itns[0]) #@ AGF
            print
            print 'Average CostFunction calculation time: (object only)',
            AGF.PrintTime(Set.cum_CG_time[0], Set.cum_CostFunction_itns[0]) 
                    #@ AGF 
        elif Set.decon_type == 'myopic' and Set.terms[0] == 0:

            print 'Total CG iterations =', Set.cum_CG_itns[1]
            print 'Total Cost Function calculations =',
            print Set.cum_CostFunction_itns[1]
            print 'Average CG iteration time: (PSF only)',
            AGF.PrintTime(Set.cum_CG_time[1], Set.cum_CG_itns[1]) #@ AGF
            print
            print 'Average CostFunction calculation time: (PSF only)',
            AGF.PrintTime(Set.cum_CG_time[1], Set.cum_CostFunction_itns[1])
                    #@ AGF
        else:

            print 'Total CG iterations =', Set.cum_CG_itns
            print 'Total Cost Function calculations =',
            print Set.cum_CostFunction_itns
            print 'Average CG iteration time:', '\n\t[',

            for i in range(len(Set.cum_CG_itns)):
            
                AGF.PrintTime(Set.cum_CG_time[i], Set.cum_CG_itns[i]) #@ AGF

                if i != len(Set.cum_CG_itns) - 1:

                    print ' ',

            print ']  (obj PSF tot)'
            print 'Average CostFunction calculation time:', '\n\t[',

            for i in range(len(Set.cum_CostFunction_itns)):

                AGF.PrintTime(Set.cum_CG_time[i], Set.cum_CostFunction_itns[i])
                        #@ AGF

                if i != len(Set.cum_CostFunction_itns) - 1:

                    print ' ',

            print ']  (obj PSF tot)'

        print '---'
        print "Total Set-up Time =",
        AGF.PrintTime(Set.setup_time) #@ AGF
        print "\t( PSF Processing Time =", 
        AGF.PrintTime(Set.PSF_processing_time) #@ AGF
        print ")"
        print 'Total Deconvolutions =', Set.cum_deconvolutions

        print 'Total AIDA Deconvolution Time =', 
        AGF.PrintTime(Set.cum_decon_time) #@ AGF
        print
        
        if Set.cum_deconvolutions > 1:
        
            print 'Average Deconvolution Time =',
            AGF.PrintTime(Set.cum_decon_time, Set.cum_deconvolutions) #@ AGF   
            print

        print 'Total Time for AIDA Run =',
        AGF.PrintTime(time.time() - Set.start_setuptime), #@ AGF
        print   

    print '===='
    print '\nAIDA End Time:', time.asctime()
####
