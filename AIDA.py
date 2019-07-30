#!/usr/bin/env priithon
################################################################################
#
#   File:       AIDA.py
#
#   Summary:    Master Python file that runs the Adaptive Image Deconvolution 
#               Algorithm (AIDA)
#
#   Authors:    Erik F.Y. Hom (Sedat Lab, UCSF) with assistance from
#               Sebastian Haase (Sedat Lab, UCSF)
#               Clement Chalumeau (SETI Institute)
#               Franck Marchis (SETI Institute)
#
#   Other:      See 'AIDA_version.py' for date and version details
#               See 'LICENSE.txt' for use, license, and reference details
#
################################################################################


#from Priithon.all import *
#Used from Priithon: N (numpy)
import numpy as N
import time, os, sys, getopt
import AIDA_Functions as Fun
import AIDA_Settings as Set
import AIDA_Check_Settings as Chk
import AIDA_Version as Ver

#20090716 N.seterr(divide="raise", over=None, under=None, invalid=None)

__author__ = "Erik F.Y. Hom, Sebastian Haase, Timothy K. Lee (Sedat Lab, UCSF)"
__license__ = "See 'LICENSE.txt' for details"
__version__ = "AIDA version: " + str(Ver.version)
__date__ = "Modification date: " + str(Ver.date)


def usage():
    print """
usage: priithon AIDA.py   -S settings.py
                          -i image *
                          -h psf *
                          -k key label
                          -r results output directory
                          -d decon type
                          -n initial psf guess
                          -3 dimension of data
                          -f output file format
                          -b background
                          -a darkImageFile
                          -s sigma detector value
                          -L lambda_object
                          -T theta
                          -O lambda_OTF
                          -H lamda_PSF
                          -l lambda_object_scaling
                          -t theta_scaling
                          -v lambda_OTF_scaling
                          -u lambda_PSF_scaling

       (read AIDA.py for additional options)
"""
							### EHom (20130611): added "n" option above to specify initial guess of PSFs
    sys.exit(1)
    

#######  FUNCTION: main()  #######
####
def main():
    ParseCommandLine()
    RunAIDA()
####

#######  FUNCTION: ParseCommandLine()  #######
####
def ParseCommandLine():
    try:

        (opts, args) = getopt.getopt(sys.argv[1:],
                "S:i:h:k:y:r:d:n:p:3f:b:a:s:L:T:O:H:l:t:v:u",
                ["info=", "terms=", "mem=", "guess=", "clean=", "backper=", 
                 "nsig=", "psfthresh=", "otfthresh=", "fill=", "vfloor=",
                 "vrad=", "originpsf=", "deriv=", "lapl=", "fullft=", "objtol=",
                 "psftol=", "pcgiterarr=", "objiterarr=", "psfiterarr=", 
                 "maxiter=", "maxseqstops=", "maxupobj=", "maxuppsf=", 
                 "fracconverge=", "objccgtol=", "psfccgtol=", "objccgiter=", 
                 "psfccgiter=", "ccgxmin=", "ccgivec=", "lgrid=", "tgrid=", 
                 "debug"]) # seb 20090714
                ### functionality of a number of these short and long options
                ### still (as of 2007/03) need to be tested
    except getopt.GetoptError:

        usage()
        sys.exit(2)


    for (opt, arg) in opts:

        if opt == '-S':

            if not os.path.isfile(arg):

                print "'Settings_File' ('-S' option) specified does not exist!"
                sys.exit(3)

            ReadSettingsFile(settings_filename=arg)


    ########################################################################
    #######  COMMAND LINE INTERFACE / PROCESS INPUT ; DON'T TOUCH!!  #######
    ########################################################################

    ### ASSIGN VALUES PROVIDED IN COMMAND LINE ###
    for (opt, arg) in opts:

        if opt == '-i':             Set.image_filenames = arg   ##<- these 2 are 
        if opt == '-h':             Set.PSF_filenames = arg     ##<- required
        if opt == '-k':             Set.dataset_label = arg
        if opt == '-y':             Set.dataset_list = arg
        if opt == '-r':             Set.results_directory = str(arg)
        if opt == '-d':             Set.decon_type = str(arg)
        if opt == '-n':             Set.initial_psf = str(arg)
        if opt == '-p':             Set.precision = str(arg)
        if opt == '-3':             Set.dimension = 3
        if opt == '-f':             Set.output_format = str(arg)
        if opt == '-b':             Set.background = arg
        if opt == '-s':             Set.sigma_det = arg
        if opt == '-a':             Set.dark_image = arg
        if opt == '-L':             Set.lambda_object_input = arg
        if opt == '-T':             Set.theta_input = arg
        if opt == '-O':             Set.lambda_OTF_input = arg
        if opt == '-H':             Set.lambda_PSF_input = arg
        if opt == '-l':             Set.lambda_object_scaling = arg
        if opt == '-t':             Set.theta_scaling = arg
        if opt == '-v':             Set.lambda_OTF_scaling = arg
        if opt == '-u':             Set.lambda_PSF_scaling = arg

        if opt == '--info':         Set.info_level = arg
        if opt == '--terms':        Set.terms = arg
        if opt == '--mem':          Set.memory_usage_level = arg
        if opt == '--guess':        Set.initial_object_guess = str(arg)
        if opt == '--clean':        Set.cleaning = arg
        if opt == '--backper':      Set.PSF_background_percent = arg
        if opt == '--nsig':         Set.nsigmas = arg
        if opt == '--psfthresh':    Set.PSF_threshold_percent = arg
        if opt == '--otfthresh':    Set.OTF_threshold_percent = arg
        if opt == '--fill':         Set.fill = arg
        if opt == '--ufloor':       Set.u_floor_per_dim = arg
        if opt == '--vfloor':       Set.v_floor_per_dim = arg
        if opt == '--vrad':         Set.vrad = arg
        if opt == '--originpsf':    Set.origin_centered_PSF_flag = arg
        if opt == '--deriv':        Set.derivative_operator = str(arg)
        if opt == '--lapl':         Set.laplacian_operator = arg
        if opt == '--objtol':       Set.object_PCG_tolerance = arg
        if opt == '--psftol':       Set.PSF_PCG_tolerance = arg
        if opt == '--pcgiterarr':   Set.PCG_iter_array = arg
        if opt == '--objiterarr':   Set.object_PCG_iter_array = arg
        if opt == '--psfiterarr':   Set.PSF_PCG_iter_array = arg
        if opt == '--maxiter':      Set.max_total_PCG_blocks = arg
        if opt == '--maxseqstops':  Set.max_sequential_PCG_stops = arg
        if opt == '--maxupobj':     Set.max_uphill_object_PCG_steps = arg
        if opt == '--maxuppsf':     Set.max_uphill_PSF_PCG_steps = arg
        if opt == '--fracconverge': Set.Nframes_fractional_convergence = arg
        if opt == '--objccgtol':    Set.object_CCG_tolerance = arg
        if opt == '--psfccgtol':    Set.PSF_CCG_tolerance = arg
        if opt == '--objccgiter':   Set.max_object_CCG_iter = arg
        if opt == '--psfccgiter':   Set.max_PSF_CCG_iter = arg
        if opt == '--ccgxmin':      Set.xmin_value = arg
        if opt == '--ccgivec':      Set.ivec_value = arg
        if opt == '--lgrid':        Set.Lgrid = arg
        if opt == '--tgrid':        Set.Tgrid = arg
        if opt == '--debug':        Set.debug = True


def ReadSettingsFile(settings_filename):
    """
    useful function if AIDA is used from PyShell
    """
    Set.Settings_File = settings_filename
    execfile(Set.Settings_File, Set.__dict__)
    


#######  FUNCTION: main()  #######
####
def RunAIDA(settings_filename=None):
    """
    this is the *real* main function that starts AIDA decon 
    - exclusing commandline handling
    if `settings_filename` is given here, it will get `execfile`d in `Set`
       (note that `ParseCommandLine()` it doing this by itself)
    """
    if settings_filename is not None:
        ReadSettingsFile(settings_filename)
    Set.start_run_time = time.time()
    Chk.Check_Settings(Set)

    # ##########################################################################
    # ##########  ADAPTIVE IMAGE DECONVOLUTION ALGORITHM ; DON'T TOUCH!! #######
    # ##########################################################################
    ## loop over DIFFERENT datasets
    for i in range(len(Set.dataset_image_prefix_list)):

        Set.results_directory = Set.results_directory_list[i]
        Set.image_filenames = Set.dataset_image_prefix_list[i]
        Set.PSF_filenames = Set.dataset_PSF_prefix_list[i]
        Set.dataset_label = Set.dataset_label_list[i]
        Set.background = Set.dataset_background_list[i]
        Set.sigma_det = Set.dataset_sigma_det_list[i]
        Set.dark_image = Set.dataset_dark_image_list[i]

        try:
            _oldErrSettings = N.seterr(**Set.debug_numpy_error_settings)

            Fun.AIDA_Setup()
            Fun.Start_Print()

            if Set.decon_type == 'classical':
                Fun.DeconSingleFrame()

            elif Set.decon_type == 'myopic':
                Fun.DeconSingleFrame()

            elif Set.decon_type == 'nobjects':
                Fun.DeconMultiFrame()

            elif Set.decon_type == 'npsfs':
                Fun.DeconMultiFrame()

            Fun.Final_Print()
        finally:
            N.seterr(**_oldErrSettings)
####

if __name__ == "__main__":
    
    if False:       # profile, for debugging purposes
    
        import profile
        profile.run("main()", "profile2")
    else:
        main()
    
