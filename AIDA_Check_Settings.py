################################################################################
#
#   File:       AIDA_Check_Settings.py
#
#   Summary:    Function to check input values for AIDA settings loaded from
#               'AIDA_Settings.py', user specified settings file, or the
#               command line
#
#   Author:     Erik F.Y. Hom (Sedat Lab, UCSF)
#
#   Other:      See 'AIDA_version.py' for date and version details
#               See 'LICENSE.txt' for use, license, and reference details
#
#   Modifications : Clement Chalumeau (SETI intern) 2016/04/06
#                Remove all priithon occurences
################################################################################


# from Priithon.all import *
import numpy as N

import os, glob, types


#######  Function: Check_Settings  #######
####    [checked]
def Check_Settings(Set):   #Fun
    
    ### if opt == '-y'
    if Set.dataset_list is not None:
    
        Set.image_filenames = []
        Set.PSF_filenames = []
        Set.dataset_label = []
    
        if type(Set.dataset_list) in (types.ListType, types.TupleType):
        
            for i in range(len(Set.dataset_list)):
            
                file_tuple = Set.dataset_list[i]
                length = len(file_tuple)
                
                if length > 3:
                
                    message = "\nThere are too many tuple elements!:\n\t" + \
                            str(file_tuple)
                    raise ValueError, message
    
                elif length < 2:
    
                    message = "\n'Each list element of 'dataset_list' " + \
                            "('-y' option) must contain at least an " + \
                            "'image_filenames' and a 'PSF_filenames'\n" + \
                            "You provided: " + str(file_tuple) + "\n\t" + \
                            "at position " + str(i)
                    raise ValueError, message
                
                Set.image_filenames.append(file_tuple[0])
                Set.PSF_filename.append(file_tuple[1])
    
                if len(glob.glob(file_tuple[0] + '*')) == 0:
    
                    message = "\n'image_filenames' value in 'dataset_list' " + \
                            "('-y' option):\n\t'" + str(file_tuple[0]) + \
                            "' does not point to a file that exists!"
                    raise ValueError, message
    
                if len(glob.glob(file_tuple[1] + '*')) == 0:
    
                    message = "\n'PSF_filenames' value in 'dataset_list' " + \
                            "('-y' option):\n\t'" + str(file_tuple[1]) + \
                            "' does not point to a file that exists!"
                    raise ValueError, message
    
                if length == 3:
                
                    Set.dataset_label.append(file_tuple[2])
                elif length >= 2:
                
                    Set.dataset_label.append(None)
        else:
    
            message = "\n'dataset_list' ('-y' option) supplied must be a " + \
                    "list of tuple string names of the form:\n\t" + \
                    "[('image_filenames', 'PSF_filenames', {'dataset_label'})"+\
                    ", ... ]\n\twhere {'dataset_label'} is optional"
            raise ValueError, message
    
    
    ### if opt == '-i'
    if Set.dataset_list is None:
    
        if Set.image_filenames == '':
    
            Set.image_filenames = raw_input("Image file to deconvolve: ")
    
            if Set.image_filenames == '':
    
                Set.image_filenames = FN() # where is GUI !?  HOW DO??
        elif Set.image_filenames is None:
        
            message = "\n'image_filenames' value ('-i' option) is 'None'!\n" + \
                    "Please check your Settings File entry and/or make " + \
                    "sure you have used the\n'-S' option to specify a " + \
                    "non-default Settings File"
            raise ValueError, message
    
        if type(Set.image_filenames) is types.StringType:

            Set.dataset_image_prefix_list = [Set.image_filenames]
                
            if len(glob.glob(Set.dataset_image_prefix_list[0] + '*')) == 0:
    
                message = "\n'image_filenames' value ('-i' option):\n\t'" + \
                        str(Set.dataset_image_prefix_list) + "' does not " + \
                        "point to a file that exists!"
                raise ValueError, message


        elif type(Set.image_filenames) in (types.ListType, types.TupleType):

            # assume a simple list
            if type(Set.image_filenames[0]) is types.StringType:
        
                # check existence of files
                for file in Set.image_filenames:
            
                    if len(glob.glob(file + '*')) == 0:
    
                        message = "\n'image_filenames' ENTRY ('-i' option):" + \
                                "\n\t'" + str(file) + "' does not point to " + \
                                "a file that exists!"
                        raise ValueError, message
    
                if Set.decon_type.lower() in ('nobjects', 'npsfs'):
                
                    Set.dataset_image_prefix_list = [Set.image_filenames[:]]
                else:
                
                    Set.dataset_image_prefix_list = Set.image_filenames[:]

                
            # assume we have a 'list of a list'
            elif type(Set.image_filenames[0]) in \
                    (types.ListType, types.TupleType):
    
                for i in range(len(Set.image_filenames)):
    
                    for file in Set.image_filenames[i]:
                
                        if len(glob.glob(file + '*')) == 0:
    
                            message = "\n'image_filenames' ENTRY ('-i' " + \
                                    "option): " + str(file) + "\nin LIST " + \
                                    "does not point to a file that exists!"
                            raise ValueError, message
    
                if Set.decon_type.lower() in ('nobjects', 'npsfs'):

                    Set.dataset_image_prefix_list = Set.image_filenames[:]
                else:       # collapse to a single list for monoframe decon

                    Set.dataset_image_prefix_list = Set.image_filenames[0]
        else:
    
            message = "\n'image_file' value ('-i' option) should be a " + \
                    "filename string or a list of filename strings!"
            raise ValueError, message

    else:
    
        message = "\n'image_file' value ('-i' option) should be a filename " + \
                "string or a list of filename strings!"
        raise ValueError, message
    
    
    ### if opt == '-h'
    if Set.dataset_list is None:
    
        if Set.PSF_filenames == '':
    
            Set.PSF_filenames = raw_input("PSF(s) to use: ")
    
            if Set.PSF_filenames == '':
    
                Set.PSF_filenames = FN() # where is GUI !?  HOW DO??
        elif Set.PSF_filenames is None:
        
            message = "\n'PSF_filename' value ('-h' option) is 'None'!\n" + \
                    "Please check your Settings File entry and/or make " + \
                    "sure you have used the\n'-S' option to specify a " + \
                    "non-default Settings File"
            raise ValueError, message
    
        if type(Set.PSF_filenames) is types.StringType:
        
            if len(glob.glob(Set.PSF_filenames + '*')) == 0:
    
                message = "\n'PSF_filenames' value ('-h' option):\n\t'" + \
                        str(Set.PSF_filenames) + "' does not point to a " + \
                        "file that exists!"
                raise ValueError, message
        
            # if image_filename was a List
            if type(Set.dataset_image_prefix_list) in \
                    (types.ListType, types.TupleType):

                # if was a list of a list
                if type(Set.dataset_image_prefix_list[0]) in \
                        (types.ListType, types.TupleType):
                
                        Set.dataset_PSF_prefix_list = \
                                [len(Set.dataset_image_prefix_list) * \
                                [Set.PSF_filenames]]
                else:   # just a monoframe decon list
                
                    Set.dataset_PSF_prefix_list = \
                            len(Set.dataset_image_prefix_list) * \
                            [Set.PSF_filenames]
            else:
        
                Set.dataset_PSF_prefix_list = [Set.PSF_filenames]
        elif type(Set.PSF_filenames) in (types.ListType, types.TupleType):
    
            #if len(Set.PSF_filenames) != len(Set.dataset_image_prefix_list) and \
			# Edited by TBM...check against image list, not image prefix list
            if len(Set.PSF_filenames) != len(Set.image_filenames) and \
                    len(Set.PSF_filenames) != 1:
                
                message = "\n'PSF_filename' entry ('-h' option) length of:" + \
                        "\n" + str(len(Set.PSF_filenames)) + "\ndoes not " + \
                        "match the 'image_filenames' entry ('-i' option) " + \
                        "length of: \n" + \
                        str(len(Set.dataset_image_prefix_list)) + "!"
                #raise ValueError, message
                #Edited by TBM...just don't throw error...
                print message
    
            elif type(Set.PSF_filenames[0]) is types.StringType:
        
                # check that files exist
                for file in Set.PSF_filenames:
    
                    if len(glob.glob(file + '*')) == 0:
    
                        message = "\n'PSF_filenames' value ('-h' option):\n" + \
                                "\t'" + str(file) + "' does not point to a " + \
                                "file that exists!"
                        raise ValueError, message
                    
                if len(Set.PSF_filenames) == 1:
                
                    Set.dataset_PSF_prefix_list = \
                            len(Set.dataset_image_prefix_list) * \
                            [Set.PSF_filenames[0]]
                            
                else:
                
                    Set.dataset_PSF_prefix_list = Set.PSF_filenames[:]
                
            elif type(Set.PSF_filenames[0]) in (types.ListType, types.TupleType):
            
                # check that files exist
                for i in range(len(Set.PSF_filenames)):
        
                    for file in Set.PSF_filenames[i]:
    
                        if len(glob.glob(file + '*')) == 0:
    
                            message = "\n'PSF_filenames' entry value ('-h' " + \
                                    "option):\n" + "\t'" + str(file) + \
                                    "' does not point to a file that exists!"
                            raise ValueError, message
                            
                Set.dataset_PSF_prefix_list = Set.PSF_filenames[:]
        else:
    
            message = "\n'PSF_filenames' ENTRY ('-h' option) should be a " + \
                    "filename string or a list of filename strings " + \
                    "(or a list of a list)!"
            raise ValueError, message

    else:
    
        message = "\n'PSF_filenames' value ('-h' option) should be a " + \
                "filename string or a list of filename strings (or a list " + \
                "of a list)!"
        raise ValueError, message
    
    
    ### if opt == '-k'
    if Set.dataset_list is None:
        
        if type(Set.dataset_label) in (types.StringType, types.NoneType):
    
            if Set.dataset_label == '':
    
                Set.dataset_label = None
            
            Set.dataset_label_list = len(Set.dataset_image_prefix_list) * \
                    [Set.dataset_label]
        elif type(Set.dataset_label) in (types.ListType, types.TupleType):
    
            if len(Set.dataset_label) != len(Set.dataset_image_prefix_list) and \
                    len(Set.dataset_label) != 1:
    
                message = "\n'dataset_label' array ('-k' option) supplied " + \
                        "of length: " + str(len(Set.dataset_label_list)) + \
                        "\n" + "does not match the 'image_file' array " + \
                        "('-i' option) supplied of length: " + \
                        str(len(Set.dataset_image_prefix_list)) + "!"
                raise ValueError, message
                
            elif len(Set.dataset_label) == len(Set.dataset_image_prefix_list):
            
                for i in range(len(Set.dataset_label)):
                    
                    if Set.dataset_label[i] == '':
    
                        Set.dataset_label[i] = None
                    elif type(Set.dataset_label[i]) in \
                            (types.ListType, types.TupleType):
                        
                        if len(Set.dataset_label[i]) == 1:
                        
                            Set.dataset_label[i] = Set.dataset_label[i][0]
                        else:
                        
                            message = "\n'dataset_label' ('-h' option) " + \
                                    "cannot be a list of a list of more " + \
                                    "than one element!\nPlease check " + \
                                    "entry #" + str(i)
                            raise ValueError, message
                                    
                Set.dataset_label_list = Set.dataset_label[:]

            # assume len(Set.dataset_label) == 1
            elif type(Set.dataset_label[0]) in (types.StringType,
                    types.NoneType):
    
                if Set.dataset_label == '':
    
                    Set.dataset_label = None

                Set.dataset_label_list = len(Set.dataset_image_prefix_list) * \
                        [Set.dataset_label[0]]
            else:
                
                message = "\n'dataset_label' entry ('-h' option) length of:" + \
                        "\n" + str(len(Set.dataset_label)) + "\ndoes not " + \
                        "match the 'image_filenames' entry ('-i' option) " + \
                        "length of: \n" + \
                        str(len(Set.dataset_image_prefix_list)) + "!"
                raise ValueError, message
        else:
                
            message = "\n'dataset_label' value ('-h' option) should be a " + \
                    "string to label the decon run \nor a list of label " + \
                    "strings for each dataset!"
            raise ValueError, message
    else:
                
        message = "\n'dataset_label' value ('-h' option) should be 'None', " + \
                "a string to label the decon run, \nor a list of label " + \
                "strings for each dataset!"
        raise ValueError, message
    
    
    ### if opt == '-r'
    if type(Set.results_directory) is types.StringType:
    
        Set.results_directory_list = len(Set.dataset_image_prefix_list) * \
                [Set.results_directory]
    elif type(Set.results_directory) in (types.ListType, types.TupleType):
    
        if len(Set.results_directory) == 1:
        
            Set.results_directory_list = len(Set.dataset_image_prefix_list) * \
                    [Set.results_directory[0]]
    
        elif len(Set.results_directory) == len(Set.dataset_image_prefix_list):
    
            Set.results_directory_list = Set.results_directory[:]
        else:
        
            message = "\n'results_directory' array ('-r' option) supplied of " + \
                    "length: " + str(len(Set.results_directory)) + "\n" + \
                    "does not match the 'image_filenames' array ('-i' option) " + \
                    "supplied of length: " + \
                    str(len(Set.dataset_image_prefix_list)) + "!"
            raise ValueError, message
    else:
    
        message = "\n'results_directory' value ('-r' option) should be a name " + \
                "of a directory (string) or a list of directory names!"
        raise ValueError, message
    
    
    ### if opt = '-d'
    if Set.decon_type.lower() not in ('classical', 'myopic', 'nobjects', 
            'npsfs', 'si', 'siclassical'):
    
        message = "\n'decon_type' specification ('-d' option) is not valid!\n" + \
                "Options are: 'classical', 'myopic', 'nobjects', or 'npsfs'"
        raise ValueError, message
    Set.decon_type = Set.decon_type.lower()
    
    
    ### if opt == '-p'
    if Set.precision.lower() == 'single':
    
        Set.dtype = N.float32
        Set.complex_type = N.complex64
    elif Set.precision.lower() == 'double':
    
        Set.dtype = N.float64
        Set.complex_type = N.complex128
    else:
    
        message = "\n'precision_type' ('-p' option) must be 'double' or 'single'!"
        raise ValueError, message
    Set.precision = Set.precision.lower()

    
    ### if opt == '-3'
    if Set.dimension not in (2,3):
    
        message = "\nImage 'dimension' supplied is not valid!\n" + \
                "AIDA can only handle images that are 2D or 3D\n"
        raise ValueError, message
    
    elif Set.dimension == 2:
    
        Set.axesFT = (-2,-1)
    else:
    
        Set.axesFT = (-3,-2,-1)
    
    
    ### if opt == '-f'
    if Set.output_format not in ('f', 'm','t', None):
    
        print "WARNING: 'output_format' value ('-f' option) must be 'f' or 'm' or 't'"
        print "for FITS or MRC or TIFF formats, respectively"
        print "Will default to using 'output_format' = 'image_format' instead"
        
        Set.output_format = None
    
    
    ### if opt == '-b'
    if type(Set.background) in (types.IntType, types.LongType, types.FloatType, 
            types.NoneType): # or isinstance(Set.background,  N.number):
    
        Set.dataset_background_list = len(Set.dataset_image_prefix_list) * \
                [Set.background]
    elif type(Set.background) in (types.ListType, types.TupleType):
    
        if len(Set.background) == 1:
        
            Set.dataset_background_list = len(Set.dataset_image_prefix_list) * \
                    [Set.background[0]]
        elif len(Set.background) == len(Set.dataset_image_prefix_list) or \
                len(Set.dataset_image_prefix_list) == 1:
    
            for i in range(len(Set.background)):
            
                if type(Set.background[i]) in (types.ListType, types.TupleType):
                
                    if len(Set.background[i]) != 1 and \
                            len(Set.background[i]) != \
                            len(Set.dataset_image_prefix_list[i]):
                            
                        message = "\n'background' array ('-b' option) " + \
                                "supplied is a list of a list that\n" + \
                                "does not match that of 'image_filenames'" + \
                                "array ('-i' option)!"
                        raise ValueError, message   
                    elif len(Set.background[i]) == 1:
                    
                        Set.background[i] = Set.background[i][0]

            # note that Set.background could be a list of a list for:
            # (1) Nframes decon; (2) when multiple images are in a single file
            Set.dataset_background_list = Set.background[:]
        else:
        
            message = "\n'background' array ('-b' option) supplied of " + \
                    "length: " + str(len(Set.background)) + "\n" + \
                    "does not match the 'image_filenames' array ('-i' option) " + \
                    "supplied of length: " + \
                    str(len(Set.dataset_image_prefix_list)) + "!"
            raise ValueError, message
    else:
    
        message = "\n'background' image value ('-b' option) entered must " + \
                "be a number or a list of numbers!"
        raise ValueError, message
    
    
    ### if opt == '-s'  
    if type(Set.sigma_det) in (types.IntType, types.LongType, types.FloatType, 
            types.NoneType):
        
            Set.dataset_sigma_det_list = len(Set.dataset_image_prefix_list) * \
                    [Set.sigma_det]
    elif type(Set.sigma_det) in (types.ListType, types.TupleType):
        
#       if type(Set.dataset_image_prefix_list[0]) in \
#               (types.ListType, types.TupleType):
#
#       if len(Set.sigma_det) == 1:
#
#           if Set.sigma_det[0] == '':
#           
#               Set.sigma_det[0] = None
#
#           Set.dataset_sigma_det_list = len(Set.dataset_image_prefix_list) * \
#                   [Set.sigma_det[0]]
        if len(Set.sigma_det) == len(Set.dataset_image_prefix_list) or \
                len(Set.dataset_image_prefix_list) == 1:
            
            for i in range(len(Set.sigma_det)):
            
                if Set.sigma_det[i] == '':
                
                    Set.sigma_det[i] = None
                elif type(Set.sigma_det[i]) in (types.ListType, types.TupleType):
                
                    if len(Set.sigma_det[i]) != 1 and \
                            len(Set.sigma_det[i]) != \
                            len(Set.dataset_image_prefix_list[i]):
                            
                        message = "\n'sigma_det' array ('-s' option) " + \
                                "supplied is a list of a list that\n" + \
                                "does not match that of 'image_filenames'" + \
                                "array ('-i' option)!"
                        raise ValueError, message   
                    elif len(Set.sigma_det[i]) == 1:
                    
                        if Set.sigma_det[i][0] == '':
                        
                            Set.sigma_det[i][0] = None
                        
                        Set.sigma_det[i] = Set.sigma_det[i][0]

            # note that Set.sigma_det could be a list of a list for:
            # (1) Nframes decon; (2) when multiple images are in a single file
            Set.dataset_sigma_det_list = Set.sigma_det[:]
        else:
    
            message = "\n'sigma_det' array ('-s' option) supplied of " + \
                    "length: " + str(len(Set.sigma_det)) + "\n" + \
                    "does not match the 'image_file' array ('-i' option) " + \
                    "supplied of length: " + \
                    str(len(Set.dataset_image_prefix_list)) + "!"
            raise ValueError, message
    else:
        
        message = "\n'sigma_det' value entered must be a number or a list " + \
                "of numbers!"
        raise ValueError, message
            
    
    ### if opt == '-a'
    if Set.dark_image == '' or Set.dark_image is None:
    
        Set.dataset_dark_image_list = len(Set.dataset_image_prefix_list) * [None]
    elif type(Set.dark_image) is types.StringType:
        
        if os.path.isfile(Set.dark_image):
    
            Set.dataset_dark_image_list = len(Set.dataset_image_prefix_list) * \
                    [os.path.abspath(Set.dark_image)]
        else:
    
            message = "\n'dark_image' file ('-a' option) specified does " + \
                    "not exist!"
            raise ValueError, message
    elif type(Set.dark_image) in (types.ListType, types.TupleType):
        
        if len(Set.dark_image) == 1 and Set.dark_image[0] in ('', None):
    
            Set.dataset_dark_image_list = len(Set.dataset_image_prefix_list) * \
                    [None]
        elif len(Set.dark_image) == len(Set.dataset_image_prefix_list) or \
                len(Set.dataset_image_prefix_list) == 1:
        
            for i in range(len(Set.dark_image)):
            
                if Set.dark_image[i] == '':
            
                    Set.dark_image[i] = None
                elif type(Set.dark_image[i]) is types.StringType:
                
                    if not os.path.isfile(Set.dark_image[i]):
        
                        message = "\n'dark_image' file ENTRY " + \
                                "('-a' option):\n" + str(Set.dark_image[i]) + \
                                " does not exist!"
                
                elif type(Set.dark_image[i]) in \
                        (types.ListType, types.TupleType):
                
                    if len(Set.dark_image[i]) != 1 and \
                            len(Set.dark_image[i]) != \
                            len(Set.dataset_image_prefix_list[i]):
                            
                        message = "\n'dark_image' array ('-a' option) " + \
                                "supplied is a list of a list that\n" + \
                                "does not match that of 'image_filenames'" + \
                                "array ('-i' option)!"
                        raise ValueError, message   
                    elif len(Set.dark_image[i]) == 1:
                    
                        if Set.dark_image[i][0] == '':
                        
                            Set.dark_image[i][0] = None
                        elif not os.path.isfile(Set.dark_image[i]):
        
                            message = "\n'dark_image' file ENTRY " + \
                                    "('-a' option):\n" + \
                                    str(Set.dark_image[i]) + " does not exist!"
                                    
                        Set.dark_image[i] = Set.dark_image[i][0]
                    else:
                    
                        for dark in Set.dark_image[i]:
                        
                            if not os.path.isfile(dark):
        
                                message = "\n'dark_image' file ENTRY " + \
                                        "('-a' option):\n" + str(dark) + \
                                        " does not exist!"

            # note that Set.sigma_det could be a list of a list for:
            # (1) Nframes decon; (2) when multiple images are in a single file
            Set.dataset_dark_image_list = Set.dark_image[:]

        else:
    
            message = "\n'dark_image' array ('-s' option) supplied of " + \
                    "length: " + str(len(Set.dark_image)) + "\n" + \
                    "does not match the 'image_file' array ('-i' option) " + \
                    "supplied of length: " + \
                    str(len(Set.dataset_image_prefix_list)) + "!"
            raise ValueError, message
    else:
    
        message = "\n'dark_file' value ('-a' option) should be a filename " + \
                "string or a list of filename strings!"
        raise ValueError, message
    
    
    ### if opt == '-L'
    if Set.lambda_object_input is not None and (type(Set.lambda_object_input) \
            not in (types.IntType, types.LongType, types.FloatType)):
    
        message = "\n'lambda_object' input value ('-L' option) must be" + \
                "a number!"
        raise ValueError, message
    
    
    ### if opt == '-T'
    if Set.theta_input is not None and (type(Set.theta_input) not in \
            (types.IntType, types.LongType, types.FloatType)):
    
        message = "\n'theta' input value ('-T' option) must be a number!"
        raise ValueError, message
        
        
    ### if opt == '-O'
    if Set.lambda_OTF_input is not None and (type(Set.lambda_OTF_input) \
            not in (types.IntType, types.LongType, types.FloatType)):
    
        message = "\n'lambda_OTF' input value ('-O' option) must be" + \
                "a number!"
        raise ValueError, message
    
    
    ### if opt == '-H'
    if Set.lambda_PSF_input is not None and (type(Set.lambda_PSF_input) \
            not in (types.IntType, types.LongType, types.FloatType)):
    
        message = "\n'lambda_PSF' input value ('-H' option) must be" + \
                "a number!"
        raise ValueError, message
    
    ### if opt == '-l'
    if type(Set.lambda_object_scaling) not in (types.IntType, types.LongType, 
            types.FloatType):
    
        message = "\n'lambda_object_scaling' input value ('-l' option) " + \
                "must be a number!"
        raise ValueError, message
    else:
    
        Set.lambda_object_scaling = float(Set.lambda_object_scaling)
    

    ### if opt == '-t'
    if type(Set.theta_scaling) not in (types.IntType, types.LongType, 
            types.FloatType):
    
        message = "\n'theta_scaling' input value ('-t' option) must be " + \
                "a number!"
        raise ValueError, message
    else:
    
        Set.theta_scaling = float(Set.theta_scaling)
    

    ### if opt == '-v'
    if type(Set.lambda_OTF_scaling) not in (types.IntType, types.LongType, 
            types.FloatType):
    
        message = "\n'lambda_OTF_scaling' input value ('-v' option) must " + \
                "be a number!"
        raise ValueError, message
    else:
    
        Set.lambda_OTF_scaling = float(Set.lambda_OTF_scaling)
    
    
    ### if opt == '-u'
    if type(Set.lambda_PSF_scaling) not in (types.IntType, types.LongType, 
            types.FloatType):
    
        message = "\n'lambda_PSF_scaling' input value ('-u' option) must " + \
                "be a number!"
        raise ValueError, message
    else:
    
        Set.lambda_PSF_scaling = float(Set.lambda_PSF_scaling)



    ######################################
    ###### === long name options === #####
    ######################################
    
    ### if opt == '--info'
    if Set.info_level not in (0,1,2,3):
    
        message = "\n'infolevel' value ('--info' option) must be [0-3]!"
        raise ValueError, message
    
    
    ### if opt == '--terms'
    if type(Set.terms) is not types.TupleType:
    
        message = "\n'terms' value ('--terms' option) entered must be a tuple " + \
                "of the form: \t(1, 1, 1, 0)"
        raise ValueError, message
        
    elif len(Set.terms) < 1 or len(Set.terms) > 4:
    
        message = "\n'terms' needs to be a tuple of length > 1 and < 4!"
        raise ValueError, message
    
    for each in Set.terms:
        
        if each not in (0,1):
            
            message = "\n'terms' entries must be '1' or '0'!"
            raise ValueError, message
                
    if len(Set.terms) >= 2 and Set.terms[1] == 1 and N.sum(Set.terms) == 1:
            
        message = "\n'terms'='" + str(Set.terms) + "' are invalid!"
        raise ValueError, message
        # cannot have terms = (0,1,0,0)
    elif N.sum(Set.terms) == 0:
    
        message = "\n'terms'='" + str(Set.terms) + "' are invalid!"
        raise ValueError, message
        # cannot have terms = (0,0,0,0)
    
    
    ### if opt == '--mem'       ##!! not fully implemented yet !!##
    if Set.memory_usage_level not in (0,1,2,3):
    
        message = "\n'memusagelevel' value ('--mem' option) must be [0-3]!"
        raise ValueError, message
    
    
    ### if opt == '--guess'
    if Set.initial_object_guess not in ('zero', 'wiener', 'image'):
    
        message = "\n'initial_guess' value ('--guess' option) must be:\n\t" + \
                "'zero', 'wiener', or 'image'!"
        raise ValueError, message
    

    ### if opt == '--clean'
    if type(Set.cleaning) not in (types.ListType, types.TupleType) or \
            len(Set.cleaning) != 4:
        
        message = "\n'clean' ('--clean' option) must be a list or tuple " + \
                "of 0's and 1's of length 4 (e.g., (1,1,1,1))!"
        raise ValueError, message
    else:
        
        for i in range(4):
            
            if Set.cleaning[i] not in (0,1):
                
                message = "\n'clean' ('--clean' option) must be a list " + \
                        "or tuple of 0's and 1's\nentry '" + str(i) + \
                        "' has a value of " + str(clean[i]) + "!"
                raise ValueError, message
    

    ### if opt == '--center'
    if Set.PSF_centering not in ('unknown', 'origin', 'array_center'):
        message = "\n'PSF_centering' ('--center' option) value must be " + \
                "'unknown', 'origin', or 'array_center'!"
        raise ValueError, message
    elif Set.PSF_centering == 'unknown':
    
        Set.center = None
    else:
    
        Set.center = Set.PSF_centering
    
    
    ### if opt == '--backper'
    if type(Set.PSF_subtract_background_percent) not in (types.IntType, 
            types.FloatType):
    
        message = "\n'PSF_subtract_background_percent' value ('--backper' " + \
                "option) must be a number!"
        raise ValueError, message
    
    
    ### if opt == '--nsig'
    if type(Set.nsigmas) not in (types.IntType, types.FloatType):
    
        message = "\n'nsigmas' value ('--nsig' option) must be a number!"
        raise ValueError, message
    
    
    ### if opt == '--psfthresh'
    if type(Set.PSF_threshold_percent) not in (types.IntType, types.FloatType):
    
        message = "\n'PSF_threshold_percent' value ('--psfthresh' " + \
                "option) must be a number!"
        raise ValueError, message
    
    
    ### if opt == '--otfthresh'
    if type(Set.OTF_threshold_percent) not in (types.IntType, types.FloatType):
    
        message = "\n'OTF_threshold_percent' value ('--otfthresh' " + \
                "option) must be a number!"
        raise ValueError, message
    
    
    ### if opt == '--fill'
    if type(Set.fill) not in (types.IntType, types.FloatType):
    
        message = "\n'fill' value ('--fill' option) must be a number!"
        raise ValueError, message


    ### if opt == '--ufloor'
    if type(Set.u_floor_per_dim) not in (types.IntType, types.FloatType):
    
        message = "\n'u_floor_per_dim' value ('--ufloor' option) " + \
                "must be a number!"
        raise ValueError, message   
    

    ### if opt == '--vfloor'
    if type(Set.v_floor_per_dim) not in (types.IntType, types.FloatType):
    
        message = "\n'v_floor_per_dim' value ('--vfloor' option) " + \
                "must be a number!"
        raise ValueError, message
    
    
    ### if opt == '--vrad'
    if Set.dimensions_to_radially_average_v not in (0,1,2,3):
    
        message = "\n'radially_average_v' value ('--vrad' option) must be [0-3]!"
        raise ValueError, message
    
    
    ### if opt == '--originpsf'
    if Set.origin_centered_PSF_flag not in (True, False):
    
        message = "\n'origin_centered_output' value ('--originpsf') must be " + \
                "'True' or 'False'!"
        raise ValueError, message
    
    
    ### if opt == '--deriv'
    if Set.derivative_operator not in ('pixel', 'symmetric', 'FC'):
    
        message = "\n'derivative_operator' value ('--deriv' option) must be:\n" + \
                "\t'pixel', 'symmetric', or 'FC'!"
        raise ValueError, message
    
    
    ### if opt == '--lapl'
    if Set.laplacian_operator not in (0,1,2,3):
    
        message = "\n'laplacian_operator' value ('--lapl' option) must be [0-3]!"
        raise ValueError, message
    
    
    ### if opt == '--objtol'
    if type(Set.object_PCG_tolerance) not in (types.IntType, types.FloatType):
    
        message = "\n'object_tolerance' value ('--objtol' option) must be " + \
                "a number!"
        raise ValueError, message
    else:
    
        Set.original_object_PCG_tolerance = Set.object_PCG_tolerance

        
    ### if opt == '--psftol'
    if type(Set.PSF_PCG_tolerance) not in (types.IntType, types.FloatType):
    
        message = "\n'PSF_tolerance' value ('--psftol' option) must be " + \
                "a number!"
        raise ValueError, message
    else:
    
        Set.original_PSF_PCG_tolerance = Set.PSF_PCG_tolerance
            
        
    ### if opt == '--pcgiterarr'
    if type(Set.PCG_iter_array) not in (types.ListType, types.TupleType):
    
        message = "\n'PCG_iter_array'input value ('--pcgiterarr') must be " + \
                "a tuple or list of the form:\n\t(x,x,...) or [x,x,...]"
        raise ValueError, message
    
    
    ### if opt == '--objiterarr'
    if Set.object_PCG_iter_array is None:
    
        Set.object_PCG_iter_array = Set.PCG_iter_array
    elif type(Set.object_PCG_iter_array) not in (types.ListType, 
            types.TupleType):
    
        message = "\n'object_PCG_iter_array'input value ('--objiterarr') \n" + \
                "must be a tuple or list of the form:\n\t" + \
                "(x,x,...) or [x,x,...]"
        raise ValueError, message
    
    
    ### if opt == '--psfiterarr'
    if Set.PSF_PCG_iter_array is None:
    
        Set.PSF_PCG_iter_array = tuple(N.array(Set.PCG_iter_array))
    if type(Set.PSF_PCG_iter_array) not in (types.ListType, types.TupleType):
    
        message = "\n'PSF_PCG_iter_array'input value ('--psfiterarr') \n" + \
                "must be a tuple or list of the form:\n\t" + \
                "(x,x,...) or [x,x,...]"
        raise ValueError, message
    
    
    ### if opt == '--maxiter'
    if type(Set.max_total_PCG_blocks) is not types.IntType:
    
        message = "\n'max_total_PCG_blocks' value ('--maxiter' option) " + \
                "must be an interger!"
        raise ValueError, message
    
    elif Set.max_total_PCG_blocks > len(Set.PCG_iter_array):
    
        print "Note: 'max_total_PCG_blocks' value ('--maxiter' option) " + \
                "reset to match length of input 'PCG_iter_array'"
        Set.max_total_PCG_blocks = len(Set.PCG_iter_array)
    
    
    ### if opt == '--maxseqstops'
    if type(Set.max_sequential_PCG_stops) is not types.IntType:
    
        message = "\n'max_sequential_stops' value ('--maxseqstops' option) " + \
                "must be an integer!"
        raise ValueError, message
    
    
    ### if opt == '--maxupobj'
    if type(Set.max_uphill_object_PCG_steps) is not types.IntType:
    
        message = "\n'max_uphill_object_PCG_steps' value ('--maxupobj' " + \
                "option) must be an integer!"
        raise ValueError, message
    
    
    ### if opt == '--maxuppsf'
    if type(Set.max_uphill_PSF_PCG_steps) is not types.IntType:
    
        message = "\n'max_uphill_PSF_PCG_steps' value ('--maxuppsf' " + \
                "option) must be an integer!"
        raise ValueError, message


    ### if opt == '--risetol'
    if type(Set.rising_tol_ratio) is not types.FloatType:
    
        message = "\n'rising_tol_ratio' value ('--risetol' " + \
                "option) must be a float!"
        raise ValueError, message
    
    
    ### if opt == '--fracconverge'
    if type(Set.Nframes_fractional_convergence) is not types.FloatType:
    
        message = "\n'fraction_Nframes_convergence' value " + \
                "('--fracconverge' option) must be a float!"
        raise ValueError, message
    
    
    ### if opt == '--objccgtol'
    if type(Set.object_CCG_tolerance) is not types.FloatType:
    
        message = "\n'object_CCG_tolerance' value ('--objccgtol' option) " + \
                "must be a float!"
        raise ValueError, message
    
    
    ### if opt == '--psfccgtol'
    if type(Set.PSF_CCG_tolerance) is not types.FloatType:
    
        message = "\n'PSF_CCG_tolerance' value ('--psfccgtol' option) " + \
                "must be a float!"
        raise ValueError, message
    
    
    ### if opt == '--objccgiter'
    if type(Set.max_object_CCG_iter) is not types.IntType:
    
        message = "\n'max_object_CCG_iter' value ('--objccgiter' option) " + \
                "must be an integer!"
        raise ValueError, message
    
    
    ### if opt == '--psfccgiter'
    if type(Set.max_PSF_CCG_iter) is not types.IntType:
    
        message = "\n'max_PSF_CCG_iter' value ('--psfccgiter' option) " + \
                "must be an integer!"
        raise ValueError, message
    
    
    ### if opt == '--ccgxmin'
    if type(Set.xmin_value) is not types.FloatType:
    
        message = "\n'xmin_value' value ('--ccgxmin' option) must be a float!"
        raise ValueError, message
    
    
    ### if opt == '--ccgivec'
    if type(Set.ivec_value) is not types.IntType:
    
        message = "\n'ivec_value' value ('--ccgivec' option) must be " + \
                "an integer!"
        raise ValueError, message
    
    
    
    ### if opt == '--lgrid'
    if hasattr(Set, 'Lgrid'):
    
        if len(Set.Lgrid) != 3:
    
            message = "\n'Lgrid' value ('--Lgrid' option) must be a " + \
                    "tuple of the following format:\n\t" + \
                    "('lambda_object_center', 'multiply_factor', 'exponent')"
            raise ValueError, message
        
        elif len(Set.Lgrid) == 2:
    
            try:
        
                Set.lambda_object_multiply_factor = float(Set.Lgrid[0])
            except ValueError:
            
                message = "\n'Lgrid[0]' value ('--Lgrid' option) for " + \
                        "'lambda_object_multiply_factor' must be a " + \
                        "positive number!"
                print message
        
            if Set.lambda_object_multiply_factor < 0.:
        
                message = "\n'Lgrid[0]' value ('--Lgrid' option) for " + \
                        "'lambda_object_multiply_factor' must be a " + \
                        "positive number!"
                raise ValueError, message
    
            try:
        
                Set.lambda_object_above_below_exp = float(Set.Lgrid[1])
            except ValueError:
            
                message = "\n'Lgrid[1]' value ('--Lgrid' option) for " + \
                        "'lambda_object_above_below_exp' must be a number!"
                print message
    
        else:
    
            if Set.Lgrid[0] == 'a':
        
                Set.orig_lambda_object_center = Set.lambda_object_input
            else:
            
                try:
            
                    Set.orig_lambda_object_center = float(Set.Lgrid[0])
                except ValueError:
            
                    message = "\n'Lgrid[0]' value ('--Lgrid' option) must " + \
                            "be 'a' for adaptive lambda_object grid search " + \
                            "centering or a number!"
                    print message
                
            try:
        
                Set.lambda_object_multiply_factor = float(Set.Lgrid[1])
            except ValueError:
            
                message = "\n'Lgrid[1]' value ('--Lgrid' option) for " + \
                        "'lambda_object_multiply_factor' must be a " + \
                        "positive number!"
                print message
        
            if Set.lambda_object_multiply_factor < 0.:
        
                message = "\n'Lgrid[1]' value ('--Lgrid' option) for " + \
                        "'lambda_object_multiply_factor' must be a " + \
                        "positive number!"
                raise ValueError, message
    
            try:
        
                Set.lambda_object_above_below_exp = float(Set.Lgrid[2])
            except ValueError:
            
                message = "\n'Lgrid[2]' value ('--Lgrid' option) for " + \
                    "'lambda_object_above_below_exp' must be a number!"
                print message
    
    
    ### if opt == '--tgrid'
    if hasattr(Set, 'Tgrid'):
    
        if len(Set.Tgrid) != 3:
    
            message = "\n'Tgrid' value ('--tgrid' option) must be a tuple " + \
                    "of the following format:\n\t('theta_center', " + \
                    "'multiply_factor', 'exponent')"
            raise ValueError, message
        
        elif len(Set.Tgrid) == 2:
    
            try:
        
                Set.theta_multiply_factor = float(Set.Tgrid[0])
            except ValueError:
            
                message = "\n'Tgrid[0]' value ('--Tgrid' option) for " + \
                        "'theta_factor' must be a positive number!"
                print message
        
            if Set.theta_multiply_factor < 0.:
        
                message = "\n'Tgrid[0]' value ('--Tgrid' option) for " + \
                        "'theta_multiply_factor' must be a positive number!"
                raise ValueError, message
    
            try:
        
                Set.theta_above_below_exp = float(Set.Tgrid[1])
            except ValueError:
            
                message = "\n'Tgrid[1]' value ('--Tgrid' option) for " + \
                        "'theta_above_below_exp' must be a number!"
                print message
    
        else:
    
            if Set.Tgrid[0] == 'a':
        
                Set.inv_theta_center = 1./Set.theta_input
            else:
            
                try:
            
                    Set.inv_theta_center = 1./float(Set.Tgrid[0])
                except ValueError:
            
                    message = "\n'Tgrid[0]' value ('--Tgrid' option) must " + \
                            "be 'a' for adaptive theta grid search " + \
                            "centering or a number!"
                    print message
                
            try:
        
                Set.theta_multiply_factor = float(Set.Tgrid[1])
            except ValueError:
            
                message = "\n'Tgrid[1]' value ('--Tgrid' option) for " + \
                        "'theta_multiply_factor' must be a positive number!"
                print message
        
            if Set.theta_multiply_factor < 0.:
        
                message = "\n'Tgrid[1]' value ('--Tgrid' option) for " + \
                        "'theta_multiply_factor' must be a positive number!"
                raise ValueError, message
    
            try:
        
                Set.theta_above_below_exp = float(Set.Tgrid[2])
            except ValueError:
            
                message = "\n'Tgrid[2]' value ('--Tgrid' option) for " + \
                    "'theta_above_below_exp' must be a number!"
                print message
    
####
