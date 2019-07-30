################################################################################
#
#   File:       AIDA_Version.py
#
#   Summary:    Provides date and vesion information for the Adaptive Image 
#               Deconvolution Algorithm (AIDA) package (AIDAP)
#
#   Authors:    Erik F.Y. Hom (Sedat Lab, UCSF)
#
#   Other:      See 'AIDA_version.py' for date and version details
#               See 'LICENSE.txt' for use, license, and reference details
#
################################################################################


version = '1.4.1'
date = 'September 15, 2016'

# Priithon has been embedded and is no longer required

header = """
================================================================================
         <<< Adaptive Image Deconvolution Algorithm (AIDA) Package >>>
                                       by
    E.F.Y. Hom, F. Marchis, C. Chalumeau, T.K. Lee, S. Haase, D.A. Agard, and J.W. Sedat
                 Copyright (c) 2006, University of California
                 
     ("AIDA: An Adaptive Image Deconvolution Algorithm with Application to
   Multi-Frame and Three-Dimensional Data"  J. Opt. Soc. Amer. A., 24(6):1580-1600)
   
================================================================================
\n\t\tAIDA version: """ + str(version) + "  [ " + date + " ]"
