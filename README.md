aida
====

Code for the Adaptive Image Deconvolution Algorithm (AIDA) of Hom et al. (2007) J. Opt. Soc. Amer. A

The Adaptive Image Deconvolution Algorithm (AIDA) package was developed mainly by Erik Hom (with invaluable help from Sebastian Haase, creator of Priithon, and Franck Marchis, Planetary Astronomer at SETI Institute & UC-Berkeley) while in John Sedat's lab at UCSF; details are described in J. Opt. Soc. Am. A 24:1580-1600 (2007).

To contribute please read section Information about algorithm and citations

# Getting Started 

In order to run AIDA, you must have some form of Python 3 installed.You must also have the following python libraries: numpy, pillow, astropy, scipy, r-fftw, cython, pip, gcc/libgcc, matplotlib

 ### or

Recommended to use anaconda environment control and import the provided environment.yml 
file using: conda env create -f environment.yml
followed by: conda activate Python36

# Example Images
![alt text](https://raw.githubusercontent.com/harshnagarkar/aida/master/example/jupiter.png? "Jupiter AIDA comparison")

![alt text](https://raw.githubusercontent.com/harshnagarkar/aida/master/example/Uranus%20and%20jupiter.png "Uranus and Jupiter Deconvolved using AIDA")
# Running the Program

### Using AIDA by Command Line

Running through AIDA.py:
Before launching AIDA.py, open AIDA_Settings.py in a text editor. Change any settings to the
value you desire. Only the image filename, and psf filename are required. Other options may be
left as the default values. (For new users, it is recommended not to adjust other settings.) Then
launch AIDA.py and wait for the deconvolution to end.


### Using AIDA by GUI

Running through AIDA_app.py (Recommended):
Simply launch AIDA_app.py to have the AIDA Graphical User Interface appear. You may use
either the text field or the browse buttons to add the image file, psf file, and results directory
name. If wishing to run AIDA with the default settings, then press run and wait for completion.
You will know when the deconvolution is finished when the Graphical User Interface
automatically closes.


# Information About this algorithm:

AIDA is a de novo implementation and extension of the MISTRAL myopic deconvolution method developed by Mugnier et al. (2004) (see J. Opt. Soc. Am. A 21:1841-1854). The MISTRAL approach, used thus far to process astronomical images, has been shown to yield object reconstructions with excellent edge preservation and photometric precision. AIDA improves upon the original MISTRAL implementation in the following ways:

• AIDA is free and intended for further open source development - please help improve it; we welcome any feedback!

• AIDA is written in open source Numerical Python(numpy) instead of the commercial Interactive Data Language (IDL, Research Systems Inc.)

• AIDA typically runs at least 15 times faster

• AIDA includes a scheme to automatically estimate the 2 hyperparameters used to control noise suppression and edge-preserving regularization in the deconvolution process; this frees the user from quite a bit of ad hoc parameter guessing and in practice, speeds up the production of satisfactory reconstructions by an additional order of magnitude

• AIDA can deconvolve multiple frame data and three-dimensional image stacks encountered in adaptive optics and light microscopic imaging.

Although development and maintenance of AIDA has been minimal the past couple of years, this project site has been set up to revive these efforts and nurture community use and development.

With colleagues Franck Marchis and Clement Chalumeau at SETI, AIDA now runs independent of the Priithon platform (AIDA v1.4)

With James Long, we have reinstituted FFTW3 capabilities, which leads to a 30% increase in speed.  This, along with a GUI interface for running the AIDA code will be forthcoming in AIDA v1.5.

There are still bugs we are ferreting out to enable image input of TIFF files, and deconvolution of 3D image data, which was broken in previous code updates.  Please be patient!

# Citations
```
For additional help please read this paper
Hom, Erik F. Y., et al. “AIDA: an Adaptive Image Deconvolution Algorithm with Application to
Multi-Frame and Three-Dimensional Data.” Journal of the Optical Society of America A,
vol. 24, no. 6, 2007, p. 1580., doi:10.1364/josaa.24.001580.
```
