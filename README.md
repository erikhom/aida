aida
====

Code for the Adaptive Image Deconvolution Algorithm (AIDA) of Hom et al. (2007) J. Opt. Soc. Amer. A

The Adaptive Image Deconvolution Algorithm (AIDA) package was developed mainly by Erik Hom (with invaluable help from Sebastian Haase, creator of Priithon, and Franck Marchis, Planetary Astronomer at SETI Institute & UC-Berkeley) while in John Sedat's lab at UCSF; details are described in J. Opt. Soc. Am. A 24:1580-1600 (2007).

AIDA is a de novo implementation and extension of the MISTRAL myopic deconvolution method developed by Mugnier et al. (2004) (see J. Opt. Soc. Am. A 21:1841-1854). The MISTRAL approach, used thus far to process astronomical images, has been shown to yield object reconstructions with excellent edge preservation and photometric precision. AIDA improves upon the original MISTRAL implementation in the following ways:

• AIDA is free and intended for further open source development - please help improve it; we welcome any feedback!

• AIDA is written in open source Numerical Python instead of the commercial Interactive Data Language (IDL, Research Systems Inc.)

• AIDA typically runs at least 15 times faster

• AIDA includes a scheme to automatically estimate the 2 hyperparameters used to control noise suppression and edge-preserving regularization in the deconvolution process; this frees the user from quite a bit of ad hoc parameter guessing and in practice, speeds up the production of satisfactory reconstructions by an additional order of magnitude

• AIDA can deconvolve multiple frame data and three-dimensional image stacks encountered in adaptive optics and light microscopic imaging.

Although development and maintenance of AIDA has been minimal the past couple of years, this project site has been set up to revive these efforts and nurture community use and development.

Stay tuned for improvements as we are able to get to them...
