#!/usr/bin/env priithon
################################################################################
#
#   File:       AIDA_app.py
#
#   Summary:    Python application that runs the Adaptive Image Deconvolution 
#               Algorithm (AIDA)
#
#   Authors:    James Long (Hom Lab, UM)
#
#   Other:      See 'AIDA_version.py' for date and version details
#               See 'LICENSE.txt' for use, license, and reference details
#
################################################################################

from tkinter import filedialog
from tkinter import ttk
from tkinter import *
from PIL import Image, ImageTk
from astropy.io import fits
import AIDA
import AIDA_Settings as Set
import AIDA_Functions as Fun

class Window(Frame):
    def __init__(self, master = None):
        Frame.__init__(self, master)
        self.master = master
        self.init_window()
        
    def init_window(self):
        self.master.title("AIDA")
        
#        quitButton = Button(self, text = "Exit", command = self.showImg)
#        quitButton.place(x=0, y=0)
##Image Input Creation
        self.imageLabel = Label(root, text = "Image File:")
        self.imageLabel.grid(row = 0, column = 0, sticky = 'W')
        self.imageEntry = Entry(root, bd = 3)
        self.imageEntry.grid(row = 1, column = 0)
        self.imageButton = Button(root,text = "Browse...", command= self.image_browse)
        self.imageButton.grid(row = 1, column = 1, sticky = 'W')

##PSF Input Creation
        self.psfLabel = Label(root, text = "PSF File:")
        self.psfLabel.grid(row = 2, column = 0, sticky = 'W')
        self.psfEntry = Entry(root, bd = 3)
        self.psfEntry.grid(row = 3, column = 0)
        self.psfButton = Button(root,text = "Browse...", command = self.psf_browse)
        self.psfButton.grid(row = 3, column = 1, sticky = 'W')
        
##Destination Creation
        self.destLabel = Label(root, text = "Destination Folder:")
        self.destLabel.grid(row = 4, column = 0, sticky = 'W')
        self.destEntry = Entry(root, bd = 3)
        self.destEntry.grid(row = 5, column = 0)
        self.destButton = Button(root,text = "Browse...", command = self.dest_browse)
        self.destButton.grid(row = 5, column = 1, sticky = 'W')
##Decon type
        self.type = StringVar()
        self.type.set('myopic')
        self.typeLabel = Label(root, text = "Deconvolution Type:")
        self.typeLabel.grid(row = 6, column = 0, sticky = 'W')
        self.radioMyopic = Radiobutton(root, text='myopic', state = ACTIVE, variable=self.type, value='myopic')
        self.radioMyopic.grid(row = 7, column = 0, sticky = 'W')
        self.radioClassical = Radiobutton(root, text='classical', variable=self.type, value='classical')
        self.radioClassical.grid(row = 7, column = 0, sticky = 'E')
        self.radioNobjects = Radiobutton(root, text='nobjects', variable=self.type, value='nobjects')
        self.radioNobjects.grid(row = 7, column = 1, sticky = 'W')
        self.radioNpsfs = Radiobutton(root, text="npsfs", variable=self.type, value='npsfs')
        self.radioNpsfs.grid(row = 7, column = 2, sticky = 'W')
        
##Background
        self.backgroundLabel = Label(root, text = "Background Value:")
        self.backgroundLabel.grid(row = 8, column = 0, sticky = 'W')
        self.backgroundEntry = Entry(root, bd = 3, width = 5)
        self.backgroundEntry.insert(END, 0)
        self.backgroundEntry.grid(row = 8, column = 0, sticky = 'E')
##Sigma_det
        self.sigmaLabel = Label(root, text = "Sigma Det:")
        self.sigmaLabel.grid(row = 9, column = 0, sticky = 'W')
        self.sigmaEntry = Entry(root, bd = 3, width = 5)
        self.sigmaEntry.insert(END, 5)
        self.sigmaEntry.grid(row = 9, column = 0, sticky = 'E')
##dark_image
        self.darkLabel = Label(root, text = "Dark Image:")
        self.darkLabel.grid(row = 10, column = 0, sticky = 'W')
        self.darkEntry = Entry(root, bd = 3)
        self.darkEntry.insert(END, 'None')
        self.darkEntry.grid(row = 11, column = 0)
        self.darkButton = Button(root,text = "Browse...", command= self.dark_browse)
        self.darkButton.grid(row = 11, column = 1, sticky = 'W')


        



##Run Button Creation
        self.runButton = Button(root,text = "Run", command = self.runAIDA)
        self.runButton.grid(row = 12, column = 0, sticky = 'W')
        
        self.optionsButton = Button(root,text = "Advanced...", command = self.advancedWindow)
        self.optionsButton.grid(row = 12, column = 1, sticky = 'W')
        

        
    def client_exit(self):
        exit()
        
    def runAIDA(self):
        Set.image_filenames = self.imageEntry.get()
        Set.PSF_filenames = self.psfEntry.get().split(' ')
        Set.results_directory = self.destEntry.get()
        Set.decon_type = self.type.get()
        Set.background = int(self.backgroundEntry.get())
        Set.sigma_det = int(self.sigmaEntry.get())
        if self.darkEntry.get().strip() == 'None':
            Set.dark_image = None;
        else:
            Set.dark_image = self.darkEntry.get()
        
        AIDA.RunAIDA()
        root.destroy()
    
    def update(self):
        line = os.read(pipein, 100)
        self.text.insert(END, line)
        self.text.after(1000, self.update)
        

    def image_browse(self):
        filename = filedialog.askopenfilename(initialdir = "/", title = "Select file", filetypes = (("fits files", "*.fits"),("tiff files", "*.tiff") ), multiple = 'true')
        self.imageEntry.delete(0,END)
        self.imageEntry.insert(0,filename)
        
    def true_image_browse(self):
        filename = filedialog.askopenfilename(initialdir = "/", title = "Select file", filetypes = (("fits files", "*.fits"),("tiff files", "*.tiff") ), multiple = 'true')
        self.trueImageEntry.delete(0,END)
        self.trueImageEntry.insert(0,filename)

    
    def dark_browse(self):
        filename = filedialog.askopenfilename(initialdir = "/", title = "Select file", filetypes = (("fits files", "*.fits"),("tiff files", "*.tiff") ), multiple = 'false')
        self.darkEntry.delete(0,END)
        self.darkEntry.insert(0,filename)


    def psf_browse(self):
        filename = filedialog.askopenfilename(initialdir = "/", title = "Select file", filetypes = (("fits files", "*.fits"),("tiff files", "*.tiff") ), multiple = 'true')
        self.psfEntry.delete(0,END)
        self.psfEntry.insert(0,filename)
        
    def true_psf_browse(self):
        filename = filedialog.askopenfilename(initialdir = "/", title = "Select file", filetypes = (("fits files", "*.fits"),("tiff files", "*.tiff") ), multiple = 'true')
        self.truePSFEntry.delete(0,END)
        self.truePSFEntry.insert(0,filename)

    
    def dest_browse(self):
        filename = filedialog.askdirectory(initialdir="/",title = "Select directory")
        self.destEntry.delete(0,END)
        self.destEntry.insert(0,filename)
        
    def advancedWindow(self):
        self.advWin = Toplevel(root)
        self.tabControl = ttk.Notebook(self.advWin)
        self.tabControl.grid(row =0,column=0)
        tab1 = Frame(self.tabControl)
        self.tabControl.add(tab1,text = 'Options 1')
        tab2 = Frame(self.tabControl)
        self.tabControl.add(tab2, text = 'Options 2')
        tab3 = Frame(self.tabControl)
        self.tabControl.add(tab3, text = 'Options 3')
        tab4 = Frame(self.tabControl)
        self.tabControl.add(tab4, text = 'Options 4')
        tab5 = Frame(self.tabControl)
        self.tabControl.add(tab5, text = 'Options 5')
        tab1.columnconfigure(0, weight=1)        
        tab2.columnconfigure(0, weight=1)
        tab3.columnconfigure(0, weight=1)
        tab4.columnconfigure(0, weight=1)
        tab5.columnconfigure(0, weight=1)



#dimension
        self.dimenLabel = Label(tab1, text = 'Image Dimensions:')
        self.dimenLabel.grid(row = 0, column = 0)
        self.dimenBox = Spinbox(tab1, values = (2,3), textvariable = 2)
        self.dimenBox.grid(row = 0, column = 1)
#zeta
        self.zetaLabel = Label(tab1, text = 'Zeta:')
        self.zetaLabel.grid(row = 1, column = 0)
        self.zetaBox = Spinbox(tab1, from_=0, to=1000, increment = 0.1)
        self.zetaBox.delete(0,END)
        self.zetaBox.insert(END,3.0)
        self.zetaBox.grid(row = 1, column = 1)

#lambda_object_scaling
        self.lamObjLabel = Label(tab1, text = 'Lambda Object Scaling:')
        self.lamObjLabel.grid(row = 2, column = 0)
        self.lamObjBox = Spinbox(tab1, from_=0, to=1000, increment = 0.1)
        self.lamObjBox.delete(0,END)
        self.lamObjBox.insert(END,1.0)
        self.lamObjBox.grid(row = 2, column = 1)

#theta_scaling
        self.thetaLabel = Label(tab1, text = 'Theta Scaling:')
        self.thetaLabel.grid(row = 3, column = 0)
        self.thetaBox = Spinbox(tab1, from_=0, to=1000, increment = 0.1)
        self.thetaBox.delete(0,END)
        self.thetaBox.insert(END,1.0)
        self.thetaBox.grid(row = 3, column = 1)

#lambda_OTF_scaling
        self.lamOtfLabel = Label(tab1, text = 'Lambda OTF Scaling:')
        self.lamOtfLabel.grid(row = 4, column = 0)
        self.lamOtfBox = Spinbox(tab1, from_=0, to=1000, increment = 0.1)
        self.lamOtfBox.delete(0,END)
        self.lamOtfBox.insert(END,1.0)
        self.lamOtfBox.grid(row = 4, column = 1)

#lambda_PSF_scaling
        self.lamPsfLabel = Label(tab1, text = 'Lambda PSF Scaling:')
        self.lamPsfLabel.grid(row = 5, column = 0)
        self.lamPsfBox = Spinbox(tab1, from_=0, to=1000, increment = 0.1)
        self.lamPsfBox.delete(0,END)
        self.lamPsfBox.insert(END,1.0)
        self.lamPsfBox.grid(row = 5, column = 1)
        
#info_level
        self.infoLabel = Label(tab1, text = 'Info Level:')
        self.infoLabel.grid(row = 6, column = 0)
        self.infoBox = Spinbox(tab1, values = (0,1,2))
        self.infoBox.delete(0,END)
        self.infoBox.insert(END,2)
        self.infoBox.grid(row = 6, column = 1)

#memory_usage_level
        self.memoryLabel = Label(tab1, text = 'Memory Usage Level:')
        self.memoryLabel.grid(row = 7, column = 0)
        self.memoryBox = Spinbox(tab1, values = (0,1,2,3))
        self.memoryBox.delete(0,END)
        self.memoryBox.insert(END,3)
        self.memoryBox.grid(row = 7, column = 1)

#terms checkboxes
        
#precision
        self.precision = StringVar()
        self.precision.set("double")
        self.precisionLabel = Label(tab1, text = "Precision:")
        self.precisionLabel.grid(row = 9, column = 0)
        self.radioSingle = Radiobutton(tab1, text="single", state = ACTIVE, variable=self.precision, value="single")
        self.radioSingle.grid(row = 9, column = 1, sticky = 'W')
        self.radioDouble = Radiobutton(tab1, text="double", variable=self.precision, value="double")
        self.radioDouble.grid(row = 9, column = 1, sticky = 'E')

#output_format
        self.outFormat = StringVar()
        self.outFormat.set('None')
        self.outFormatLabel = Label(tab1, text = "Output Format:")
        self.outFormatLabel.grid(row = 10, column = 0)
        self.outFormatNone = Radiobutton(tab1, text="None", state = ACTIVE, variable=self.outFormat, value='None')
        self.outFormatNone.grid(row = 10, column = 1, sticky = 'W')
        self.outFormatF= Radiobutton(tab1, text="f", variable=self.outFormat, value="f")
        self.outFormatF.grid(row = 10, column = 1)
        self.outFormatM= Radiobutton(tab1, text="m", variable=self.outFormat, value="m")
        self.outFormatM.grid(row = 10, column = 1, sticky = 'E')
        
#output_intermediate_files_level
        self.outLevelLabel = Label(tab1, text = 'Output IM Files Level:')
        self.outLevelLabel.grid(row = 11, column = 0)
        self.outLevelBox = Spinbox(tab1, from_=0, to=1000)
        self.outLevelBox.delete(0,END)
        self.outLevelBox.insert(END,0)
        self.outLevelBox.grid(row = 11, column = 1)

#resizeimage_flag = True                # default = True
        self.resizeImage = BooleanVar()
        self.resizeImage.set(True)
        self.resizeImageLabel = Label(tab2, text = "Resize Image:")
        self.resizeImageLabel.grid(row = 0, column = 0)
        self.radioresizeImageTrue = Radiobutton(tab2, text="True", state = ACTIVE, variable=self.resizeImage, value=True)
        self.radioresizeImageTrue.grid(row = 0, column = 1, sticky = 'W')
        self.radioresizeImageFalse = Radiobutton(tab2, text="False", variable=self.resizeImage, value=False)
        self.radioresizeImageFalse.grid(row = 0, column = 1, sticky = 'E')

        
#zresize_flag = False                   # default = False
        self.zresize = BooleanVar()
        self.zresize.set(False)
        self.zresizeLabel = Label(tab2, text = "Z Resize:")
        self.zresizeLabel.grid(row = 1, column = 0)
        self.radiozresizeTrue = Radiobutton(tab2, text="True", variable=self.zresize, value=True)
        self.radiozresizeTrue.grid(row = 1, column = 1, sticky = 'W')
        self.radiozresizeFalse = Radiobutton(tab2, text="False", state = ACTIVE, variable=self.zresize, value=False)
        self.radiozresizeFalse.grid(row = 1, column = 1, sticky = 'E')

        
#initial_object_guess = 'zero'          # default = 'zero'  ['image', 'wiener']
        self.initGuess = StringVar()
        self.initGuess.set('zero')
        self.initGuessLabel = Label(tab2, text = "Initial Object Guess:")
        self.initGuessLabel.grid(row = 2, column = 0)
        self.initGuessZero = Radiobutton(tab2, text="zero", state = ACTIVE, variable=self.initGuess, value='zero')
        self.initGuessZero.grid(row = 2, column = 1, sticky = 'W')
        self.initGuessImage= Radiobutton(tab2, text="image", variable=self.initGuess, value="image")
        self.initGuessImage.grid(row = 2, column = 1)
        self.initGuessWiener= Radiobutton(tab2, text="wiener", variable=self.initGuess, value="wiener")
        self.initGuessWiener.grid(row = 2, column = 1, sticky = 'E')

#exclude_radius = 5                     # default = 5
        self.excludeRadiusLabel = Label(tab2, text = 'Exclude Radius:')
        self.excludeRadiusLabel.grid(row = 3, column = 0)
        self.excludeRadiusBox = Spinbox(tab2, from_=0, to=1000)
        self.excludeRadiusBox.delete(0,END)
        self.excludeRadiusBox.insert(END,5)
        self.excludeRadiusBox.grid(row = 3, column = 1)

#cleaning = (1,1,1,1)                   # default = (1,1,1,1)
        
#PSF_centering = 'unknown'              # default = 'unknown' ['origin', 'array_center']
        self.psfCentering = StringVar()
        self.psfCentering.set('unknown')
        self.psfCenteringLabel = Label(tab2, text = "PSF Centering:")
        self.psfCenteringLabel.grid(row = 5, column = 0)
        self.psfCenteringUnknown = Radiobutton(tab2, text="unknown", state = ACTIVE, variable=self.psfCentering, value='unknown')
        self.psfCenteringUnknown.grid(row = 5, column = 1, sticky = 'W')
        self.psfCenteringOrigin= Radiobutton(tab2, text="origin", variable=self.psfCentering, value="origin")
        self.psfCenteringOrigin.grid(row = 5, column = 1, sticky = 'E')
        self.psfCenteringArray= Radiobutton(tab2, text="array center", variable=self.psfCentering, value="array_center")
        self.psfCenteringArray.grid(row = 5, column = 2, sticky = 'W')

        
#PSF_subtract_background_percent = 0.   # default = 0.
        self.subtractBackLabel = Label(tab2, text = 'PSF Subtract Background %:')
        self.subtractBackLabel.grid(row = 6, column = 0)
        self.subtractBackBox = Spinbox(tab2, from_=0, to=100)
        self.subtractBackBox.delete(0,END)
        self.subtractBackBox.insert(END,0.0)
        self.subtractBackBox.grid(row = 6, column = 1)

#nsigmas = 2.                           # default = 2
        self.NsigmaLabel = Label(tab2, text = 'N Sigmas:')
        self.NsigmaLabel.grid(row = 7, column = 0)
        self.NsigmaBox = Spinbox(tab2, from_=0, to=1000)
        self.NsigmaBox.delete(0,END)
        self.NsigmaBox.insert(END,2)
        self.NsigmaBox.grid(row = 7, column = 1)

#PSF_threshold_percent = 1e-7           # default = 1e-7
        self.psfThreshLabel = Label(tab2, text = 'PSF Threshold %:')
        self.psfThreshLabel.grid(row = 8, column = 0)
        self.psfThreshEntry = Entry(tab2)
        self.psfThreshEntry.delete(0,END)
        self.psfThreshEntry.insert(END,1e-7)
        self.psfThreshEntry.grid(row = 8, column = 1)
        
#OTF_threshold_percent = 1e-7           # default = 1e-7
        self.otfThreshLabel = Label(tab2, text = 'OTF Threshold %:')
        self.otfThreshLabel.grid(row = 9, column = 0)
        self.otfThreshEntry = Entry(tab2)
        self.otfThreshEntry.delete(0,END)
        self.otfThreshEntry.insert(END,1e-7)
        self.otfThreshEntry.grid(row = 9, column = 1)

#fill = 0.                              # default = 0.
        
#u_floor_per_dim = 1e-6                 # default = 1e-6
        self.uFloorLabel = Label(tab2, text = 'UFloor/dim:')
        self.uFloorLabel.grid(row = 11, column = 0)
        self.uFloorEntry = Entry(tab2)
        self.uFloorEntry.delete(0,END)
        self.uFloorEntry.insert(END,1e-6)
        self.uFloorEntry.grid(row = 11, column = 1)
        
#v_floor_per_dim = 1e-6                 # default = 1e-6
        self.vFloorLabel = Label(tab2, text = 'VFloor/dim:')
        self.vFloorLabel.grid(row = 12, column = 0)
        self.vFloorEntry = Entry(tab2)
        self.vFloorEntry.delete(0,END)
        self.vFloorEntry.insert(END,1e-6)
        self.vFloorEntry.grid(row = 12, column = 1)

#v_threshold = 0.01                     # default = 0.01
        self.vThreshLabel = Label(tab3, text = 'V Threshold:')
        self.vThreshLabel.grid(row = 0, column = 0)
        self.vThreshBox = Spinbox(tab3, from_=0, to=1000, increment = 0.01)
        self.vThreshBox.delete(0,END)
        self.vThreshBox.insert(END,0.01)
        self.vThreshBox.grid(row = 0, column = 1)

#dimensions_to_radially_average_v = 0   # default = 0  [2, 3]
        self.d2RadLabel = Label(tab3, text = 'Dimensions to Radially Avg V:')
        self.d2RadLabel.grid(row = 1, column = 0)
        self.d2RadBox = Spinbox(tab3, values = (0,2,3))
        self.d2RadBox.delete(0,END)
        self.d2RadBox.insert(END,0)
        self.d2RadBox.grid(row = 1, column = 1)

#origin_centered_PSF_flag = False       # default = False
        self.originCenter = BooleanVar()
        self.originCenter.set(False)
        self.originCenterLabel = Label(tab3, text = "Origin Centered PSF:")
        self.originCenterLabel.grid(row = 2, column = 0)
        self.radiooriginCenterTrue = Radiobutton(tab3, text="True", variable=self.originCenter, value=True)
        self.radiooriginCenterTrue.grid(row = 2, column = 1, sticky = 'W')
        self.radiooriginCenterFalse = Radiobutton(tab3, text="False", state = ACTIVE, variable=self.originCenter, value=False)
        self.radiooriginCenterFalse.grid(row = 2, column = 1, sticky = 'E')
        
#derivative_operator = 'FC'             # default = 'FC  ['pixel', 'symmetric']
        self.derivativeOp = StringVar()
        self.derivativeOp.set('FC')
        self.derivativeOpLabel = Label(tab3, text = "Derivative Operator:")
        self.derivativeOpLabel.grid(row = 3, column = 0)
        self.derivativeOpUnknown = Radiobutton(tab3, text="FC", state = ACTIVE, variable=self.derivativeOp, value='FC')
        self.derivativeOpUnknown.grid(row = 3, column = 1, sticky = 'W')
        self.derivativeOpOrigin= Radiobutton(tab3, text="pixel", variable=self.derivativeOp, value='pixel')
        self.derivativeOpOrigin.grid(row = 3, column = 1)
        self.derivativeOpArray= Radiobutton(tab3, text="symmetric", variable=self.derivativeOp, value='symmetric')
        self.derivativeOpArray.grid(row = 3, column = 2, sticky = 'W')

#laplacian_operator = 3                 # default = 3  [0, 1, 2]
        self.laplacianLabel = Label(tab3, text = 'Laplacian Operator:')
        self.laplacianLabel.grid(row = 4, column = 0)
        self.laplacianBox = Spinbox(tab3, values = (0,1,2,3))
        self.laplacianBox.delete(0,END)
        self.laplacianBox.insert(END,3)
        self.laplacianBox.grid(row = 4, column = 1)

#rho = 1.                               # default = 1.
        self.rhoLabel = Label(tab3, text = 'Rho:')
        self.rhoLabel.grid(row = 5, column = 0)
        self.rhoBox = Spinbox(tab3, from_=0, to=1000)
        self.rhoBox.delete(0,END)
        self.rhoBox.insert(END,1.0)
        self.rhoBox.grid(row = 5, column = 1)

#object_PCG_tolerance = 0.1             # default = 0.1
        self.objToleranceLabel = Label(tab3, text = 'Object PCG Tolerance:')
        self.objToleranceLabel.grid(row = 6, column = 0)
        self.objToleranceBox = Spinbox(tab3, from_=0, to=1000, increment = 0.01)
        self.objToleranceBox.delete(0,END)
        self.objToleranceBox.insert(END,0.1)
        self.objToleranceBox.grid(row = 6, column = 1)

#PSF_PCG_tolerance = 0.1                # default = 0.1
        self.psfToleranceLabel = Label(tab3, text = 'PSF PCG Tolerance:')
        self.psfToleranceLabel.grid(row = 7, column = 0)
        self.psfToleranceBox = Spinbox(tab3, from_=0, to=1000, increment = 0.01)
        self.psfToleranceBox.delete(0,END)
        self.psfToleranceBox.insert(END,0.1)
        self.psfToleranceBox.grid(row = 7, column = 1)

#PCG_iter_array = (1,3,5,7,9,11,13,15)  # default = (1,3,5,7,9,11,13,15)
        self.pcgArrayLabel = Label(tab3, text = 'PCG Iteration Array:')
        self.pcgArrayLabel.grid(row = 8, column = 0)
        self.pcgArrayEntry = Entry(tab3)
        self.pcgArrayEntry.delete(0,END)
        self.pcgArrayEntry.insert(END,(1,3,5,7,9,11,13,15))
        self.pcgArrayEntry.grid(row = 8, column = 1)

#object_PCG_iter_array = None           # default = None; set to PCG_iter_array
        self.objPCG = BooleanVar()
        self.objPCG.set(False)
        self.objPCGLabel = Label(tab3, text = "Object PCG Arry:")
        self.objPCGLabel.grid(row = 9, column = 0)
        self.radioobjPCGTrue = Radiobutton(tab3, text="Set to PCG Iteration", variable=self.objPCG, value=True)
        self.radioobjPCGTrue.grid(row = 9, column = 1, sticky = 'W')
        self.radioobjPCGFalse = Radiobutton(tab3, text="None", state = ACTIVE, variable=self.objPCG, value=False)
        self.radioobjPCGFalse.grid(row = 9, column = 2, sticky = 'W')

#PSF_PCG_iter_array = None              # default = None; set to PCG_iter_array
        self.psfPCG = BooleanVar()
        self.psfPCG.set(False)
        self.psfPCGLabel = Label(tab3, text = "PSF PCG Array:")
        self.psfPCGLabel.grid(row = 10, column = 0)
        self.radiopsfPCGTrue = Radiobutton(tab3, text="Set to PCG Iteration", variable=self.psfPCG, value=True)
        self.radiopsfPCGTrue.grid(row = 10, column = 1, sticky = 'W')
        self.radiopsfPCGFalse = Radiobutton(tab3, text="None", state = ACTIVE, variable=self.psfPCG, value=False)
        self.radiopsfPCGFalse.grid(row = 10, column = 2, sticky = 'W')

#max_classical_PCG_iter = (24,)                 # default = (24,)
        
        
        
#max_total_PCG_blocks = len(Set.PCG_iter_array) # default = len(PCG_iter_array)
        
        
        
#max_sequential_PCG_stops = 2                   # default = 2
        
        self.maxSeqPCGStopsLabel = Label(tab4, text = 'Max Sequential PCG Stops:')
        self.maxSeqPCGStopsLabel.grid(row = 0, column = 0)
        self.maxSeqPCGStopsBox = Spinbox(tab4, from_=0, to=1000)
        self.maxSeqPCGStopsBox.delete(0,END)
        self.maxSeqPCGStopsBox.insert(END,2)
        self.maxSeqPCGStopsBox.grid(row = 0, column = 1)
        
#max_uphill_object_PCG_steps = 2                # default = 2
        
        self.maxUpObjPCGStepsLabel = Label(tab4, text = 'Max Uphill Object PCG Steps:')
        self.maxUpObjPCGStepsLabel.grid(row = 1, column = 0)
        self.maxUpObjPCGStepsBox = Spinbox(tab4, from_=0, to=1000)
        self.maxUpObjPCGStepsBox.delete(0,END)
        self.maxUpObjPCGStepsBox.insert(END,2)
        self.maxUpObjPCGStepsBox.grid(row = 1, column = 1)
        
#max_uphill_PSF_PCG_steps = 2                   # default = 2
        
        self.maxUpPSFPCGStepsLabel = Label(tab4, text = 'Max Uphill PSF PCG Steps:')
        self.maxUpPSFPCGStepsLabel.grid(row = 2, column = 0)
        self.maxUpPSFPCGStepsBox = Spinbox(tab4, from_=0, to=1000)
        self.maxUpPSFPCGStepsBox.delete(0,END)
        self.maxUpPSFPCGStepsBox.insert(END,2)
        self.maxUpPSFPCGStepsBox.grid(row = 2, column = 1)
        
#rising_tol_ratio = 1.03                        # default = 1.03
        
        self.risingTolRatLabel = Label(tab4, text = 'Rising Tolerance Ratio:')
        self.risingTolRatLabel.grid(row = 3, column = 0)
        self.risingTolRatBox = Spinbox(tab4, from_=0, to=1000, increment = 0.01)
        self.risingTolRatBox.delete(0,END)
        self.risingTolRatBox.insert(END,1.03)
        self.risingTolRatBox.grid(row = 3, column = 1)
        
#max_optimization_stops = 3                     # default = 3
        
        self.maxOptimStopsLabel = Label(tab4, text = 'Max Optimization Stops:')
        self.maxOptimStopsLabel.grid(row = 4, column = 0)
        self.maxOptimStopsBox = Spinbox(tab4, from_=0, to=1000)
        self.maxOptimStopsBox.delete(0,END)
        self.maxOptimStopsBox.insert(END,3)
        self.maxOptimStopsBox.grid(row = 4, column = 1)
        
#max_rising_stops = 3                           # default = 3
        
        self.maxRisingStopsLabel = Label(tab4, text = 'Max Rising Stops:')
        self.maxRisingStopsLabel.grid(row = 5, column = 0)
        self.maxRisingStopsBox = Spinbox(tab4, from_=0, to=1000)
        self.maxRisingStopsBox.delete(0,END)
        self.maxRisingStopsBox.insert(END,3)
        self.maxRisingStopsBox.grid(row = 5, column = 1)
        
#Nframes_fractional_convergence = 0.99          # default = 0.99
        
        self.NFramesFractConvLabel = Label(tab4, text = 'NFrames Fractional Convergence:')
        self.NFramesFractConvLabel.grid(row = 6, column = 0)
        self.NFramesFractConvBox = Spinbox(tab4, from_=0, to=1000, increment = 0.01)
        self.NFramesFractConvBox.delete(0,END)
        self.NFramesFractConvBox.insert(END,0.99)
        self.NFramesFractConvBox.grid(row = 6, column = 1)
        
#object_CCG_tolerance = 1e-7        # default = 1e-7
        
        self.objCCGTolerLabel = Label(tab4, text = 'Object CCG Tolerance:')
        self.objCCGTolerLabel.grid(row = 7, column = 0)
        self.objCCGTolerEntry = Entry(tab4)
        self.objCCGTolerEntry.delete(0,END)
        self.objCCGTolerEntry.insert(END,1e-7)
        self.objCCGTolerEntry.grid(row = 7, column = 1)
        
#PSF_CCG_tolerance = 1e-7           # default = 1e-7
        
        self.PSFCCGTolerLabel = Label(tab4, text = 'PSF CCG Tolerance:')
        self.PSFCCGTolerLabel.grid(row = 8, column = 0)
        self.PSFCCGTolerEntry = Entry(tab4)
        self.PSFCCGTolerEntry.delete(0,END)
        self.PSFCCGTolerEntry.insert(END,1e-7)
        self.PSFCCGTolerEntry.grid(row = 8, column = 1)
        
#max_object_CCG_iter = 24           # default = 24
        
        self.maxObjCCGIterLabel = Label(tab4, text = 'Max Object CCG Iteration:')
        self.maxObjCCGIterLabel.grid(row = 9, column = 0)
        self.maxObjCCGIterBox = Spinbox(tab4, from_=0, to=1000)
        self.maxObjCCGIterBox.delete(0,END)
        self.maxObjCCGIterBox.insert(END,24)
        self.maxObjCCGIterBox.grid(row = 9, column = 1)

        
#max_PSF_CCG_iter = 24              # default = 24
        
        self.maxPSFCCGIterLabel = Label(tab4, text = 'Max PSF CCG Iteration:')
        self.maxPSFCCGIterLabel.grid(row = 10, column = 0)
        self.maxPSFCCGIterBox = Spinbox(tab4, from_=0, to=1000)
        self.maxPSFCCGIterBox.delete(0,END)
        self.maxPSFCCGIterBox.insert(END,24)
        self.maxPSFCCGIterBox.grid(row = 10, column = 1)
        
#xmin_value = 0.                    # default = 0.
        
        self.xMinLabel = Label(tab4, text = 'Min X Value:')
        self.xMinLabel.grid(row = 11, column = 0)
        self.xMinBox = Spinbox(tab4)
        self.xMinBox.delete(0,END)
        self.xMinBox.insert(END,0.0)
        self.xMinBox.grid(row = 11, column = 1)


#lambda_object_input = None     # default = None positive float
        
        self.lambdaObjInputLabel = Label(tab4, text = 'Lambda Object Input:')
        self.lambdaObjInputLabel.grid(row = 12, column = 0)
        self.lambdaObjInputEntry = Entry(tab4)
        self.lambdaObjInputEntry.grid(row = 12, column = 1)

        
#theta_input = None             # default = None positive float
        
        self.thetaInputLabel = Label(tab5, text = 'Theta Input:')
        self.thetaInputLabel.grid(row = 0, column = 0)
        self.thetaInputEntry = Entry(tab5)
        self.thetaInputEntry.grid(row = 0, column = 1)
        
#lambda_OTF_input = None        # default = None; set to 1/N_total_pixels
        
        self.lambdaOTFInputLabel = Label(tab5, text = 'Lambda OTF Input:')
        self.lambdaOTFInputLabel.grid(row = 1, column = 0)
        self.lambdaOTFInputEntry = Entry(tab5)
        self.lambdaOTFInputEntry.grid(row = 1, column = 1)
        
#lambda_PSF_input = None        # default = None; set to 1.
        
        self.lambdaPSFInputLabel = Label(tab5, text = 'Lambda PSF Input:')
        self.lambdaPSFInputLabel.grid(row = 2, column = 0)
        self.lambdaPSFInputEntry = Entry(tab5)
        self.lambdaPSFInputEntry.grid(row = 2, column = 1)

#orig_lambda_object_center = 0.01           # default = 0.01
        
        self.origLambdaObjCenLabel = Label(tab5, text = 'Orig Lambda Object Center:')
        self.origLambdaObjCenLabel.grid(row = 3, column = 0)
        self.origLambdaObjCenBox = Spinbox(tab5, from_=0, to=1000, increment = 0.01)
        self.origLambdaObjCenBox.delete(0,END)
        self.origLambdaObjCenBox.insert(END,0.01)
        self.origLambdaObjCenBox.grid(row = 3, column = 1)
        
#lambda_object_multiply_factor = 2.         # default = 2.
        
        self.lamObjMultFactorLabel = Label(tab5, text = 'Lambda Object Mult Factor:')
        self.lamObjMultFactorLabel.grid(row = 4, column = 0)
        self.lamObjMultFactorBox = Spinbox(tab5, from_=0, to=1000, increment = 0.1)
        self.lamObjMultFactorBox.delete(0,END)
        self.lamObjMultFactorBox.insert(END,2.0)
        self.lamObjMultFactorBox.grid(row = 4, column = 1)
        
#lambda_object_above_below_exp = 0.         # default = 0.; i.e. no grid search
        
        self.lamObjAboBelExpLabel = Label(tab5, text = 'Lambda Object Above Below Exp:')
        self.lamObjAboBelExpLabel.grid(row = 5, column = 0)
        self.lamObjAboBelExpBox = Spinbox(tab5, increment = 0.1)
        self.lamObjAboBelExpBox.delete(0,END)
        self.lamObjAboBelExpBox.insert(END,0.0)
        self.lamObjAboBelExpBox.grid(row = 5, column = 1)
        
#theta_multiply_factor = 2.                 # default = 2.
        
        self.thetaMultFactorLabel = Label(tab5, text = 'Theta Mult Factor:')
        self.thetaMultFactorLabel.grid(row = 6, column = 0)
        self.thetaMultFactorBox = Spinbox(tab5, from_=0, to=1000, increment = 0.1)
        self.thetaMultFactorBox.delete(0,END)
        self.thetaMultFactorBox.insert(END,2.0)
        self.thetaMultFactorBox.grid(row = 6, column = 1)
        
#theta_above_below_exp = 0.                 # default = 0.
        
        self.thetaAboBelExpLabel = Label(tab5, text = 'Theta Above Below Exp:')
        self.thetaAboBelExpLabel.grid(row = 7, column = 0)
        self.thetaAboBelExpBox = Spinbox(tab5, increment = 0.1)
        self.thetaAboBelExpBox.delete(0,END)
        self.thetaAboBelExpBox.insert(END,0.0)
        self.thetaAboBelExpBox.grid(row = 7, column = 1)
        
#lambda_OTF_multiply_factor = 2.            # default = 2.
        
        self.lamOTFMultFactorLabel = Label(tab5, text = 'Lambda OTF Mult Factor:')
        self.lamOTFMultFactorLabel.grid(row = 8, column = 0)
        self.lamOTFMultFactorBox = Spinbox(tab5, from_=0, to=1000, increment = 0.1)
        self.lamOTFMultFactorBox.delete(0,END)
        self.lamOTFMultFactorBox.insert(END,2.0)
        self.lamOTFMultFactorBox.grid(row = 8, column = 1)
        
#lambda_OTF_above_below_exp = 0.            # default = 0.; i.e. no grid search
        
        self.lamOTFAboBelExpLabel = Label(tab5, text = 'Lambda OTF Above Below Exp:')
        self.lamOTFAboBelExpLabel.grid(row = 9, column = 0)
        self.lamOTFAboBelExpBox = Spinbox(tab5, increment = 0.1)
        self.lamOTFAboBelExpBox.delete(0,END)
        self.lamOTFAboBelExpBox.insert(END,0.0)
        self.lamOTFAboBelExpBox.grid(row = 9, column = 1)
        
#grid_search_lambda_object_estimate_flag = False    # default = False
        
        self.gridSearchLamObjEst = BooleanVar()
        self.gridSearchLamObjEst.set(False)
        self.gridSearchLamObjEstLabel = Label(tab5, text = "Grid Search Lambda Object Est:")
        self.gridSearchLamObjEstLabel.grid(row = 10, column = 0)
        self.radiogridSearchLamObjEstTrue = Radiobutton(tab5, text="True", variable=self.gridSearchLamObjEst, value=True)
        self.radiogridSearchLamObjEstTrue.grid(row = 10, column = 1, sticky = 'W')
        self.radiogridSearchLamObjEstFalse = Radiobutton(tab5, text="False", state = ACTIVE, variable=self.gridSearchLamObjEst, value=False)
        self.radiogridSearchLamObjEstFalse.grid(row = 10, column = 1, sticky = 'E')

#band_limited_constraint = False            # default = False
        
        self.bandLimitConstraint = BooleanVar()
        self.bandLimitConstraint.set(False)
        self.bandLimitConstraintLabel = Label(tab5, text = "Band Limited Constraint:")
        self.bandLimitConstraintLabel.grid(row = 11, column = 0)
        self.radiobandLimitConstraintTrue = Radiobutton(tab5, text="True", variable=self.bandLimitConstraint, value=True)
        self.radiobandLimitConstraintTrue.grid(row = 11, column = 1, sticky = 'W')
        self.radiobandLimitConstraintFalse = Radiobutton(tab5, text="False", state = ACTIVE, variable=self.bandLimitConstraint, value=False)
        self.radiobandLimitConstraintFalse.grid(row = 11, column = 1, sticky = 'E')
        
#true_object_file = None                    # default = None    
        
        self.trueImageLabel = Label(tab5, text = "True Image File:")
        self.trueImageLabel.grid(row = 12, column = 0, sticky = 'W')
        self.trueImageEntry = Entry(tab5)
        self.trueImageEntry.grid(row = 12, column = 1, sticky = 'W')
        self.trueImageButton = Button(tab5,text = "Browse...", command= self.true_image_browse)
        self.trueImageButton.grid(row = 12, column = 2, sticky = 'W')
        
#true_PSF_file = None                       # default = None
        
        self.truePSFLabel = Label(tab5, text = "True PSF File:")
        self.truePSFLabel.grid(row = 13, column = 0, sticky = 'W')
        self.truePSFEntry = Entry(tab5)
        self.truePSFEntry.grid(row = 13, column = 1, sticky = 'W')
        self.truePSFButton = Button(tab5, text = "Browse...", command= self.true_psf_browse)
        self.truePSFButton.grid(row = 13, column = 2, sticky = 'W')

#debug = False             # default = False

        self.debug = BooleanVar()
        self.debug.set(False)
        self.debugLabel = Label(tab5, text = "Resize Image:")
        self.debugLabel.grid(row = 14, column = 0)
        self.radiodebugTrue = Radiobutton(tab5, text="True", variable=self.debug, value=True)
        self.radiodebugTrue.grid(row = 14, column = 1, sticky = 'W')
        self.radiodebugFalse = Radiobutton(tab5, text="False", state = ACTIVE, variable=self.debug, value=False)
        self.radiodebugFalse.grid(row = 14, column = 1, sticky = 'E')
        
#End of Tab Widgets
        
        
        
        
        
        
        
        
        
        
        
        self.tabControl.select(tab1)
        self.tabControl.enable_traversal()
        
        self.confirmButton = Button(self.advWin, text ="Confirm", command = self.cancel)
        self.confirmButton.grid(row = 1, column = 0, sticky = 'W')
        self.cancelButton = Button(self.advWin, text="Cancel", command = self.cancel)
        self.cancelButton.grid(row = 1, column = 0, sticky = 'E')
        
    def cancel(self):
        self.advWin.destroy()
        
    def confirm(self):
        Set.dimension = int(self.dimenBox.get())
        Set.zeta = float(self.zetaBox.get())
        Set.lambda_object_scaling = float(self.lamObjBox.get())
        Set.theta_scaling = float(self.thetaBox.get())
        Set.lambda_OTF_scaling = float(self.lamOtfBox.get())
        Set.lambda_PSF_scaling = float(self.lamPsfBox.get())
        Set.info_level = int(self.infoBox.get())
        Set.memory_usage_level = int(self.memoryBox.get())
        
        Set.precision = self.precision.get()
        if self.outFormat == 'None':
            Set.output_format = None;
        else:
            Set.output_format = self.outFormat.get()
        Set.output_intermediate_files_level = int(self.outLevelBox.get())
        Set.resizeimage_flag = self.resizeImage.get()
        Set.zresize_flag = self.zresize.get()
        Set.initial_object_guess = self.initGuess.get()
        Set.exclude_radius = int(self.excludeRadiusBox.get())
        
        Set.PSF_centering = self.psfCentering.get()
        Set.PSF_subtract_background_percent = float(self.subtractBackBox.get())
        Set.nsigmas = float(self.NsigmaBox.get())
        Set.PSF_threshold_percent = float(self.psfThreshEntry.get())
        Set.OTF_threshold_percent = float(self.otfThreshEntry.get())
        
        Set.u_floor_per_dim = float(self.uFloorEntry.get())
        Set.v_floor_per_dim = float(self.vFloorEntry.get())
        Set.v_threshold = float(self.vThreshBox.get())
        Set.dimensions_to_radially_average_v = int(self.d2RadBox.get())
        Set.origin_centered_PSF_flag = self.originCenter.get()
        Set.derivative_operator = self.derivativeOp.get()
        Set.laplacian_operator = int(self.laplacianBox.get())
        Set.rho = float(self.rhoBox.get())
        Set.object_PCG_tolerance = float(self.objToleranceBox.get())
        Set.PSF_PCG_tolerance = float(self.psfToleranceBox.get())
        Set.PCG_iter_array = self.pcgArrayEntry.get()
        if self.objPCG:
            Set.object_PCG_iter_array = Set.PCG_iter_array;
        else:
            Set.object_PCG_iter_array = None
        if self.psfPCG:
            Set.PSF_PCG_iter_array = Set.PCG_iter_array;
        else:
            Set.PSF_PCG_iter_array = None
        Set.max_sequential_PCG_stops = int(self.maxSeqPCGStopsBox.get())
        Set.max_uphill_object_PCG_steps = int(self.maxUpObjPCGStepsBox.get())
        Set.max_uphill_PSF_PCG_steps = int(self.maxUpPSFPCGStepsBox.get())
        Set.rising_tol_ratio = float(self.risingTolRatBox.get())
        Set.max_optimization_stops = int(self.maxOptimStopsBox.get())
        Set.max_rising_stops = int(self.maxRisingStopsBox.get())
        Set.Nframes_fractional_convergence = float(self.NFramesFractConvBox.get())
        Set.object_CCG_tolerance = float(self.objCCGTolerEntry.get())
        Set.PSF_CCG_tolerance = float(self.PSFCCGTolerEntry.get())
        Set.max_object_CCG_iter = int(self.maxPSFCCGIterBox.get())
        Set.xmin_value = float(self.xMinBox.get())
        
        if self.lambdaObjInputEntry.get() == '':
            Set.lambda_object_input = None;
        else:
            Set.lambda_object_input = float(self.lambdaObjInputEntry.get())
            
        if self.thetaInputEntry.get() == '':
            Set.theta_input = None;
        else:
            Set.theta_input = float(self.thetaInputEntry.get())
            
        if self.lambdaOTFInputEntry.get() == '':
            Set.lambda_OTF_input = None;
        else:
            Set.lambda_OTF_input = float(self.lambdaOTFInputEntry.get())
            
        if self.lambdaPSFInputEntry.get() == '':
            Set.lambda_PSF_input = None;
        else:
            Set.lambda_PSF_input = float(self.lambdaPSFInputEntry.get())
        
        Set.orig_lambda_object_center = float(self.origLambdaObjCenBox.get())
        Set.lambda_object_multiply_factor = float(self.lamObjMultFactorBox.get())
        Set.lambda_object_above_below_exp = float(self.lamObjAboBelExpBox.get())
        Set.theta_multiply_factor = float(self.thetaMultFactorBox.get())
        Set.theta_above_below_exp = float(self.thetaAboBelExpBox.get())
        Set.lambda_OTF_multiply_factor = float(self.lamOTFAboBelExpBox.get())
        Set.lambda_OTF_above_below_exp = float(self.lamOTFAboBelExpBox.get())
        Set.grid_search_lambda_object_estimate_flag = self.gridSearchLamObjEst.get()
        Set.band_limited_constraint = self.bandLimitConstraint.get()
        
        if self.trueImageEntry.get() == '':
            Set.true_object_file = None;
        else:
            Set.true_object_file = self.trueImageEntry.get()
            
        if self.truePSFEntry.get() == '':
            Set.true_PSF_file = None;
        else:
            Set.true_PSF_file = self.truePSFEntry.get()
        
        Set.debug = self.debug.get()
        self.advWin.destroy()










root = Tk()
root.resizable(False,False)
app = Window(root)
root.mainloop()

#top = Tk(className="AIDA")
#
#
#images = ""
#psf = ""
#
#top.mainloop()