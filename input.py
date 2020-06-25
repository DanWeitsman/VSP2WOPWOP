#%%        VSP2WOPWOP User Input File

#   Author: Daniel Weitsman
'''''
This input file is used to set up the cases for which to generate the patch and functional data files.
OpenVSP2WOPWOP supports any number of DegenGeom files and loading conditions, which are specified in the comma-delimitted lists below (dataFileNames and CT).
A folder that corresponds to each DegenGeom file will be automatically generated in the user specified directory.
This folder will contain the geometry and compact geometry files as well as all the loading functional data files,
which correspond to the number of set operational conditions (CT).
'''''

#%%
import os
dirPyScript = os.getcwd()

#%%
#   Directory where the DegenGeom and airfoil polars are located, by default this is set to the location of the TestCase folder.
# dirDataFile= os.path.abspath(os.path.join(dirPyScript,'TestCase/BoeingModel360'))
dirDataFile ='/Users/danielweitsman/Desktop/Masters_Research/ACG/AS365 FF'
#%%
#   Names of DegenGeom files add as many geometry cases as you wish separated by commas.
dataFileNames = ["AS365N_DegenGeom.csv"]
#   Operational mode set to one for single geometry under multiple loading conditions, set equal to two for multiple geometries at single loading condition
OperMode = 2
#   Set equal to one in order to save the main dictionary in a .pkl file, which contains all the computed geometric and loading parameters for the case. This file can be read using the following command: pickle.load(open("file.pkl", "rb"))
savePickle = 1

#%% Airfoil polars setup
#   Name of XFoil files containing airfoil polars. The number of files should correspond to the number of airfoil sections.
airfoilPolarFileName = [["oa212Re24E5.dat","oa209Re24E5.dat"],["oa212Re24E5.dat","oa209Re24E5.dat"],["oa212Re30E5.dat","oa209Re30E5.dat"],
                         ["oa212Re36E5.dat","oa209Re36E5.dat"],["oa212Re39E5.dat","oa209Re39E5.dat"]]
#   Non-dimensional starting radial location of each airfoil cross section
XsecLocation = [0.29,0.88]
#   Starting angle of attack over which to evaluate the lift curve slop (degrees)
aStart = 2
#   Duration over which to evaluate the lift curve slope (degrees)
aLength = 3
#   Set to 1 in order for the airfoil polars to be plotted, along with the interval that the lift curve slope is evaluated over
check = 0

#%%     Input for loading module
# Number of Blades
Nb = 4
#   Forward velocity, set to zero for hover
Vx = [220*0.514444]
#   Vertical climb velocity
Vz = [0]
#   Shaft tilt angle (degrees), forward tilt is designated by a negative tilt angle
alphaShaft = [-12]
# Set desired operational thrust condition, (CT is the trim variable). Add as many comma-delimited loading conditions as you wish into the list. Separate functional data files will be written for each case.
T = [element * 4.44822*0.06 for element in [5400,6000,8000]]
#  Populate this list with the rotational rates (rpm) that were used to compute CT. The length of this comma-delimited list should be the same as that of CT.
omega = [298,308]
#   Initial collective pitch setting (degrees) used for trimming the rotor. This is an arbitrary value that should be adjusted to achieve convergence by ensuring that the computed angles of attack along the blade span lie within the limits of the Xfoil polars in the loading module.
thetaInit = 1
# Position of loading line from the blade's leading edge as a percentage of the chord
loadPos = 0.25
#   Set equal to 1 to include Pardtl's tip loss formulation
tipLoss = 1
# Density (kg/m^3)
rho = 1.225
# Speed of sound (m/s)
c = 340
#   Blade moment of inertia, only required for forward flight.
Ib = 1/3*217.423*(5.965-1.72985)**2
#   Nondimensionalized rotating flap frequency, this quantity is dependent on the non-rotating natural frequency and hinge offset of the blade, only required for forward flight.
nuBeta = 1.03

#%%     Namelist file write
#   Set equal to one in order to edit and write a namelist file for each case
nmlWrite = 1
#   Name of actual namelist file
NmlFileName = 'AS365N_FF.nam'
#   Duration of the case, expressed as shaft revolutions
nRev = 1
#   Set sample rate (Hz), preferably as a power of two.
nt = 2**14
#   Radial position of the observers from the rotor hub
radius = [304.8]
#   Number of in-plane observers
nbtheta = 60
#   Minimum angle of in-plane observer (degrees)
thetamin = 0
#   Maximum angle of in-plane observer (degrees)
thetamax = 360
#   Number of observers outside of the rotor plane (degrees)
nbpsi = 15
#   Minimum inclination angle of the out-of-plane observers (degrees)
psimin = -30
#   Maximum inclination angle of the out-of-plane observers (degrees)
psimax = 30


#%%     Patch file write parameters
# Patch file write directory, this is currently configured to create a directory titled 'PatchFiles' in dirDataFile, however it can be set to whatever path the user desires.
dirPatchFile = os.path.abspath(os.path.join(dirDataFile,'220knts'))
#   Lifting line compact patch file name
compactGeomFileName = "CompactGeom"
#   Blade geometry patch file name
geomFileName = "Geom"
#   Loading functional data file name
loadingFileName = "Load"

#%%     Pegg broadband noise set up

BBNoise = 0
#   Pegg Broadband data file name
bbFileName = "PeggBB"

#%% Packs user input parameters into a dictionary (no need to edit)
UserIn = {'dirDataFile':dirDataFile,'OperMode':OperMode,'savePickle':savePickle,'dataFileName':dataFileNames,'airfoilPolarFileName':airfoilPolarFileName,'XsecLocation':XsecLocation,
          'aStart':aStart,'aLength':aLength,'check':check,'Nb':Nb,'Vx': Vx,'Vz':Vz,'alphaShaft':alphaShaft,'omega':omega,'T':T,'thetaInit':thetaInit,'loadPos':loadPos,'tipLoss':tipLoss,'rho':rho,'c':c,'Ib':Ib,'nuBeta':nuBeta,
          'nmlWrite':nmlWrite,'dirPatchFile':dirPatchFile,'compactGeomFileName':compactGeomFileName,'geomFileName':geomFileName,'loadingFileName':loadingFileName,'bbFileName':bbFileName,
          'BBNoise':BBNoise,'NmlFileName':NmlFileName,'nRev':nRev,'nt':nt,'radius':radius,'nbtheta':nbtheta,'thetamin':thetamin,'thetamax':thetamax,'nbpsi':nbpsi,'psimin':psimin,'psimax':psimax}

#todo better estimate gamma, disable radial force component, run prediction without aircraft motion