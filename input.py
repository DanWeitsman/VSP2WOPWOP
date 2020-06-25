# %%        VSP2WOPWOP User Input File

#      Author: Daniel Weitsman
"""'' This input file is used to set up the cases for which to generate the patch and functional data files.
OpenVSP2WOPWOP supports any number of DegenGeom files and loading conditions, which are specified in the
comma-delimitted lists below (dataFileNames and CT). A folder that corresponds to each DegenGeom file will be
automatically generated in the user specified directory. This folder will contain the geometry and compact geometry
files as well as all the loading functional data files, which correspond to the number of set operational conditions
(CT). '' """

# %% Package import
import os

dirPyScript = os.getcwd()

# Directory where the DegenGeom and airfoil polars are located, by default this is set to the location of the
# test case folder.
dirDataFile = '/Users/danielweitsman/Desktop/Masters_Research/ACG/AS365 FF'

# %%
#   Names of DegenGeom files add as many geometry cases as you wish separated by commas.
dataFileNames = ["AS365N_DegenGeom.csv"]
# Operational mode: set equal to one for design mode, which is ideal for geometric parametric studies, where each
# variant of the blade geometry geometry corresponds to a different operating condition. Set this quantity equal to
# two to run in analysis mode. The analysis mode is tailored for running operational condition sweeps on each blade
# geometry.
OperMode = 2
# Set equal to one in order to save the main dictionary, which contains all the computed geometric and loading
# parameters for the case to a .pkl file. This file can be read using the following command: pickle.load(open(
# "file.pkl", "rb"))
savePickle = 1

# %%
''''Patch file write directory setup'''

# Directory were you would like to write out the patch and functional data flies. By default this path is configured
# to be a seperate folder in the 'dirDataFile'.
dirPatchFile = os.path.abspath(os.path.join(dirDataFile, '220knts'))
#   Lifting line compact patch file name (without the extension).
compactGeomFileName = "CompactGeom"
#   Blade geometry patch file name (without the extension).
geomFileName = "Geom"
#   Loading functional data file name (without the extension).
loadingFileName = "Load"

# %%
'''Airfoil Cross Section Configuration'''

# Name of XFoil files containing the airfoil cross section polars. The number of files should correspond to the
# number of 'XsecLocation'.
airfoilPolarFileName = [["oa212Re24E5.dat", "oa209Re24E5.dat"], ["oa212Re24E5.dat", "oa209Re24E5.dat"],
                        ["oa212Re30E5.dat", "oa209Re30E5.dat"],
                        ["oa212Re36E5.dat", "oa209Re36E5.dat"], ["oa212Re39E5.dat", "oa209Re39E5.dat"]]
#   Non-dimensional radial location of each airfoil cross section
XsecLocation = [0.29, 0.88]
# Starting angle of attack over which to evaluate the lift curve slop (degrees). The airfoil properties must have
# been evaluated at this angle of attack.
aStart = 2
#   Range over which to evaluate to evaluate the lift curve slope (degrees).
aLength = 3
# Set equal to one for the airfoil polars to be plotted, along with the interval that the lift curve slope is evaluated.
check = 0

# %%
''' Operating conditions configuration '''
# Number of Blades
Nb = 4
#   Forward velocity, set to zero for hover, input as many comma-delimited forward flight velocities as you wish to be
# evaluated.
Vx = [220 * 0.514444]
#   Vertical climb velocity
Vz = [0]
#   Shaft tilt angle (degrees), forward tilt is designated by a negative tilt angle
alphaShaft = [-12]
#   Set the gross weight of the aircraft (N). Add as many comma-delimited thrust conditions as you wish into the list.
# Separate functional data files will be written for each case. If you are running in the design mode ('OperMode' = 1),
# the length of this list should be equal to the number of geometric cases ('dataFileNames').
T = [element * 4.44822 * 0.06 for element in [5400, 6000, 8000]]
# Populate this list with the rotational rates (rpm) of the rotor. If you are running in the design mode ('OperMode'
# = 1), the length of this list should be equal to the number of geometric cases ('dataFileNames').
omega = [298, 308]
#   Initial collective pitch setting (degrees) used for trimming the rotor. This is an arbitrary value that should be
# adjusted to achieve convergence by ensuring that the computed angles of attack along the blade span lie within the
# limits of the Xfoil polars.
thetaInit = 1
#   Position of loading line from the blade's leading edge as a percentage of the chord.
loadPos = 0.25
#   Set equal to 1 to include Pardtl's tip loss formulation
tipLoss = 1
# Density (kg/m^3)
rho = 1.225
# Speed of sound (m/s)
c = 340
#   Blade moment of inertia (only required for forward flight).
Ib = 1 / 3 * 217.423 * (5.965 - 1.72985) ** 2
#   Nondimensionalized rotating flap frequency, this quantity is dependent on the non-rotating natural frequency and
# hinge offset of the blade (only required for forward flight).
nuBeta = 1.03

#%%
'''Namelist file configuration'''
#   Set equal to one in order to write out a namelist file for each case.
nmlWrite = 1
#   Name of the namelist file that will be the same for all cases.
NmlFileName = 'AS365N_FF.nam'
#   Duration of the case, expressed as shaft revolutions
nRev = 1
#   Set sample rate (Hz), preferably as a power of two.
nt = 2**14
# Observer type set equal to one for a single observer, two for a rectangular observer grid, and three for a
# spherical grid. Then complete the respective section below.
obsType = 3

# %%
'''Single Observer'''
#   x-position (ahead of the rotor) of observer
xLoc = 0
#   y-position (to the side of the rotor) of observer
yLoc = 0
#   z-position (vertical) of observer
zLoc = 0

# %%
'''Rectangular Observer Grid'''
#   Number of observers along the x-direction (ahead of the rotor).
nbx = 0
#   Minimum x coordinate
xMin = 0
#   Maximum x coordinate
xMax = 0
#   Number of observers along the y-direction (to the side of the rotor).
nby = 0
#   Minimum y coordinate
yMin = 0
#   Maximum y coordinate
yMax = 0
#   Number of observers along the z-direction (normal to the rotor plane).
nbz = 0
#   Minimum z coordinate
zMin = 0
#   Maximum z coordinate
zMax = 0

# %%
'''Spherical Observer Grid'''
#   Radial position of the observers from the rotor hub.
radius = [304.8]
#   Number of in-plane observers
nbtheta = 60
#   Minimum angle of in-plane observer (degrees)
thetamin = 0
#   Maximum angle of in-plane observer (degrees)
thetamax = 360
#   Number of observers outside of the rotor plane (degrees)
nbpsi = 11
#   Minimum inclination angle of the out-of-plane observers (degrees)
psimin = -30
#   Maximum inclination angle of the out-of-plane observers (degrees)
psimax = 30

# %%
'''Pegg broadband noise set up '''

BBNoise = 0
#   Pegg Broadband data file name
bbFileName = "PeggBB"

# %% Packs user input parameters into a dictionary (no need to edit)
UserIn = {'dirDataFile': dirDataFile, 'OperMode': OperMode, 'savePickle': savePickle, 'dataFileName': dataFileNames,
          'airfoilPolarFileName': airfoilPolarFileName, 'XsecLocation': XsecLocation,
          'aStart': aStart, 'aLength': aLength, 'check': check, 'Nb': Nb, 'Vx': Vx, 'Vz': Vz, 'alphaShaft': alphaShaft,
          'omega': omega, 'T': T, 'thetaInit': thetaInit, 'loadPos': loadPos, 'tipLoss': tipLoss, 'rho': rho, 'c': c,
          'Ib': Ib, 'nuBeta': nuBeta,'nmlWrite': nmlWrite, 'dirPatchFile': dirPatchFile, 'compactGeomFileName':compactGeomFileName,
          'geomFileName': geomFileName, 'loadingFileName': loadingFileName, 'bbFileName': bbFileName,
          'BBNoise': BBNoise, 'NmlFileName': NmlFileName, 'nRev': nRev, 'nt': nt, 'obsType': obsType, 'xLoc': xLoc,
          'yLoc': yLoc, 'zLoc': zLoc, 'nbx': nbx, 'nby': nby, 'nbz': nbz, 'xMin': xMin, 'yMin': yMin, 'zMin': zMin,
          'xMax': xMax, 'yMax': yMax, 'zMax': zMax,'radius': radius, 'nbtheta': nbtheta,
          'thetamin': thetamin, 'thetamax': thetamax, 'nbpsi': nbpsi, 'psimin': psimin, 'psimax': psimax}
