'''
VSP2WOPWOP User Input File

Author: Daniel Weitsman

This input file is used to configure the cases for which to generate the patch and functional data files. VSP2WOPWOP
supports any number of DegenGeom files operating at any number of loading conditions, which are specified in the
comma-delimitted lists below (dataFileNames, Vx, W, omega). A folder that corresponds to each DegenGeom file,
containing the patch and namelist files for each loading condition, will be generated in the user specified directory.
'''

# %%
# Names of DegenGeom files add as many geometry cases as you wish separated by commas.
dataFileNames = ["agilitygofly_DegenGeom.csv"]

#%%
'''Airfoil Cross Section Configuration'''

# Name of XFoil files containing the airfoil cross section polars. The number of files should correspond to the
# number of 'XsecLocation'.
airfoilPolarFileName = [['test_polar_v1.dat']]

#   Non-dimensional radial location of each airfoil cross section
XsecLocation = [0.268]

# Starting angle of attack over which to evaluate the lift curve slop (degrees). The airfoil properties must have
# been evaluated at this angle of attack.
aStart = 2

#   Range over which to evaluate to evaluate the lift curve slope (degrees).
aLength = 3

# Set equal to one for the airfoil polars to be plotted, along with the interval that the lift curve slope is evaluated.
check = 0

# %%
'''Vehicle Configuration'''
#   Number of rotors
Nrotor = 2
#   location of each rotor with respect to the origin
rotorLoc = [[0,0,0.3048],[0,0,-0.3048]]
#   Direction of rotation of each rotor set equal to 1 for CCW and -1 for CW, when viewed from above.
rotation = [1,-1]

#%%

# Number of Blades
Nb = 4

#   Forward velocity, set to zero for hover, input as many comma-delimited forward flight velocities as you wish to be
# evaluated.
Vx = 0

#   Vertical climb velocity (m/s), set to positive for ascent and negative for descent
Vz = 0

#   Shaft tilt angle (degrees), forward tilt is designated by a negative tilt angle
alphaShaft = 0

#   Position of loading line from the blade's leading edge as a percentage of the chord.
loadPos = 0.25

# Density (kg/m^3)
rho = 1.225

# Speed of sound (m/s)
c = 340

# %%
'''Broadband noise analysis configuration '''

# Set equal to '1' to conduct a broadband noise prediction,'0' otherwise. The Brooks, Pope, and Marcolini (BPM)
# semi-empirical method is used for the prediction. Does not work yet with coaxial op data - 3/18/21
BBNoiseFlag = 0

# %%
'''Namelist file configuration'''
#   Set equal to one in order to write out a namelist file for each case.
nmlWrite = 1

#   Name of the namelist file that will be the same for all cases.
NmlFileName = 'op.nam'

#%%
'''Observer namelist configuration'''
#   Duration of the case, expressed as shaft revolutions
nRev = 1

#   Set sample rate (Hz), preferably as a power of two.
nt = 2 ** 14

# Observer type set equal to '1' for a single observer, '2' for a rectangular observer grid, and '3' for a
# spherical grid. Then complete the respective section below.
obsType = 2

# %%
'''Single Observer'''
#   x-position (ahead of the rotor) of observer
xLoc = 0

#   y-position (to the side of the rotor) of observer
yLoc = 14.165447920000007

#   z-position (vertical) of observer
zLoc = 14.165447920000007

# %%
'''Rectangular Observer Grid'''
#   Number of observers along the x-direction (ahead of the rotor).
nbx = 1

#   Minimum x coordinate
xMin = 0

#   Maximum x coordinate
xMax = [0]

#   Number of observers along the y-direction (to the side of the rotor).
nby = 1

#   Minimum y coordinate
yMin = 14.165447920000007

#   Maximum y coordinate
yMax = [14.165447920000007]

#   Number of observers along the z-direction (normal to the rotor plane).
nbz = 3

#   Minimum z coordinate
zMin = -1.4888485708273682

#   Maximum z coordinate
zMax = [1.4888485708273682]

# %%
'''Spherical Observer Grid'''
#   Radial position of the observers from the rotor hub.
radius = [3.2]

#   Number of in-plane observers
nbTheta = 1

#   Minimum angle of in-plane observer (degrees)
thetaMin = 0

#   Maximum angle of in-plane observer (degrees)
thetaMax = 0

#   Number of observers outside of the rotor plane (degrees)
nbPsi = 4

#   Minimum inclination angle of the out-of-plane observers (degrees)
psiMin = -45

#   Maximum inclination angle of the out-of-plane observers (degrees)
psiMax = 0

# %% Packs user input parameters into a dictionary (no need to edit)
UserIn = {
          'dataFileName':dataFileNames,'airfoilPolarFileName': airfoilPolarFileName, 'XsecLocation': XsecLocation,
          'aStart': aStart, 'aLength': aLength, 'check': check,'Nrotor':Nrotor,'rotorLoc':rotorLoc, 'Nb': Nb,'rotation':rotation, 'Vx': Vx, 'Vz': Vz, 'alphaShaft': alphaShaft,
          'loadPos': loadPos,'rho': rho, 'c': c,
          'nmlWrite': nmlWrite,
          'BBNoiseFlag': BBNoiseFlag, 'NmlFileName': NmlFileName, 'nRev': nRev, 'nt': nt, 'obsType': obsType, 'xLoc': xLoc,
          'yLoc': yLoc, 'zLoc': zLoc, 'nbx': nbx, 'nby': nby, 'nbz': nbz, 'xMin': xMin, 'yMin': yMin, 'zMin': zMin,
          'xMax': xMax, 'yMax': yMax, 'zMax': zMax, 'radius': radius, 'nbtheta': nbTheta,
          'thetamin': thetaMin, 'thetamax': thetaMax, 'nbpsi': nbPsi, 'psimin': psiMin, 'psimax': psiMax}
