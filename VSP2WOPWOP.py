'''

VSP2WOPWOP Main Script

Author: Daniel Weitsman

This is the main script which calls on and coordinates the executes of all the functions comprising the program.
There really shouldn't be a need to edit this script. If the code is ran in a python IDE three dictionaries,
titled geomParams,XsecPolar,loadParams will be returned these contain the analyzed geometric parameters of the
blades, the lift curve characteristics, and the aerodynamic loading/performance information, respectively.
'''
# %%

from input import UserIn
import os
from shutil import rmtree
from functions.DegenGeom import ParseDegenGeom
from functions.GeomProcess import geomProcess
from functions.polarRead import polarRead
from functions.loadingHover import loadingAxialHover
from functions.loadingFF import loadingFF
from functions.GeomPatchFileWrite import GeomPatchFileWrite
from functions.ConstantLoadingPatchFileWrite import ConstantLoadingPatchFileWrite
from functions.PeriodicLoadingPatchFileWrite import PeriodicLoadingPatchFileWrite
from functions.CompactGeomPatchFileWrite import CompactGeomPatchFileWrite
from functions.nmlWrite import nml_write
from functions.CaseFileWrite import caseFile_write
from functions.PeggWrite import PeggBBDataFileWrite
from functions.ConstantBPMWrite import ConstantBPMWrite
from functions.FullGeomPatchFileWrite import FullGeomPatchFileWrite

#todo figure out lifting line normal orientation and corresponding cb, run BPM case with single geom patch file, update nml_write module
#todo write error handaling module for input file, attempt to replace the conditional statments in main

# %%
def main():
    #   Creates a parent directory for the case files to be written out to
    if os.path.exists(UserIn['dirPatchFile']) == 0:
        os.mkdir(UserIn['dirPatchFile'])

    #   Initializes empty dictionary and lists
    MainDict = {}
    loadParams = {}
    XsecPolar = {}
    globalFolder = []

    #   Iterates over each DegenGeom geometry file
    for iter_geom, dataFileName in enumerate(UserIn['dataFileName']):

        #   Parses and returns data contained in the DegenGeom file
        [dataSorted, indHeader] = ParseDegenGeom(UserIn['dirDataFile'], dataFileName)

        # Processes the DegenGeom to extract the blade geometric properties (e.g. radius, root cut-out, local solidity,
        # radial pitch and chord distributions)
        geomParams = geomProcess(dataSorted, indHeader, UserIn['loadPos'], UserIn['Nb'])

        # Reads and evaluates the XFoil polars, this is only done once during the first iteration of the outer for
        # loop.
        if iter_geom == 0:
            for i, n in enumerate(UserIn['omega']):
                polarReadOut = polarRead(UserIn, i)
                XsecPolar = {**XsecPolar, **{str(round(n)) + 'RPM': polarReadOut}}

        # Creates a directory for each geometry where the respective loading, patch, and namelist files will be
        # written.
        dirSaveFile = os.path.abspath(os.path.expanduser(UserIn['dirPatchFile'] + os.path.sep + dataFileName[:-4]))
        if os.path.exists(dirSaveFile) == 1:
            rmtree(dirSaveFile)
        os.mkdir(dirSaveFile)

        #   Design Mode: Single loading condition/XFoil polar per DegenGeom geometry
        if UserIn['OperMode'] == 1:

            #   Writes out the blade geometry and lifting line compact geometry patch files
            FullGeomPatchFileWrite(UserIn['geomFileName'], geomParams, dirSaveFile)
            CompactGeomPatchFileWrite(UserIn['compactGeomFileName'], geomParams['nXsecs'], geomParams['liftLineCoord'],geomParams['liftLineNorm'],
                                      dirSaveFile)

            #todo change conditional statments to assertions in perhaps a seperate module that handles input module errors

            # In the design mode each DegenGeom variant can be trimmed to the same or have its own respective thrust
            # condition. This if statement selects the correct thrust condition before passing it on to the loading
            # module.
            if len(UserIn['T']) > 1:
                T = UserIn['T'][iter_geom]
            else:
                T = UserIn['T'][0]

            if len(UserIn['Vz']) > 1:
                Vz = UserIn['Vz'][iter_geom]
            else:
                Vz = UserIn['Vz'][0]

            if len(UserIn['omega']) > 1:
                omega = UserIn['omega'][iter_geom]
            else:
                omega = UserIn['omega'][0]

            if len(UserIn['Vx']) > 1:
                Vx = UserIn['Vx'][iter_geom]
            else:
                Vx = UserIn['Vx'][0]

            if len(UserIn['alphaShaft']) > 1:
                alphaShaft = UserIn['alphaShaft'][iter_geom]
            else:
                alphaShaft = UserIn['alphaShaft'][0]

            # This section of code determines whether to run the hover/axial or forward flight module and writes out
            # the corresponding constant or periodic functional data file, respectively.
            if Vx == 0:
                loadParams = loadingAxialHover(UserIn, geomParams, XsecPolar[list(XsecPolar.keys())[iter_geom]], T, omega, Vz)
                ConstantLoadingPatchFileWrite(UserIn['loadingFileName'], loadParams, geomParams['nXsecs'], dirSaveFile)
            else:

                loadParams = loadingFF(UserIn, geomParams, XsecPolar[list(XsecPolar.keys())[iter_geom]], T, omega, Vx, Vz, alphaShaft)
                PeriodicLoadingPatchFileWrite(UserIn['loadingFileName'], loadParams, geomParams['nXsecs'], omega, dirSaveFile)

            if UserIn['BBNoiseFlag'] == 1:
                if UserIn['BBNoiseModel'] == 1:
                    PeggBBDataFileWrite(geomParams, loadParams,dirSaveFile)
                if UserIn['BBNoiseModel'] == 2:
                    ConstantBPMWrite(geomParams, loadParams,dirSaveFile)

            if UserIn['nmlWrite'] == 1:
                nml_write(UserIn, loadParams, dirSaveFile, Vx, Vz, omega, alphaShaft, iter_geom,geomParams['nXsecs'])

            globalFolder.append(dataFileName[:-4])

            if iter_geom == len(UserIn['dataFileName'])-1:
                caseFile_write(globalFolder, UserIn['NmlFileName'], UserIn['dirPatchFile'])

        #   Analysis Mode: Multiple loading condition per geometry
        if UserIn['OperMode'] == 2:

            for iter_thrust, nThrust in enumerate(UserIn['T']):
                for iter_Vx, nVx in enumerate(UserIn['Vx']):
                    for iter_Vz, nVz in enumerate(UserIn['Vz']):
                        for iter_omega, nOmega in enumerate(UserIn['omega']):

                            globalFolderName = 'T_'+'{:.2e}'.format(nThrust) + 'N_Vx_' + str(round(nVx * 1.944))\
                                               +'Kts_Vz_' + str(round(nVz)) + 'ms_Nr_'+ str(round(nOmega)) + 'RPM'
                            dirCaseFile = os.path.abspath(os.path.expanduser(dirSaveFile + os.path.sep + globalFolderName))

                            if os.path.exists(dirCaseFile) == 1:
                                rmtree(dirCaseFile)
                            os.mkdir(dirCaseFile)

                            FullGeomPatchFileWrite(UserIn['geomFileName'], geomParams, dirCaseFile)
                            CompactGeomPatchFileWrite(UserIn['compactGeomFileName'], geomParams['nXsecs'],
                                                      geomParams['liftLineCoord'],geomParams['liftLineNorm'], dirCaseFile)

                            if len(UserIn['alphaShaft']) > 1:
                                alphaShaft = UserIn['alphaShaft'][iter_Vx]
                            else:
                                alphaShaft = UserIn['alphaShaft'][0]

                            if nVx == 0:
                                loadingOut = loadingAxialHover(UserIn, geomParams, XsecPolar[list(XsecPolar.keys())[iter_omega]],
                                                               nThrust, nOmega,nVz)
                                ConstantLoadingPatchFileWrite(UserIn['loadingFileName'], loadingOut, geomParams['nXsecs'],
                                                              dirCaseFile)
                            else:
                                loadingOut = loadingFF(UserIn, geomParams, XsecPolar[list(XsecPolar.keys())[iter_omega]], nThrust,
                                                       nOmega, nVx, nVz, alphaShaft)
                                PeriodicLoadingPatchFileWrite(UserIn['loadingFileName'], loadingOut, geomParams['nXsecs'],
                                                              nOmega, dirCaseFile)

                            if UserIn['BBNoiseFlag'] == 1:
                                if UserIn['BBNoiseModel'] ==1:
                                    PeggBBDataFileWrite(geomParams, loadingOut,dirCaseFile)
                                if UserIn['BBNoiseModel'] == 2:
                                    ConstantBPMWrite(geomParams, loadingOut,dirCaseFile)

                            if UserIn['nmlWrite'] == 1:
                                nml_write(UserIn, loadingOut, dirCaseFile,nVx, nVz, nOmega, alphaShaft, iter_geom,geomParams['nXsecs'])

                            loadParams = {**loadParams, **{globalFolderName: loadingOut}}

                            globalFolder.append(globalFolderName)

            caseFile_write(globalFolder, UserIn['NmlFileName'], dirSaveFile)

        MainDict = {**MainDict, **{'UserIn': UserIn,
                                   dataFileName[:-4]: {'geomParams': geomParams, 'XsecPolar': XsecPolar,
                                                       'loadParams': loadParams}}}

        if UserIn['savePickle'] == 1:
            import pickle
            with open(os.path.abspath(os.path.expanduser(UserIn['dirPatchFile'] + os.path.sep + 'MainDict.pkl')),
                      "wb") as f:
                pickle.dump(MainDict, f)

        # Use the following command to load data
        # pickle.load(open("file.pkl", "rb"))

    return MainDict


if __name__ == '__main__':
    print(__name__)
    MainDict = main()
