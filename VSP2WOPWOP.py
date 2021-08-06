
'''

VSP2WOPWOP Main Script

Author: Daniel Weitsman

This is the main script which calls on and coordinates the execution of all the functions comprising the program.
There really shouldn't be a need to edit this script. If the code is ran in a python IDE three dictionaries,
titled geomParams,XsecPolar,loadParams will be returned these contain the analyzed geometric parameters of the
blades, the lift curve characteristics, and the aerodynamic loading/performance information, respectively.
'''

# %% Adds current command line path to the sys.path list so that the input.py module can be imported from the users
# current directory

import os
from sys import path

path.insert(0, os.getcwd())
from input import UserIn

from shutil import rmtree
from AnalyzeDegenGeom import AnalyzeDegenGeom
from ProcessGeom import ProcessGeom
from polarRead import polarRead
from loadingHover import loadingHover
from loadingFF import loadingFF
from ConstantLoadingPatchFileWrite import ConstantLoadingPatchFileWrite
from PeriodicLoadingPatchFileWrite import PeriodicLoadingPatchFileWrite
from nmlWrite import nml_write
from CaseFileWrite import caseFile_write
from ConstantBPMWrite import ConstantBPMWrite
from PeriodicBPMWrite import PeriodicBPMWrite
from GeomPatchFileWrite import GeomPatchFileWrite
from ErrorHandles import ErrorHandles
from writeHDF5 import writeHDF5

# %%

def main():

    def fwrite(UserIn, geomParams, XsecPolar, dirSaveFile,T, omega, Vx,Vz,alphaShaft,iter_geom):

        #   Writes out the blade geometry and lifting line compact geometry patch files
        GeomPatchFileWrite(UserIn['geomFileName'], geomParams, dirSaveFile)

        # This section of code determines whether to run the hover/axial or forward flight module and writes out
        # the corresponding constant or periodic functional data file, respectively.
        if Vx == 0:
            loadParams = loadingHover(UserIn, geomParams, XsecPolar, T, omega, Vz)
            ConstantLoadingPatchFileWrite(UserIn['loadingFileName'], loadParams, geomParams['nXsecs'], dirSaveFile)
        else:
            loadParams = loadingFF(UserIn, geomParams, XsecPolar, T, omega, Vx, Vz, alphaShaft)
            PeriodicLoadingPatchFileWrite(UserIn['loadingFileName'], loadParams, geomParams['nXsecs'], dirSaveFile)

        #   If the broadband noise flag is enabled the following modules are invoked to write out the BPM binary  data files.
        if UserIn['BBNoiseFlag'] == 1:
            if Vx == 0:
                ConstantBPMWrite(geomParams, loadParams, dirSaveFile)
            else:
                PeriodicBPMWrite(geomParams, loadParams, UserIn['nRev'], dirSaveFile)

        #   writes out namelist file
        if UserIn['nmlWrite'] == 1:
            nml_write(UserIn, loadParams, dirSaveFile, Vx, Vz, omega, alphaShaft, iter_geom, geomParams['nXsecs'])

        return loadParams

    #   Checks that the user inputs were provided correctly
    param_len = ErrorHandles(UserIn)

    #   Creates a parent directory for the case files to be written out to
    if not os.path.exists(os.path.join(os.getcwd(),UserIn['outputFolderName'])):
        os.mkdir(os.path.join(os.getcwd(),UserIn['outputFolderName']))

    #   Initializes empty dictionary and lists
    MainDict = {}
    loadParams = {}
    XsecPolar = {}
    globalFolder = []

    #   Iterates over each DegenGeom geometry file
    for iter_geom, dataFileName in enumerate(UserIn['dataFileName']):
        if isinstance(dataFileName, bytes):
            dataFileName=dataFileName.decode('utf-8')
        #   Parses and returns data contained in the DegenGeom file
        [dataSorted, indHeader] = AnalyzeDegenGeom(dataFileName)

        # Processes the DegenGeom to extract the blade geometric properties (e.g. radius, root cut-out, local solidity,
        # radial pitch and chord distributions)
        geomParams = ProcessGeom(dataSorted, indHeader, UserIn['loadPos'], UserIn['Nb'], UserIn['rotation'])

        # Reads and evaluates the XFoil polars, this is only done once during the first iteration of the outer for loop.
        if iter_geom == 0:
            for i, n in enumerate(UserIn['omega']):
                polarReadOut = polarRead(UserIn, i)
                XsecPolar = {**XsecPolar, **{f"{round(n)}RPM": polarReadOut}}

        # Creates a directory for each geometry where the respective loading, patch, and namelist files will be
        # written.
        dirSaveFile = os.path.abspath(os.path.join(os.getcwd(),UserIn['outputFolderName'], dataFileName[:-4]))
        if os.path.exists(dirSaveFile):
            rmtree(dirSaveFile)
        os.mkdir(dirSaveFile)

        #   Design Mode: Single loading condition/XFoil polar per DegenGeom geometry
        if UserIn['operMode'] == 1:

            loadParams = fwrite(UserIn, geomParams, XsecPolar[list(XsecPolar.keys())[iter_geom]], dirSaveFile, UserIn['T'][iter_geom], UserIn['omega'][iter_geom], UserIn['Vx'][iter_geom], UserIn['Vz'][iter_geom], UserIn['alphaShaft'][iter_geom],iter_geom)

            globalFolder.append(dataFileName[:-4])

            if iter_geom == len(UserIn['dataFileName']) - 1:
                caseFile_write(globalFolder, UserIn['NmlFileName'], os.path.join(os.getcwd(),UserIn['outputFolderName']))

        if UserIn['operMode'] == 2:

            # loops through the UserIn parameter with the greatest amount of values.
            for i in range(max(param_len)):
            # formats the path and name of the folder corresponding to each case
                globalFolderName = f"T_{round(UserIn['T'][i],1)}N_Nr_{round(UserIn['omega'][i],1)}RPM"
                print(globalFolderName)
                dirCaseFile = os.path.abspath(os.path.join(dirSaveFile, globalFolderName))

            # checks whether the case directory exists, and overwrites if true.
                if os.path.exists(dirCaseFile):
                    rmtree(dirCaseFile)
                os.mkdir(dirCaseFile)

            # computes blade loads and writes out the geometry, loading, BB noise, and namelist files.
                loadOut = fwrite(UserIn, geomParams, XsecPolar[list(XsecPolar.keys())[i]], dirCaseFile, UserIn['T'][i], UserIn['omega'][i], UserIn['Vx'][i], UserIn['Vz'][i], UserIn['alphaShaft'][i] , iter_geom)

            # compiles computed loads for each case into a dictionary
                loadParams = {**loadParams, **{globalFolderName: loadOut}}

            # adds the name of each case to the 'globalFolder' list
                globalFolder.append(globalFolderName)

            # writes out the cases.nam file
            caseFile_write([os.path.join(os.path.basename(dirSaveFile),x) for x in globalFolder], UserIn['NmlFileName'], os.path.dirname(dirSaveFile))


        #   Analysis Mode: Multiple loading condition per geometry
        if UserIn['operMode'] == 3:

            for itr_T, T in enumerate(UserIn['T']):
                for itr_Vx, Vx in enumerate(UserIn['Vx']):
                    for itr_Vz, Vz in enumerate(UserIn['Vz']):
                        for itr_omega, omega in enumerate(UserIn['omega']):

                            # formats the path and name of the folder corresponding to each case
                            globalFolderName = f"T_{round(T,1)}N_Vx_{round(Vx*1.944,1)}Kts_Vz_{round(Vz,1)}ms_Nr_{round(omega,1)}RPM"
                            print(globalFolderName)
                            dirCaseFile = os.path.abspath(os.path.join(dirSaveFile , globalFolderName))

                            # checks whether the case directory exists, and overwrites if true.
                            if os.path.exists(dirCaseFile):
                                rmtree(dirCaseFile)
                            os.mkdir(dirCaseFile)

                            # computes blade loads and writes out the geometry, loading, BB noise, and namelist files.
                            loadOut = fwrite(UserIn, geomParams, XsecPolar[list(XsecPolar.keys())[itr_omega]], dirSaveFile, T, omega, Vx, Vz, UserIn['alphaShaft'][0], iter_geom)

                            loadParams = {**loadParams, **{globalFolderName: loadOut}}

                            globalFolder.append(globalFolderName)

            caseFile_write([os.path.join(os.path.basename(dirSaveFile),x) for x in globalFolder], UserIn['NmlFileName'], dirSaveFile)

        MainDict = {**MainDict, **{'UserIn': UserIn,
                                   dataFileName[:-4]: {'geomParams': geomParams, 'XsecPolar': XsecPolar,
                                                       'loadParams': loadParams}}}

        if UserIn['saveHDF5'] == 1:
            writeHDF5(MainDict, os.path.join(os.getcwd(),UserIn['outputFolderName']))

    return MainDict


if __name__ == '__main__':
    MainDict = main()
