'''

VSP2WOPWOP Main Script

Author: Daniel Weitsman

This is the main script which calls on and coordinates the executes of all the functions comprising the program.
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
from designModeVal import designModeVal
from writeHDF5 import writeHDF5
import argparse

import h5py
# %%
def main():
    #   Checks that the user inputs were provided correctly
    ErrorHandles(UserIn)

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
        [dataSorted, indHeader] = AnalyzeDegenGeom(UserIn['dirDataFile'], dataFileName)

        # Processes the DegenGeom to extract the blade geometric properties (e.g. radius, root cut-out, local solidity,
        # radial pitch and chord distributions)
        geomParams = ProcessGeom(dataSorted, indHeader, UserIn['loadPos'], UserIn['Nb'], UserIn['rotation'])

        # Reads and evaluates the XFoil polars, this is only done once during the first iteration of the outer for
        # loop.
        if iter_geom == 0:
            for i, n in enumerate(UserIn['omega']):
                polarReadOut = polarRead(UserIn, i)
                XsecPolar = {**XsecPolar, **{str(round(n)) + 'RPM': polarReadOut}}

        # Creates a directory for each geometry where the respective loading, patch, and namelist files will be
        # written.
        dirSaveFile = os.path.abspath(os.path.expanduser(os.getcwd() + os.path.sep + dataFileName[:-4]))
        if os.path.exists(dirSaveFile) == 1:
            rmtree(dirSaveFile)
        os.mkdir(dirSaveFile)

        #   Design Mode: Single loading condition/XFoil polar per DegenGeom geometry
        if UserIn['OperMode'] == 1:

            #   Writes out the blade geometry and lifting line compact geometry patch files
            GeomPatchFileWrite(UserIn['geomFileName'], geomParams, dirSaveFile)

            # This function returns the values, which are specified as lists in the input module, corresponding to the
            # current DegenGeom index
            T, Vz, Vx, omega, alphaShaft, XsecPolar_select = designModeVal(UserIn, XsecPolar, iter_geom)

            # This section of code determines whether to run the hover/axial or forward flight module and writes out
            # the corresponding constant or periodic functional data file, respectively.
            if Vx == 0:

                if isinstance(args.mat_loading_file,str):
                    import numpy as np
                    mat_read = h5py.File(os.path.join(os.getcwd(), args.mat_loading_file), 'r')

                    if args.load_type:
                        loadParams = {'dFx':mat_read[args.rotor]['inplane'][args.azimuth,:],'dFy':np.zeros(np.shape(mat_read[args.rotor]['inplane'][args.azimuth,:])),'dFz':mat_read[args.rotor]['outplane'][args.azimuth,:],'th':[0,0,0]}
                        ConstantLoadingPatchFileWrite(UserIn['loadingFileName'], loadParams, geomParams['nXsecs'], dirSaveFile)
                    else:
                        loadParams = {'dFx':mat_read[args.rotor]['inplane'][:],'dFy':np.zeros(np.shape(mat_read[args.rotor]['inplane'][:])),'dFz':mat_read[args.rotor]['outplane'][:],'th':[0,0,0],'phi': np.linspace(0,2*np.pi,np.shape(mat_read[args.rotor]['inplane'][:])[0])}
                        PeriodicLoadingPatchFileWrite(UserIn['loadingFileName'], loadParams, geomParams['nXsecs'], omega,dirSaveFile)

                else:
                    loadParams = loadingHover(UserIn, geomParams, XsecPolar_select, T, omega, Vz)
                    ConstantLoadingPatchFileWrite(UserIn['loadingFileName'], loadParams, geomParams['nXsecs'], dirSaveFile)
            else:
                loadParams = loadingFF(UserIn, geomParams, XsecPolar_select, T, omega, Vx, Vz, alphaShaft)
                PeriodicLoadingPatchFileWrite(UserIn['loadingFileName'], loadParams, geomParams['nXsecs'], omega,
                                              dirSaveFile)

            if UserIn['BBNoiseFlag'] == 1:
                if Vx == 0:
                    ConstantBPMWrite(geomParams, loadParams, dirSaveFile)
                else:
                    PeriodicBPMWrite(geomParams, loadParams, UserIn['nRev'], omega, dirSaveFile)

            if UserIn['nmlWrite'] == 1:
                nml_write(UserIn, loadParams, dirSaveFile, Vx, Vz, omega, alphaShaft, iter_geom, geomParams['nXsecs'])

            globalFolder.append(dataFileName[:-4])

            if iter_geom == len(UserIn['dataFileName']) - 1:
                caseFile_write(globalFolder, UserIn['NmlFileName'], UserIn['dirPatchFile'])

        #   Analysis Mode: Multiple loading condition per geometry
        if UserIn['OperMode'] == 2:

            for iter_thrust, nThrust in enumerate(UserIn['T']):
                for iter_Vx, nVx in enumerate(UserIn['Vx']):
                    for iter_Vz, nVz in enumerate(UserIn['Vz']):
                        for iter_omega, nOmega in enumerate(UserIn['omega']):

                            globalFolderName = 'T_' + '{:.2e}'.format(nThrust) + 'N_Vx_' + str(round(nVx * 1.944)) \
                                               + 'Kts_Vz_' + str(round(nVz)) + 'ms_Nr_' + str(round(nOmega)) + 'RPM'
                            dirCaseFile = os.path.abspath(
                                os.path.expanduser(dirSaveFile + os.path.sep + globalFolderName))

                            if os.path.exists(dirCaseFile) == 1:
                                rmtree(dirCaseFile)
                            os.mkdir(dirCaseFile)

                            GeomPatchFileWrite(UserIn['geomFileName'], geomParams, dirCaseFile)

                            if len(UserIn['alphaShaft']) > 1:
                                alphaShaft = UserIn['alphaShaft'][iter_Vx]
                            else:
                                alphaShaft = UserIn['alphaShaft'][0]

                            if nVx == 0:
                                loadingOut = loadingHover(UserIn, geomParams,
                                                          XsecPolar[list(XsecPolar.keys())[iter_omega]],
                                                          nThrust, nOmega, nVz)
                                ConstantLoadingPatchFileWrite(UserIn['loadingFileName'], loadingOut,
                                                              geomParams['nXsecs'],
                                                              dirCaseFile)
                            else:
                                loadingOut = loadingFF(UserIn, geomParams,
                                                       XsecPolar[list(XsecPolar.keys())[iter_omega]], nThrust, nOmega,
                                                       nVx, nVz, alphaShaft)
                                PeriodicLoadingPatchFileWrite(UserIn['loadingFileName'], loadingOut,
                                                              geomParams['nXsecs'],
                                                              nOmega, dirCaseFile)

                            if UserIn['BBNoiseFlag'] == 1:
                                if nVx == 0:
                                    ConstantBPMWrite(geomParams, loadingOut, dirSaveFile)
                                else:
                                    PeriodicBPMWrite(geomParams, loadingOut, UserIn['nRev'], nOmega, dirSaveFile)

                            if UserIn['nmlWrite'] == 1:
                                nml_write(UserIn, loadingOut, dirCaseFile, nVx, nVz, nOmega, alphaShaft, iter_geom,
                                          geomParams['nXsecs'])

                            loadParams = {**loadParams, **{globalFolderName: loadingOut}}

                            globalFolder.append(globalFolderName)

            caseFile_write(globalFolder, UserIn['NmlFileName'], dirSaveFile)

        MainDict = {**MainDict, **{'UserIn': UserIn,
                                   dataFileName[:-4]: {'geomParams': geomParams, 'XsecPolar': XsecPolar,
                                                       'loadParams': loadParams}}}

        if UserIn['saveHDF5'] == 1:
            writeHDF5(MainDict, UserIn)

    return MainDict


parser = argparse.ArgumentParser()
parser.add_argument('mat_loading_file', help="Name of mat file containing blade loading information.",type = str)
parser.add_argument('-op','--rotor',help = "Rotor from which to reference the loads (op1 or op2), defaults to op1.",type=str,default = 'op1')

group = parser.add_mutually_exclusive_group()
group.add_argument('-p',"--load_type",help = "Specifies if the loads are constant or periodic. Include -p as a command line argument if loads are periodic and vary around the rotor disk but repeat every rotor revolution, exclude otherwise.",action='store_false')
group.add_argument('-phi','--azimuth',help = 'If the loads are specified as constant, select a corresponding azimuthal index from which to reference the loads.',type = int)

# parser.add_argument("-l","--load type",help = "Are the loads constant or periodic",action=store_true)

args = parser.parse_args()


if __name__ == '__main__':
    MainDict = main()
