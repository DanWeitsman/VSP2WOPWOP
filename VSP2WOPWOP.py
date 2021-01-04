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
import numpy as np
import os
from shutil import rmtree
from functions.DegenGeom import ParseDegenGeom
from functions.GeomProcess import geomProcess
from functions.polarRead import polarRead
from functions.loadingHover import loadingAxialHover
from functions.ConstantLoadingPatchFileWrite import ConstantLoadingPatchFileWrite
from functions.PeriodicLoadingPatchFileWrite import PeriodicLoadingPatchFileWrite
from functions.nmlWrite import nml_write
from functions.CaseFileWrite import caseFile_write
from functions.PeggWrite import PeggBBDataFileWrite
from functions.ConstantBPMWrite import ConstantBPMWrite
from functions.PeriodicBPMWrite import PeriodicBPMWrite
from functions.GeomPatchFileWrite import GeomPatchFileWrite
from functions.ErrorHandles import ErrorHandles
from functions.loadingFF import  loadingFF
from functions.designModeVal import designModeVal
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
        [dataSorted, indHeader] = ParseDegenGeom(UserIn['dirDataFile'], dataFileName)

        # Processes the DegenGeom to extract the blade geometric properties (e.g. radius, root cut-out, local solidity,
        # radial pitch and chord distributions)
        geomParams = geomProcess(dataSorted, indHeader, UserIn['loadPos'], UserIn['Nb'],UserIn['rotation'])

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
            GeomPatchFileWrite(UserIn['geomFileName'], geomParams, dirSaveFile)

            # This function returns the values, which are specified as lists in the input module, corresponding to the
            # current DegenGeom index
            T, Vz, Vx, omega, alphaShaft,XsecPolar_select = designModeVal(UserIn,XsecPolar, iter_geom)

            # This section of code determines whether to run the hover/axial or forward flight module and writes out
            # the corresponding constant or periodic functional data file, respectively.
            if Vx == 0:
                loadParams = loadingAxialHover(UserIn, geomParams, XsecPolar_select, T, omega, Vz)
                ConstantLoadingPatchFileWrite(UserIn['loadingFileName'], loadParams, geomParams['nXsecs'], dirSaveFile)
            else:
                loadParams = loadingFF(UserIn, geomParams, XsecPolar_select, T, omega, Vx, Vz, alphaShaft)
                PeriodicLoadingPatchFileWrite(UserIn['loadingFileName'], loadParams, geomParams['nXsecs'], omega, dirSaveFile)

            if UserIn['BBNoiseFlag'] == 1:
                if UserIn['BBNoiseModel'] == 1:
                    PeggBBDataFileWrite(geomParams, loadParams,dirSaveFile)
                if UserIn['BBNoiseModel'] == 2:
                    if Vx == 0:
                        ConstantBPMWrite(geomParams, loadParams,dirSaveFile)
                    else:
                        PeriodicBPMWrite(geomParams,loadParams,UserIn['nRev'],omega,dirSaveFile)

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

                            GeomPatchFileWrite(UserIn['geomFileName'], geomParams, dirCaseFile)

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
                                loadingOut = loadingFF(UserIn, geomParams, XsecPolar[list(XsecPolar.keys())[iter_omega]], nThrust, nOmega, nVx, nVz, alphaShaft)
                                PeriodicLoadingPatchFileWrite(UserIn['loadingFileName'], loadingOut, geomParams['nXsecs'],
                                                              nOmega, dirCaseFile)

                            if UserIn['BBNoiseFlag'] == 1:
                                if UserIn['BBNoiseModel'] == 1:
                                    PeggBBDataFileWrite(geomParams, loadingOut,dirCaseFile)
                                if UserIn['BBNoiseModel'] == 2:
                                    if nVx == 0:
                                        ConstantBPMWrite(geomParams, loadingOut, dirSaveFile)
                                    else:
                                        PeriodicBPMWrite(geomParams, loadingOut, UserIn['nRev'], nOmega, dirSaveFile)

                            if UserIn['nmlWrite'] == 1:
                                nml_write(UserIn, loadingOut, dirCaseFile,nVx, nVz, nOmega, alphaShaft, iter_geom,geomParams['nXsecs'])

                            loadParams = {**loadParams, **{globalFolderName: loadingOut}}

                            globalFolder.append(globalFolderName)

            caseFile_write(globalFolder, UserIn['NmlFileName'], dirSaveFile)

        MainDict = {**MainDict, **{'UserIn': UserIn,
                                   dataFileName[:-4]: {'geomParams': geomParams, 'XsecPolar': XsecPolar,
                                                       'loadParams': loadParams}}}

        if UserIn['saveHDF5'] == 1:
            import h5py
            with h5py.File(os.path.abspath(os.path.join(UserIn['dirPatchFile'], 'MainDict.h5')), 'w') as f_write:

                def encode_str_list(str_list):
                    for i, elem in enumerate(str_list):
                        if isinstance(elem, list):
                            encode_str_list(elem)
                        elif isinstance(elem, str):
                            str_list[i] = elem.encode()
                        else:
                            break
                    return str_list

                def write_dict_hdf5(d,parent=''):
                    for key , value in d.items():
                        if isinstance(value, dict):
                            write_dict_hdf5(value, parent+'/'+key)
                        else:
                            # print(key)
                            if isinstance(value, list):
                                value = encode_str_list(value)
                            elif isinstance(value, str):
                                value = value.encode()
                            f_write.create_dataset(parent+'/'+key, shape=np.shape(value), data=value)

                write_dict_hdf5(MainDict)

    return MainDict


if __name__ == '__main__':
    print(__name__)
    MainDict = main()

    # f = h5py.File(os.path.abspath(os.path.expanduser(UserIn['dirPatchFile'] + os.path.sep + 'MainDict.hdf5')),'w')