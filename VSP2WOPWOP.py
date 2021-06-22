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
from ConstantLoadingPatchFileWrite import ConstantLoadingPatchFileWrite
from PeriodicLoadingPatchFileWrite import PeriodicLoadingPatchFileWrite
from nmlWrite import nml_write
from CaseFileWrite import caseFile_write
from ConstantBPMWrite import ConstantBPMWrite
from PeriodicBPMWrite import PeriodicBPMWrite
from GeomPatchFileWrite import GeomPatchFileWrite
from ErrorHandles import ErrorHandles
import argparse

import h5py
# %%
def main():
    #   Checks that the user inputs were provided correctly
    ErrorHandles(UserIn)

    #   Initializes empty dictionary and lists
    MainDict = {}
    XsecPolar = {}
    globalFolder = []

    #   Iterates over each DegenGeom geometry file
    for iter_geom, dataFileName in enumerate(UserIn['dataFileName']):

        # Reads and evaluates the XFoil polars, this is only done once during the first iteration of the outer for
        # loop.
        if iter_geom == 0:
            for i, n in enumerate(UserIn['airfoilPolarFileName']):
                polarReadOut = polarRead(UserIn, n)
                XsecPolar = {**XsecPolar, **{i: polarReadOut}}

        # Creates a directory for each geometry where the respective loading, patch, and namelist files will be
        # written.
        dirSaveFile = os.path.abspath(os.path.expanduser(os.getcwd() + os.path.sep + dataFileName[:-4]))
        if os.path.exists(dirSaveFile) == 1:
            rmtree(dirSaveFile)
        os.mkdir(dirSaveFile)

        #   Parses and returns data contained in the DegenGeom file
        [dataSorted, indHeader] = AnalyzeDegenGeom(dataFileName)

        # This function returns the values, which are specified as lists in the input module, corresponding to the
        # current DegenGeom index
        # T, Vz, Vx, omega, alphaShaft, XsecPolar_select = designModeVal(UserIn, XsecPolar, iter_geom)

        # This section of code determines whether to run the hover/axial or forward flight module and writes out
        # the corresponding constant or periodic functional data file, respectively.

        if isinstance(args.mat_loading_file,str):
            import numpy as np
            mat_read = h5py.File(os.path.join(os.getcwd(), args.mat_loading_file), 'r')

        # Processes the DegenGeom to extract the blade geometric properties (e.g. radius, root cut-out, local solidity,
        # radial pitch and chord distributions)
        for i in range(UserIn['Nrotor']):
            geomParams = ProcessGeom(dataSorted, indHeader, UserIn['loadPos'], UserIn['Nb'], UserIn['rotation'][i])
            GeomPatchFileWrite('op'+ str(i+1) +'_geom', geomParams, dirSaveFile)

            if args.load_type:
                loadParams = {'dFx':mat_read['op'+ str(i+1)]['inplane'][args.azimuth,:],'dFy':np.zeros(np.shape(mat_read['op'+ str(i+1)]['inplane'][args.azimuth,:])),'dFz':mat_read['op'+ str(i+1)]['outplane'][args.azimuth,:],'th':[0,0,0]}
                ConstantLoadingPatchFileWrite('op'+ str(i+1) +'_load', loadParams, geomParams['nXsecs'], dirSaveFile)
            else:
                loadParams = {'dFx':mat_read['op'+ str(i+1)]['inplane'],'dFy':np.zeros(np.shape(mat_read['op'+ str(i+1)]['inplane'])),'dFz':mat_read['op'+ str(i+1)]['outplane'],'th':[0,0,0],'phi': np.linspace(0,2*np.pi,np.shape(mat_read['op'+ str(i+1)]['inplane'])[0])}
                PeriodicLoadingPatchFileWrite('op'+ str(i+1) +'_load', loadParams, geomParams['nXsecs'], np.squeeze(mat_read['op'+ str(i+1)]['rpm'][()]) ,dirSaveFile)

        # if UserIn['BBNoiseFlag'] == 1:
        #     if Vx == 0:
        #         ConstantBPMWrite(geomParams, loadParams, dirSaveFile)
        #     else:
        #         PeriodicBPMWrite(geomParams, loadParams, UserIn['nRev'], omega, dirSaveFile)

        if UserIn['nmlWrite'] == 1:
            nml_write(UserIn, dirSaveFile,np.squeeze(mat_read['op'+ str(i+1)]['rpm'][()]), iter_geom, geomParams['nXsecs'],mat_read)

        globalFolder.append(dataFileName[:-4])

        if iter_geom == len(UserIn['dataFileName']) - 1:
            caseFile_write(globalFolder, UserIn['NmlFileName'])

        caseFile_write(globalFolder, UserIn['NmlFileName'])


    return MainDict


parser = argparse.ArgumentParser()
parser.add_argument('mat_loading_file', help="Name of mat file containing blade loading information.",type = str)

group = parser.add_mutually_exclusive_group()
group.add_argument('-p',"--load_type",help = "Specifies if the loads are constant or periodic. Include -p as a command line argument if loads are periodic and vary around the rotor disk but repeat every rotor revolution, exclude otherwise.",action='store_false')
group.add_argument('-phi','--azimuth',help = 'If the loads are specified as constant, select a corresponding azimuthal index from which to reference the loads.',type = int)

args = parser.parse_args()


if __name__ == '__main__':
    MainDict = main()
