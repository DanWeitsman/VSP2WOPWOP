'''
VSP2WOPWOP Namelist File Write

Author: Daniel Weitsman

This function write out a namelist file for each case to its respective directory. If the users wishes to add any
other entrees to the namelist file or edit any of the default parameters it can be done here. A few things to note:
by default the 'sigmaflag' is set to '.false.', rho and c are the only environmental constants that the user can vary
from the input module. If you would like to set other environmental constants to values that differ from their
default, simply add them to that environmental constants namelist. The same goes for change of bases, if you would
like to add or omit a 'cb' namelist, simply edit the respective dictionary. If multiple of the same namelists are
used in a single container the corresponding dictionaries must be placed into a list.
'''

# %%
import os
import numpy as np


# %%
def nml_write(UserIn, loadParams, dirSaveFile, nVx, nVz, nOmega, alphaShaft,iter_geom):

    # The position of the observers can be maintained constant or varied for each geometry variant by
    # populating the list of any dimension of the observer grid in the input file with multiple comma
    # delimited values. This set of if statements checks if multiple observer position have been
    # specified in the input file.

    if UserIn['obsType'] == 2:
        if len(UserIn['xMax']) > 1:
            xMax = UserIn['xMax'][iter_geom]
        else:
            xMax = UserIn['xMax'][0]

        if len(UserIn['yMax']) > 1:
            yMax = UserIn['xMax'][iter_geom]
        else:
            yMax = UserIn['xMax'][0]

        if len(UserIn['zMax']) > 1:
            zMax = UserIn['xMax'][iter_geom]
        else:
            zMax = UserIn['xMax'][0]

    elif UserIn['obsType'] == 3:
        if len(UserIn['radius']) > 1:
            radius = UserIn['radius'][iter_geom]
        else:
            radius = UserIn['radius'][0]

    # Determines the sampling rate, as a power of 2, based on the desired sampling rate and the duration of the run
    # specified in the input module.

    nt = 2 ** np.arange(4, 25)
    ntIdeal = UserIn['nt'] * ((nOmega / 60) ** -1 * UserIn['nRev'])
    ind = [np.squeeze(np.where(nt <= ntIdeal))[-1], np.squeeze(np.where(nt >= ntIdeal))[0]]
    indSelect = np.squeeze(
        [iter for iter, n in enumerate(abs(1 - ntIdeal / nt[ind])) if n == np.min(abs(1 - ntIdeal / nt[ind]))])
    nt = nt[ind[indSelect]]

    # %%
    # Configuration of each namelist.
    environmentin = {
        'environmentin':
            {
                'nbSourceContainers': 1,
                'nbObserverContainers': 1,
                'spectrumFlag': '.false.',
                'SPLdBFLAG': '.false.',
                'SPLdBAFlag': '.false.',
                'OASPLdBAFlag': '.false.',
                'OASPLdBFlag': '.false.',
                'thicknessNoiseFlag': '.true.',
                'loadingNoiseFlag': '.true.',
                'acousticPressureFlag': '.true.',
                'totalNoiseFlag': '.true.',
                'debugLevel': 1,
                'ASCIIOutputFlag': '.true.',
                'sigmaflag': '.false.',
                'loadingSigmaFlag': '.true.',
                'pressureSigmaFlag': '.true.',
                'normalSigmaFlag': '.true.',
                'observerSigmaflag': '.true.',
                'loadingNoisesigmaFlag': '.true.',
                'thicknessNoiseSigmaFlag': '.true.',
                'totalNoiseSigmaFlag': '.true.',
                'accelerationSigmaFlag': '.true.',
                'MdotrSigmaFlag': '.true.',
                'velocitySigmaFlag': '.true.',
            }
    }

    environmentconstants = {
        'environmentconstants':
            {
                'rho': 1.225,
                'c': 340.2939
            }
    }

    if UserIn['obsType'] == 1:
        observerin = {
            'observerin':
                {
                    'title': "'Observers'",
                    'attachedto': "'Aircraft'",
                    'nt': nt,
                    'tmin': 0.0,
                    'tmax': (nOmega / 60) ** -1 * UserIn['nRev'],
                    'highpassfrequency': 10,
                    'lowpassfrequency': 20000,
                    'xLoc': UserIn['xLoc'],
                    'yLoc': UserIn['yLoc'],
                    'zLoc': UserIn['zLoc'],
                    'nbBase': 0
                },
            'cb':
                {
                    'Title': "'Observer motion'",
                    'TranslationType': "'KnownFunction'",
                    'VH': [nVx, 0, 0],
                }
        }
    elif UserIn['obsType'] == 2:
        observerin = {
            'observerin':
                {
                    'title': "'Observers'",
                    'attachedto': "'Aircraft'",
                    'nt': nt,
                    'tmin': 0.0,
                    'tmax': (nOmega / 60) ** -1 * UserIn['nRev'],
                    'highpassfrequency': 10,
                    'lowpassfrequency': 20000,
                    'nbx': UserIn['nbx'],
                    'xMin': UserIn['xMin'],
                    'xMax': xMax,
                    'nby': UserIn['nby'],
                    'yMin': UserIn['yMin'],
                    'yMax': yMax,
                    'nbz': UserIn['nbz'],
                    'zMin': UserIn['zMin'],
                    'zMax': zMax,
                    'nbBase': 0
                },
            'cb':
                {
                    'Title': "'Observer motion'",
                    'TranslationType': "'KnownFunction'",
                    'VH': [nVx, 0, 0],
                }
        }

    elif UserIn['obsType'] == 3:
        observerin = {
            'observerin':
                {
                    'title': "'Observers'",
                    'attachedto': "'Aircraft'",
                    'nt': nt,
                    'tmin': 0.0,
                    'tmax': (nOmega / 60) ** -1 * UserIn['nRev'],
                    'highpassfrequency': 10,
                    'lowpassfrequency': 20000,
                    'radius': radius,
                    'nbtheta': UserIn['nbtheta'],
                    'nbpsi': UserIn['nbpsi'],
                    'thetamin': UserIn['thetamin'] * (np.pi / 180),
                    'thetamax': UserIn['thetamax'] * (np.pi / 180),
                    'psimin': UserIn['psimin'] * (np.pi / 180),
                    'psimax': UserIn['psimax'] * (np.pi / 180),
                    'nbBase': 0
                },
            'cb':
                {
                    'Title': "'Observer motion'",
                    'TranslationType': "'KnownFunction'",
                    'VH': [nVx, 0, 0],
                }
        }

    aircraft = {
        'containerin':
            {
                'Title': "'Aircraft'",
                'nbContainer': 1,
                'nbBase': 1,
                'dTau': (nOmega / 60) ** -1 / 120,
            },
        'cb':
            {
                'Title': "'Aircraft motion'",
                'TranslationType': "'KnownFunction'",
                'VH': [nVx, 0, nVz],
            }
    }

    rotor = {
        'containerin':
            {
                'Title': "'Rotor'",
                'nbContainer': 2 * UserIn['Nb'],
                'nbBase': 3,
            },
        'cb':
            [
                {
                    'Title': "'Rotation'",
                    'rotation': '.true.',
                    'AngleType': "'KnownFunction'",
                    'Omega': nOmega / 60 * 2 * np.pi,
                    'Psi0': 0,
                    'AxisValue': [0, 0, 1],
                },
                {
                    'Title': "'Rotate to align blades with zero azimuth'",
                    'AngleType': "'Timeindependent'",
                    'AxisType': "'Timeindependent'",
                    'AxisValue': [0, 0, 1],
                    'angleValue': -np.pi / 2,
                },
                {
                    'Title': "'Rotor disk AoA tilt'",
                    'AngleType': "'Timeindependent'",
                    'AxisType': "'Timeindependent'",
                    'AxisValue': [0, 1, 0],
                    'angleValue': - alphaShaft * np.pi / 180,
                }
            ]
    }
    nml = [environmentin, environmentconstants, observerin, aircraft, rotor]

    # Loop over the blade count and appends the container of each blade to the nml list.
    for Nb in range(0, UserIn['Nb']):
        blade_geom = {
            'containerin':
                {
                    'Title': "'" + 'Thickness - Blade ' + str(Nb + 1) + "'",
                    'patchGeometryFile': "'" + UserIn['geomFileName'] + '.dat' + "'",
                    'nbBase': 3,
                },
            'cb':
                [
                    {
                        'Title': "'Constant Rotation'",
                        'AngleType': "'Timeindependent'",
                        'AxisType': "'Timeindependent'",
                        'AxisValue': [0, 0, 1],
                        'angleValue': 2 * np.pi / UserIn['Nb'] * Nb,
                    },
                    {
                        'Title': "'Pitch'",
                        'AngleType': "'Periodic'",
                        'A0': loadParams['th'][0],
                        'A1': -loadParams['th'][1],
                        'B1': -loadParams['th'][2],
                        'AxisValue': [0, 1, 0],
                    },
                    {
                        'Title': "'Coning'",
                        'AngleType': "'Timeindependent'",
                        'AxisType': "'Timeindependent'",
                        'angleValue': loadParams['beta'][0],
                        'AxisValue': [1, 0, 0],
                    }
                ]
        }

        nml.append(blade_geom)
        blade_load = {
            'containerin':
                {
                    'Title': "'" + 'Loading-Blade ' + str(Nb + 1) + "'",
                    'patchGeometryFile': "'" + UserIn['compactGeomFileName'] + '.dat' + "'",
                    'patchLoadingFile': "'" + UserIn['loadingFileName'] + '.dat' + "'",
                    'periodicKeyOffset': 2 * np.pi / UserIn['Nb'] * Nb,
                    'nbBase': 1,
                },
            'cb':
                {
                    'Title': "'Constant Rotation'",
                    'AngleType': "'Timeindependent'",
                    'AxisType': "'Timeindependent'",
                    'AxisValue': [0, 0, 1],
                    'angleValue': 2 * np.pi / UserIn['Nb'] * Nb,
                }
        }
        nml.append(blade_load)
#%%
    #   writes out each namelist to the file
    with open(os.path.expanduser(dirSaveFile + os.path.sep + UserIn['NmlFileName']), 'w') as f:
        #   iterates over each container in the nml list
        for nml_iter in nml:
            # iterates over each namelist in every container
            for i, contents in enumerate(nml_iter.items()):
                f.write('&' + contents[0] + '\n')
                #   if there is a single 'cb' namelist per container:
                if type(contents[1]) is not list:
                    #   iterates and writes out each entree in each namelist
                    for entry in contents[1].items():
                        if type(entry[1]) is not list:
                            f.write(entry[0] + ' = ' + str(entry[1]) + '\n')
                        else:
                            f.write(entry[0] + ' = ' + str(entry[1])[1:-1] + '\n')
                    f.write('/' + '\n' + '\n')
                #   if there are multiple 'cb' namelist per container:
                else:
                    for ii, list_nml in enumerate(contents[1]):
                        if ii > 0:
                            f.write('&' + contents[0] + '\n')
                        for entry in list_nml.items():
                            if type(entry[1]) is not list:
                                f.write(entry[0] + ' = ' + str(entry[1]) + '\n')
                            else:
                                f.write(entry[0] + ' = ' + str(entry[1])[1:-1] + '\n')
                        f.write('/' + '\n' + '\n')
