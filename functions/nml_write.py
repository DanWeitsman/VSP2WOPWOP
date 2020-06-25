import os
import numpy as np
# %%
def nml_write(UserIn,loadParams,dirSaveFile,i):

    if UserIn['radius'] is list:
        radius = UserIn['radius'][i]
    else:
        radius = UserIn['radius'][0]

    nOmega = UserIn['omega'][i]

    if UserIn['Vx'] is list:
        nVx = UserIn['Vx'][i]
    else:
        nVx = UserIn['Vx'][0]

    if UserIn['alphaShaft'] is list:
        alphaShaft = UserIn['alphaShaft'][i]
    else:
        alphaShaft = UserIn['alphaShaft'][0]

    nt = 2 ** np.arange(4, 25)
    ntIdeal = UserIn['nt'] * ((nOmega / 60) ** -1 * UserIn['nRev'])
    ind = [np.squeeze(np.where(nt <= ntIdeal))[-1], np.squeeze(np.where(nt >= ntIdeal))[0]]
    indSelect = np.squeeze([iter for iter, n in enumerate(abs(1 - ntIdeal / nt[ind])) if n == np.min(abs(1 - ntIdeal / nt[ind]))])
    nt = nt[ind[indSelect]]

    #%%
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
                'loadingSigmaFlag': '.false.',
                'pressureSigmaFlag': '.false.',
                'normalSigmaFlag': '.false.',
                'observerSigmaflag': '.false.',
                'loadingNoisesigmaFlag': '.false.',
                'thicknessNoiseSigmaFlag': '.false.',
                'totalNoiseSigmaFlag': '.false.',
                'accelerationSigmaFlag': '.false.',
                'MdotrSigmaFlag': '.false.',
                'velocitySigmaFlag': '.false.',
            }
    }

    environmentconstants = {
        'environmentconstants':
            {
                'rho': 1.225,
                'c': 340.2939
            }
    }

    observerin = {
        'observerin':
            {
                'title': "'Observers'",
                'attachedto': "'Aircraft'",
                'nt': nt,
                'tmin': 0.0,
                'tmax': (nOmega/60)**-1*UserIn['nRev'],
                'highpassfrequency': 2,
                'lowpassfrequency': 20000,
                'radius': radius,
                'nbtheta': UserIn['nbtheta'],
                'nbpsi': UserIn['nbpsi'],
                'thetamin': UserIn['thetamin']*(np.pi/180),
                'thetamax': UserIn['thetamax']*(np.pi/180),
                'psimin': UserIn['psimin']*(np.pi/180),
                'psimax': UserIn['psimax']*(np.pi/180),
                'nbBase': 0
            },
        'cb':
            {
                'Title': "'Observer motion'",
                'TranslationType': "'KnownFunction'",
                'VH': [nVx, 0, 0],
            }
    }

    # todo add other observer options and make user select and fill in values input module

    aircraft = {
        'containerin':
            {
                'Title': "'Aircraft'",
                'nbContainer': 1,
                'nbBase': 1,
                'dTau': (nOmega / 60)**-1/180,
            },
        'cb':
                {
                    'Title': "'Aircraft motion'",
                    'TranslationType': "'KnownFunction'",
                    'VH': [nVx, 0, 0],
                }
    }

    rotor = {
        'containerin':
            {
                'Title': "'Rotor'",
                'nbContainer': 2*UserIn['Nb'],
                'nbBase': 3,
            },
        'cb':
        [
            {
                'Title': "'Rotation'",
                'rotation': '.true.',
                'AngleType': "'KnownFunction'",
                'Omega': nOmega/60*2*np.pi,
                'Psi0': 0,
                'AxisValue': [0, 0, 1],
            },
            {
                'Title': "'Rotate to align blades with zero azimuth'",
                'AngleType': "'Timeindependent'",
                'AxisType': "'Timeindependent'",
                'AxisValue': [0, 0, 1],
                'angleValue': -np.pi/2,
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

    for Nb in range(0,UserIn['Nb']):
        blade_geom = {
            'containerin':
                {
                    'Title': "'"+'Thickness - Blade '+ str(Nb+1)+"'",
                    'patchGeometryFile': "'"+UserIn['geomFileName']+'.dat'+"'",
                    'nbBase': 3,
                },
            'cb':
                [
                    {
                        'Title': "'Constant Rotation'",
                        'AngleType': "'Timeindependent'",
                        'AxisType': "'Timeindependent'",
                        'AxisValue': [0, 0, 1],
                        'angleValue': 2*np.pi/UserIn['Nb']*Nb,
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

        #todo create zero padded dictionary entree for flap and cyclic pitch for hover cases

        nml.append(blade_geom)
        blade_load = {
            'containerin':
                {
                    'Title': "'"+'Loading-Blade '+ str(Nb+1)+"'",
                    'patchGeometryFile': "'"+UserIn['compactGeomFileName']+'.dat'+"'",
                    'patchLoadingFile': "'"+UserIn['loadingFileName']+'.dat'+"'",
                    'periodicKeyOffset': 2*np.pi/UserIn['Nb']*Nb,
                    'nbBase': 1,
                },
            'cb':
                {
                    'Title': "'Constant Rotation'",
                    'AngleType': "'Timeindependent'",
                    'AxisType': "'Timeindependent'",
                    'AxisValue': [0, 0, 1],
                    'angleValue': 2*np.pi/UserIn['Nb']*Nb,
                }
        }
        nml.append(blade_load)


    with open(os.path.expanduser(dirSaveFile+os.path.sep+UserIn['NmlFileName']), 'w') as f:
        for nml_iter in nml:
            for i, contents in enumerate(nml_iter.items()):
                f.write('&' + contents[0] + '\n')
                if type(contents[1]) is not list:
                    for entry in contents[1].items():
                            if type(entry[1]) is not list:
                                f.write(entry[0] + ' = ' + str(entry[1]) + '\n')
                            else:
                                f.write(entry[0] + ' = ' + str(entry[1])[1:-1] + '\n')
                    f.write('/' + '\n' + '\n')
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


