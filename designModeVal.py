#       VSP2WOPWOP Operational Mode 1 Value Return

#   Author: Daniel Weitsman
#   10/28/20

# When using the design mode (operMode = 1), distinct operating conditions can be defined for each blade geometry.
# This module returns the correct quantity corresponding to each geometry variant.

def designModeVal(UserIn,XsecPolar, i):
    # In the design mode each  can be trimmed to a distinct operating condition. This set of if
    # statement checks the length of the quantities specified as lists in the input module. If the list consists of a
    # single element, this value is retained for all subsequent DegenGeom variant. However, if its length is greater
    # than one the value corresponding to the current DegenGeom iteration is returned.

    if len(UserIn['T']) > 1:
        T = UserIn['T'][i]
    else:
        T = UserIn['T'][0]

    if len(UserIn['Vz']) > 1:
        Vz = UserIn['Vz'][i]
    else:
        Vz = UserIn['Vz'][0]

    if len(UserIn['omega']) > 1:
        omega = UserIn['omega'][i]
    else:
        omega = UserIn['omega'][0]

    if len(UserIn['Vx']) > 1:
        Vx = UserIn['Vx'][i]
    else:
        Vx = UserIn['Vx'][0]

    if len(UserIn['alphaShaft']) > 1:
        alphaShaft = UserIn['alphaShaft'][i]
    else:
        alphaShaft = UserIn['alphaShaft'][0]

    if len(UserIn['airfoilPolarFileName']) > 1:
        XsecPolar = XsecPolar[list(XsecPolar.keys())[i]]
    else:
        XsecPolar = XsecPolar[list(XsecPolar.keys())[0]]

    return T,Vz,Vx,omega,alphaShaft,XsecPolar