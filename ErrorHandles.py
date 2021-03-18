'''

VSP2WOPWOP Input Module Error Handles

Author: Daniel Weitsman
8/6/20

This module performs all the error handaling on the input module to ensure that all the user input data is correct
'''

def ErrorHandles(UserIn):


    assert type(UserIn['dataFileName']) is list, "Ensure that 'dataFileNames' is specified as a comma-delimited list"
    assert type(UserIn['airfoilPolarFileName'][0]) is list, "Ensure that 'airfoilPolarFileName' is specified as a nested comma-delimited list"
    assert type(UserIn['XsecLocation']) is list, "Ensure that 'XsecLocation' is specified as a comma-delimited list"
    # assert type(UserIn['Vx']) is list, "Ensure that 'Vx' is specified as a comma-delimited list"
    # assert type(UserIn['Vz']) is list, "Ensure that 'Vz' is specified as a comma-delimited list"
    # assert type(UserIn['alphaShaft']) is list, "Ensure that 'alphaShaft' is specified as a comma-delimited list"
    # assert type(UserIn['T']) is list, "Ensure that 'T' is specified as a comma-delimited list"
    # assert type(UserIn['omega']) is list, "Ensure that 'omega' is specified as a comma-delimited list"

    if UserIn['obsType']==1:
        assert type(UserIn['xMax']) is list, "Ensure that 'xMax' is specified as a comma-delimited list"
        assert type(UserIn['yMax']) is list, "Ensure that 'yMax' is specified as a comma-delimited list"
        assert type(UserIn['zMax']) is list, "Ensure that 'zMax' is specified as a comma-delimited list"
    elif UserIn['obsType']==2:
        assert type(UserIn['radius']) is list, "Ensure that 'radius' is specified as a comma-delimited list"
