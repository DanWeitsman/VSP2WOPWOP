'''

VSP2WOPWOP Input Module Error Handles

Author: Daniel Weitsman
8/6/20

This module performs all the error handling on the input module to ensure that all the user input data is correct
'''

def ErrorHandles(UserIn):

    import numpy as np

    #   asserts the data type of the input parameters
    assert type(UserIn['dataFileName']) is list, "Ensure that 'dataFileNames' is specified as a comma-delimited list"
    assert type(UserIn['airfoilPolarFileName'][0]) is list, "Ensure that 'airfoilPolarFileName' is specified as a nested comma-delimited list"
    assert type(UserIn['XsecLocation']) is list, "Ensure that 'XsecLocation' is specified as a comma-delimited list"
    assert type(UserIn['Vx']) is list, "Ensure that 'Vx' is specified as a comma-delimited list"
    assert type(UserIn['Vz']) is list, "Ensure that 'Vz' is specified as a comma-delimited list"
    assert type(UserIn['alphaShaft']) is list, "Ensure that 'alphaShaft' is specified as a comma-delimited list"
    assert type(UserIn['T']) is list, "Ensure that 'T' is specified as a comma-delimited list"
    assert type(UserIn['omega']) is list, "Ensure that 'omega' is specified as a comma-delimited list"

    if UserIn['obsType']==1:
        assert type(UserIn['xMax']) is list, "Ensure that 'xMax' is specified as a comma-delimited list"
        assert type(UserIn['yMax']) is list, "Ensure that 'yMax' is specified as a comma-delimited list"
        assert type(UserIn['zMax']) is list, "Ensure that 'zMax' is specified as a comma-delimited list"

    elif UserIn['obsType']==2:
        assert type(UserIn['radius']) is list, "Ensure that 'radius' is specified as a comma-delimited list"

    #   creates a list from the UserIn parameters that could be of variable lengths
    params = [UserIn['T'],UserIn['omega'],UserIn['alphaShaft'],UserIn['Vz'], UserIn['Vx']]
    #   determines maximum length of the input parameters
    param_len = np.array([len(x) for x in params])
    #   tests to see if any of the input parameters were specified to have an intermediate length if in OperMode 1 or 2
    if UserIn['operMode'] == 1 or UserIn['operMode'] == 2:
        if not np.all(param_len==1):
            assert not np.any((param_len > 1) == (param_len < np.max(param_len))), "In 'operMode' 1 or 3 the number of values for 'T', 'omega, 'Vx', 'Vz', 'alphaShaft' and 'airfoilPolarFileName' should be equal, ensure that this is the case if more than one value is specified for any of these parameters."
        #   duplicates list if the length of any of the input parameters is 1.
            UserIn['T'], UserIn['omega'], UserIn['alphaShaft'], UserIn['Vz'], UserIn['Vx'] = [list(np.full(np.max(param_len), x[0])) if len(x) == 1 else x for x in params]

            if len(UserIn['airfoilPolarFileName'])==1:
                [UserIn['airfoilPolarFileName'].append(UserIn['airfoilPolarFileName'][0]) for x in range(len(UserIn['omega']))]
                
    return param_len