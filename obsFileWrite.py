'''
This module writes an ASCII single-zoned Plot3D file containing the custom observer grid
Author: Daniel Weitsman
11/17/21
'''

def obsFileWrite(UserIn,dirSaveFile):
    import os
    import numpy as np
    x = np.array(UserIn['x'])
    y = np.array(UserIn['y'])
    z = np.array(UserIn['z'])


    header = np.ones(3)
    if np.any(np.diff( x)) != 0:
        header[0] = len(x)
    elif np.any(np.diff( y)) != 0:
        header[1] = len(y)
    elif np.any(np.diff( z)) != 0:
        header[2] = len(z)

    with open(os.path.abspath(os.path.join(dirSaveFile,'obs_grid.ascii')),'w') as f:
        f.write(str(header.astype(int))[1:-1]+'\n')
        f.write(str(x.astype(float))[1:-1]+'\n')
        f.write(str(y.astype(float))[1:-1]+ '\n')
        f.write(str(z.astype(float))[1:-1]+'\n')
