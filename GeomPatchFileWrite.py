#       VSP2WOPWOP Blade Geometry and Compact Lifting Line Binary Patch File Write

#   Author: Daniel Weitsman

# This function writes a binary patch file that specifies the constant structured blade geometry and the compact
# lifting line upon which the blade loads are prescribed.


#%%
def GeomPatchFileWrite(geomFileName, geomParams, dirSaveFile):

    #%% imports necessary modules
    import struct
    import numpy as np
    import os
    #%%

    fileName = geomFileName
    magic_number = 42                                               # 4-byte signed
    version_number = [1,0]                                          # Version number 1.0
    units = 'Pa'                                                    #32-byte
    comments = 'Blade and lifting line geometry patch file'
    geometryFile = 1                                                # If this is a geometry file = 1, -1 for psuedo geometry
    Nzones = 2                                                      # Number of zones (blade geometry and compact lifting line)
    grid_type = 1                                                   # Structured (1) or unstructured (2) grid
    geometry_type = 1                                               # Constant (1), periodic (2), aperiodic (3), multi-time aperiodic (4), quasiperiodic (5), multi-time quasiperiodic (6)
    vector_centering = 1                                            # Normal vectors are node centered (1), face centered (2) (Only node centered normals are supported by v3.4.4)
    precision = 1                                                   # Floating points are single (1) or double (2) precision. WOPWOP only supports single
    iblank = 0                                                      # iblank values are included (1) or not included (0)
    zoneName = ['Blade','Lifting Line']                             # Zone names
    iMax = [geomParams['pntsPerXsec'], 1]                           # number of chordwise elements
    jMax = [geomParams['nXsecs'],geomParams['nXsecs']-1]              # number of spanwise elements

    #%%

    with open(os.path.abspath(os.path.join(dirSaveFile,fileName + '.dat')),'bw') as f_bin:

        f_bin.write(struct.pack('<i', magic_number))
        f_bin.write(struct.pack('<i', version_number[0]))
        f_bin.write(struct.pack('<i', version_number[1]))
        f_bin.write(struct.pack('<32s', bytes(units, encoding='ascii')))
        comments_bin = struct.pack('<' + str(len(comments)) + 's', bytes(comments, encoding='ascii'))
        f_bin.write(comments_bin)
        f_bin.write(struct.pack(str(1024 - len(comments_bin)) + 'x'))
        f_bin.write(struct.pack('<i', geometryFile))
        f_bin.write(struct.pack('<i', Nzones))
        f_bin.write(struct.pack('<i', grid_type))
        f_bin.write(struct.pack('<i', geometry_type))
        f_bin.write(struct.pack('<i', vector_centering))
        f_bin.write(struct.pack('<i', precision))
        f_bin.write(struct.pack('<i', iblank))
        f_bin.write(struct.pack('<i', 0))

        #structured header
        for i in range(Nzones):
            f_bin.write(struct.pack('<32s', bytes(zoneName[i], encoding='ascii')))
            f_bin.write(struct.pack('<i', iMax[i]))
            f_bin.write(struct.pack('<i', jMax[i]))

        for i in range(np.size(geomParams['surfNodes'],1)):
            f_bin.write(struct.pack('<' + str(np.size(geomParams['surfNodes'],0)) + 'f',*geomParams['surfNodes'][:,i]))

        for i in range(np.size(geomParams['surfNorms'],1)):
            f_bin.write(struct.pack('<' + str(np.size(geomParams['surfNorms'],0)) + 'f',*geomParams['surfNorms'][:,i]))

        for i in range(np.size(geomParams['liftLineCoord'],1)):
            f_bin.write(struct.pack('<' + str(np.size(geomParams['liftLineCoord'],0)) + 'f',*geomParams['liftLineCoord'][:,i]))

        for i in range(np.size(geomParams['liftLineNorm'],1)):
            f_bin.write(struct.pack('<' + str(np.size(geomParams['liftLineNorm'],0)) + 'f',*geomParams['liftLineNorm'][:,i]))
