#       VSP2WOPWOP Compact Geometry Binary Patch File Write

#   Author: Daniel Weitsman

#   This function writes a compact geometry binary patch file upon which the loading vectors are prescribed. The units, comments and zoneName variables
#   are most likely the only parameters that the user may wish to vary.


#%%
import numpy as np
import struct
import os
#%%
def CompactGeomPatchFileWrite(liftLineFileName, nXsecs, liftLineCoord, liftLineNorm,dirSaveFile):
    # data = liftLineCoord
    # norm = np.append(np.zeros((np.size(data,0),2)),-np.ones((np.size(data,0),1)),1)

    #%%

    fileName = liftLineFileName
    magic_number = 42                #4-byte signed
    version_number = [1,0]
    units = 'Pa'                  #32-byte
    comments = 'Compact lifting line geometry'
    geometryFile = 1                 # If this is a geometry file = 1, -1 for psuedo geometry
    Nzones = 1                       # number of zones
    grid_type = 1                    # structured (1) or unstructured (2) grid
    geometry_type = 1                # constant (1), periodic (2), aperiodic (3), multi-time aperiodic (4), quasiperiodic (5), multi-time quasiperiodic (6)
    vector_centering = 1             # normal vectors are node centered (1), face centered (2)
    precision = 1                    # Floating points are single (1) or double (2) precision. WOPWOP only supports single
    iblank = 0                       # iblank values are included (1) or not included (0)
    zoneName = ['Lifting Line']
    iMax = 1                          # number of chordwise elements
    jMax = nXsecs            # number of spanwise elements

    #%%

    # os.chdir(dirPatchFile)

    with open(os.path.abspath(os.path.expanduser(dirSaveFile+ '/'+ fileName + '.dat')),'bw') as f_bin:

        magic_number_bin = struct.pack('<i', magic_number)
        f_bin.write(magic_number_bin)
        version_number_bin1 = struct.pack('<i', version_number[0])
        f_bin.write(version_number_bin1)
        version_number_bin2 = struct.pack('<i', version_number[1])
        f_bin.write(version_number_bin2)
        units_bin = struct.pack('<32s', bytes(units, encoding='ascii'))
        f_bin.write(units_bin)
        comments_bin = struct.pack('<'+str(len(comments))+'s',bytes(comments, encoding='ascii'))
        f_bin.write(comments_bin)
        comments_pad_bin = struct.pack(str(1024-len(comments_bin))+'x')
        f_bin.write(comments_pad_bin)
        geometryFile_bin = struct.pack('<i', geometryFile)
        f_bin.write(geometryFile_bin)
        Nzones_bin = struct.pack('<i', Nzones)
        f_bin.write(Nzones_bin)
        grid_type_bin = struct.pack('<i', grid_type)
        f_bin.write(grid_type_bin)
        geometry_type_bin = struct.pack('<i', geometry_type)
        f_bin.write(geometry_type_bin)
        vector_centering_bin = struct.pack('<i', vector_centering)
        f_bin.write(vector_centering_bin)
        precision_bin = struct.pack('<i', precision)
        f_bin.write(precision_bin)
        iblank_bin = struct.pack('<i', iblank)
        f_bin.write(iblank_bin)
        f_bin.write(struct.pack('<i', 0))

        zoneName_bin = struct.pack('<32s', bytes(zoneName[0], encoding='ascii'))
        f_bin.write(zoneName_bin)
        iMax_bin = struct.pack('<i', iMax)
        f_bin.write(iMax_bin)
        jMax_bin = struct.pack('<i', jMax)
        f_bin.write(jMax_bin)

        for i in range(np.size(liftLineCoord,1)):
            data_bin = struct.pack('<' + str(np.size(liftLineCoord,0)) + 'f',*liftLineCoord[:,i])
            f_bin.write(data_bin)

        for i in range(np.size(liftLineNorm,1)):
            norm_bin = struct.pack('<' + str(np.size(liftLineNorm,0)) + 'f',*liftLineNorm[:,i])
            f_bin.write(norm_bin)

