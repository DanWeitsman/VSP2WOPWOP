#       VSP2WOPWOP Compact Loading Functional Binary Data File Write

#   Author: Daniel Weitsman

#   This function writes a compact loading binary functional data file for a constant structured blade geometry. The comments and zoneName variables
#   are most likely the only parameters that the user may wish to vary.

#%%
import numpy as np
import struct
import os
#%%
def ConstantLoadingPatchFileWrite(loadingFileName, loadParams, nXsecs,dirSaveFile):

    aeroLoads = np.transpose([loadParams['dFx'],loadParams['dFy'],loadParams['dFz']])
    # aeroLoads = aeroLoads/np.expand_dims(loadParams['compactArea'],axis = 1)

    #%%
    fileName = loadingFileName
    magic_number = 42                #4-byte signed
    version_number = [1,0]
    comments = "Constant loading patch file"
    Nzones = 1                       # number of zones
    grid_type = 1                    # structured (1) or unstructured (2) grid
    geom_type = 1                    # constant (1), periodic (2), aperiodic (3), multi-time aperiodic (4), quasiperiodic (5), multi-time quasiperiodic (6)
    vector_centering = 1             # normal vectors are node centered (1), face centered (2)
    data_type = 2                    # data is surface pressure (1), surface loading vector (2), flow parameters (3)
    ref_frame = 3                    # reference frame is a stationary (1), rotating ground-fixed frame (2), patch-fixed frame (3). Note that this has no effect on pressure data, and that “2” and “3” are equivalent for load vectors.
    precision = 1                    # Floating points are single (1) or double (2) precision. WOPWOP only supports single
    dataZones = [1,-1]               # number of zones with data, zone designation (negative to skip thickness calc).

    zoneName = "LiftingLine"
    iMax = 1                         # number of chordwise elements
    jMax = nXsecs                    # number of spanwise elements


    #%%
    # os.chdir(dirPatchFile)

    with open(os.path.abspath(os.path.expanduser(dirSaveFile+ '/'+ fileName + '.dat')),'bw') as f_bin:

        magic_number_bin = struct.pack('<i', magic_number)
        f_bin.write(magic_number_bin)
        version_number_bin1 = struct.pack('<i', version_number[0])
        f_bin.write(version_number_bin1)
        version_number_bin2 = struct.pack('<i', version_number[1])
        f_bin.write(version_number_bin2)
        comments_bin = struct.pack('<' + str(len(comments)) + 's', bytes(comments, encoding='ascii'))
        f_bin.write(comments_bin)
        comments_pad_bin = struct.pack(str(1024 - len(comments_bin)) + 'x')
        f_bin.write(comments_pad_bin)
        func_type_bin = struct.pack('<i', 2)
        f_bin.write(func_type_bin)
        Nzones_bin = struct.pack('<i', Nzones)
        f_bin.write(Nzones_bin)
        grid_type_bin = struct.pack('<i', grid_type)
        f_bin.write(grid_type_bin)
        geom_type_bin = struct.pack('<i', geom_type)
        f_bin.write(geom_type_bin)
        vector_centering_bin = struct.pack('<i', vector_centering)
        f_bin.write(vector_centering_bin)
        data_type_bin = struct.pack('<i', data_type)
        f_bin.write(data_type_bin)
        ref_frame_bin = struct.pack('<i', ref_frame)
        f_bin.write(ref_frame_bin)
        precision_bin = struct.pack('<i', precision)
        f_bin.write(precision_bin)
        f_bin.write(struct.pack('<i', 0))
        f_bin.write(struct.pack('<i', 0))

        for i in range(len(dataZones)):
            dataZones_bin = struct.pack('<i', dataZones[i])
            f_bin.write(dataZones_bin)

        zoneName_bin = struct.pack('<32s', bytes(zoneName, encoding='ascii'))
        f_bin.write(zoneName_bin)
        iMax_bin = struct.pack('<i', iMax)
        f_bin.write(iMax_bin)
        jMax_bin = struct.pack('<i', jMax)
        f_bin.write(jMax_bin)

        for i in range(np.size(aeroLoads,1)):
            data_bin = struct.pack('<' + str(np.size(aeroLoads,0)) + 'f', *aeroLoads[:,i])
            f_bin.write(data_bin)

