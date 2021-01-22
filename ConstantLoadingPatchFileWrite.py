#       VSP2WOPWOP Compact Loading Functional Binary Data File Write

#   Author: Daniel Weitsman

# This function writes a compact loading binary functional data file for a constant structured blade geometry. The
# comments and zoneName variables are most likely the only parameters that the user may wish to vary.

#%%
def ConstantLoadingPatchFileWrite(loadingFileName, loadParams, nXsecs,dirSaveFile):
    #%% imports necessary modules
    import numpy as np
    import struct
    import os
    #%%
    aeroLoads = np.transpose([loadParams['dFx'],loadParams['dFy'],loadParams['dFz']])
    # aeroLoads = aeroLoads/np.expand_dims(loadParams['compactArea'],axis = 1)

    #%%
    # 4-byte signed
    magic_number = 42
    version_number = [1,0]
    comments = "Constant loading patch file"
    # number of zones, must correspond to the same number of zones as specified in the geometry patch file
    Nzones = 2
    # structured (1) or unstructured (2) grid
    grid_type = 1
    # constant (1), periodic (2), aperiodic (3), multi-time aperiodic (4), quasiperiodic (5), multi-time quasiperiodic (6)
    geom_type = 1
    # normal vectors are node centered (1), face centered (2) (Only node centered normals are supported in v3.4.4)
    vector_centering = 1
    # data is surface pressure (1), surface loading vector (2), flow parameters (3)
    data_type = 2
    # reference frame is a stationary (1), rotating ground-fixed frame (2), patch-fixed frame (3). Note that this has
    # no effect on pressure data, and that “2” and “3” are equivalent for load vectors.
    ref_frame = 3
    # Floating points are single (1) or double (2) precision  (only single precision is supported in v3.4.4)
    precision = 1
    # number of zones with data, zone designation (negative to skip thickness noise calculation)
    dataZones = [1,-2]

    zoneName = "Blade Loads"
    # number of chordwise elements
    iMax = 1
    # number of spanwise elements
    jMax = nXsecs


    #%%
    # os.chdir(dirPatchFile)

    with open(os.path.abspath(os.path.join(dirSaveFile, loadingFileName + '.dat')),'bw') as f_bin:

        f_bin.write(struct.pack('<i', magic_number))
        f_bin.write(struct.pack('<i', version_number[0]))
        f_bin.write(struct.pack('<i', version_number[1]))
        comments_bin = struct.pack('<' + str(len(comments)) + 's', bytes(comments, encoding='ascii'))
        f_bin.write(comments_bin)
        f_bin.write(struct.pack(str(1024 - len(comments_bin)) + 'x'))
        f_bin.write(struct.pack('<i', 2))
        f_bin.write(struct.pack('<i', Nzones))
        f_bin.write(struct.pack('<i', grid_type))
        f_bin.write(struct.pack('<i', geom_type))
        f_bin.write(struct.pack('<i', vector_centering))
        f_bin.write(struct.pack('<i', data_type))
        f_bin.write(struct.pack('<i', ref_frame))
        f_bin.write(struct.pack('<i', precision))
        f_bin.write(struct.pack('<i', 0))
        f_bin.write(struct.pack('<i', 0))

        for i in range(len(dataZones)):
            f_bin.write(struct.pack('<i', dataZones[i]))

        f_bin.write(struct.pack('<32s', bytes(zoneName, encoding='ascii')))
        f_bin.write(struct.pack('<i', iMax))
        f_bin.write(struct.pack('<i', jMax))

        for i in range(np.size(aeroLoads,1)):
            f_bin.write(struct.pack('<' + str(np.size(aeroLoads,0)) + 'f', *aeroLoads[:,i]))

