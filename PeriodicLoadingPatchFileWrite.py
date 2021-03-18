#       VSP2WOPWOP Compact Period Loading Functional Binary Data File Write

#   Author: Daniel Weitsman

#   This function writes a compact periodic loading binary functional data file for a structured blade geometry. The comments and zoneName variables
#   are most likely the only parameters that the user may wish to vary.


#%%
def PeriodicLoadingPatchFileWrite(loadingFileName, loadParams, nXsecs, omega,dirSaveFile):
    #%% imports necessary modules
    import numpy as np
    import struct
    import os

    #%%
    loads = [loadParams['dFx'],loadParams['dFy'],loadParams['dFz']]
    # aeroLoads = aeroLoads/np.expand_dims(loadParams['compactArea'],axis = 1)

    #%%
    magic_number = 42                #4-byte signed
    version_number = [1,0]
    comments = "Constant loading patch file"
    Nzones = 2                       # number of zones
    grid_type = 1                    # structured (1) or unstructured (2) grid
    geom_type = 2                    # constant (1), periodic (2), aperiodic (3), multi-time aperiodic (4), quasiperiodic (5), multi-time quasiperiodic (6)
    vector_centering = 1             # normal vectors are node centered (1), face centered (2)
    data_type = 2                    # data is surface pressure (1), surface loading vector (2), flow parameters (3)
    ref_frame = 3                    # reference frame is a stationary (1), rotating ground-fixed frame (2), patch-fixed frame (3). Note that this has no effect on pressure data, and that “2” and “3” are equivalent for load vectors.
    precision = 1                    # Floating points are single (1) or double (2) precision. WOPWOP only supports single
    dataZones = [1,-2]               # number of zones with data, zone designation (negative to skip thickness calc).

    zoneName = "LiftingLine"
    period = (omega/60) ** -1        # period [sec]
    nkey = len(loadParams['phi'])    # number of keys specified in radians
    iMax = 1                         # number of chordwise elements
    jMax = nXsecs                    # number of spanwise elements
    keys = loadParams['phi']         # Array of keys
    # keys = np.linspace(0, period, nkey)

#%%

    with open(os.path.abspath(os.path.expanduser(dirSaveFile+ '/'+ loadingFileName + '.dat')),'bw') as f_bin:

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
        f_bin.write(struct.pack('<f', period))
        f_bin.write(struct.pack('<i', nkey))
        f_bin.write(struct.pack('<i', iMax))
        f_bin.write(struct.pack('<i', jMax))

        for i,key in enumerate(keys):
            f_bin.write(struct.pack('<f', key))
            for ii, df in enumerate(loads):
                f_bin.write(struct.pack('<'+str(np.shape(df)[1])+'f', *df[i]))


