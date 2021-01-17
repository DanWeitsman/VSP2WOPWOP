#       VSP2WOPWOP Periodic BPM Binary Broadband Data File

#   Author: Daniel Weitsman

#   This function writes out the binary BPM broadband data file for a non-uniform blade with periodic loading data.


#%%


def PeriodicBPMWrite(geomParams,loadParams,nRev,omega,dirSaveFile):
    #%% imports necessary modules
    import os
    import struct
    import numpy as np
    #%% TEflowAngle hardcoded to 1 degree
    TEflowAngle = np.zeros(len(geomParams['chordDist']))*(np.pi/180)

    #%%

    magic_number = 42                        # 4-byte signed
    nSect = geomParams['nXsecs']             # Number of bade sections
    uniformBlade = 0                         # Set equal to '1' if each blade section is uniform and '0' otherwise.
    sectChordFlag = 1                        # Set equal to '1' if the blade section chord is included,'0' otherwise
    sectLengthFlag = 1                       # Set equal to '1' if the blade section length is included,'0' otherwise
    TEthicknessFlag = 1                      # Set equal to '1' if the TE thickness is included,'0' otherwise
    TEflowAngleFlag = 1                      # Set equal to '1' if the TE flow angle (radians) is included,'0' otherwise
    TipLCSFlag = 0                           # Set equal to '1' if the blade section lift-curve slope is included,'0' otherwise
    SectAOAFlag = 1                          # Set equal to '1' if the blade section angle of attack (radians) is included,'0' otherwise
    UFlag = 1                                # Set equal to '1' if the blade section freestream speed is included,'0' otherwise
    timeType = 2                             # Set equal to '1' if the data is constant, 2 if periodic, 3 if aperoodic.

    period = (omega/60) ** -1                #  Period of revolution
    Nsteps = loadParams['phiRes']         #  Azimuthal loading resolution
    time = np.linspace(0,nRev*period,nRev*Nsteps)
    rev_count = 0

#%%

    with open(os.path.abspath(os.path.expanduser(dirSaveFile + '/BPM.dat')), 'bw') as f_bin:

        f_bin.write(struct.pack('>i', magic_number))
        f_bin.write(struct.pack('>i', nSect))
        f_bin.write(struct.pack('>i', uniformBlade))
        f_bin.write(struct.pack('>i', sectChordFlag))
        f_bin.write(struct.pack('>i', sectLengthFlag))
        f_bin.write(struct.pack('>i', TEthicknessFlag))
        f_bin.write(struct.pack('>i', TEflowAngleFlag))
        f_bin.write(struct.pack('>i', TipLCSFlag))
        f_bin.write(struct.pack('>i', SectAOAFlag))
        f_bin.write(struct.pack('>i', UFlag))
        f_bin.write(struct.pack('>i', timeType))

        f_bin.write(struct.pack('>i', nRev))
        f_bin.write(struct.pack('>f', period))
        f_bin.write(struct.pack('>i', Nsteps))

        #   Constant blade section data
        for i in range(0, nSect):
            #   sectional chord
            f_bin.write(struct.pack('>f', geomParams['chordDist'][i]))
            #   sectional length
            f_bin.write(struct.pack('>f', geomParams['sectLen'][i]))
            #   TE Thickness
            f_bin.write(struct.pack('>f',  geomParams['TE_thick'][i]))
            #   TE Flow Angle
            f_bin.write(struct.pack('>f', TEflowAngle[i]))

        for i,t in enumerate(time):
            f_bin.write(struct.pack('>f', t))
            if i > 0 and i % Nsteps == 0:
                rev_count += 1
            #   Constant blade section data
            # for ii in range(0, nSect):
            #     f_bin.write(struct.pack('>f', loadParams['AoA'][i-rev_count*Nsteps][ii]))
            #     f_bin.write(struct.pack('>f', loadParams['U'][i-rev_count*Nsteps][ii]))

            f_bin.write(struct.pack('>' + str(np.shape(loadParams['AoA'])[1]) + 'f', *loadParams['AoA'][i-rev_count*Nsteps]))
            f_bin.write(struct.pack('>' + str(np.shape(loadParams['U'])[1]) + 'f', *loadParams['U'][i-rev_count*Nsteps]))
