#       VSP2WOPWOP Constant BPM Binary Broadband Data File

#   Author: Daniel Weitsman

#   This function writes out the binary BPM broadband data file for constant cases.

import os
import struct
import numpy as np

#%%


def ConstantBPMWrite(geomParams,loadParams,dirSaveFile):

    TEThickness = np.ones(len(geomParams['chordDist']))*0.001
    TEflowAngle = np.zeros(len(geomParams['chordDist']))*1*(np.pi/180)

    #%%

    magic_number = 42                   # 4-byte signed
    nSect = geomParams['nXsecs']       # Number of bade sections
    uniformBlade = 0                        # Set equal to '1' if each blade section is uniform and '0' otherwise.
    sectChordFlag = 1                       # Set equal to '1' if the blade section chord is included,'0' otherwise
    sectLengthFlag = 1                         # Set equal to '1' if the blade section length is included,'0' otherwise
    TEthicknessFlag = 1                         # Set equal to '1' if the TE thickness is included,'0' otherwise
    TEflowAngleFlag = 1                         # Set equal to '1' if the TE flow angle (radians) is included,'0' otherwise
    TipLCSFlag = 0                   # Set equal to '1' if the blade section lift-curve slope is included,'0' otherwise
    SectAOAFlag = 1                         # Set equal to '1' if the blade section angle of attack (radians) is included,'0' otherwise
    UFlag = 1                           # Set equal to '1' if the blade section freestream speed is included,'0' otherwise
    timeType = 1                        # Set equal to '1' if the data is constant, 2 if periodic, 3 if aperoodic.

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


    #   Constant blade section data
        for i in range(0,nSect):
            #   sectional chord
            f_bin.write(struct.pack('>f', geomParams['chordDist'][i]))
            #   sectional length
            f_bin.write(struct.pack('>f', geomParams['sectLen'][i]))
            #   TE Thickness
            f_bin.write(struct.pack('>f', TEThickness[i]))
            #   TE Flow Angle
            f_bin.write(struct.pack('>f', TEflowAngle[i]))

        # for i in range(0,nSect):
            #   Effective sectional AoA
        f_bin.write(struct.pack('>'+str(len(loadParams['AoA']))+'f', *loadParams['AoA']))
            # f_bin.write(struct.pack('>f', loadParams['AoA'][i]))
            #   Sectional tip lift-curve slope
        # for i in range(nXsecs):
        #     f_bin.write(struct.pack('<f', loadParams['ClaDist'][i]))
        # for i in range(0,nSect):
            #   Sectional freestream velocity
        f_bin.write(struct.pack('>'+str(len(loadParams['UP']))+'f', *loadParams['UP']))

