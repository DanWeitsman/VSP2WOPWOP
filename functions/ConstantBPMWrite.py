#       VSP2WOPWOP BPM Binary Broadband Data File

#   Author: Daniel Weitsman

#   This function writes out the binary BPM broadband data file for constant cases.

#todo edit nml_write to include pegg/bpm flag if write is enabled in the input module
import os
import struct


#%%


def ConstantBPMWrite(geomParams,loadParams,dirSaveFile):

    #%%

    magic_number = 42                   # 4-byte signed
    nXsecs = geomParams['nXsecs']       # Number of bade sections
    XsecType = 0                        # Set equal to '1' if each blade section is uniform and '0' otherwise.
    XsecChord = 1                       # Set equal to '1' if the blade section chord is included,'0' otherwise
    XsecLen = 1                         # Set equal to '1' if the blade section length is included,'0' otherwise
    thickTE = 0                         # Set equal to '1' if the TE thickness is included,'0' otherwise
    alphaTE = 0                         # Set equal to '1' if the TE flow angle (radians) is included,'0' otherwise
    XsecAoA = 1                         # Set equal to '1' if the blade section angle of attack (radians) is included,'0' otherwise
    XsecLiftSlope = 1                   # Set equal to '1' if the blade section lift-curve slope is included,'0' otherwise
    XsecU = 1                           # Set equal to '1' if the blade section freestream speed is included,'0' otherwise
    dataType = 1                        # Set equal to '1' if the data is constant, 2 if periodic, 3 if aperoodic.

#%%

    with open(os.path.abspath(os.path.expanduser(dirSaveFile + '/BPM.dat')), 'bw') as f_bin:

        f_bin.write(struct.pack('<i', magic_number))
        f_bin.write(struct.pack('<i', nXsecs))
        f_bin.write(struct.pack('<i', XsecType))
        f_bin.write(struct.pack('<i', XsecChord))
        f_bin.write(struct.pack('<i', XsecLen))
        f_bin.write(struct.pack('<i', thickTE))
        f_bin.write(struct.pack('<i', alphaTE))
        f_bin.write(struct.pack('<i', XsecAoA))
        f_bin.write(struct.pack('<i', XsecLiftSlope))
        f_bin.write(struct.pack('<i', XsecU))
        f_bin.write(struct.pack('<i', dataType))

    #todo if does not work, create two for loops to seperate constant and varying variables

    #   Constant blade section data
        for i in range(nXsecs):
            #   sectional chord
            f_bin.write(struct.pack('<f', geomParams['chordDist'][i]))
            #   sectional length
            f_bin.write(struct.pack('<f', geomParams['sectLen'][i]))
            #   TE Thickness
            #f_bin.write(struct.pack('<f', geomParams['sectLen'][i]))
            #   TE Flow Angle
            #f_bin.write(struct.pack('<f', geomParams['sectLen'][i]))
            #   Effective sectional AoA
            f_bin.write(struct.pack('<f', loadParams['AoA'][i]))
            #   Sectional tip lift-curve slope
            f_bin.write(struct.pack('<f', loadParams['ClaDist'][i]))
            #   Sectional freestream velocity
            f_bin.write(struct.pack('<f', loadParams['UP'][i]))

