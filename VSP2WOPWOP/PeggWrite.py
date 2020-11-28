#       VSP2WOPWOP Pegg Method Broadband Data File Write

#   Author: Daniel Weitsman

#   This function writes a binary data file for that is needed by PSU-WOPWOP to predict the broadband noise contributions.
#%%
import struct
import os
#%%

def PeggBBDataFileWrite(geomParams,loadParams,dirSaveFile):


#%%
    magic_number = 42   # 4-byte signed
    IncBladeArea = 1    # Set equal to one in order to include blade area (4-byte signed)
    IncBladeR = 1       # Set equal to one in order to include blade radius (4-byte signed)
    IncOmega = 1        # Set equal to one in order to include the rotational rate expressed in rad/sec (4-byte signed)
    IncCL = 1           # Set equal to one in order to include the rotor's average lift coefficient (4-byte signed)
    IncHubAxis = 0      # Set equal to one in order to specify the rotor's hub axis as a vector, default is [0,0,1] (4-byte signed)
    IncThrust = 1       # Set equal to one in order to include the total thrust expressed as a vector, default is [Tx,Ty,Tz] (4-byte signed)
    DataType = 1        # Set equal to one if data is constant, two if periodic and three if aperiodic (4-byte signed)

#%%

    with open(os.path.abspath(os.path.expanduser(dirSaveFile + '/Pegg.dat')), 'bw') as f_bin:

        magic_number_bin = struct.pack('<i', magic_number)
        f_bin.write(magic_number_bin)
        IncBladeArea_bin = struct.pack('<i', IncBladeArea)
        f_bin.write(IncBladeArea_bin)
        IncBladeR_bin = struct.pack('<i', IncBladeR)
        f_bin.write(IncBladeR_bin)
        IncOmega_bin = struct.pack('<i', IncOmega)
        f_bin.write(IncOmega_bin)
        IncCL_bin = struct.pack('<i', IncCL)
        f_bin.write(IncCL_bin)
        IncHubAxis_bin = struct.pack('<i', IncHubAxis)
        f_bin.write(IncHubAxis_bin)
        IncThrust_bin = struct.pack('<i', IncThrust)
        f_bin.write(IncThrust_bin)
        DataType_bin = struct.pack('<i', DataType)
        f_bin.write(DataType_bin)

        bladeArea_bin = struct.pack('<f', geomParams['bladeArea'])
        f_bin.write(bladeArea_bin)
        bladeR_bin = struct.pack('<f', geomParams['R'])
        f_bin.write(bladeR_bin)
        bladeR_bin = struct.pack('<f', geomParams['R'])
        f_bin.write(bladeR_bin)
        omega_bin = struct.pack('<f', loadParams['omega'])
        f_bin.write(omega_bin)
        CL_bin = struct.pack('<f', loadParams['CL'])
        f_bin.write(CL_bin)
        T_bin = struct.pack('<f', loadParams['T'])
        f_bin.write(T_bin)
