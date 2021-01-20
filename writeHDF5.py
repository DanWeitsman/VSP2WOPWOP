'''VSP2WOPWOP HDF5 Write Module
Author: Daniel Weitsman
1/16/21
This function writes the main dictionary (MainDict) to an HDF5 file. HDF5 is a universal file format that is structured
like MATLAB .mat files. The content of HDF5 files can be read using the following function in the h5py package:
h5py.File("path+file"), 'r').
'''
def writeHDF5(MainDict,UserIn):

    #%% import necessary modules
    import os
    import h5py
    import numpy as np
    #%%
    with h5py.File(os.path.abspath(os.path.join(UserIn['dirPatchFile'], 'MainDict.h5')), 'w') as f_write:

        def encode_str_list(str_list):
            for i, elem in enumerate(str_list):
                if isinstance(elem, list):
                    encode_str_list(elem)
                elif isinstance(elem, str):
                    str_list[i] = elem.encode()
                else:
                    break
            return str_list

        def write_dict_hdf5(d, parent=''):
            for key, value in d.items():
                if isinstance(value, dict):
                    write_dict_hdf5(value, parent + '/' + key)
                else:
                    # print(key)
                    if isinstance(value, list):
                        value = encode_str_list(value)
                    elif isinstance(value, str):
                        value = value.encode()
                    f_write.create_dataset(parent + '/' + key, shape=np.shape(value), data=value)

        write_dict_hdf5(MainDict)