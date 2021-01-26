import os
import f90nml
import matplotlib.pyplot as plt
import numpy as np
import h5py

def apply_to_namelist(functions, cases_directory='', cases='cases.nam'):
    """
    Apply functions in list to each case folder name in the cases namelist
    :param functions: list of functions to be applied to each folder
    :param cases_directory: directory containing cases namelist, default ''
    :param cases: name of cases namelist, default 'cases.nam'
    """
    namelist = f90nml.read(os.path.abspath(os.path.expanduser(cases_directory + '/' + cases)))
    casenames = namelist['casename']
    if type(casenames) is not list:
        casenames = [casenames]
    for casename in casenames:
        globalfoldername = casename['globalfoldername']
        if os.path.isabs(globalfoldername):
            casefolder = os.path.abspath(os.path.expanduser(casename['globalfoldername']))
        else:
            casefolder = os.path.abspath(os.path.expanduser(cases_directory + casename['globalfoldername']))
        for f in functions:
            f(casefolder)

def import_from_namelist(file, cases_directory='', cases='cases.nam'):
    """
    Apply functions in list to each case folder
    :param file: name of saved variables in npz file to be imported
    :param cases_directory: directory containing cases namelist, default ''
    :param cases: name of cases namelist, default 'cases.nam'
    """
    namelist = f90nml.read(os.path.abspath(os.path.expanduser(cases_directory + '/' + cases)))
    casenames = namelist['casename']

    extract = {}
    if type(casenames) is not list:
        casenames = [casenames]
    for casename in casenames:
        globalfoldername = casename['globalfoldername']
        if os.path.isabs(globalfoldername):
            casefolder = os.path.abspath(os.path.expanduser(casename['globalfoldername']))
        else:
            casefolder = os.path.abspath(os.path.expanduser(cases_directory + casename['globalfoldername']))

        if type(file) is list:
            temp_extract = {}
            for i,f in enumerate(file):
                data_extract = {}
                data = np.load(os.path.abspath(os.path.expanduser(casefolder + '/' + f)))
                keys = ['names', 'function_values', 'geometry_values']
                for ii, arrg in enumerate(data.files):
                    data_extract = {**data_extract, **{keys[ii]: data[arrg]}}
                temp_extract = {**temp_extract, **{f[:-4]: data_extract}}
            extract = {**extract, **{globalfoldername[2:-1]: temp_extract}}

        else:
            data_extract = {}
            data = np.load(os.path.abspath(os.path.expanduser(casefolder+'/'+file)))
            keys = ['names','function_values','geometry_values']
            for i, arrg in enumerate(data.files):
                data_extract = {**data_extract,**{keys[i]:data[arrg]}}
            extract = {**extract,**{globalfoldername[2:-1]:data_extract}}

    return extract


def read_wopwop_output(prefix, name_postfix='.nam', function_postfix='.fn', geometry_postfix='.x', save_npz= '' , save_h5= '' ):
    """
    Read ASCII PSU-WOPWOP data from files
    :param save_npz: set = 1 in order to write the data to an npz file
    :param prefix: File path + prefix, e.g., /Path/to/case/pressure
    :param name_postfix: postfix for namelist (default .nam)
    :param function_postfix: postfix for function file (default .fn)
    :param geometry_postfix: postfix for geometry file (default .x)
    :return: Tuple containing: List of output function names,
    function values (observers1 x observers2 x samples x function types), and
    observer (X,Y,Z) positions (observers1 x observers 2 x samples x 3)
    """
    # Read in data names from namelist
    with open(os.path.abspath(os.path.expanduser(prefix  +name_postfix)), 'r') as f:
        names = f.read().splitlines()
    # Read in values from function file
    with open(os.path.abspath(os.path.expanduser(prefix  +function_postfix)), 'r') as f:
        function_data = f.read().splitlines()
        header = np.array(function_data[0].split()).astype(int)
        assert (len(header) == 4)
        nf_observers1 = header[0]
        nf_observers2 = header[1]
        nf_samples = header[2]
        nf_types = header[3]
        assert (len(names) == nf_types)
        function_values = np.squeeze(np.array([x.split() for x in function_data[1:]]).astype(float))
        # Reshape into numpy array (observers1,observers2,samples,types) using Fortran index ordering
        function_values = function_values.reshape(nf_observers1, nf_observers2, nf_samples, nf_types, order='F')
    with open(os.path.abspath(os.path.expanduser(prefix +geometry_postfix)), 'r') as f:
        geometry_data = f.read().splitlines()
        header = np.array(geometry_data[0].split()).astype(int)
        assert (len(header) == 3)
        ng_observers1 = header[0]
        ng_observers2 = header[1]
        ng_samples = header[2]
        assert (ng_observers1 == nf_observers1)
        assert (ng_observers2 == nf_observers2)
        assert (ng_samples == nf_samples)
        geometry_values = np.squeeze(np.array([x.split() for x in geometry_data[1:]]).astype(float))
        # Reshape into numpy array ((observers1,observers2,samples,3), where last dimension is X,Y,Z coordinates
        geometry_values = geometry_values.reshape(ng_observers1, ng_observers2, ng_samples, 3, order='F')
    if save_npz:
        np.savez(os.path.abspath(os.path.expanduser(prefix+name_postfix[:-4])) ,names, function_values, geometry_values)
    if save_h5:
        data = {'names':names,'geometry_values':geometry_values,'function_values':function_values}
        names = [elem.encode() for elem in names]
        with h5py.File(os.path.abspath(os.path.expanduser(prefix+name_postfix[:-4]+'.h5')),'w') as h5_file:
            # grp = h5_file.create_group('data')
            for k,v in data.items():
                h5_file.create_dataset(k,shape=np.shape(v), data=v)
    return names, function_values, geometry_values

def extract_wopwop_quant(case_directory,prefix):

    """
    Generates a single .npz file for the desired output files
    :param case_directory: Folder containing OASPL output files
    :param prefix: Metric name to be appended to namelist and function files
    """

    read_wopwop_output(prefix = os.path.abspath(os.path.expanduser(case_directory + '/')),
                                                                         name_postfix='/'+prefix +'.nam',
                                                                         function_postfix='/'+prefix + '.fn',
                                                                         geometry_postfix = '/'+prefix + '.x',save_h5=True)


def plot_pressure_time_histories(case_directory, remove_dc=True, subplots=False, savename=None):
    """
    Plot pressure time history data
    :param case_directory: Folder containing pressure.x, pressure.fn, pressure.nam output
    :param remove_dc: Set to remove DC component of signals (default True)
    :param subplots: Set to plot each time history in a subplot of the same main plot (default False)
    :param savename: Optional filename to save images.  Image file type inferred from extension.
    When subplots=false, observer indicies are appended to filename
    """
    names, function_values, geometry_values = read_wopwop_output(prefix=case_directory + '/pressure')
    shape = function_values.shape
    # Remove DC component from all channels
    if remove_dc:
        means = np.mean(function_values[:, :, :, 1:], axis=2)
        mean_shape = list(shape)
        mean_shape[2] = 1
        mean_shape[3] = mean_shape[3] - 1
        function_values[:, :, :, 1:] = function_values[:, :, :, 1:] - np.reshape(means, mean_shape)
    # Set up subplot axes (still need to work on formatting)
    if subplots:
        ni = shape[0]
        nj = shape[1]
        fig, axs = plt.subplots(nj, ni)
    # Set "nice" y-scale values at around the most significant digit of the max value
    scale_value = float('%.0e' % float(1.5 * np.max(np.abs(function_values[:, :, :, 1:]))))

    # Iterate through observers
    for oi, obs_row_functions in enumerate(function_values):
        for oj, obs_functions in enumerate(obs_row_functions):
            if not subplots:
                fig, ax = plt.subplots(1, 1)
            else:
                ax = axs[oj * ni + oi]
            if remove_dc:
                obs_functions = obs_functions - np.mean(obs_functions, axis=0)
            for n, name in enumerate(names[1:]):
                ax.plot(obs_functions[:, 0], obs_functions[:, n + 1], label=name)
                ax.set_xlabel('Time, s')
                ax.set_ylabel('Pressure, Pa')
                ax.set_ylim(-scale_value, scale_value)
                ax.legend(loc='lower center')
            if not subplots and savename is not None:
                (root,ext) = os.path.splitext(savename)
                fig.savefig(root + '_' + str(oi) + '_' + str(oj) + ext)
    if subplots and savename is not None:
        fig.savefig(savename)


def plot_oaspl_polar(case_directory, metric='dB', axis=0, savename=None):
    """
    Plot polar plot of sound pressure levels
    :param case_directory: Folder containing OASPL output files
    :param metric: Metric name to be appended to namelist and function files
    :param axis: Axis about which to calculate polars
    :param savename: Optional filename to save images.  Image file type inferred from extension.
    """

    names, function_values, geometry_values = read_wopwop_output(case_directory + '/OASPL',
                                                                 name_postfix=metric + '.nam',
                                                                 function_postfix=metric + '.fn',save_npz =1)

    # Calculate angle about selected axis
    sweep_axes = np.delete(range(3), axis)
    angles = np.squeeze(np.arctan2(geometry_values[:, :, :, sweep_axes[1]], geometry_values[:, :, :, sweep_axes[0]]))
    # Setup polar plot
    fig, ax = plt.subplots(1, 1, subplot_kw=dict(polar=True))
    for n, name in enumerate(names[1:]):
        ax.plot(angles, np.squeeze(function_values[:, :, :, n + 1]), label=name)
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.25))
    if savename is not None:
        fig.savefig(savename, bbox_inches='tight')


def plot_spectra(case_directory, subplots=False, savename=None, semilogx=False):
    names, function_values, geometry_values = read_wopwop_output(prefix=case_directory + '/spl_spectrum')
    shape = function_values.shape
    """
    Plot sound pressure level spectra
    :param case_directory: Folder containing spl_spectrum.x, spl_spectrum.fn, spl_spectrum.nam output
    :param subplots: Set to plot each time history in a subplot of the same main plot (default False)
    :param savename: Optional filename to save images.  Image file type inferred from extension.  
     When subplots=false, observer indicies are appended to filename
    :param semilogx: Set to plot x-axis on log scale (default False)
    """
    # Set up subplot axes (still need to work on formatting)
    if subplots:
        ni = shape[0]
        nj = shape[1]
        fig, axs = plt.subplots(nj, ni)
    # Set "nice" y-scale values at around the most significant digit of the max value
    scale_value = float('%.0e' % float(1.5 * np.max(function_values[:, :, :, 1:])))

    # Iterate through observers
    for oi, obs_row_functions in enumerate(function_values):
        for oj, obs_functions in enumerate(obs_row_functions):
            if not subplots:
                fig, ax = plt.subplots(1, 1)
            else:
                ax = axs[oj * ni + oi]
            # Only plot the total
            if semilogx:
                ax.semilogx(obs_functions[:, 0], obs_functions[:, -1])
            else:
                ax.plot(obs_functions[:, 0], obs_functions[:, -1])
            ax.set_xlabel('Frequency, Hz')
            ax.set_ylabel('Sound Pressure Level, dB')
            ax.set_ylim(0, scale_value)

            if not subplots and savename is not None:
                (root, ext) = os.path.splitext(savename)
                fig.savefig(root + '_' + str(oi) + '_' + str(oj) + ext)
    if subplots and savename is not None:
        fig.savefig(savename)
