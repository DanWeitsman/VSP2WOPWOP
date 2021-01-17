#       VSP2WOPWOP OpenVSP DegenGeom File Parse

#   Author: Daniel Weitsman

#   This function parses through the DegenGeom csv file and places the degenerate geometry components in a dictionary.
#   This is an independent function that supports any arbitrary number of components, which do not necessarily need to be rotor blades.

#%%

def AnalyzeDegenGeom(dirDataFile, dataFileName):

#%% imports necessary modules
    import os
    import csv
    import numpy as np
#%%
    data = {}
    n = 0
    # Opens and reads csv data
    with open(os.path.abspath(os.path.expanduser(dirDataFile+'/'+dataFileName)), 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        # Detects and deletes empty rows
        for row in csvreader:
            if len(row) != 0:
                data[n] = row
                n = n + 1

    # %%
    # Initializes a list for each degenerate geometry type
    indHeader = []
    # Names of the degenerate geometries
    DegenGeom = ["SURFACE_NODE", "SURFACE_FACE", "PLATE", "STICK_NODE", "STICK_FACE", "POINT"]

    #   Iterates through the data and locates the indices corresponding to the header of each degenerate geometry type
    for i, n in enumerate(data.items()):
        if any(x == n[1][0] for x in DegenGeom):
            indHeader.append(n)

    # %%
    #   This function extracts and returns a segment of data between two specified indices
    #   dataset: the data set
    #   startInd: starting index
    #   endInd: ending index

    def dictExtractasArray(dataset, startInd, endInd):
        temp = {}
        for i in range(startInd, endInd):
            temp[i] = dataset[i]
        temp = np.array(list(temp.values()))
        return temp

    # %%
    #   This section of the code extracts the degenerate geometries corresponding to each component, ...
    #   converts the degenerate geometries parameters into a numpy array, which are then organized ...
    #   in a nested dictionary. Each component has its own dictionary, which can be accessed by: ...
    #   dataSorted['"Component Name(numbered sequentially by default)"']['"Degenerate Geometry Name"']

    #   Number if components
    nComponents = int(data[2][0])
    #   Initializes two empty dictionaries, dataSorted is the dictionary that contains the final sorted data...
    #   dataTemp is a temporary dictionary used to store the degenerate geometry data of each individual component.
    dataTemp = {}
    dataSorted = {}

    for ii in range(nComponents):
        for i, n in enumerate(indHeader[0 + len(DegenGeom) * ii:len(DegenGeom) + len(DegenGeom) * ii]):

            if n[1][0] == DegenGeom[1]:
                dataTemp = {**dataTemp, **{
                    n[1][0]: dictExtractasArray(data, n[0] + 1, indHeader[i + len(DegenGeom) * ii + 1][0] - 1)}}

            elif n[1][0] == DegenGeom[2]:
                dataTemp = {**dataTemp, **{"PLATE_NORM": dictExtractasArray(data, n[0] + 1, n[0] + int(n[1][1]) + 2)}}
                dataTemp = {**dataTemp, **{"PLATE_NODE": dictExtractasArray(data, n[0] + int(n[1][1]) + 2,
                                                                            indHeader[i + len(DegenGeom) * ii + 1][
                                                                                0] - 1)}}

            elif n[1][0] == DegenGeom[3]:
                data[n[0] + 1] = data[n[0] + 1][0:-1]
                dataTemp = {**dataTemp, **{
                    n[1][0]: dictExtractasArray(data, n[0] + 1, indHeader[i + len(DegenGeom) * ii + 1][0] - 1)}}

            elif n[1][0] == DegenGeom[4]:
                dataTemp = {**dataTemp, **{
                    n[1][0]: dictExtractasArray(data, n[0] + 1, indHeader[i + len(DegenGeom) * ii + 1][0] - 1)}}

            elif n[1][0] == DegenGeom[5]:
                dataTemp = {**dataTemp, **{n[1][0]: dictExtractasArray(data, n[0] + 1, n[0] + 3)}}

            else:
                dataTemp = {**dataTemp,
                            **{n[1][0]: dictExtractasArray(data, n[0] + 1, indHeader[i + len(DegenGeom) * ii + 1][0])}}

        dataSorted = {**dataSorted, **{'Component ' + str(ii + 1): dataTemp}}

    return dataSorted, indHeader
