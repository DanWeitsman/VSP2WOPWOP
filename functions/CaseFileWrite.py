import f90nml
import os
import numpy as np

def CaseFileWrite(caseFile):

    dictlist = []
    for i in range(0,len(caseFile)):
        temp = list(caseFile[list(caseFile.keys())[i]].values())
        dictlist.append(temp)

    keys = list(caseFile[list(caseFile.keys())[0]].keys()),
    casename = caseFile

    keys = list(caseFile['0'].keys())

    nml = {'casename': {keys[0]:dictlist[0][0],
                       keys[1]:dictlist[0][1],
                       keys[0]: dictlist[1][0],
                       keys[1]: dictlist[1][1],
                        keys[0]: dictlist[2][0],
                        keys[1]: dictlist[2][1]
                       }
           }
    # fill = {}
    # for i,n in enumerate(caseFile):
    #     fill = {**fill,**{n}}
    # nml = {'casename': {caseFile['0'],caseFile['1'],caseFile['2']}}

    with open('cases.nam', 'w') as nml_file:
       f90nml.write(nml, nml_file)