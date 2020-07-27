import os

def caseFile_write(globalFolderName, caseFileName, dirSaveFile):
    with open(os.path.expanduser(dirSaveFile + os.path.sep + 'cases.nam'), 'w') as f:

        if type(globalFolderName) is list:
            for globFold in globalFolderName:
                f.write('&casename\n' + 'globalFolderName = ' + "'./" + globFold + "/'\n" + "caseNameFile = '" + caseFileName + "'\n/\n")
        else:
            f.write('&casename\n'+ 'globalFolderName = ' + "'./" + globalFolderName[0] + "/'\n" + "caseNameFile = '" + caseFileName + "'\n/\n")