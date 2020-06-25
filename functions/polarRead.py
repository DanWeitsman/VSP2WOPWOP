#       VSP2WOPWOP Airfoil Polar Analysis

#   Author: Daniel Weitsman

#   This function parses through the Xfoil polar files and detects the maximum coefficient of lift along with the corresponding angle of attack
#   as well as the angle of attack that corresponds to the nominal lift coefficient. The lift curve slope is also computed based on the
#   angle of attack interval specified in the input file. These quantities are assembled into a dictionary which is returned to the user.

#%%
import os
import numpy as np

#%%
def polarRead(UserIn,iii):

    dirDataFile = UserIn['dirDataFile']
    airfoilPolarFileName = UserIn['airfoilPolarFileName'][iii]
    # airfoilName = UserIn['airfoilName']
    aStart =  UserIn['aStart']
    aLength = UserIn['aLength']

    #%%
    XsecPolar = {}
    os.chdir(dirDataFile)

    #%%
    def SeekAbsMin(PolarData):
        data = PolarData[:,1]
        Min = np.min(abs(data))
        MinInd = np.squeeze(np.where(data == Min))
        if np.size(MinInd) == 0:
            MinInd = np.squeeze(np.where(data == -Min))
            MinInd = int(MinInd)
            Min = -Min
        else:
            MinInd =int(MinInd)
        return PolarData[MinInd,0], Min, PolarData[MinInd,2]

    def SeekClMax(PolarData,intBounds,searchBounds):
        MaxInd = []
        for iii, n in enumerate(PolarData[intBounds:-intBounds,1]):
            if PolarData[iii,1] < n and n > PolarData[iii+2*intBounds,1]:
                MaxInd.append(iii)

        MaxInd = np.asarray(MaxInd)
        MaxInd = MaxInd[np.squeeze(np.where((searchBounds[1] > PolarData[MaxInd,0]) == (searchBounds[0] < PolarData[MaxInd,0])))]
        if np.size(MaxInd) > 1:
            MaxInd = MaxInd[np.squeeze(np.where(PolarData[MaxInd, 1] == np.max(PolarData[MaxInd, 1])))]
        return PolarData[MaxInd,0], PolarData[MaxInd,1],PolarData[MaxInd,2]

    #%%
    for i,file in enumerate(airfoilPolarFileName):
        with open(file) as f:
            data = f.read()
        data = data.split("\n")
        dataParse = [x.split(" ")[1:] for x in data[12:-1]]

        #%%
        airfoilName = str(data[3].split()[-1])
        dataSort = np.zeros((len(dataParse),7))

        for ii,n in enumerate(dataParse):
            n = [x for x in n if x != ""]
            dataSort[ii] = np.asarray(n,dtype = np.float64)

        polar = dataSort[:,0:3]

        [alpha0, ClMin,CdMin] = SeekAbsMin(polar)
        [alpha1,ClMax,CdMax] = SeekClMax(polar,4,[1,20])

        ind = np.squeeze(np.where(polar[:,0] == aStart) + np.where(polar[:, 0] == aStart + aLength))
        polar[:,0]= polar[:,0]*(np.pi/180)
        Cla = (polar[ind[1],1]-polar[ind[0],1])/(polar[ind[1],0]-polar[ind[0],0])

        b = polar[ind[0],1]-Cla*(polar[ind[0],0])
        y = Cla*polar[ind[0]:ind[1],0]+b

        XsecPolar = {**XsecPolar, **{airfoilName:{"Polar": polar, "Alpha0": alpha0, "ClMaxAlpha": alpha1, "ClMin": ClMin,"CdMin": CdMin,"ClMax": ClMax, "CdMax": CdMax,"Lift Slope": Cla}}}

        if UserIn['check'] == 1:
            import matplotlib.pyplot as plt

            plt.plot(XsecPolar[airfoilName]['Polar'][:, 0] * (180 / np.pi), XsecPolar[airfoilName]['Polar'][:, 1],label=airfoilName+', '+str(round(UserIn['omega'][iii]))+'RPM')
            plt.plot(polar[ind[0]:ind[1],0]*(180/np.pi),y,color='r')
            plt.scatter(alpha0*(np.pi/180)*(180/np.pi), ClMin,color='g')
            plt.scatter(alpha1*(np.pi/180*(180/np.pi)), ClMax,color='g')

            plt.ylabel('Lift Coefficient')
            plt.xlabel('Angle of Attack')
            plt.xlim(-10, 20)
            plt.ylim(-0.5,1.9)
            plt.grid('on')
            plt.legend()
            plt.show()

    return XsecPolar