#       VSP2WOPWOP Airfoil Polar Analysis

#   Author: Daniel Weitsman

#   This function parses through the Xfoil polar files and detects the maximum coefficient of lift along with the corresponding angle of attack
#   as well as the angle of attack that corresponds to the nominal lift coefficient. The lift curve slope is also computed based on the
#   angle of attack interval specified in the input file. These quantities are assembled into a dictionary which is returned to the user.

#%%
import os
import numpy as np
import bisect
import time
#%%
def polarRead(UserIn,iii):

    dirDataFile = UserIn['dirDataFile']
    airfoilPolarFileName = UserIn['airfoilPolarFileName'][iii]
    # airfoilName = UserIn['airfoilName']
    aStart =  UserIn['aStart']
    aLength = UserIn['aLength']
    XsecPolar = {}

    #%%
    def min_max_CL_search(PolarData):
        minInd = bisect.bisect_left(PolarData[:,1], np.min(abs(PolarData[:,1])))
        maxInd = minInd + np.squeeze(np.where(np.sign(np.diff(PolarData[minInd:, 1])) == -1))[0]
        return minInd,PolarData[minInd,0], PolarData[minInd,1], PolarData[minInd,2],maxInd,PolarData[maxInd,0], PolarData[maxInd,1],PolarData[maxInd,2]

    # def SeekClMax(PolarData,minInd):
    #
    #     maxInd = minInd+np.squeeze(np.where(np.sign(np.diff(PolarData[minInd:,1])) == -1))[0]
    #
    #     # maxInd = []
    #     # for iii, n in enumerate(PolarData[intBounds:-intBounds,1]):
    #     #     if PolarData[iii,1] < n and n > PolarData[iii+2*intBounds,1]:
    #     #         maxInd.append(iii)
    #     #
    #     # maxInd = np.asarray(maxInd)
    #     # maxInd = maxInd[np.squeeze(np.where((searchBounds[1] > PolarData[maxInd,0]) == (searchBounds[0] < PolarData[maxInd,0])))]
    #     # if np.size(maxInd) > 1:
    #     #     maxInd = maxInd[np.squeeze(np.where(PolarData[maxInd, 1] == np.max(PolarData[maxInd, 1])))]
    #     return PolarData[maxInd,0], PolarData[maxInd,1],PolarData[maxInd,2],maxInd

    #%%
    for i, file in enumerate(airfoilPolarFileName):
        with open(os.path.expanduser(dirDataFile+os.path.sep+file)) as f:
            data = f.read()
        data = data.split("\n")
        dataParse = [x.split(" ")[1:] for x in data[12:-1]]

        #%%
        airfoilName = str(data[3].split()[-1])
        dataSort = np.zeros((len(dataParse),7))

        for ii,n in enumerate(dataParse):
            n = [x for x in n if x != ""]
            dataSort[ii] = np.asarray(n,dtype = np.float64)

        polar = dataSort[:,:3]

        minInd, alpha0, ClMin, CdMin, maxInd, alpha1, ClMax, CdMax = min_max_CL_search(polar)

        dx = polar[1, 0] - polar[0, 0]
        padInd = 5
        polar= polar[minInd-padInd:maxInd+padInd,:]

        polar[:,0] =  polar[:,0]%360
        polar2 = np.zeros((int((polar[0,0]-polar[-1,0]-dx)/dx),3))
        polar2[:,0] = np.arange(polar[-1,0]+dx, polar[0,0], dx)
        polar2[:, 1] = np.interp(polar2[:,0],[polar[-1,0],polar[0,0]],[polar[-1,1],polar[0,1]])
        polar2[:, 2] = CdMax
        polar = np.concatenate((polar, polar2))

        ind = [bisect.bisect_left(polar[:,0],aStart),bisect.bisect_left(polar[:, 0], aStart + aLength)]

        polar[:,0] = polar[:,0]*np.pi/180
        Cla = (polar[ind[1],1]-polar[ind[0],1])/(polar[ind[1],0]-polar[ind[0],0])
        b = polar[ind[0],1]-Cla*(polar[ind[0],0])
        y = Cla*polar[ind[0]:ind[1],0]+b

        XsecPolar = {**XsecPolar, **{airfoilName:{"Polar": polar, "Alpha0": alpha0*np.pi/180, "alphaMax": alpha1*np.pi/180, "ClMin": ClMin,"CdMin": CdMin,"ClMax": ClMax, "CdMax": CdMax,"Lift Slope": Cla}}}

        if UserIn['check'] == 1:
            import matplotlib.pyplot as plt
            polar[:,0] = polar[:,0]*180/np.pi
            zeroInd = np.squeeze(np.where(polar[:,0]==alpha0))
            maxInd = np.squeeze(np.where(polar[:,0]==alpha1))
            plt.plot(polar[zeroInd:maxInd+2*padInd, 0], polar[zeroInd:maxInd+2*padInd, 1],label=airfoilName+', '+str(round(UserIn['omega'][iii]))+'RPM')
            plt.plot(polar[ind[0]:ind[1],0],y,color='r')
            plt.scatter(alpha0, ClMin,color='g')
            plt.scatter(alpha1, ClMax,color='g')

            plt.ylabel('Lift Coefficient')
            plt.xlabel('Angle of Attack')
            plt.xlim(-10, 20)
            plt.ylim(-0.5,1.9)
            plt.grid('on')
            plt.legend()
            plt.show()

    return XsecPolar