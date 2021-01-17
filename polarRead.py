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

        ind = [bisect.bisect_left(polar[:,0],aStart),bisect.bisect_left(polar[:, 0], aStart + aLength)]

        polar[:,0] = polar[:,0]*np.pi/180
        Cla = (polar[ind[1],1]-polar[ind[0],1])/(polar[ind[1],0]-polar[ind[0],0])
        b = polar[ind[0],1]-Cla*(polar[ind[0],0])
        y = Cla*polar[ind[0]:ind[1],0]+b

        XsecPolar = {**XsecPolar, **{airfoilName:{"Polar": polar, "Alpha0": alpha0*np.pi/180, "alphaMax": alpha1*np.pi/180, "ClMin": ClMin,"CdMin": CdMin,"ClMax": ClMax, "CdMax": CdMax,"Lift Slope": Cla}}}

        if UserIn['check'] == 1:
            import matplotlib.pyplot as plt

            color = ['tab:blue','tab:red']
            fig,ax1 = plt.subplots(1,1,figsize = (6.4,4.5))
            #,label=airfoilName+', '+str(round(UserIn['omega'][iii]
            ax1.plot(polar[:,0]*180/np.pi, polar[:, 1],color = color[0],label=airfoilName+', '+str(round(UserIn['omega'][iii]))+' rpm')
            ax1.plot(polar[ind[0]:ind[1],0]*180/np.pi,y,color='k')
            ax1.tick_params(axis='y', labelcolor=color[0])
            ax1.scatter(alpha0, ClMin,color='g')
            ax1.scatter(alpha1, ClMax,color='g')
            ax1.set_ylabel('Lift Coefficient',color = color[0])
            ax1.set_xlabel('Angle of Attack')
            ax1.set_xlim(-10, 20)
            ax1.set_ylim(-0.5,1.9)
            ax2 = ax1.twinx()
            ax2.plot(polar[:,0]*180/np.pi, polar[:, 2],color = color[1],label=airfoilName+', '+str(round(UserIn['omega'][iii]))+' rpm')
            ax2.set_ylabel('Drag Coefficient',color = color[1])
            ax2.set_ylim(0,0.125)
            ax2.tick_params(axis='y', labelcolor=color[1])
            plt.grid('on')
            plt.legend()
            plt.show()

    return XsecPolar