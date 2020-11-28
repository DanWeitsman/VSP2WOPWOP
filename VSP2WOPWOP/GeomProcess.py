#       VSP2WOPWOP Blade Geometry Analysis

#   Author: Daniel Weitsman

#   This function analyzes the blade geometry to determine various geometric parameters, such as the solidity, twist, tapper distributionss.
#   These quantities are then assembled into a dictionary.

# %%
import numpy as np
from VSP2WOPWOP.NodeCenteredNorms import NodeCenteredNorms
#%%

def geomProcess(dataSorted, indHeader, loadPos, Nb,rotation):

    #   number of spanwise airfoil sections
    nXsecs = int(indHeader[0][1][1])
    #   number of points per section
    pntsPerXsec = int(indHeader[0][1][2])

    #   extracts surface nodes
    surfNodes = np.float64(dataSorted['Component 1']['SURFACE_NODE'][1:, :3])

    #   extracts LE and TE Nodes
    LENodes = np.float64(dataSorted['Component 1']['STICK_NODE'][1:, :3])
    TENodes = np.float64(dataSorted['Component 1']['STICK_NODE'][1:, 3:6])

    TE_thick = surfNodes[1::pntsPerXsec,0]-surfNodes[(pntsPerXsec-2)::pntsPerXsec,0]

    #   computes the node centered surface normals that are scaled by the area of each surface element
    ScaledNodeCenteredSurfNorms = NodeCenteredNorms(surfNodes,pntsPerXsec,nXsecs)


    #   reflects blade over y-z plane if the blade is rotating CW
    if rotation == 2:
        surfNodes[:,0] = -surfNodes[:,0]
        LENodes[:,0] = -LENodes[:,0]
        TENodes[:,0] = -TENodes[:,0]
        ScaledNodeCenteredSurfNorms[:,0] = -ScaledNodeCenteredSurfNorms[:,0]

    #   computes the chord distribution along the blade span
    chordDist = np.linalg.norm(abs(LENodes - TENodes), axis=1)

    # precone = np.arctan(LENodes[:, 2][-1]/LENodes[:, 1][-1])
    # preconez = LENodes[:, 1]*np.tan(precone)

    #   computes the twist distribution along the blade span [rad] (this quantitiy varies based on the orientation of the blade in OpenVSP)
    twistDist =np.arctan(-LENodes[:,2]/LENodes[:,0])

    # stickNodes = np.float64(comp1Data['SURFACE_FACE'][1:, :])
    # plateNodes = np.float64(comp1Data['PLATE_NODE'][1:, :3])

    # blade radius and root cutout
    R = LENodes[:,1][-1]
    e = LENodes[:,1][0]

    # Rotor disk area
    A = np.pi*R**2

    #   Blade planform chordwise distribution (used to more accurately compute the local blade solidity)
    # planformChord = TENodes[:, 0] - LENodes[:, 0]

    #   Dimensional radial vector
    rdim = np.linspace(e, R, len(chordDist))
    #   Non-dimensional radial vector
    r = np.linspace(e / R, 1, len(chordDist))

    sectLen = np.insert(np.diff(rdim), 0, np.diff(rdim)[0])

     # Local solidity computed between cross-sections
    # solDist = np.cumsum(Nb * (planformChord[:-1] + np.diff(planformChord) / 2) * np.diff(rdim)) / \
    #           (np.pi * ((rdim[:-1] + np.diff(rdim) / 2) ** 2))
    # r = np.linspace(e / R, 1, len(solDist))

    #   Local solidity computed at cross-sections
    solDist = Nb*chordDist/(np.pi*R)
    sol = np.mean(solDist)

    #   Coordinates of the lifting line
    liftLineCoord = LENodes - (LENodes - TENodes) * loadPos
    liftLineNorm = np.transpose((np.sin(twistDist),np.zeros(len(twistDist)),np.cos(twistDist)))
    # liftLineNorm = np.transpose((np.zeros(len(twistDist)), np.zeros(len(twistDist)), np.ones(len(twistDist))))
    geomParams = {'liftLineCoord':liftLineCoord,'liftLineNorm':liftLineNorm,'R':R,'e':e,'diskArea':A,'sectLen':sectLen,'chordDist':chordDist,'twistDist':twistDist,'solDist':solDist,
                  'solidity':sol,'surfNodes':surfNodes,'surfNorms':ScaledNodeCenteredSurfNorms,'nXsecs':nXsecs,'pntsPerXsec':pntsPerXsec,'rdim':rdim,'r':r,'TE_thick':TE_thick}


    return geomParams

#
    # import matplotlib.pyplot as plt
    # fig = plt.figure()
    # ax = fig.gca(projection = '3d')
    # ax.auto_scale_xyz([-2, 2], [10, 60], [-1, 1])
    # ax.pbaspect = [.09, 1, .05]
    # ax.set(xlabel = 'x',ylabel = 'y',zlabel = 'z')
    # # ax.scatter3D(surfNodes[:,0], surfNodes[:,1], surfNodes[:,2],c = 'red',linewidths = 1)
    # ax.scatter3D(LENodes[:,0], LENodes[:,1], LENodes[:,2],c = 'green',linewidths = 1)
    # ax.scatter3D(TENodes[:,0], TENodes[:,1], TENodes[:,2],c = 'blue',linewidths = 1)
    # ax.scatter3D(surfNodes[:,0], surfNodes[:,1], surfNodes[:,2], c='yellow', linewidths=.2)
    # ax.scatter3D(surfNodes[1::pntsPerXsec,0], surfNodes[1::pntsPerXsec,1], surfNodes[1::pntsPerXsec,2], c='red', linewidths=1)
    # ax.scatter3D(surfNodes[(pntsPerXsec-2)::pntsPerXsec,0], surfNodes[(pntsPerXsec-2)::pntsPerXsec,1], surfNodes[(pntsPerXsec-2)::pntsPerXsec,2], c='red', linewidths=1)
    #
    # plt.show()
    # # ax.quiver(surfNodes[:,0], surfNodes[:,1], surfNodes[:,2],ScaledNodeCenteredSurfNorms[:,0]*0.005,ScaledNodeCenteredSurfNorms[:,1]*0.005,ScaledNodeCenteredSurfNorms[:,2]*0.005)
    # ax.quiver(liftLineCoord[:,0], liftLineCoord[:,1], liftLineCoord[:,2],liftLineNorm[:,0]*0.005,liftLineNorm[:,1]*0.005,liftLineNorm[:,2]*0.005)
