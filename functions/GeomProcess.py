#       VSP2WOPWOP Blade Geometry Analysis

#   Author: Daniel Weitsman

#   This function analyzes the blade geometry to determine various geometric parameters, such as the solidity, twist, tapper distributionss.
#   These quantities are then assembled into a dictionary.

# %%
import numpy as np
from functions.NodeCenteredNorms import NodeCenteredNorms
#%%

def geomProcess(dataSorted, indHeader, loadPos, Nb):


    comp1Data = dataSorted['Component 1']

    nXsecs = int(indHeader[0][1][1])
    pntsPerXsec = int(indHeader[0][1][2])

    surfNodes = np.float64(comp1Data['SURFACE_NODE'][1:, :3])

    surfNorms = np.float64(comp1Data['SURFACE_FACE'][1:, :])
    #   Duplicates two most inboard sets of normal vectors and areas to node-center vectors (required for sigma surfaces)
    FaceCenteredSurfNorms = surfNorms[:, :3] / np.expand_dims(surfNorms[:, 3], axis=1)

    #   Extracts node centered surface normal vectors that are normalized by the face areas
    ScaledNodeCenteredSurfNorms = NodeCenteredNorms(surfNodes,pntsPerXsec,nXsecs)

    #   LE and TE Nodes
    LENodes = np.float64(comp1Data['STICK_NODE'][1:, :3])
    TENodes = np.float64(comp1Data['STICK_NODE'][1:, 3:6])

    #   Chord distribution along the blade span
    chordDist = np.linalg.norm(LENodes - TENodes, axis=1)

    # precone = np.arctan(LENodes[:, 2][-1]/LENodes[:, 1][-1])
    # preconez = LENodes[:, 1]*np.tan(precone)

    #   Twist distribution along the blade span in radians (this quantitiy varies based on the orientation of the blade in OpenVSP)
    twistDist =np.arctan(-LENodes[:,2]/LENodes[:,0])
    # stickNodes = np.float64(comp1Data['SURFACE_FACE'][1:, :])
    # plateNodes = np.float64(comp1Data['PLATE_NODE'][1:, :3])

    # Blade radius and root cutout
    R = np.linalg.norm(LENodes[-1])
    e = np.linalg.norm(LENodes[0])

    # Rotor disk area
    A = np.pi*R**2

    #   Blade planform chordwise distribution (used to more accurately compute the local blade solidity)
    planformChord = TENodes[:, 0] - LENodes[:, 0]

    #   Dimensional radial vector
    rdim = np.linspace(e, R, len(planformChord))
    #   Non-dimensional radial vector
    r = np.linspace(e / R, 1, len(planformChord))

    sectLen = np.insert(np.diff(rdim), 0, np.diff(rdim)[0])

     # Local solidity computed between cross-sections
    # solDist = np.cumsum(Nb * (planformChord[:-1] + np.diff(planformChord) / 2) * np.diff(rdim)) / \
    #           (np.pi * ((rdim[:-1] + np.diff(rdim) / 2) ** 2))
    # r = np.linspace(e / R, 1, len(solDist))

    #   Local solidity computed at cross-sections
    solDist = np.cumsum(Nb * chordDist * sectLen) / (np.pi *rdim ** 2)
    solDist[0] = 0

    bladeArea = np.sum(Nb * chordDist * np.insert(np.diff(rdim), 0, np.diff(rdim)[0]))
    sol = bladeArea/(np.pi*R**2)
    #   Coordinates of the lifting line
    liftLineCoord = LENodes - (LENodes - TENodes) * loadPos
    liftLineNorm = np.transpose(np.array((np.sin(twistDist),np.zeros(len(twistDist)),np.cos(twistDist))))

    geomParams = {'liftLineCoord':liftLineCoord,'liftLineNorm':liftLineNorm,'R':R,'e':e,'diskArea':A,'bladeArea':bladeArea,'sectLen':sectLen,'chordDist':chordDist,'twistDist':twistDist,'solDist':solDist,
                  'solidity':sol,'surfNodes':surfNodes,'surfNorms':ScaledNodeCenteredSurfNorms,'nXsecs':nXsecs,'pntsPerXsec':pntsPerXsec,'rdim':rdim,'r':r}


    return geomParams


# n = 49
# fig = plt.figure()
# ax = fig.gca(projection = '3d')
# ax.auto_scale_xyz([-2, 2], [10, 60], [-1, 1])
# ax.pbaspect = [.09, 1, .05]
# ax.set(xlabel = 'x',ylabel = 'y',zlabel = 'z')
# ax.scatter3D(surfNodes[:,0], surfNodes[:,1], surfNodes[:,2],c = 'red',linewidths = 1)
# ax.scatter3D(surfNodes[24,0], surfNodes[24,1], surfNodes[24,2],c = 'yellow',linewidths = 5)
# plt.show()
# ax.quiver(surfNodes[n*pntsPerXsec:n*pntsPerXsec+pntsPerXsec,0], surfNodes[n*pntsPerXsec:n*pntsPerXsec+pntsPerXsec,1], surfNodes[n*pntsPerXsec:n*pntsPerXsec+pntsPerXsec,2],ScaledNodeCenteredSurfNorms[n*pntsPerXsec:n*pntsPerXsec+pntsPerXsec,0]*0.05,ScaledNodeCenteredSurfNorms[n*pntsPerXsec:n*pntsPerXsec+pntsPerXsec,1]*0.05,ScaledNodeCenteredSurfNorms[n*pntsPerXsec:n*pntsPerXsec+pntsPerXsec,2]*0.05)
