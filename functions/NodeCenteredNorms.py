#       VSP2WOPWOP Scaled Node-Centered Normal Vectors

#   Author: Daniel Weitsman

#   This function computes the face-area scaled node-centered normal vectors needed for the patch files.
#   Even though the DegenGeom file contains the face centered normals, the version of PSU-WOPWOP that I was running was
#   throwing an error when I used the face-centered normals, so I wrote this function to bypass the glitch.

#%%
import numpy as np
#%%
#   This function computes the node centered normal vectors
#   Input:
#       surfNodes: Coordinates of the surface nodes
#       pntsPerXsec: Number of chordwise nodes [i]
#       nXsecs: Number of spanwise nodes [j]
#   Output:
#       NormNodeCenteredNorms: Node centered surface normals, which are normalized by the corresponding face areas

    #%%
def NodeCenteredNorms(surfNodes,pntsPerXsec,nXsecs):

    # Parses and sorts nodes to correspond to the i and j direction
    sortedNodes = np.zeros((pntsPerXsec, nXsecs, 3))
    for i in range(0, nXsecs):
        sortedNodes[:, i, :] = surfNodes[pntsPerXsec * i:pntsPerXsec * i + pntsPerXsec]

    #%%     Separates out nodes based on blade features

    TENodes = sortedNodes[0, :, :]
    LENodes = sortedNodes[int((pntsPerXsec - 1) / 2), :, :]

    RootNodes = sortedNodes[:, 0, :]
    TipNodes = sortedNodes[:, -1, :]

    BottomSurf = sortedNodes[1:int((pntsPerXsec - 1) / 2), 1:-1, :]
    TopSurf = sortedNodes[int((pntsPerXsec - 1) / 2):pntsPerXsec, 1:-1, :]

    #%%    Computes difference between adjacent chordwise and spanwise nodes needed for calculating the cross product
    # chordwise forward difference
    ChordSortedNodesFwdDiff = np.diff(sortedNodes, axis=0)
    # spanwise forward difference
    SpanSortedNodesFwdDiff = np.diff(sortedNodes, axis=1)
    # chordwise backwards difference
    ChordSortedNodesBackDiff = np.diff(sortedNodes[::-1, :, :], axis=0)[::-1, :, :]
    # spanwise backwards difference
    SpanSortedNodesBackDiff = np.diff(sortedNodes[:, ::-1, :], axis=1)[:, ::-1, :]

    #%%     Mid-section nodes cross product computation (excluding root, tip, trailing edge upper, and trailing edge lower nodes (pntsPerXsec-2 x nXsecs-2))
    # computes the cross product between every mid-section node and its four adjecent nodes
    a = np.cross(ChordSortedNodesFwdDiff[1:,1:-1,:],SpanSortedNodesFwdDiff[1:-1,1:,:])
    b = np.cross(SpanSortedNodesFwdDiff[1:-1,1:,:],ChordSortedNodesBackDiff[:-1,1:-1,:])
    c = np.cross(ChordSortedNodesBackDiff[:-1,1:-1,:],SpanSortedNodesBackDiff[1:-1,:-1,:])
    d = np.cross(SpanSortedNodesBackDiff[1:-1,:-1,:],ChordSortedNodesFwdDiff[1:,1:-1,:])
    # average of the four cross product calculations
    midSectCrossp = np.mean([a, b, c, d], axis=0)
    # average of the magnitudes of the cross product calculations to attain face areas
    midSectArea = np.mean([np.linalg.norm(a, axis=2), np.linalg.norm(b, axis=2), np.linalg.norm(c, axis=2), np.linalg.norm(d, axis=2)],axis=0)

    #%% TE cross product calculation (1 x nXsecs-2))
    #  TE upper surface normal calculation
    aTE = np.cross(SpanSortedNodesFwdDiff[0,1:,:],ChordSortedNodesBackDiff[-1,:,:][1:-1,:])
    bTE = np.cross(ChordSortedNodesBackDiff[-1,:,:][1:-1,:],SpanSortedNodesBackDiff[0,:-1,:])
    # average of the two cross product calculations
    UpperTECrossp = np.mean([aTE,bTE],axis = 0)
    # average of the magnitudes of the two cross product calculation to attain face areas
    upperTEArea = np.mean([np.linalg.norm(aTE,axis = 1),np.linalg.norm(bTE,axis = 1)],axis =0)

    #  TE lower surface normal calculation
    cTE = np.cross(ChordSortedNodesFwdDiff[0,:,:][1:-1,:],SpanSortedNodesFwdDiff[0,1:,:])
    dTE = np.cross(SpanSortedNodesBackDiff[0,:-1,:],ChordSortedNodesFwdDiff[0,:,:][1:-1,:])
    # average of the two cross product calculations
    LowerTECrossp = np.mean([cTE,dTE],axis = 0)
    # average of the magnitudes of the two cross product calculation to attain face areas
    lowerTEArea = np.mean([np.linalg.norm(cTE,axis = 1),np.linalg.norm(dTE,axis = 1)],axis =0)

    #%% Root and tip nodes cross product calculation (pntsPerXsec-2 x 1))

    aRoot = np.cross(ChordSortedNodesFwdDiff[1:, 0, :], SpanSortedNodesFwdDiff[1:-1, 0, :])
    bRoot = np.cross(SpanSortedNodesFwdDiff[1:-1, 0, :], ChordSortedNodesBackDiff[:-1, 0, :])
    RootCrossp = np.mean([aRoot, bRoot], axis=0)
    RootArea = np.mean([np.linalg.norm(aRoot, axis=1), np.linalg.norm(bRoot, axis=1)], axis=0)

    aTip = np.cross(ChordSortedNodesBackDiff[:-1, -1, :], SpanSortedNodesBackDiff[1:-1, -1, :])
    bTip = np.cross(SpanSortedNodesBackDiff[1:-1, -1, :], ChordSortedNodesFwdDiff[1:, -1, :])
    TipCrossp = np.mean([aTip, bTip], axis=0)
    TipArea = np.mean([np.linalg.norm(aTip, axis=1), np.linalg.norm(bTip, axis=1)], axis=0)

    #%% Corner cross product calculation

    RootLower = np.cross(ChordSortedNodesFwdDiff[0, 0, :], SpanSortedNodesFwdDiff[0, 0, :])
    RootUpper = np.cross(SpanSortedNodesFwdDiff[0, 0, :], ChordSortedNodesBackDiff[-1, 0, :])
    TipLower = np.cross(ChordSortedNodesFwdDiff[0, -1, :], SpanSortedNodesFwdDiff[0, -1, :])
    TipUpper = np.cross(SpanSortedNodesFwdDiff[0, -1, :], ChordSortedNodesBackDiff[-1, -1, :])

    RootLowerArea = np.linalg.norm(RootLower)
    RootUpperArea = np.linalg.norm(RootUpper)
    TipLowerArea = np.linalg.norm(TipLower)
    TipUpperArea = np.linalg.norm(TipUpper)

    #%% Assembles all cross products and face areas back into a single sorted array

    rootNorms = np.concatenate((np.expand_dims(RootLower, axis=0), RootCrossp, np.expand_dims(RootUpper, axis=0)))
    tipNorms = np.concatenate((np.expand_dims(TipLower, axis=0), TipCrossp, np.expand_dims(TipUpper, axis=0)))
    midSectNorms = np.concatenate((np.expand_dims(LowerTECrossp, axis=0), midSectCrossp, np.expand_dims(UpperTECrossp, axis=0)))
    sortedNorms = np.concatenate((np.expand_dims(rootNorms, axis=1), midSectNorms, np.expand_dims(tipNorms, axis=1)),axis=1)

    rootArea = np.concatenate((np.expand_dims(RootLowerArea, axis=0), RootArea, np.expand_dims(RootUpperArea, axis=0)))
    tipArea = np.concatenate((np.expand_dims(TipLowerArea, axis=0), TipArea, np.expand_dims(TipUpperArea, axis=0)))
    midSectArea = np.concatenate((np.expand_dims(lowerTEArea, axis=0), midSectArea, np.expand_dims(upperTEArea, axis=0)))
    sortedArea = np.concatenate((np.expand_dims(rootArea, axis=1), midSectArea, np.expand_dims(tipArea, axis=1)),axis=1)

    #%% Normalized normal vectors by face area

    NormsortedNorms = sortedNorms/np.expand_dims(sortedArea,axis = 2)

    #%% Reshapes the assembled array to original input dimensions

    NodeCenteredNorms = np.zeros((pntsPerXsec * nXsecs, 3))
    NormNodeCenteredNorms = np.zeros((pntsPerXsec * nXsecs, 3))
    for i in range(0, nXsecs):
        NodeCenteredNorms[pntsPerXsec * i:pntsPerXsec * i + pntsPerXsec, :] = sortedNodes[:, i, :]
        NormNodeCenteredNorms[pntsPerXsec * i:pntsPerXsec * i + pntsPerXsec, :] = NormsortedNorms[:, i, :]

    return NormNodeCenteredNorms