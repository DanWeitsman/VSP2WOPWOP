#       VSP2WOPWOP Blade Loading Analysis

#   Author: Daniel Weitsman

#   This function trims the rotor to the desired thrust condition, which is specified in the input file, and
#   computes the aerodynamic loads using BEMT. These quantities are then assembled into a dictionary, which is returned to the user.

#%%
import numpy as np

# %%
def loadingHover(UserIn, geomParams, XsecPolar, T, omega):

    #   Unpacks needed variables from dictionaries
    # T = UserIn['T'][ii]
    # omega = UserIn['omega']
    Nb = UserIn['Nb']
    R = geomParams['R']
    e = geomParams['e']
    chordDist = geomParams['chordDist']
    twistDist = geomParams['twistDist']
    solDist = geomParams['solDist']
    XsecLocation = UserIn['XsecLocation']
    airfoilName = list(XsecPolar.keys())
    rho = UserIn['rho']
    tipLoss = UserIn['tipLoss']
    rdim = geomParams['rdim']
    r = geomParams['r']
    Adisk = geomParams['diskArea']
    sol = geomParams['solidity']

    #%%
    twist = twistDist
    omega = omega / 60 * 2 * np.pi

    TargetCT = T/(rho * Adisk * (omega * R) ** 2)

    Th0 = UserIn['thetaInit'] * (np.pi / 180)
    thetaDist = Th0 + twist
    InitLam = np.full(len(chordDist), np.sqrt(TargetCT / 2))

    if len(XsecLocation) > 1:
        PolarInd = []
        for ii in range(1, len(XsecLocation)):
            PolarInd.append(np.squeeze(np.where((XsecLocation[ii - 1] <= np.round(r, 5)) == (XsecLocation[ii] >= np.round(r, 5)))))
        PolarInd.append(np.arange(PolarInd[-1][-1] + 1, len(r)))
    else:
        PolarInd = np.arange(0,len(r))

    a = np.ones((len(r)))
    a0 = np.ones((len(r)))

    if len(XsecLocation) > 1:
        for ii in range(0, len(PolarInd)):
            a[PolarInd[ii]] = XsecPolar[airfoilName[ii]]['Lift Slope']
            a0[PolarInd[ii]] = XsecPolar[airfoilName[ii]]['Alpha0'] * (np.pi / 180)
    else:
        a = a*XsecPolar[airfoilName[0]]['Lift Slope']
        a0 = a0*XsecPolar[airfoilName[0]]['Alpha0'] * (np.pi / 180)

    # %%
    def TipLoss(lambdaInit, ThetaDist):

        if tipLoss == 1:
            iter = 0
            err = np.ones(len(lambdaInit))

            while np.any(err > 0.005):
                # froot = 0.5*Nb*(r/((1 - r)*lam/r))
                f = 0.5 * Nb * ((1 - r) / lambdaInit)
                F = (2 / np.pi) * np.arccos(np.e ** (-f))
                lam = (solDist * a / (16 * F)) * (np.sqrt(1 + 32 * F / (solDist * a) * ThetaDist * r) - 1)
                # lam = np.sqrt(1/(8*F)*CL*r**2)
                err = np.abs((lam - lambdaInit) / lam)
                err[np.where(np.isnan(err) == 1)] = 0
                lambdaInit = lam
                iter = iter + 1

        else:
            lam = (solDist * a / 16) * (np.sqrt(1 + 32 / (solDist * a) * ThetaDist * r) - 1)

        lam[np.where(np.isnan(lam) == 1)] = 0
        lam[0] = lam[1]
        lam[-1] = lam[-2]
        # return (lam, F, err, iter)
        return (lam)

    def EvalCT(lam, th):
        AoA = th - lam / r
        [CL, CD] = PolarLookup(AoA)
        dCT = 0.5 * solDist * CL * r ** 2
        CT = np.trapz(dCT, r)

        return CT, dCT, CL, CD, AoA

    def PolarLookup(alpha):
        CL = np.zeros(len(alpha))
        CD = np.zeros(len(alpha))

        for i, n in enumerate(alpha):
            if len(XsecLocation) > 1:
                for iii, p in enumerate(PolarInd):
                    if np.any(i == p):
                        polar = (XsecPolar[airfoilName[iii]]['Polar'])
                        leftInd = np.squeeze(np.where(polar[:, 0] < n))[-1]
                        rightInd = np.squeeze(np.where(polar[:, 0] > n))[0]

                        CL[i] = polar[leftInd, 1] + (n - polar[leftInd, 0]) * (polar[rightInd, 1] - polar[leftInd, 1]) / \
                                (polar[rightInd, 0] - polar[leftInd, 0])
                        CD[i] = polar[leftInd, 2] + (n - polar[leftInd, 0]) * (polar[rightInd, 2] - polar[leftInd, 2]) / \
                                (polar[rightInd, 0] - polar[leftInd, 0])
            else:
                polar = (XsecPolar[airfoilName[0]]['Polar'])

                leftInd = np.squeeze(np.where(polar[:, 0] < n))[-1]
                rightInd = np.squeeze(np.where(polar[:, 0] > n))[0]

                CL[i] = polar[leftInd, 1] + (n - polar[leftInd, 0]) * (polar[rightInd, 1] - polar[leftInd, 1]) / \
                        (polar[rightInd, 0] - polar[leftInd, 0])
                CD[i] = polar[leftInd, 2] + (n - polar[leftInd, 0]) * (polar[rightInd, 2] - polar[leftInd, 2]) / \
                        (polar[rightInd, 0] - polar[leftInd, 0])

        return (CL, CD)


    # %%

    errCT = 1
    ii = 0
    delta = 0.0001

    while errCT > 0.0005:

        delTh0 = Th0 + delta * Th0
        delThetaDist = delTh0 + twist

        #[Lambda, F, Err, i]  = TipLoss(InitLam, thetaDist)
        Lambda= TipLoss(InitLam, thetaDist)
        # [delLambda, delF, delErr, deli] = TipLoss(InitLam, delThetaDist)
        delLambda = TipLoss(InitLam, delThetaDist)

        [CT, dCT, dCL, dCD, AoA] = EvalCT(Lambda, thetaDist)
        [delCT, deldCT, delCL, delCD, delAoA] = EvalCT(delLambda, delThetaDist)

        J = (delCT - CT) / (delta * Th0)
        errNR = TargetCT - CT
        errCT = abs(errNR) / TargetCT

        Th0 = delTh0 + 0.5 * J ** -1 * errNR
        thetaDist = Th0 + twist
        InitLam = delLambda

        ii = ii + 1

# %%
    CL = np.trapz(dCL,r)
    CD = np.trapz(dCD,r)

    dCP = 0.5 * solDist * (Lambda / r * dCL + dCD) * r ** 3
    CP = np.trapz(dCP,r)
    P = CP* rho * Adisk * (omega * R) ** 3

    dT = dCT * rho * Adisk * (omega * R) ** 2
    T = np.trapz(dT,r)
    dFz = dT / Nb

    dQ = dCP * rho * Adisk * (omega * R) ** 2 * R
    Q = np.trapz(dQ,r)
    dFx = dQ / (Nb * r * R)

    FM = CP/(1.15*CP+sol/8*CD)


    dFx[np.where(np.isnan(dFx) == 1)] = 0
    dFz[np.where(np.isnan(dFz) == 1)] = 0
    dFy = np.squeeze(np.zeros((len(r))))

    #%%
    dL = 0.5*rho*(omega * R) ** 2*chordDist*dCL
    dD = 0.5 * rho * (omega * R) ** 2 * chordDist * dCD
    dFz2 = dL*np.cos(Lambda / r)-dD*np.sin(Lambda / r)
#todo resolve lift coefficient from blade normal to allign with the TPP
    #%%

    loads = [dFx, dFy, dFz]
    # compactArea = np.insert(np.diff(rdim),0,np.diff(rdim)[0])
    compactArea = chordDist

    loadParams = {'th':[Th0,0,0],'beta':[0,0,0],'CT':CT,'T':T, 'dCT':dCT,'dT':dT, 'CP':CP,'P':P,'Q':Q, 'dCP':dCP,
                  'dQ':dQ,'dCL':dCL, 'dCD':dCD, 'CL':CL,'CD':CD,'FM':FM, 'AoA':AoA,'Lambda':Lambda,'dFx':dFx,'dFy':dFy,'dFz':dFz,'dFz2':dFz2,'compactArea':compactArea,'omega':omega}
    return loadParams







# # %%
# # figdir = os.path.abspath(os.path.join(input.dirDataFile,'Figures/CL.png'))
# # with cbook.get_sample_data(figdir) as image_file:
#
# # fig, ax = plt.subplots()
# # image = plt.imread(figdir)
# # im = ax.imshow(image)
#
#
# fig = plt.figure(figsize=[6.4, 4.5], )
# ax = fig.gca()
# plt.plot(r, CL)
# ax.set_ylabel('Lift Coefficient')
# ax.set_xlabel('Nondimensional radial position, r/R')
# ax.set_title('CT/$\sigma$=0.01')
# plt.grid()
#
# fig = plt.figure(figsize=[6.4, 4.5], )
# ax = fig.gca()
# plt.plot(r, CD)
# ax.set_ylabel('Drag Coefficient')
# ax.set_xlabel('Nondimensional radial position, r/R')
# ax.set_title('CT/$\sigma$=0.01')
# plt.grid()
# # plt.axes([0.25 ,1, 0.4 ,0.9])
# # [lam, F, err, i] = TipLoss(lamInit,thInit,r)
# # lam[np.where(np.isnan(lam) == 1)] = 0
# # AoA = th_init-lam/r
# # [Cl,Cd]=PolarLookup(AoA)
# # dCT = 0.5*sig*Cl*r**2*dr
# # CT_temp = np.trapz(dCT)
# #
#
# # errCT= abs((CT_temp-CT)/CT_temp)
# # th_init = 6*CT_temp/(np.mean(chord)*a)+3/2*np.sqrt(CT_temp/2)
# # ii = ii+1
