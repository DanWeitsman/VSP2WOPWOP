#       VSP2WOPWOP Blade Loading Analysis for Hover and Axial Flight

#   Author: Daniel Weitsman

# This function trims the rotor to the desired thrust condition, which is specified in the input file, and computes
# the aerodynamic loads using BEMT. These quantities are then assembled into a dictionary, which is returned to the
# user.

# %%
import numpy as np
from scipy.optimize import least_squares,minimize


# %%
def loadingAxialHover(UserIn, geomParams, XsecPolar, T, omega, Vz):

    def TipLoss(lambdaInit, ThetaDist):
        """
        This function applies the fixed point itteration method to compute the inflow distribution and applies
        Prandtl's tip loss formulation, if specified for in the input module
        :param lambdaInit: Initial guess for the inflow ratio
        :param ThetaDist: Blade twist distribution (rad)
        :return:
        :param:  lam: radial inflow distribution
        """
        if tipLoss == 1:

            iter = 0
            err = np.ones(len(lambdaInit))

            while np.any(err > 0.005):
                # froot = 0.5*Nb*(r/((1 - r)*lam/r))
                f = 0.5 * Nb * ((1 - r) / lambdaInit)
                F = (2 / np.pi) * np.arccos(np.e ** (-f))
                #lam = (solDist * a / (16 * F)) * (np.sqrt(1 + 32 * F / (solDist * a) * ThetaDist * r) - 1)
                lam = np.sqrt(1/4*(solDist*a/(8*F)-lam_c)**2+solDist*a*ThetaDist*r/(8*F))-(solDist*a/(16*F)-lam_c/2)
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
        return lam

    def EvalCT(lamInit, th):
        """
        This function computes the radial loading distribution based on the thrust coefficient
        :param lam: Initial guess for the radial inflow distribution
        :param th: Radial twist distribution
        :return:
        :param CT: Radially integrated thrust coefficient
        :param dCT: Incremental thrust coefficient
        :param dCL:  Radial distribution of the lift coefficient
        :param dCD:  Radial distribution of the drag coefficient
        :param AoA: Radial angle of attack distribution
        """

        lam = TipLoss(lamInit, th)
        AoA = th - lam / r
        [dCL, dCD] = PolarLookup(AoA)
        dCT = 0.5 * solDist * dCL * r ** 2
        CT = np.trapz(dCT, r)

        return CT, dCT, dCL, dCD, lam,AoA

    def PolarLookup(alpha):
        """
        This function looks up the sectional blade loads from the XFoil polar based on the computed angle of attack
        distribution
        :param alpha: Angle of attack distribution from EvalCT
        :return:
        :param dCL:  Radial distribution of the lift coefficient
        :param dCD:  Radial distribution of the drag coefficient
        """
        dCL = np.zeros(len(alpha))
        dCD = np.zeros(len(alpha))

        for i, n in enumerate(alpha):
            if len(XsecLocation) > 1:
                for iii, p in enumerate(PolarInd):
                    if np.any(i == p):
                        polar = (XsecPolar[airfoilName[iii]]['Polar'])
                        leftInd = np.squeeze(np.where(polar[:, 0] < n))[-1]
                        rightInd = np.squeeze(np.where(polar[:, 0] > n))[0]

                        dCL[i] = polar[leftInd, 1] + (n - polar[leftInd, 0]) * (
                                    polar[rightInd, 1] - polar[leftInd, 1]) / \
                                 (polar[rightInd, 0] - polar[leftInd, 0])
                        dCD[i] = polar[leftInd, 2] + (n - polar[leftInd, 0]) * (
                                    polar[rightInd, 2] - polar[leftInd, 2]) / \
                                 (polar[rightInd, 0] - polar[leftInd, 0])
            else:
                polar = (XsecPolar[airfoilName[0]]['Polar'])

                leftInd = np.squeeze(np.where(polar[:, 0] < n))[-1]
                rightInd = np.squeeze(np.where(polar[:, 0] > n))[0]

                dCL[i] = polar[leftInd, 1] + (n - polar[leftInd, 0]) * (polar[rightInd, 1] - polar[leftInd, 1]) / \
                         (polar[rightInd, 0] - polar[leftInd, 0])
                dCD[i] = polar[leftInd, 2] + (n - polar[leftInd, 0]) * (polar[rightInd, 2] - polar[leftInd, 2]) / \
                         (polar[rightInd, 0] - polar[leftInd, 0])

        return dCL, dCD

    def residuals( Th0, twistDist, targCT):
        '''
        This function computes the residuals from the target trim variables.
        :param lamInit: initial estimate for the inflow distribution
        :param Th0: collective pitch setting
        :param twistDist: radial twist distribution
        :param targCT: target thrust coefficient
        :return:
        :param res: residual in the percentage error between the target and computed thrust coefficient
        '''

        th = Th0 + twistDist
        trim_out = EvalCT(lamInit, th)
        trimCT = trim_out[0]
        res = np.abs((targCT - trimCT) / targCT)
        #res = targCT - trimCT

        return res

    # %%
    #   This block of code defines parameters that are used throughout the remainder of the module
    Nb = UserIn['Nb']
    R = geomParams['R']
    chordDist = geomParams['chordDist']
    twistDist = geomParams['twistDist']
    solDist = geomParams['solDist']
    XsecLocation = UserIn['XsecLocation']
    airfoilName = list(XsecPolar.keys())
    rho = UserIn['rho']
    tipLoss = UserIn['tipLoss']
    r = geomParams['r']
    Adisk = geomParams['diskArea']
    sol = geomParams['solidity']

    # %%
    #   converts rotational rate from degrees to radians
    omega = omega / 60 * 2 * np.pi
    #   Target thrust coefficient
    targCT = T / (rho * Adisk * (omega * R) ** 2)
    #   Converts initial guess for the collective pitch setting from degrees to radians
    Th0 = UserIn['thetaInit'] * (np.pi / 180)
    #   Initial guess for the radial inflow distribution
    lamInit = np.full(len(chordDist), np.sqrt(targCT / 2))
    #   Axial climb/descent inflow ratio
    lam_c = Vz / (omega * R)

    assert -2 <= Vz/np.sqrt(T/(2*rho*Adisk)) <= 0,'Non-physical solution, 1D assumption of momentum theory is violated'
        # raise ValueError('Non-physical solution, 1D assumption of momentum theory is violated')

    #   Populates a list with indices that correspond to each airfoil cross-section.
    if len(XsecLocation) > 1:
        PolarInd = []
        for ii in range(1, len(XsecLocation)):
            PolarInd.append(
                np.squeeze(np.where((XsecLocation[ii - 1] <= np.round(r, 5)) == (XsecLocation[ii] >= np.round(r, 5)))))
        PolarInd.append(np.arange(PolarInd[-1][-1] + 1, len(r)))
    else:
        PolarInd = np.arange(0, len(r))

    #   initializes and populates an array with the sectional lift curve slope.
    a = np.ones((len(r)))
    a0 = np.ones((len(r)))

    if len(XsecLocation) > 1:
        for ii in range(0, len(PolarInd)):
            a[PolarInd[ii]] = XsecPolar[airfoilName[ii]]['Lift Slope']
            a0[PolarInd[ii]] = XsecPolar[airfoilName[ii]]['Alpha0'] * (np.pi / 180)
    else:
        a = a * XsecPolar[airfoilName[0]]['Lift Slope']
        a0 = a0 * XsecPolar[airfoilName[0]]['Alpha0'] * (np.pi / 180)

    # %%

    # This function employs the non-linear least square optimization method (LM) to compute the necessary collective
    # pitch settings to minimize the residuals, which is the percentage error between the target and computed thrust
    # coefficient.
    trim_sol = least_squares(residuals, Th0, args=[ twistDist, targCT], method='lm')

    # Evaluates the core function, 'EvalCT' once more with the trimmed collective pitch to attain the remainder of
    # load parameters needed to compute the integrated and pre-integraated blade load.
    CT, dCT, dCL, dCD, lam, AoA = EvalCT(lamInit, trim_sol.x+twistDist)

    # Assembles the pitch setting into a list, which is referenced in the 'nml_write' module
    th = list(trim_sol.x)+[0,0]

    # errCT = 1
    # ii = 0
    # delta = 0.0001
    #
    # while errCT > 0.0005:
    #
    #     thetaDist = Th0+twistDist
    #     delTh0 = Th0 + delta * Th0
    #     delThetaDist = delTh0 + twistDist
    #
    #     [CT, dCT, dCL, dCD, lam,AoA] = EvalCT(lamInit, thetaDist)
    #     [delCT, deldCT, delCL, delCD,dellam, delAoA] = EvalCT(lamInit, delThetaDist)
    #
    #     J = (delCT - CT) / (delta * Th0)
    #     errNR = targCT - CT
    #     errCT = abs(errNR) / targCT
    #
    #     Th0 = delTh0 + 0.5 * J ** -1 * errNR
    #     thetaDist = Th0 + twistDist
    #     lamInit = dellam
    #
    #     ii = ii + 1

    # %%

    #   Integrated lift and drag coefficients
    CL = np.trapz(dCL, r)
    CD = np.trapz(dCD, r)

    #   Distribution and integrated of the power/torque coefficient
    dCP = 0.5 * solDist * (lam / r * dCL + dCD) * r ** 3
    CP = np.trapz(dCP, r)

    #   Power required by the rotor
    P = CP * rho * Adisk * (omega * R) ** 3

    #   Distribution and integrated thrust
    dT = dCT * rho * Adisk * (omega * R) ** 2
    T = np.trapz(dT, r)

    #   Resolved force component that is normal to the rotor plane
    dFz = dT / Nb

    #   Distribution and integrated torque
    dQ = dCP * rho * Adisk * (omega * R) ** 2 * R
    Q = np.trapz(dQ, r)

    #   Resolved in-plane force component
    dFx = dQ / (Nb * r * R)

    #   Figure of merit, induced power factor = 1.15
    FM = CP / (1.15 * CP + sol / 8 * CD)

    #   Sets any infinite values of the computed force components (primarily near the blade root) equal to zero.
    dFx[np.where(np.isnan(dFx) == 1)] = 0
    dFz[np.where(np.isnan(dFz) == 1)] = 0
    dFy = np.squeeze(np.zeros((len(r))))

    #   Assembles loads into a list that is referenced in the 'ConstantLoadingPatchFileWrite' module
    loads = [dFx, dFy, dFz]

#%%
    # Assembles all computed load parameters into a dictionary
    loadParams = {'residuals':trim_sol.fun,'th': th, 'beta': [0, 0, 0], 'CT': CT, 'T': T, 'dCT': dCT, 'dT': dT, 'CP': CP, 'P': P,
                  'Q': Q, 'dCP': dCP, 'dQ': dQ, 'dCL': dCL, 'dCD': dCD, 'CL': CL, 'CD': CD, 'FM': FM, 'AoA': AoA, 'lambda': lam,
                  'dFx': dFx, 'dFy': dFy, 'dFz': dFz, 'omega': omega}
    return loadParams

# todo implement axial flight and ensure the inflow distribution does not change. Then compare the results attained
#  from the hover and FF modules and determine which one to use for axial flight conditions.
#
# %% # figdir = os.path.abspath(os.path.join(input.dirDataFile,'Figures/CL.png')) # with cbook.get_sample_data(figdir) as
#  image_file:
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
