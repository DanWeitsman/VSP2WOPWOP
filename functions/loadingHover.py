#       VSP2WOPWOP Blade Loading Analysis for Hover and Axial Flight

#   Author: Daniel Weitsman

# This function trims the rotor to the desired thrust condition, which is specified in the input file, and computes
# the aerodynamic loads using BEMT. These quantities are then assembled into a dictionary, which is returned to the
# user.

# %%
import numpy as np
from scipy.optimize import least_squares,minimize
import bisect

# %%
def loadingAxialHover(UserIn, geomParams, XsecPolar, T, omega, Vz):

    def rpm_residuals(omega):
        '''
        This function computes the residuals based on the percentage difference between the computed and target thrust.
        The rotational rate is adjusted until this residual is minimized.
        :param omega: rotational rate [rad/s]
        :return:
        :param res:  percentage error between the target and computed T
        '''

        trim_out = rpm_trim(omega)
        res = np.abs((T - trim_out[0]*rho*np.pi*R**2*(omega*R)**2) / T)
        return res

    def coll_residuals(th0):
        '''
        This function computes the residuals based on the percentage difference between the computed and target CT.
        The collective pitch is adjusted until this residual is minimized.
        :param th0: collective pitch setting
        :return:
        :param res:  percentage error between the target and computed CT
        '''
        th = th0+twistDist
        trim_out = coll_trim(th)
        # res = targCT - trim_out[0]
        res = np.abs((targCT - trim_out[0]) / targCT)
        print(res)
        return res

    def rpm_trim(omega):
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
        CT_init = T / (rho * np.pi * R ** 2 * (omega * R) ** 2)
        lam_init = np.sqrt(CT_init / 2)
        err = 1
        while np.any(err > 0.0005):
            lam = TipLoss(lam_init, twistDist)
            AoA = twistDist - lam / r
            dCL, dCD = PolarLookup(AoA)
            dCT = 0.5 * solDist * dCL * r ** 2
            CT = np.trapz(dCT, r)
            err = np.abs((lam - lam_init) / lam)
            lam_init = lam

        return CT, dCT, dCL, dCD, lam, AoA

    def coll_trim(th):
        """
        This function computes the radial loading distribution based on the thrust coefficient
        :param th: collective pitch setting+twist distribution [rad]
        :return:
        :param CT:  integrated thrust coefficient
        :param dCT:  radial thrust coefficient distribution
        :param dCL:  radial lift coefficient distribution
        :param dCD:  radial drag coefficient distribution
        :param lam:  radial inflow coefficient distribution
        :param AoA: radial angle of attack distribution [rad]
        """
        lam = TipLoss(lamInit, th)
        AoA = th - lam / r
        dCL, dCD = PolarLookup(AoA)
        dCT = 0.5 * solDist * (dCL*np.cos(lam/r)-dCD*np.sin(lam/r))* r ** 2
        CT = np.trapz(dCT, r)

        return CT, dCT, dCL, dCD, lam, AoA

    def TipLoss(lambdaInit, ThetaDist):
        """
        This function applies the fixed point iteration method to compute the inflow distribution and applies
        Prandtl's tip loss formulation, if specified for in the input module
        :param lambdaInit: Initial guess for the inflow ratio
        :param ThetaDist: collective pitch angle + twist distribution (rad)
        :return:
        :param:  lam: radial inflow distribution
        """
        if tipLoss == 1:
            iter = 0
            err = np.ones(len(r))
            while np.any(err > 0.005):
                # froot = 0.5*Nb*(r/((1 - r)*lam/r))
                f = 0.5 * Nb * ((1 - r) / lambdaInit)
                F = (2 / np.pi) * np.arccos(np.e ** (-f))
                # lam = np.sqrt((solDist*XsecPolarExp['Lift Slope']/16-lam_c/2)**2+solDist*XsecPolarExp['Lift Slope']/8*ThetaDist*r)-(solDist*XsecPolarExp['Lift Slope']/16-lam_c/2)
                # lam = (solDist * XsecPolarExp['Lift Slope']/ (16 * F)) * (np.sqrt(1 + 32 * F / (solDist *XsecPolarExp['Lift Slope']) * ThetaDist * r) - 1)
                lam = np.sqrt(1/4*(solDist*XsecPolarExp['Lift Slope']/(8*F)-lam_c)**2+solDist*XsecPolarExp['Lift Slope']*ThetaDist*r/(8*F))-(solDist*XsecPolarExp['Lift Slope']/(16*F)-lam_c/2)
                err = np.abs((lam - lambdaInit) / lam)
                err[np.where(np.isnan(err) == 1)] = 0
                lambdaInit = lam
                iter = iter + 1
        else:
            lam = (solDist * XsecPolarExp['Lift Slope'] / 16) * (np.sqrt(1 + 32 / (solDist * XsecPolarExp['Lift Slope']) * ThetaDist * r) - 1)

        lam[np.where(np.isnan(lam) == 1)] = 0
        # lam[0] = lam[1]
        # lam[-1] = lam[-2]
        return lam

    def PolarLookup(AoA):
        """
        This function linearly interpolates the sectional blade load coefficients from the XFoil polar based on the
        computed angle of attack distribution. If the blade section is stalled CL at that section is linearly
        interpolated between the maximum and minimum CL, while CD is simply set to its maximum value for the
        respective airfoil.
        :param alpha: angle of attack distribution
        return:
        :param dCL:  radial lift coefficient distribution
        :param dCD:  radial drag coefficient distribution
        """
        dCL = np.zeros(len(AoA))
        dCD = np.zeros(len(AoA))
        for i, alpha in enumerate(AoA):
            if alpha > XsecPolar[polarInd[i]]['alphaMax']:
                dCL[i] = np.interp(alpha ,xp = [XsecPolar[polarInd[i]]['alphaMax'],XsecPolar[polarInd[i]]['Alpha0']%(2*np.pi)],fp = [XsecPolar[polarInd[i]]['ClMax'],XsecPolar[polarInd[i]]['ClMin']])
                dCD[i] = XsecPolar[polarInd[i]]['CdMax']
            else:
                 AoA_ind = bisect.bisect(XsecPolar[polarInd[i]]['Polar'][:, 0], alpha)
                 dCL[i] = np.interp(alpha ,xp = [XsecPolar[polarInd[i]]['Polar'][AoA_ind-1,0],XsecPolar[polarInd[i]]['Polar'][AoA_ind,0]],fp = [XsecPolar[polarInd[i]]['Polar'][AoA_ind-1,1],XsecPolar[polarInd[i]]['Polar'][AoA_ind,1]])
                 dCD[i] = np.interp(alpha ,xp = [XsecPolar[polarInd[i]]['Polar'][AoA_ind-1,0],XsecPolar[polarInd[i]]['Polar'][AoA_ind,0]],fp = [XsecPolar[polarInd[i]]['Polar'][AoA_ind-1,2],XsecPolar[polarInd[i]]['Polar'][AoA_ind,2]])
        return dCL, dCD

    # %%
    #   This block of code defines parameters that are used throughout the remainder of the module
    Nb = UserIn['Nb']
    R = geomParams['R']
    chordDist = geomParams['chordDist']
    twistDist = geomParams['twistDist']
    solDist = geomParams['solDist']
    XsecLocation = UserIn['XsecLocation']
    rho = UserIn['rho']
    tipLoss = UserIn['tipLoss']
    r = geomParams['r']
    Adisk = geomParams['diskArea']
    sol = geomParams['solidity']
    #   converts rotational rate from degrees to radians
    omega = omega / 60 * 2 * np.pi
    #   Sectional free stream velocity
    UP = omega*geomParams['rdim']
    #   Target thrust coefficient
    targCT = T / (rho * Adisk * (omega * R) ** 2)
    #   Converts initial guess for the collective pitch setting from degrees to radians
    th0 = UserIn['thetaInit'] * (np.pi / 180)
    #   Initial guess for the radial inflow distribution
    lamInit = np.ones(len(r))*np.sqrt(targCT / 2)
    #   Axial climb/descent inflow ratio
    lam_c = Vz / (omega * R)

    if -2 < Vz/np.sqrt(T/(2*rho*Adisk)) < 0:
         raise ValueError('Non-physical solution, 1D assumption of momentum theory is violated')

#%% This section of the code populates an array of the airfoil names based on their radial location along the blade span

# initializes the expanded Xsect polar dictionary, which will store all the airfoil parameters corresponding to their
    # radial location
    XsecPolarExp = {}
    polarInd = []
    #   if multiple airfoil sections are used along the blade span are used this section of the code would be executed
    if len(XsecLocation) > 1:
        ind = np.zeros((len(XsecLocation) + 1))
    #   creates an array of size r that is filled with the indices corresponding to the radial location of each airfoil
        for i, Xsec in enumerate(XsecLocation):
            ind[i] = bisect.bisect(r, Xsec)
        ind[0] = 0
        ind[-1] = len(r)
    # loops through each airfoil section and their parameters, populating an array of size r, with these parameters.
        # These arrays are then written to the XsecPolarExp dictionary.
        for i, Xsec in enumerate(XsecPolar.keys()):
            polarInd.extend([Xsec] * int(ind[i + 1] - ind[i]))
            for ii, param in enumerate(list(XsecPolar[Xsec].keys())[1:]):
                if i == 0:
                    XsecPolarExp = {**XsecPolarExp, **{param:XsecPolar[Xsec][param]*np.ones(len(r))}}
                else:
                    XsecPolarExp[param][int(ind[i]):] = XsecPolar[Xsec][param]
        # if only a single airfoil section is used along the blade span the section's parameters are looped over,
        # expanded to correspond to each blade section, and assembled into the XsecPolarExp dictionary.
    else:
        polarInd = list(XsecPolar.keys())*len(r)
        for i,key in enumerate(list(XsecPolar[list(XsecPolar.keys())[0]].keys())[1:]):
            XsecPolarExp[key] = np.ones(len(r))*XsecPolar[list(XsecPolar.keys())[0]][key]
    XsecPolarExp['Lift Slope'] = XsecPolarExp['Lift Slope']*np.pi/180

    # %%

    # This function employs the non-linear least square optimization method (LM) to compute the necessary rotational rate or collective
    # pitch angle to meet the target thrust or thrust coefficient, respectively.
    if UserIn['trim'] == 1:
        trim_sol = least_squares(rpm_residuals, omega, method='lm')
        CT, dCT, dCL, dCD, lam, AoA = rpm_trim(trim_sol.x)
        omega = trim_sol.x
    else:
        trim_sol = least_squares(coll_residuals, th0, method='lm')
        CT, dCT, dCL, dCD, lam, AoA = coll_trim(trim_sol.x+twistDist)
        th = np.array([np.squeeze(trim_sol.x), 0, 0])

#%%
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
    dFy = np.zeros(len(r))


#%%
    # Assembles all computed load parameters into a dictionary
    loadParams = {'coll_residuals':trim_sol.fun,'th': th, 'beta': [0, 0, 0], 'CT': CT, 'T': T, 'dCT': dCT, 'dT': dT, 'CP': CP, 'P': P,
                  'Q': Q, 'dCP': dCP, 'dQ': dQ, 'dCL': dCL, 'dCD': dCD, 'CL': CL, 'CD': CD, 'FM': FM, 'AoA': AoA,'ClaDist':XsecPolarExp['Lift Slope'], 'lambda': lam,
                  'dFx': dFx, 'dFy': dFy, 'dFz': dFz, 'omega': omega,'UP':UP}
    return loadParams

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
