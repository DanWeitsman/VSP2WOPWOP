#       VSP2WOPWOP Blade Loading Analysis for Hover and Axial Flight

#   Author: Daniel Weitsman

# This function trims the rotor to the desired thrust condition, which is specified in the input file, and computes
# the aerodynamic loads using BEMT. These quantities are then assembled into a dictionary, which is returned to the
# user.


# %%
def loadingHover(UserIn, geomParams, XsecPolar, T, omega, Vz):

    import numpy as np
    from scipy.optimize import least_squares
    import bisect
    from scipy.fft import fft, ifft
    from scipy.interpolate import interp1d

    def rpm_residuals(omega):
        '''
        This function computes the residuals based on the percentage difference between the computed and target thrust.
        The rotational rate is adjusted until this residual is minimized.
        :param omega: rotational rate [rad/s]
        :return:
        :param res:  percentage error between the target and computed T
        '''

        trim_out = rpm_trim(omega)
        res = np.abs((T - trim_out[0]*rho*np.pi*R**2*(omega*R)**2)/T)
        print(res)
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
        th = th0+twistDist
        err = 1
        while np.any(err > 0.0005):
            lam = TipLoss(lam_init, th)
            AoA = th - lam / r
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
        # dCT = 0.5 * solDist * dCL * r ** 2
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
                F[np.where(F == 0)] = 1e-10
                lam = np.sqrt(1/4*(solDist*XsecPolarExp['Lift Slope']/(8*F)-lam_c)**2+solDist*XsecPolarExp['Lift Slope']*ThetaDist*r/(8*F))-(solDist*XsecPolarExp['Lift Slope']/(16*F)-lam_c/2)
                err = np.abs((lam - lambdaInit) / lam)
                err[np.where(np.isnan(err) == 1)] = 0
                lambdaInit = lam
                iter = iter + 1
        else:
            F = 1
            lam = np.sqrt(1/4*(solDist*XsecPolarExp['Lift Slope']/(8*F)-lam_c)**2+solDist*XsecPolarExp['Lift Slope']*ThetaDist*r/(8*F))-(solDist*XsecPolarExp['Lift Slope']/(16*F)-lam_c/2)

        # lam[np.where(np.isnan(lam) == 1)] = 1
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
        # dCL2 = np.zeros(len(AoA))
        # dCD2 = np.zeros(len(AoA))

        for k,v in XsecPolar.items():

            ind = np.where(np.array(polarInd) == k)

            if np.any(AoA[ind]>v['alphaMax']) or np.any(AoA[ind]< v['Alpha0']):
                ind2 = np.squeeze(np.where((AoA[ind]< v['Alpha0']) | (AoA[ind]>v['alphaMax'])))
                dCL[ind2] = np.interp(AoA[ind2]%(2*np.pi) ,xp = [v['alphaMax'],2*np.pi],fp = [v['ClMax'],v['ClMin']])
                dCD[ind2] = v['CdMax']
                ind = np.delete(ind, ind2)

            dCL[ind] = np.interp(AoA[ind],xp = v['Polar'][:,0],fp =v['Polar'][:,1])
            dCD[ind] = np.interp(AoA[ind],xp = v['Polar'][:,0],fp =v['Polar'][:,2])

            # x1_ind = np.squeeze(np.where(v['Polar'][:, 0] > AoA[ind][-3]))[0]
            # x1_ind = [np.squeeze(np.where(v['Polar'][:, 0] > x))[0] for x in AoA[ind]]
            # x0_ind = [np.squeeze(np.where(v['Polar'][:, 0] < x))[-1] for x in AoA[ind]]
            #
            # # x0_ind = np.squeeze(np.where(v['Polar'][:, 0] < AoA[ind][-3]))[-1]
            # dCL[ind] = (v['Polar'][x0_ind,1]*(v['Polar'][x1_ind,0]-AoA[ind])+v['Polar'][x1_ind,1]*(AoA[ind]-v['Polar'][x0_ind,0]))/(v['Polar'][x1_ind,0]-v['Polar'][x0_ind,0])
            # dCD[ind] = (v['Polar'][x0_ind,2]*(v['Polar'][x1_ind,0]-AoA[ind])+v['Polar'][x1_ind,2]*(AoA[ind]-v['Polar'][x0_ind,0]))/(v['Polar'][x1_ind,0]-v['Polar'][x0_ind,0])
            # dCL[ind] = AoA[ind]*2*np.pi
            # dCD[ind] =  0.1*dCL[ind]

        return dCL, dCD

    def fsmooth(quant):
        '''
        Frequency domain smoothing filter that retains only the 1/rev frequency components.
        :param quant: Quantitiy to filter should be of size [phi x r]
        :return:
        '''
        Xm = fft(quant, axis=0)
        Xm[10:-9] = 0
        quant_smooth = np.real(ifft(Xm, axis=0))

        return quant_smooth
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
    th0 = UserIn['thetaInit'] * np.pi / 180
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

    # %%

    # This function employs the non-linear least square optimization method (LM) to compute the necessary rotational rate or collective
    # pitch angle to meet the target thrust or thrust coefficient, respectively.
    if UserIn['trim'] == 1:
        trim_sol = least_squares(rpm_residuals, omega, method='lm')
        CT, dCT, dCL, dCD, lam, AoA = rpm_trim(trim_sol.x)
        omega = trim_sol.x
        th = np.array([th0, 0, 0])
    else:
        trim_sol = least_squares(coll_residuals, th0, method='lm')
        CT, dCT, dCL, dCD, lam, AoA = coll_trim(trim_sol.x+twistDist)
        th = np.array([np.squeeze(trim_sol.x), 0, 0])

#%%

    th1c = 0*np.pi/180
    th1s = 0*np.pi/180
    phiRes = 361
    phi = np.linspace(0,2*np.pi,phiRes)
    theta_expanded = geomParams['twistDist'] + th[0] + np.expand_dims(th1c * np.cos(phi), axis=1) + np.expand_dims(th1s * np.sin(phi), axis=1)
    dCT2, dCL2, dCD2, lam2, AoA2 =   np.array([coll_trim(theta)[1:] for theta in theta_expanded]).transpose((1,0,-1))
    dCL2, dCD2 = [fsmooth(x) for x in [dCL2,dCD2]]
    dCT2 = 0.5 * solDist * (dCL2 * np.cos(lam2 / r) - dCD2 * np.sin(lam2 / r)) * r ** 2
    dT2 = dCT2 * rho * Adisk * (omega * R) ** 2
    dFz2 = dT2 / Nb
    dCP2 = 0.5 * solDist * (lam2 / r * dCL2 + dCD2) * r ** 3
    dQ2 = dCP2 * rho * Adisk * (omega * R) ** 2 * R
    dFx2 = dQ2 / (Nb * r * R)

    if UserIn['rotation'] == 2:
        dFz2 = np.flip(dFz2,axis = 0)
        dFx2 = np.flip(dFx2, axis=0)
    # # up = lam
    # # ut = r
    # # AoA2 = theta_expanded - up/ut
    # # CL,CD = np.array([PolarLookup(alpha) for alpha in AoA2]).transpose((1,0,2))
    # # dT = rho*np.pi*R**2*(omega*R)**2*(1 / 2 * solDist * r ** 2 * (CL * np.cos(up / ut) - CD * np.sin(up / ut)))
    # # dQ = rho*np.pi*R**3*(omega*R)**2*(0.5*solDist*r**3*(CL*np.sin(up/ut)+CD*np.cos(up/ut)))
    # # dFz2 = dT / Nb
    # # dFx2 = dQ/(Nb*r*R)
    #
    # if UserIn['rotation'] == 2:
    #     dFz2 = np.flip(dFz2, axis=0)
    #     dFx2 = np.flip(dFx2, axis=0)

    #%%

    U =np.sqrt((omega*geomParams['rdim'])**2+(omega*R*lam)**2)
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

    #   Distribution and integrated torque
    dQ = dCP * rho * Adisk * (omega * R) ** 2 * R
    Q = np.trapz(dQ, r)

    # Rotates the normal force component by the collective pitch setting, so that a single change of base (CB) can be
    # applied to the blade geometry and loading vector in the namelist file. If the collective pitch CB is
    # unnecessary, then dFz = dT/Nb.
    # dFz = dT/Nb*np.cos(-th[0])-dQ/(Nb*r*R)*np.sin(-th[0])
    dFz = dT/Nb
    # Rotates the inplane force component by the collective pitch setting, so that a single change of base (CB) can be
    # applied to the blade geometry and loading vector in the namelist file. If the collective pitch CB is
    # unnecessary, then dFx =dQ/(Nb*r*R).
    # dFx = dT/Nb*np.sin(-th[0])+dQ/(Nb*r*R)*np.cos(-th[0])
    dFx = dQ/(Nb*r*R)
    #   Figure of merit, induced power factor = 1.15
    FM = CP / (1.15 * CP + sol / 8 * CD)

    #   Sets any infinite values of the computed force components (primarily near the blade root) equal to zero.
    dFx[np.where(np.isnan(dFx) == 1)] = 0
    dFz[np.where(np.isnan(dFz) == 1)] = 0
    dFy = np.zeros(len(r))

    #   if the rotor is rotating CW the force distributions are flipped along the longitudinal axis of the rotor disk.
    if UserIn['rotation'] == 2:
        dFx = -dFx

    dFz = dFz2
    # dFz = np.zeros(np.shape(dFz2))
    dFx = dFx2
    dFy = np.zeros(np.shape(dFz2))

    #%%
    # Assembles all computed load parameters into a dictionary
    loadParams = {'coll_residuals':trim_sol.fun,'th': th, 'beta': [0, 0, 0], 'CT': CT, 'T': T, 'dCT': dCT, 'dT': dT, 'CP': CP, 'P': P,
                  'Q': Q, 'dCP': dCP, 'dQ': dQ, 'dCL': dCL, 'dCD': dCD, 'CL': CL, 'CD': CD, 'FM': FM, 'AoA': AoA,'ClaDist':XsecPolarExp['Lift Slope'], 'lambda': lam,
                  'dFx': dFx, 'dFy': dFy, 'dFz': dFz, 'omega':omega*60/(2*np.pi),'U':U,'dFz2':dFz2,'dFx2':dFx2,'phi':phi}
    return loadParams

#
# %% # figdir = os.path.abspath(os.path.join(input.dirDataFile,'Figures/CL.png')) # with cbook.get_sample_data(figdir) as
#
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    # quant = np.diff(out[1],axis = 0)
    # # levels = np.linspace(-.005, .005, 50)
    # levels = np.linspace(np.min(quant),np.max(quant),50)
    # dist = ax.contourf(phi[:-1], geomParams['rdim'], quant.transpose(),levels = levels)
    # ax.set_ylim(geomParams['rdim'][0],geomParams['rdim'][-1])
    # cbar = fig.colorbar(dist)
    # cbar.ax.set_ylabel('$dFx \: [N]$')
#
# #%%
#     fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
#     quant = dFz2
#     # levels = np.linspace(-.005, .005, 50)
#     levels = np.linspace(np.min(quant),np.max(quant),50)
#     dist = ax.contourf(phi, geomParams['rdim'], quant.transpose(),levels = levels)
#     ax.set_ylim(geomParams['rdim'][0],geomParams['rdim'][-1])
#     cbar = fig.colorbar(dist)
#     cbar.ax.set_ylabel('$dFx \: [N]$')
# #%%
#     fig, ax = plt.subplots(1,1)
#     ax.plot(phi[:-1],np.diff(dCL2,axis = 0)[:,40])
#     ax.plot(phi[:-1],np.diff(dCL2_smooth,axis = 0)[:,40])
#     fig, ax = plt.subplots(1,1)
#     ax.plot(phi,dCL2[:,40])
#     ax.plot(phi,dCL2_smooth[:,40])

#
# #%%
#     fig, ax = plt.subplots(1,1)
#     ax.plot(phi[:-1],np.diff(dFz2,axis = 0)[:,40])
#     fig, ax = plt.subplots(1,1)
#     ax.plot(phi,dFz2[:,40])
# # #
#     import matplotlib.pyplot as plt
#     fig = plt.figure(figsize=[6.4, 4.5], )
#     ax = fig.gca()
#     plt.plot(r, dCL)
#     ax.set_ylabel('Lift Coefficient')
#     ax.set_xlabel('Nondimensional radial position, r/R')
#     ax.set_title('CT/$\sigma$=0.01')
#     plt.grid()
#
#     fig = plt.figure(figsize=[6.4, 4.5], )
#     ax = fig.gca()
#     plt.plot(r, dCD)
#     ax.set_ylabel('Drag Coefficient')
#     ax.set_xlabel('Nondimensional radial position, r/R')
#     ax.set_title('CT/$\sigma$=0.01')
#     plt.grid()
#     plt.axes([0.25 ,1, 0.4 ,0.9])
#     [lam, F, err, i] = TipLoss(lamInit,thInit,r)
#     lam[np.where(np.isnan(lam) == 1)] = 0
#     AoA = th_init-lam/r
#     [Cl,Cd]=PolarLookup(AoA)
#     dCT = 0.5*sig*Cl*r**2*dr
#     CT_temp = np.trapz(dCT)

#
# # errCT= abs((CT_temp-CT)/CT_temp)
# # th_init = 6*CT_temp/(np.mean(chord)*a)+3/2*np.sqrt(CT_temp/2)
# # ii = ii+1
