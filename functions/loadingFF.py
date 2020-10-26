#       VSP2WOPWOP Forward Flight Loading Calculation

#   Author: Daniel Weitsman

#   This function trims a rotor for operating in forward flight by equating the drag and side force to zero.
#   The periodic blade loads are also computed.
#%%
import bisect
import time
import numpy as np
from scipy.optimize import least_squares,minimize
#%%

def loadingFF(UserIn,geomParams,XsecPolar,W, omega,Vx,Vz,alphaShaft):

    def searchPolar(AoA):
        '''
        This function uses the bisection search algorithm to locate the  lift and drag coefficients from the XFoil
        polar for a given angle of attack (AoA). If the angle of attack exceeds that corresponding to the maximum
        AoA, the maximum lift and drag coefficients would be set.

        :param AoA: matrix filled with the angles of attack
        [rad] for each azimuthal and radial station
        :return CL: matrix of corresponding lift coefficients, having the same dimensions as AoA.
        :return CD: matrix of corresponding drag coefficients, having the same
        dimensions as AoA
        '''

        CL = np.zeros(np.shape(AoA))
        CD = np.zeros(np.shape(AoA))
        for i, azimuth in enumerate(AoA):
            for ii, radial in enumerate(azimuth):
                if radial >= XsecPolar['b540ols']['alphaMax']:
                    # CL[i, ii] = XsecPolar['b540ols']['ClMax']
                    CL[i, ii] = 0
                    CD[i, ii] = XsecPolar['b540ols']['CdMax']
                else:
                    ind = bisect.bisect_right(XsecPolar['b540ols']['Polar'][:, 0], radial)
                    CL[i,ii] = XsecPolar['b540ols']['Polar'][ind,1]
                    CD[i, ii] = XsecPolar['b540ols']['Polar'][ind, 2]
        return CL, CD

    def beta_solve(th,mu_x,lamTPP):
        '''
        This function solves for the flapping angles by evaluating the aerodynamic flap moment
        about the hinge of the blade, while retaining only the first harmonics of the flap and pitch response (HELICOPTER DYNAMICS by Chopra et al, (Chapter 1.3.4)).

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes, expressed in radians.
        :param mu_x: the advance ratio
        :param lamTPP: a constant inflow ratio with respect to the tip path plane (TTP)
        :return: an array consisting of the coning, longitudinal, and lateral flapping angles, expressed in radians.
        '''

        beta0 = gamma/nuBeta**2*(th[0]/8*(1-(e/R)**4)-geomParams['twistDist']/8*(1-(e/R)**4)+th[0]/8*mu_x**2*(1-(e/R)**2)
                                 -geomParams['twistDist']/8*mu_x**2*(1-(e/R)**2)+th[2]/6*mu_x*(1-(e/R)**3)-lamTPP/6*(1-(e/R)**3))

        A = np.array([[(nuBeta**2-1),(gamma/16*mu_x**2*(1-(e/R)**2)+gamma/8*(1-(e/R)**4))],
             [(gamma/16*mu_x**2*(1-(e/R)**2)-gamma/8*(1-(e/R)**4)),(nuBeta**2-1)]])

        B = np.array([gamma*(th[1]/8*(1-(e/R)**4)+th[1]/16*mu_x**2*(1-(e/R)**2)-beta0/6*mu_x*(1-(e/R)**3)),
             gamma*(th[0]/3*mu_x*(1-(e/R)**3)-geomParams['twistDist']/3*mu_x*(1-(e/R)**3)+th[2]/8*(1-(e/R)**4)
                    +3/16*th[2]*mu_x**2*(1-(e/R)**2)-lamTPP/4*mu_x*(1-(e/R)**2))])

        beta1c,beta1s = np.linalg.solve(A,B)

        return np.array([beta0[0],beta1c[0],beta1s[0]])

    def lam_fixed_pnt(lam,mu,alpha,CT):
        '''
        This function applies the fixed point iteration method to converge the inflow ratio. At this stage only uniform inflow
         is supported with the intent of incorporating an inflow model in the future.

        :param lam: the estimate of the inflow ratio
        :param mu: the advance ratio
        :param alpha: angle of attack of the rotor disk
        :param CT: thrust/weight coefficient
        :return: converged inflow ratio
        '''
        errFP = 1
        iii = 0
        while np.any(errFP > 0.0005):

            lam_temp = mu * np.tan(alpha) + CT / (2 * np.sqrt(mu ** 2 + lam ** 2))
            errFP = np.abs((lam_temp - lam) / lam_temp)
            lam = lam_temp
            iii = iii+1

        return lam

    def WT_trim(th,mu_x,lamTPP_init):
        '''
        This function performs a wind tunnel trim on the rotor, whereby the trim targets are the thrust coefficient, longitudinal, and lateral flapping angles.
        Since the flapping angles are dependent on the advance ratio and inflow ratio the trim procedure is solved itterativley until the inflow ratio is converged.

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes, expressed in radians.
        :param mu_x: the advance ratio
        :param lamTPP_init: a constant inflow ratio with respect to the tip path plane (TTP)

        :return:
        :param beta: converged array consisting of the coning, longitudinal, and lateral flapping angles, expressed in radians.
        :param alpha: convered rotor disk angle of attack, expressed in radians
        :param mu: converged advance ratio
        :param CT: converged thrust coefficient
        :param lamTTP_temp: converged inflow ratio
        :param theta_expanded: expanded form of the pitch variations, accounting for first harmonic fluctuations in cyclic pitch (len(phi)xlen(r)).
        :param beta_expanded: expanded form of the flap variations, accounting for first harmonic fluctuations in longitudinal and lateral flapping
        :param ut: nondimensionalized tangential velocity component, evaluated with respect to the hub plane.
        :param up: nondimensionalized normal velocity component, evaluated with respect to the hub plane.
        '''

        err = 1
        while np.any(err > 0.0001):

            beta = beta_solve(th, mu_x, lamTPP_init)
            alpha = alphaShaft+beta[1]+thFP

            # if alpha <-20*np.pi/180:
            #     alpha =-20*np.pi/180

            mu = U / (omega * R) * np.cos(alpha)
            lam = lamTPP_init-mu*beta[1]

            theta_expanded = geomParams['twistDist']+np.expand_dims(th[0],axis=0)+np.expand_dims(th[1]*np.cos(phi),axis=1)+np.expand_dims(th[2]*np.sin(phi),axis = 1)
            beta_expanded = beta[0]+beta[1]*np.cos(phi)+beta[2]*np.sin(phi)

            ut = r + mu * np.expand_dims(np.sin(phi), axis = 1)
            up = lam + r * np.expand_dims(beta[2] * np.cos(phi) - beta[1] * np.sin(phi), axis = 1) + np.expand_dims(beta_expanded * mu * np.cos(phi), axis = 1)

            # todo  insert condition for determining stalled sections if AoA > AoA_max Cl = 0, Cd = cd_max
            AoA = theta_expanded-up/ut
            # CL,CD = searchPolar(AoA)
            CL = a*(theta_expanded-up/ut)
            # ct_dist = 1/2*solDist*ut**2*(CL*np.cos(up/ut)-CD*np.sin(up/ut))*np.expand_dims(np.cos(beta_expanded),axis = 1)
            ct_dist = 1/2*solDist*ut**2*(CL*np.cos(up/ut)-0.1*CL*np.sin(up/ut))

            # ct_dist = 1/(4*np.pi) * solDist * ut ** 2 * a * (theta_expanded - up / ut)
            # ct_dist = 1/(4*np.pi)*solDist*a*(ut**2*theta_expanded-ut*up)
            CT = 1/(2*np.pi)*np.trapz(np.trapz(ct_dist,r),phi)

            lamTTP_temp = lam_fixed_pnt(lamTPP_init,mu,alpha,CT)

            err = np.abs((lamTTP_temp - lamTPP_init) / lamTTP_temp)
            lamTPP_init = lamTTP_temp
            mu_x = mu

        return beta,alpha,mu,CT,ct_dist,lamTTP_temp,theta_expanded,beta_expanded,ut,up,CL

    def fixed_pitch_trim(omega, lamTPP_init):

        mu = U / (omega * R) * np.cos(alphaInit)
        CT_temp =  W/(rho * np.pi * R ** 2 * (omega * R) ** 2)
        i = 0
        err = 1
        while err > 0.0001:

            up = lam_fixed_pnt(lamTPP_init, mu, alphaInit, CT_temp)
            ut = r + mu * np.expand_dims(np.sin(phi), axis=1)
            CT_temp_dist = 1/2*solDist*r**2*(a*(geomParams['twistDist']-up/ut)*np.cos(up/ut)-cd0*np.sin(up/ut))
            CT_temp = 1 / (2 * np.pi) * np.trapz(np.trapz(CT_temp_dist, r), phi)
            err = np.abs((lamTPP_init - up) / lamTPP_init)
            lamTPP_init = up
            i =+1

        T = CT_temp * rho * np.pi * R ** 2 * (omega * R) ** 2
        # print(omega)
        return T,up,ut,CT_temp_dist,CT_temp

    def variable_pitch_residuals(th,*args):
        '''
        This function computes the residuals from the target trim variables.

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes.
        :param mu_x: advance ratio
        :param lamTPP_init: initial estimate for the inflow ratio
        :return: difference between the trim targets and computes CT, beta1c, and beta1s.
        '''
        weights = np.array([1, 1, 1])
        trimOut = WT_trim(th,mu_x,lamTPP_init)
        res = (trimTargs - np.array([trimOut[3], trimOut[0][1], trimOut[0][2]]))*weights
        print(res)
        return res

    def fixed_pitch_residuals(omega, lamTPP_init):
        '''
        This function computes the residuals from the target trim variables.

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes.
        :param mu_x: advance ratio
        :param lamTPP_init: initial estimate for the inflow ratio
        :return: difference between the trim targets and computes CT, beta1c, and beta1s.
        '''

        trimOut = fixed_pitch_trim(omega, lamTPP_init)
        lamTPP_init = trimOut[1]
        res = trimTargs - trimOut[0]
        print(res)
        return res

    def loads_moments(ut, up, beta, theta_expanded, beta_expanded):
        '''
        This function computes the averaged hub loads and moments.

        :param ut: nondimensionalized tangential velocity component, evaluated with respect to the hub plane.
        :param up: nondimensionalized normal velocity component, evaluated with respect to the hub plane.
        :param beta: array consisting of the coning, longitudinal, and lateral flapping angles.
        :param theta_expanded: expanded form of the pitch variations, accounting for first harmonic fluctuations in cyclic pitch (len(phi)xlen(r)).
        :param beta_expanded: expanded form of the flap variations, accounting for first harmonic fluctuations in longitudinal and lateral flapping
        :return:
        :param CT: averaged thrust coefficient
        :param CH: averaged rotor drag force coefficient
        :param CY: averaged side force coefficient
        :param CQ: averaged torque coefficient
        :param CMX: averaged roll moment coefficient
        :param CMY: averaged pitch moment coefficient
        '''

        gam = np.transpose(mu_x*np.cos(phi)/(np.expand_dims(r,axis = 1)+mu_x*np.sin(phi)))
        distCH = 0.5*solDist*ut**2*(a*(theta_expanded-up/ut)*(np.expand_dims(np.sin(phi),axis = 1)*np.sin(up/ut)-np.sin(beta_expanded)*np.expand_dims(np.cos(phi),axis = 1))+cd0*(np.cos(up/ut)*np.cos(gam)*np.expand_dims(np.sin(phi),axis = 1)+np.sin(gam)*np.expand_dims(np.cos(phi),axis = 1)))
        CH = 1/(2*np.pi)*np.trapz(np.trapz(distCH,r),phi)
        CH_TTP = CH+beta[1]*CT
        H = rho * np.pi * R ** 2 * (omega * R) ** 2 * CH

        # distCY = 0.5*solDist*a*(-(up * ut * theta_expanded - up ** 2 + cd0 / a * ut ** 2) * np.expand_dims(np.cos(phi), axis = 1) - np.expand_dims(beta_expanded * np.sin(phi), axis = 1) * (ut ** 2 * theta_expanded - up * ut))
        distCY = 0.5*solDist*ut**2*(-a*(theta_expanded-up/ut)*(np.expand_dims(np.cos(phi),axis = 1)*np.sin(up/ut)+np.sin(beta_expanded)*np.expand_dims(np.sin(phi),axis = 1))+cd0*(np.cos(up/ut)*np.cos(gam)*np.expand_dims(np.cos(phi),axis = 1)+np.sin(gam)*np.expand_dims(np.sin(phi),axis = 1)))
        CY = 1/(2*np.pi)*np.trapz(np.trapz(distCY,r),phi)
        CY_TTP = CY+beta[2]*CT
        Y = rho * np.pi * R ** 2 * (omega * R) ** 2 * CY

        distCMX = solDist * a / (2*gamma) * (nuBeta**2-1-3/2*e/R) * beta[2] + e / R * 1 / (4*np.pi) * solDist * a * (ut ** 2 * theta_expanded - up * ut) * np.expand_dims(np.cos(phi), axis = 1)
        CMX = np.trapz(np.trapz(distCMX, r), phi)
        MX = rho * np.pi * R ** 3 * (omega * R) ** 2 * CMX

        distCMY = -solDist * a / (2*gamma) * (nuBeta**2-1-3/2*e/R) * beta[1] + e / R * 1 / (4*np.pi) * solDist * a * (ut ** 2 * theta_expanded - up * ut) * np.expand_dims(np.sin(phi), axis = 1)
        CMY = np.trapz(np.trapz(distCMY, r), phi)
        MY = rho * np.pi * R ** 3 * (omega * R) ** 2 * CMY

        return H,CH,Y,CY,MX,CMX,MY,CMY

#%%
    '''
    This block of code processes the input parameters and determines the initial estimates for the trim procedure 
    '''

    omega = omega/60*2*np.pi
    rho = UserIn['rho']
    Nb = UserIn['Nb']
    R = geomParams['R']
    e = geomParams['e']
    r = geomParams['r']
    sig = geomParams['solidity']
    solDist = geomParams['solDist']
    XsecLocation = UserIn['XsecLocation']
    airfoilName = list(XsecPolar.keys())

    alphaShaft = alphaShaft*(np.pi/180)
    thFP = np.arctan(Vz / Vx)

    alphaInit = alphaShaft+thFP
    U = np.linalg.norm((Vx,Vz))
    mu_z = U/(omega*R) * np.sin(alphaInit)
    mu_x = U/(omega*R) * np.cos(alphaInit)

    phiRes = 361
    phi = np.linspace(0,2*np.pi,phiRes)

    th0 = UserIn['thetaInit'] * (np.pi / 180)
    th1c = 1* (np.pi / 180)
    th1s = 1* (np.pi / 180)
    thInit = np.array([th0,th1c,th1s])

    # nuBeta = np.sqrt(1+3/2*(e/R))
    nuBeta = UserIn['nuBeta']

    targ_CT = W/(rho*np.pi*R**2*(omega*R)**2)
    targ_beta1c = 0
    targ_beta1s = 0
    if UserIn['trim'] == 3:
        trimTargs = [targ_CT,targ_beta1c,targ_beta1s]
    else:
        trimTargs = W

    # if -2 < Vz/np.sqrt(W/(2*rho*np.pi*R**2)) < 0:
    #     raise ValueError('Non-physical solution, 1D assumption of momentum theory is violated')

#%%
    '''
    This block of code populates an array of indices corresponding to the variety of airfoil sections used along the blade span.   
    '''
    if len(XsecLocation) > 1:
        polarInd = []
        for ii in range(1, len(XsecLocation)):
            polarInd.append(np.squeeze(np.where((XsecLocation[ii - 1] <= np.round(r, 5)) == (XsecLocation[ii] >= np.round(r, 5)))))
        polarInd.append(np.arange(polarInd[-1][-1] + 1, len(r)))
    else:
        polarInd = np.arange(0,len(r))

    a = np.ones((len(r)))
    a0 = np.ones((len(r)))
    cd0 = np.ones((len(r)))

    if len(XsecLocation) > 1:
        for ii in range(0, len(polarInd)):
            a[polarInd[ii]] = XsecPolar[airfoilName[ii]]['Lift Slope']
            a0[polarInd[ii]] = XsecPolar[airfoilName[ii]]['Alpha0'] * (np.pi / 180)
            cd0[polarInd[ii]] = XsecPolar[airfoilName[ii]]['CdMin']
    else:
        a = a*XsecPolar[airfoilName[0]]['Lift Slope']
        a0 = a0*XsecPolar[airfoilName[0]]['Alpha0'] * (np.pi / 180)
        cd0 = cd0*XsecPolar[airfoilName[0]]['CdMin']

    gamma = np.mean((rho*a*geomParams['chordDist']*R**4)/UserIn['Ib'])

    #%%
    #   Initial estimate  of the inflow ratio after completing the fixed point iteration procedure
    lamTPP_init = lam_fixed_pnt(mu_x * np.tan(alphaInit), mu_x, alphaInit, targ_CT)

    # Employs non-linear least square optimization method (LM) to compute the necessary rpm (fixed pitch) or
    # collective and cyclic pitch settings (variable pitch) to minimize the residuals of the trim targets. Separate
    # functions were written for each variety of pitch since the independent variable (th and omega),
    # which is optimized must be the first argument of each of these functions in order for scipy.least_squares to
    # run properly.
    if UserIn['trim'] == 3:
        trim_sol = least_squares(variable_pitch_residuals, thInit, args = [targ_CT,omega,mu_x,lamTPP_init] ,method='lm',diff_step=0.5)
        th = trim_sol.x
        #   Run WT_trim once again with the trimmed values to return quantities necessary in computing the blade loads
        beta, alpha, mu, CT, distCT, lamTPP, theta_expanded, beta_expanded, ut, up,CL = WT_trim(th, mu_x, lamTPP_init)

    else:
        trim_sol = least_squares(fixed_pitch_residuals,omega, args = [lamTPP_init],method = 'lm')
        omega = trim_sol.x
        theta_expanded = geomParams['twistDist']
        beta_expanded = np.zeros(1)
        beta = np.zeros(3)
        T,up,ut,distCT,CT = fixed_pitch_trim(omega,lamTPP_init)

    #%%

    #   Run to return hub loads
    # hubLM = loads_moments(ut, up, beta, theta_expanded, beta_expanded)

    #   Dimensionalized tangential, normal, and radial velocities
    UT = ut*(omega*R)
    UP = up*(omega*R)
    UR = mu_x*np.cos(phi)*(omega*R)

    U = np.sqrt(UT**2+UP**2)
    AoA = theta_expanded-UP/UT

    #   Computes the distributed and integrated thrust as well as the vertical force component
    dT = rho*np.pi*R**2*(omega*R)**2*distCT
    T = 1/(2*np.pi)*np.trapz(np.trapz(dT,r),phi)
    dFz = dT/Nb

    #   Computes distributed and integrated thrust as well as the in-plane force component
    distCQ = 0.5*solDist*r**3*(CL*np.sin(up/ut)+cd0*np.cos(up/ut))
    CQ = 1/(2*np.pi)*np.trapz(np.trapz(distCQ,r),phi)
    dQ = rho*np.pi*R**3*(omega*R)**2*distCQ
    Q = 1/(2*np.pi)*np.trapz(np.trapz(dQ,r),phi)

    #   Power required
    P = Q * omega
    dFx = dQ / (Nb * r * R)
    #   Radial force component, computed using eq. 1.31 of HELICOPTER DYNAMICS (2011, Chopra et al.)
    # dFy = 0.5*rho*geomParams['chordDist']*U**2*(-a*(theta_expanded-up/ut)*np.expand_dims(np.sin(beta_expanded),axis = 1)+cd0*np.sin(np.expand_dims(UR,axis=1)/UT))
    dFy = np.zeros((np.shape(dFz)))


#%%
    # assembles a dictionary with the computed parameters that is returned to the user and is referenced in other segments of the program
    loadParams = {'residuals':trim_sol.fun,'phiRes':phiRes,'ClaDist':a,'AoA':AoA,'lamTPP': lamTPP ,'alpha':alpha,'gamma':gamma,'mu_x':mu_x,'phi':phi,'th':th,'beta':beta,'CT':CT,'T':T,'CQ':CQ,'Q':Q,'P':P,
                  'UP':UP,'UT':UT,'U':U,'dFx':dFx,'dFy':dFy,'dFz':dFz}

    return loadParams
    #
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    # dist = ax.contourf(phi, r, np.expand_dims(geomParams['twistDist'],axis = 1)-np.transpose(up/ut))
    levels = np.linspace(np.min(dFz),np.max(dFz),200)
    dist = ax.contourf(phi, r, np.transpose((dFz)),levels = levels)
    plt.colorbar(dist)
