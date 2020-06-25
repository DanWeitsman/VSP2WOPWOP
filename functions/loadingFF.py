#       VSP2WOPWOP Forward Flight Loading Calculation

#   Author: Daniel Weitsman

#   This function trims a rotor for operating in forward flight by equating the drag and side force to zero.
#   The periodic blade loads are also computed.
#%%

import numpy as np
from scipy.optimize import least_squares
#%%

def loadingFF(UserIn,geomParams,XsecPolar,W, omega,Vx,Vz,alphaShaft):

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
            mu = U / (omega * R) * np.cos(alpha)
            lam = lamTPP_init-mu*beta[1]

            theta_expanded = np.expand_dims(th[0],axis=0)+np.expand_dims(th[1]*np.cos(phi),axis=1)+np.expand_dims(th[2]*np.sin(phi),axis = 1)
            beta_expanded = beta[0]+beta[1]*np.cos(phi)+beta[2]*np.sin(phi)

            ut = r + mu * np.expand_dims(np.sin(phi), axis = 1)
            up = lam + r * np.expand_dims(beta[2] * np.cos(phi) - beta[1] * np.sin(phi), axis = 1) + np.expand_dims(beta_expanded * mu * np.cos(phi), axis = 1)

            ct_dist = 1/(4*np.pi)*solDist*a*(ut**2*theta_expanded-ut*up)
            CT = np.trapz(np.trapz(ct_dist,r,axis=1),phi)

            lamTTP_temp = lam_fixed_pnt(lamTPP_init,mu,alpha,CT)

            err = np.abs((lamTTP_temp - lamTPP_init) / lamTTP_temp)
            lamTPP_init = lamTTP_temp
            mu_x = mu

        return beta,alpha,mu,CT,lamTTP_temp,theta_expanded,beta_expanded,ut,up

    def residuals(th,mu_x,lamTPP_init):
        '''
        This function computes the residuals from the target trim variables.

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes.
        :param mu_x: advance ratio
        :param lamTPP_init: initial estimate for the inflow ratio
        :return: difference between the trim targets and computes CT, beta1c, and beta1s.
        '''

        trim_out = WT_trim(th,mu_x,lamTPP_init)
        res = trimTargs - np.array([trim_out[3], trim_out[0][1], trim_out[0][2]])
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

        distCT = 1 / (4 * np.pi) * solDist * a * (ut ** 2 * theta_expanded - ut * up)
        CT = np.trapz(np.trapz(distCT, r, axis=1), phi)

        distCH = 1/(4*np.pi)*solDist*a*((up * ut * theta_expanded - up ** 2 + cd0 / a * ut ** 2) * np.expand_dims(np.sin(phi), axis = 1)
                                        - np.expand_dims(beta_expanded * np.cos(phi), axis = 1) * (ut ** 2 * theta_expanded - up * ut))
        CH = np.trapz(np.trapz(distCH, r, axis=1), phi)
        CH_TTP = CH+beta[1]*CT

        distCY = 1/(4*np.pi)*solDist*a*(-(up * ut * theta_expanded - up ** 2 + cd0 / a * ut ** 2) * np.expand_dims(np.sin(phi), axis = 1) - np.expand_dims(beta_expanded * np.sin(phi), axis = 1) * (ut ** 2 * theta_expanded - up * ut))
        CY = np.trapz(np.trapz(distCY, r, axis=1), phi)
        CY_TTP = CY+beta[2]*CT

        distCQ = 1/(4*np.pi)*solDist*a*r*(up * ut * theta_expanded - up ** 2 + cd0 / a * ut ** 2)
        CQ = np.trapz(np.trapz(distCQ, r, axis=1), phi)

        distCMX = solDist * a / (2*gamma) * (nuBeta**2-1-3/2*e/R) * beta[2] + e / R * 1 / (4*np.pi) * solDist * a * (ut ** 2 * theta_expanded - up * ut) * np.expand_dims(np.cos(phi), axis = 1)
        CMX = np.trapz(np.trapz(distCMX, r, axis=1), phi)

        distCMY = -solDist * a / (2*gamma) * (nuBeta**2-1-3/2*e/R) * beta[1] + e / R * 1 / (4*np.pi) * solDist * a * (ut ** 2 * theta_expanded - up * ut) * np.expand_dims(np.sin(phi), axis = 1)
        CMY = np.trapz(np.trapz(distCMY, r, axis=1), phi)

        return np.array([CT,CH,CY,CQ,CMX,CMY])

#%%
    '''
    This block of code processes the input parameters and determines the initial estimates for the trim procedure 
    '''

    thFP = np.arctan(Vz/Vx)

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

    alphaInit = alphaShaft+thFP
    U = np.linalg.norm((Vx,Vz))
    mu_z = U/(omega*R) * np.sin(alphaInit)
    mu_x = U/(omega*R) * np.cos(alphaInit)

    phiRes = 361
    phi = np.linspace(0,2*np.pi,phiRes)

    np.zeros((len(phi),len(r)))

    th0 = UserIn['thetaInit'] * (np.pi / 180)
    th1c = 1* (np.pi / 180)
    th1s = 1* (np.pi / 180)
    thInit = np.array([th0,th1c,th1s])
    th = thInit
    # nuBeta = np.sqrt(1+3/2*(e/R))
    nuBeta = UserIn['nuBeta']

    targ_CT = W/(rho*np.pi*R**2*(omega*R)**2)
    targ_beta1c = 0
    targ_beta1s = 0
    trimTargs = [targ_CT,targ_beta1c,targ_beta1s]

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
    #  initial estimate of the inflow ratio
    lamTPP_init = mu_x * np.tan(alphaInit)
    #   Improved estimate of the inflow ratio after going through the fixed point procedure
    lamTPP_init = lam_fixed_pnt(lamTPP_init, mu_x, alphaInit, targ_CT)
    #   employs non-linear least square optimization method (LM) to compute the necessary collective and cyclic pitch settings to mimimize the residuals.
    trim_sol = least_squares(residuals, th, args = [mu_x,lamTPP_init] ,method='lm')
    #   overwrites the initial guessed pitch settings with the computed.
    th = trim_sol.x

    if np.any(trim_sol.fun > 1e-8):
        raise NameError('Caution: large residuals, solution is not converged!')

    #%%
    #todo resolve ut and up in TPP, incorporate inflow model

    beta,alpha,mu,CT, lamTPP, theta_expanded,beta_expanded,ut,up = WT_trim(th,mu_x,lamTPP_init)
    aeroloads = loads_moments(ut, up, beta, theta_expanded, beta_expanded)

    #   dimensionalized tangential, normal, and radial velocity components in order to compute the periodic blade load distribution
    UT = ut*(omega*R)
    UP = up*(omega*R)
    UR = mu_x*np.cos(phi)*(omega*R)

    dFz = (0.5*rho*a*geomParams['chordDist']*((UT ** 2 * theta_expanded - UT * UP)*np.cos(UP/UT)-UT**2*cd0/a*np.sin(UP/UT)))*np.expand_dims(np.cos(beta_expanded),axis = 1)
    T = Nb / (2 * np.pi) * np.trapz(np.trapz(dFz, geomParams['rdim'], axis=1), phi)

    dFy = (0.5*rho*a*geomParams['chordDist']*(-(UT ** 2 * theta_expanded - UT * UP)*np.expand_dims(np.sin(beta_expanded),axis = 1)+UT**2*cd0/a*np.sin(np.expand_dims(UR,axis=1)/UT)))

    dFx = (0.5*rho*a*geomParams['chordDist']*((UT ** 2 * theta_expanded - UT * UP)*np.sin(UP/UT)+UT**2*cd0/a*np.cos(UP/UT)*np.cos(np.expand_dims(UR,axis=1)/UT)))

    Q = Nb / (2 * np.pi) * np.trapz(np.trapz(geomParams['rdim']*dFx, geomParams['rdim'], axis=1), phi)
    P = Q*omega

#%%
    # assembles a dictionary with the computed parameters that is returned to the user and is referenced in other segments of the program
    loadParams = {'residuals':trim_sol.fun,'phiRes':phiRes,'lamTPP': lamTPP ,'gamma':gamma,'mu_x':mu_x,'phi':phi,'th':th,'beta':beta,'CT':aeroloads[0],'T':T,'CH':aeroloads[1],'CY':aeroloads[2],'CQ':aeroloads[3],'Q':Q,'P':P,
                  'CMX':aeroloads[4],'CMY':aeroloads[5],'UP':UP,'UT':UT,'dFx':dFx,'dFy':dFy,'dFz':dFz}

    return loadParams
