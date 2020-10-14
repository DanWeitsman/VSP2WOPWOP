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

def fixPitchLoadFF(UserIn,geomParams,XsecPolar,W, omega,Vx,Vz,alphaShaft):


    def fixed_pitch_residuals(omega):
        '''
        This function computes the residuals from the target trim variables.

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes.
        :param mu_x: advance ratio
        :param lamTPP_init: initial estimate for the inflow ratio
        :return: difference between the trim targets and computes CT, beta1c, and beta1s.
        '''
        trimOut = fixed_pitch_trim(omega)
        res = trimTargs - trimOut[0]
        print(res)
        return res


    def fixed_pitch_trim(omega):

        mu = U / (omega * R) * np.cos(alphaInit)
        CT_temp =  W/(rho * np.pi * R ** 2 * (omega * R) ** 2)
        lamInit = lam_fixed_pnt(mu * np.tan(alphaInit), mu, alphaInit, CT_temp)
        i = 0
        err = 1
        while err > 0.0005:

            up = lam_fixed_pnt(lamInit, mu, alphaInit, CT_temp)
            ut = r + mu * np.expand_dims(np.sin(phi), axis=1)
            AoA = (th0+geomParams['twistDist']-up/ut)%(2*np.pi)
            CL,CD = searchPolar(AoA)
            CT_temp_dist = 1/2*solDist*r**2*(CL*np.cos(up/ut)-CD*AoA*np.sin(up/ut))
            CT_temp = 1 / (2 * np.pi) * np.trapz(np.trapz(CT_temp_dist, r), phi)
            err = np.abs((up - lamInit) / up)
            lamInit = up
            i =+1

        T = CT_temp * rho * np.pi * R ** 2 * (omega * R) ** 2

        return T,up,ut,CT_temp_dist,CT_temp,AoA,CL,CD, mu

    def cyclic_pitch_residuals(th, mu_x, lamTPP_init):
        '''
        This function computes the residuals from the target trim variables.

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes.
        :param mu_x: advance ratio
        :param lamTPP_init: initial estimate for the inflow ratio
        :return: difference between the trim targets and computes CT, beta1c, and beta1s.
        '''
        # weights = np.ones(len(trimTargs))
        # trimOut = cyclic_pitch_trim(th,mu_x, lamTPP_init)

        if UserIn['trim']==1:
            trimOut = cyclic_pitch_trim([th,0,0], mu_x, lamTPP_init)
            res = trimTargs - trimOut[3]
        else:
            trimOut = cyclic_pitch_trim(th, mu_x, lamTPP_init)
            res = trimTargs - np.array([trimOut[3], trimOut[0][1], trimOut[0][2]])
        print(res)
        return res

    def cyclic_pitch_trim(th, mu_x, lamTPP_init):
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
        while np.any(err > 0.0005):

            beta = beta_solve(th, mu_x, lamTPP_init)
            alpha = alphaShaft+beta[1]+thFP

            mu = U / (omega * R) * np.cos(alpha)
            lam = lamTPP_init-mu*beta[1]

            theta_expanded = geomParams['twistDist']+th[0]+np.expand_dims(th[1]*np.cos(phi),axis=1)+np.expand_dims(th[2]*np.sin(phi),axis = 1)
            beta_expanded = beta[0]+beta[1]*np.cos(phi)+beta[2]*np.sin(phi)

            ut = r + mu * np.expand_dims(np.sin(phi), axis = 1)
            up = lam + r * np.expand_dims(beta[2] * np.cos(phi) - beta[1] * np.sin(phi), axis = 1) + np.expand_dims(beta_expanded * mu * np.cos(phi), axis = 1)

            AoA = (theta_expanded-up/ut)%(2*np.pi)
            CL,CD = searchPolar(AoA)

            # ct_dist = 1/2*solDist*ut**2*(CL*np.cos(up/ut)-CD*np.sin(up/ut))*np.expand_dims(np.cos(beta_expanded),axis = 1)
            ct_dist = 1/2*solDist*r**2*(CL*np.cos(up/ut)-CD*np.sin(up/ut))
            CT = 1/(2*np.pi)*np.trapz(np.trapz(ct_dist,r),phi)

            lamTTP_temp = lam_fixed_pnt(lamTPP_init,mu,alpha,CT)
            err = np.abs((lamTTP_temp - lamTPP_init) / lamTTP_temp)

            lamTPP_init = lamTTP_temp
            mu_x = mu

        return beta,alpha,mu,CT,ct_dist,lamTTP_temp,theta_expanded,beta_expanded,ut,up,CL,CD,AoA

    def beta_solve(th,mu_x,lamTPP):
        '''
        This function solves for the flapping angles by evaluating the aerodynamic flap moment
        about the hinge of the blade, while retaining only the first harmonics of the flap and pitch response (HELICOPTER DYNAMICS by Chopra et al, (Chapter 1.3.4)).

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes, expressed in radians.
        :param mu_x: the advance ratio
        :param lamTPP: a constant inflow ratio with respect to the tip path plane (TTP)
        :return: an array consisting of the coning, longitudinal, and lateral flapping angles, expressed in radians.
        '''

        beta0 = gamma/UserIn['nuBeta']**2*(th[0]/8*(1-(e/R)**4)-geomParams['twistDist']/8*(1-(e/R)**4)+th[0]/8*mu_x**2*(1-(e/R)**2)
                                 -geomParams['twistDist']/8*mu_x**2*(1-(e/R)**2)+th[2]/6*mu_x*(1-(e/R)**3)-lamTPP/6*(1-(e/R)**3))

        A = np.array([[(UserIn['nuBeta']**2-1),(gamma/16*mu_x**2*(1-(e/R)**2)+gamma/8*(1-(e/R)**4))],
             [(gamma/16*mu_x**2*(1-(e/R)**2)-gamma/8*(1-(e/R)**4)),(UserIn['nuBeta']**2-1)]])

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
        while np.any(errFP > 0.0005):
            lam_temp = mu * np.tan(alpha) + CT / (2 * np.sqrt(mu ** 2 + lam ** 2))
            errFP = np.abs((lam_temp - lam) / lam_temp)
            lam = lam_temp
        return lam

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
                    ind = bisect.bisect_left(XsecPolar['b540ols']['Polar'][:, 0], radial)
                    if ind == len(XsecPolar['b540ols']['Polar'][:,0]):
                        CL[i,ii] = XsecPolar['b540ols']['Polar'][0,1]
                        CD[i,ii] = XsecPolar['b540ols']['Polar'][0, 2]
                    else:
                        CL[i,ii] = XsecPolar['b540ols']['Polar'][ind,1]
                        CD[i,ii] = XsecPolar['b540ols']['Polar'][ind, 2]
        return CL, CD

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
        distCH = 0.5*solDist*ut**2*(CL*(np.expand_dims(np.sin(phi),axis = 1)*np.sin(up/ut)-np.sin(beta_expanded)*np.expand_dims(np.cos(phi),axis = 1))+CD*(np.cos(up/ut)*np.cos(gam)*np.expand_dims(np.sin(phi),axis = 1)+np.sin(gam)*np.expand_dims(np.cos(phi),axis = 1)))
        CH = 1/(2*np.pi)*np.trapz(np.trapz(distCH,r),phi)
        CH_TTP = CH+beta[1]*CT
        H = rho * np.pi * R ** 2 * (omega * R) ** 2 * CH

        # distCY = 0.5*solDist*a*(-(up * ut * theta_expanded - up ** 2 + cd0 / a * ut ** 2) * np.expand_dims(np.cos(phi), axis = 1) - np.expand_dims(beta_expanded * np.sin(phi), axis = 1) * (ut ** 2 * theta_expanded - up * ut))
        distCY = 0.5*solDist*ut**2*(-CL*(np.expand_dims(np.cos(phi),axis = 1)*np.sin(up/ut)+np.sin(beta_expanded)*np.expand_dims(np.sin(phi),axis = 1))+CD*(np.cos(up/ut)*np.cos(gam)*np.expand_dims(np.cos(phi),axis = 1)+np.sin(gam)*np.expand_dims(np.sin(phi),axis = 1)))
        CY = 1/(2*np.pi)*np.trapz(np.trapz(distCY,r),phi)
        CY_TTP = CY+beta[2]*CT
        Y = rho * np.pi * R ** 2 * (omega * R) ** 2 * CY

        distCMX = solDist * a / (2*gamma) * (UserIn['nuBeta']**2-1-3/2*e/R) * beta[2] + e / R * 1 / (4*np.pi) * solDist * a * (ut ** 2 * theta_expanded - up * ut) * np.expand_dims(np.cos(phi), axis = 1)
        CMX = np.trapz(np.trapz(distCMX, r), phi)
        MX = rho * np.pi * R ** 3 * (omega * R) ** 2 * CMX

        distCMY = -solDist * a / (2*gamma) * (UserIn['nuBeta']**2-1-3/2*e/R) * beta[1] + e / R * 1 / (4*np.pi) * solDist * a * (ut ** 2 * theta_expanded - up * ut) * np.expand_dims(np.sin(phi), axis = 1)
        CMY = np.trapz(np.trapz(distCMY, r), phi)
        MY = rho * np.pi * R ** 3 * (omega * R) ** 2 * CMY

        return H,CH,Y,CY,MX,CMX,MY,CMY

#%%

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
    a = np.ones((len(r)))*XsecPolar[list(XsecPolar.keys())[0]]['Lift Slope']
    th0 = UserIn['thetaInit']*np.pi/180

    gamma = np.mean((rho * a * geomParams['chordDist'] * R ** 4) / UserIn['Ib'])

    if UserIn['trim']==0:
        trimTargs = W
        trim_sol = least_squares(fixed_pitch_residuals, omega, method = 'lm',diff_step = 0.5)
        omega = trim_sol.x
        T, up, ut, dCT, CT, AoA, CL, CD, mu_x = fixed_pitch_trim(omega)
        beta = [0,0,0]

    elif UserIn['trim']==1:
        trimTargs = W/(rho*np.pi*R**2*(omega*R)**2)
        lamTPP_init =  lam_fixed_pnt(mu_x*np.tan(alphaInit), mu_x, alphaInit, trimTargs)
        th = np.array([th0, 0, 0])
        trim_sol = least_squares(cyclic_pitch_residuals, th[0], args=[mu_x, lamTPP_init], method='lm', diff_step=0.5)
        th = trim_sol.x
        beta, alpha, mu, CT, dCT, lamTTP, theta_expanded, beta_expanded, ut, up, CL, CD, AoA = cyclic_pitch_trim([th,0,0], mu_x, lamTPP_init)

    else:
        trimTargs = [W/(rho*np.pi*R**2*(omega*R)**2),0,0]
        th = np.array([th0,np.pi/180,np.pi/180])
        lamTPP_init =  lam_fixed_pnt(mu_x*np.tan(alphaInit), mu_x, alphaInit, trimTargs[0])
        trim_sol = least_squares(cyclic_pitch_residuals, th ,args = [mu_x, lamTPP_init],method = 'lm',diff_step = 0.5)
        th = trim_sol.x
        beta,alpha,mu,CT,dCT,lamTTP,theta_expanded,beta_expanded,ut,up,CL,CD,AoA = cyclic_pitch_trim(th,mu_x,lamTPP_init)

    #%%
    beta_exp = beta[0]+beta[1]*np.cos(phi)+beta[2]*np.sin(phi)

    UT = ut*omega*R
    UP = up * omega * R
    U = np.sqrt(UT**2+UP**2)

    dT = rho*np.pi*R**2*(omega*R)**2*dCT
    T = 1/(2*np.pi)*np.trapz(np.trapz(dT,r),phi)

    dCQ = 0.5*solDist*r**3*(CL*np.sin(up/ut)+CD*np.cos(up/ut))
    CQ = 1/(2*np.pi)*np.trapz(np.trapz(dCQ,r),phi)
    dQ = rho*np.pi*R**3*(omega*R)**2*dCQ
    Q = 1/(2*np.pi)*np.trapz(np.trapz(dQ,r),phi)
    P = Q * omega

    dFz = dT / Nb
    dFx = dQ / (Nb * r * R)
    dFr = rho*np.pi*R**2*(omega*R)**2*(1/2*solDist*r**2*(-CL*np.expand_dims(np.sin(beta_exp),axis = 1)+CD*np.sin(np.expand_dims(mu_x*np.cos(phi),axis = 1)/ut)))

    H = Nb/(2*np.pi)*np.trapz(np.trapz((dFr*np.expand_dims(np.cos(phi),axis = 1)+dFx*np.expand_dims(np.sin(phi),axis = 1)),r),phi)
    Y = Nb/(2*np.pi)*np.trapz(np.trapz((dFr*np.expand_dims(np.sin(phi),axis = 1)-dFx*np.expand_dims(np.cos(phi),axis = 1)),r),phi)
    Mx = Nb/(2*np.pi)*np.trapz((UserIn['nuBeta']**2-1-3/2*e/R)*UserIn['Ib']*omega**2*beta_exp*np.sin(phi),phi)+Nb/(2*np.pi)*np.trapz(np.trapz(e*dFz*np.expand_dims(np.sin(phi),axis = 1),r),phi)
    My = -Nb/(2*np.pi)*np.trapz((UserIn['nuBeta']**2-1-3/2*e/R)*UserIn['Ib']*omega**2*beta_exp*np.cos(phi),phi)-Nb/(2*np.pi)*np.trapz(np.trapz(e*dFz*np.expand_dims(np.cos(phi),axis = 1),r),phi)

    hubLM = [H,Y,Mx,My]

    #   assembles a dictionary with the computed parameters that is returned to the user and is referenced in other segments of the program
    loadParams = {'residuals':trim_sol.fun,'phiRes':phiRes,'ClaDist':a,'AoA':AoA,'alpha':alpha,'gamma':gamma,'mu_x':mu_x,'phi':phi,'th':th,'beta':beta,'CT':CT,'T':T,'CQ':CQ,'Q':Q,'P':P,
                  'UP':UP,'UT':UT,'U':U,'dFx':dFx,'dFy':dFr,'dFz':dFz,'hubLM':hubLM}
    #
    return loadParams

    import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    # # dist = ax.contourf(phi, r, AOA)
    # # levels = np.linspace(-1,5,20)
    # dist = ax.contourf(phi, r, np.transpose(dFy))
    # plt.colorbar(dist)

