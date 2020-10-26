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

def loadingFFv2(UserIn,geomParams,XsecPolar,W, omega,Vx,Vz,alphaShaft):


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


    def variable_pitch_residuals(th, mu_x, lamTPP_init):
        '''
        This function computes the residuals from the target trim variables.

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes.
        :param mu_x: advance ratio
        :param lamTPP_init: initial estimate for the inflow ratio
        :return: difference between the trim targets and computes CT, beta1c, and beta1s.
        '''

        if UserIn['trim'] == 2:
            trimOut = cyclic_pitch_trim([th,0,0], mu_x, lamTPP_init)
            res = trimTargs - trimOut[0]
        else:
            trimOut = cyclic_pitch_trim(th, mu_x, lamTPP_init)
            res = trimTargs - np.array([trimOut[0], trimOut[3], trimOut[4]])
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
            AoA = (geomParams['twistDist']-up/ut)%(2*np.pi)
            CL,CD = searchPolar(AoA)
            CT_temp_dist = 1/2*solDist*r**2*(CL*np.cos(up/ut)-CD*AoA*np.sin(up/ut))
            CT_temp = 1 / (2 * np.pi) * np.trapz(np.trapz(CT_temp_dist, r), phi)
            err = np.abs((up - lamInit) / up)
            lamInit = up
            i =+1

        T = CT_temp * rho * np.pi * R ** 2 * (omega * R) ** 2

        return T,up,ut,CT_temp_dist,CT_temp,AoA,CL,CD, mu

    def coll_pitch_trim(th, mu, lamInit,CT_init):

        i = 0
        err = 1
        while err > 0.0005:
            up = lam_fixed_pnt(lamInit, mu, alphaInit, CT_init)
            ut = r + mu * np.expand_dims(np.sin(phi), axis=1)
            AoA = (th+geomParams['twistDist'] - up / ut) % (2 * np.pi)
            CL, CD = aeroParams(AoA)
            CT_temp_dist = 1 / 2 * solDist * r ** 2 * (CL * np.cos(up / ut) - CD * AoA * np.sin(up / ut))
            CT_temp = 1 / (2 * np.pi) * np.trapz(np.trapz(CT_temp_dist, r), phi)
            err = np.abs((CT_temp - CT_init) / CT_temp)
            CT_init = CT_temp
            lamInit = up
            i = +1

        T = CT_temp * rho * np.pi * R ** 2 * (omega * R) ** 2

        return T, up, ut, CT_temp_dist, CT_temp, AoA, CL, CD, mu

    def cyclic_pitch_trim(th, mu_x, lamTPP_init):
        '''
        This function performs a wind tunnel trim on the rotor, whereby the thrust coefficient, roll, and pitching moments are the trim targets.

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

            theta_expanded = geomParams['twistDist']+th[0]+np.expand_dims(th[1]*np.cos(phi),axis=1)+np.expand_dims(th[2]*np.sin(phi),axis = 1)
            ut = r + mu_x * np.expand_dims(np.sin(phi), axis = 1)
            up = lamTPP_init

            AoA = (theta_expanded-up/ut)%(2*np.pi)
            CL,CD = aeroParams(AoA)

            dCT = 1/2*solDist*r**2*(CL*np.cos(up/ut)-CD*np.sin(up/ut))
            CT = 1/(2*np.pi)*np.trapz(np.trapz(dCT,r),phi)

            lamTTP_temp = lam_fixed_pnt(lamTPP_init,mu_x,alphaInit,CT)
            err = np.abs((lamTTP_temp - lamTPP_init) / lamTTP_temp)
            lamTPP_init = lamTTP_temp

        dFz = rho * np.pi * R ** 2 * (omega * R) ** 2 * dCT/Nb
        Mx = Nb/(2*np.pi)*np.trapz(np.trapz(geomParams['rdim']*dFz*np.expand_dims(np.sin(phi),axis = 1),r),phi)
        My = -Nb/(2*np.pi)*np.trapz(np.trapz(geomParams['rdim']*dFz*np.expand_dims(np.cos(phi),axis = 1),r),phi)

        return CT,dCT,dFz,Mx,My,lamTTP_temp,theta_expanded,ut,up,CL,CD,AoA

#%%
    # def cyclic_pitch_trim(th, mu_x, lamTPP_init):
    #     '''
    #     This function performs a wind tunnel trim on the rotor, whereby the trim targets are the thrust coefficient, longitudinal, and lateral flapping angles.
    #     Since the flapping angles are dependent on the advance ratio and inflow ratio the trim procedure is solved itterativley until the inflow ratio is converged.
    #
    #     :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes, expressed in radians.
    #     :param mu_x: the advance ratio
    #     :param lamTPP_init: a constant inflow ratio with respect to the tip path plane (TTP)
    #
    #     :return:
    #     :param beta: converged array consisting of the coning, longitudinal, and lateral flapping angles, expressed in radians.
    #     :param alpha: convered rotor disk angle of attack, expressed in radians
    #     :param mu: converged advance ratio
    #     :param CT: converged thrust coefficient
    #     :param lamTTP_temp: converged inflow ratio
    #     :param theta_expanded: expanded form of the pitch variations, accounting for first harmonic fluctuations in cyclic pitch (len(phi)xlen(r)).
    #     :param beta_expanded: expanded form of the flap variations, accounting for first harmonic fluctuations in longitudinal and lateral flapping
    #     :param ut: nondimensionalized tangential velocity component, evaluated with respect to the hub plane.
    #     :param up: nondimensionalized normal velocity component, evaluated with respect to the hub plane.
    #     '''
    #
    #     err = 1
    #     while np.any(err > 0.0005):
    #
    #         beta = beta_solve(th, mu_x, lamTPP_init)
    #         alpha = alphaShaft+beta[1]+thFP
    #
    #         mu = U / (omega * R) * np.cos(alpha)
    #         lam = lamTPP_init-mu*beta[1]
    #
    #         theta_expanded = geomParams['twistDist']+th[0]+np.expand_dims(th[1]*np.cos(phi),axis=1)+np.expand_dims(th[2]*np.sin(phi),axis = 1)
    #         beta_expanded = beta[0]+beta[1]*np.cos(phi)+beta[2]*np.sin(phi)
    #
    #         ut = r + mu * np.expand_dims(np.sin(phi), axis = 1)
    #         up = lam + r * np.expand_dims(beta[2] * np.cos(phi) - beta[1] * np.sin(phi), axis = 1) + np.expand_dims(beta_expanded * mu * np.cos(phi), axis = 1)
    #
    #         AoA = (theta_expanded-up/ut)%(2*np.pi)
    #         CL,CD = searchPolar(AoA)
    #
    #         # ct_dist = 1/2*solDist*ut**2*(CL*np.cos(up/ut)-CD*np.sin(up/ut))*np.expand_dims(np.cos(beta_expanded),axis = 1)
    #         ct_dist = 1/2*solDist*r**2*(CL*np.cos(up/ut)-CD*np.sin(up/ut))
    #         CT = 1/(2*np.pi)*np.trapz(np.trapz(ct_dist,r),phi)
    #
    #         lamTTP_temp = lam_fixed_pnt(lamTPP_init,mu,alpha,CT)
    #         err = np.abs((lamTTP_temp - lamTPP_init) / lamTTP_temp)
    #
    #         lamTPP_init = lamTTP_temp
    #         mu_x = mu
    #
    #     return beta,alpha,mu,CT,ct_dist,lamTTP_temp,theta_expanded,beta_expanded,ut,up,CL,CD,AoA

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

#%%
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
                    ind = bisect.bisect_left(XsecPolar[polarInd[ii]]['Polar'][:, 0], radial)
                    if ind == len(XsecPolar[polarInd[ii]]['Polar'][:,0]):
                        ind = 0
                    CL[i,ii] = XsecPolar[polarInd[ii]]['Polar'][ind,1]
                    CD[i,ii] = XsecPolar[polarInd[ii]]['Polar'][ind, 2]

        return CL, CD

    def aeroParams(AoA):
        '''
        This function returns the lift and drag coefficients corresponding to a radial and azimuthal distribution of
        angles of attack. The lift coefficient for stalled blade sections is linearly interpolated between the
        section's airfoil minimum and maximum lift coefficients. The drag coefficient is assumed to be 10% of the
        lift coefficient tunless the blade section is stalled. In that case, the sectional drag coefficient is set to
        the airfoil's drag coefficient at the angle of attack corresponding to the maximum lift coefficient
        :param AoA: array of size [phiRes x len(r)] filled with the computed angles of attack
        :return:
        :param CL: lift coefficient, linearly interpolated for the stalled blade sections
        :param CD:  drag coefficient, set equal to its value at the angle of attack corresponding to themaximum lift
        coefficient for the stalled blade sections
        '''

        #   assume that the airfoil is symmetric and therefore the CL can be estimated by the product of the
        # lift-curve slope and the angle of attack
        CL = XsecPolarExp['Lift Slope'] * AoA
        #   CD is assumed to be 10% of CL
        CD = 0.1 * CL

        #   reruns the indices of stalled blade sections
        azInd, rInd = np.where(AoA > XsecPolarExp['alphaMax'])
        #   sets the CD of these sections equal to the CD @ CLmax
        CD[azInd, rInd] = XsecPolarExp['CdMax'][azInd, rInd]
        #   linearly interpolates CL between CLmin and CL
        CL[azInd, rInd] = np.interp(CL[azInd, rInd], XsecPolarExp['Alpha0'][azInd, rInd],
                                    XsecPolarExp['ClMax'][azInd, rInd])

        return CL,CD

    #%%

    omega = omega/60*2*np.pi
    rho = UserIn['rho']
    Nb = UserIn['Nb']
    R = geomParams['R']
    e = geomParams['e']
    r = geomParams['r']
    solDist = geomParams['solDist']
    XsecLocation = UserIn['XsecLocation']

    alphaShaft = alphaShaft*(np.pi/180)
    thFP = np.arctan(Vz / Vx)
    alphaInit = alphaShaft+thFP
    U = np.linalg.norm((Vx,Vz))

    mu_x = U/(omega*R) * np.cos(alphaInit)

    phiRes = 361
    phi = np.linspace(0,2*np.pi,phiRes)
    a = np.ones((len(r)))*XsecPolar[list(XsecPolar.keys())[0]]['Lift Slope']
    th0 = UserIn['thetaInit']*np.pi/180
    gamma = np.mean((rho * a * geomParams['chordDist'] * R ** 4) / UserIn['Ib'])

# %% This section of code assigns the airfoil parameters from the XFoil polar to the corresponding radial section
    # along the blade span
    XsecPolarExp = {}
    if len(XsecLocation) > 1:
        polarInd = []
        ind = np.zeros((len(XsecLocation)+1))
        for i,Xsec in enumerate(XsecLocation):
            ind[i] = bisect.bisect(r, Xsec)
        ind[0] = 0
        ind[-1] = len(r)

        for i,Xsec in enumerate(XsecPolar.keys()):
            polarInd.extend([Xsec]*int(ind[i+1]-ind[i]))
            for ii, param in enumerate(list(XsecPolar[Xsec].keys())[1:]):
                if i ==0:
                    XsecPolarExp = {**XsecPolarExp, **{param:XsecPolar[Xsec][param]*np.ones(len(r))}}
                else:
                    XsecPolarExp[param][int(ind[i]):] = XsecPolar[Xsec][param]
    else:
        for i,key in enumerate(list(XsecPolar[list(XsecPolar.keys())[0]].keys())[1:]):
            XsecPolarExp[key] = np.ones((phiRes,len(r)))*XsecPolar[list(XsecPolar.keys())[0]][key]


#%%
    if UserIn['trim']==1:
        trimTargs = W
        trim_sol = least_squares(fixed_pitch_residuals, omega, method = 'lm',diff_step = 0.5)
        omega = trim_sol.x
        T, up, ut, dCT, CT, AoA, CL, CD, mu_x = fixed_pitch_trim(omega)
        beta = [0,0,0]
        alpha = alphaInit
        th = [th0,0,0]

    elif UserIn['trim'] == 2:
        trimTargs = W/(rho*np.pi*R**2*(omega*R)**2)
        lamTPP_init =  lam_fixed_pnt(mu_x*np.tan(alphaInit), mu_x, alphaInit, trimTargs)
        trim_sol = least_squares(variable_pitch_residuals, th0, args=[mu_x, lamTPP_init], method='lm', diff_step=0.5)
        th = [trim_sol.x,0,0]
        CT,dCT,dFz,Mx,My,lamTTP_temp,theta_expanded,ut,up,CL,CD,AoA = cyclic_pitch_trim(th,mu_x,lamTPP_init)

    else:
        trimTargs = [W/(rho*np.pi*R**2*(omega*R)**2),0,0]
        th = np.array([th0,np.pi/180,np.pi/180])
        lamTPP_init =  lam_fixed_pnt(mu_x*np.tan(alphaInit), mu_x, alphaInit, trimTargs[0])
        trim_sol = least_squares(variable_pitch_residuals, th ,args = [mu_x, lamTPP_init],method = 'lm',diff_step=0.5)
        th = trim_sol.x
        CT,dCT,dFz,Mx,My,lamTTP_temp,theta_expanded,ut,up,CL,CD,AoA = cyclic_pitch_trim(th,mu_x,lamTPP_init)

    #%%

    # if np.sum(trim_sol.fun/trimTargs > 1e-4):
    #     raise ValueError('caution, trim solution is not converged')


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
    # dFr = rho*np.pi*R**2*(omega*R)**2*(1/2*solDist*r**2*(-CL*np.expand_dims(np.sin(beta_exp),axis = 1)+CD*np.sin(np.expand_dims(mu_x*np.cos(phi),axis = 1)/ut)))
    dFr = np.zeros((np.shape(dFz)))

    if UserIn['rotation'] == 2:
        dFz = np.flip(dFz,axis = 0)
        dFx = np.flip(dFx, axis=0)

    #   hub force
    H = Nb/(2*np.pi)*np.trapz(np.trapz((dFr*np.expand_dims(np.cos(phi),axis = 1)+dFx*np.expand_dims(np.sin(phi),axis = 1)),r),phi)
    #   side force
    Y = Nb/(2*np.pi)*np.trapz(np.trapz((dFr*np.expand_dims(np.sin(phi),axis = 1)-dFx*np.expand_dims(np.cos(phi),axis = 1)),r),phi)
    #   roll moment
    Mx = Nb/(2*np.pi)*np.trapz((UserIn['nuBeta']**2-1-3/2*e/R)*UserIn['Ib']*omega**2*beta_exp*np.sin(phi),phi)+Nb/(2*np.pi)*np.trapz(np.trapz(geomParams['rdim']*dFz*np.expand_dims(np.sin(phi),axis = 1),r),phi)
    #    pitching moment
    My = -Nb/(2*np.pi)*np.trapz((UserIn['nuBeta']**2-1-3/2*e/R)*UserIn['Ib']*omega**2*beta_exp*np.cos(phi),phi)-Nb/(2*np.pi)*np.trapz(np.trapz(geomParams['rdim']*dFz*np.expand_dims(np.cos(phi),axis = 1),r),phi)
    hubLM = [H,Y,Mx,My]

    #   assembles a dictionary with the computed parameters that is returned to the user and is referenced in other segments of the program
    loadParams = {'residuals':trim_sol.fun,'phiRes':phiRes,'ClaDist':a,'AoA':AoA,'alpha':alpha,'gamma':gamma,'mu_x':mu_x,'phi':phi,'th':th,'beta':beta,'CT':CT,'T':T,'CQ':CQ,'Q':Q,'P':P,
                  'UP':UP,'UT':UT,'U':U,'dFx':dFx,'dFy':dFr,'dFz':dFz,'hubLM':hubLM}
    #
    return loadParams


    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    quant = dFz
    levels = np.linspace(np.min(quant),np.max(quant),50)
    dist = ax.contourf(phi, r, np.transpose(quant),levels = levels)
    cbar = fig.colorbar(dist)
    cbar.ax.set_ylabel('dFz')
    #
