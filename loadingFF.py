#       VSP2WOPWOP Forward Flight Loading Calculation

#   Author: Daniel Weitsman

#   This function trims a rotor for operating in forward flight by equating the drag and side force to zero.
#   The periodic blade loads are also computed.

#%%

def loadingFF(UserIn, geomParams, XsecPolar, W, omega, Vx, Vz, alphaShaft):

    import bisect
    import numpy as np
    from scipy.optimize import least_squares
    from scipy.interpolate import interp1d

    def fixed_pitch_residuals(omega):
        '''
        This function computes the residuals between the trim targets and computed trim variables for rpm trim.

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes.
        :param mu_x: advance ratio
        :param lamTPP_init: initial estimate for the inflow ratio
        :return: difference between the trim targets and computes CT, beta1c, and beta1s.
        '''
        trimOut = fixed_pitch_trim(omega)
        res = trimTargs - trimOut[0]
        print(f'Trim residuals: T = {round(res,6)}N')
        return res


    def variable_pitch_residuals(th,mu, lamTPP_init):
        '''
        This function computes the residuals between the trim targets and computed trim variables for collective and cyclic pitch trim.

        :param th: an array of three elements the first being the collective pitch setting, followed by the lateral and longituinal cyclic pitch amplitudes.
        :param mu_x: advance ratio
        :param lamTPP_init: initial estimate for the inflow ratio
        :return: difference between the trim targets and computes CT, beta1c, and beta1s.
        '''

        if UserIn['trim'] == 2:
            print(th*180/np.pi)
            trimOut = variable_pitch_trim([th,0,0], mu, lamTPP_init)
            res = trimTargs - trimOut[0]
            print(f'Trim residuals: dCT = {round(res, 6)}')
        else:
            trimOut = variable_pitch_trim(th, mu, lamTPP_init)
            res = trimTargs - np.array([trimOut[0], trimOut[2], trimOut[3]])
            print(f'Trim residuals: dCT = {round(res[0], 6)}, CMx = {round(res[1], 6)}, CMy = {round(res[2], 6)}')
        return res

    def fixed_pitch_trim(omega):
        '''
        This function performs an rpm trim, whereby the rotational rate is varied until the desired thrust is achieved.
        :param omega: rotational rate [rad/s]
        :return:
        '''

        mu = U / (omega * R)
        CT =  W/(rho * np.pi * R ** 2 * (omega * R) ** 2)
        lamTPP_init = inflowModSelect(2, mu*np.tan(alphaInit), mu, CT)
        dCT = CT*np.ones(np.shape(lamTPP_init))

        err = 1
        while np.any(err > 0.0005):

            up = inflowModSelect(UserIn['inflowMod'],lamTPP_init, mu, CT,dCT)
            ut = r + mu * np.expand_dims(np.sin(phi), axis=1)
            AoA = th0+geomParams['twistDist']-up/ut
            CL,CD = np.array([aeroParams(x) for x in AoA]).transpose(1,0,2)
            AoA = th0+geomParams['twistDist']-up/ut
            CL,CD = aeroParams(AoA)
            dCT = 1/2*solDist*r**2*(CL*np.cos(up/ut)-CD*AoA*np.sin(up/ut))
            CT = 1 / (2 * np.pi) * np.trapz(np.trapz(dCT, r), phi)
            err = np.abs((up - lamTPP_init) / up)
            lamTPP_init = up

        T = CT * rho * np.pi * R ** 2 * (omega * R) ** 2

        return T,CT,dCT,lamTPP_init,ut,up,CL,CD,AoA,mu

    def variable_pitch_trim(th,mu, lamTPP_init):
        '''
        This function performs a collective/cyclic pitch trim, whereby the thrust coefficient, roll, and pitching moments are the trim targets.

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
        i = 0
        while np.any(err > 0.0005):

            theta_expanded = geomParams['twistDist']+th[0]+np.expand_dims(th[1]*np.cos(phi),axis=1)+np.expand_dims(th[2]*np.sin(phi),axis = 1)
            ut = r + mu*np.cos(alphaInit) * np.expand_dims(np.sin(phi), axis = 1)

            AoA = theta_expanded-up/ut

            CL,CD = np.array([aeroParams(x) for x in AoA]).transpose(1,0,2)

            dCT = 0.5*solDist*r**2*CL
            CT = 1/(2*np.pi)*np.trapz(np.trapz(dCT,r),phi)

            # dCT = 1/2*solDist*r**2*(CL*np.cos(up/ut)-CD*np.sin(up/ut))
            # CT = 1/(2*np.pi)*np.trapz(np.trapz(dCT,r),phi)

            lamTTP_temp = inflowModSelect(UserIn['inflowMod'], lamTPP_init, mu, CT, dCT)
            err = np.abs((lamTTP_temp - lamTPP_init) / lamTTP_temp)
            lamTPP_init = lamTTP_temp

        Mx = 1/(2*np.pi)*np.trapz(np.trapz(r*dCT*np.expand_dims(np.sin(phi),axis = 1),r),phi)*rho*(omega*R)**2*np.pi*R**3
        My = -1/(2*np.pi)*np.trapz(np.trapz(r*dCT*np.expand_dims(np.cos(phi),axis = 1),r),phi)*rho*(omega*R)**2*np.pi*R**3

        return CT,dCT,Mx,My,lamTTP_temp,theta_expanded,ut,up,CL,AoA

    def inflowModSelect(model, lam, mu,CT, *args):
        '''
        This function selects and returns the converged inflow distribution based on the model specified in the user input module.
        :param model: integer specifing the model selected in the input module (UserIn['inflowMod'])
        :param lam: initial guess for the inflow distribution, can be an arbitrary sized array
        :param mu: = standard advance ratio (V/(omega*R)), unresolved into parallel and perpendicular components to the rotor disk.
        :param CT: thrust coefficient
        :return:
        '''
        if model ==1:
            lam = constant_inflow(lam, mu, CT)
        elif model == 2:
            lam = linear_inflow(lam, mu, CT)
        elif model == 3:
            lam = drees_inflow(lam, mu, CT)
        elif model == 4:
            lam = pitt_peters_inflow(lam, mu, args[0])
        return lam

    def constant_inflow(lam,mu,CT):
        '''
        This function applies the fixed point iteration method to converge the constant inflow ratio.

        :param lam: the estimate of the inflow ratio
        :param mu: the advance ratio
        :param alpha: angle of attack of the rotor disk
        :param CT: thrust/weight coefficient
        :return: converged inflow ratio
        '''
        errFP = 1
        # mu = mu * np.cos(alphaInit)
        while np.any(errFP > 1e-5):
            lam_temp =  CT / (2 * np.sqrt(mu * np.cos(alphaInit) ** 2 + (lam +mu*np.sin(alphaInit))** 2))
            errFP = np.abs((lam_temp - lam) / lam_temp)
            lam = lam_temp
        return lam

    def linear_inflow(lam,mu,CT):
        '''
        This function utilizes the fixed point itteration method to converge the Glauert's linear inflow model.

        :param lam: the estimate of the inflow ratio
        :param mu: the advance ratio
        :param alpha: angle of attack of the rotor disk
        :param CT: thrust/weight coefficient
        :return: converged inflow ratio
        '''
        err = 1
        # mu = mu*np.cos(alphaInit)
        while np.any(err > 0.0005):
            lam_temp = CT / (2 * np.sqrt(mu*np.cos(alphaInit) ** 2 + lam ** 2))*(1+1.2*r*np.expand_dims(np.cos(phi),axis = 1))
            err = np.abs((lam_temp - lam) / lam_temp)
            lam = lam_temp
        return lam

    def drees_inflow(lam,mu,CT):
        '''
        This function utilizes the fixed point itteration method to converge the Drees's inflow model.

        :param lam: the estimate of the inflow ratio
        :param mu: the advance ratio
        :param alpha: angle of attack of the rotor disk
        :param CT: thrust/weight coefficient
        :return: converged inflow ratio
        '''
        err = 1
        while np.any(err > 0.0005):
            wake_skew = np.arctan(mu*np.cos(alphaInit)/lam)
            kx = 4/3*((1-np.cos(wake_skew)-1.8*mu**2)/np.sin(wake_skew))
            ky = -2*mu
            lam_temp = CT / (2 * np.sqrt(mu ** 2 + lam ** 2))*(1+kx*r*np.expand_dims(np.cos(phi),axis = 1)+ky*r*np.expand_dims(np.sin(phi),axis = 1))
            err = np.abs((lam_temp - lam) / lam_temp)
            lam = lam_temp
        return lam

    def pitt_peters_inflow(lam, mu,CT,dCT):
        '''
        This function computes the inflow distribution based on the steady component of the Pitt-Peters model. This
        formulation was originally presented in, Pitt, Dale M., and David A. Peters. "Theoretical prediction of
        dynamic-inflow derivatives." (1980) and then again in Chen, Robert TN. "A survey of nonuniform inflow models
        for rotorcraft flight dynamics and control applications." (1989). This model takes into account the effects
        that the hub moments have on the steady (1st harmonic) induced velocity distribution. This model should be
        used when performing the fixed/collective pitch trim since the inflow distribution would inevitably vary in
        order to produce the necessary reaction to counteract the hub moments.
        :param lam: initial guess for the inflow distribution
        :param mu: advance ratio
        :param dCT: radial and azimuthal distribution of the thrust coefficient
        '''

        CMX = 1/(2 * np.pi) * np.trapz(np.trapz(r * dCT  * np.expand_dims(np.sin(phi), axis=1), r), phi)
        CMY = -1/(2 * np.pi) * np.trapz(np.trapz(r * dCT * np.expand_dims(np.cos(phi), axis=1), r), pphisi)

        lam = constant_inflow(np.sqrt(CT/2), mu, CT)
        wake_skew = np.arctan(mu*np.cos(alphaInit)/lam)
        vt = np.sqrt((mu*np.cos(alphaInit))**2+lam**2)
        vm = ((mu*np.cos(alphaInit))**2+lam*(lam+CT/(2*np.sqrt(mu**2 + lam**2))))/vt

        L = np.array([[0.5/vt,0,15*np.pi/(64*vm)*np.tan(wake_skew/2)],[0,-4/(vm*(1+np.cos(wake_skew))),0],[15*np.pi/(64*vt)*np.tan(wake_skew/2),0,-4*np.cos(wake_skew)/(vm*(1+np.cos(wake_skew)))]])
        lam_0,lam_1c,lam_1s = np.dot(L,[CT,CMX,CMY])
        lam = lam_0 + lam_1c*r*np.expand_dims(np.cos(phi),axis = 1)+ lam_1s*r*np.expand_dims(np.sin(phi),axis = 1)

        return lam

    def aeroParams(AoA):
        """
        This function linearly interpolates the sectional airfoil lift and drag coefficients from the XFoil polar based on the
        computed angle of attack distribution. It is assumed that the airfoils are not stalled in the region of angle of attack ranging from -alpha_max+alpha0 to alpha_max. 
        In this region, the lift coefficient variation is reflected and negated so that it varies from -CL_max to CL_max. While the drag coefficient is only reflected.
        
        For stalled blade sections the coefficient of lift is linearly interpolated from CL_max to -CL_max, while the drag coefficient is set to its value corresponding to the
        the stall angle of attack (alpha_max).

        :param alpha: angle of attack distribution [rad]
        return:
        :param CL:  lift coefficient distribution
        :param CD:  drag coefficient distribution
        """

        # allocate empty arrays for CL and CD that have the same dimensions as AoA. 
        CL = np.empty(AoA.shape)
        CD = np.empty(AoA.shape)

        # loops through each airfoil section that makes up the blade.
        for k,v in XsecPolar.items():
            
            # turns polarInd list into a 2d array that has the same dimensions as AoA so that it can be used as an index. 
            airfoil_key_expanded= np.repeat(np.expand_dims(np.array(polarInd),axis = 0),len(AoA),axis = 0)
            
            # determines the prestall region ranging from alpha0 to alphaMax. 
            interp_range = np.array([np.squeeze(np.where(v['Polar'][:,0]==v['Alpha0'])),np.squeeze(np.where(v['Polar'][:,0]==v['alphaMax']))])
            # concatenates the negated version of interp_range to account for negative values of the angles of attack. 
            x = (np.concatenate((-(v['Polar'][interp_range[0]:interp_range[-1],0]-v['Alpha0'])[::-1][:-1],v['Polar'][interp_range[0]:interp_range[-1],0]-v['Alpha0']))+v['Alpha0'])
            # negates, flips and concatenates the variation in the lift cofficients that correspond to negative values of the angles of attack. 
            CL_y = np.concatenate((-(v['Polar'][interp_range[0]:interp_range[-1],1])[::-1][:-1],v['Polar'][interp_range[0]:interp_range[-1],1]))
            # flips and concatenates the variations in the drag cofficients that correspond to negative values of the angles of attack. 
            CD_y = np.concatenate((v['Polar'][interp_range[0]:interp_range[-1],2][::-1][:-1],v['Polar'][interp_range[0]:interp_range[-1],2]))

            # creates the linear interpolation function for the pre-stalled blade sections. 
            CL_prestall_interp = interp1d(x,CL_y)
            CD_prestall_interp = interp1d(x,CD_y)

            # indecies corresponding to the pre-stalled blade sections. Criteria: AoA >= -alpha_max+alpha0 and AoA <= alpha_max
            prestall_ind = np.array((AoA >= x[0]) & (AoA <= x[-1]) & (k == airfoil_key_expanded))

            #   performs the linear interpolation on the blade sections that are not stalled
            CL[prestall_ind] = CL_prestall_interp(AoA[prestall_ind])
            CD[prestall_ind] = CD_prestall_interp(AoA[prestall_ind])

            # creates the linear interpolation function for the stalled blade sections, where CL ranges from CL_max to -CL_max from alpha_max to -alpha_max+alpha0. 
            CL_stall_interp = interp1d(np.array([x[-1],x[0]%(2*np.pi)]),np.array([CL_y[-1],CL_y[0]]))

            dCL[ind] = np.interp(AoA[ind],xp = v['Polar'][:,0],fp =v['Polar'][:,1])
            # dCL[ind] = 2*np.pi*AoA[ind]
            dCD[ind] = np.interp(AoA[ind],xp = v['Polar'][:,0],fp =v['Polar'][:,2])
        return dCL, dCD

#%%

    omega = omega/60*2*np.pi
    rho = UserIn['rho']
    Nb = UserIn['Nb']
    R = geomParams['R']
    r = geomParams['r']
    solDist = geomParams['solDist']
    XsecLocation = UserIn['XsecLocation']
    CT = W/(rho*np.pi*R**2*(omega*R)**2)


    alphaShaft = alphaShaft*(np.pi/180)
    thFP = np.arctan(Vz / Vx)
    alphaInit = alphaShaft+thFP
    U = np.linalg.norm((Vx,Vz))

    mu = U/(omega*R)

    phiRes = 361
    phi = np.linspace(0,2*np.pi,phiRes)
    a = np.ones((len(r)))*XsecPolar[list(XsecPolar.keys())[0]]['Lift Slope']
    th0 = UserIn['thetaInit']*np.pi/180

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
        polarInd = list(XsecPolar.keys()) * len(r)
        for i,key in enumerate(list(XsecPolar[list(XsecPolar.keys())[0]].keys())[1:]):
            XsecPolarExp[key] = np.ones((phiRes,len(r)))*XsecPolar[list(XsecPolar.keys())[0]][key]
            polarInd = np.full(len(r),list(XsecPolar.keys())[0])
    # polarInd = np.repeat(np.expand_dims(polarInd,axis = 0),361,axis = 0)


#%%
    if UserIn['trim']==1:
        trimTargs = W
        trim_sol = least_squares(fixed_pitch_residuals, omega, method = 'lm',diff_step = 0.5)
        omega = trim_sol.x
        th = np.array([th0,0 ,0 ])
        T,CT,dCT,lam,ut,up,CL,CD,AoA,mu = fixed_pitch_trim(omega)
        theta_expanded = np.empty((np.shape(AoA)))*th0

    elif UserIn['trim'] == 2:
        trimTargs = CT
        lamTPP_init =  inflowModSelect(1, mu*np.tan(alphaInit), mu, trimTargs)
        trim_sol = least_squares(variable_pitch_residuals, th0, args=[mu, lamTPP_init], method='lm')
        th = np.array([np.squeeze(trim_sol.x),0 ,0 ])
        CT,dCT,CMX,CMY,lam,th_exp,ut,up,CL,AoA = variable_pitch_trim(th,mu, lamTPP_init)

    else:
        trimTargs = [CT,0,0]
        th = np.array([th0,np.pi/180,np.pi/180])
        lamTPP_init =  inflowModSelect(1, mu*np.tan(alphaInit), mu, trimTargs[0])
        trim_sol = least_squares(variable_pitch_residuals, th ,args = [mu, lamTPP_init],method = 'lm')
        th = trim_sol.x
        CT,dCT,CMX,CMY,lam,th_exp,ut,up,CL,AoA,AoA_eff = variable_pitch_trim(th,mu, lamTPP_init)


#%%
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

    # resolves loading vectors to vertical and horizontal directions so that a change of base can be applied to the
    # blade geometry account for the pitching motion in the namelist file - 1/18/21
    dFz =  dT/Nb*np.cos(-theta_expanded)-dQ/(Nb*r*R)*np.sin(-theta_expanded)
    dFx = dT/Nb*np.sin(-theta_expanded)+dQ/(Nb*r*R)*np.cos(-theta_expanded)
    # dFr = rho*np.pi*R**2*(omega*R)**2*(1/2*solDist*r**2*(-CL*np.expand_dims(np.sin(beta_exp),axis = 1)+CD*np.sin(np.expand_dims(mu_x*np.cos(phi),axis = 1)/ut)))
    dFr = np.zeros((np.shape(dFz)))

    #   if the rotor is rotating CW the force distributions are flipped along the longitudinal axis of the rotor disk.
    if UserIn['rotation'] == 2:
        dFz = np.flip(dFz,axis = 0)
        dFx = np.flip(dFx, axis=0)
        AoA = np.flip(AoA, axis=0)
        U = np.flip(U, axis=0)

        if UserIn['inflowMod'] !=1:
            lam = np.flip(lam, axis=0)
        th[2] = -th[2]

    #   hub force
    H = Nb/(2*np.pi)*np.trapz(np.trapz((dFr*np.expand_dims(np.cos(phi),axis = 1)+dFx*np.expand_dims(np.sin(phi),axis = 1)),r),phi)
    #   side force
    Y = Nb/(2*np.pi)*np.trapz(np.trapz((dFr*np.expand_dims(np.sin(phi),axis = 1)-dFx*np.expand_dims(np.cos(phi),axis = 1)),r),phi)
    #   roll moment
    Mx = Nb / (2 * np.pi) * np.trapz(np.trapz(geomParams['rdim'] * dFz * np.expand_dims(np.sin(phi), axis=1), r), phi)
    #   pitch moment
    My = -Nb / (2 * np.pi) * np.trapz(np.trapz(geomParams['rdim'] * dFz * np.expand_dims(np.cos(phi), axis=1), r), phi)
    hubLM = [H,Y,Mx,My]


    #   assembles a dictionary with the computed parameters that is returned to the user and is referenced in other segments of the program
    loadParams = {'residuals':trim_sol.fun,'phiRes':phiRes,'ClaDist':a,'AoA':AoA,'alpha':alphaInit,'mu':mu,'phi':phi,'th':th,'CT':CT,'T':T,'CQ':CQ,'Q':Q,'P':P,
                  'UP':UP,'UT':UT,'U':U,'dFx':dFx,'dFy':dFr,'dFz':dFz,'hubLM':hubLM,'omega':omega*60/(2*np.pi)}
    #
    return loadParams

    #
    # import matplotlib.pyplot as plt
    #
    # fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    # quant = AoA*180/np.pi
    # levels = np.linspace(np.min(quant),np.max(quant),50)
    # dist = ax.contourf(phi, geomParams['rdim'], np.transpose(quant),levels = levels)
    # ax.set_ylim(geomParams['rdim'][0],geomParams['rdim'][-1])
    # cbar = fig.colorbar(dist)
    # # cbar.ax.set_ylabel(r'$\lambda$')
    # cbar.ax.set_ylabel(r'$\alpha \ [\circ]$')
    # plt.show()
    
    #
    # N = 3
    # for i in range(N):
    #     print(th[1]*np.cos(2*np.pi/N*i)+th[2]*np.sin(2*np.pi/N*i))
    # for i in range(N):
    #     print(-th[1]*np.sin(2*np.pi/N*i)+th[2]*np.cos(2*np.pi/N*i))