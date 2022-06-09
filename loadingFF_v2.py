#       VSP2WOPWOP Forward Flight Loading Calculation

#   Author: Daniel Weitsman

#   This function trims a rotor for operating in forward flight by equating the drag and side force to zero.
#   The periodic blade loads are also computed.

#%%

# from curses.ascii import US


def loadingFF(UserIn, geomParams, XsecPolar, W, omega, Vx, Vz, alphaShaft):

    import bisect
    import numpy as np
    from scipy.optimize import least_squares

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
            trimOut = variable_pitch_trim_2([th,0,0], mu)
            res = trimTargs - trimOut[0]
            print(f'Trim residuals: T = {round(res, 6)}N')
        else:
            trimOut = variable_pitch_trim_2(th, mu,lamTPP_init)
            res = trimTargs - np.array([trimOut[0], trimOut[2], trimOut[3]])
            print(f'Trim residuals: T = {round(res[0], 6)}N, Mx = {round(res[1], 6)}Nm, My = {round(res[2], 6)}Nm')
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
            ut = r + mu * np.expand_dims(np.sin(psi), axis=1)
            AoA = (th0+geomParams['twistDist']-up/ut)%(2*np.pi)
            CL,CD = aeroParams(AoA)
            dCT = 1/2*solDist*r**2*(CL*np.cos(up/ut)-CD*AoA*np.sin(up/ut))
            CT = 1 / (2 * np.pi) * np.trapz(np.trapz(dCT, r), psi)
            err = np.abs((up - lamTPP_init) / up)
            lamTPP_init = up

        T = CT * rho * np.pi * R ** 2 * (omega * R) ** 2

        return T,CT,dCT,lamTPP_init,ut,up,CL,CD,AoA,mu

    def variable_pitch_trim(th,mu, lamTPP_init,alphaInit):
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
        AoA_err = 1
        lam_err = 1
        i = 0

        theta_expanded = geomParams['twistDist']+th[0]+np.expand_dims(th[1]*np.cos(psi),axis=1)+np.expand_dims(th[2]*np.sin(psi),axis = 1)

        if UserIn['flap']:

            beta_exp = beta_init[0]+beta_init[1]*np.cos(psi)+beta_init[2]*np.sin(psi)
            ut_init = r + mu*np.cos(alphaInit) * np.expand_dims(np.sin(psi), axis = 1)
            up_init = mu * np.tan(alphaInit)+lamTPP_init + r * np.expand_dims(np.gradient(beta_exp,edge_order = 2), axis=1) / omega + np.expand_dims(mu * np.sin(beta_exp) * np.cos(psi), axis=1)
            AoA_init = (theta_expanded - up_init / ut_init)

            # lamTPP_init = np.ones(np.shape(theta_expanded)) * lamTPP_init
            # up = lamTPP_init + r * np.expand_dims(beta_exp,axis = 1) / omega+ np.expand_dims(mu * np.sin(beta_exp) * np.cos(phi),axis = 1)

            while np.any(AoA_err > 0.0005):

                while np.any(lam_err > 0.0005):
                    up_init = lamTPP_init + r * np.expand_dims(beta_exp, axis=1) / omega + np.expand_dims(mu * np.sin(beta_exp) * np.cos(psi), axis=1)
                    AoA_init = theta_expanded - up_init / ut_init
                    CL, CD = np.array([PolarLookup(x) for x in AoA_init]).transpose(1, 0, 2)
                    dCT = 1 / 2 * solDist * r ** 2 * (CL * np.cos(up_init / ut_init) - CD * np.sin(up_init / ut_init))
                    CT = 1 / (2 * np.pi) * np.trapz(np.trapz(dCT, r), psi)
                    lamTTP_temp = inflowModSelect(UserIn['inflowMod'], lamTPP_init, mu, CT,dCT)
                    lam_err = abs((lamTTP_temp - lamTPP_init) / lamTTP_temp)
                    lamTPP_init = lamTTP_temp

                lam_const = constant_inflow(np.sqrt(CT/2), mu, CT)
                beta0 = gamma*(th[0]/8*(1+mu**2)+th_tw/10*(1+5/6*mu**2)-mu/6*th[2]-lam_const/6)
                beta1s = -(4/3*mu*beta0)/(1+1/2*mu**2)+th[1]
                beta1c = -(8/3*mu*(th[0]-3/4*lam_const+3/4*mu*th[2]+3/4*th_tw)/(1-1/2*mu**2))-th[2]
                beta_exp = beta0 + beta1c * np.cos(psi) + beta1s * np.sin(psi)
                alphaInit = alphaShaft+beta1c+thFP+beta0
                ut = r + mu*np.cos(alphaInit) * np.expand_dims(np.sin(psi), axis = 1)
                up = lamTPP_init + r * np.expand_dims(beta_exp, axis=1) / omega + np.expand_dims(mu * np.sin(beta_exp) * np.cos(psi), axis=1)
                AoA = (theta_expanded - up / ut) % (2 * np.pi)
                
                if UserIn['UnsteadyAero']:
                    da_ds = [(-1/12*AoA[(i+2)%360]+2/3*AoA[(i+1)%360]-2/3*AoA[(i-1)%360]+1/12*AoA[(i-2)%360])/ds[i] for i in range(len(psi)-1)]

                    X_err = 1
                    X = np.empty(np.shape(s))
                    Y = np.empty(np.shape(s))

                    X[0] = A[0]*da_ds[0]
                    Y[0]= A[1]*da_ds[0]

                    while np.any(X_err > 1e-5) or np.any(Y_err > 1e-5):
                        for i in range(len(psi)-1):
                            X[i+1] = X[i]*np.exp(-B[0]*ds[i])+A[0]*da_ds[i]*np.exp(-B[0]*ds[i]/2)
                            Y[i+1] = Y[i]*np.exp(-B[1]*ds[i])+A[1]*da_ds[i]*np.exp(-B[1]*ds[i]/2)
                        X_err = X[-1]-X[0]
                        Y_err= Y[-1]-Y[0]
                        X[0] = X[-1]
                        Y[0] = Y[-1]

                    AoA_eff = AoA-X-Y

                    AoA_err = (AoA_eff - AoA_init) / AoA_eff
                    AoA_init = AoA_eff
  
                else:
                    AoA_err = (AoA - AoA_init) / AoA
                    AoA_init = AoA

                up_init = up
                ut_init = ut

        else:
            while np.any(lam_err > 0.0005):

                # if UserIn['flap']:

                    # beta_exp = beta_init[0]+np.expand_dims(beta_init[1]*np.cos(phi),axis=1)+np.expand_dims(beta_init[2]*np.sin(phi),axis = 1)
                    # ut = r + mu*np.cos(alphaInit) * np.expand_dims(np.sin(phi), axis = 1)
                    # up = lamTPP_init+r*beta_exp/omega+mu*np.sin(beta_exp)*np.cos(phi)
                    # AoA = (theta_expanded - up / ut) % (2 * np.pi)
                    # CL, CD = aeroParams(AoA)
                    # dCT = 1 / 2 * solDist * r ** 2 * (CL * np.cos(up / ut) - CD * np.sin(up / ut))
                    # CT = 1 / (2 * np.pi) * np.trapz(np.trapz(dCT, r), phi)
                    # lamTTP_temp = inflowModSelect(UserIn['inflowMod'], lamTPP_init, mu, CT, dCT)
                    # beta0 = gamma*(th[0]/8*(1+mu**2)+th_twist/10*(1+5/6*mu**2)-mu/6*th[2]-lamTTP_temp/6)
                    # beta1s = -(4/3*mu*beta0)/(1+1/2*mu**2)+th[1]
                    # beta1c = -(8/3*mu*(th[0]-3/4*lamTTP_temp+3/4*mu*th[2]+3/4*th_twist)/(1-1/2*mu**2))-th[2]
                    # beta_exp = beta0+np.expand_dims(beta1c*np.cos(phi),axis=1)+np.expand_dims(beta1s*np.sin(phi),axis = 1)
                    # alphaInit = alphaShaft+beta1c
                    #
                    # err = np.abs((lamTTP_temp - lamTPP_init) / lamTTP_temp)
                    # lamTPP_init = lamTTP_temp

                # else:
                # theta_expanded = geomParams['twistDist']+th[0]+np.expand_dims(th[1]*np.cos(phi),axis=1)+np.expand_dims(th[2]*np.sin(phi),axis = 1)
                ut = r + mu*np.cos(alphaInit) * np.expand_dims(np.sin(psi), axis = 1)
                up = mu * np.tan(alphaInit)+lamTPP_init
                AoA =theta_expanded-up/ut
                CL,CD = np.array([PolarLookup(x) for x in AoA]).transpose(1,0,2)
                dCT = 1/2*solDist*r**2*(CL*np.cos(up/ut)-CD*np.sin(up/ut))
                CT = 1/(2*np.pi)*np.trapz(np.trapz(dCT,r),psi)
                lamTTP_temp = inflowModSelect(UserIn['inflowMod'],lamTPP_init, mu, CT, dCT)
                lam_err = abs((lamTTP_temp - lamTPP_init) / lamTTP_temp)
                lamTPP_init = lamTTP_temp

        # Mx = 1/(2*np.pi)*np.trapz(np.trapz(r*dCT*np.expand_dims(np.sin(phi),axis = 1),r),phi)*rho*(omega*R)**2*np.pi*R**3
        # My = -1/(2*np.pi)*np.trapz(np.trapz(r*dCT*np.expand_dims(np.cos(phi),axis = 1),r),phi)*rho*(omega*R)**2*np.pi*R**3
        Mx = 1/(2*np.pi)*np.trapz(np.trapz(r*dCT*np.expand_dims(np.sin(psi),axis = 1),r),psi)
        My = -1/(2*np.pi)*np.trapz(np.trapz(r*dCT*np.expand_dims(np.cos(psi),axis = 1),r),psi)

        # return CT,dCT,Mx,My,lamTTP_temp,theta_expanded,ut,up,CL,CD,AoA,beta0,beta1c,beta1s,beta_exp,alphaInit
        if UserIn['UnsteadyAero']:
            return CT,dCT,Mx,My,lamTTP_temp,theta_expanded,ut,up,CL,CD,AoA,AoA_eff
        else:
            return CT,dCT,Mx,My,lamTTP_temp,theta_expanded,ut,up,CL,CD,AoA


    def variable_pitch_trim_2(th,mu,lam_init):

        lam_err = 1
        i = 0
        CT = trimTargs[0]
        th_exp = twist+th[0]+np.expand_dims(th[1]*np.cos(psi),axis=1)+np.expand_dims(th[2]*np.sin(psi),axis = 1)
        ut = r+mu*np.expand_dims(np.sin(psi), axis = 1)

        if UserIn['flap']:
            lam_const = constant_inflow(np.sqrt(CT/2),mu,CT)
            beta0 = gamma*(th[0]/8*(1+mu**2)+th_tw/10*(1+5/6*mu**2)-mu/6*th[2]-lam_const/6)
            beta1s = -(4/3*mu*beta0)/(1+1/2*mu**2)+th[1]
            beta1c = -(8/3*mu*(th[0]-3/4*lam_const+3/4*mu*th[2]+3/4*th_tw)/(1-1/2*mu**2))-th[2]
            beta_exp = beta0 + beta1c * np.cos(psi) + beta1s * np.sin(psi)
            up_beta =  r * np.expand_dims(np.gradient(beta_exp,edge_order = 2), axis=1) / omega + np.expand_dims(mu * np.sin(beta_exp) * np.cos(psi), axis=1)

        while np.any(lam_err > 1e-4):

            if UserIn['pwake']:
                up = up_beta+mu*np.tan(alphaInit)+lam_init
                v = np.sqrt(up**2+ut**2)*omega*R
                AoA = th_exp-np.arctan2(up,ut)
                if UserIn['UnsteadyAero']:
                    AoA_eff = unsteady_aero(AoA)
                    gamma_b = np.array([np.matmul(I_tot_inv[:,:,ind],(-v*AoA_eff)[ind]) for ind in range(len(psi))])
                else:
                    gamma_b = np.array([np.matmul(I_tot_inv[:,:,ind],(-v*AoA)[ind]) for ind in range(len(psi))])
                    AoA_eff =0
                # up_temp = (I_fw_tot*np.expand_dims(np.expand_dims(np.amax(gamma_b,axis = 1),axis = -1),axis = -1)).transpose()
                # dv_fw =  I_fw*np.expand_dims(np.expand_dims(np.expand_dims(np.expand_dims(np.amax(gamma_b[:,int(.5*len(r)):],axis = -1),axis = -1),axis = -1),axis = -1),axis = -1)
                # lam = np.sum(np.sum(dv_fw[:,-1]*(omega*R)**-1,axis = -1),axis = -1)
                # dv_fw =  I_fw*np.expand_dims(gamma_b[:,int(.75*len(r))],axis = -1)
                lam = I_fw_tot[:,-1]*np.expand_dims(np.amax(gamma_b[:,int(.25*len(r)):],axis = -1),axis = -1)*(omega*R)**-1
                # lam = np.array([np.sum(np.sum(np.roll(dv_fw[-1,:,:-1],-ind,axis = 1)[:,::int(360/Nb)],axis = -1),axis = -1) for ind in range(len(psi))])/(omega*R)
                # lam = np.array([np.sum(np.sum(dv_fw[-1,:,:-1][:,ind-int(np.floor(ind/90)*360/Nb)::int(360/Nb)],axis = -1),axis = -1) for ind in range(len(psi))])/(omega*R)

            else:
                up = up_beta+mu*np.tan(alphaInit)+lam_init
                AoA = th_exp-np.arctan2(up,ut)
                if UserIn['UnsteadyAero']:
                    AoA_eff = unsteady_aero(AoA)
                    CL = 2*np.pi*AoA_eff
                else:
                    AoA_eff =0
                    CL = 2*np.pi*AoA
                dCT = 0.5*solDist*r**2*CL
                CT = 1/(2*np.pi)*np.trapz(np.trapz(dCT,r),psi)
                lam = inflowModSelect(UserIn['inflowMod'],lam_init, mu, CT, dCT)

            lam_err =abs((lam - lam_init) / lam)
            # lam_err[np.isnan(lam_err)] = 1
            lam_init = lam
            i+=1
            print(i)
            print(np.max(lam_err))

        if UserIn['pwake']:

            CL = 2*gamma_b/(v*chordDist)
            # stall_ind = AoA>(XsecPolarExp['alphaMax']+2*np.pi/180)
            # if np.any(stall_ind):
            #     # dCL = np.interp(AoA[stall_ind]%(2*np.pi),xp = [XsecPolar[list(XsecPolar.keys())[0]]['alphaMax'],2*np.pi],fp = [XsecPolar[list(XsecPolar.keys())[0]]['ClMax'],XsecPolar[list(XsecPolar.keys())[0]]['ClMin']])
            #     CL[stall_ind] = 1.175*np.sin(2*AoA[stall_ind]*180/np.pi)
                
            # CL = XsecPolarExp['Lift Slope']*(AoA-XsecPolarExp['Alpha0'])
            dCT = 0.5*solDist*r**2*CL
            CT = 1/(2*np.pi)*np.trapz(np.trapz(dCT,r),psi)

        # print(f'{np.array([beta0,beta1c,beta1s])*180/np.pi}')
        CMX = 1/(2*np.pi)*np.trapz(np.trapz(r*dCT*np.expand_dims(np.sin(psi),axis = 1),r),psi)
        CMY = -1/(2*np.pi)*np.trapz(np.trapz(r*dCT*np.expand_dims(np.cos(psi),axis = 1),r),psi)


        return CT,dCT,CMX,CMY,lam,th_exp,ut,up,CL,AoA,AoA_eff


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
            lam = pitt_peters_inflow(lam, mu,CT,args[0])
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
            lam_temp =  CT / (2 * np.sqrt(mu ** 2 + (lam +mu_z)** 2))
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
            lam_temp = CT / (2 * np.sqrt(mu ** 2 + lam ** 2))*(1+1.2*r*np.expand_dims(np.cos(psi),axis = 1))
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
            wake_skew = np.arctan(mu/lam)
            kx = 4/3*((1-np.cos(wake_skew)-1.8*mu**2)/np.sin(wake_skew))
            ky = -2*mu
            lam_temp = CT / (2 * np.sqrt(mu ** 2 + lam ** 2))*(1+kx*r*np.expand_dims(np.cos(psi),axis = 1)+ky*r*np.expand_dims(np.sin(psi),axis = 1))
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

        CMX = 1/(2 * np.pi) * np.trapz(np.trapz(r * dCT  * np.expand_dims(np.sin(psi), axis=1), r), psi)
        CMY = -1/(2 * np.pi) * np.trapz(np.trapz(r * dCT * np.expand_dims(np.cos(psi), axis=1), r), psi)

        lam = constant_inflow(np.sqrt(CT/2), mu, CT)
        wake_skew = np.arctan(mu*np.cos(alphaInit)/lam)
        vt = np.sqrt((mu*np.cos(alphaInit))**2+lam**2)
        vm = ((mu*np.cos(alphaInit))**2+lam*(lam+CT/(2*np.sqrt(mu**2 + lam**2))))/vt

        L = np.array([[0.5/vt,0,15*np.pi/(64*vm)*np.tan(wake_skew/2)],[0,-4/(vm*(1+np.cos(wake_skew))),0],[15*np.pi/(64*vt)*np.tan(wake_skew/2),0,-4*np.cos(wake_skew)/(vm*(1+np.cos(wake_skew)))]])
        lam_0,lam_1c,lam_1s = np.dot(L,[CT,CMX,CMY])
        lam = lam_0 + lam_1c*r*np.expand_dims(np.cos(psi),axis = 1)+ lam_1s*r*np.expand_dims(np.sin(psi),axis = 1)

        return lam

    def aeroParams(AoA):
        '''
        This function returns the lift and drag coefficients corresponding to a radial and azimuthal distribution of the
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

        # #   assume that the airfoil is symmetric and therefore the CL can be estimated by the product of the
        # # lift-curve slope and the angle of attack
        CL = XsecPolarExp['Lift Slope'] * AoA
        # #   CD is assumed to be 10% of CL
        CD = 0.1 * CL

        #   reruns the indices of stalled blade sections
        azInd, rInd = np.where(AoA > XsecPolarExp['alphaMax'])
        #   sets the CD of these sections equal to the CD @ CLmax
        CD[azInd, rInd] = XsecPolarExp['CdMax'][azInd, rInd]
        # CL[azInd, rInd] = XsecPolarExp['ClMin'][azInd, rInd]
        #   linearly interpolates CL between CLmin and CL
        CL[azInd, rInd] = XsecPolarExp['ClMax'][azInd, rInd]+(AoA[azInd, rInd]-XsecPolarExp['alphaMax'][azInd, rInd])*(XsecPolarExp['ClMin'][azInd, rInd]-XsecPolarExp['ClMax'][azInd, rInd])/(XsecPolarExp['Alpha0'][azInd, rInd]+2*np.pi-XsecPolarExp['alphaMax'][azInd, rInd])
        return CL, CD

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
            ind = np.squeeze(np.where(np.array(polarInd) == k))

            if np.any(AoA[ind]>v['alphaMax']) or np.any(AoA[ind] < v['Alpha0']):

                ind2 = np.where((AoA[ind]< v['Alpha0']) | (AoA[ind]>v['alphaMax']))
                dCL[ind2] = np.interp(AoA[ind2]%(2*np.pi) ,xp = [v['alphaMax'],2*np.pi],fp = [v['ClMax'],v['ClMin']])

                if np.any(AoA[ind2]%(2*np.pi)>-v['alphaMax']%(2*np.pi)):
                    ind3 = np.squeeze(np.where(AoA[ind2]%(2*np.pi)>-v['alphaMax']%(2*np.pi)))
                    dCD[ind3] = np.interp(AoA[ind2][ind3]%(2*np.pi) ,xp = [-v['alphaMax']%(2*np.pi),2*np.pi],fp = [v['CdMax'],v['CdMin']])
                    ind2 = np.delete(ind2, ind3)
                else:
                    dCD[ind2] = v['CdMax']

                ind = np.delete(ind, ind2)

            dCL[ind] = np.interp(AoA[ind],xp = v['Polar'][:,0],fp =v['Polar'][:,1])
            # dCL[ind] = 2*np.pi*AoA[ind]
            dCD[ind] = np.interp(AoA[ind],xp = v['Polar'][:,0],fp =v['Polar'][:,2])

        return dCL, dCD

    def beddoes_fwake(wake_age):

        psi_w = np.arange(0,360*wake_age)*np.pi/180
        psi_diff = np.expand_dims(psi_b,axis = 1)-psi_w

        x_v =np.cos(psi_diff)+mu*psi_w
        y_v =np.sin(psi_diff)
        
        # rigid wake:
        # z_v = -mu*np.tan(wake_skew)*np.ones(np.shape(x_v))

        z_v = np.zeros(np.shape(x_v))

        if np.any(x_v < -np.cos(psi_diff)):
            z_v[x_v < -np.cos(psi_diff)] = (mu_z*psi_w-lam_const*(1+E*(x_v -0.5*mu*psi_w-abs(y_v**3)))*psi_w)[x_v < -np.cos(psi_diff)]


        if np.any(0 < np.cos(psi_diff)):
            z_v[0 < np.cos(psi_diff)] = (mu_z*psi_w-2*lam_const*(1-E*abs(y_v**3))*psi_w)[0 < np.cos(psi_diff)]

        if np.any(z_v == 0):
            z_v[((0 > np.cos(psi_diff)) & (x_v > -np.cos(psi_diff)))] = (mu_z*psi_w - 2*lam_const*x_v/mu*(1-E*abs(y_v**3)))[((0 > np.cos(psi_diff)) & (x_v > -np.cos(psi_diff)))]
    
        return np.array((x_v,y_v,z_v))
        
    def beddoes_fwake_2(wake_age):

        psi_w = np.arange(0,360*wake_age)*np.pi/180
        psi_diff = np.expand_dims(psi_b,axis = 1)-psi_w

        x_v = np.cos(psi_diff)+mu*psi_w
        y_v = np.sin(psi_diff)
        x_e =  np.array([np.repeat(np.expand_dims(np.roll(np.cos(psi_b),i+1)[:-1],axis = -1),wake_age,axis = -1).flatten(order = 'F') for i in range(len(psi_b))])
        # x_vr = np.expand_dims(lift_line_exp[0,-1],axis = -1)
        # y_vr = np.expand_dims(lift_line_exp[1,-1],axis = -1)


        z_v = np.zeros(np.shape(x_v))

        z_v_x_e = -lam_const/mu*((k0-kx*np.abs(y_v)**3+ky*y_v)*(x_e - np.expand_dims(np.cos(psi_b),axis = -1)) +kx/2*(x_e**2 - np.expand_dims(np.cos(psi_b),axis = -1)**2))-mu_z/mu*(x_e - np.expand_dims(np.cos(psi_b),axis = -1))

        ind = np.abs(x_v) >= np.abs(x_e)
        for i in range(len(psi_b)):
            ind[i,np.squeeze(np.where(ind[i]))[1]:] = True 

        if np.any(ind != 1):
            z_v[ind != 1] = (-lam_const/mu*((k0-kx*np.abs(y_v)**3+ky*y_v)*(x_v - np.expand_dims(np.cos(psi_b),axis = -1)) +kx/2*(x_v**2 - np.expand_dims(np.cos(psi_b),axis = -1)**2))-mu_z/mu*(x_v - np.expand_dims(np.cos(psi_b),axis = -1)))[ind != 1]


        if np.any(ind):
            z_v[ind] = (z_v_x_e-2*lam_const/mu*(k0-kx*np.abs(y_v)**3+ky*y_v)*(x_v-x_e) -mu_z/mu*(x_v-x_e))[ind]

        # if np.any(z_v == 0):
        #     z_v[((0 > np.cos(psi_diff)) & (x_v > -np.cos(psi_diff)))] = (mu_z*psi_w - 2*lam_const*x_v/mu*(1-E*abs(y_v**3)))[((0 > np.cos(psi_diff)) & (x_v > -np.cos(psi_diff)))]
        # import matplotlib.pyplot as plt
        # for i in range(int(360/20)):
        #     plt.plot(z_v[i*10,:360])
        # plt.show()

        return np.array((x_v,y_v,z_v))
        # out = np.array([wake_coord(Nb_ind) for Nb_ind in range(Nb)])

    def beddoes_fwake_3(wake_age):

        psi_w = np.arange(0,360*wake_age)*np.pi/180
        psi_diff = np.expand_dims(psi_b,axis = 1)-psi_w

        x_v_2 = np.sin(psi_diff)+mu*psi_w
        y_v = np.sin(psi_diff)

        ind = 150
        x_e =  np.array([np.repeat(np.expand_dims(np.roll(np.cos(psi_b),i+1)[:-1],axis = -1),wake_age,axis = -1).flatten(order = 'F') for i in range(len(psi_b))])
        x_vr = np.expand_dims(lift_line_exp[0,-1],axis = -1)
        y_vr = np.expand_dims(lift_line_exp[1,-1],axis = -1)


        z_v = np.zeros(np.shape(x_v))

        z_v_x_e = -lam_const/mu*((k0-kx*np.abs(y_v)**3+ky*y_v)*(x_e - np.expand_dims(np.cos(psi_b),axis = -1)) +kx/2*(x_e**2 - np.expand_dims(np.cos(psi_b),axis = -1)**2))-mu_z/mu*(x_e - np.expand_dims(np.cos(psi_b),axis = -1))

        ind = np.abs(x_v) >= np.abs(x_e)
        for i in range(len(psi_b)):
            ind[i,np.squeeze(np.where(ind[i]))[1]:] = True 

        if np.any(ind != 1):
            z_v[ind != 1] = (-lam_const/mu*((k0-kx*np.abs(y_v)**3+ky*y_v)*(x_v - np.expand_dims(np.cos(psi_b),axis = -1)) +kx/2*(x_v**2 - np.expand_dims(np.cos(psi_b),axis = -1)**2))-mu_z/mu*(x_v - np.expand_dims(np.cos(psi_b),axis = -1)))[ind != 1]


        if np.any(ind):
            z_v[ind] = (z_v_x_e-2*lam_const/mu*(k0-kx*np.abs(y_v)**3+ky*y_v)*(x_v-x_e) -mu_z/mu*(x_v-x_e))[ind]

        # if np.any(z_v == 0):
        #     z_v[((0 > np.cos(psi_diff)) & (x_v > -np.cos(psi_diff)))] = (mu_z*psi_w - 2*lam_const*x_v/mu*(1-E*abs(y_v**3)))[((0 > np.cos(psi_diff)) & (x_v > -np.cos(psi_diff)))]
        # import matplotlib.pyplot as plt
        # for i in range(int(360/20)):
        #     plt.plot(z_v[i*10,:360])
        # plt.show()

        return np.array((x_v,y_v,z_v))

    def beddoes_nwake(wake_age):

        psi_w = np.arange(0,360*wake_age)*np.pi/180
        psi_diff = np.expand_dims(psi_b,axis = 1)-psi_w

        x_v =np.expand_dims(np.expand_dims(geomParams['r'],axis = -1),axis = -1)*np.cos(psi_diff)+mu*psi_w
        y_v = np.expand_dims(np.expand_dims(geomParams['r'],axis = -1),axis = -1)*np.sin(psi_diff)
        z_v = np.zeros(np.shape(x_v))


        if np.any(x_v <= -np.cos(psi_diff)):
            z_v[x_v < -np.cos(psi_diff)] = (mu_z*psi_w-lam_const*(1+E*(x_v -0.5*mu*psi_w-abs(y_v**3)))*psi_w)[x_v < -np.cos(psi_diff)]

        if np.any(0 < np.cos(psi_diff)):
            z_v[np.zeros(np.shape(x_v)) < np.cos(psi_diff)] = (mu_z*psi_w-2*lam_const*(1-E*abs(y_v**3))*psi_w)[np.zeros(np.shape(x_v)) < np.cos(psi_diff)]

        if np.any(z_v == 0):
            # z_v[z_v == 0] = (mu_z*psi_w - 2*lam_const*x_v/mu*(1-E*abs(y_v)**3))[z_v == 0]

            z_v[((0 > np.cos(psi_diff)) & (x_v > -np.cos(psi_diff)))] = (mu_z*psi_w - 2*lam_const*x_v/mu*(1-E*abs(y_v**3)))[((0 > np.cos(psi_diff)) & (x_v > -np.cos(psi_diff)))]
    
        
        # out = np.array([wake_coord(Nb_ind) for Nb_ind in range(Nb)])

        return np.array((x_v,y_v,z_v))
   
    def beddoes_nwake_2(wake_age):

        psi_w = np.arange(0,360*wake_age)*np.pi/180
        psi_diff = np.expand_dims(psi_b,axis = 1)-psi_w

        x_v =np.expand_dims(np.expand_dims(geomParams['r'],axis = -1),axis = -1)*np.cos(psi_diff)+mu*psi_w
        y_v = np.expand_dims(np.expand_dims(geomParams['r'],axis = -1),axis = -1)*np.sin(psi_diff)
        z_v = np.zeros(np.shape(x_v))

        x_e =  np.array([np.roll(np.cos(psi_b),i+1)[:int(360*wake_age)] for i in range(len(psi_b))])
        z_v_x_e = -lam_const/mu*((k0-kx*np.abs(y_v)**3+ky*y_v)*(x_e - np.expand_dims(np.cos(psi_b),axis = -1)) +kx/2*(x_e**2 - np.expand_dims(np.cos(psi_b),axis = -1)**2))-mu_z/mu*(x_e - np.expand_dims(np.cos(psi_b),axis = -1))


        ind = np.abs(x_v) >= np.abs(x_e)
        # for i in range(len(psi_b)):
        #     ind[:,i,np.squeeze(np.where(ind[i]))[1]:] = True 


        if np.any(ind != 1):
            z_v[ind != 1] = (-lam_const/mu*((k0-kx*np.abs(y_v)**3+ky*y_v)*(x_v - np.expand_dims(np.cos(psi_b),axis = -1)) +kx/2*(x_v**2 - np.expand_dims(np.cos(psi_b),axis = -1)**2))-mu_z/mu*(x_v - np.expand_dims(np.cos(psi_b),axis = -1)))[ind != 1]

        if np.any(ind):
            z_v[ind] = (z_v_x_e-2*lam_const/mu*(k0-kx*np.abs(y_v)**3+ky*y_v)*(x_v-x_e) -mu_z/mu*(x_v-x_e))[ind]
    
        
        # out = np.array([wake_coord(Nb_ind) for Nb_ind in range(Nb)])

        return np.array((x_v,y_v,z_v))

    def unsteady_aero(AoA):
    
        # da_ds = np.diff(AoA,axis = 0)
        da_ds = [(-1/12*AoA[(i+2)%360]+2/3*AoA[(i+1)%360]-2/3*AoA[(i-1)%360]+1/12*AoA[(i-2)%360]) for i in range(len(psi)-1)]
        # da_ds = [3/2*AoA[i]-2*AoA[i-1]+1/2*AoA[i-2] for i in range(len(psi)-1)]

        X_err = 1
        X = np.empty(np.shape(s))
        Y = np.empty(np.shape(s))

        X[0] = A[0]*da_ds[0]
        Y[0]= A[1]*da_ds[0]

        while np.any(X_err > 1e-5) or np.any(Y_err > 1e-5):
            for i in range(len(psi)-1):
                X[i+1] = X[i]*np.exp(-B[0]*ds[i])+A[0]*da_ds[i]*np.exp(-B[0]*ds[i]/2)
                Y[i+1] = Y[i]*np.exp(-B[1]*ds[i])+A[1]*da_ds[i]*np.exp(-B[1]*ds[i]/2)
            X_err = X[-1]-X[0]
            Y_err= Y[-1]-Y[0]
            X[0] = X[-1]
            Y[0] = Y[-1]

        AoA_eff = AoA-X-Y

        return AoA_eff


    #%%
    omega = omega/60*2*np.pi
    rho = UserIn['rho']
    Nb = UserIn['Nb']
    R = geomParams['R']
    r = geomParams['r']
    r = (r[1:]+r[:-1])/2
    twist = (geomParams['twistDist'][:-1]+geomParams['twistDist'][1:])/2
    r_TE = np.linalg.norm(geomParams['TENodes'],axis =1)
    r_TE = (r_TE[:-1]+r_TE[1:])/2


    solDist = (geomParams['solDist'][:-1]+geomParams['solDist'][1:])/2
    chordDist = (geomParams['chordDist'][:-1]+geomParams['chordDist'][1:])/2

    XsecLocation = UserIn['XsecLocation']
    th_tw = geomParams['twistDist'][-1]-geomParams['twistDist'][0]
    gamma = 6

    alphaShaft = alphaShaft*(np.pi/180)
    thFP = np.arctan(Vz / Vx)
    alphaInit = alphaShaft+thFP
    U = np.linalg.norm((Vx,Vz))

    mu = U*np.cos(alphaInit)/(omega*R)
    mu_z = U*np.sin(alphaInit)/(omega*R)

    beta_init = np.array([1,1,1])*np.pi/180
    phiRes = 361
    psi = np.linspace(0,2*np.pi,phiRes)
    a = np.ones((len(r)))*XsecPolar[list(XsecPolar.keys())[0]]['Lift Slope']
    th0 = UserIn['thetaInit']*np.pi/180
    CT = W/(rho*np.pi*R**2*(omega*R)**2)
    lam_init =  constant_inflow(np.sqrt(CT/2),mu,CT)

    if UserIn['UnsteadyAero']:
                #%%
        s = np.expand_dims(psi,axis = 1)*r*R/(chordDist/2)
        ds = np.diff(s,axis = 0)
        # Beddoes indicial linear response approximation
        # A = np.array([0.3,0.7])
        # B = np.array([.1,.53])
        #  W.P. Jones indicial response approximation 
        A = np.array([0.156,0.335])
        B = np.array([.041,.335])
        # da_ds = np.gradient(AoA,axis = 0, edge_order = 2)

        # da_ds = [(3*AoA[i]-4*AoA[i-1]+AoA[i-2])/(ds[i]**2) for i in range(len(psi)-1)]



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


    if UserIn['pwake']:

        lam_const = constant_inflow(np.sqrt(CT/2),mu,CT)
        # wake_skew = np.arctan((mu_z-lam_const)/mu)
        wake_skew = np.arctan(mu/(mu_z+lam_const))
        E = wake_skew/2
        # kx = 4/3*(1-np.cos(wake_skew)-1.8*mu**2)/np.sin(wake_skew)
        # ky = -2*mu
        # ky = 0
        # kx = E
        # k0 = 1+8*kx/(15*np.pi)
        psi_b = np.arange(0,361)*np.pi/180


        alphaShaft_dcm = np.array([[np.cos(alphaShaft),0,np.sin(alphaShaft)],[0,1,0],[-np.sin(alphaShaft),0,np.cos(alphaShaft)]])
        rotor_frm_dcm = np.array([[np.sin(psi),np.cos(psi),np.zeros(len(psi))],[-np.cos(psi),np.sin(psi),np.zeros(len(psi))],[np.zeros(len(psi)),np.zeros(len(psi)),np.ones(len(psi))]])
        tot_dcm = np.matmul(alphaShaft_dcm,rotor_frm_dcm)
        lift_line = geomParams['LENodes'] - (geomParams['LENodes']-geomParams['TENodes'])*0.25
        control_pnt = geomParams['LENodes'] - (geomParams['LENodes']-geomParams['TENodes'])*0.75
        control_pnt = (control_pnt[1:]+control_pnt[:-1])/2
        lift_line_exp = np.matmul(lift_line/R,rotor_frm_dcm)
        # lift_line_exp = np.array([np.matmul(lift_line_exp[:,:,i].transpose(),alphaShaft_dcm) for i in range(len(psi))]).transpose(2,1,0)
        control_pnt_exp = np.matmul(control_pnt,rotor_frm_dcm)/R
        # lift_line_exp = np.matmul(lift_line_exp,alphaShaft_dcm)



        surfNodes = geomParams['surfNodes'].reshape(geomParams['pntsPerXsec'],geomParams['nXsecs'],3,order = 'F')
        
        # y = np.array((np.expand_dims(surfNodes[:, :, 1], axis=-1) * np.cos(psi) + np.expand_dims(surfNodes[:, :, 0],axis=-1) * np.sin(psi),
        #             np.expand_dims(-surfNodes[:, :, 0], axis=-1) * np.cos(psi) + np.expand_dims(surfNodes[:, :, 1], axis=-1) * np.sin(psi),
        #             np.expand_dims(surfNodes[:, :, -1], axis=-1) * np.ones(len(psi))))/R

        # control_pnt = geomParams['LENodes'] - (geomParams['LENodes'] - geomParams['TENodes']) * 0.75
        # control_pnt = (control_pnt[1:]+control_pnt[:-1])/2
        # control_pnt_exp = np.array((np.expand_dims(control_pnt[:, 1], axis=-1) * np.cos(psi) + np.expand_dims(control_pnt[:, 0],axis=-1) * np.sin(psi),
        #         np.expand_dims(-control_pnt[:, 0], axis=-1) * np.cos(psi) + np.expand_dims(control_pnt[:, 1], axis=-1) * np.sin(psi),
        #         np.expand_dims(control_pnt[:, -1], axis=-1) * np.ones(len(psi))))/R

        # # lift_line = geomParams['LENodes'] - (geomParams['LENodes'] - geomParams['TENodes']) * 0.25
        # lift_line_exp = np.array((np.expand_dims(lift_line[:, 1], axis=-1) * np.cos(psi) + np.expand_dims(lift_line[:, 0],axis=-1) * np.sin(psi),
        #         np.expand_dims(-lift_line[:, 0], axis=-1) * np.cos(psi) + np.expand_dims(lift_line[:, 1], axis=-1) * np.sin(psi),
        #         np.expand_dims(lift_line[:, -1], axis=-1) * np.ones(len(psi))))/R
        
        control_pnt_norm = np.zeros(np.shape(control_pnt.transpose()))
        control_pnt_norm[-1] = 1

        wake_age_fw = 10
        fwake = beddoes_fwake(wake_age = wake_age_fw)
        wake_age_nw = 0.1
        nwake = beddoes_nwake(wake_age = wake_age_nw)


    #%% Near-wake influence coefficients
        r_nw = np.expand_dims(np.expand_dims(control_pnt_exp,axis  = 2),axis = -1) - np.expand_dims(nwake,axis = 1)
        l_v_nw = np.expand_dims(np.diff(nwake,axis = -1),axis = 1)
        r_nw_unit = r_nw/np.expand_dims(np.linalg.norm(r_nw,axis =0),axis = 0)
        # l_v_nw_unit = np.expand_dims(l_v_nw/np.expand_dims(np.linalg.norm(l_v_nw,axis =0),axis = 0),axis = 1)
        x_p_nw = np.cross(r_nw[:,:,:,:,:-1], r_nw[:,:,:,:,1:],axis = 0)
        I_nw = (4*np.pi)**-1*np.sum(l_v_nw*(r_nw_unit[:,:,:,:,:-1] - r_nw_unit[:,:,:,:,1:]),axis = 0)*x_p_nw/np.expand_dims(np.linalg.norm(x_p_nw,axis = 0)**2,axis = 0)
        I_nw_tot = np.sum(I_nw,axis = -1)


        # #%%
        # r_nw_1 = np.flip(r_nw[:,0,0,0,:],axis = 1)
        # r_nw_2 = r_nw[:,0,1,0,:]
        # r_nw_1_unit = r_nw_1/np.expand_dims(np.linalg.norm(r_nw_1,axis = 0),axis = 0)
        # r_nw_2_unit = r_nw_2/np.expand_dims(np.linalg.norm(r_nw_2,axis = 0),axis = 0)
        # l_v_nw_1 = np.diff(np.flip(nwake[:,0,0],axis = -1)) 
        # l_v_nw_2 = np.diff(nwake[:,1,0],axis = 1)         
        # x_p_nw_1 = np.cross(r_nw_1[:,:-1], r_nw_1[:,1:],axis = 0)
        # x_p_nw_2 = np.cross(r_nw_2[:,:-1], r_nw_2[:,1:],axis = 0)
        # I_nw_1 = (4*np.pi)**-1*np.sum(l_v_nw_1*(r_nw_1_unit[:,:-1] - r_nw_1_unit[:,1:]),axis = 0)*x_p_nw_1/np.expand_dims(np.linalg.norm(x_p_nw_1,axis = 0)**2,axis = 0)
        # I_nw_2 = (4*np.pi)**-1*np.sum(l_v_nw_2*(r_nw_2_unit[:,:-1] - r_nw_2_unit[:,1:]),axis = 0)*x_p_nw_2/np.expand_dims(np.linalg.norm(x_p_nw_2,axis = 0)**2,axis = 0)




                # h_nw = np.linalg.norm(np.cross(r_nw[:,:,:,:,:-1],l_v_nw_unit,axis = 0),axis = 0)
        # I_nw = (4*np.pi*h_nw)**-1*np.sum(l_v_nw_unit*(r_nw_unit[:,:,:,:,:-1] - r_nw_unit[:,:,:,:,1:]),axis = 0)*np.cross(r_nw[:,:,:,:,:-1], r_nw[:,:,:,:,1:],axis = 0)/np.expand_dims(np.linalg.norm(np.cross(r_nw[:,:,:,:,:-1], r_nw[:,:,:,:,1:],axis = 0),axis = 0),axis = 0)
        # cross_nw = np.cross(l_v_nw, r_nw[:,:,:,:,:-1],axis = 0)
        # I_nw_3 = (4*np.pi)**-1*np.sum(l_v_nw*(r_nw_unit[:,:,:,:,:-1] - r_nw_unit[:,:,:,:,1:]),axis = 0)*cross_nw/np.expand_dims(np.linalg.norm(cross_nw,axis = 0)**2,axis = 0)

        #np.linalg.norm(np.cross(r_nw[:,:,:,:,:-1],l_v_nw_unit,axis = 0),axis = 0)-np.linalg.norm(np.cross(r_nw[:,:,:,:,:-1],r_nw[:,:,:,:,1:],axis = 0)/np.linalg.norm(l_v_nw,axis = 0),axis = 0)
        # np.linalg.norm(np.cross(l_v_nw,r_nw[:,:,:,:,:-1],axis = 0),axis = 0)/np.linalg.norm(l_v_nw,axis = 0)-np.linalg.norm(np.cross(r_nw[:,:,:,:,:-1],r_nw[:,:,:,:,1:],axis = 0)/np.linalg.norm(l_v_nw,axis = 0),axis = 0)
    # I_b[:,0,0,0]+I_nw_tot[:,0,0,0]+I_nw_tot[:,0,1,0]

        
        
    #%% Bound circulation influence coefficients

        # r_b = np.expand_dims(control_pnt_exp,axis = 2)- np.expand_dims(lift_line_exp,axis = 1)
        r_b = np.expand_dims(control_pnt_exp,axis = -1)-np.expand_dims(lift_line_exp.transpose(0,2,1),axis = 1)
        r_b_unit = r_b/np.expand_dims(np.linalg.norm(r_b,axis =0),axis = 0)
        l_v_b = np.expand_dims(np.diff(lift_line_exp.transpose(0,2,1),axis = -1),axis = 1)
        x_p_b = np.cross(r_b[:,:,:,:-1], r_b[:,:,:,1:],axis = 0)
        I_b = (4*np.pi)**-1*np.sum(l_v_b*(r_b_unit[:,:,:,:-1] - r_b_unit[:,:,:,1:]),axis = 0)*x_p_b/np.expand_dims(np.linalg.norm(x_p_b,axis = 0)**2,axis = 0)

        # l_v_b_unit = np.expand_dims(l_v_b/np.expand_dims(np.linalg.norm(l_v_b,axis =0),axis = 0),axis = 1)
        # h_b = np.linalg.norm(np.cross(r_b[:,:,:,:-1],l_v_b_unit,axis = 0),axis = 0)
        # I_b_2 = (4*np.pi*h_b)**-1*np.sum(l_v_b_unit*(r_b_unit[:,:,:,:-1] - r_b_unit[:,:,:,1:]),axis = 0)*np.cross(r_b_unit[:,:,:,:-1], r_b_unit[:,:,:,1:],axis = 0)
        # I_b = (4*np.pi)**-1*np.sum(l_v_b*(r_b_unit[:,:,:,:-1]-r_b_unit[:,:,:,1:]),axis = 0)*np.sum(cross_b/np.expand_dims(np.linalg.norm(cross_b,axis = 0),axis = 0)**2*np.expand_dims(np.expand_dims(control_pnt_norm,axis =-1),axis = -1),axis = 0)
        
        # I_b_tot = np.sum(I_b,axis = -1)

        # I_tot = np.expand_dims(I_b_tot,axis = 1) +I_nw_tot
        # I_tot_2 = np.sum(I_tot,axis = 1)
        # I_tot = I_nw_tot+I_b.transpose(0,1,-1,2)

        I_tot = np.array([(I_b[:,:,:,j]-I_nw_tot[:,:,j,:]+I_nw_tot[:,:,j+1,:]) for j in range(len(geomParams['r'])-1)]).transpose(1,2,0,-1)

        I_tot_inv = np.array([np.linalg.inv(I) for I in I_tot[-1].transpose(-1,0,1)]).transpose(1,2,0)

    # #%% far-wake influence coefficients
        t_wake = (np.arange(wake_age_nw*360,wake_age_fw*360)*np.pi/180)/omega
    #     #  kinematic viscosity (m^2/s)
        nu = 14.88e-6
        # Reynolds number of tip vortex
        Re_v = (2*omega*np.mean(chordDist)*R)/nu*CT/geomParams['solidity']
        #   radius of vortex core (Squire-1965)
        r_c = np.sqrt(4*1.256*nu*(1+6.5e-5*Re_v)*t_wake)
        # r_c = np.sqrt(4*1.256*nu*t_wake)

        I_fw = np.empty((len(psi),3,len(r),Nb,len(t_wake)))

        for ind in range(len(psi)):
            r_fw = np.expand_dims(np.expand_dims(control_pnt_exp[:,:,ind],axis = -1),axis = -1)-np.expand_dims(np.roll(fwake[:,:-1,int(wake_age_nw*360)-1:],-ind,axis = 1)[:,::int(360/Nb)],axis = 1)
            # r_fw = np.expand_dims(np.expand_dims(np.roll(control_pnt_exp[:,:,ind],axis = -1),axis = -1)-np.expand_dims(np.expand_dims(np.roll(fwake[:,:-1,int(wake_age_nw*360)-1:],-ind,axis = 1)[:,::int(360/Nb)],axis = 1),axis = 1)
            # r_fw = np.expand_dims(np.expand_dims(control_pnt_exp[:,:,ind-int(np.floor(ind/90)*360/Nb)::int(360/Nb)],axis = -1),axis = -1)-np.expand_dims(fwake[:,:,int(wake_age_nw*360):],axis = 1)
            # r_fw = np.expand_dims(control_pnt_exp,axis  = -1) - np.expand_dims(fwake[:,:,int(wake_age_nw*360)-1:],axis = 1)
            r_fw_unit = r_fw/np.linalg.norm(r_fw,axis =0)
            l_v_fw = np.diff(r_fw,axis = -1)
            l_v_fw_unit = l_v_fw/np.linalg.norm(l_v_fw,axis =0)
            x_p_fw = np.cross(r_fw[:,:,:,:-1],r_fw[:,:,:,1:],axis = 0)
            h_fw = np.linalg.norm(x_p_fw,axis = 0)/np.linalg.norm(l_v_fw,axis = 0)
            # n = 2
            I_fw[ind] = (4*np.pi*r_c)**-1*(h_fw/r_c)/((1+(h_fw/r_c)**4)**0.5)*np.sum(l_v_fw_unit*(r_fw_unit[:,:,:,:-1]-r_fw_unit[:,:,:,1:]),axis = 0)*x_p_fw/np.linalg.norm(x_p_fw,axis = 0)
        I_fw_tot = np.sum(np.sum(I_fw,axis= -1),axis = -1)
        
        # I_fw = (4*np.pi)**-1*(h_fw/(np.expand_dims(np.expand_dims(r_c,axis =0),axis = 0)**4+h_fw**4)**0.5)*np.sum(l_v_fw_unit*(r_fw_unit[:,:,:,:-1]-r_fw_unit[:,:,:,1:]),axis = 0)*np.cross(r_fw_unit[:,:,:,:-1],r_fw_unit[:,:,:,1:],axis = 0)

        # I_fw_tot = np.array([np.sum(np.sum(np.roll(I_fw[:,:,:-1],-ind,axis = 2)[:,:,::int(360/Nb)],axis = -1),axis = -1) for ind in range(len(psi))])
        # I_fw_tot = np.array([np.sum(np.sum(I_fw[-1,:,ind::int(360/Nb)][:,:,:Nb],axis = -1),axis = -1) for ind in range(len(psi))])


#     up_err = 1
#     i = 0
#     while np.any(up_err > 0.00001):
#         AoA = twist-np.arctan2(up,ut)
#         V = np.sqrt(ut**2+up**2)*omega*R
#         gamma_v_2 = V*AoA/I_tot_2.transpose()
#         up_temp = (I_fw_tot*np.amax(gamma_v_2,axis = 1)).transpose()
#         up_err =  np.abs((up_temp - up) / up_temp)
#         up = up_temp
#         i =+1
#             gamma = np.dot(I_tot_inv[:,:,90],V[90])*AoA[90]
#     gamma_v_tot = np.sum(gamma_v,axis = 0)
#     up = I_fw_tot*np.max(gamma_v_tot)
# #%%
    if UserIn['trim']==1:
        trimTargs = W
        trim_sol = least_squares(fixed_pitch_residuals, omega, method = 'lm',diff_step = 0.5)
        omega = trim_sol.x
        th = np.array([th0,0 ,0 ])
        T,CT,dCT,lam,ut,up,CL,CD,AoA,mu = fixed_pitch_trim(omega)
        theta_expanded = np.empty((np.shape(AoA)))*th0

    elif UserIn['trim'] == 2:
        trimTargs = W/(rho*np.pi*R**2*(omega*R)**2)
        lamTPP_init =  inflowModSelect(1, mu*np.tan(alphaInit), mu, trimTargs,alphaInit)
        trim_sol = least_squares(variable_pitch_residuals, th0, args=[mu, lamTPP_init], method='lm')
        th = np.array([np.squeeze(trim_sol.x),0 ,0 ])
        CT,dCT,Mx,My,lam,theta_expanded,ut,up,CL,CD,AoA = variable_pitch_trim(th,mu, lamTPP_init)

    else:
        trimTargs = [W/(rho*np.pi*R**2*(omega*R)**2),0*W/(rho*np.pi*R**2*(omega*R)**2*R),0*W/(rho*np.pi*R**2*(omega*R)**2*R)]
        th = np.array([th0,0*np.pi/180,0*np.pi/180])
        trim_sol = least_squares(variable_pitch_residuals, th ,args = [mu, lam_init],method = 'lm')
        th = trim_sol.x
        # th[2] = th[2]+5*np.pi/180
        # CT,dCT,Mx,My,lam,theta_expanded,ut,up,CL,CD,AoA,beta0,beta1c,beta1s,beta_exp,alphaInit = variable_pitch_trim(th,mu, lamTPP_init,alphaInit)
        CT,dCT,CMX,CMY,lam,th_exp,ut,up,CL,AoA,AoA_eff = variable_pitch_trim_2(th,mu, lam_init)




#%



    # # gamma_v = np.multiply(I_tot_inv*np.expand_dims(V.transpose(),axis = 0),np.expand_dims(AoA.transpose(),axis = 0))
    # gamma = np.dot(I_tot_inv[:,:,90],V[90])*AoA[90]
    # gamma_v_tot = np.sum(gamma_v,axis = 0)
    # up = I_fw_tot*np.max(gamma_v_tot)

    # #%%

#%%

    # import matplotlib.pyplot as plt
    # ind = 210 
    # fig = plt.figure(figsize=(6.4, 4.5))
    # fig.subplots_adjust(top=1.3, bottom=-.3)
    # ax = fig.add_subplot(111, projection='3d')
    # ax.set_box_aspect((1, 1, 1))
    # for Nb_ind in range(Nb):
    #     ax.plot3D(fwake[0,int(ind+Nb_ind*360/Nb)%360,35:],fwake[1,int(ind+Nb_ind*360/Nb)%360,35:],fwake[2,int(ind+Nb_ind*360/Nb)%360,35:])
    #     # ax.plot_surface(nwake[0,:,int(ind+Nb_ind*360/Nb)%360], nwake[1,:,int(ind+Nb_ind*360/Nb)%360], nwake[2,:,int(ind+Nb_ind*360/Nb)%360])
    #     # ax.plot_wireframe(y[0,:,:,int(ind+360/Nb*Nb_ind)], y[1,:,:,int(ind+360/Nb*Nb_ind)], y[2,:,:,int(ind+360/Nb*Nb_ind)])
    # # ax.plot_wireframe(lift_line_exp[0,:,::180]*np.cos(-alphaShaft)+lift_line_exp[2,:,::180]*np.sin(-alphaShaft), lift_line_exp[1,:,::180], -lift_line_exp[0,:,::180]*np.sin(-alphaShaft)+lift_line_exp[2,:,::180]*np.cos(-alphaShaft))
    # # ax.plot_wireframe(lift_line_exp[0], lift_line_exp[1], lift_line_exp[2])

    # ax.set(xlabel ='x/R',ylabel = 'y/R',zlabel = 'z/R')
    # ax.view_init(0,90)
    # ax.axis([-1.5,1.5,-1.5,1.5])
    # ax.set_zlim([-1.5,1.5])
    # plt.show()

    # #%%
    # AoA[25:45] = AoA[25:45] + 2*np.pi/180

    #%%
    s = np.expand_dims(psi,axis = 1)*r*R/(chordDist/2)
    ds = np.diff(s,axis = 0)
    # Beddoes indicial linear response approximation
    # A = np.array([0.3,0.7])
    # B = np.array([.1,.53])
    #  W.P. Jones indicial response approximation 
    A = np.array([0.156,0.335])
    B = np.array([.041,.335])
    # da_ds = np.gradient(AoA,axis = 0, edge_order = 2)
    da_ds = np.diff(AoA,axis = 0)
    # da_ds = [(3*AoA[i]-4*AoA[i-1]+AoA[i-2])/(ds[i]**2) for i in range(len(psi)-1)]
    # da_ds = [(3/2*AoA[i]-2*AoA[i-1]+1/2*AoA[i-2]) for i in range(len(psi)-1)]
    # da_ds = [(1/2*AoA[i+1]-1/2*AoA[i-1]) for i in range(len(psi)-1)]
    # da_ds = [(-1/12*AoA[(i+2)%360]+2/3*AoA[(i+1)%360]-2/3*AoA[(i-1)%360]+1/12*AoA[(i-2)%360]) for i in range(len(psi)-1)]

    
    # X_err = 1
    # X = np.empty(np.shape(s))
    # Y = np.empty(np.shape(s))

    # X[0] = A[0]*da_ds[0]
    # Y[0]= A[1]*da_ds[0]

    # while np.any(X_err > 1e-5) or np.any(Y_err > 1e-5):
    #     # for i in range(len(psi)-1):
    #     #     X[i+1] = X[i]*np.exp(-B[0]*ds[i])+A[0]*da_ds[i]
    #     #     Y[i+1] = Y[i]*np.exp(-B[1]*ds[i])+A[1]*da_ds[i]
    #     # X_err = X[-1]-X[0]
    #     # Y_err= Y[-1]-Y[0]
    #     # X[0] = X[-1]
    #     # Y[0] = Y[-1]

    #     for i in range(len(psi)-1):
    #         X[i+1] = X[i]*np.exp(-B[0]*ds[i])+A[0]*da_ds[i]*np.exp(-B[0]*ds[i]/2)
    #         Y[i+1] = Y[i]*np.exp(-B[1]*ds[i])+A[1]*da_ds[i]*np.exp(-B[1]*ds[i]/2)
    #     X_err = X[-1]-X[0]
    #     Y_err= Y[-1]-Y[0]
    #     X[0] = X[-1]
    #     Y[0] = Y[-1]

    # AoA_eff = AoA-X-Y

    



#%%
    # CD = 0.1*CL
    CD = 0.01-(XsecPolarExp['Lift Slope']*(1-0.95)*XsecPolarExp['Alpha0'])*AoA+(XsecPolarExp['Lift Slope']*(1-0.95))*AoA**2
    UT = ut*omega*R
    UP = up * omega * R
    U = np.sqrt(UT**2+UP**2)

    dT = rho*np.pi*R**2*(omega*R)**2*dCT
    T = 1/(2*np.pi)*np.trapz(np.trapz(dT,r),psi)

    dCQ = 0.5*solDist*r**3*(CL*np.sin(up/ut)+CD*np.cos(up/ut))
    CQ = 1/(2*np.pi)*np.trapz(np.trapz(dCQ,r),psi)
    dQ = rho*np.pi*R**3*(omega*R)**2*dCQ
    Q = 1/(2*np.pi)*np.trapz(np.trapz(dQ,r),psi)
    P = Q * omega

    # resolves loading vectors to vertical and horizontal directions so that a change of base can be applied to the
    # blade geometry account for the pitching motion in the namelist file - 1/18/21
    # dFz =  dT/Nb*np.cos(-theta_expanded)-dQ/(Nb*r*R)*np.sin(-theta_expanded)
    # dFx = dT/Nb*np.sin(-theta_expanded)+dQ/(Nb*r*R)*np.cos(-theta_expanded)
    dFz = dT / Nb
    dFx = dQ / (Nb * r * R)
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
    H = Nb/(2*np.pi)*np.trapz(np.trapz((dFr*np.expand_dims(np.cos(psi),axis = 1)+dFx*np.expand_dims(np.sin(psi),axis = 1)),r),psi)
    #   side force
    Y = Nb/(2*np.pi)*np.trapz(np.trapz((dFr*np.expand_dims(np.sin(psi),axis = 1)-dFx*np.expand_dims(np.cos(psi),axis = 1)),r),psi)
    # #   roll moment
    # Mx = Nb / (2 * np.pi) * np.trapz(np.trapz(geomParams['rdim'] * dFz * np.expand_dims(np.sin(psi), axis=1), r), psi)
    # #   pitch moment
    # My = -Nb / (2 * np.pi) * np.trapz(np.trapz(geomParams['rdim'] * dFz * np.expand_dims(np.cos(psi), axis=1), r), psi)

    # Mx = 1/(2*np.pi)*np.trapz(np.trapz(r*dCT*np.expand_dims(np.sin(psi),axis = 1),r),psi)*rho*(omega*R)**2*np.pi*R**3
    # My = -1/(2*np.pi)*np.trapz(np.trapz(r*dCT*np.expand_dims(np.cos(psi),axis = 1),r),psi)*rho*(omega*R)**2*np.pi*R**3
    # print(f'{Mx2},{My2}')
    # hubLM = [H,Y,Mx,My]

    # if UserIn['pwake']:
    #     loadParams = {'residuals':trim_sol.fun,'phiRes':phiRes,'ClaDist':a,'AoA':AoA,'alpha':alphaInit,'mu':mu,'phi':psi,'th':th,'CT':CT,'T':T,'CQ':CQ,'Q':Q,'P':P,
    #                 'UP':UP,'UT':UT,'U':U,'dFx':dFx,'dFy':dFr,'dFz':dFz,'dCT':dCT,'lam':lam,'omega':omega*60/(2*np.pi),'fwake':fwake,'nwake':nwake,'y':y}

    # else:
    #   assembles a dictionary with the computed parameters that is returned to the user and is referenced in other segments of the program
    loadParams = {'residuals':trim_sol.fun,'phiRes':phiRes,'ClaDist':a,'AoA':AoA,'alpha':alphaInit,'mu':mu,'phi':psi,'th':th,'CT':CT,'T':T,'CQ':CQ,'Q':Q,'P':P,
                'UP':UP,'UT':UT,'U':U,'dFx':dFx,'dFy':dFr,'dFz':dFz,'dCT':dCT,'lam':lam,'omega':omega*60/(2*np.pi)}
    #
    return loadParams

    #
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    # AoA[AoA>np.pi] = AoA[AoA>np.pi]-2*np.pi
    quant = AoA
    # levels = np.linspace(-.75, 1.2, 50)
    levels = np.linspace(np.min(quant),np.max(quant),50)
    dist = ax.contourf(psi, r, quant.transpose(),levels = levels)
    # ax.set_ylim(geomParams['rdim'][0],geomParams['rdim'][-1])
    cbar = fig.colorbar(dist)
    # cbar.ax.set_ylabel(r'$\lambda$')
    cbar.ax.set_ylabel(r'$\alpha \ [\circ]$')
    plt.show()
    

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,1)
    ax.plot(180/np.pi*AoA[:,int(0.75*len(r))])
    ax.plot(180/np.pi*AoA_eff[:,int(0.75*len(r))])
    ax.set_ylabel(r'$\alpha \ [\circ]$')
    ax.set_xlabel(r'$\psi \ [\circ]$')
    ax.set_xlim([0,360])
    ax.grid('on')
    ax.legend(['without unsteady effects','with unsteady effects'])
    plt.show()



    fig, ax = plt.subplots(1,1)
    for i in range(8):
        ax.plot(gamma_b[int(i*(360/8))])
    plt.show()

    N = 3
    for i in range(N):
        print(th[1]*np.cos(2*np.pi/N*i)+th[2]*np.sin(2*np.pi/N*i))
    for i in range(N):
        print(-th[1]*np.sin(2*np.pi/N*i)+th[2]*np.cos(2*np.pi/N*i))