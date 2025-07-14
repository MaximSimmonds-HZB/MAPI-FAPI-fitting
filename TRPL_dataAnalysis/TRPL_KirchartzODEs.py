import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


def TRPL_function3(params, xdata, N_c2, N_v2, kT, G0, Eg2, t_off = False):
    
    """
    Function to solve ODEs using xdata as time axis. t_off will define a time at which laser is off.
    
    Parameters
    ----------
    params : params, lmfit.Parameters like. 
    
    xdata [s]: time axis, array like.
    
    mode: slection of the ODEs to be solved. string: "Decay", "Rise-Decay", "k2".
     
    t_off [s]: time at which the pulse ends. float.
    
    N0 [cm-3, cm-3, cm-3]: initial concnetration for decay. 
        
    Returns
    -------
    n, p, ntr at time = t 
    
    """
    
    krad_c = 10**params['krad_c'] 
    
    Et1 = params['Et1'] 
    Nt1_c = 10**params['Nt1_c']
    b1n_c = 10**params['b1n_c'] 
    b1p_c = 10**params['b1p_c']
    
    G0_1sun = G0*10**params['G0_1sun_c']
    
    ni2 = N_c2*N_v2*np.exp(-Eg2/(kT))

    #Solvig ODE, calculation of n(t), p(t) and ntr(t)
    a = 1
    
    arg = [G0_1sun, krad_c, Et1, Nt1_c, b1n_c, b1p_c, N_c2, N_v2, ni2, Eg2, kT]
    t = xdata
    sol_eq = solve_ivp(One_Shallow_SteadyState3, [0, 0.1], [ni2, 1e18*10**Nt1_c*(np.exp(-(Et1)/(2*kT)))], args = arg, method='LSODA', t_eval = t)
    arg = [0, krad_c, Et1, Nt1_c, b1n_c, b1p_c, N_c2, N_v2, ni2, Eg2, kT]

    sol = solve_ivp(One_Shallow_SteadyState3, [0, 1000e-6], [sol_eq.y[0][-1], sol_eq.y[1][-1]], args = arg, method='LSODA', t_eval = t)
    n = sol.y[0].T
    p = sol.y[0].T + sol.y[1].T
    ntr = sol.y[1].T
    
    y = krad_c*1e-11*1e6*(n*p-ni2) #in m3

    print(y[0]/np.amax(y))
        
    return np.log10(1e6*(y/np.amax(y)))


def One_Shallow_SteadyState3(t, z, G, krad_c, Et1, Nt1_c, b1n_c, b1p_c, N_c2, N_v2, ni2, Eg, kT):
    """
    Function to propagate ODE with initial condition in the dark, constant generation G. 
    Rise of PL during an infinite pulse. 
    Designed to be used by the scipy.optimize.solve_ivp() function.

    Parameters
    ----------
    t [s]: time axis, vector like. 
    
    z [cm3, cm3, cm3]: ODE variables, array like. z = [n, p, ntr1]
    
    G [cm-3.s-1]: Generation rate. Float.
    krad_c [1e11*cm3.s-1]: Coefficient for radiative recombination. Float.
    Et1 [eV]: Trap level energy, referenced from the valence band. Float.
    Nt1_c [1e-18*cm-3]: Coefficient for trap density. Float.
    b1n_c [1e11*cm3.s-1]: Electron capture coefficient. Float.
    b1n_c [1e11*cm3.s-1]: Hole capture coefficient. Float.
        
    Returns
    -------
    n, p, ntr at time = t 
    
    """
    
    n, ntr1 = z
    #print(str(t))
    
    p = n + ntr1
    
    krad = krad_c*1e-11
    
    Nt1 = Nt1_c*1e18
    b1n = b1n_c*1e-11 
    b1p = b1p_c*1e-11
    
    if (ntr1 > Nt1):
        ntr1 = Nt1
    
    #Defining of dn/dt; dp/dt; dptr/dt
    e_capture1 = b1n*n*(Nt1-ntr1)
    e_emission1 = N_c2*np.exp(-(Eg-Et1)/kT)*b1n*ntr1
    
    h_capture1 = b1p*p*ntr1
    h_emission1 = N_v2*np.exp(-(Et1)/kT)*b1p*(Nt1-ntr1)
        
    f_n =  G - krad*(n*p-ni2) - e_capture1 + e_emission1
    #f_p =  G - krad*(n*p-ni2) - h_capture1 + h_emission1    
    f_ntr1 =  e_capture1 - e_emission1 - h_capture1 + h_emission1
    
    return f_n, f_ntr1

