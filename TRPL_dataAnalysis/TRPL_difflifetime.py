import matplotlib.pyplot as plt
import scipy
import numpy as np
import pandas as pd
#from scipy.signal import savgol_filter

def fitfunc(x, *args):
        """
        
        Used for fitting of a multi (varying amount defined by len(args)) of exponential decay.

        Parameters
        ----------
        x: evaluation points of function.
        
        *args : arguments for multiexponential decay. len(args) has to be 2n + 1.
            
        Returns
        -------
        f = arg[0] + arg[1]*exp(-arg[2]*x) + ... + arg[i]*exp(-arg[i+1]*x)
        
        """
        #print(args)
        params = np.array([arg for arg in args])
        
        if ((len(params) == 1) or (len(params) > 20)):
            print("Number of params is wrong n = "+str(len(params))+"\n")
        
        s = params[0]
        for i, a in enumerate(params):
            if (i == 0):
                continue
            if not(i % 2):
                continue
            else:
                s = s + params[i]*np.exp(-params[i+1]*x)
                
        return (s)

def fitfunc2(x, *args):

    params = np.array([arg for arg in args])

    if ((len(params) < 2) or (len(params) > 20)):
            print("Number of params in fitfunc is wrong n = "+str(len(params))+"\n")

    numerator = 0
    denum = 0
    for i in range(int(len(params)/2)):
        numerator = numerator + params[i]*(x**i)
        denum = denum + params[int(len(params)/2)+i]*(x**i)
    

    return numerator/denum


def fit_difflifetimes(time, TRPL_raw, n0, noise_level, n_exp = 3, l2 = None):
    """
    Calculates differential lifetime values, given a raw TRPL table, using a arbritrary amount of exponentials to fit the data.

    Parameters
    ----------
    time: time values, array like.
    
    TRPL_denoised : denoised TRPL values, array-like
    
    powers : Powers corresponding to the second dimension of TRPL_raw, array-like
    
    thickness : thickness of the samples, in cm.
    
    l2: number of data points considered for fitting.
        
    Returns
    -------
    time_fit: x_axis considered for fit
    densities2: the calculated carrier densities following n0*sqrt(data)
    diff_taus: diffenrential lifetimes, tau_diff = -2*(dt/d(log(fit))). following Thomas Kirchartzs' work
    
    """
    ns_raw = TRPL_raw.transpose()
    
    time = time.transpose()
    
    diff_taus = []
    densities2 = []
    time_fit = []
    
    print("Number of exponentials for fit used is = "+str(n_exp)+"\n")
    
    #previous_ps = [1, 1, 1, 1, 1, 1, 1]
    f, ax = plt.subplots(3, len(n0), figsize=(25,12))
    for i in range(len(n0)):
        t = time[i, (ns_raw[i, :] > 2*noise_level[i]) & (time[i,:] > 0)]
        pl = ns_raw[i, (ns_raw[i, :] > 2*noise_level[i]) & (time[i,:] > 0)]
        #initial guess for fitting
        p = [1]*(2*n_exp[i]+1)
        #p = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        if (l2 == None):
            #previous_ps, pcov = scipy.optimize.curve_fit(fitfunc, t[(pl >= 0)], np.log(pl[(pl >= 0)]), maxfev = 150000, p0 = p) 
            previous_ps, pcov = scipy.optimize.curve_fit(fitfunc, t, (pl), maxfev = 150000, p0 = p) 
            #np.exp(np.log(start = 1, stop = 1+len(t[t < 1e-9*l2[i]]))))#method = 'lm', loss)
            fit = (fitfunc(t, *previous_ps))       
            tau_diff = -1*(np.diff(t)/np.diff(np.log(fit)))
            print("L2 is None")
        else:
            #previous_ps, pcov = scipy.optimize.curve_fit(fitfunc, t[(t < 1e-9*l2[i]) & (pl >= 0)], np.log(pl[(t < 1e-9*l2[i]) & (pl >= 0)]), maxfev = 150000, p0 = p)
            previous_ps, pcov = scipy.optimize.curve_fit(fitfunc, t[(t < 1e-9*l2[i])], pl[(t < 1e-9*l2[i])], maxfev = 1500000, p0 = p)
            #np.exp(np.log(start = 1, stop = 1+len(t[t < 1e-9*l2[i]]))))#method = 'lm', loss)
            fit = (fitfunc(t[t < 1e-9*l2[i]], *previous_ps))
            tau_diff = -2*(np.diff(t[t < 1e-9*l2[i]])/np.diff(np.log(fit)))
    

        carrier_densities_fit = np.sqrt(fit/np.max(fit))*(n0[i])
        
        #Plotting
        ax[0, i].scatter(1e9*t, pl, marker = 'x')
        ax[0, i].plot(1e9*t[:len(fit)], (fit), color = 'orange')
        ax[0, i].set_yscale("log")
        ax[0, i].set_xlim([min(1e9*t), l2[i]*2])
        ax[0, i].set_xlabel("time [ns]")
        ax[0, i].set_ylabel("PL counts [#]")

        ax[1, i].plot(1e9*t[:len(tau_diff)], tau_diff)
        ax[1, i].set_xlim([min(1e9*t), max(1e9*t[:len(tau_diff)])])
        ax[1, i].set_xlabel("time [ns]")
        ax[1, i].set_ylabel("Differential lifetime [s]")
        
        
        ax[2, i].plot(carrier_densities_fit[1:], tau_diff)
        ax[2, i].set_xlabel("Carrier Concentration [cm-3]")
        ax[2, i].set_ylabel("Differential lifetime [s]")
        ax[2, i].set_xscale("log")
        ax[2, i].set_yscale("log")
    
        densities2.append(carrier_densities_fit)
        diff_taus.append(tau_diff)
        time_fit.append(t[:len(fit)])
        
    f = plt.figure()
    for i, (b, c) in enumerate(zip(diff_taus, densities2)):
        plt.plot(c[:-1], b)
        #plt.plot(c[0:cut], savgol_filter(b[0:cut], 20, 7, mode = "nearest"))
        
    ax = f.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel("Carrier Density [cm-3]")
    ax.set_ylabel("Differential lifetime [s]")
    plt.legend(range(len(ns_raw)), loc = 'upper left')
    
    return time_fit, densities2, diff_taus

def diff_savedata(files_selection, times, densities, tau_diffs, filename, BG, Nc, Nv, kT = 27.7*1e-3):
    """
    Saves the differential lifetime values. 

    Parameters - They come from the output of the fit_difflifetimes function
    ----------
    files_selection: file names, array like.

    times: time values, array like.
    
    densities : densities, array-like
    
    tau_diffs : differential lifetimes, array-like
    
    filenames : name of the .csv file to be saved, path-like
    
    BG, Nc, Nv: parameters used for calculating qfls, floats.
        
    Returns
    -------
    data: Dataframe that was used to save the .csv file.
    
    """

    data = pd.DataFrame()
    kT = 27.7*1e-3 #eV
    ni = np.sqrt(Nc*Nv*np.exp(-BG/(kT)))
    l = len(times[0])
    for i, (time, density, tau_d) in enumerate(zip(times, densities, tau_diffs)): 
        qfls = kT*np.log((density[:-1]*density[:-1])/(ni*ni))
        if (i == 0):
            how = 'right'
        elif(len(time) > l):
            how = 'right'
        else:
            how = 'left'
        
        data = data.join(pd.Series(time[:-1]).rename("time: "+files_selection[i]), how = how)
        data = data.join(pd.Series(density[:-1]).rename("density: "+files_selection[i]), how = how)
        data = data.join(pd.Series(qfls).rename("qfls: "+files_selection[i]), how = how)
        data = data.join(pd.Series(tau_d).rename("tau_diff: "+files_selection[i]), how = how)

        l = len(time)

    data.to_csv(filename)

    return data
