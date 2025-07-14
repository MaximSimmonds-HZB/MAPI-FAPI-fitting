import pandas as pd
import numpy as np
from os.path import isfile, join
import os

def readout_function(filepath, integration_time_seconds, reprate, denoise = False, retime = False, mode = "manual"):
    """
    Imports a TRPL data file.  

    Parameters
    ----------
    filepath : Path to the TRPL data file. 
    
    integration_time_seconds : Integration time of the measurement. 
    
    reprate: repetition rate of all the experiments. 
    
    denoise: denoise values? True: Automatic denoising (careful, check the procedure). False: No denoising. 
                            Value: The value is removed from the counts.
                            
    
    retime: shift time axis so that max(counts) is at t = 0.
    
    
    Returns
    -------
    time vector, TRPL normalised vector, TRPL (denoised) vector, TRPL raw data vector, Noise level
    
    """

    if(mode == "manual"):
        df = pd.read_csv(filepath, skiprows=8, header = None, encoding='latin-1', delimiter = '\t')
        binsize = 1e-9*df.iloc[0].astype('float64')[0]
        print(binsize)
        TRPL = df.iloc[2:,0].to_numpy(dtype = None).astype('int')
    elif(mode == "auto"):
        df = pd.read_csv(filepath, delimiter = '\t')
        df = df.drop(df.columns[0:1], axis=1)
        t = 1e-12*df["bins [ps]"].to_numpy(dtype = None).astype('int')
        binsize = t[1]-t[0]
        print(binsize)
        TRPL = df["counts [#]"].to_numpy(dtype = None).astype('int')
        t = np.transpose(t)
    elif(mode == "wannsee"):
        df = pd.read_csv(filepath, skiprows=11, header = None, on_bad_lines='skip')
        binsize = 1e-12*(df.iloc[1, 0].astype('float64') - df.iloc[0, 0].astype('float64'))
        TRPL = df.iloc[0:,1].to_numpy(dtype = None).astype('int')
        t = 1e-12*(df.iloc[0:,0].to_numpy(dtype = None).astype('float64'))
        print(binsize)

    TRPL = np.transpose(TRPL)
    #Remove Uniform noise
    if(not(denoise)):
        noise = 0
        #noise = np.mean(TRPL[0][(np.argmax(TRPL[0][:])-noise_start_i):(np.argmax(TRPL[0][:])-noise_end_i)])
    elif(denoise < 0):  
        noise = np.mean(np.trim_zeros(TRPL, trim='b')[denoise:]) #takes the zeros out of the end of the data
        #noise = np.mean(TRPL[denoise:])
    elif(denoise > 0):
        noise = np.mean(TRPL[:denoise])
         
    TRPL_denoise = TRPL[:] - noise
    TRPL_denoise = rate_calculation_function(TRPL_denoise, integration_time_seconds, binsize, reprate)
    #TRPL = rate_calculation_function(TRPL, integration_time_seconds, binsize, reprate)
        
    #Normalize
    TRPL_n = TRPL_denoise/np.amax(TRPL_denoise)
    
    if (retime):
        if(mode == "manual"):
            t_TRPL = binsize*(np.arange(len(TRPL_n))-np.argmax(TRPL_n))
        elif(mode == "auto" or mode == "wannsee"):
            t_TRPL = t-t[np.argmax(TRPL)]
    else:
        if(mode == "manual"):
            t_TRPL = binsize*(np.arange(len(TRPL_n)))
        else:
            t_TRPL = t
    
    print(noise)
    return t_TRPL, TRPL_n, TRPL_denoise, TRPL, noise

def TRPL_folder_read(folderpath, integration_time_seconds, reprate, denoise = False, retime = False, mode = "manual"):
    """
    Imports a folder of TRPL files.  

    Parameters
    ----------
    folderpath : Path to the folder containing data. 
    
    integration_time_seconds : Integration time of the measurement. 
    
    reprate: repetition rate of all the experiments. 
    
    denoise: denoise values? True: Automatic denoising (careful, check the procedure). False: No denoising. 
                            Value: The value is removed from the counts.
                            
    
    retime: shift time axis so that max(counts) is at t = 0.
    
    
    Returns
    -------
    file names, time array, TRPL normalised array, TRPL (denoised), TRPL raw data, Noise levels
    
    """
    
    if(mode == "wannsee"):
        files = [f for f in os.listdir(folderpath) if (isfile(join(folderpath, f)) and f.endswith(".csv") and (not(f.startswith("._"))))]
        files.sort()
    else:
        files = [f for f in os.listdir(folderpath) if (isfile(join(folderpath, f)) and f.endswith(".dat") and (not(f.startswith("._"))))]
        files.sort()
    
    ts, TRPLs_n, TRPL_subsMean, TRPL_raw, noise = readout_function(join(folderpath, files[0]), integration_time_seconds[0], reprate[0], denoise = denoise, retime = retime, mode = mode)
    
    Noise = [noise]
    for i, f in enumerate(files):
        if(i == 0):
            continue
        
        t, trpl_n, trpl_subsMean, trpl_raw, noise = readout_function(join(folderpath, f), integration_time_seconds[i], reprate[i], denoise = denoise, retime = retime, mode = mode)
        
        ts = np.c_[ts, t]
        TRPLs_n = np.c_[TRPLs_n, trpl_n]
        TRPL_subsMean = np.c_[TRPL_subsMean, trpl_subsMean]
        TRPL_raw = np.c_[TRPL_raw, trpl_raw]
        Noise.append(noise)
        
    return files, ts, TRPLs_n, TRPL_subsMean, TRPL_raw, Noise

def rate_calculation_function(count, integration_time_seconds, binsize_seconds, reprate, tau_COUNT_APD = 45e-9):
    """
    Corrects nonlinearities of counts for higher Rates. 

    Parameters
    ----------
    count : Measured Counts, array like. 
    
    integration_time_seconds : Integration time of the measurement. 
    
    binsize: binsize in ns. 
    
    tau_COUNT_APD: Value for the dead time of the APD, found on the Laser Components COUNT manual. 
        
    Returns
    -------
    Measured Rate
    
    """
    rate_measured = count/(binsize_seconds*integration_time_seconds*reprate)
    
    return rate_measured

def correction_function(count, integration_time_seconds, binsize_nanoseconds, reprate, tau_COUNT_APD = 45e-9):
    """
    Corrects nonlinearities of counts for higher Rates. 

    Parameters
    ----------
    count : Measured Counts, array like. 
    
    integration_time_seconds : Integration time of the measurement. 
    
    binsize: binsize in ns. 
    
    tau_COUNT_APD: Value for the dead time of the APD, found on the Laser Components COUNT manual. 
        
    Returns
    -------
    Corrected Counts
    
    """
    rate_measured = count/(1e-9*binsize_nanoseconds*integration_time_seconds*reprate)
    rate_actual = rate_measured/(1-rate_measured*tau_COUNT_APD)
    
    return rate_actual*(binsize_nanoseconds*integration_time_seconds*reprate)

def SPV_folder_read(folderpath, setup = "Dittrich"):

    if (setup == "Dittrich"):
        files = [f for f in os.listdir(folderpath) if (isfile(join(folderpath, f)) and f.endswith("txt"))]
    elif(setup == "Wannsee"):
        files = [f for f in os.listdir(folderpath) if (isfile(join(folderpath, f)) and f.endswith(".Wfm.csv"))]
    
    files.sort()
    
    dfs = import_SPV(join(folderpath, files[0]), setup = setup)
    #SPVs = SPVs.transpose()
    
    for i, f in enumerate(files):
        if(i == 0):
            continue
        
        df = import_SPV(join(folderpath, f), setup = setup)
        
        #ts = np.c_[ts, t]
        dfs = pd.concat([dfs, df], axis = 1)
        #SPVs = SPVs.set_index(spv.index).join(spv.set_index(spv.index), lsuffix = '_'+str(i-1), rsuffix = '_'+str(i))
        #SPVs = np.c_[SPVs, spv]
        #SPVs = SPVs.join(spv)
        
    return files, dfs

def import_SPV(filepath, setup):
    
    if (setup == "Dittrich"):
        df = pd.read_csv(filepath, skiprows=0, sep = "\t", comment='#')
        names = [os.path.basename(filepath)[:-4]+": "+str(n) for n in df.columns]
        df.columns = names
    elif(setup == "Wannsee"):
        df = pd.read_csv(filepath, skiprows=0, header = None, sep = ",")
        cols = ["time [s]", "voltage [V]"]
        names = [os.path.basename(filepath)+": "+cols[n] for n in range(len(df.columns))]
        df.columns = names
    
    #t = df.iloc[:,0].to_numpy().astype("float64")
    #SPVs = df.iloc[:,0].to_numpy().astype("float64")

    #print(df)
  
    return df

def TRPL_savedata(*arrays, times, selection, colnames, filename, files, sep=',', col_filenames = True):

    """
    Saves the trPL arrays in the specified order. 

    Parameters - They come from the output of the TRPL_folder_read function
    ----------
    *arrays: number of arrays to be saved next to time values, array, array, array, like.

    times: time values, array like.
    
    selection : slection of measurements to be saved, array-like
    
    colnames : prefix for column names corresponding to the *arrays, array-like
    
    filename : name of the .csv file to be saved, path-like
    
    files: file names corresponding to the times, array-like.
        
    Returns
    -------
    data: Dataframe that was used to save the .csv file.
    
    """

    data = pd.DataFrame()
    l = len(times[0])
    for i, time in enumerate((times[:, selection].transpose())):
        if (i == 0):
            how = 'right'        
        elif (len(time) > l):
            how = 'right'
        else:
            how = 'left'
        
        if (col_filenames):
            data = data.join(pd.Series(time).rename("time: "+files[selection[i]]), how = how)
        else:
            data = data.join(pd.Series(time).rename("time"), how = how)
    
        for j, array in enumerate(arrays):
            if (col_filenames):
                #data[] = pd.Series()
                data = data.join(pd.Series(array[:, selection].transpose()[i]).rename(colnames[j]+": "+files[selection[i]]), how = how)
            else:
                data = data.join(pd.Series(array[:, selection].transpose()[i]).rename(colnames[j]), how = how)
                #data[colnames[j]] = pd.Series(array[:, selection].transpose()[i])
        #print(elements)

        l = len(time)

    data.to_csv(filename, index=None, sep = sep)
    return data

