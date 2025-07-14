import sys
sys.path.append(r"C://PyDAQmx-1.4.6")
import time
import ni_pulser_loop_new as DAQpulser
from snAPI.Main import *
import matplotlib
matplotlib.use('TkAgg',force=True)
from matplotlib import pyplot as plt
import pandas as pd
import os
from pylablib.devices import Thorlabs

print("Switched to:",matplotlib.get_backend())

def launch_trPL(sn: snAPI, time_s: int, save_directory: float, syncEdge: bool = 1, syncLevel:int = -300, syncOffset:int = 0, chan1Edge:bool = 0, chan1Level:int = -120, chan1Offset:int = 0, resolution_s:float = 1e-9, chanStop:int = 5):
    #sn = snAPI()
    # get first available device
    sn.getDevice()
    sn.setLogLevel(logLevel=LogLevel.DataFile, onOff=True)
    #initialize the device
    sn.initDevice(MeasMode.Histogram)
    # set the configuration for your device type
    binning_factor = int(np.log2(resolution_s/(250*1e-12)))
    sn.device.setHistoLength(chanStop) #2^{(10 + \mathrm{lengthCode})}
    sn.device.setBinning(binning_factor) # binsize = 250ps * 2**n
    #sync
    sn.device.setSyncEdgeTrig(syncLevel, syncEdge)
    sn.device.setSyncChannelOffset(syncOffset)
    #input channel
    sn.device.setInputEdgeTrig(0,chan1Level,chan1Edge)
    sn.device.setInputChannelOffset(0, chan1Offset)
    cntRt = sn.getCountRates()
    print(cntRt)
    if(cntRt[0] > 8e6):
       print("Count Rates are too high, counts = ", cntRt)
    else:
        sn.histogram.measure(acqTime = int(1e3*time_s), savePTU=True, waitFinished=False)

    while True:
        finished = sn.histogram.isFinished()
        data, bins = sn.histogram.getData()
        plt.figure(1)
        if len(data):
            plt.clf()
            #plt.semilogy(bins, data[0], linewidth=2.0, label='sync')
            for c in range(1, 1+sn.deviceConfig["NumChans"]):
                plt.semilogy(bins, data[c], linewidth=2.0, label=f'chan{c}', alpha = 0.5)
            plt.xlabel('Time [ps]')
            plt.ylabel('Counts')
            plt.xlim([0,0.5e6])
            plt.ylim([0.1,1e5])
            plt.legend()
            plt.title("Counts / Time")
            plt.pause(1)

        if finished:
            measurement = np.array([bins, data[1]])
            df = pd.DataFrame(data = measurement.transpose(), columns = ["bins [ps]", "counts [#]"])
            newname = list(os.path.splitext(os.path.basename(save_directory)))
            newname.append(newname[-1])
            #read the thorlabs filter wheel
            wheel = Thorlabs.FW("COM8")
            add_on = "_"+str(time_s)+"s_"+str(int(1e-3*cntRt[0]))+"kHz_ND"+str(wheel.get_position())+"_"+str(int(cntRt[1]))+"cps"
            wheel.close()
            newname[1] = add_on
            a = ""
            for n in newname:
                a = str(a)+str(n)

            #Save file with measurement parameters
            save_directory = os.path.join(os.path.dirname(save_directory), a)
            df.to_csv(save_directory, sep="\t")
            wheel.close()
            break

    return

def launch_LPtrPL(sn: snAPI, p: DAQpulser, time_s: int, pulseLength_s: int, repRate_Hz: int, save_directory: float, syncEdge: bool = 0, syncLevel:int = -250, syncOffset_ps:int = -99999, chan1Edge:bool = 1, chan1Level:int = -120, chan1Offset_ps:int = 99999, resolution_s:float = 1e-9, chanStop:int = 5):
    #sn = snAPI()
    p.start_pulses(pulseLength_s, repRate_Hz)
    # get first available device
    sn.getDevice()
    sn.setLogLevel(logLevel=LogLevel.DataFile, onOff=True)
    #initialize the device
    sn.initDevice(MeasMode.Histogram)
    # set the configuration for your device type
    binning_factor = int(np.log2(resolution_s/(250*1e-12)))
    sn.device.setHistoLength(chanStop) #2^{(10 + \mathrm{lengthCode})}
    sn.device.setBinning(binning_factor) # binsize = 250ps * 2**n
    #sync
    sn.device.setSyncEdgeTrig(syncLevel, syncEdge)
    sn.device.setSyncChannelOffset(syncOffset_ps)
    #input channel
    sn.device.setInputEdgeTrig(0,chan1Level,chan1Edge)
    sn.device.setInputChannelOffset(0, chan1Offset_ps)
    cntRt = sn.getCountRates()
    time.sleep(10.0)
    if(cntRt[1] > 8e6):
       print("Count Rates are too high, ND wheel is changed = ", cntRt)
    else:
        sn.histogram.measure(acqTime = int(1e3*time_s),savePTU=True, waitFinished=False)

    while True:
        finished = sn.histogram.isFinished()
        data, bins = sn.histogram.getData()
        plt.figure(10)
        if len(data):
            plt.clf()
            #plt.semilogy(bins, data[0], linewidth=2.0, label='sync')
            for c in range(1, 1+sn.deviceConfig["NumChans"]):
                plt.semilogy(bins, data[c], linewidth=2.0, label=f'chan{c}', alpha = 0.5)
            plt.xlabel('Time [ps]')
            plt.ylabel('Counts')
            plt.xlim([0,5e6])
            plt.ylim([0.1,1e5])
            plt.legend()
            plt.title("Counts / Time")
            plt.pause(1)

        if finished:
            measurement = np.array([bins, data[1]])
            df = pd.DataFrame(data = measurement.transpose(), columns = ["bins [ps]", "counts [#]"])
            newname = list(os.path.splitext(os.path.basename(save_directory)))
            newname.append(newname[-1])
            #read the thorlabs filter wheel
            wheel = Thorlabs.FW("COM8")
            add_on = "_"+str(time_s)+"s_"+str(int(cntRt[0]))+"Hz_ND"+str(wheel.get_position())+"_PW-"+str(1e6*pulseLength_s)+"us_"+str(int(cntRt[1]))+"cps"
            wheel.close()
            newname[1] = add_on
            a = ""
            for n in newname:
                a = str(a)+str(n)

            #Save file with measurement parameters
            save_directory = os.path.join(os.path.dirname(save_directory), a)
            df.to_csv(save_directory, sep="\t")
            wheel.close()
            p.stop_tasks()
            p.clear_tasks()
            break

    return

if(__name__ == "__main__"):
    #Set Devices
    sn = snAPI()
    p = DAQpulser.NationalInstrumentsPulser(clock_channel="ctr_ch1", device='Dev2', high_time=25e-9, low_time=1e-3)
    t = 7200 #s
    directory = r"D:\\Maxim\\20250128-LPtrPL-FAPI\\FAPI27-Spot3\\"
    reprate = 1000
    PLs = [10e-6, 100e-6, 700e-6]
    
    # file = "M04-MAPI5-20uW.dat"
    # filepath = directory+file
    # print(filepath)
    # launch_LPtrPL(sn, p, t, float(700e-6), reprate, filepath, resolution_s = 1000e-12, syncEdge = 1, syncLevel = 300)
    name = "-LPtrPL-FAPI27-Spot3-0p2uW.dat"
    
    for i, pulse in enumerate(PLs):
        file = "M0"+str(i+5)+name
        filepath = directory+file
        print(filepath)
        launch_LPtrPL(sn, p, t, float(pulse), reprate, filepath, resolution_s = 1000e-12, syncEdge = 0, syncLevel = -300)