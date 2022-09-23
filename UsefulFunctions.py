import numpy as np
from scipy.optimize import curve_fit
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

def chi2(y,yfit):
    res = 0
    for i in range(len(yfit)):
        res = res + (y[i]-yfit[i])**2/(0.9724)**2 #The denominator should be the unbias Sipm voltage in mV
    return (res / len(yfit))

def waveform(x, Ci, start, m, end, d, A):
    condlist = [x < start, (x >= start) & (x <= end), x > end]
    funclist = [lambda x: Ci, lambda x: m*x+(Ci - m*start), lambda x: (A*(np.exp(x/d)-np.exp(end/d)) + m*(end-start) + Ci)]
    #(A*(np.exp(x/d)-np.exp(end/d)) + m*(end-start) + Ci)
    #(m*(end-start) + Ci)*np.exp((x-end)/d) + Cf
    
    return np.piecewise(x, condlist, funclist)

def gaus(x, normalization, location, stdev):
    return normalization* np.exp(-0.5*((x - location)/stdev)**2)

def Landau(x, normalization,location,stdev):
    #print(type(x))
    u=(x-location)*3.591/stdev/2.355
    renormalization = 1.64872*normalization
    return renormalization * np.exp(-u/2 - np.exp(-u)/2)

def LandauLinear(xs,normalization,location,stdev):
    if "array" in str(type(xs)) or "list" in str(type(xs)):
        linapprox=[]
        for x in xs:
            u=(x-location-0.05165303306066127)*0.42462845010615713/stdev
            if u < -3.8 or u > 500: linapp = 0
            else:
                for i in range(len(LandauXs)):
                    if LandauXs[i] < u: index = i
                    else: break
                linapp = normalization*(u-LandauXs[index])*(LandauYs[index+1]-LandauYs[index]) \
                /(LandauXs[index+1]-LandauXs[index])+renormalization*LandauYs[index]
            linapprox.append(linapp)
    else:
        u=(xs-location-0.05165303306066127)*0.42462845010615713/stdev
        if u < -3.8 or u > 500: linapprox = 0
        else:
            for i in range(len(LandauXs)):
                if LandauXs[i] < u: index = i
                else: break

            linapprox = normalization*(u-LandauXs[index])*(LandauYs[index+1]-LandauYs[index]) \
            /(LandauXs[index+1]-LandauXs[index])+renormalization*LandauYs[index]
    return linapprox

def findindex(xvals,xval): 
    minima = 66e6
    for i,elem in enumerate(xvals):
        # Update minimum difference from xval and index of the minimum difference
        if abs(xval-elem) < minima: 
            minima = abs(xval-elem) 
            index = i
    return index

def get_chi2(voltages):
    avg = np.mean(voltages[0:150])
    chisq = chi2(voltages[0:150],[avg]*150)
    return chisq

def get_amplitude_raw(voltages):
    v_max=np.max(voltages)
    v_min=np.min(voltages[0:findindex(voltages,v_max)+1])
    return v_max-v_min

def get_amplitude_base(voltages):
    v_min = np.mean(voltages[0:150])
    v_max=np.max(voltages)
    return v_max-v_min

def get_amplitude_smooth(voltages,degree):
    voltagesmooth = voltage_smooth(voltages,degree)
    v_min = np.mean(voltagesmooth[0:150])
    v_max=np.max(voltagesmooth)
    return v_max-v_min

def get_amplitude_fit(fitparams):
    linearrise = fitparams[2]*(fitparams[3]-fitparams[1])
    return linearrise

def get_time_raw(times,voltages):
    prevoltages = voltages[0:findindex(voltages,np.max(voltages))+1]
    halfamp = 0.5*(np.max(voltages)+np.min(prevoltages))
    halfindex = findindex(voltages,np.max(voltages)) #initially halfindex = maxindex
    haltime = 0
    while voltages[halfindex] - halfamp > 0:
        halfindex -= 1
    if abs(voltages[halfindex] - halfamp) < abs(voltages[halfindex + 1] - halfamp):
        halftime = times[halfindex]
    else: 
        halftime = times[halfindex + 1]
    return halftime

def get_time_base(times,voltages):
    halfamp = 0.5*(np.max(voltages)+np.mean(voltages[0:150]))
    halfindex = findindex(voltages,np.max(voltages)+1)
    haltime = 0
    while voltages[halfindex] - halfamp > 0 and halfindex > 0:
        halfindex -= 1
    if abs(voltages[halfindex] - halfamp) < abs(voltages[halfindex + 1] - halfamp):
        halftime = times[halfindex]
    else: 
        halftime = times[halfindex + 1]
    return halftime

def get_time_fit(fitparams):
    pulse_time = fitparams[1]/2+fitparams[3]/2
    return pulse_time

def get_time_fit_CDF(fitparams):
    pulse_time = fitparams[1]/2+fitparams[3]/2
    return pulse_time

def get_time_smooth(times,voltages,degree): # Mid point smoothing function
    voltagesmooth = voltage_smooth(voltages, degree) 

    halfamp = 0.5*(np.max(voltagesmooth)+np.sum(voltagesmooth[0:150])/151)
    halfindex = np.min([findindex(voltagesmooth,np.max(voltagesmooth))+1,len(voltagesmooth)-2])
    haltime = 0
    while ((voltagesmooth[halfindex] - halfamp) > 0 and (halfindex > 1)):
        halfindex -= 1
    if abs(voltagesmooth[halfindex] - halfamp) < abs(voltagesmooth[halfindex + 1] - halfamp):
        halftime = times[halfindex]
    else: 
        halftime = times[halfindex + 1]
    return halftime

def voltage_smooth(voltages, degree):
    degree = (int(degree/2))*2+1
    newarray = voltages[:-degree+1]/degree
    for deg in range(degree):
        if deg == 0: continue
        newarray += voltages[deg:np.size(voltages)-degree+deg+1]/degree
    newarray = np.insert(newarray,0,voltages[0:int((degree - 1) / 2)])
    newarray = np.append(newarray,voltages[-int((degree - 1) / 2):])
    
    return newarray

def get_dataframe(inputfiles, whichstats, channels=[1,2,3,4], rmscut=1.5, residualcut=5, p0=[(0,100,1,110,-100,100),(0,100,1,110,-100,100),(0,100,1,110,-100,100),(0,100,1,110,-100,100),(0,100,1,110,-100,100)],eventstart=1,verbose=False,show_residual=False,vieweventstart=1,viewevents=10):
    #error handling
    error = False
    if "list" not in str(type(inputfiles)):
        print("ERROR: inputfiles (first input) must be a list of strings.")
        error = True
    else:
        for inputfile in inputfiles:
            if "str" not in str(type(inputfile)): 
                print("ERROR: inputfiles (first input) must be a list of strings.")
                error = True
    if "list" not in str(type(whichstats)):
        print("ERROR: whichstats (second input) must be a list of six boolean values.")
        error = True
    elif len(whichstats) != 6: 
        print("ERROR: whichstats (second input) must be a list of six boolean values.")
        error = True
    else:
        for whichstat in whichstats:
            if not (whichstat == 0 or whichstat == 1): 
                print("ERROR: whichstats (second input) must be a list of ten boolean values.")
                error = True
    if "list" not in str(type(channels)):
        print("ERROR: channels must be a list of integers (1 through 4).")
        error = True
    else:
        for chan in channels:
            if not (chan == 1 or chan ==2 or chan ==3 or chan == 4): 
                print("ERROR: channels must be a list of integers (1 through 4).")
                error = True
    
    
    
    
    do_chi2 = whichstats[0]
    do_raw = whichstats[1]
    do_base = whichstats[2]
    do_smooth = whichstats[3]
    do_fit = whichstats[4]
    do_time_CDF = whichstats[5]
    
    
    #build the collumns of the dataframe
    #toc = ttime.perf_counter_ns()

    channelnames=[]
    for ch in channels:
        channelnames.append(f'ch{ch}')

    stats=[]
    if do_chi2: stats.append("chisq")
    if do_raw: 
        stats.append("P2P_raw")
        stats.append("time_raw")
    if do_base: 
        stats.append("P2P_base")
        stats.append("time_base")
    if do_smooth: 
        stats.append("P2P_smooth")
        stats.append("time_smooth")
    if do_fit: 
        stats.append("P2P_fit")
        stats.append("time_fit")
    if do_time_CDF: stats.append("time_CDF")

    filenumber = len(inputfiles)

    #Add each stat for each channel to the dataframe
    din = {a+"_"+b:[] for a in channelnames for b in stats}

    
    for i in range(0, filenumber):                             #iterate through files
        if error: continue
        with open(inputfiles[i]) as f:
            print(f"File: {inputfiles[i]}")
            current_file = (f.read().split('-- Event'))
        
        for j in range(eventstart, len(current_file)):                  #iterate through events len(current_file)
            if (verbose and (j < vieweventstart or j > viewevents+vieweventstart) and j%10==0) or (not verbose and j%10==0): print(f"Event: {j}",end="\r")
            
            #grab the data from each channel
            time = np.array([])
            voltage = [np.array([])]*4
            lines = current_file[j].split('\n')
            
            if verbose and (j >= vieweventstart) and (j <= viewevents+vieweventstart): #show the waveform fit line
                print(f"Event Number {j}")
                fig,ax = plt.subplots(1,4,figsize=(32,4))
                ax[0].set_title("Ch. 1")
                ax[1].set_title("Ch. 2")
                ax[2].set_title("Ch. 3")
                ax[3].set_title("Ch. 4")
                
            for line in lines[6:1028]:                         #iterate through data points
                values = line.split()
                time = np.append(time, float(values[2]))

                for channel in channels:
                    voltage[channel-1] = np.append(voltage[channel-1], float(values[channel+2]))
                
            
            #calculate stats for each channel
            for channel in channels:
                totalrms = sum((voltage[channel-1]-np.mean(voltage[channel-1]))**2)/len(voltage[channel-1])
                if totalrms < rmscut:
                    popt = (np.mean(voltage[channel-1]),0,0,0,1,0)
                    fit_voltage = waveform(time,*popt)
                    #remove "blips"
                    residual = np.abs(voltage[channel-1] - fit_voltage)

                    if do_chi2:
                        chisq = 0
                        din[f'ch{channel}_chisq'].append(chisq)
                            

                    #calculate amplitude & time
                    if do_raw: 
                        amplitude = 0
                        din[f'ch{channel}_P2P_raw'].append(amplitude)
                        pulse_time = 0
                        din[f'ch{channel}_time_raw'].append(pulse_time)

                    if do_base: 
                        amplitude = 0
                        din[f'ch{channel}_P2P_base'].append(amplitude)
                        pulse_time = 0
                        din[f'ch{channel}_time_base'].append(pulse_time)
                           
                    if do_smooth: 
                        amplitude = 0
                        din[f'ch{channel}_P2P_smooth'].append(amplitude)
                        pulse_time = 0
                        din[f'ch{channel}_time_smooth'].append(pulse_time)

                    if do_fit: 
                        amplitude = 0
                        din[f'ch{channel}_P2P_fit'].append(amplitude)
                        pulse_time = 0
                        din[f'ch{channel}_time_fit'].append(pulse_time)                      

                else:
                    #remove "blips"
                    smoothvoltage = voltage_smooth(voltage[channel-1],5)
                    residual = np.abs(voltage[channel-1] - smoothvoltage)
                    mask = residual < residualcut
                    for i, ele in enumerate(mask):
                        if ele == 0:# and np.abs(time[i] - get_time_fit(popt)) > 20:
                            voltage[channel-1][i] = smoothvoltage[i]
                            
                    if do_fit: 
                        popt, pcov = curve_fit(waveform, time, voltage[channel-1],p0=p0[channel-1],maxfev = 100000)#,bounds ([-10,10,1,60,-1000,0],[10,100,30,140,0,1000]))
                    
                    #calculate chi^2
                    if do_chi2:
                        chisq = get_chi2(voltage[channel-1])
                        din[f'ch{channel}_chisq'].append(chisq)
                            

                    #calculate amplitude and time
                    if do_raw: 
                        amplitude_raw = get_amplitude_raw(voltage[channel-1])
                        din[f'ch{channel}_P2P_raw'].append(amplitude_raw)
                        pulse_time_raw = get_time_raw(time,voltage[channel-1])
                        din[f'ch{channel}_time_raw'].append(pulse_time_raw)
                            
                    if do_base: 
                        amplitude_base = get_amplitude_base(voltage[channel-1])
                        din[f'ch{channel}_P2P_base'].append(amplitude_base)
                        pulse_time_base = get_time_base(time, voltage[channel-1])
                        din[f'ch{channel}_time_base'].append(pulse_time_base)
                            
                    if do_smooth: 
                        amplitude_smooth = get_amplitude_smooth(voltage[channel-1],5)
                        din[f'ch{channel}_P2P_smooth'].append(amplitude_smooth)
                        pulse_time_smooth = get_time_smooth(time, voltage[channel-1], 5)
                        din[f'ch{channel}_time_smooth'].append(pulse_time_smooth)

                    if do_fit: 
                        amplitude_fit = get_amplitude_fit(popt)
                        din[f'ch{channel}_P2P_fit'].append(amplitude_fit)
                        pulse_time_fit = get_time_fit(popt)
                        din[f'ch{channel}_time_fit'].append(pulse_time_fit)

                        
                if verbose and (j >= vieweventstart) and (j <= viewevents+vieweventstart): #show the waveform fit line
                    print(f"Channel {channel} RMS: {totalrms:.2f}")
                    if do_fit: print(f"fit params: {popt[0]:.2f}, {popt[1]:.1f}, {popt[2]:.2f}, {popt[3]:.1f}, {popt[4]:.2f}, {popt[5]:.3f}")
                    ts = np.linspace(0,np.max(time),501)
                    if totalrms > rmscut and do_fit: 
                        fits = waveform(ts,*popt)
                    else:
                        fits = np.mean(voltage[channel-1])*501
                    
                    #plot the data, smooth if smoothing, fit if fitting
                    ax[channel-1].plot(time,voltage[channel-1],label="raw")
                    if do_smooth: ax[channel-1].plot(time,smoothvoltage,label="smooth_5", color='orange')
                    if do_fit: ax[channel-1].plot(ts,fits,label="fit",color="black")
                    
                    
                    if show_residual: ax[channel-1].plot(time,residual,label="residual")

                    #draw the P2P and time 
                    if do_raw: 
                        baseline_raw = np.min(voltage[channel-1][0:findindex(voltage[channel-1],np.max(voltage[channel-1]))+1])
                        ax[channel-1].vlines(pulse_time_raw, ymin=baseline_raw, ymax=amplitude_raw+baseline_raw, color='r', label="t_raw")
                        ax[channel-1].hlines(baseline_raw, xmin=pulse_time_raw-30, xmax=pulse_time_raw+30, color='r')
                        ax[channel-1].hlines(amplitude_raw+baseline_raw, xmin=pulse_time_raw-30, xmax=pulse_time_raw+30, color='r', label="A_raw")
                        
                    if do_base: 
                        baseline_base = np.mean(voltage[channel-1][0:150])
                        ax[channel-1].vlines(pulse_time_base, ymin=baseline_base, ymax=amplitude_base+baseline_base, color='b',label="t_base")
                        ax[channel-1].hlines(baseline_base, xmin=pulse_time_base-30, xmax=pulse_time_base+30, color='b')
                        ax[channel-1].hlines(amplitude_base+baseline_base, xmin=pulse_time_base-30, xmax=pulse_time_base+30, color='b',label="A_base")
                        
                    if do_smooth: 
                        baseline_smooth = np.mean(smoothvoltage[0:151])
                        ax[channel-1].vlines(pulse_time_smooth, ymin=baseline_smooth, ymax=amplitude_smooth+baseline_smooth, color='purple',label="t_smooth_5")
                        ax[channel-1].hlines(baseline_smooth, xmin=pulse_time_smooth-30, xmax=pulse_time_smooth+30, color='purple')
                        ax[channel-1].hlines(amplitude_smooth+baseline_smooth, xmin=pulse_time_smooth-30, xmax=pulse_time_smooth+30, color='purple', label="A_smooth")
                        
                    if do_fit: 
                        ax[channel-1].vlines(pulse_time_fit, ymin=popt[0], ymax=popt[0]+amplitude_fit, color='g',label="t_fit")
                        ax[channel-1].hlines(popt[0], xmin=pulse_time_fit-30, xmax=pulse_time_fit+30, color='g')
                        ax[channel-1].hlines(amplitude_fit+popt[0], xmin=pulse_time_fit-30, xmax=pulse_time_fit+30, color='g',label="A_fit")

            if verbose and (j >= vieweventstart) and (j <= viewevents+vieweventstart):
                plt.legend()        
                plt.show()
        print("")
    df = pd.DataFrame(din)                         
    print(f"Done! Total events analyzed: {len(df)}")

    return df