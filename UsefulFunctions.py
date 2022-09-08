import numpy as np
from scipy.optimize import curve_fit
import numpy as np, pandas as pd

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

def get_chi2(fitparams,times,voltages):
    startindex=1
    for i in range(len(times)): 
        if times[i] < fitparams[1]: startindex = i + 1
        else: break
    chisq = chi2(voltages[0:startindex],waveform(times[0:startindex], *fitparams))
    return chisq

def get_amplitude_raw(voltages):
    v_max=max(voltages)
    v_min=min(voltages[0:findindex(voltages,v_max)+1])
    return v_max-v_min

def get_amplitude_base(fitparams,voltages):
    v_min = fitparams[0]
    v_max=max(voltages)
    return v_max-v_min

def get_amplitude_smooth(fitparams,voltages,degree):
    voltagesmooth = voltage_smooth(voltages,degree)
    v_min = fitparams[0]
    v_max=max(voltagesmooth)
    return v_max-v_min

def get_amplitude_fit(fitparams):
    linearrise = fitparams[2]*(fitparams[3]-fitparams[1])
    return linearrise

def get_time_raw(times,voltages):
    prevoltages = voltages[0:findindex(voltages,max(voltages))+1]
    halfamp = 0.5*(max(voltages)+min(prevoltages))
    halfindex = findindex(voltages,max(voltages)) #initially halfindex = maxindex
    haltime = 0
    while voltages[halfindex] - halfamp > 0:
        halfindex -= 1
    if abs(voltages[halfindex] - halfamp) < abs(voltages[halfindex + 1] - halfamp):
        halftime = times[halfindex]
    else: 
        halftime = times[halfindex + 1]
    return halftime

def get_time_base(fitparams,times,voltages):
    halfamp = 0.5*(max(voltages)-fitparams[0]) + fitparams[0]
    halfindex = findindex(voltages,max(voltages)+1)
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

def get_time_smooth(fitparams,times,voltages,degree): # Mid point smoothing function
    voltagesmooth = voltage_smooth(voltages, degree) 

    halfamp = 0.5*(max(voltagesmooth)-fitparams[0]) + fitparams[0]
    halfindex = findindex(voltagesmooth,max(voltagesmooth)+1)
    haltime = 0
    while voltagesmooth[halfindex] - halfamp > 0:
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
        newarray += myarray[deg:np.size(myarray)-degree+deg+1]/degree
    newarray = np.insert(newarray,0,myarray[0:int((degree - 1) / 2)])
    newarray = np.append(newarray,myarray[-int((degree - 1) / 2):])
    
    return newarray

def get_dataframe(inputfiles, channelnum, rmscut, residualcut, whichstats, p0):
    do_chi2 = whichstats[0]
    do_amplitude_raw = whichstats[1]
    do_amplitude_base = whichstats[2]
    do_amplitude_smooth = whichstats[3]
    do_amplitude_fit = whichstats[4]
    do_time_raw = whichstats[5]
    do_time_base = whichstats[6]
    do_time_fit = whichstats[7]
    do_time_CDF = whichstats[8]
    do_time_smooth = whichstats[9]
    
    
    #build the collumns of the dataframe
    #toc = ttime.perf_counter_ns()

    channelnames=[]
    for ch in channelnum:
        channelnames.append(f'ch{ch}')

    stats=[]
    if do_chi2: stats.append("chisq")
    if do_amplitude_raw: stats.append("P2P_raw")
    if do_amplitude_base: stats.append("P2P_base")
    if do_amplitude_smooth: stats.append("P2P_smooth")
    if do_amplitude_fit: stats.append("P2P_fit")
    if do_time_raw: stats.append("time_raw")
    if do_time_base: stats.append("time_base")
    if do_time_fit: stats.append("time_fit")
    if do_time_CDF: stats.append("time_CDF")
    if do_time_smooth: stats.append("time_smooth")

    filenumber = len(inputfiles)

    #Add each stat for each channel to the dataframe
    din = {a+"_"+b:[] for a in channelnames for b in stats}

    
    for i in range(0, filenumber):                             #iterate through files
            
        with open(inputfiles[i]) as f:
            current_file = (f.read().split('-- Event'))
        
        for j in range(1, len(current_file)):                  #iterate through events len(current_file)
                #grab the data from each channel
            time = np.array([])
            voltage = [np.array([])]*4
            lines = current_file[j].split('\n')
            for line in lines[6:1028]:                         #iterate through data points
                values = line.split()
                time = np.append(time, float(values[2]))

                for channel in channelnum:
                    voltage[channel-1] = np.append(voltage[channel-1], float(values[channel+2]))
                
            
            #calculate stats for each channel
            for channel in channelnum:
                totalrms = sum((voltage[channel-1]-np.mean(voltage[channel-1]))**2)/len(voltage[channel-1])
                if totalrms < rmscut:
                    popt = (np.mean(voltage[channel-1]),0,0,0,1,0)
                    fit_voltage = uf.waveform(time,*popt)
                    #remove "blips"
                    residual = np.abs(voltage[channel-1] - fit_voltage)

                    if do_chi2:
                        chisq = 0
                        din[f'ch{channel}_chisq'].append(chisq)
                            

                    #calculate amplitude
                    if do_amplitude_raw: 
                        amplitude = 0
                        din[f'ch{channel}_P2P_raw'].append(amplitude)

                    if do_amplitude_base: 
                        amplitude = 0
                        din[f'ch{channel}_P2P_base'].append(amplitude)
                           
                    if do_amplitude_smooth: 
                        amplitude = 0
                        din[f'ch{channel}_P2P_smooth'].append(amplitude)

                    if do_amplitude_fit: 
                        amplitude = 0
                        din[f'ch{channel}_P2P_fit'].append(amplitude)
                      

                    #calculate time
                    if do_time_raw: 
                        pulse_time = 0
                        din[f'ch{channel}_time_raw'].append(pulse_time)
                            
                        
                    if do_time_base:
                        pulse_time = 0
                        din[f'ch{channel}_time_base'].append(pulse_time)
                      

                    if do_time_fit: 
                        pulse_time = 0
                        din[f'ch{channel}_time_fit'].append(pulse_time)
                            
                        
                    if do_time_smooth:
                        pulse_time = 0
                        din[f'ch{channel}_time_smooth'].append(pulse_time)
                      

                else:
                    popt, pcov = curve_fit(waveform, time, voltage[channel-1],p0=p0[channel-1],
                                           maxfev = 100000)#,bounds=([-10,60,0,60,0,-1],[10,140,100,140,3000,0])
                    
                    fit_voltage = waveform(time,*popt)
                    #remove "blips"
                    residual = np.abs(voltage[channel-1] - fit_voltage)
                    mask = residual < residualcut
                    for i, ele in enumerate(mask):
                        if ele == 0:# and np.abs(time[i] - get_time_fit(popt)) > 20:
                            voltage[channel-1][i] = waveform(time[i],*popt)
                    
                    #calculate chi^2
                    if do_chi2:
                        chisq = get_chi2(popt,time,voltage[channel-1])
                        din[f'ch{channel}_chisq'].append(chisq)
                            

                    #calculate amplitude
                    if do_amplitude_raw: 
                        amplitude = get_amplitude_raw(voltage[channel-1])
                        din[f'ch{channel}_P2P_raw'].append(amplitude)
                            

                    if do_amplitude_base: 
                        amplitude = get_amplitude_base(popt,voltage[channel-1])
                        din[f'ch{channel}_P2P_base'].append(amplitude)
                            
                    if do_amplitude_smooth: 
                        amplitude = get_amplitude_smooth(popt,voltage[channel-1],5)
                        din[f'ch{channel}_P2P_smooth'].append(amplitude)

                    if do_amplitude_fit: 
                        amplitude = get_amplitude_fit(popt)
                        din[f'ch{channel}_P2P_fit'].append(amplitude)
                            

                    #calculate time
                    if do_time_raw: 
                        pulse_time = get_time_raw(time,voltage[channel-1])
                        din[f'ch{channel}_time_raw'].append(pulse_time)
                            

                    if do_time_base:
                        pulse_time = get_time_base(popt, time, voltage[channel-1])
                        din[f'ch{channel}_time_base'].append(pulse_time)
                            

                    if do_time_fit: 
                        pulse_time = get_time_fit(popt)
                        din[f'ch{channel}_time_fit'].append(pulse_time)
                           
                        
                    if do_time_smooth:
                        pulse_time = get_time_smooth(popt, time, voltage[channel-1], 3)
                        din[f'ch{channel}_time_smooth'].append(pulse_time)
                                            

    return pd.DataFrame(din)