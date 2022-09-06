import numpy as np

def chi2(y,yfit):
    res = 0
    for i in range(len(yfit)):
        res = res + (y[i]-yfit[i])**2/(0.9724)**2 #The denominator should be the unbias Sipm voltage in mV
    return (res / len(yfit))

def waveform(x, Ci, start, m, end, d, Cf):
    condlist = [x < start, (x >= start) & (x <= end), x > end]
    funclist = [lambda x: Ci, lambda x: m*x+(Ci - m*start), lambda x: (m*(end-start) + Ci)*np.exp((x-end)/d) + Cf]
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
        if times[i] < fitparams[1]: startindex = i
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
    voltage_smooth = []
    for i in range(len(voltages) - degree - 1):
        vol = 0
        for j in range(degree):
            vol += voltages[j+i]
        voltage_smooth.append(vol / degree)
    for i in range(int((degree - 1) / 2)):
        voltage_smooth.append(voltages[-(int((degree - 1) / 2) - i)])
        voltage_smooth.insert(0,voltages[int((degree - 1) / 2) - i - 1])

    halfamp = 0.5*(max(voltage_smooth)-fitparams[0]) + fitparams[0]
    halfindex = findindex(voltage_smooth,max(voltage_smooth)+1)
    haltime = 0
    while voltage_smooth[halfindex] - halfamp > 0:
        halfindex -= 1
    if abs(voltage_smooth[halfindex] - halfamp) < abs(voltage_smooth[halfindex + 1] - halfamp):
        halftime = times[halfindex]
    else: 
        halftime = times[halfindex + 1]
    return halftime

def voltage_smooth(voltages, degree):
    voltage = []
    for i in range(len(voltages) - degree + 1):
        vol = 0
        for j in range(degree):
            vol += voltages[j+i]
        voltage.append(vol / degree)
    for i in range(int((degree - 1) / 2)):
        voltage.append(voltages[-(int((degree - 1) / 2) - i)])
        voltage.insert(0,voltages[int((degree - 1) / 2) - i - 1])
    return voltage

