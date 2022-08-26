import matplotlib.pyplot as plt
import awkward as ak
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import numpy as np, pandas as pd

def chi2(y,yfit):
    res = 0
    for i in range(len(yfit)):
        res = res + (y[i]-yfit[i])**2/(0.9724)**2 #The denominator should be the unbias Sipm voltage in mV
    return (res / len(yfit))

def waveform(x, C, start, m, end, A, d):
    condlist = [x < start, (x >= start) & (x < end), x >= end]
    funclist = [lambda x: C, lambda x: m*x+(C - m*start), lambda x: (A*(np.exp(d*x)-np.exp(d*end)) + m*(end-start) + C)]
    return np.piecewise(x, condlist, funclist)

def gaus(x, normalization, location, stdev):
    return normalization* np.exp(-0.5*((x - location)/stdev)**2)

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
    chisq = chi2(voltages[0:startindex],waveform(time[0:startindex], *fitparams))
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