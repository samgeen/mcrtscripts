'''
Plot the equilibrium temperature of the gas
Sam Geen, April 2020
'''

import os, sys
import numpy as np
import matplotlib.pyplot as plt

import pickle as pik

import scipy.interpolate

frig=False

if frig:
    sys.path.append("../Cooling/old")
    import cooling_module as cooling
else:
    sys.path.append("../Cooling/")
    import cooling


Zsolar = 0.014
kB = 1.380649e-16
X = 0.74
mH = 1.66e-24
gamma = 5.0/3.0

def runforZ(metal):
    ns = np.logspace(-1,10,110)
    #print(ns)
    T2 = ns*0.0 + 1e2 # Initial temperature
    dT2 = T2+0.0
    zsolars = ns*0.0 + metal/Zsolar
    dt = 1e9
    ncell = len(ns)
    print("Doing...")
    loops = 0
    while np.abs(dT2/T2).max() > 0.001:
        if frig:
            dT2 = cooling.solve_cooling_frig(ns,T2,zsolars,dt,gamma,ncell)
        else:
            dT2 = cooling.solve_cooling_ramses(ns,T2,zsolars,dt,gamma,ncell)
        T2 += dT2
        #print(dT2, T2)
        loops += 1
    print("Done!", loops)
    return ns, T2

def runforZCached(metal):
    cachefilename = "../cache/runforZ_Z"+str(metal)+".pik"
    if not os.path.exists(cachefilename):
        ns, T2 = runforZ(metal)
        f = open(cachefilename,"wb")
        pik.dump(ns, f)
        pik.dump(T2, f)
        f.close()
    else:
        f = open(cachefilename,"rb")
        ns = pik.load(f)
        T2 = pik.load(f)
        f.close()
    return ns, T2


_coolingdict = {}
_coolingbypressuredict = {}

def NeutralTemperature(nH,metal):
    if not metal in _coolingdict:
        ns, T2 = runforZCached(metal)    
        finterp = scipy.interpolate.interp1d(nH,T2)
        _coolingdict[metal] = finterp
    return _coolingdict[metal](nH)

def NeutralTemperatureFromPressure(Pressure,metal):
    if not metal in _coolingbypressuredict:
        ns, T2 = runforZCached(metal)   
        # Calculate sound speed squared
        cs2 = gamma * kB * T2 * X / mH
        Pshell = mH/X*ns*cs2
        finterp = scipy.interpolate.interp1d(Pshell,T2,bounds_error=False,fill_value=(T2[0],T2[-1]))
        _coolingbypressuredict[metal] = finterp
    return _coolingbypressuredict[metal](Pressure)

def run(plotpressure=False):
    plt.clf()
    for metal in [0.001,0.01,0.1,1.0]:
        n, T = runforZ(metal*Zsolar)
        #print(metal, n, T)
        if not plotpressure:
            y = T
            ylabel = "T / K"
            plotname = "temperature"
        else:
            y = 1.5*n*kB*T
            ylabel = "Thermal pressure / erg / cm$^{-3}$"
            plotname = "pressure"
        plt.plot(n,y,label="Z/Z$_{\odot}$ = "+str(metal))
    plt.legend(frameon=False)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("$n_H$ / cm$^{-3}$")
    plt.ylabel(ylabel)
    if plotname == "temperature":
        plt.ylim([5.0,2000.0])
    coolmod = "ramses"
    if frig:
        coolmod = "frig"
    plt.savefig("../plots/equilibrium"+plotname+"_"+coolmod+".pdf")

if __name__=="__main__":
    run()
    run(True)