'''
Plot the stall radius 
Sam Geen, May 2019
'''

import numpy as np
import scipy.interpolate

from consts import *

import customplot
import matplotlib.pyplot as plt 


def finddensitiesOLD(radii,ninner,Arp,QH):
    # NOTE: ninner balances the wind pressure
    # Make up some numbers
    nwind = 0.01 * ninner
    alphaB = 2e-13
    rw = 0.5 * radii[-1]
    drw = 0.25 * radii[-1]
    # Solve the equation
    for i in range(0,100):
        densities = ninner / (1 - Arp*(radii-rw)*ninner)
        dr = radii[1] - radii[0]
        recomb = np.sum(4.0*np.pi*radii[radii > rw]**2.0 * dr * densities[radii > rw]**2 * alphaB)
        if recomb < QH:
            rw += drw
        else:
            rw -= drw
        rw = max(rw,radii[0])
        rw = min(rw,radii[-1])
        drw *= 0.5
        print recomb/QH, rw/radii[-1]
    densities[radii <= rw] = nwind
    return densities, rw

def finddensities(radii,ninner,Arp,QH,rwvsri = 0.25):
    # NOTE: ninner balances the wind pressure
    # Make up some numbers
    nwind = 0.01 * ninner
    alphaB = 2e-13
    rw = rwvsri * radii[-1]
    # Solve the equation
    densities = ninner / (1 - Arp*(radii-rw)*ninner)
    dr = radii[1] - radii[0]
    recomb = np.sum(4.0*np.pi*radii[radii > rw]**2.0 * dr * densities[radii > rw]**2 * alphaB)
    print recomb/QH, rw/radii[-1]
    densities[radii <= rw] = nwind
    return densities, rw

def plot():
    plt.clf()
    radii = np.linspace(0.0,1.0,100)*pc
    ninner = 100.0
    Arp = 0.5 / ninner / radii[50]
    QH = 5e47
    #densities, rw = finddensities(radii,ninner,Arp,QH)
    #print rw, densities
    # Plot uniform case (no RP)
    densities, rw = finddensities(radii,ninner,0.0,QH)
    plt.plot(radii/pc,densities,"k--",label="No radiation pressure")
    densities, rw = finddensities(radii,0.7*ninner,Arp,QH,rwvsri=0.4)
    plt.plot(radii/pc,densities,"k-",label="Radiation pressure")
    # Plot stuff
    plt.xlabel("$r$ / pc")
    plt.ylabel("$n_i(r)$ / cm$^{-3}$")
    plt.legend(frameon=False,loc="lower right",fontsize="medium")
    #plt.yscale("log")
    plt.savefig("plots/densitygradient_with_rp.pdf")

if __name__=="__main__":
    plot()
