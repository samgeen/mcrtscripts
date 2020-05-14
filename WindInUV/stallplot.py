'''
Plot the stall radius 
Sam Geen, May 2019
'''

import numpy as np
import scipy.interpolate

import customplot
import matplotlib.pyplot as plt 


def integrate(r0, w):
    nvals = 10000
    dt = 0.001
    rs = np.zeros(nvals+1)
    ts = np.arange(0.0,dt*nvals+dt/2.0,dt)
    r = r0
    rs[0] = r
    index = (2.0*w-3.0) * 0.25
    for i in range(0,nvals):
        v0 = r**(1.0-0.5*w)
        drdt = (r**index - 1.0)*v0
        r += drdt * dt
        if r < 0:
            rs = rs[0:i]
            ts = ts[0:i]
            break
        rs[i+1] = r
    return ts, rs

def plot(noiso = False):
    plt.clf()
    # Plot uniform case
    ts, rs = integrate(r0=0.1,w=0.0)
    plt.plot(ts,rs,"k--",label="Uniform (w=0)")
    ts, rs = integrate(r0=2.0,w=0.0)
    plt.plot(ts,rs,"k--")
    # Plot isothermal case
    col = "k"
    if noiso:
        col = "w"
    ts, rs = integrate(r0=1.1,w=2.0)
    plt.plot(ts,rs,col+"-",label="Isothermal (w=2)")
    ts, rs = integrate(r0=0.9,w=2.0)
    plt.plot(ts,rs,col+"-")
    plt.plot(ts,rs*0.0+1.0,"k-",alpha=0.25)
    # Plot stuff
    plt.xlabel("t")
    plt.ylabel("R")
    plt.legend(frameon=False,loc="lower left")
    if not noiso:  
        plt.savefig("plots/stallplot.pdf")
    else:
        plt.savefig("plots/stallplot_noiso.pdf")

if __name__=="__main__":
    plot()
    plot(noiso=True)
