'''
Plot the stall radius 
Sam Geen, May 2019
'''

import numpy as np
import scipy.interpolate

import customplot
import matplotlib.pyplot as plt 


def rwri(Cw):
    # Solve the equation
    const = 2**(1.0/3.0) * Cw**(-4.0/3.0)
    roots = np.roots([1.0,0.0,0.0,-1.0,-const])
    # This will be a 4-element array with complex numbers in it
    # rw is the one positive, wholly real value in this array
    roots = roots[np.isreal(roots)]
    root = roots[roots > 0.0][0]
    rirw = root.real
    return 1.0/rirw

def plot():
    plt.clf()
    # Plot uniform case
    npts = 50
    Cws = np.logspace(-1.0,1.0,npts)
    rwris = np.zeros(npts)
    for i, Cw in enumerate(Cws):
        rwris[i] = rwri(Cw)
    plt.plot(Cws,rwris,"k-")
    # Plot stuff
    plt.xlabel("$C_{w}$")
    plt.ylabel("$r_w / r_i$")
    plt.xscale("log")
    plt.savefig("plots/rirw_vs_cw.pdf")

if __name__=="__main__":
    plot()
