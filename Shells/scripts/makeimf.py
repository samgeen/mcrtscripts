'''
Make an IMF
Sam Geen, March 2016
'''

import random
import numpy as np

import pickle as pik
import os, glob

import customplot
import matplotlib.pyplot as plt

def SampleIMF(Mcluster,seed=None):
    # Cravenly stolen from:
    # https://python4mpia.github.io/fitting_data/MC-sampling-from-Salpeter.html
    # Convert limits from M to logM.
    cachename = None
    if seed is not None:
        # Cache format "caches/cache_mass_?????.dat")
        mstr = str(int(Mcluster)).zfill(7)
        cachestr = str(seed).zfill(5)
        cachename = 'cache/imfs/imf_'+mstr+'_'+cachestr+'.dat' 
        if os.path.exists(cachename):
            f = open(cachename,"rb")
            Mcluster = pik.load(f)
            masses = pik.load(f)
            f.close()
            return np.array(masses)
    print("Making IMF...",)
    mMin = 1.0
    mMax = 120.0
    log_M_Min = np.log(mMin)
    log_M_Max = np.log(mMax)
    alpha = 2.35
    # Since Salpeter SMF decays, maximum likelihood occurs at M_min
    maxlik = mMin**(1.0 - alpha)
    # Prepare array for output masses.
    masses = []
    # Fill in array.
    lastM = 0.0
    while np.sum(masses) < Mcluster:
        # Draw candidate from logM interval.
        logM = random.uniform(log_M_Min,log_M_Max)
        M    = np.exp(logM)
        lastM = M
        # Compute likelihood of candidate from Salpeter SMF.
        likelihood = M**(1.0 - alpha)
        # Accept randomly.
        u = random.uniform(0.0,maxlik)
        if (u < likelihood):
            masses.append(M)
    # Find closest answer to target mass
    if np.abs(np.sum(masses) - Mcluster) > \
       np.abs(np.sum(masses) - M - Mcluster):
        masses.remove(M)
    masses.sort()
    # Save (if a seed is specified)
    if cachename is not None:
        f = open(cachename,"wb")
        pik.dump(Mcluster,f)
        pik.dump(masses,f)
        f.close()
    print("Done")
    return np.array(masses)

def testplot(mcluster):
    masses = SampleIMF(mcluster)
    n = np.arange(0,len(masses))+1.0
    sample = np.arange(0,len(masses),10)
    sample = np.concatenate([sample,[len(masses)-1]])
    
    print(sample)
    dndm = np.diff(n) / np.diff(masses)
    mcent = 0.5*(masses[1:] + masses[:-1])
    sx = dndm
    hist,edges = np.histogram(masses,bins=masses[sample])
    ecent = 0.5*(edges[1:] + edges[:-1])
    x = []
    y = []
    for i in range(0,len(hist)):
        dx = edges[i+1] - edges[i]
        if dx == 0.0:
            dx = 1e-5
            hist[i] = 0.0
        x.append(edges[i])
        x.append(edges[i+1])
        y.append(hist[i]/dx)
        y.append(hist[i]/dx)
    plt.plot(mcent,dndm)
    plt.plot(x,y)
    plt.xlim([1.0,120.0])
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("M / Msun")
    plt.ylabel("dN/dM")
    plt.savefig("testimf.pdf")

if __name__=="__main__":
    testplot(1e3)
