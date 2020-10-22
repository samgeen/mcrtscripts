'''
Sample from the (Chabrier) IMF
Sam Geen, January 2019
'''


import matplotlib as mpl
mpl.use("Agg")

from . import chabrier
import numpy as np
import scipy.interpolate
import scipy.stats

mcen, imfcen, N, M = None,None,None,None

def _SampleChunk(fimf,nsample=10000):
    # METHOD 1 - Monte carlo sample 2D space (issue with low probabilities)
    # Make 2d sample of normalised mass and probability
    #Mmonte = np.random.random(nsample)
    #Pmonte = np.random.random(nsample)
    # If the probability is below the IMF, reject it
    #Pimf = fimf(Mmonte)
    #return Mmonte[Pmonte > Pimf]
    # METHOD 2 - reverse cumulative distribution
    return fimf(np.random.random(nsample))

def SampleIMF(mcluster,mmin=5.0,mmax=120.0,seed=None):
    # mcen - mass in centre of each bin
    # imfcen - dN/dM for each bin
    # N - number-averaged fraction of mass in each bin
    # M - mass-averaged fraction in each bin
    global mcen, imfcen, N, M
    # Seed the random generator
    if seed is not None:
        np.random.seed(seed)
    # Read the IMF and store it in globals
    if mcen is None:
        mcen, imfcen, N, M = chabrier.makeimf()
        imfcen /= imfcen.max() # Normalise to allow probability sampling from 0 to 1
    # Get indices of bins between mmin and mmax
    inds = (mcen >= mmin) * (mcen <= mmax)
    mcenmin = mcen[inds].min()
    mcenmax = mcen[inds].max()
    # Mass of stars to sample
    msample = mcluster * M[inds].sum()
    # Sample stars until we reach msample
    stars = []
    keepsampling = True
    mcennorm = (mcen[inds] - mcenmin)/(mcenmax - mcenmin)
    #fimf = scipy.interpolate.interp1d(mcennorm,imfcen[inds])
    Ncumul = np.cumsum(N)
    Ncumul /= Ncumul.max()
    Ncumul = np.concatenate(([0.0],Ncumul))
    Mtomap = np.concatenate(([0.0],mcen))
    fimf = scipy.interpolate.interp1d(Ncumul,Mtomap)
    #return disn.rvs() * (mcenmax - mcenmin) + mcenmin
    nsample = int(msample)
    while keepsampling:
        # Sample the IMF in chunks
        chunk = _SampleChunk(fimf,nsample)
        #stars = np.concatenate((stars,chunk * (mcenmax - mcenmin) + mcenmin))
        stars = np.concatenate((stars,chunk))
        # Check for end of sampling
        if np.sum(stars) >= mcluster:
            keepsampling = False
            mcumul = np.cumsum(stars)
            stars = np.array(stars).flatten()
            stars = stars[mcumul < mcluster]
    return stars

if __name__=="__main__":
    import matplotlib.pyplot as plt
    # Make sampling test
    mtot = 1e6
    stars = SampleIMF(mtot,5.0,1000.0)
    stars = np.sort(stars)
    counts = np.arange(0,len(stars))
    counts = counts / float(counts.max())
    # Compare to IMF generated
    Ncumul = np.cumsum(N)
    Ncumul /= Ncumul.max()
    # Plot
    plt.clf()
    plt.plot(stars,counts)
    plt.plot(mcen,Ncumul)
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("testsample.pdf")
