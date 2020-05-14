# Implements Chabrier IMF
# Sam Geen, December 2014

import numpy as np
import matplotlib.pyplot as plt

log0p08 = np.log10(0.08)
log0p22 = np.log10(0.22)

def smallfunc(m):
    #rf = 0.158/m * np.exp(-0.5*((np.log10(m) - log0p08)/0.69)**2)
    def rffunc(m):
        return 0.086/m * np.exp(-0.5*((np.log10(m) - log0p22)/0.57)**2)
    k = rffunc(1.0) # Normalise to 1 at M=1
    return rffunc(m)/k

def bigfunc(m):
    return m**-2.3

def makeimf():
    # CHECK NORMALISATION!!!!
    # mcen - mass in centre of each bin
    # imfcen - dN/dM for each bin
    # N - number-averaged fraction of mass in each bin
    # M - mass-averaged fraction in each bin
    mmax = np.log10(120.0)
    small = 10.0**np.arange(-3,0,0.001)
    big = 10.0**np.arange(0,mmax,0.001)
    imf = np.concatenate((smallfunc(small),bigfunc(big)))
    masses = np.concatenate((small,big))
    mcen = 0.5*(masses[1:]+masses[:-1])
    imfcen = 0.5*(imf[1:]+imf[:-1])
    N = np.diff(masses)*0.5*(imfcen)
    k = 1.0/np.trapz(imf*masses,masses)
    #print "k =",k
    M = mcen * N
    M /= np.sum(M)
    N /= np.sum(N)
    return mcen, imfcen, N, M

def findnorm():
    masses,imf,d1,d2 = makeimf()
    k = 1.0/np.trapz(imf*masses,masses)
    #print "knorm =",k
    return k

def mmaxes():
    # CHECK NORMALISATION!!!
    mcen, imf, N, M = makeimf()
    mclusters = 10**np.arange(3,6,0.01)
    mmaxes = mclusters*0
    for i in range(0,len(mclusters)):
        mmaxes[i] = mcen[N*mclusters[i] > 0.5].max()
    plt.plot(mclusters, mmaxes)
    plt.xlabel ("Cluster Mass / M$_{\odot}$")
    plt.ylabel ("Maximum stellar mass in cluster / M$_{\odot}$")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("mmax.png")

def mmaxes_ana():
    k = findnorm()
    mcen, imf, N, M = makeimf()
    mclusters = 10**np.arange(1,4,0.01)
    mlim = 120.0
    mmax = (k/0.3 * mclusters * (2.0**0.3-1))**(1.0/1.3)
    mmax[mmax > mlim] = mlim
    return mclusters, mmax

def plotmmaxes():
    mclusters, mmax = mmaxes_ana()
    plt.clf()
    plt.plot(mclusters,mmax)
    plt.xlabel ("Cluster Mass / M$_{\odot}$")
    plt.ylabel ("Maximum stellar mass in cluster / M$_{\odot}$")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.savefig("mmax_ana.png")
    

def plot():
    mcluster = 1
    mcen, imf, N, M = makeimf()
    N *= mcluster
    #mmax = mcen[N > 0.5].max()
    #print mmax
    plt.clf()
    #plt.plot(mcen,N)
    #plt.plot(mcen,mcen*0+0.5,"k--")
    #plt.plot(N*0+mmax,N,"k--")
    plt.plot(mcen,imf)
    plt.xscale("log")
    plt.yscale("log")
    #plt.ylabel("N*"+str(mcluster))
    plt.ylabel("dN/dM")
    plt.xlabel("M/M$_{\odot}$")
    plt.savefig("chabrierimf.png")


if __name__=="__main__":
    plot()
    plotmmaxes()
    mmaxes_ana()
