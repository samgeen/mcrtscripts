'''
Make stellar sources and cache them
Sam Geen, August 2018
'''

import os, sys
sys.path.append("/home/stgeen0/StellarSources")
sys.path.append("/home/stgeen0/StellarSources/Fortran/f2py")

import numpy as np
import scipy.interpolate

import cPickle as pik

from Input import readgeneva, derivedgeneva, readstarburst, readsbwinds

import singlestar, makeimf
from IMFUtils import sampleimf

from consts import *

# Verbose mode?
verbose = True

# Stellar track stuff
windtracks = None
trackloc = "/home/stgeen0/StellarSources/"
starmasses = np.arange(5,121,5)
starmetals = [0.002,0.014]

# Set up single star module
startableloc = "/home/stgeen0/StellarSources/Compressed/singlestar_z0.014"
print "Reading star tables in", startableloc
singlestar.star_setup(startableloc)  

# Cache folder
cachefolder = "/home/stgeen0/MCRT/mcrtscripts/WindInUV/"

# IMF function for sampling stellar clusters
# Runs stars through a nearest 5 Msun sampling
def IMF(clustermass,seed=None):
    masses = sampleimf.SampleIMF(clustermass,mmin=5.0,mmax=120.0,seed=seed)
    # Make sure we don't include small stars
    masses = masses[masses > 5.0]
    # Round to the nearest 5
    masses = (((masses/5.0 + 0.5).astype(int))*5.0).astype(float)
    return masses

# IMF function for sampling stellar clusters
# Runs stars through a nearest 5 Msun sampling
def IMF_OLD(clustermass,seed=None):
    if seed is not None:
        mstr = str(int(clustermass)).zfill(7)
        cachestr = str(seed).zfill(5)
        cachename = cachefolder+"/cache/imfs_massiveonly/imfmo_"+mstr+"_"+cachestr+".dat"
        if os.path.exists(cachename):
            f = open(cachename,"rb")
            masses = pik.load(f)
            f.close()
            return np.array(masses)
    # Is clustermass an iterable or a single value?
    try:
        dummy = iter(clustermass)
    except TypeError:
        # Single value: make an IMF up to clustermass
        masses = makeimf.SampleIMF(clustermass,seed=seed)
    else:
        # Multiple values: you got your cluster there already
        masses = np.array(clustermass)
    masses = masses[masses > 5.0]
    # Round to the nearest 5
    masses = (((masses/5.0 + 0.5).astype(int))*5.0).astype(float)
    # Save (if a seed is specified)
    if cachename is not None:
        f = open(cachename,"wb")
        pik.dump(masses,f)
        f.close()
    return masses

# Object to be cached (DON'T INSTANTIATE OUTSIDE THIS MODULE! USE Star BELOW)
class _Star(object):
    def __init__(self,mass,metal):
        global windtracks
        if windtracks is None:
            windtracks = readgeneva.Tracks(trackloc)
        self._mass = mass
        self._metal = metal
        self._iscluster = False
        # Wind stuff
        windtrack = windtracks.FindTrack(mass,metal,rotating=True)
        self._windtimes = windtrack["Time"]*yrins
        self._windlum = np.array(derivedgeneva.WindLuminosity(windtrack))
        self._windmom = np.array(derivedgeneva.WindMomentumRate(windtrack))
        self._Teff = 10.0**np.array(windtrack["logTe"])
        self._windlumFunc = None
        self._windmomFunc = None
        self._TeffFunc = None
        # Spectral stuff
        spectrack = readstarburst.Track(mass,metal,trackloc)
        self._spectimes = spectrack.Times()*yrins # in s
        self._nphotons = np.array(spectrack.NPhotonsVsTime(13.6, 1000.0))
        self._ephotons = np.array(spectrack.EPhotonsVsTime(13.6, 1000.0))
        self._nphotonFunc = None
        self._ephotonFunc = None

    @property
    def mass(self):
        return self._mass

    @property
    def metal(self):
        return self._metal

    def Copy(self):
        return Star(self._mass, self._metal)

    def WindLuminosity(self, t):
        if self._windlumFunc is None:
            self._windlumFunc = scipy.interpolate.interp1d(self._windtimes,self._windlum,fill_value="extrapolate")
        return self._windlumFunc(t)

    def WindMomentumRate(self, t):
        if self._windmomFunc is None:
            self._windmomFunc = scipy.interpolate.interp1d(self._windtimes,self._windmom,fill_value="extrapolate")
        return self._windmomFunc(t)
    
    def Teff(self,t):
        if self._TeffFunc is None:
            self._TeffFunc = scipy.interpolate.interp1d(self._windtimes,self._Teff,fill_value="extrapolate")
        return self._TeffFunc(t)

    def EPhoton(self, t):
        if self._ephotonFunc is None:
            self._ephotonFunc = scipy.interpolate.interp1d(self._spectimes,self._ephotons/self._nphotons,fill_value="extrapolate")
        return self._ephotonFunc(t)

    def NPhotons(self, t):
        if self._nphotonFunc is None:        
            self._nphotonFunc = scipy.interpolate.interp1d(self._spectimes,self._nphotons,fill_value="extrapolate")
        return self._nphotonFunc(t)

    def Lifetime(self):
        # Use the f2py module that reads the processed tables for this
        return singlestar.star_lifetime(self._mass)
        
class _ClusterStarEmission(object):
    '''
    Intermediate cluster object that collects stellar emissions for faster cluster output
    '''        
    def __init__(self,metal):
        self._metal = metal
        # Dictionary of stars with mass as the key and star object as the value
        self._stars = AllStars(metal)
        if len(self._stars) > 0:
            self._lifetime = np.max([singlestar.star_lifetime(star.mass) for star in self._stars.values()])
        else:
            self._lifetime = 0.0
        # Cache cluster information
        self._currcluster = None
        self._counts = {}

    @property
    def metal(self):
        return self._metal

    def _CheckCluster(self,cluster):
        if self._currcluster is not cluster:
            self._counts = {}
            for starmass in self._stars.iterkeys():
                self._counts[starmass] = np.count_nonzero(cluster.starmasses==starmass)
            self._currcluster = cluster

    def _SumForFunc(self,cluster,t,func):
        self._CheckCluster(cluster)


    def WindLuminosity(self, cluster, t):
        self._CheckCluster(cluster)
        val = 0.0
        for starmass, star in self._stars.iteritems():
            val += self._counts[starmass]*star.WindLuminosity(t)
        return val

    def WindMomentumRate(self, cluster, t):
        self._CheckCluster(cluster)
        val = 0.0
        for starmass, star in self._stars.iteritems():
            val += self._counts[starmass]*star.WindMomentumRate(t)
        return val

    def Teff(self,cluster,t):
        self._CheckCluster(cluster)
        val = 0.0
        # Find the maximum Teff, which will set the effective temperature of the cluster
        for starmass, star in self._stars.iteritems():
            val = max(val,star.Teff(t))
        return val

    def EPhoton(self, cluster, t):
        # NOTE!!! WE NEED THE AVERAGE PHOTON ENERGY HERE, NOT THE TOTAL
        self._CheckCluster(cluster)
        eval = 0.0
        nval = 0.0
        for starmass, star in self._stars.iteritems():
            nphoton = self._counts[starmass]*star.NPhotons(t)
            eval += star.EPhoton(t)*nphoton
            nval += nphoton
        if nval > 0:
            return eval/nval
        else:
            return 0.0

    def NPhotons(self, cluster, t):
        self._CheckCluster(cluster)
        val = 0.0
        for starmass, star in self._stars.iteritems():
            val += self._counts[starmass]*star.NPhotons(t)
        return val

_clusterStarEmissions = {}
def ClusterStarEmission(metal):
    if not metal in _clusterStarEmissions:
        _clusterStarEmissions[metal] = _ClusterStarEmission(metal)
    return _clusterStarEmissions[metal]


class ClusterOnTheFly(object):
    # Better to use if there are many clusters to generate
    def __init__(self,mass,metal,seed=1,maxtime_Myr = 130.0,dt_Myr=0.01):
        # NOTE: mass can be one mass or a list
        # if mass is one value, populate a cluster with these masses
        # if mass is multiple values, use this directly as a cluster
        self._mass = mass
        self._metal = metal
        self._starmasses = np.atleast_1d(IMF(mass,seed))

    @property
    def mass(self):
        return self._mass

    @property
    def metal(self):
        return self._metal

    @property
    def starmasses(self):
        return self._starmasses

    def WindLuminosity(self, t):
        return ClusterStarEmission(self._metal).WindLuminosity(self,t)

    def WindMomentumRate(self, t):
        return ClusterStarEmission(self._metal).WindMomentumRate(self,t)

    def Teff(self,t):
        return ClusterStarEmission(self._metal).Teff(self,t)

    def EPhoton(self, t):
        return ClusterStarEmission(self._metal).EPhoton(self,t)

    def NPhotons(self, t):
        return ClusterStarEmission(self._metal).NPhotons(self,t)

    def Lifetime(self):
        # Return maximum stellar lifetime
        #return self._lifetime
        # TODO: Figure out a clean, fast way to do this
        raise NotImplementedError

class ClusterPrecomputed(object):
    # Better to use if few clusters are needed and they are sampled many times
    def __init__(self,mass,metal,seed=1,maxtime_Myr = 130.0,dt_Myr=0.01):
        # NOTE: mass can be one mass or a list
        # if mass is one value, populate a cluster with these masses
        # if mass is multiple values, use this directly as a cluster
        self._mass = mass
        self._metal = metal
        self._masses = np.atleast_1d(IMF(mass,seed))
        #self._stars = [Star(mass,metal) for mass in self._masses]
        if len(self._stars) > 0:
            self._lifetime = np.max([singlestar.star_lifetime(star.mass) for star in self._stars])
        else:
            self._lifetime = 0.0
        # Construct new tables by summing over the old tables
        self._times = np.linspace(0.0,maxtime_Myr,int(maxtime_Myr/dt_Myr)+1)*yrins*1e6
        self._windlum = np.zeros(len(self._times))
        self._windmom = np.zeros(len(self._times))
        self._Teff = np.zeros(len(self._times))
        self._nphotons = np.zeros(len(self._times))
        self._ephotons = np.zeros(len(self._times))
        for star in self._stars:
            self._windlum += star.WindLuminosity(self._times)
            self._windmom += star.WindMomentumRate(self._times)
            self._Teff = np.maximum(self._Teff, star.Teff(self._times))
            nphotons = star.NPhotons(self._times)
            self._nphotons += nphotons
            self._ephotons += star.EPhoton(self._times)*nphotons
        self._windlumFunc = None
        self._windmomFunc = None
        self._TeffFunc = None
        self._nphotonFunc = None
        self._ephotonFunc = None   

    @property
    def mass(self):
        return self._mass

    @property
    def metal(self):
        return self._metal

    @property
    def starmasses(self):
        return self._starmasses

    def WindLuminosity(self, t):
        if self._windlumFunc is None:
            self._windlumFunc = scipy.interpolate.interp1d(self._times,self._windlum,fill_value="extrapolate")
        return self._windlumFunc(t)

    def WindMomentumRate(self, t):
        if self._windmomFunc is None:
            self._windmomFunc = scipy.interpolate.interp1d(self._times,self._windmom,fill_value="extrapolate")
        return self._windmomFunc(t)

    def Teff(self,t):
        if self._TeffFunc is None:
            self._TeffFunc = scipy.interpolate.interp1d(self._windtimes,self._Teff,fill_value="extrapolate")
        return self._TeffFunc(t)

    def EPhoton(self, t):
        if self._ephotonFunc is None:
            self._ephotonFunc = scipy.interpolate.interp1d(self._times,self._ephotons/self._nphotons,fill_value="extrapolate")
        return self._ephotonFunc(t)

    def NPhotons(self, t):
        if self._nphotonFunc is None:        
            self._nphotonFunc = scipy.interpolate.interp1d(self._times,self._nphotons,fill_value="extrapolate")
        return self._nphotonFunc(t)

    def Lifetime(self):
        # Return maximum stellar lifetime
        #return self._lifetime
        # TODO: Figure out a clean, fast way to do this
        raise NotImplementedError

class Binary(object):
    # Two equal mass stars
    def __init__(self,mass,metal):
        self._mass = mass
        self._metal = metal
        self._star = Star(mass,metal)

    def WindLuminosity(self, t):
        return 2.0*self._star.WindLuminosity(t)

    def WindMomentumRate(self, t):
        return 2.0*self._star.WindMomentumRate(t)

    def Teff(self,t):
        return self._star.Teff(t)

    def EPhoton(self, t):
        return 2.0*self._star.EPhoton(t)

    def NPhotons(self, t):
        return 2.0*self._star.NPhotons(t)

    def Lifetime(self):
        # Return maximum stellar lifetime
        return self._star.Lifetime(t)

# A factory method with a caching wrapper
starcaches = {}
def Star(mass, metal):
    # Set up star cache
    if not metal in starcaches:
        starcaches[metal] = {}
    # Make star or load from cache
    if not mass in starcaches[metal]:
        cachename = cachefolder+"/cache/stars/star_M"+str(int(mass))+"_Z"+str(metal)+".pik"
        if os.path.exists(cachename):
            if verbose:
                print "Loading star, M =", mass, "Z =",metal, "from cache"
            f = open(cachename,"rb")
            star = pik.load(f)
            f.close()
        else:
            if verbose:
                print "Making and saving star, M =", mass, "Z =",metal, "to cache"
            star = _Star(mass, metal)
            f = open(cachename,"wb")
            pik.dump(star,f)
            f.close()
        starcaches[metal][mass] = star
    # Return the star object saved in the cache
    return starcaches[metal][mass]

starlist = None
def AllStars(metal):
    # some
    global starlist
    if starlist is None:
        starlist = {}
        for mass in starmasses:
            starlist[mass] = Star(mass,metal)
    # BODY
    return starlist

def test():
    # Brief test of the caching system
    star = Star(60,0.014)
    star2 = Star(60,0.014)
    print star2.WindLuminosity(1e6*yrins)

def populate():
    # Run through the parameter space and make cached stars
    for mass in starmasses:
        for metal in starmetals:
            star = Star(mass,metal)

if __name__=="__main__":
    populate()
