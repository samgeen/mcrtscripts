'''
Find dense clumps in the simulation
Sam Geen, March 2015
'''

import customplot
import matplotlib.pyplot as plt

import Hamu

import numpy as np

import os, errno

import pymses
from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import outflowmodel, linestyles

class ClumpFinder(object):
    def __init__(self, snap, nthresh=1e5):
        self._snap = snap
        self._nthresh = nthresh

    def Run(self):
        amr = self._snap.amr_source(["rho","B-left","B-right"])
        cell_source = CellsToPoints(amr)
        cells = cell_source.flatten()
        clumpcells = self._FilterByDensity(cells)
        clumps = self._FindClumps(clumpcells)
        return clumps
        
    def _FilterByDensity(self, cells):
        rhos = cells["rho"]*self._snap.info["unit_density"].express(C.H_cc)
        return cells.filtered_by_mask(rhos > self._nthresh)

    def _FindClumps(self, cells):
        '''
        Finds the nearest neighbours of the clumps
        '''
        clumps = []
        while len(cells["rho"]) > 0:
            print "Number of cells left:", len(cells["rho"])
            print "CLUMP FINDING..."
            # Find any clump in the list of cells
            clump = ClumpMaker(self._snap, cells)
            clumps.append(clump)
            # Remove cells not in clump
            cells = cells.filtered_by_mask(clump.inclump == 0)
            print "CLUMP FOUND!"
        print "DONE BUT NOW FIND CLUMP PROPERTIES"
        return clumps

class ClumpMaker(object):
    '''
    Identifies a single clump and returns its stats
    '''
    def __init__(self, snap, cells):
        self._snap = snap
        self._ncells = len(cells["rho"])
        self._cells = cells
        self._inclump = np.zeros((self._ncells),dtype="int32")
        self._radii = cells.get_sizes()*2.0 # Make the radii "big enough"
        self._setup = False
        self._mass = 0.0
        self._peakdens = 0.0
        self._centre = 0.0
        self._radius = 0.0
        self._distance = 0.0
        self._bfactor = 0.0

    def _Setup(self):
        if not self._setup:
            self._FindClump()
            self._FindClumpProperties()
            self._setup = True
            

    def _FindClump(self):
        newfound = np.array([0]) # Pick the first cell
        self._inclump[0] = 1
        # Loop while newly found cells exist
        icount = 0
        while len(newfound) > 0:
            icount += 1
            # Make new newfound list and search through previous list
            oldfound = newfound
            newfound = np.array([],dtype="int32")
            # Find new neighbours
            for icell in oldfound:
                found = self._FindNeighbours(icell)
                newfound = np.concatenate((newfound,found))
            # Uniquify newfound
            newfound = np.unique(newfound)
            # Identify these as being in the clump
            self._inclump[newfound] = 1
        print "Num cells in clump:", len(self._inclump[self._inclump == 1])

    def _FindClumpProperties(self):
        cells = self._cells.filtered_by_mask(self._inclump == 1)
        munit = self._snap.info["unit_mass"].express(C.Msun)
        dunit = self._snap.info["unit_density"].express(C.H_cc)
        runit = self._snap.info["unit_length"].express(C.pc)
        magunit = self._snap.info["unit_mag"].express(C.Gauss)
        rhos = cells["rho"]*dunit
        self._mass = np.sum(cells["rho"]*cells.get_sizes()**3.0*munit)
        self._peakdens = rhos.max()
        central = rhos == self._peakdens
        cells.points -= 0.5 # centre of the volume
        cells.points *= runit
        self._centre = cells.points[central,:]+0.0
        self._distance = np.sqrt(np.sum(self._centre**2))
        print "CENTRE", self._centre
        points = cells.points - self._centre
        self._radius = np.sqrt(np.sum(points**2,1).max())
        # 1e6 because Gauss -> microgauss
        magnetic = 0.5*(cells["B-left"]+cells["B-right"])*magunit*1e6
        # b from Bertoldi & McKee 1990
        self._bfactor = (np.sqrt(np.sum(magnetic**2,1))*rhos**(-0.5))[central]
        
    @property
    def inclump(self):
        self._Setup()
        return self._inclump

    @property
    def centre(self):
        self._Setup()
        return self._centre

    @property
    def distance(self):
        self._Setup()
        return self._distance

    @property
    def peakdensity(self):
        self._Setup()
        return self._peakdens

    @property
    def radius(self):
        self._Setup()
        return self._radius

    @property
    def mass(self):
        self._Setup()
        return self._mass

    @property
    def bfactor(self):
        self._Setup()
        return self._bfactor
            
    def _FindNeighbours(self,icell):
        '''
        Find neighbouring cells
        TODO: TIDY UP INTERFACE
        '''
        cells = self._cells
        radius = self._radii[icell]
        points = cells.points - cells.points[icell,:]
        # Find points that are inside the radius and not already in the clump
        nbors = (np.sum(points**2,1) < radius**2) * (self._inclump == 0)
        return np.where(nbors)[0]
    
def ClumpProperties(clumps,simname):
    path = "../plots/clumps/"+simname
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
    print "Number of clumps:", len(clumps)
    masses = np.array([clump.mass for clump in clumps])
    PlotProperty(masses,"Mass","Msun",simname)
    radii = np.array([clump.radius for clump in clumps])
    PlotProperty(radii,"R","pc",simname)
    distance = np.array([clump.distance for clump in clumps])
    PlotProperty(distance,"Distance","pc",simname)
    Plot2Props(masses,distance,"Mass","Msun","Distance","pc",simname)

def PlotProperty(prop,label,unit,simname):
    prop.sort()
    nums = range(0,len(prop))
    plt.clf()
    plt.plot(prop,nums)
    plt.xlabel(label+" / "+unit)
    plt.ylabel("N(>"+label+")")
    plt.xscale("log")
    plt.savefig("../plots/clumps/"+simname+"/"+label+".pdf") 

def Plot2Props(prop1,prop2,label1,unit1,label2,unit2,simname):
    prop1.sort()
    prop2.sort()
    plt.clf()
    plt.plot(prop1,prop2)
    plt.xlabel(label1+" / "+unit1)
    plt.ylabel(label2+" / "+unit2)
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("../plots/clumps/"+simname+"/"+label1+"vs"+label2+".pdf")     

def FindClumps(snap, nthresh=1e6):
    clumps = ClumpFinder(snap,nthresh)
    return clumps.Run()

FindClumpsHamu = Hamu.Algorithm(FindClumps)

def CompareMassDistance(simnames,time=1.25,folder="photons",suffix=""):
    myr = None
    masses = {}
    radii = {}
    distances = {}
    bfactors = {}
    plt.clf()
    sims = {}
    colins = ["k","r","m","b"]
    cols = {}
    icol = 0
    try:
        os.mkdir("../plots/clumps/"+folder)
    except:
        pass
    for simname in simnames:
        cols[simname] = colins[icol]
        icol += 1
    noclumps = True
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        sims[simname] = sim
        if myr is None:
            snap = sim.Snapshots()[0]
            myr = snap.RawData().info["unit_time"].express(C.Myr)
        snap = sim.FindAtTime(time/myr)
        clumps = FindClumpsHamu(snap)
        masses[simname] = np.array([clump.mass for clump in clumps])
        radii[simname] = np.array([clump.radius for clump in clumps])
        distances[simname] = np.array([clump.distance for clump in clumps])
        bfactors[simname] = np.array([clump.bfactor for clump in clumps])
        sorter = np.argsort(distances[simname])
        masses[simname] = masses[simname][sorter]
        radii[simname] = radii[simname][sorter]
        distances[simname] = distances[simname][sorter]
        bfactors[simname] = bfactors[simname][sorter]
        if len(clumps) > 0:
            noclumps = False
    if noclumps:
        print "NO CLUMPS AT THIS TIME! RETURNING"
        return
    # Plot mass vs distance
    mdsyms = []
    labels = []
    for simname in simnames:
        di = np.log10(distances[simname])
        ma = np.log10(masses[simname])
        fit = np.polyfit(di,ma,1)
        print "mass, distance fit (coeff, power):", fit[1], fit[0]
        #plt.plot(np.log10(distances[simname]),
        #         np.log10(masses[simname]),
        #         label=sims[simname].Label())
        rscale = np.log10(radii[simname])
        rscale = 19*(rscale - rscale.min())/(rscale.max() - rscale.min())+1
        areas = np.pi*(rscale)**2
        mdsym = plt.scatter(di, ma,
                    #s=areas,
                    alpha=0.7,
                    label=sims[simname].Label(),
                    color=cols[simname])
        labels.append(sims[simname].Label())
        mdsyms.append(mdsym)
    if len(simnames) > 1:
        plt.legend(mdsyms,labels,scatterpoints=1,
                   fontsize="x-small",loc="lower right")
    plt.xlabel("log(Distance to Source / pc)")
    plt.ylabel("log(Clump Mass / M$_{\odot}$)")
    plt.savefig("../plots/clumps/"+folder+"/massdistance"+suffix+".pdf")
    # Plot radius vs distance
    plt.clf()
    for simname in simnames:
        #plt.plot(np.log10(distances[simname]),
        #         np.log10(radii[simname]),
        #         label=sims[simname].Label())
        mdsym = plt.scatter(np.log10(distances[simname]),
                    np.log10(radii[simname]),
                    alpha=0.7,
                    label=sims[simname].Label(),
                    color=cols[simname])
        labels.append(sims[simname].Label())
        mdsyms.append(mdsym)
    if len(simnames) > 1:
        plt.legend(mdsyms,labels,scatterpoints=1,
                   fontsize="x-small",loc="lower right")

    #plt.legend(fontsize="xx-small",loc="lower right")
    plt.xlabel("log(Distance to Source / pc)")
    plt.ylabel("log(Clump Radius / pc)")
    plt.savefig("../plots/clumps/"+folder+"/radiusdistance"+suffix+".pdf")
    # Plot mass vs radius
    plt.clf()
    for simname in simnames:
        sorter = np.argsort(masses[simname])
        masses[simname] = masses[simname][sorter]
        radii[simname] = radii[simname][sorter]
        bfactors[simname] = bfactors[simname][sorter]
        #plt.plot(np.log10(masses[simname]),
        #         np.log10(radii[simname]),
        #         label=sims[simname].Label())
        mdsym = plt.scatter(np.log10(radii[simname]),
                            np.log10(masses[simname]),
                            alpha=0.7,
                            label=sims[simname].Label(),
                            color=cols[simname])
        labels.append(sims[simname].Label())
        mdsyms.append(mdsym)
    if len(simnames) > 1:
        plt.legend(mdsyms,labels,scatterpoints=1,
                   fontsize="x-small",loc="lower right")
    #plt.legend(fontsize="xx-small",loc="lower right")
    plt.xlabel("log(Clump Radius / pc)")
    plt.ylabel("log(Clump Mass / M$_{\odot}$)")
    plt.savefig("../plots/clumps/"+folder+"/radiusmass"+suffix+".pdf")
    # Plot mass vs bfactor
    plt.clf()
    for simname in simnames:
        plt.plot(np.log10(masses[simname]),
                 np.log10(bfactors[simname]),
                 label=sims[simname].Label())
    plt.legend(fontsize="x-small",loc="lower right")
    plt.xlabel("Clump Mass / M$_{\odot}$")
    plt.ylabel("b $\mu G$")
    plt.savefig("../plots/clumps/"+folder+"/massbfactor"+suffix+".pdf")

def EvaporationRate(simname,time=1.25,mag=False,tevout=False):
    myr = None
    masses = {}
    radii = {}
    distances = {}
    bfactors = {}
    sims = {}
    sim = Hamu.Simulation(simname)
    nthresh = 1e6 #3e5
    if "_C2" in simname:
        nthresh = 1e7
    if "_C2" in simname:
        nthresh = 1e8
    if myr is None:
        snap = sim.Snapshots()[0]
        myr = snap.RawData().info["unit_time"].express(C.Myr)
    snap = sim.FindAtTime(time/myr)
    clumps = FindClumpsHamu(snap,nthresh=nthresh)
    times = np.arange(0,5,0.1)
    mtots = np.zeros((len(times)))
    flux = 10.0**float(simname[1:3]) # "N48_M4_B02"
    S49 = flux/1e49
    phi1 = 1.0
    c5 = 0.545
    tevs = []
    for clump in clumps:
        m1 = clump.mass
        R1 = clump.distance
        n03 = clump.peakdensity/1e3
        b = clump.bfactor
        # Equation 4.10a
        if not mag:
            tev = 4.48e5 * phi1 * c5 ** (-1.2) * S49 ** (-0.2) * \
                R1 ** 0.2 * m1 ** 0.4 / 1e6 # in Myr
        else:
            tev = 1.30e5 * phi1 * (n03**(1.0/6.0)/b)**(6.0/7.0) * \
                S49 ** (-2.0/7.0) * R1 ** (4.0/7.0) * m1 ** (3.0/7.0) / 1e6
        tevs.append(tev)
        mclump = m1 * (1.0 - times / tev)
        mclump[mclump < 0.0] = 0.0
        mtots += mclump
    if tevout:
        return clumps, tevs
    return times, mtots

def AccretionEvaporation(sim,tff,tstart=1.25,mag=False,nofluxsim=None):
    '''
    Use slow accretion model + evaporation
    '''
    myr = None
    masses = {}
    radii = {}
    distances = {}
    bfactors = {}
    sims = {}
    simname = sim.Name()
    outflowmodel.tstart = 1.25 # HACK!!
    dens = outflowmodel.FindnH(nofluxsim)
    rs = outflowmodel.Findrstromgren(sim,dens=dens)
    cii = outflowmodel.cs
    def Findtion(radius):
        # Try spitzer
        tion = 4.0/7.0 * rs / cii * \
            ((radius/rs)**(7.0/4.0) - 1.0)
        #print "TION FOUND:", tion/outflowmodel.Myrins
        #print "USING DENS", dens
        #print radius, rs
        return tion
    def Findrspitzer(ts):
        ts = ts*outflowmodel.Myrins
        return rs * (1.0 + 7.0/4.0 * cii / rs * ts)**(4.0/7.0)
    # Set thresholds
    nthresh = 3e5
    if "_C2" in simname:
        nthresh = 1e6
    elif "_C" in simname:
        nthresh = 1e6
    if myr is None:
        snap = sim.Snapshots()[0]
        myr = snap.RawData().info["unit_time"].express(C.Myr)
    snap = sim.FindAtTime(tstart/myr)
    clumps = FindClumpsHamu(snap,nthresh=nthresh)
    tstart = 0.0
    tend = 5.0
    ntimes = 10000
    times = np.arange(tstart,tend,(tend - tstart)/ntimes)
    ntimes = len(times)
    mtots = np.zeros((ntimes))
    flux = 10.0**float(simname[1:3]) # "N48_M4_B02"

    S49 = flux/1e49
    phi1 = 1.0
    c5 = 0.545
    kmsinpcMyr = 1.02269032 # 1 km/s / (pc / Myr)
    #c0 = 0.545*1e5
    #vturb = 10*1e5
    #spd = np.sqrt(vturb**2 + c0**2)
    #c5 = spd/1e5
    # HACK
    #c5 = 1.0
    nclump = 1e5

    def dmdt(mfrac):
        mf0p4 = mfrac**0.4
        rfrac = 1.5 * betaR * ((0.4 * v1 / vR + 1) * (1 - mf0p4) + \
                                   mf0p4 * np.log(mf0p4)) + 1
        rfrac = rfrac**(5.0/3.0)
        return -2.5 * mfrac**0.6 / tev * rfrac ** (-0.4)

    for clump in clumps:
        m1 = clump.mass
        R1 = clump.distance
        rc = clump.radius
        nc = clump.peakdensity
        v1 = 0.0 # Try this, TODO: replace with actual velocity?
        print "rc/pc", rc#*outflowmodel.pcincm
        print "DISTANCE", R1
        n03 = clump.peakdensity/1e3
        b = clump.bfactor
        # Accretion phase
        tion = Findtion(R1*outflowmodel.pcincm)/outflowmodel.Myrins#+tstart
        mclump = m1 * (1.0 + times/tff)
        m1 *= 1.0 + tion/tff
        vR = 0.8 * c5*1e5
        if not mag:
            # Equation 4.10a
            tev = 4.48e5 * phi1 * c5 ** (-1.2) * S49 ** (-0.2) * \
                R1 ** 0.2 * m1 ** 0.4 / 1e6 # in Myr
            # Equation 4.19
            betaR = 5.22 * 0.8 * (vR / (c5*1e5))* c5**(-1.2) * \
                S49 ** (-0.2) * R1 ** (-0.6) * m1 ** (0.4)
        else:
            # Equation 4.10b
            tev = 1.30e5 * phi1 * (n03**(1.0/6.0)/b)**(6.0/7.0) * \
                S49 ** (-2.0/7.0) * R1 ** (4.0/7.0) * m1 ** (3.0/7.0) / 1e6
        # KH instability
        # diff of Spitzer solution
        if flux > 1e2:
            # Follow shell model
            #vwind = cii ** (4.0/7.0) * \
            #    (4.0/7.0 * rs / ((times-tion)*outflowmodel.Myrins)) ** (3.0/7.0)
            # Constant radius model
            vwind = 2.0/7.0 * R1*outflowmodel.pcincm / \
                ((times-tion)*outflowmodel.Myrins)
        else:
            vwind = 1e30
        print "vwind", vwind/1e5
        # Use initial density 1e3
        rspitz = Findrspitzer(times-tion)
        next = 1e2*(rspitz/rs)**(-3.0/2.0)
        contrast = nc/next
        print "contrast", contrast
        tKH = 2.0*rc * outflowmodel.pcincm * contrast**0.5 / vwind
        tKH /= outflowmodel.Myrins
        # Rayleigh Taylor model
        #aRT = 1.66e-24/0.76 *
        print "tKH, tev", tKH, tev
        # Do we actually have a flux?
        if flux > 1e2:
            # Linear model
            #tfact = (1.0 - (times - tion) / tev)
            #tfact[tfact < 0.0] = 0.0
            #mevap = m1 * tfact** (5.0/3.0)
            # Integration model
            nint = ntimes
            dt = times[1] - times[0] # TIMES SHOULD HAVE SAME DT OVER ARRAY
            mevap = np.zeros((ntimes))
            mevap += 1.0
            istart = np.where(times == times[times > tion].min())[0]
            for iint in range(istart,nint-1):
                mevap[iint+1] = mevap[iint]+dmdt(mevap[iint])*dt
            mevap[np.isnan(mevap)] = 0.0
            mevap *= m1
            # Exponential evaporation model
            #mevap = m1 * np.exp(-(times-tion)/(tev))
            # KH model
            #mevap *= np.exp(-(times-tion)/(tKH))
            mclump[times > tion] = mevap[times > tion]
        mclump[mclump < 0.0] = 0.0
        mtots += mclump
    return times, mtots

def BondiAccretion(snap,times):
    '''
    Calculate a Bondi accretion rate
    '''
    clumps = FindClumpsHamu(snap,nthresh=1e6)
    #times = np.arange(0,5,0.1)
    mtots = np.zeros((len(times)))
    #flux = 10.0**float(simname[1:3]) # "N48_M4_B02"
    myrins = 3.15569e13
    msolaring = 1.9891e33
    #S49 = flux/1e49
    phi1 = 1.0
    c0 = 0.545*1e5
    vturb = 10*1e5
    spd = np.sqrt(vturb**2 + c0**2)
    G = 6.674e-8
    nH = 3e4 # Eh, dunno
    dens = 1.66e-24 / 0.76 * nH
    coeff = 4.0 * np.pi * G**2 / spd**3 * dens
    for clump in clumps:
        m1 = clump.mass
        #print times * myrins * coeff * m1 * msolaring
        mclump = m1 / (1.0 - times * myrins * coeff * m1 * msolaring)
        mclump[mclump < 0.0] = 0.0
        mtots += mclump
    #print mtots
    print times, mtots
    return times, mtots
    

def PlotEvaporation(simnames):
    plt.clf()
    cols = ["k","b","m","r"]
    icol = 0
    for simname in simnames:
        col = cols[icol]
        sim = Hamu.Simulation(simname)
        times, masses = EvaporationRate(simname,time=1.25,mag=False)
        plt.plot(times, masses,col,
                 label=sim.Label()+" Thermal")
        times, masses = EvaporationRate(simname,time=1.25,mag=True)
        plt.plot(times, masses,col+"--",
                 label=sim.Label()+" Magnetic")
        icol += 1
    plt.legend(fontsize="x-small",loc="upper right")
    plt.yscale("log")
    plt.ylabel("Total Mass / M$_{\odot}$")
    plt.xlabel("Time / Myr")
    plt.savefig("../plots/clumps/photons/evaporation.pdf")
    print simnames


def PlotTEvap(simnames,tfact=1.0):
    plt.clf()
    cols = ["b","m","r","k"]
    icol = 0
    mdsyms = []
    labels = []
    for simname in simnames:
        tff = 1.25
        if "_C2" in simname:
            tff = 1.25*0.75**3
        elif "_C" in simname:
            tff = 1.25*0.5**3
        col = cols[icol]
        sim = Hamu.Simulation(simname)
        clumps, tevs = EvaporationRate(simname,time=tff*tfact,mag=False,tevout=True)
        masses = [clump.mass for clump in clumps]
        print len(masses), len(tevs)
        #plt.plot(masses,tevs,col,
        #         label=sim.Label())

        mdsym = plt.scatter(masses, tevs,
                    #s=areas,
                    alpha=0.7,
                    label=sim.Label(),
                    color=col)
        labels.append(sim.Label())
        mdsyms.append(mdsym)
        icol += 1
    if len(simnames) > 1:
        plt.legend(mdsyms,labels,scatterpoints=1,
                   fontsize="x-small",loc="lower right")
        icol += 1
    #plt.legend(fontsize="xx-small",loc="upper right")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel("Evaporation Time / Myr")
    plt.xlabel("Clump Mass / Msolar")
    tfacttxt = ""
    if tfact != 1.0:
        tfacttxt = "_"+str(tfact)+"t_ff"
        tfacttxt = tfacttxt.replace(".","p")
    plt.savefig("../plots/clumps/photons/tevap"+tfacttxt+".pdf")
    print simnames

def RunForSim(simname):
    sim = Hamu.Simulation(simname)
    snap = sim.Snapshots()[-1]
    clumps = FindClumps(snap)
    ClumpProperties(clumps,sim.Name())

def CompareForTimes(simname):
    folder = simname+""
    sim = Hamu.Simulation(simname)
    for snap in sim.Snapshots():
        suffix = r"_"+str(snap.OutputNumber()).zfill(5)
        myr = snap.RawData().info["unit_time"].express(C.Myr)
        time = snap.Time()*myr
        CompareMassDistance([simname],time=time,
                            folder=folder,suffix=suffix)

def FindDist(snap):
    # Find mean distance from clumps to centre of mass
    clumps = FindClumpsHamu(snap.hamusnap)
    masses = np.array([clump.mass for clump in clumps])
    distances = np.array([clump.distance for clump in clumps])
    dist = np.sum(distances*masses) / np.sum(masses)
    if len(distances) > 0:
        print "MINMAX DIST", distances.min(), distances.max()
        print "MASS WEIGHTED MEAN DIST", dist
    return dist
FindDistHamu = Hamu.Algorithm(FindDist)


def FindCoMDist(snap):
    # Find distance between source and centre of mass
    clumps = FindClumpsHamu(snap.hamusnap)
    masses = np.array([clump.mass for clump in clumps])
    x = np.array([clump.centre[0][0] for clump in clumps])
    y = np.array([clump.centre[0][1] for clump in clumps])
    z = np.array([clump.centre[0][2] for clump in clumps])
    dist = np.sqrt(np.sum(x*masses)**2 + \
                       np.sum(y*masses)**2 + \
                       np.sum(z*masses)**2) \
                       / np.sum(masses)
    return dist
FindCoMDistHamu = Hamu.Algorithm(FindCoMDist)

def FindDistToCoM(snap):
    # Find mean clump distance to centre of mass
    clumps = FindClumpsHamu(snap.hamusnap)
    masses = np.array([clump.mass for clump in clumps])
    if len(masses) == 0:
        return 0.0
    mtot = np.sum(masses)
    # Find com
    x = np.array([clump.centre[0][0] for clump in clumps])
    y = np.array([clump.centre[0][1] for clump in clumps])
    z = np.array([clump.centre[0][2] for clump in clumps])
    com = np.array([np.sum(x*masses), np.sum(y*masses),np.sum(z*masses)])
    com /= mtot
    # Find distance to com
    x -= com[0]
    y -= com[1]
    z -= com[2]
    dist = np.sqrt(x**2 + y**2 + z**2)
    dist = np.sum(dist*masses) / mtot
    return dist
FindDistToCoMHamu = Hamu.Algorithm(FindDistToCoM)

def FindMinDist(snap):
    # Find minimum distance between clump and source
    clumps = FindClumpsHamu(snap.hamusnap)
    # NOTE: include radius
    x = np.array([clump.centre[0][0] - clump.radius for clump in clumps])
    if len(x) == 0:
        return 0.0
    y = np.array([clump.centre[0][1] - clump.radius for clump in clumps])
    z = np.array([clump.centre[0][2] - clump.radius for clump in clumps])
    mindist = np.sqrt(np.min(x**2 + y**2 + z**2))
    return mindist
FindMinDistHamu = Hamu.Algorithm(FindMinDist)

def ClumpDistanceOverTime(simnames,cols,lines,name="",com=False,mindist=False):
    '''
    Calculate the mean mass-weighted clump distance over time
    '''
    plt.clf()
    if mindist:
        com = False
    for col, line, simname in zip(cols,lines,simnames):
        sim = Hamu.Simulation(simname)
        print "Running ClumpDistanceOverTime for", simname
        times = []
        dists = []
        # Gather data
        for snap in sim.Snapshots():
            suffix = r"_"+str(snap.OutputNumber()).zfill(5)
            myr = snap.RawData().info["unit_time"].express(C.Myr)
            time = snap.Time()*myr
            if not mindist:
                if not com:
                    dist = FindDistToCoMHamu(snap)
                else:
                    dist = FindCoMDistHamu(snap)
            else:
                dist = FindMinDistHamu(snap)
            dists.append(dist)
            times.append(time)
        # Plot
        plt.plot(times,dists,color=col,linestyle=line,label=sim.Label())
    # Plot radius of massive clump as example
    if not mindist:
        # (Mindist already subtracts the clump radius)
        plt.plot([0,5],[0.25,0.25],"k:")
    plt.legend(fontsize="x-small",loc="upper left",ncol=2,
               frameon=False)
    plt.yscale("log")
    plt.xlabel("Time / Myr")
    plt.ylabel("Mass-weighted mean distance to centre of mass / pc")
    plt.xlim([0,5])
    try:
        os.mkdir("../plots/clumps/")
    except:
        pass
    plt.ylim([0.1,10])
    if len(name) > 0:
        name = "disttocom_"+name
    else:
        name = "disttocom"
    if com:
        plt.ylabel("Distance of Centre of Mass to Source / pc")
        plt.ylim([0.1,4])
        if len(name) > 0:
            name = "centreofmass_"+name
        else:
            name = "centreofmass"
    if mindist:
        plt.ylabel("Minimum clump distance to source / pc")
        plt.ylim([0.1,10])
        if len(name) > 0:
            name = "mindist_"+name
        else:
            name = "mindist"
        
    if len(name) > 0:
        name = "_"+name
    plt.savefig("../plots/clumps/distanceovertime"+name+".pdf")

if __name__=="__main__":
    #RunForSim("N49_M4_B02")
    simnames = ["N"+num+"_M4_B02" for num in ["00","47","48","49"]]
    #CompareMassDistance(simnames)
    #CompareMassDistance(["N00_M4_B02"],time=1.45)
    #CompareForTimes("N48_M4_B02_C")
    #PlotEvaporation(simnames)
    simnames = ["N"+num+"_M4_B02" for num in ["47","48","49"]]
    simnames.append("N48_M4_B02_C")
    #PlotTEvap(simnames,tfact=2.0)
    simnames = ["N00_M4_B02","N47_M4_B02","N48_M4_B02","N49_M4_B02",
                "N00_M4_B02_C", "N48_M4_B02_C",
                "N00_M4_B02_C2","N48_M4_B02_C2"]
    cols = [linestyles.col(simname) for simname in simnames]
    lines = [linestyles.line(simname) for simname in simnames]
    #lines = ["k--","r","k","m",
    #         "c--","c",
    #         "b--","b"]
    #ClumpDistanceOverTime(simnames,cols,lines,"all",com=False)
    #ClumpDistanceOverTime(simnames,cols,lines,"all",com=True)
    ClumpDistanceOverTime(simnames,cols,lines,"all",mindist=True)

