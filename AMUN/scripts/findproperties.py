'''
Plot properties in simulations
Sam Geen, November 2015
'''

from startup import *

from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import matplotlib.patheffects as pe

import triaxfinder

def massinsnap(snap,nHthresh=0.0):
    print "Finding mass in snap", snap.iout
    amr = snap.amr_source(["rho","P"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    vols = (cells.get_sizes())**3.0
    dens = cells["rho"]
    mass = dens*vols*snap.info["unit_mass"].express(C.Msun)
    if nHthresh > 0.0:
        mass = mass[dens > nHthresh]
    return np.sum(mass)

def Pnontherminsnap(snap,nHthresh=0.0):
    print "Finding Pnontherm in snap", snap.iout
    amr = snap.amr_source(["Pnontherm"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    vols = (cells.get_sizes())**3.0
    return np.max(cells["Pnontherm"])

def meandensinsnap(snap,nHthresh=0.0):
    print "Finding mass in snap", snap.iout
    amr = snap.amr_source(["rho","P"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    vols = (cells.get_sizes()*snap.info["unit_length"].express(C.cm))**3.0
    dens = cells["rho"]*snap.info["unit_density"].express(C.H_cc)
    mass = dens*vols
    if nHthresh > 0.0:
        mass = mass[dens > nHthresh]
        vols = vols[dens > nHthresh]
    return np.sum(mass)/np.sum(vols)

def etherminsnap(snap,wind=False):
    print "Finding max T in snap", snap.iout
    amr = snap.amr_source(["rho","P"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    vols = (cells.get_sizes())**3.0
    mass = cells["rho"]*vols*\
        snap.info["unit_mass"].express(C.g)
    temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    etherm = 1.0/(gamma - 1.0) * (mass / mH) * kB * temp
    if wind:
        mask = np.where(temp > 1e5)
        etherm = etherm[mask]
    return np.sum(etherm)

def ekininsnap(snap,wind=False):
    print "Finding kinetic energy in snap", snap.iout
    amr = snap.amr_source(["rho","vel"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    umass = snap.info["unit_mass"].express(C.g) 
    uvel = snap.info["unit_velocity"].express(C.cm/C.s)
    ue = umass*uvel**2
    rhos = cells["rho"]
    vels = cells["vel"]
    vols = (cells.get_sizes())**3.0
    spds = np.sqrt(np.sum(vels**2.0,1))
    ekin = 0.5*rhos*vols*spds**2
    if wind:
        # Over 100 km/s
        mask = np.where(spds*uvel/1e5 > 100.0)
        ekin = ekin[mask]
    ekin =  np.sum(ekin)*ue
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    print "KINETIC ENERGY FOUND @ ",time,"Myr:", ekin
    return ekin

def windenergyinsnap(snap):
    return energyinsnap(snap,wind=True)

def energyinsnap(snap,wind=False):
    etherm = etherminsnap(snap,wind)
    ekin = ekininsnap(snap,wind)
    etot = ekin+etherm
    return (etherm,ekin,etot)

def maxTinsnap(snap):
    print "Finding max T in snap", snap.iout
    amr = snap.amr_source(["rho","P"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    maxT = temp.max()
    return maxT

def momentuminsnap(snap,nHlow=0.0,nHhigh=0.0):
    print "Finding momentum in snap", snap.iout
    amr = snap.amr_source(["rho","vel"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    uvel = snap.info["unit_velocity"].express(C.cm/C.s)
    rhos = cells["rho"]
    vels = cells["vel"]
    vols = (cells.get_sizes())**3.0
    spds = np.sqrt(np.sum(vels**2.0,1))
    pos = cells.points-0.5
    rads = pos+0.0
    dist = np.sqrt(np.sum(pos**2,1))
    for i in range(0,3):
        rads[:,i] /= dist
    rvel = np.sum(rads*vels,1)
    rvelp = rvel+0.0
    rvelp[rvelp < 0.0] = 0.0
    moms = rhos*rvelp*vols
    # Apply threshold?
    if nHlow > 0.0:
        moms = moms[rhos > nHlow]
        rhos = rhos[rhos > nHlow]
    if nHhigh > 0.0:
        moms = moms[rhos < nHhigh]
        rhos = rhos[rhos < nHhigh]
    mom =  np.sum(moms)
    umass = snap.info["unit_mass"].express(C.g)
    mom *= umass*uvel
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    print "MOMENTUM FOUND @ ",time,"Myr:", mom
    return mom 

def totalmomentuminsnap(snap,nHlow=0.0,nHhigh=0.0):
    print "Finding *total* momentum in snap", snap.iout
    amr = snap.amr_source(["rho","vel"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    uvel = snap.info["unit_velocity"].express(C.cm/C.s)
    rhos = cells["rho"]
    vels = cells["vel"]
    vols = (cells.get_sizes())**3.0
    spds = np.sqrt(np.sum(vels**2.0,1))
    moms = rhos*spds*vols
    # Apply threshold?
    if nHlow > 0.0:
        moms = moms[rhos > nHlow]
        rhos = rhos[rhos > nHlow]
    if nHhigh > 0.0:
        moms = moms[rhos < nHhigh]
        rhos = rhos[rhos < nHhigh]
    mom =  np.sum(moms)
    umass = snap.info["unit_mass"].express(C.g)
    mom *= umass*uvel
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    print "MOMENTUM FOUND @ ",time,"Myr:", mom
    return mom

def windradiusinsnap(snap):
    return radiusinsnap(snap,wind=True)

def radiusinsnap(snap,wind=False):
    print "Finding radius of HII region in snap", snap.iout
    boxlen = snap.info["boxlen"]
    amr = snap.amr_source(["rho","P","xHII"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    ion = cells["xHII"]
    vols = (cells.get_sizes())**3.0
    pos = cells.points+0.0
    rhos = cells["rho"]
    mcode = vols*rhos
    msum = np.sum(mcode)
    thresh = 0.1 # HACK!
    if not wind:
        mask = np.where(ion > thresh)
        weight = mcode[mask]*ion[mask]
    else:
        mask = np.where(temp > 1e5)
        weight = mcode[mask]
    x = pos[mask,0]*boxlen
    y = pos[mask,1]*boxlen
    z = pos[mask,2]*boxlen
    ncells = weight.shape[0]
    if ncells == 0:
        return 0.0
    com = np.array([np.sum(x*weight),np.sum(y*weight),np.sum(z*weight)]) / np.sum(weight)
    x -= com[0]
    y -= com[1]
    z -= com[2]
    dist = np.sqrt(x**2+y**2+z**2)
    radius = np.sum(dist*weight) / np.sum(weight)
    #import pdb; pdb.set_trace()
    #print x.shape, y.shape, z.shape
    #a,b,c = triaxfinder.findaxes(weight,x,y,z,ncells)
    #radius = np.sqrt(a**2+b**2+c**2)
    #print "RADIUS FOUND", radius
    return radius

def radiusinsnap3(snap):
    '''
    2nd implementation: measure volume of ionised gas and sphericise
    '''
    print "Finding radius of HII region in snap", snap.iout
    boxlen = snap.info["boxlen"]
    amr = snap.amr_source(["rho","xHII"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    #temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    ion = cells["xHII"]
    vols = (cells.get_sizes())**3.0
    pos = cells.points+0.0
    rhos = cells["rho"]
    mcode = vols*rhos
    msum = np.sum(mcode)
    # Get volume of ionised gas
    # Assume gas is either fully ionised or neutral on a sub-grid scale
    # Also include a basic threshold
    thresh = 0.1 # HACK!
    mask = np.where(ion > thresh)
    ionvol = np.sum(vols[mask]*ion[mask])
    ionrad = ionvol**(1.0/3.0) * (3.0 / 4.0 / np.pi)
    return ionrad*boxlen

def maxBfieldinsnap(snap):
    amr = snap.amr_source(["B-left","B-right"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    b = 0.5*(cells["B-left"]+cells["B-right"])  
    ud = snap.info["unit_density"].express(C.g_cc)
    ul = snap.info["unit_length"].express(C.cm)
    ut = snap.info["unit_time"].express(C.s)
    unit = np.sqrt(4.0*np.pi*ud*(ul/ut)**2)*1e6 # microGauss
    Bmag = np.sqrt((b**2).sum(axis=1))*unit
    return np.max(Bmag)
    
def run(func,simnames,plotname):
    name = func.__name__
    plt.clf()
    hamufunc = Hamu.Algorithm(func)
    # Get and plot lines
    for simname in simnames:
        print "Running function", name, "for simulation", simname
        sim = hamusims[simname]
        times, Ts = timefuncs.timefunc(sim,hamufunc)
        plt.plot(times,Ts,
                 color=linestyles.Colour(simname),label=linestyles.Label(simname),
                 path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()]) 
    # Do plot stuff
    plt.xlabel("Time / Myr") 
    plt.legend(frameon=False,fontsize="x-small") 
    # Quantity-specific stuff
    if name == "totalmomentuminsnap":
        plt.ylabel("Momentum / g cm / s")
    else:
        plt.ylabel(name)
        
    #if name == "energyinsnap":
    #    xlims = np.array(plt.gca().get_xlim())
    #    plt.plot(xlims,xlims*0.0+1e51,"k:")
    #    plt.ylim([1e49,2e51])
    #if tlog:
    #    plt.xscale("log")
    #    plt.xlim([1e-6,1])
    plt.yscale("log")
    #if "maxT" in name:
    #    plt.ylabel("Temperature / K")
    #if "momentum" in name:
    #    if nsnsub:
    #        plt.ylabel("Added Momentum / g cm/s")
    #    else:
    #        plt.ylabel("Total Momentum / g cm/s")
    #    if Hamu.CurrentWorkspace() == "HIISN" and not nsnsub:
    #        plt.ylim([4e43,2e45])
    #    if Hamu.CurrentWorkspace() == "HIISN4":
    #        plt.ylim([2e42,2e43])
    # Check for user-specified ylim
    #if ylim is not None:
    #    plt.ylim(ylim)
    if "radius" in name:
        boxlen = 0.121622418993404E+03
        xlim = plt.gca().get_xlim()
        plt.plot(xlim,[0.5*boxlen,0.5*boxlen],"k--")
        plt.ylabel("Mean HII region radius / pc")
        plt.yscale("linear")
    #plt.xlabel("Time After SN / Myr")
    #ws = Hamu.CurrentWorkspace()
    #suffix = filesuffix
    #if suffix == "":
    #    if nsnsub:
    #        suffix += "_nsnsub"
    #    if tlog:
    #        suffix += "_tlog"
    plt.savefig("../plots/"+name+"_"+plotname+".pdf")

if __name__=="__main__":
    for func in [totalmomentuminsnap,radiusinsnap][::-1]:
        run(func,imfsims,"imf")
        run(func,icsims,"ic")
