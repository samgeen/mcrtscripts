'''
Plot properties in simulations
Sam Geen, November 2015
'''

from startup import *

from pymses.filters import CellsToPoints
from pymses.utils import constants as C
from pymses.analysis import sample_points, bin_cylindrical
from pymses.utils.regions import Cylinder
from pymses.filters import CellsToPoints,RegionFilter

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
    # Physical conversions
    X = 0.76
    mH = 1.6735326381e-24
    kB = 1.38062e-16
    G = 6.674e-8
    gamma = 1.4 # RAMSES hard-coded
    pcincm = 3.086e18
    Msuning = 1.9891e33
    mHing = 1.66e-24
    Myrins = 3.1556926e13
    print "Finding max T in snap", snap.iout
    amr = snap.amr_source(["rho","P","vel","xHII","xHeII","xHeIII"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    vols = (cells.get_sizes())**3.0
    mass = cells["rho"]*vols*\
        snap.info["unit_mass"].express(C.g)
    mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                              0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))  
    temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)*mufunc(cells)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    etherm = 1.0/(gamma - 1.0) * (mass / mH) * kB * temp
    if wind:
        vels = cells["vel"]
        uvel = snap.info["unit_velocity"].express(C.cm/C.s)
        spds = np.sqrt(np.sum(vels**2.0,1))
        mask = np.logical_or(spds*uvel/1e5 > 100.0,temp > 1e5)
        etherm = etherm[mask]
    return np.sum(etherm)

def ekininsnap(snap,wind=False):
    print "Finding kinetic energy in snap", snap.iout
    amr = snap.amr_source(["rho","P","vel","xHII","xHeII","xHeIII"])
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
        mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                                  0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))  
        temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)*mufunc(cells)
        # Over 100 km/s or T > 1e5 K
        mask = np.logical_or(spds*uvel/1e5 > 100.0,temp > 1e5)
        ekin = ekin[mask]
    ekin =  np.sum(ekin)*ue
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    print "KINETIC ENERGY FOUND @ ",time,"Myr:", ekin
    return ekin

def Lcoolinsnap(snap,wind=False):
    print "Finding cooling luminosity in snap", snap.iout
    import rtcooling
    amr = snap.amr_source(["rho","P","vel","xHII","xHeII","xHeIII"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    umass = snap.info["unit_mass"].express(C.g) 
    uvel = snap.info["unit_velocity"].express(C.cm/C.s)
    ue = umass*uvel**2
    rhos = cells["rho"]
    vels = cells["vel"]
    vols = (cells.get_sizes()*snap.info["unit_length"].express(C.cm))**3.0
    spds = np.sqrt(np.sum(vels**2.0,1))
    mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                                0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))  
    temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)*mufunc(cells)
    dedt = rtcooling.Finddedt(T2,nH,xion,Zsolar,Np=None,Fp=None,p_gas=None,a_exp=np.array([1.0])):
    Lcool = dedt * vols
    if wind:
        # Over 100 km/s or T > 1e5 K
        windmask = np.logical_or(spds*uvel/1e5 > 100.0,temp > 1e5)
        Lcool = Lcool[mask]
    Lcool =  np.sum(Lcool)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    print "COOLING LUMINOSITY IN erg/s FOUND @ ",time,"Myr:", Lcool
    return Lcool

# Version 1: ekin only measures fast winds, etherm only hot winds
# Version 2: Modified functions above to count both fast and hot winds for both ekin and etherm
def windenergyinsnap2(snap):
    return energyinsnap(snap,wind=True)

def windLcoolinsnap(snap):
    return Lcoolinsnap(snap,wind=True)

def energyinsnap(snap,wind=False):
    etherm = etherminsnap(snap,wind)
    ekin = ekininsnap(snap,wind)
    etot = ekin+etherm
    return (etherm,ekin,etot)

def energytotinsnap(snap,wind=False):
    etherm = etherminsnap(snap,wind)
    ekin = ekininsnap(snap,wind)
    etot = ekin+etherm
    return tot

def maxTinsnap(snap):
    print "Finding max T in snap", snap.iout
    amr = snap.amr_source(["rho","P"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    maxT = temp.max()
    return maxT

def totalmomentuminsnapS(snap,nHlow=0.0,nHhigh=0.0):
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

    # Add momentum of sinks
    currsinks = sinks.FindSinks(snap.hamusnap)
    sinkvx = currsinks.vx
    sinkvy = currsinks.vy
    sinkvz = currsinks.vz
    sinkm  = currsinks.mass
    spdsS = np.sqrt(sinkvx**2+sinkvy**2+sinkvz**2)
    momsS = sinkm*spdsS

    # Apply threshold?
    if nHlow > 0.0:
        moms = moms[rhos > nHlow]
        rhos = rhos[rhos > nHlow]
    if nHhigh > 0.0:
        moms = moms[rhos < nHhigh]
        rhos = rhos[rhos < nHhigh]
    mom   = np.sum(moms) 
    momS  = np.sum(momsS)
    umass = snap.info["unit_mass"].express(C.g)
    mom  *= umass*uvel
    momS *= umass*uvel
    momT  =  mom + momS
    time  = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    print "MOMENTUM FOUND @ ",time,"Myr:", momT
    return momT

def momentuminsnap(snap,centre=(0.5,0.5,0.5),nHlow=0.0,nHhigh=0.0):
    print "Finding momentum in snap", snap.iout
    amr = snap.amr_source(["rho","vel"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    uvel = snap.info["unit_velocity"].express(C.cm/C.s)
    rhos = cells["rho"]
    vels = cells["vel"]
    vols = (cells.get_sizes())**3.0
    spds = np.sqrt(np.sum(vels**2.0,1))
    pos = cells.points+0.0
    pos[:,0] -= centre[0]
    pos[:,1] -= centre[1]
    pos[:,2] -= centre[2]
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

def radialmomentumatstarinsnap(snap):
    stellar = stellars.FindStellar(snap.hamusnap)
    if len(stellar.mass) == 0:
        return 0.0
    imax = np.where(stellar.mass == np.max(stellar.mass))[0][0]
    sinkid = stellar.sinkid[imax]-1
    sink = sinks.FindSinks(snap.hamusnap)
    boxlen = snap.info["boxlen"]
    starpos = np.array([sink.x[sinkid],sink.y[sinkid],sink.z[sinkid]])/boxlen
    return momentuminsnap(snap,starpos)

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
    return radiusinsnap3(snap,wind=True)

def densityprofileinsnap(snap):
    print "Finding density profile in snap", snap.iout

    boxlen = snap.info["boxlen"] * 3.086e+16   # m

    nbins = 200
    center = [ 0.5, 0.5, 0.5 ]        # in box units
    radius = 1.                       # in box units
    thickn = 0.2 		      # in box units 
    normal = [ 0, 0, 1 ]              # Norm = 1

    factorL   = snap.info["unit_length"].express(C.pc)          # pc  
    factorD   = snap.info["unit_density"].express(C.g_cc)       # g.cm^-3
    factorD2  = snap.info["unit_density"].express(C.H_cc)       # H.cm^-3
    factorM   = snap.info["unit_mass"].express(C.g)             # g
    factorM   = factorM / (1.9889 * 1e33)                       # Msun

    cyl       = Cylinder(center, normal, radius, thickn)

    # AMR density field point sampling
    source = snap.amr_source(["rho"])
    filt_source = RegionFilter(cyl, source)
    filt_cell_source = CellsToPoints(filt_source)
    dset = filt_cell_source.flatten()

    mass_weight_func  = lambda dset: dset["rho"] *dset["size"]**3 * factorM
    rho_weight_func   = lambda dset: dset["rho"]                  * factorD	

    r_bins = np.linspace(0, radius, nbins)            # in box units 
    # Geometrical midpoint of the bins
    bins_centers = (r_bins[1:]+r_bins[:-1])/2.    # in box units
    bin_width    = dtR                            # in box units

    # Profile computation
    rho_profile   = bin_cylindrical(dset, center, normal, rho_weight_func, r_bins, divide_by_counts=True)  # for density profile

    return rho_profile
    

def surfacedensityprofileinsnap(snap):
    print "Finding surface density profile in snap", snap.iout

    boxlen = snap.info["boxlen"]      # pc

    nbins = 200
    center = [ 0.5, 0.5, 0.5 ]        # in box units
    radius = 1.                       # in box units
    thickn = 0.2 		      # in box units 
    normal = [ 0, 0, 1 ]              # Norm = 1

    factorL   = snap.info["unit_length"].express(C.pc)          # pc  
    factorD   = snap.info["unit_density"].express(C.g_cc)       # g.cm^-3
    factorD2  = snap.info["unit_density"].express(C.H_cc)       # H.cm^-3
    factorM   = snap.info["unit_mass"].express(C.g)             # g
    factorM   = factorM / (1.9889 * 1e33)                       # Msun

    cyl       = Cylinder(center, normal, radius, thickn)

    # AMR density field point sampling
    source = snap.amr_source(["rho"])
    filt_source = RegionFilter(cyl, source)
    filt_cell_source = CellsToPoints(filt_source)
    dset = filt_cell_source.flatten()

    mass_weight_func  = lambda dset: dset["rho"] *dset["size"]**3 * factorM
    rho_weight_func   = lambda dset: dset["rho"]                  * factorD	

    r_bins = np.linspace(0, radius, nbins)            # in box units 
    # Geometrical midpoint of the bins
    bins_centers = (r_bins[1:]+r_bins[:-1])/2.    # in box units
    bin_width    = dtR                            # in box units

    # Profile computation
    rho_profile   = bin_cylindrical(dset, center, normal, rho_weight_func, r_bins, divide_by_counts=False)
    mass_profile  = bin_cylindrical(dset, center, normal, mass_weight_func, r_bins, divide_by_counts=False)

    surf_dens_prof = []
    for j, rr in enumerate(r_bins):
        if (j > 0):
            Apc  = np.pi * ( (rr * boxlen)**2 - (r_bins[j-1] * boxlen)**2 )  # pc^2
            surf_dens_prof.append(mass_profile[j-1]/Apc)                     # Msun/pc^2
 
    exit() 
    
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

def radiusinsnap3(snap,wind=False):
    '''
    2nd implementation: measure volume of ionised gas and sphericise
    '''
    print "Finding radius of HII region in snap", snap.iout
    boxlen = snap.info["boxlen"]
    amr = snap.amr_source(["rho","xHII","xHeII","xHeIII","P","vel"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    #temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    ion = cells["xHII"]
    vols = (cells.get_sizes())**3.0
    pos = cells.points+0.0
    rhos = cells["rho"]
    mcode = vols*rhos
    msum = np.sum(mcode)
    # Speeds for winds
    if wind:
        uvel = snap.info["unit_velocity"].express(C.cm/C.s)
        vels = cells["vel"]
        spds = np.sqrt(np.sum(vels**2.0,1))
        mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                                  0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))  
        temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)*mufunc(cells)
        # Over 100 km/s or T > 1e5 K
        mask = np.logical_or(spds*uvel/1e5 > 100.0,temp > 1e5)
    else:
        thresh = 0.1 # fiducial "non-ionised" threshold
        mask = np.where(ion > thresh)
    # Get volume of ionised gas
    # Assume gas is either fully ionised or neutral on a sub-grid scale
    # Also include a basic threshold
    ionvol = np.sum(vols[mask]*ion[mask])
    ionrad = ionvol**(1.0/3.0) * (3.0 / 4.0 / np.pi)
    return ionrad*boxlen

def photodensinsnap(snap):
    '''
    Average density of photoionised gas (not including wind-shocked gas)
    '''
    print "Finding photoionised gas density of HII region in snap", snap.iout
    boxlen = snap.info["boxlen"]
    amr = snap.amr_source(["rho","xHII","xHeII","xHeIII","P","vel"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    #temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    ion = cells["xHII"]
    vols = (cells.get_sizes())**3.0
    pos = cells.points+0.0
    rhos = cells["rho"]
    mcode = vols*rhos
    msum = np.sum(mcode)
    # Speeds for winds
    uvel = snap.info["unit_velocity"].express(C.cm/C.s)
    vels = cells["vel"]
    spds = np.sqrt(np.sum(vels**2.0,1))
    mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                              0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))  
    temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)*mufunc(cells)
    # Over 100 km/s or T > 1e5 K
    notwindmask = np.logical_and(spds*uvel/1e5 < 100.0,temp < 1e5)
    thresh = 0.1 # fiducial "non-ionised" threshold
    mask = np.logical_and(ion > thresh,notwindmask)
    # Get volume-weighted average density of ionised gas
    ionvols = vols[mask]*ion[mask]
    ionmass = rhos[mask] * ionvols
    iondens = np.sum(ionmass) / np.sum(ionvols)
    return iondens

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
