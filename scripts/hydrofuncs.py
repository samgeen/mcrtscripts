'''
Functions defining hydro variables in RAMSES for use by pymses
Sam Geen, November 2014
'''

import numpy as np
from pymses.utils import constants as C
import rtcooling

def hydro_label(hydro):
    '''
    Returns a label for the given hydro variable in some default units
    '''
    if hydro == "rho":
        return "Density / atoms/cm$^{3}$"
    if hydro == "NH":
        return "$N_{\mathrm{H}}$ / cm$^{-2}$"
    if hydro == "P":
        return "Pressure / ergs/cm$^3$"
    if hydro == "T":
        return "Temperature / K"
    if "xH" in hydro:
        return "Ionisation Fraction "+hydro
    if hydro == "gpe":
        return "Gravitational Potential Energy"
    if hydro == "Bmag":
        return "|B| / $\mu G$"
    if hydro == "Bx":
        return "B.$x$ / $\mu G$"
    if hydro == "By":
        return "B.$y$ / $\mu G$"
    if hydro == "Bz":
        return "B.$z$ / $\mu G$"
    if hydro == "vrad":
        return "V$_{radial}$ / km/s"
    if hydro == "spd":
        return "Gas Speed / km/s"
    if hydro == "Pram":
        return "Ram Pressure / ergs/cm$^3$"
    if hydro == "vx":
        return "Velocity X / km/s"
    if hydro == "vy":
        return "Velocity Y / km/s"
    if hydro == "vz":
        return "Velocity Z / km/s"
    if hydro == "Lcool":
        return "L$_{cool}$ / erg/s"
    
def scale_by_units(ro, hydro):
    '''
    Returns a lambda function that scales a dataset by some default units
    '''
    def bunit(dset):
        ud = ro.info["unit_density"].express(C.g_cc)
        ul = ro.info["unit_length"].express(C.cm)
        ut = ro.info["unit_time"].express(C.s)
        unit = np.sqrt(4.0*np.pi*ud*(ul/ut)**2)*1e6 # microGauss
        return unit
    def bfunc(dset,axis):
        b = 0.5*(dset["B-left"]+dset["B-right"])
        return b[:,axis]*bunit(dset)
    if hydro == "rho":
        unit = ro.info["unit_density"].express(C.H_cc)
        return lambda dset: dset["rho"]*unit
    if hydro == "NH":
        unit = ro.info["unit_density"].express(C.g_cc) * \
               ro.info["unit_length"].express(C.cm)
        return scop(lambda dset: dset["rho"]*unit)
    if hydro == "P":
        unit = ro.info["unit_pressure"].express(C.barye)
        return lambda dset: dset["P"]*unit
    if hydro == "T":
        mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                     0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
        unit = ro.info["unit_temperature"].express(C.K)
        return lambda dset: dset["P"]/dset["rho"]*unit*mufunc(dset)
    if "xH" in hydro:
        unit = 1.0
        return lambda dset: dset[hydro]*unit
    if hydro == "Bmag":
        def bmagfunc(dset):
            #b = 0.5*(dset["B-left"]+dset["B-right"])
            #up = ro.info["unit_pressure"].express(C.barye)
            #unit = np.sqrt(up/(4.0*np.pi))*1e-6 # microGauss
            # Magnitude of the 3-vector for each cell
            b = 0.5*(dset["B-left"]+dset["B-right"])
            return np.sqrt((b**2).sum(axis=1))*bunit(dset)
        return bmagfunc
    if hydro == "Bx":
        return lambda dset: bfunc(dset,0)
    if hydro == "By":
        return lambda dset: bfunc(dset,1)
    if hydro == "Bz":
        return lambda dset: bfunc(dset,2)
    if hydro == "vrad":
        unit = ro.info["unit_velocity"].express(C.km/C.s)
        def findvrad(dset):
            pos = dset.points-0.5
            rad = pos # Be careful here! Reference, not value copy
            dist = np.sqrt(np.sum(pos**2,1))
            for i in range(0,3):
                rad[:,i] /= dist
            return np.sum(rad*dset["vel"],1)*unit
        return findvrad
    if hydro == "vx":
        unit = ro.info["unit_velocity"].express(C.km/C.s)
        return lambda dset: dset["vel"][:,0]*unit
    if hydro == "vy":
        unit = ro.info["unit_velocity"].express(C.km/C.s)
        return lambda dset: dset["vel"][:,1]*unit
    if hydro == "vz":
        unit = ro.info["unit_velocity"].express(C.km/C.s)
        return lambda dset: dset["vel"][:,2]*unit
    if hydro == "spd":
        unit = ro.info["unit_velocity"].express(C.km/C.s)
        return lambda dset: np.sqrt(np.sum(dset["vel"]**2,1))*unit
    if hydro == "Pram":
        unit = ro.info["unit_pressure"].express(C.barye)
        return lambda dset: dset["rho"]*np.sum(dset["vel"]**2,1)*unit
    if hydro == "Lcool":
        return rtcooling.dedtOnCells
    # None of those? Return unitless
    return lambda dset: dset[hydro]

def amr_source(ro, hydro,extra=[]):
    '''
    Load the necessary data from the hydro files
    '''
    if type(hydro) == type("duck"):
        hydros = [hydro] # default is just the variable we need
    else:
        hydros = hydro
    pymsesvars = []
    # Run through hydro variables
    for hydro in hydros:
        if hydro == "NH":
            pymsesvars += ["rho"]
        elif hydro == "rho":
            pymsesvars += ["rho"]
        elif hydro == "T":
            pymsesvars += ["rho","P","xHII","xHeII","xHeIII"]
        elif hydro[0] == "B":
            pymsesvars += ["B-left","B-right"]
        elif hydro[0] == "v": # Any velocity
            pymsesvars += ["vel"]
        elif hydro == "spd":
            pymsesvars += ["vel"]
        elif hydro == "Pram":
            pymsesvars += ["rho","vel"]
        elif hydro == "Lcool":
            pymsesvars += ["rho","P","vel","xHII","xHeII","xHeIII","NpHII","NpHeII","NpHeIII"]
        else:
            pymsesvars += [hydro]
        # Add extras
    for e in extra:
        if not e in pymsesvars:
            pymsesvars.append(e)
    # Make unique list
    pymsesvars = list(set(pymsesvars))
    print "Making AMR source with the variables", pymsesvars
    # Output pymses AMR source
    return ro.amr_source(pymsesvars)
    
def cmap(hydro):
    '''
    Returns cmap for this hydro variable 
    '''
    if hydro == "rho":
        return "RdPu_r"
    if hydro == "NH":
        return "RdPu_r"
    if hydro == "P":
        return "jet"
    if hydro == "T":
        return "jet"
    if "xH" in hydro:
        return "RdPu_r"
    if hydro == "Bmag":
        return "RdPu_r"
    if hydro[0] == "v": # any velocity
        return "RdPu_r"
    if hydro == "spd":
        return "RdPu_r"
    if hydro == "Pram":
        return "RdPu_r"
    if hydro == "Lcool":
        return "YlOrRd"
    # None of those? Return false
    return "linear"

def yscale(hydro):
    '''
    Returns whether this function takes a log scale or not
    '''
    if hydro == "rho":
        return "log"
    if hydro == "NH":
        return "log"
    if hydro == "P":
        return "log"
    if hydro == "T":
        return "log"
    if "xH" in hydro:
        return "linear"
    if hydro == "Bmag":
        return "log"
    if hydro[0] == "v": # any velocity
        return "linear"
    if hydro == "spd":
        return "linear"
    if hydro == "Pram":
        return "log"
    if hydro == "Lcool":
        return "log"
    # None of those? Return false
    return "linear"

def hydro_range(hydro):
    if hydro == "rho":
        return (0, 6)
    if hydro == "NH":
        return (20.5,23.0)
    if hydro == "P":
        return (None, None) # No idea
    if hydro == "T":
        return (3,6)
    if "xH" in hydro:
        return (-6,0)
    if hydro == "gpe":
        return (None, None)
    if hydro == "Bmag":
        return (None, None)
    print "Using default hydro_range for ", hydro
    return (None,None)
