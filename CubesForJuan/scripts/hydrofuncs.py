'''
Functions defining hydro variables in RAMSES for use by pymses
Sam Geen, November 2014
'''

import numpy as np
from abc import ABCMeta, abstractmethod

from pymses.utils import constants as C

#import sys
#sys.path.append("../../Cooling/f2py") 
#import rtcooling

# Base hydro variable class used to define the interface
# e.g. density, temperature, magnetic field strength, etc
class AbstractHydro(ABCMeta):
    __metaclass__ = ABCMeta
    
    def __init__(self):
        pass

    @abstractmethod
    def Label(self):
        '''
        Returns a string corresponding to a plot axis label (e.g. "Temperature / K")
        '''
        return None

    @abstractmethod
    def ScaleByUnits(self, snap):
        '''
        Returns a function that takes pymses cells and outputs scalar values
        '''
        return None

    @abstractmethod
    def BaseVariables(self):
        '''
        Returns a list of all the base hydro fields for this hydro object
        '''
        return None

    @abstractmethod
    def ColourMap(self):
        '''
        Pyplot / matplotlib colour map to use
        '''
        return None

    @abstractmethod
    def AxisScale(self):
        '''
        Should this be displayed in log or linear space?
        '''
        return None

    def AxisRange(self):
        '''
        What range of values should be shown? i.e. min/max on colour maps / axes
        NOTE: If AxisScale is log, these should be log units
        TODO: Change this, it introduces an extra dependency we don't need
        '''
        return (None, None)

    def SurfaceQuantity(self):
        '''
        Is this a surface quantity? i.e. column density, etc
        Not used by everything so defaults to false
        '''
        return False

# Allows compact creation of hydro variables
# Can also define hydro variables by inheriting from AbstractHydro
# TODO: Relink to AbstractHydro
# TODO: Add setter structures like in cmap to all fields
class Hydro(object):
    
    def __init__(self,label,scalebyunits,basevariables,colourmap,axisscale,axisrange=(None, None),surfacequantity=False):
        self._label = label
        self._scalebyunits = scalebyunits
        self._basevariables = basevariables
        self._colourmap = colourmap
        self._axisscale = axisscale
        self._axisrange = axisrange
        self._surfacequantity = surfacequantity

    def Label(self):
        return self._label
    
    def ScaleByUnits(self, snap):
        return self._scalebyunits(snap)
    
    def BaseVariables(self):
        return self._basevariables

    def ColourMap(self,cmap=None):
        if cmap is not None:
            self._colourmap = cmap
        return self._colourmap

    def AxisScale(self):
        return self._axisscale

    def AxisRange(self):
        return self._axisrange

    def SurfaceQuantity(self):
        return self._surfacequantity

class AllHydros(object):
    def __init__(self):
        self._hydros = None
        self._Setup()

    def __setitem__(self,name,newHydro):
        self._hydros[name] = newHydro

    def __getitem__(self,name):
        if not name in self._hydros:
            return self._hydros["DEFAULTEMPTY"]
        return self._hydros[name]

    def __contains__(self,name):
        return name in self._hydros

    def _Setup(self):
        h = {}
        # Cheat Sheet for Hydro's constructor:
        # Hydro(label,scalebyunits,basevariables,colourmap,axisscale,axisrange=(None, None),surfacequantity=False)
        # About the lambda lambdas: pymses needs a function that takes dset. But we also need ro for units.
        #   So we make a first lambda that fills in ro, then we have a lambda dset we can pass
        #   I know, I know...
        # Mass density
        func = lambda ro: lambda dset: dset["rho"]*ro.info["unit_density"].express(C.g_cc)
        h["rho"] = Hydro("Density / g/cm$^{3}$",func,["rho"],"RdPu_r","log",(-28, -19))
        # Hydrogen number density
        func = lambda ro: lambda dset: dset["rho"]*ro.info["unit_density"].express(C.H_cc)
        h["nH"] = Hydro("Density / atoms/cm$^{3}$",func,["rho"],"RdPu_r","log",(0, 6))
        # Column (number) density
        func =  lambda ro: lambda dset: dset["rho"]*ro.info["unit_density"].express(C.H_cc)*ro.info["unit_length"].express(C.cm)
        h["NH"] = Hydro("$N_{\mathrm{H}}$ / cm$^{-2}$",func,["rho"],"RdPu_r","log",(20.5,23.0),surfacequantity=True)
        # Thermal pressure
        func = lambda ro: lambda dset: dset["P"]*ro.info["unit_pressure"].express(C.barye)
        h["P"] = Hydro("Pressure / ergs/cm$^3$",func,["P"],"YlOrRd","log",(None, None))
        # Temperature
        mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                     0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
        func = lambda ro: lambda dset: dset["P"]/dset["rho"]*ro.info["unit_temperature"].express(C.K)*mufunc(dset)
        h["T"] = Hydro("Temperature / K",func,["rho","P","xHII","xHeII","xHeIII"],"YlOrRd","log",(0,8))
        # Ionisation Fraction(s)
        h["xHII"] = Hydro("$x_{\mathrm{HII}}$",lambda ro: lambda dset: dset["xHII"],["xHII"],"Reds","linear",(0,1))
        h["xHeII"] = Hydro("$x_{\mathrm{HeII}}$",lambda ro: lambda dset: dset["xHeII"],["xHeII"],"Reds","linear",(0,1))
        h["xHeIII"] = Hydro("$x_{\mathrm{HeIII}}$",lambda ro: lambda dset: dset["xHeIII"],["xHeIII"],"Reds","linear",(0,1))
        h["NpHII"] = Hydro("$Np_{\mathrm{HII}}$",lambda ro: lambda dset: dset["NpHII"],["NpHII"],"Reds","linear",(None,None))
        h["NpHeII"] = Hydro("$Np_{\mathrm{HeII}}$",lambda ro: lambda dset: dset["NpHeII"],["NpHeII"],"Reds","linear",(None,None))
        h["NpHeIII"] = Hydro("$Np_{\mathrm{HeIII}}$",lambda ro: lambda dset: dset["NpHeIII"],["NpHeIII"],"Reds","linear",(None,None))
        # Gravitational energy - NOTE: doesn't really work in Pymses, including for completeness
        # h["gpe"] = Hydro("Gravitational Potential Energy",func???,[???],"YlOrRd","log",(None,None))
        # Magnetic field strength
        def bfuncforsnap(axis):
            def howmanylayersareyouonmydude(ro):
                ud = ro.info["unit_density"].express(C.g_cc)
                ul = ro.info["unit_length"].express(C.cm)/ro.info["boxlen"]
                ut = ro.info["unit_time"].express(C.s)
                unit = np.sqrt(4.0*np.pi*ud*(ul/ut)**2)*1e6 # microGauss
                def bfunc(dset):
                    b = 0.5*(dset["B-left"]+dset["B-right"])
                    if axis == "all":
                        # NOTE!!! FOR VISUALISATION, axis=2, FOR OTHER ANALYSIS axis=1
                        # FIGURE OUT NEAT WAY TO SWITCH THIS
                        return np.sqrt((b**2).sum(axis=1))*unit
                    else:
                        return b[:,axis]*unit
                return bfunc
            return howmanylayersareyouonmydude
        h["Bmag"] = Hydro("$|\mathrm{B}|~/~\mu G$",bfuncforsnap("all"),["B-left","B-right"],"PuOr","log",(None, None))
        # B-field vector (x)
        h["Bx"] = Hydro("B.$x$ / $\mu G$",bfuncforsnap(0),["B-left","B-right"],"PuBuGn","log",(None, None))
        # B-field vector (y)
        h["By"] = Hydro("B.$y$ / $\mu G$",bfuncforsnap(1),["B-left","B-right"],"PuBuGn","log",(None, None))
        # B-field vector (z)
        h["Bz"] = Hydro("B.$z$ / $\mu G$",bfuncforsnap(2),["B-left","B-right"],"PuBuGn","log",(None, None))
        # Radial velocity
        def vradfunc(ro):
            unit = ro.info["unit_velocity"].express(C.km/C.s)
            def findvrad(dset):
                pos = dset.points-0.5
                rad = pos # Be careful here! Reference, not value copy
                dist = np.sqrt(np.sum(pos**2,1))
                for i in range(0,3):
                    rad[:,i] /= dist
                return np.sum(rad*dset["vel"],1)*unit
            return findvrad
        h["vrad"] = Hydro("$v_{\mathrm{radial}}$ / km/s",vradfunc,["vel"],"BuGn","log",(None, None))
        # Gas speed
        func = lambda ro: lambda dset: np.sqrt(np.sum(dset["vel"]**2,1))*ro.info["unit_velocity"].express(C.km/C.s)
        h["spd"] = Hydro("Gas Speed / km/s",func,["vel"],"BuGn","log",(None, None))
        # Velocity vector (x)
        func = lambda ro: lambda dset: dset["vel"][:,0]*ro.info["unit_velocity"].express(C.km/C.s)
        h["vx"] = Hydro("$v_{\mathrm{x}}$ / km/s",func,["vel"],"BuGn","log",(None, None))
        # Velocity vector (y)
        func = lambda ro: lambda dset: dset["vel"][:,1]*ro.info["unit_velocity"].express(C.km/C.s)
        h["vy"] = Hydro("$v_{\mathrm{y}}$ / km/s",func,["vel"],"BuGn","log",(None, None))
        # Velocity vector (z)
        func = lambda ro: lambda dset: dset["vel"][:,2]*ro.info["unit_velocity"].express(C.km/C.s)
        h["vz"] = Hydro("$v_{\mathrm{z}}$ / km/s",func,["vel"],"BuGn","log",(None, None))
        # Ram pressure
        func = lambda ro: lambda dset: dset["rho"]*np.sum(dset["vel"]**2,1)*ro.info["unit_pressure"].express(C.barye)
        h["Pram"] = Hydro("Ram Pressure / ergs/cm$^3$",func,["rho","vel"],"YlOrRd","log",(None, None))
        # Cooling rate (energy density change per unit time)
        #func = lambda ro: rtcooling.dedtOnCells(ro)
        #coolvars = ["rho","P","vel","xHII","xHeII","xHeIII","NpHII","NpHeII","NpHeIII"]
        #h["Lcool"] = Hydro("L$_{cool}$ / erg/s/cm$^{-3}$",func,coolvars,"YlOrRd","log",(None, None))
        # EMPTY DEFAULT HYDRO TO PREVENT ERRORS
        h["DEFAULTEMPTY"] = Hydro("",lambda dset: dset[hydro],[],"linear","linear",(None, None))
        # Done!
        self._hydros = h
        
allhydros = AllHydros()

def hydro_label(hydro):
    '''
    Returns a label for the given hydro variable in some default units
    '''
    global allhydros
    return allhydros[hydro.replace("max","")].Label()

def scale_by_units(ro, hydro):
    '''
    Returns a lambda function that scales a dataset by some default units
    '''
    global allhydros
    return allhydros[hydro.replace("max","")].ScaleByUnits(ro)

def amr_source(ro, hydro,extra=[]):
    '''
    Load the necessary data from the hydro files
    '''
    global allhydros
    if type(hydro) == type("duck"):
        hydros = [hydro] # default is just the variable we need
    else:
        hydros = hydro
    hydros = [h.replace("max","") for h in hydros]
    pymsesvars = []
    # Run through hydro variables
    for hydro in hydros:
        if hydro in allhydros:
            pymsesvars += allhydros[hydro].BaseVariables()
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

def surface_quantity(hydro):
    '''
    Used by pymses - if true, divides each pixel by the pixel surface area
    '''
    global allhydros
    return allhydros[hydro.replace("max","")].SurfaceQuantity()
    
def cmap(hydro):
    '''
    Returns cmap for this hydro variable 
    '''
    global allhydros
    return allhydros[hydro.replace("max","")].ColourMap()

def yscale(hydro):
    '''
    Returns whether this function takes a log scale or not
    '''
    global allhydros
    return allhydros[hydro.replace("max","")].AxisScale()

def hydro_range(hydro):
    '''
    Returns the min/max values to plot/show as a tuple
    '''
    global allhydros
    return allhydros[hydro.replace("max","")].AxisRange()
