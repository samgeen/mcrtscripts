'''
Functions defining hydro variables in RAMSES for use by pymses
Sam Geen, November 2014
'''

import numpy as np
from abc import ABCMeta, abstractmethod

from pymses.utils import constants as C

import sys
sys.path.append("../../Cooling/f2py")
usecooling = True
try:
    import rtcooling
except ImportError:
    print("rtcooling MODULE FAILED TO IMPORT, USING WITHOUT COOLING IN HYDROFUNCS")
    usecooling = False
    
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
    
    def __init__(self,label,scalebyunits,basevariables,colourmap,axisscale,axisrange=(None, None),surfacequantity=False,axisrescalefunc=None):
        self._label = label
        self._scalebyunits = scalebyunits
        self._basevariables = basevariables
        self._colourmap = colourmap
        self._axisscale = axisscale
        self._axisrange = axisrange
        self._surfacequantity = surfacequantity
        if axisrescalefunc is None:
            axisrescalefunc = lambda x: x
        self._axisrescalefunc = axisrescalefunc

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

    def AxisRescaleFunc(self):
        return self._axisrescalefunc
    
    def AxisRange(self):
        return self._axisrange

    def SurfaceQuantity(self):
        return self._surfacequantity

class AllHydros(object):
    def __init__(self):
        self._hydros = None
        # Global variables for the functions to use
        self._globals = {}
        self._Setup()

    def __setitem__(self,name,newHydro):
        self._hydros[name] = newHydro

    def __getitem__(self,name):
        if not name in self._hydros:
            return self._hydros["DEFAULTEMPTY"](name)
        return self._hydros[name]

    def __contains__(self,name):
        return name in self._hydros

    def AddGlobals(self,globalsdict):
        for key, value in globalsdict.items():
            self._globals[key] = value
    
    def _Setup(self):
        h = {}
        # Cheat Sheet for Hydro's constructor:
        # Hydro(label,scalebyunits,basevariables,colourmap,axisscale,axisrange=(None, None),surfacequantity=False)
        # About the lambda lambdas: pymses needs a function that takes dset. But we also need ro for units.
        #   So we make a first lambda that fills in ro, then we have a lambda dset we can pass
        #   I know, I know...
        # Cell size
        func = lambda ro: lambda dset: dset.get_sizes()*ro.info["unit_length"].express(C.cm)
        h["dx"] = Hydro("Cell size / cm",func,["rho"],"RdPu_r","linear",(None,None))
        # Mass density
        func = lambda ro: lambda dset: dset["rho"]*ro.info["unit_density"].express(C.g_cc)
        h["rho"] = Hydro("Density / g/cm$^{3}$",func,["rho"],"RdPu_r","log",(-25, -17))
        # Hydrogen number density
        func = lambda ro: lambda dset: dset["rho"]*ro.info["unit_density"].express(C.H_cc)
        h["nH"] = Hydro("Density / atoms/cm$^{3}$",func,["rho"],"RdPu_r","log",(0, 6))
        # Column (number) density
        func =  lambda ro: lambda dset: dset["rho"]*ro.info["unit_density"].express(C.H_cc)*ro.info["unit_length"].express(C.cm)
        h["NH"] = Hydro("$N_{\mathrm{H}}$ / cm$^{-2}$",func,["rho"],"RdPu_r","log",(20.5,24.0),surfacequantity=True)
        # Thermal pressure
        Pfunc = lambda ro: lambda dset: dset["P"]*ro.info["unit_pressure"].express(C.barye)
        h["P"] = Hydro("Pressure / ergs/cm$^3$",Pfunc,["P"],"YlOrRd","log",(None, None))
        # Thermal energy density (same as thermal pressure)
        h["Etherm"] = Hydro("Thermal energy density / ergs/cm$^3$",Pfunc,["P"],"YlOrRd","log",(-16.0,-7.0))
        # Kinetic energy density
        EKfunc = lambda ro: lambda dset: dset["rho"]*ro.info["unit_density"].express(C.g_cc)* \
                  np.sum(dset["vel"]**2,1)*(ro.info["unit_velocity"].express(C.cm/C.s))**2
        h["Ekin"] = Hydro("Kinetic energy density / ergs/cm$^3$",EKfunc,["rho","vel"],"PuBuGn","log",(-16.0, -7.0))
        # Kinetic energy / Thermal Energy
        EKperEThermfunc = lambda ro: lambda dset: dset["rho"]*ro.info["unit_density"].express(C.g_cc)* \
                  np.sum(dset["vel"]**2,1)*(ro.info["unit_velocity"].express(C.cm/C.s))**2 / \
                  (dset["P"]*ro.info["unit_pressure"].express(C.barye))
        h["EkinperEtherm"] = Hydro("$E_{kin} / E_{therm}$",EKperEThermfunc,["rho","vel","P"],"RdYlBu","log",(-3.0,3.0))
        # Temperature
        mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                     0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
        func = lambda ro: lambda dset: dset["P"]/dset["rho"]*ro.info["unit_temperature"].express(C.K)*mufunc(dset)
        h["T"] = Hydro("Temperature / K",func,["rho","P","xHII","xHeII","xHeIII"],"YlOrRd","log",(0,8))
        # Ionisation Fraction(s)
        h["xHII"] = Hydro("$x_{\mathrm{HII}}$",lambda ro: lambda dset: dset["xHII"],["xHII"],"PRGn","linear",(0,1))
        h["xHeII"] = Hydro("$x_{\mathrm{HeII}}$",lambda ro: lambda dset: dset["xHeII"],["xHeII"],"Reds","linear",(0,1))
        h["xHeIII"] = Hydro("$x_{\mathrm{HeIII}}$",lambda ro: lambda dset: dset["xHeIII"],["xHeIII"],"Reds","linear",(0,1))
        # Photon flux in photons / cm^2 / s
        h["NpFUV"] = Hydro("$F_{\mathrm{FUV}}$",lambda ro: lambda dset: dset["NpFUV"]*ro.info["unit_velocity"].express(C.cm/C.s),["NpFUV"],"Reds","log",(None,None))
        h["NpHII"] = Hydro("$F_{\mathrm{HII}}$",lambda ro: lambda dset: dset["NpHII"]*ro.info["unit_velocity"].express(C.cm/C.s),["NpHII"],"PRGn","log",(None,None))
        h["NpHeII"] = Hydro("$F_{\mathrm{HeII}}$",lambda ro: lambda dset: dset["NpHeII"]*ro.info["unit_velocity"].express(C.cm/C.s),["NpHeII"],"Reds","log",(None,None))
        h["NpHeIII"] = Hydro("$F_{\mathrm{HeIII}}$",lambda ro: lambda dset: dset["NpHeIII"]*ro.info["unit_velocity"].express(C.cm/C.s),["NpHeIII"],"Reds","log",(None,None))
        # Gravitational energy - NOTE: doesn't really work in Pymses, including for completeness
        # h["gpe"] = Hydro("Gravitational Potential Energy",func???,[???],"YlOrRd","log",(None,None))
        # Magnetic field strength
        def bfuncforsnap(axis):
            def howmanylayersareyouonmydude(ro):
                ud = ro.info["unit_density"].express(C.g_cc)
                ul = ro.info["unit_length"].express(C.cm) / ro.info["boxlen"]
                ut = ro.info["unit_time"].express(C.s)
                unit = np.sqrt(4.0*np.pi*ud*(ul/ut)**2)*1e6 # microGauss
                def bfunc(dset):
                    b = 0.5*(dset["B-left"]+dset["B-right"])
                    if axis == "all":
                        Bout = np.sqrt((b**2).sum(axis=1))*unit
                        return Bout
                    return b[:,axis]*unit
                return bfunc
            return howmanylayersareyouonmydude
        h["Bmag"] = Hydro("|B| / $\mu G$",bfuncforsnap("all"),["B-left","B-right"],"PuOr","log",(None, None))
        # B-fild vector (x)
        h["Bx"] = Hydro("B.$x$ / $\mu G$",bfuncforsnap(0),["B-left","B-right"],"PuBuGn","linear",(None, None))
        # B-field vector (y)
        h["By"] = Hydro("B.$y$ / $\mu G$",bfuncforsnap(1),["B-left","B-right"],"PuBuGn","linear",(None, None))
        # B-field vector (z)
        h["Bz"] = Hydro("B.$z$ / $\mu G$",bfuncforsnap(2),["B-left","B-right"],"PuBuGn","linear",(None, None))
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
        def vnonradfunc(ro):
            unit = ro.info["unit_velocity"].express(C.km/C.s)
            def findvnonrad(dset):
                try:
                    pts = dset.points
                except:
                    pts = dset.get_cell_centers()
                pos = pts+0.0
                pos[:,0] = pts[:,0]-self._globals["starpos"][0]
                pos[:,1] = pts[:,1]-self._globals["starpos"][1]
                pos[:,2] = pts[:,2]-self._globals["starpos"][2]
                rad = pos # Be careful here! Reference, not value copy
                dist = np.sqrt(np.sum(pos**2,1))
                for i in range(0,3):
                    rad[:,i] /= dist
                vrad = np.sum(rad*dset["vel"],1)*unit
                spd = np.sum(dset["vel"]**2,1)*unit
                vnonrad = spd - vrad
                return vnonrad
            return findvnonrad
        def vradfracfunc(ro):
            unit = ro.info["unit_velocity"].express(C.km/C.s)
            def findvrad(dset):
                pos = dset.points-0.5
                rad = pos # Be careful here! Reference, not value copy
                dist = np.sqrt(np.sum(pos**2,1))
                for i in range(0,3):
                    rad[:,i] /= dist
                vrad = np.sum(rad*dset["vel"],1)
                spd = np.sqrt(np.sum(dset["vel"]**2,1))
                return vrad / spd
            return findvrad
        h["vrad"] = Hydro("$v_{\mathrm{radial}}$ / km/s",vradfunc,["vel"],"BuGn","log",(None, None))
        h["vnonrad"] = Hydro("$v_{\mathrm{non-radial}}$ / km/s",vnonradfunc,["vel"],"BuGn","log",(None, None))
        # Gas speed
        spdfunc = lambda ro: lambda dset: np.sqrt(np.sum(dset["vel"]**2,1))*ro.info["unit_velocity"].express(C.km/C.s)
        h["spd"] = Hydro("Gas Speed / km/s",spdfunc,["vel"],"BuGn","log",(None, None))
        # Radial fraction
        h["vradfrac3"] = Hydro("$v_{rad} / |v|$",vradfracfunc,["vel"],"BuGn","linear",(None, None))
        # Velocity vector (x)
        func = lambda ro: lambda dset: dset["vel"][:,0]*ro.info["unit_velocity"].express(C.km/C.s)
        h["vx"] = Hydro("$v_{\mathrm{x}}$ / km/s",func,["vel"],"BuGn","linear",(None, None))
        # Velocity vector (y)
        func = lambda ro: lambda dset: dset["vel"][:,1]*ro.info["unit_velocity"].express(C.km/C.s)
        h["vy"] = Hydro("$v_{\mathrm{y}}$ / km/s",func,["vel"],"BuGn","linear",(None, None))
        # Velocity vector (z)
        func = lambda ro: lambda dset: dset["vel"][:,2]*ro.info["unit_velocity"].express(C.km/C.s)
        h["vz"] = Hydro("$v_{\mathrm{z}}$ / km/s",func,["vel"],"BuGn","linear",(None, None))
        # Vorticity - NOTE: THIS IS HANDLED BY SPECIAL CODE IN sliceMap !!!
        # Number refers to number of pixels to sample across
        h["vorticity1px"] = Hydro("$| \\nabla \\times \mathbf{v}|(1 px) / Myr^{-1}$",vradfracfunc,["vel"],"BuGn","log",(None, None))
        h["vorticity2px"] = Hydro("$| \\nabla \\times \mathbf{v}|(2 px) / Myr^{-1}$",vradfracfunc,["vel"],"BuGn","log",(None, None))
        h["vorticity4px"] = Hydro("$| \\nabla \\times \mathbf{v}|(4 px) / Myr^{-1}$",vradfracfunc,["vel"],"BuGn","log",(None, None))
        h["vorticity4px_timescale"] = Hydro("$1 / | \\nabla \\times \mathbf{v}|(4 px) / Myr$",vradfracfunc,["vel"],"BuGn_r","log",(None, None))
        # Ram pressure
        func = lambda ro: lambda dset: dset["rho"]*np.sum(dset["vel"]**2,1)*ro.info["unit_pressure"].express(C.barye)
        h["Pram"] = Hydro("Ram Pressure / ergs/cm$^3$",func,["rho","vel"],"YlOrRd","log",(None, None))
        # Cooling-related functions
        coolvars = ["rho","P","vel","xHII","xHeII","xHeIII","NpHII","NpHeII","NpHeIII"]
        if usecooling:
            # Cooling rate (energy density change per unit time)
            coolfunc = lambda ro: rtcooling.dedtOnCells(ro)
            h["Lcool"] = Hydro("L$_{cool}$ / erg/s/cm$^{-3}$",coolfunc,coolvars,"YlOrRd","log",(None, None))
            # Damkoehler number
            def damkoehlerfunc(ro):
                def finddamkoehlerremapped(dset):
                    etherm = Pfunc(ro)(dset)
                    dedt = coolfunc(ro)(dset)
                    dedt[dedt < 0.0] = 0.0
                    tcool = etherm / dedt
                    L = self._globals["Ldamkoehler"]
                    vnonrad = np.abs(vnonradfunc(ro)(dset) * np.sqrt(3.0/2.0))
                    D = L / vnonrad
                    D = D / tcool
                    # Map to a sigmoid for plotting
                    remap = D#np.arctan(10.0*(D-1))*2.0/np.pi
                    #print("DEBUG",L, vnonrad.min(),vnonrad.max(), tcool.min(),tcool.max())
                    return remap
                return finddamkoehlerremapped
            def damkoehlerrescalefunc(D):
                return np.arctan(5.0*(D-1))*2.0/np.pi
            h["Damkoehler4"] = Hydro("Negligible Turbulent Mixing $\leftrightarrow$ Strong Turbulent Mixing",damkoehlerfunc,coolvars,"RdBu","linear",(-1.0, 1.0),axisrescalefunc=damkoehlerrescalefunc)
        # EMPTY DEFAULT HYDRO TO PREVENT ERRORS
        h["DEFAULTEMPTY"] = lambda hydro: Hydro("",lambda ro: lambda dset: dset[hydro],coolvars,"Blues_r","linear",(None, None))
        # Done!
        self._hydros = h
        
allhydros = AllHydros()

def hydro_label(hydro):
    '''
    Returns a label for the given hydro variable in some default units
    '''
    global allhydros
    return allhydros[hydro].Label()

def scale_by_units(ro, hydro):
    '''
    Returns a lambda function that scales a dataset by some default units
    '''
    global allhydros
    return allhydros[hydro].ScaleByUnits(ro)

def axis_rescale_func(hydro):
    global allhydros
    return allhydros[hydro].AxisRescaleFunc()

def amr_source(ro, hydro,extra=[]):
    '''
    Load the necessary data from the hydro files
    '''
    global allhydros
    if type(hydro) == type("duck"):
        hydros = [hydro] # default is just the variable we need
    else:
        hydros = hydro
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
    print("Making AMR source with the variables", pymsesvars)
    # Output pymses AMR source
    return ro.amr_source(pymsesvars)

def surface_quantity(hydro):
    '''
    Used by pymses - if true, divides each pixel by the pixel surface area
    '''
    global allhydros
    return allhydros[hydro].SurfaceQuantity()
    
def cmap(hydro):
    '''
    Returns cmap for this hydro variable 
    '''
    global allhydros
    return allhydros[hydro].ColourMap()

def yscale(hydro):
    '''
    Returns whether this function takes a log scale or not
    '''
    global allhydros
    return allhydros[hydro].AxisScale()

def hydro_range(hydro):
    '''
    Returns the min/max values to plot/show as a tuple
    '''
    global allhydros
    return allhydros[hydro].AxisRange()
