'''
Calculates ray-traced maps with pymses
RB I guess but heavily taken from Sam Geen 2018
'''

# Get all the startup modules for the project
from startup import *

# Various extra pymses visualisation imports
from pymses.utils import constants as C
import pymses.analysis.visualization as v
scop = v.ScalarOperator

from pymses.sources.ramses.output import *
from pymses.analysis.visualization import *

# Axes up and across for each line of sight
ups = {'x':'z','y':'x','z':'y'}
acrosses = {'x':'y','y':'z','z':'x'}

IMSIZE = 1024

class MaxOperator(Operator):
    def __init__(self, ro, hydro_op):
        #self._unit = ro.info["unit_temperature"].express(C.K)
        #def Tfunc(dset):
        #    mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
        #               0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
        #    T = dset["P"]/dset["rho"]*self._unit*mufunc(dset)
        #    return T
        d = {"h":  hydro_op}
        Operator.__init__(self, d, is_max_alos=True)

    def operation(self, int_dict):
        mapHydro = int_dict.values()[0]
        mask = (mapHydro <= 0.0)
        mapHydro[mask] = 0.0
        return mapHydro

    '''
def pymses_func(ro, hydro):
    if hydro == "rho":
        unit = ro.info["unit_density"].express(C.H_cc)
        #return scop(lambda dset: dset["rho"]*unit)
        return v.FractionOperator(lambda dset: dset["rho"]**2*unit**2,
                                  lambda dset: dset["rho"]*unit)
    if hydro == "P":
        unit = ro.info["unit_pressure"].express(C.barye)
        return scop(lambda dset: dset["P"]*unit)
    if hydro == "T":
        return MaxTempOperator(ro) 
    if "xH" in hydro:
        unit = 1.0
        return scop(lambda dset: dset[hydro]*unit)
    if hydro == "Bmag":
        def bmagfunc(dset):
            b = 0.5*(dset["B-left"]+dset["B-right"])
            # Magnitude of the 3-vector for each cell
            return np.sqrt((b**2).sum(axis=1))
        return scop(bmagfunc)
    # None of those? Return unitless
    sco = scop(lambda dset: dset[hydro])
    return sco

def hydro_range(hydro):
    if hydro == "rho":
        return (0,8)
    if hydro == "P":
        return (None, None) # No idea
    if hydro == "T":
        return (0,5)
    if "xH" in hydro:
        return (-6,0)
    if hydro == "gpe":
        return (None, None)
    if hydro == "Bmag":
        return (None, None)

def hydro_label(hydro):
    if hydro == "rho":
        return "Density / atoms/cm$^{3}$"
    if hydro == "P":
        return "Pressure / dyne"
    if hydro == "T":
        return "Temperature / K"
    if "xH" in hydro:
        return "Ionisation Fraction "+hydro
    if hydro == "gpe":
        return "Gravitational Potential Energy"
    if hydro == "Bmag":
        return "|B| (code units)"
'''

def _MapRayTrace(snap,hydro='rho',los='z',zoom=1.0,starC=False):
    #allhydros = ["rho","P","xHII","B-left","B-right","xHeII","xHeIII"]
    domax = False
    if "max" in hydro:
        hydro = hydro.replace("max","")
        domax = True
    amr = hydrofuncs.amr_source(snap,hydro)
    size = np.zeros(2)+zoom
    centre = np.zeros(3)+0.5
    if starC:
        stars = stellars.FindStellar(snap)
        boxlen = snap.info["boxlen"]
        centre = np.array([stars.x[0], stars.y[0], stars.z[0]])/boxlen
    up = ups[los]

    cam  = v.Camera(center=centre, line_of_sight_axis=los, 
                    region_size=size, up_vector=up, 
                    map_max_size=IMSIZE, log_sensitive=True)
    rt = pymses.analysis.visualization.raytracing.RayTracer(snap,amr.fields_to_read)
    def makeray(snap,hydro):
        if domax:
            op = hydrofuncs.scale_by_units(snap,hydro)
            hydro_op = MaxOperator(snap,op)
        else:
            if hydro == "rho" or hydro == "nH":
                func = hydrofuncs.scale_by_units(snap,hydro)
                hydro_op = v.FractionOperator(lambda dset: func(dset)*func(dset),
                                              lambda dset: func(dset))
            else:
                hydro_op = scop(hydrofuncs.scale_by_units(snap,hydro))
        im = rt.process(hydro_op,cam,surf_qty = hydrofuncs.surface_quantity(hydro))
        print "Made ray trace map for "+hydro+" (with min/max:", im.min(), im.max(), ")"
        return im
    im = makeray(snap,hydro)
    return im

_MapRayTraceHamu = Hamu.Algorithm(_MapRayTrace)
#_MapRayTraceHamu._force_replace_cache = True

class RayTraceMap(object):
    '''
    Ray density map 
    '''
    def __init__(self,snap,hydro,los='z',pixlength=None,zoom=1.0,starC=False):
        '''
        los                         - Line of sight (strings 'x','y','z')
        pixlength                   - Length of pixels in parsecs
        zoom                        - Factor to zoom (<1 = zoom, 1 = full box)
        '''
        self._snap  = snap.RawData()
        self._los   = los
        self._hydro = hydro
        self._zoom  = zoom
        if pixlength is None:
            # NOTE: boxlen should be in pc!!
            pixlength = self._snap.info["boxlen"] * zoom / float(IMSIZE)
        self._pixlength = pixlength
        self._raytrace = None
        self._starC = starC

    def getRaytraceMap(self):
        if self._raytrace is None:
            self._raytrace = _MapRayTraceHamu(self._snap.hamusnap,self._hydro,self._los,self._zoom,self._starC)
        return self._raytrace

    # Maybe in the future helpful
    #def returnMap(self):
    #    # Return a copy to prevent editing the main array
    #    im = self._coldens+0.0
    #    # Apply thresholds to image
    #    # Do this on-the-fly to allow remapping of limits
    #    if self._NHlow is not None:
    #        im[im < NHtoColDens(self._NHlow)] = 0.0
    #    if self._NHhigh is not None:
    #        coldenshigh = NHtoColDens(self._NHhigh)
    #        im[im > coldenshigh] = coldenshigh
    #    return im
