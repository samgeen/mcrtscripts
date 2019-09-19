'''
Calculates slice maps with pymses
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

import stellars

# Axes up and across for each line of sight
ups = {'x':'z','y':'x','z':'y'}
acrosses = {'x':'y','y':'z','z':'x'}

IMSIZE = 1024

class MaxTempOperator(Operator):
    def __init__(self, ro):
        self._unit = ro.info["unit_temperature"].express(C.K)
        def Tfunc(dset):
            mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                       0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
            T = dset["P"]/dset["rho"]*self._unit*mufunc(dset)
            return T
        d = {"T": Tfunc}
        Operator.__init__(self, d, is_max_alos=True)

    def operation(self, int_dict):
        mapT = int_dict.values()[0]
        return mapT

def pymses_func(ro, hydro):
    if hydro == "rho":
        unit = ro.info["unit_density"].express(C.H_cc)
        return scop(lambda dset: dset["rho"]*unit)
    if hydro == "P":
        unit = ro.info["unit_pressure"].express(C.barye)
        return scop(lambda dset: dset["P"]*unit)
    if hydro == "T":
        mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                     0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
        unit = ro.info["unit_temperature"].express(C.K)
        return scop(lambda dset: dset["P"]/dset["rho"]*unit*mufunc(dset))
    if "xH" in hydro:
        unit = 1.0
        return scop(lambda dset: dset[hydro]*unit)
    if hydro == "Bmag":
        def bmagfunc(dset):
            b = 0.5*(dset["B-left"]+dset["B-right"])
            # Magnitude of the 3-vector for each cell
            return np.sqrt((b**2).sum(axis=1))
        return scop(bmagfunc)
    if hydro == "IRflux":
        #def irfluxfunc(dset):
        #    print 'in HERE'
        #    exit()
        #    TIR_Trap_op = dset["Pnontherm"]
        #    return TIR_Trap_op
        print 'in HERE [1]'
        #exit()
        return scop(lambda dset: dset["Pnontherm"])


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
    return (None,None)
    
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

def _MapSlice(snap,hydro='rho',los='z',zoom=1.0,starC=False):
    amr = hydrofuncs.amr_source(snap,hydro)

    if not starC:
        centre = np.zeros(3)+0.5
    else:
        stars = stellars.FindStellar(snap)
        boxlen = snap.info["boxlen"]
        centre = np.array([stars.x[0], stars.y[0], stars.z[0]])
        centre = centre / boxlen
    up = ups[los]

    size = np.zeros(2)+zoom

    cam  = v.Camera(center=centre, line_of_sight_axis=los, 
                    region_size=size, up_vector=up, 
                    map_max_size=IMSIZE, log_sensitive=True)

    def makeslice(snap,hydro):
        hydro_op = scop(hydrofuncs.scale_by_units(snap,hydro))
        slc = pymses.analysis.visualization.SliceMap(amr, cam, hydro_op, z=0.0)
        print "Made slice (min/max:", slc.min(), slc.max(), ")"
        return slc

    slc = makeslice(snap,hydro)
    return slc

_MapSliceHamu = Hamu.Algorithm(_MapSlice)

class SliceMap(object):
    '''
    Slice map 
    '''
    def __init__(self,snap,hydro,los='z',pixlength=None,zoom=1.0,starC=False):
        '''
        los                         - Line of sight (strings 'x','y','z')
        pixlength                   - Length of pixels in parsecs
        zoom                        - Factor to zoom (<1 = zoom, 1 = full box)
        '''
        self._snap  = snap
        self._los   = los
        self._hydro = hydro
        self._zoom  = zoom
        self._starC = starC
        if pixlength is None:
            # NOTE: boxlen should be in pc!!
            #pixlength = snap.info["boxlen"] * zoom / float(IMSIZE)
            pixlength = snap.info["boxlen"] / float(IMSIZE)
        self._pixlength = pixlength
        self._slice = _MapSliceHamu(self._snap.hamusnap,self._hydro,self._los,self._zoom, self._starC)

    def getSliceMap(self):
        return self._slice
