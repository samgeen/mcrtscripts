'''
Calculates slice maps with pymses
Rebekka Bieri, Sam Geen 2018
'''

# Get all the startup modules for the project
from startup import *

# Various extra pymses visualisation imports
from pymses.utils import constants as C
import pymses.analysis.visualization as v
scop = v.ScalarOperator

from pymses.sources.ramses.output import *
from pymses.analysis.visualization import *

import skimage.transform

import stellars

# Axes up and across for each line of sight
ups = {'x':'z','y':'x','z':'y'}
acrosses = {'x':'y','y':'z','z':'x'}
lostoi = {"x":0, "y":1, "z":2}

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
        print('in HERE [1]')
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

    try:
        snap = snap.RawData()
    except:
        pass

    centre = np.zeros(3)+0.5
    if starC:
        stars = stellars.FindStellar(snap)
        boxlen = snap.info["boxlen"]
        centre[lostoi[los]] = np.array([stars.x[0], stars.y[0], stars.z[0]])[lostoi[los]]/boxlen
    up = ups[los]
    across = acrosses[los]

    size = np.zeros(2)+zoom


    def makeslice(snap,hydro):
        hydro_op = scop(hydrofuncs.scale_by_units(snap,hydro))
        slc = pymses.analysis.visualization.SliceMap(amr, cam, hydro_op, z=0.0)
        print("Made slice (min/max:", slc.min(), slc.max(), ")")
        return slc


    if not "vorticity" in hydro:
        cam  = v.Camera(center=centre, line_of_sight_axis=los, 
                        region_size=size, up_vector=up, 
                        map_max_size=IMSIZE, log_sensitive=True)

        slc = makeslice(snap,hydro)
    else:
        # MAKE VORTICITY MAPS YES THIS IS HORRIBLE
        # Get pixel size
        NEWIMSIZE = IMSIZE
        # Step across multiple pixels / cells?
        if hydro == "vorticity2px":
            NEWIMSIZE /= 2
        if hydro == "vorticity4px":
            NEWIMSIZE /= 4
        dxcam = zoom / float(NEWIMSIZE)
        # Make camera again in case
        cam  = v.Camera(center=centre, line_of_sight_axis=los, 
                    region_size=size, up_vector=up, 
                    map_max_size=NEWIMSIZE, log_sensitive=True)
        # dx in km (to match velocity units)
        # We pre-divide every slice map by this to make calculations easier later
        dxphys = dxcam * boxlen * pcincm / 1000.0
        # Get xyz in frame of image (ensure right-handed coordinate system)
        # We need this because we calculate d/dx etc in frame of image
        vx0 = makeslice(snap,"v"+across) / dxphys
        vy0 = makeslice(snap,"v"+up) / dxphys
        vz0 = makeslice(snap,"v"+los) / dxphys
        # Make new slice + dx
        cx = centre+0.0
        cx[lostoi[across]] += dxcam
        cam = v.Camera(center=cx, line_of_sight_axis=los, 
                    region_size=size, up_vector=up, 
                    map_max_size=NEWIMSIZE, log_sensitive=True)
        vxx = makeslice(snap,"v"+across) / dxphys
        vyx = makeslice(snap,"v"+up) / dxphys
        vzx = makeslice(snap,"v"+los) / dxphys
        # Make new slice + dy
        cy = centre+0.0
        cy[lostoi[up]] += dxcam
        cam = v.Camera(center=cy, line_of_sight_axis=los, 
                    region_size=size, up_vector=up, 
                    map_max_size=NEWIMSIZE, log_sensitive=True)
        vxy = makeslice(snap,"v"+across) / dxphys
        vyy = makeslice(snap,"v"+up) / dxphys
        vzy = makeslice(snap,"v"+los) / dxphys
        # Make new slice + dz
        cz = centre+0.0
        # HACK TEST
        cz[lostoi[los]] -= dxcam
        cam = v.Camera(center=cz, line_of_sight_axis=los, 
                    region_size=size, up_vector=up, 
                    map_max_size=NEWIMSIZE, log_sensitive=True)
        vxz = makeslice(snap,"v"+across) / dxphys
        vyz = makeslice(snap,"v"+up) / dxphys
        vzz = makeslice(snap,"v"+los) / dxphys
        # Make vorticity map
        vortx = (vzy - vz0) - (vyz - vy0) 
        vorty = (vxz - vx0) - (vzx - vz0) 
        vortz = (vyx - vy0) - (vxy - vx0) 
        # Find magnitude in Myr^-1
        slc = np.sqrt(vortx**2 + vorty**2 + vortz**2) * Myrins
        # Find turnover timescale?
        if "timescale" in hydro:
            slc = 1.0 / slc
        # Compare the eddy turnover speed (dxphys / curl V) to the bulk gas speed
        if "speedcompare" in hydro:
            spd = np.sqrt(vx0**2 + vy0**2 + vz0**2) * dxphys
            slc = dxphys * slc / Myrins / spd
        # Resize the output image if needed
        if NEWIMSIZE != IMSIZE:
            slc = skimage.transform.resize(slc, (IMSIZE, IMSIZE))
    return slc

_MapSliceHamu = Hamu.Algorithm(_MapSlice)
#_MapSliceHamu._force_replace_cache = True

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
        self._snap  = snap.RawData()
        self._los   = los
        self._hydro = hydro
        self._zoom  = zoom
        self._starC = starC
        if pixlength is None:
            # NOTE: boxlen should be in pc!!
            #pixlength = snap.info["boxlen"] * zoom / float(IMSIZE)
            pixlength = self._snap.info["boxlen"] / float(IMSIZE)
        self._pixlength = pixlength
        self._slice = None

    def getSliceMap(self):
        if self._slice is None:
            self._slice = _MapSliceHamu(self._snap.hamusnap,self._hydro,self._los,self._zoom, self._starC)
        return self._slice
