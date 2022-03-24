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
    boxlen = snap.info["boxlen"]
    levelmax = snap.info["levelmax"]
    dx = 2.0**(-levelmax)
    up = ups[los]
    across = acrosses[los]

    # Centre the image on the first star to form?
    if starC:
        stars = stellars.FindStellar(snap)
        centre[lostoi[across]] = np.array([stars.x[0], stars.y[0], stars.z[0]])[lostoi[across]]/boxlen
        centre[lostoi[up]] = np.array([stars.x[0], stars.y[0], stars.z[0]])[lostoi[up]]/boxlen
        centre[lostoi[los]] = np.array([stars.x[0], stars.y[0], stars.z[0]])[lostoi[los]]/boxlen

    size = np.zeros(2)+zoom


    def makeslice(snap,hydro):
        hydro_op = scop(hydrofuncs.scale_by_units(snap,hydro))
        slc = pymses.analysis.visualization.SliceMap(amr, cam, hydro_op, z=0.0)
        print("Made slice (min/max:", slc.min(), slc.max(), ")")
        return slc


    if not "vorticity" in hydro or "vdispersion" in hydro:
        cam  = v.Camera(center=centre, line_of_sight_axis=los, 
                        region_size=size, up_vector=up, 
                        map_max_size=IMSIZE, log_sensitive=True)

        slc = makeslice(snap,hydro)
    else:
        # MAKE VORTICITY MAPS YES THIS IS HORRIBLE
        # Get pixel size
        NEWIMSIZE = IMSIZE
        # Step across multiple pixels / cells?
        if "2px" in hydro:
            NEWIMSIZE = IMSIZE / 2
        if "4px" in hydro:
            NEWIMSIZE = IMSIZE/ 4
        #dxcam = zoom / float(NEWIMSIZE) * 0.5 # Undersample to prevent cell skipping effects
        dxcam = dx * float(IMSIZE) / float(NEWIMSIZE)
        #centre = centre+0.5*dxcam
        # Make camera again in case
        cam  = v.Camera(center=centre, line_of_sight_axis=los, 
                    region_size=size, up_vector=up, 
                    map_max_size=NEWIMSIZE, log_sensitive=True)
        # dx in km (to match velocity units)
        # We pre-divide every slice map by this to make calculations easier later
        dxphys = dxcam * boxlen * pcincm
        if "vorticity" in hydro:
            # Get xyz in frame of image (ensure right-handed coordinate system)
            # We need this because we calculate d/dx etc in frame of image
            vx0 = makeslice(snap,"v"+across)
            vy0 = makeslice(snap,"v"+up) 
            vz0 = makeslice(snap,"v"+los) 
            # Make new slice + dx
            cx = centre+0.0
            cx[lostoi[across]] += dxcam
            cam = v.Camera(center=cx, line_of_sight_axis=los, 
                           region_size=size, up_vector=up, 
                           map_max_size=NEWIMSIZE, log_sensitive=True)
            vxx = makeslice(snap,"v"+across) 
            vyx = makeslice(snap,"v"+up)
            vzx = makeslice(snap,"v"+los)
            # Make new slice + dy
            cy = centre+0.0
            cy[lostoi[up]] += dxcam
            cam = v.Camera(center=cy, line_of_sight_axis=los, 
                           region_size=size, up_vector=up, 
                           map_max_size=NEWIMSIZE, log_sensitive=True)
            vxy = makeslice(snap,"v"+across) 
            vyy = makeslice(snap,"v"+up) 
            vzy = makeslice(snap,"v"+los) 
            # Make new slice + dz
            cz = centre+0.0
            # HACK TEST
            cz[lostoi[los]] += dxcam
            cam = v.Camera(center=cz, line_of_sight_axis=los, 
                           region_size=size, up_vector=up, 
                           map_max_size=NEWIMSIZE, log_sensitive=True)
            vxz = makeslice(snap,"v"+across)
            vyz = makeslice(snap,"v"+up) 
            vzz = makeslice(snap,"v"+los) 
            # Get vorticity components in s^{-1}
            # x1000 to convert from km/s to cgs
            vortx = ((vzy - vz0) - (vyz - vy0)) / dxphys * 1000.0
            vorty = ((vxz - vx0) - (vzx - vz0)) / dxphys * 1000.0
            vortz = ((vyx - vy0) - (vxy - vx0)) / dxphys * 1000.0
            # N = 4 here for the mean and velocity dispersion
            vxmean = (vx0 + vxx + vxy + vxz) / 4.0
            vymean = (vy0 + vyx + vyy + vyz) / 4.0
            vzmean = (vz0 + vzx + vzy + vzz) / 4.0
            #vdisp = 0.0
            #vdisp += (vx0 - vxmean)**2 + (vxx - vxmean)**2 + (vxy - vxmean)**2 + (vxz - vxmean)**2
            #vdisp += (vy0 - vymean)**2 + (vyx - vymean)**2 + (vyy - vymean)**2 + (vyz - vymean)**2
            #vdisp += (vz0 - vzmean)**2 + (vzx - vzmean)**2 + (vzy - vzmean)**2 + (vzz - vzmean)**2
            #vdisp = np.sqrt(vdisp / 4.0) * 1000.0 # --> cm/s
            # Speed in cgs from km/s
            spd = np.sqrt(vx0**2 + vy0**2 + vz0**2) * 1000.0
            # Make vorticity map
            # Find magnitude in Myr^{-1}
            slc = np.sqrt(vortx**2 + vorty**2 + vortz**2) * Myrins 
            # Find turnover timescale?
            if "timescale" in hydro:
                slc = 1.0 / slc
            # Compare the eddy turnover speed (dxphys / curl V) to the bulk gas speed
            if "speedcompare" in hydro:
                slc = dxphys * slc / Myrins / spd
        # Make velocity dispersion
        if "vdispersion" in hydro:
            # Camera plane slice
            vx0 = makeslice(snap,"v"+across)
            vy0 = makeslice(snap,"v"+up) 
            vz0 = makeslice(snap,"v"+los)
            # +los slice
            cplus = centre+0.0
            cplus[lostoi[los]] += dxcam
            cam  = v.Camera(center=cplus, line_of_sight_axis=los, 
                            region_size=size, up_vector=up, 
                            map_max_size=NEWIMSIZE, log_sensitive=True)
            vxp = makeslice(snap,"v"+across)
            vyp = makeslice(snap,"v"+up) 
            vzp = makeslice(snap,"v"+los)
            # -los slice
            cminus = centre+0.0
            cminus[lostoi[los]] -= dxcam
            cam  = v.Camera(center=cminus, line_of_sight_axis=los, 
                            region_size=size, up_vector=up, 
                            map_max_size=NEWIMSIZE, log_sensitive=True)
            vxm = makeslice(snap,"v"+across)
            vym = makeslice(snap,"v"+up) 
            vzm = makeslice(snap,"v"+los)
            # Make thin lasagna of 3 slices around the middle image
            # Will have shape 3, NEWIMSIZE, NEWIMSIZE
            vxgrid = np.array([vxm,vx0,vxp])
            # Make statistics
            vxmean = (vx0 + vxp + vxm) / 3.0
            vymean = (vy0 + vyp + vym) / 3.0
            vzmean = (vz0 + vzp + vzm) / 3.0
            nim = NEWIMSIZE
            vdisp = np.zeros((nim,nim))
            nsamples = 27 # 3x3x3 around each pixel
            # i1 is for the 3-deep slice lasagna
            for i1 in [0,1,2]:
                # i2 is the slice x axis
                for i2 in [-1,0,1]:
                    # i3 is the slice y axis
                    for i3 in [-1,0,1]:
                        smid = (slice(1,nim-1),slice(1,nim-1))
                        sgrid = (i1,slice(1+i2,nim-1+i2),slice(1+i3,nim-1+i3))
                        vdisp[smid] += (vxgrid[sgrid] - vxmean[smid])**2
                        vdisp[smid] += (vygrid[sgrid] - vymean[smid])**2
                        vdisp[smid] += (vzgrid[sgrid] - vzmean[smid])**2
            # Note: this is the 3D velocity dispersion
            vdisp = np.sqrt(vdisp / float(nsamples))
            # Bulk speed
            spd = np.sqrt(vx0**2 + vy0**2 + vz0**2)
            # Make images
            slc = vdisp +0.0
            if "speedcompare" in hydro:
                slc = vdisp / spd
        # Resize the output image if needed
        if NEWIMSIZE != IMSIZE:
            slc = skimage.transform.resize(slc, (IMSIZE, IMSIZE))
    return centre[lostoi[across]], centre[lostoi[up]], slc

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
        self._cx = None
        self._cy = None
        if pixlength is None:
            # NOTE: boxlen should be in pc!!
            #pixlength = snap.info["boxlen"] * zoom / float(IMSIZE)
            pixlength = self._snap.info["boxlen"] / float(IMSIZE)
        self._pixlength = pixlength
        self._slice = None

    def getSliceMap(self):
        #if "vdispersion" in self._hydro or "vorticity" in self._hydro:
        #Hamu.GLOBALFORCEREPLACECACHE = True
        if self._slice is None:
            self._cx, self._cy, self._slice = _MapSliceHamu(self._snap.hamusnap,self._hydro,self._los,self._zoom, self._starC)
        #Hamu.GLOBALFORCEREPLACECACHE = False
        return self._cx, self._cy, self._slice
