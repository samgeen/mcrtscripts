'''
Find maps of column density from a snapshot
Sam Geen, June 2016
'''

# Get all the startup modules for the project
from startup import *

# Various extra pymses visualisation imports
from pymses.utils import constants as C
import pymses.analysis.visualization as v
scop = v.ScalarOperator

# Astropy stuff for PSF
import astropy.convolution as conv

# Axes up and across for each line of sight
ups = {'x':'z','y':'x','z':'y'}
acrosses = {'x':'y','y':'z','z':'x'}

IMSIZE = 1024

def NH_op(snap):
    unit = snap.info["unit_density"].express(C.g_cc) * \
           snap.info["unit_length"].express(C.cm)
    return scop(lambda dset: dset["rho"]*unit)



def _MapColDens(snap,los='z'):
    # NOTE!! FINDS COLUMN DENSITY IN g cm^2
    hydro = "rho"
    amr = hydrofuncs.amr_source(snap,hydro)
    hydro_op = NH_op(snap)
    size = np.zeros(2)+1.0
    centre = np.zeros(3)+0.5
    up = ups[los]

    cam  = v.Camera(center=centre, line_of_sight_axis=los, 
                    region_size=size, up_vector=up, 
                    map_max_size=IMSIZE, log_sensitive=True)
    rt = pymses.analysis.visualization.raytracing.RayTracer(snap,["rho"])
    def makeray(snap,hydro,dosurf=True):
        # Note: dosurf gives cumulative (column) values rather than maxima
        hydro_op = NH_op(snap)
        im = rt.process(hydro_op,cam,surf_qty=dosurf)
        print "Made map with (min/max:", im.min(), im.max(), ")"        
        return im
    im = makeray(snap,hydro)
    return im

_MapColDensHamu = Hamu.Algorithm(_MapColDens)

class DensityMap(object):
    '''
    Map of densities with ability to convert & assign thresholds
    '''
    def __init__(self,snap,los='z',NHlow=None,NHhigh=None,
                 columndensitymap=None,pixlength=None,PSFsize=None):
        '''
        los                         - Line of sight (strings 'x','y','z')
        NHlow                       - Values below this are set to zero
        NHhigh                      - Values above this are set to NHhigh
        columndensitymap (optional) - Image to input
                 Note: columndensitymap SHOULD BE IN g cm^2, i.e. MASS UNITS
                 USE FUNCTIONS IN startup.py TO REMAP TO THESE UNITS
        pixlength                   - Length of pixels in parsecs
        PSFsize                     - Size of the PSF in parsecs (consistent with pixlen)
        '''
        self._snap = snap
        self._los = los
        if pixlength is None:
            # NOTE: boxlen should be in pc!!
            pixlength = snap.info["boxlen"] / float(IMSIZE)
        self._pixlength = pixlength
        self._PSFsize = PSFsize
        self._NHlow = NHlow
        self._NHhigh = NHhigh
        self._coldens = None
        if columndensitymap is None:
            self._SetupMap()
        else:
            self._coldens = columndensitymap

    def _Clone(self):
        newmap = self._coldens+0.0 # Good old numpy copy
        return DensityMap(self._snap,self._los,
                          self._NHlow,self._NHhigh,
                          columndensitymap=newmap,
                          pixlength=self._pixlength)

    def Flatten(self):
        newmap = self._Clone()
        newmap._coldens = newmap._coldens.flatten()
        newmap._coldens = np.sort(newmap._coldens)
        return newmap

    def _SetupMap(self):
        self._coldens = _MapColDensHamu(self._snap.hamusnap,self._los)
        # Apply PSF to image
        if self._PSFsize is not None:
            width = self._PSFsize / self._pixlength
            print "Convolving with a PSF of length", width, "pixels"
            kernel = conv.Gaussian2DKernel(width)
            oldcoldens = self._coldens
            self._coldens = conv.convolve(oldcoldens,kernel)

    def _PixelLength(self):
        return self._pixlength

    def NHlow(self,newlow):
        # Assign new lower NH limit
        self._NHlow = newlow

    def NHhigh(self,newhigh):
        # Assing new upper NH limit
        self._NHhigh = newhigh

    def ColDens(self):
        # Return a copy to prevent editing the main array
        im = self._coldens+0.0
        # Apply thresholds to image
        # Do this on-the-fly to allow remapping of limits
        if self._NHlow is not None:
            im[im < NHtoColDens(self._NHlow)] = 0.0
        if self._NHhigh is not None:
            coldenshigh = NHtoColDens(self._NHhigh)
            im[im > coldenshigh] = coldenshigh
        return im
    
    def NH(self):
        return ColDenstoNH(self.ColDens())
    
    def NH2(self):
        return ColDenstoNH2(self.ColDens())
    
    def Ak(self):
        return ColDenstoAk(self.ColDens())
    
    def Mass(self):
        dens = self.ColDens()
        pixlen = self._PixelLength()
        umass = (pcincm)**2 / Msuning
        masses = pixlen**2 * dens * umass
        return masses

    def Area(self):
        coldens = self.ColDens().flatten()
        coldens = coldens[coldens > 0.0]
        areainpix = float(len(coldens))
        # Physical length of one pixel
        pixlen = self._PixelLength()
        # Physical area
        cloudarea = np.sum(areainpix * pixlen**2)
        return cloudarea

    def Radius(self):
        return np.sqrt(self.Area()/np.pi)
    
    # Plotting routines
    def CumulativeMass(self,xfunc="NH"):
        im = self.Flatten()
        x = eval("im."+xfunc+"()")
        y = im.Mass()
        x = x[::-1]
        y = np.cumsum(y[::-1])
        return x, y

    def MassRadius(self):
        '''
        Mass vs radius for various Ak limits as per Lombardi+ 2010
        '''
        Aklow = 0.1
        Akhigh = 5.0
        aks = np.arange(Aklow,Akhigh+0.00001,0.2)
        lak = len(aks)
        radii = np.zeros(lak)
        masses = np.zeros(lak)
        for i in range(0,lak):
            self.NHlow(AktoNH(aks[i]))
            radii[i] = self.Radius()
            masses[i] = np.sum(self.Mass())
        # Done!
        return radii, masses

def CumulativeMass(snap,los,func,NHlow=None,NHhigh=None,imageIn=None,
                   pixlen=None,PSFsize=None):
    im = DensityMap(snap,los,NHlow,NHhigh,imageIn,pixlen,PSFsize)
    return im.CumulativeMass(func)
CumulativeMassHamu = Hamu.Algorithm(CumulativeMass)

def MassRadius(snap,los,NHlow=None,NHhigh=None,imageIn=None,
               pixlen=None,PSFsize=None):
    im = DensityMap(snap,los,NHlow,NHhigh,imageIn,pixlen,PSFsize)
    return im.MassRadius()
MassRadiusHamu = Hamu.Algorithm(MassRadius)

        
