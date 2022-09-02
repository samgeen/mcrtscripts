'''
Make searchable octree data for RAMSES outputs
Sam Geen, August 2022
'''

from startup import *

import h5py

import rdmfile, sinks, stellars

from pymses.utils import constants as C 
from pymses.filters import CellsToPoints

MEMSAFE = False


class Octree:
    # Octree class that takes points referencing data at indices in a flat array
    # Simple implementation that allows top-down search and writing to file
    # Intended to be used by a GLSL shader to visualise the octree data
    # GLSL code for scanning octs
    """
    vec3 temppos = pos;
    int ioct = 1;
    int ix = 0;
    int iy = 0;
    int iz = 0;
    uint ic = 0;
    int ioct = 0;
    while (ioct > 0)
    {
        temppos *= 2;
        // Integer cell positions in the oct, should each be either 0 or 1
        ix = int(temppos.x);
        iy = int(temppos.y);
        iz = int(temppos.z);
        // Get linear coordinate in octree data array
        ic = ioct + 4*ix = 2*iy + iz;
        // Find either next position in the octree or a leaf cell array position
        ioct = octree[ic];
        // Check whether we have leaf cell data:
        // -ve leaf cell data to output
        // +ve location of next level oct in octree, keep searching
        if (ioct < 0)
        {
            return -ioct;
        }
        // Zoom in on the next oct by moving to the appropriate branch cell corner
        temppos -= vec3(ix,iy,iz);
    }
    """
    def __init__(self, xs,ys,zs, indices=None):
        self._children = None
        self._indices = None
        if indices is None:
            indices = np.arange(len(xs))
        self._Populate(xs,ys,zs,indices)

    def _Populate(self,xs,ys,zs,indices):
        # Populate octs with new indices until we reach a leaf oct
        xs *= 2
        ys *= 2
        zs *= 2
        ixs = xs.astype(np.int32)
        iys = ys.astype(np.int32)
        izs = zs.astype(np.int32)
        ics = 4*ixs + 2*iys + izs
        #print(ics)
        # Check if leafs reached using basic method (do we only have one sub-oct left?)
        #print(xs, ys,zs, len(xs))
        #if len(xs) < 8:
        #    print("Wrong number of leaf cells, something went wrong")
        #    raise ValueError
        #if ixs.max() > 1:
        #    print("ixs bad")
        #    raise ValueError
        if len(xs) == 1:
            # Set list indices in correct position order in oct
            self._indices = np.atleast_1d(indices)
        else:
            # Zoom on new octs
            #xs -= 0.5*(xs.min()+xs.max())
            #ys -= 0.5*(ys.min()+ys.max())
            #zs -= 0.5*(zs.min()+zs.max())
            xs -= ixs
            ys -= iys
            zs -= izs
            # Make and fill new octs
            self._children = []
            for ic in range(0,8):
                mask = ics == ic
                self._children.append(Octree(xs[mask],ys[mask],zs[mask],indices[mask]))

    def WriteToArray(self,startpos=0):
        # Write the octree to a linear array
        if self._children is None:
            # Write -ve indices of cell data to indicate it's a leaf
            return -self._indices
        else:
            # Create array for cell indices
            nchildren = 8
            newarr = np.zeros(nchildren)
            # Make a new start position for arrays to begin at
            newstart = startpos+nchildren
            for ic, child in enumerate(self._children):
                # Write the child's octree cell indices
                newarr[ic] = newstart
                # Iteratively populate the child array
                childarr = child.WriteToArray(newstart)
                newarr = np.concatenate((newarr,childarr))
                newstart += len(childarr)
            # Return the result
            return newarr


def makeoctree(snap,folder,hydros,pos,radius):
    # NOTE: Ignore pos, radius for a pure list
    # If MEMSAFE is set, try to run per hydro variable to save memory
    # This is of course slower and uses more disk reads
    if MEMSAFE and len(hydros) > 1:
        for hydro in hydros:
            run(snap,folder,[hydro],pos,radius)
    # Open snapshot
    ro = snap.RawData()
    print("Sampling grid ...",)
    amr = hydrofuncs.amr_source(ro,hydros)
    # Make grid
    lmin = ro.info["levelmin"]
    lsize = 512
    boxlen = ro.info["boxlen"]
    unitcm = ro.info["unit_length"].express(C.cm)/ro.info["boxlen"]
    boxlencm = ro.info["boxlen"]*unitcm
    # Sample for each hydro variable
    hydros = ["octarray","stars"]+hydros+["x","y","z","amrlevel"]
    cell_source = CellsToPoints(amr)
    samples = cell_source.flatten()
    cells = samples
    octree = Octree(cells.points[:,0],cells.points[:,1],cells.points[:,2])
    octarray = octree.WriteToArray()
    dtype = "f"
    filename = folder+"/"+"celloctree_"+str(snap.OutputNumber()).zfill(5)+".hdf5"
    print("Making", filename,"...",)
    f = h5py.File(filename,"w")
    for hydro in hydros:
        print(hydro,"...")
        if hydro not in "xyz amrlevel stars octarray":
            scale = hydrofuncs.scale_by_units(ro,hydro)
            hydrocube = scale(samples)
            #hydrocube = hydrocube.reshape([lsize,lsize,lsize])
        else:
            if hydro == "x":
                hydrocube = cells.points[:,0]*boxlencm
            if hydro == "y":
                hydrocube = cells.points[:,1]*boxlencm
            if hydro == "z":
                hydrocube = cells.points[:,2]*boxlencm
            if hydro == "octarray":
                dtype = "i"
                hydrocube = octarray
            if hydro == "amrlevel":
                dtype = "i"
                hydrocube = cells.get_sizes()
                # Get level from cell sizes
                # 1.0001 limits issues with int casting 0.9999->0
                hydrocube = -np.log2(hydrocube)*1.0001
                hydrocube = hydrocube.astype(np.int32)
                minlevel = hydrocube.min()
                maxlevel = hydrocube.max()
                # Make info file
                infoname = folder+"/infos/"+\
                    "/info_"+str(snap.OutputNumber()).zfill(5)+".yaml"
                makeinfofile(infoname,boxlencm,minlevel,maxlevel)
        if hydro == "stars":
            # Make star file
            makestarfiles(snap,folder)
        if hydro != "stars":
            minmax = (hydrocube.min(),hydrocube.max())
            print(hydro,"type",dtype," min/max=",minmax,"of",hydrocube.shape,"cells")
            dset = f.create_dataset(hydro,hydrocube.shape,dtype=dtype)
            dset[:] = hydrocube
    f.close()
    print("Done!")


def makeinfofile(infofilename,boxlencm,minref,maxref):
    # boxlen in cm, minimum level, maximum level
    template = '''boxlen: BOXLEN
    minref: MINREF
    maxref: MAXREF
    '''
    txt = template+""
    txt = txt.replace("BOXLEN",str(boxlencm))
    txt = txt.replace("MINREF",str(minref))
    txt = txt.replace("MAXREF",str(maxref))
    f = open(infofilename,"w")
    f.write(txt)
    f.close()
    
def makedir(folder):
    try:
        os.mkdir(folder)
    except:
        pass

def makestarfiles(snap,folder):
    # Massive star files
    stellar = stellars.FindStellar(snap.RawData())
    ro = snap.RawData()
    sink = sinks.FindSinks(ro)
    snaptimeMyr = snap.Time()*ro.info["unit_time"].express(C.Myr)
    # Remove mass from sinks in massive stars and make stellar positions
    stellarx = []
    stellary = []
    stellarz = []
    stellarvx = []
    stellarvy = []
    stellarvz = []
    scale_v = ro.info["unit_velocity"].express(C.cm/C.s)/1e5
    for sinkid, mass in zip(stellar.sinkid, stellar.mass):
        mask = sink.id == sinkid
        sink.mass[mask] -= mass
        stellarx.append(sink.x[mask])
        stellary.append(sink.y[mask])
        stellarz.append(sink.z[mask])
        stellarvx.append(sink.vx[mask])
        stellarvy.append(sink.vy[mask])
        stellarvz.append(sink.vz[mask])
    stellarx = np.array(stellarx)
    stellary = np.array(stellary)
    stellarz = np.array(stellarz)
    stellarvx = np.array(stellarvx)
    stellarvy = np.array(stellarvy)
    stellarvz = np.array(stellarvz)
    # Write the sink files
    filename = folder+"/sinks/sinks_"+str(snap.OutputNumber()).zfill(5)+".fits"
    rdm = rdmfile.RDMFile("sinks")
    rdm.AddPoints(sink.x,sink.y,sink.z,"position/pc")
    rdm.AddPoints(sink.vx*scale_v,sink.vy*scale_v,sink.vz*scale_v,"velocity")
    rdm.AddArray(sink.mass,"mass/Msun")
    rdm.Write(filename) 
    # Write the massive star files
    stellarage = snaptimeMyr - stellar.tcreated
    filename = folder+"/massivestars/massivestars_"+str(snap.OutputNumber()).zfill(5)+".fits"
    rdm = rdmfile.RDMFile("massivestars")
    rdm.AddPoints(stellarx,stellary,stellarz,"position/pc")
    rdm.AddPoints(stellarvx*scale_v,stellarvy*scale_v,stellarvz*scale_v,"velocity")
    rdm.AddArray(stellar.mass,"mass/Msun")
    rdm.AddArray(stellarage,"Age/Myr")
    # Get derived stellar properties
    nstars = len(stellarage)
    radii = np.zeros(nstars)
    Teffs = np.zeros(nstars)
    currmasses = np.zeros(nstars)
    surfacegs = np.zeros(nstars)

    for i, mass, age in zip(range(0,nstars), stellar.mass, stellarage):
        # Get stellar radius / cm
        ageins = age * Myrins
        radius = singlestar.star_radius(mass,ageins)
        radii[i] = radius
        # Get stellar effective temperature / K
        Teff = singlestar.star_teff(mass,ageins)
        Teffs[i] = Teff
        # Get stellar current mass / g
        energy, massloss = singlestar.star_winds(mass,0.0,ageins)
        currmassing = mass*Msuning - massloss
        currmasses[i] = currmassing
        # Get stellar surface gravity
        surfaceg = G * currmassing / radius**2
        surfacegs[i] = surfaceg
    rdm.AddArray(radii,"radius/cm")
    rdm.AddArray(Teffs,"Teff/K")
    rdm.AddArray(currmasses,"Current Mass/g")
    rdm.AddArray(surfacegs,"Surface Gravity/cm/s/s")
    rdm.Write(filename) 

def runforsim(sim,nouts=None,times=None,pos=None,radius=None,makecubes=True):
    #hydros = ["vx","vy","vz","rho","T","xHII","Bx","By","Bz","Bmag"]
    hydros = ["T","rho","xHII","P","xHeII","xHeIII","Bx","By","Bz","vx","vy","vz",
              "NpFUV","NpHII","NpHeII","NpHeIII"][::-1]
    simname = sim.Name()
    # Make folders
    simfolder = "../cubes/fits/"+simname
    makedir(simfolder)
    for hydro in hydros:
        makedir(simfolder+"/"+hydro)
    for cpos in "xyz":
        makedir(simfolder+"/"+cpos)
    makedir(simfolder+"/amrlevel")
    makedir(simfolder+"/infos")
    makedir(simfolder+"/sinks")
    makedir(simfolder+"/massivestars")
    run = makeoctree
    if times is None:
        for snap in sim.Snapshots():  
            run(snap,folder=simfolder,hydros=hydros,pos=pos,radius=radius)
    else:
        for time in times:
            snap = timefuncs.findsnapshot(sim,time)
            run(snap,folder=simfolder,hydros=hydros,pos=pos,radius=radius)
 
def testrun():
    # Use IMF2, winds + UV
    sim = hamusims["SEED1_35MSUN_CDMASK_WINDUV"]
    # Pick the last output - TODO, CHANGE TO SPECIFIC OUTPUT!!!
    #nout = snap.OutputNumber()
    #nouts = [nout]
    #snaps = sim.Snapshots()
    #outnums = [snap.OutputNumber() for snap in snaps]
    #print(outnums)
    #snap = snaps[np.where(np.array(outnums) == 53)]
    snap = sim.Snapshots()[50]
    myr = snap.RawData().info["unit_time"].express(C.Myr)
    #times = [snap.Time()*myr for snap in sim.Snapshots()] # np.linspace(0.0,1.0,11)+3.2 # Myr
    times = [(snap.Time()*myr,"Myr"),] # np.linspace(0.0,1.0,11)+3.2 # Myr
    #snap = sim.FindAtTime(times[0]/myr)
    #boxlen = snap.RawData().info["boxlen"]
    # Pick the first star
    #stars = stellars.FindStellar(snap)
    #pos = [stars.x[0],stars.y[0],stars.z[0]]
    #pos = np.array(pos)/boxlen
    #print(pos)
    #radius = 25.0
    #pos = np.zeros(3)+0.5
    radius = 0.25
    runforsim(sim,times=times,pos=pos,radius=radius,makecubes=False)

if __name__=="__main__":
    simnames = seedset
    for simname in simnames:
        sim = hamusims[simname]
        times = [(t,"MyrFirstStar") for t in [0.1,0.2,0.3]]
        runforsim(sim,times=times,makecubes=False)
