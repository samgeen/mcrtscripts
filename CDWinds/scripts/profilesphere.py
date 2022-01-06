
'''
Find the spherically averaged profiles of the cloud
Sam Geen, November 2014
'''

from startup import *

from pymses.utils import constants as C

rextent = None

def mkdir(path):
    try:
        os.mkdir(path)
    except OSError:
        pass # Directory already exists, whatevs
    
def powerlaw(snap,hydro,npoints=1000000,centre=[0.5,0.5,0.5],rcut=0.5,rinner=0.0,nfit=1.0):
    '''
    Calculates a power law fit r = factor * prof ** power
    '''
    r, prof = profileHamu(snap,hydro,npoints,centre,rcut)
    prof = prof[r > rinner]
    r = r[r > rinner]
    notnan = np.logical_not(np.isnan(r*prof))
    r = r[notnan]
    prof = prof[notnan]
    prof[prof == 0.0] = prof[prof > 0.0].min() # Safety net
    fit = np.polyfit(np.log10(r),np.log10(prof),deg=1)
    power,factor = fit
    return power, ((10.0**factor)/nfit)**(-1/power)

def profile(snap, hydro,npoints=1000000,centre=[0.5,0.5,0.5],rcut=0.5,hcut=None):
    '''
    snap - pymses snapshot
    hydro - name of hydro variable
    npoints - number of points to sample
    rcut - radius to cut at
    hcut - range of hydro variables to accept (None = don't cut)
    '''
    amr = hydrofuncs.amr_source(snap,hydro,extra=["rho"])
    radius = rcut
    nbins = 200
    # TODO: Figure out sampling on centre better?
    hfunc = hydrofuncs.scale_by_units(snap, hydro)
    sphere = pymses.utils.regions.Sphere(centre, radius)
    points = sphere.random_points(int(npoints))
    point_dset = pymses.analysis.sample_points(amr,points)
    # TODO: IMPLEMENT HCUT HERE
    weight_func = hfunc
    bins = np.linspace(0.0,1.0,nbins) * radius
    prof = pymses.analysis.bin_spherical(point_dset,centre, weight_func,
                                                bins, divide_by_counts=True)
    bins *= snap.info["unit_length"].express(C.pc)
    binsout = (bins[1:]+bins[:-1])*0.5
    return binsout,prof

profileHamu = Hamu.Algorithm(profile)

def column_density(snap,npoints=1000000,rcut=0.5):
    '''
    Calculate the spherically-averaged column density in H atoms/cm^2
    '''
    r,p = profileHamu(snap,"rho",npoints,rcut)
    # Rescale from pc to cm
    r *= snap.RawData().info["unit_length"].express(C.cm) / \
        snap.RawData().info["unit_length"].express(C.pc)
    # Get bin boundaries
    diffs = r[1:]-r[:-1]
    p = p[:-1] # Last bin probably doesn't matter, eh
    return np.sum(p*diffs)

def central_density(snap,npoints=1000000,rcut=0.5):
    r,p = profileHamu(snap,"rho",npoints,rcut)
    # Rescale from pc to cm
    r *= snap.RawData().info["unit_length"].express(C.cm) / \
        snap.RawData().info["unit_length"].express(C.pc)
    # Sample inner cells
    return np.nanmean(p[0:10])

def mean_density(snap,dummy=None,rinner=None):
    '''
    Find mean cloud density
    NOTE - rinner, rcut are dummy parameters
    '''
    global rextent
    r,p = profileHamu(snap,"rho")
    # Rescale from pc to cm
    pctocm = snap.RawData().info["unit_length"].express(C.cm) / \
        snap.RawData().info["unit_length"].express(C.pc)
    r *= pctocm
    # Sample inner cells
    boxlen = snap.RawData().info["boxlen"]
    if rextent is None:
        extent = 10
        p = p[0:extent]
        r = r[0:extent]

    else:
        print("R REXTENT", r, rextent)
        p = p[r <= rextent*pcincm]
        r = r[r <= rextent*pcincm]
        #extent = int(len(p)*rextent/(boxlen/2.0)+0.5)
    #print "USING EXTENT", extent, "FOR REXTENT:", rextent, "BOXLEN:", boxlen,
    #print "LEN ARRAY:", len(p)
    dr = np.diff(r)
    p = 0.5 * (p[1:] + p[:-1]) # mean at midpoint of dr
    r = 0.5 * (r[1:] + r[:-1])
    dm = p*4.0*np.pi*r**2*dr
    mass = np.sum(dm)
    dens = mass / (4.0/3.0 * np.pi * r[-1]**3)
    return dens


def plot(simname, hydro, npoints=1000000,firstonly=False):
    # Set up
    sim = Hamu.Simulation(simname)
    for snap in sim.Snapshots():
        r,p = profileHamu(snap,hydro,npoints)
        print(r.shape, p.shape)
        mkdir("../plots/profiles")
        mkdir("../plots/profiles/"+simname)
        path = "../plots/profiles/"+simname+"/"+hydro
        mkdir(path)
        # Plot
        plt.clf()
        plt.plot(r,p)
        plt.yscale(hydrofuncs.yscale(hydro))
        plt.xlabel("Radius / pc")
        plt.ylabel(hydrofuncs.hydro_label(hydro))
        stub = path+"/profile_"+str(snap.OutputNumber()).zfill(5)
        plt.savefig(stub+".png")
        # Save data
        np.savetxt(stub+".dat",np.array((r,p)).T,delimiter=", ")
        if firstonly:
            break

if __name__ == "__main__":
    sims = ["N48_M4_B02_F2","N48_M4_B02_F3",
            "N47_M4_B02","N48_M4_B00","N00_M4_B02",
            "N48_M4_B02","N48_M4_B02_HR","N49_M4_B02"]
    sims = ["N00_M4_B02_C8"]#,"N48_M4_B02_C2"]
    hydros = ["rho","vrad"]#["rho","vrad","spd","P","T","xHII"]
    for sim in sims:
        for hydro in hydros:
            plot(sim,hydro,firstonly=False)
