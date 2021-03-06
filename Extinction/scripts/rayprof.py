'''
Make profiles in rays positioned uniformly on a sphere
Sam Geen, November 2014
'''

from startup import *

import isomesh
from pymses.utils import constants as C
import profilesphere

rays = None
rads = None
nprofs = 100
nray = 1000
rextent = None

def mkdir(path):
    try:
        os.mkdir(path)
    except OSError:
        pass # Directory already exists, whatevs

def makerays(centre):
    '''
    Make rays
    Run once only automatically!
    '''
    global rays, rads, nprofs, nray
    # Get rays arranged uniformly on a sphere (using a spiral pattern)
    mesh = isomesh.Mesh(nprofs)
    rays = mesh.Vertices()
    rays = np.repeat(rays,nray,axis=0)
    # Scale rays so that they shine out from the centre of the box
    dr = 1.0/nray
    arrscale = np.arange(0,1.0,dr)*0.5
    arrscale = np.tile(arrscale, nprofs)
    rays[0:nprofs*nray,0] *= arrscale
    rays[0:nprofs*nray,1] *= arrscale
    rays[0:nprofs*nray,2] *= arrscale
    rads = np.sqrt(np.sum(rays**2,1))
    rays[:,0] += centre[0]
    rays[:,1] += centre[1]
    rays[:,2] += centre[2]

def makeprofs(snap, hydro,centres,amr=None):
    '''
    Make profiles
    '''
    global rays, rads, nprofs, nray
    if rays is None or rads is None:
        makerays([0,0,0])
    allrays = np.zeros((rays.shape[0]*len(centres),3))
    iray = 0
    lenray = nprofs*nray
    # Gather all rays
    for centre in centres:
        allrays[iray:iray+lenray,0] = rays[:,0]+centre[0]
        allrays[iray:iray+lenray,1] = rays[:,1]+centre[1]
        allrays[iray:iray+lenray,2] = rays[:,2]+centre[2]
        iray += lenray
    if amr is None:
        amr = hydrofuncs.amr_source(snap,hydro)
    dset = pymses.analysis.sample_points(amr,allrays)
    allp = hydrofuncs.scale_by_units(snap, hydro)(dset)
    rs = []
    profs = []
    iray = 0
    for centre in centres:
        r = rads[iray:iray+lenray]*snap.info["unit_length"].express(C.pc)
        p = allp[iray:iray+lenray]
        iray += lenray
        rs.append(r)
        profs.append(p)
    return rs, profs

makeprofsHamu = Hamu.Algorithm(makeprofs)

def sampleray(rads, profs):
    '''
    Sample a ray in the array of rays hurray
    '''
    ind = np.random.randint(0,nprofs)*nray
    rad = rads[ind:ind+nray]
    prof= profs[ind:ind+nray]
    return rad, prof

def medianprofile(snap,hydro,centre):
    global nray, nprofs,rays,rads
    #sim = Hamu.Simulation(simname)
    r,p = makeprofsHamu(snap.hamusnap,hydro,centre)
    radii = r[0:nray]
    profs = np.reshape(p,(nprofs, nray)).T # No idea
    percs = np.percentile(profs,[0,25,50,75,100],axis=1).T
    pmed = percs[0:nprofs,2]
    return radii, pmed

medianprofileHamu = Hamu.Algorithm(medianprofile)
# Compatability hack for other scripts
profileHamu = Hamu.Algorithm(medianprofile)

def makecolumn(snap,centres,amr):
    global nray, nprofs,rays,rads
    # Get column density around centre along each ray
    rs,profs = makeprofs(snap,"nH",centres,amr)
    radii = (rs[0])[0:nray]*pcincm # Convert radius from pc to cm
    dr = radii[1] - radii[0]
    columns = []
    for prof in profs:
        p = np.reshape(prof,(nprofs, nray)).T # No idea
        cols = np.sum(p,axis=0)*dr
        columns.append(cols)
        # Debug check to make sure I didn't mess up the index in the sum above
        if len(cols) != nprofs:
            print(len(cols), nprofs)
            print("OOPS NOT NPROFS")
            raise ValueError
    # Return NH
    return columns

    
def medianradius(snap,centre,xlim=0.1):
    '''
    Find the median ionisation front radius
    '''
    r,p = makeprofsHamu(snap,"xHII",centre)
    profs = np.reshape(p,(nprofs, nray)).T # No idea
    rlims = np.zeros((nprofs))
    for iprof in range(0,nprofs):
        prof = profs[iprof,:]
        if prof.max() > xlim:
            rlims[iprof] = r[profs < xlim]
    return np.percentile(rlims,[50])[0]


def cutprofile(snap,hydro,hcut,centre):
    '''
    Cut any rays where the profile is above a given value
    '''
    global nray, nprofs,rays,rads
    #sim = Hamu.Simulation(simname)
    r,p = makeprofsHamu(snap.hamusnap,hydrom,centre)
    radii = r[0:nray]
    profs = np.reshape(p,(nprofs, nray)).T # Reshape to 2D array
    # Find the maximum along each profile
    maxprofs = np.max(profs,axis=0)
    # Now select the profiles where we have only values below the cut
    belowcut = np.where(maxprofs < hcut)[0]
    profs = profs[:,belowcut]
    # Now find the mean of these profiles
    pmed = np.median(profs,axis=1)
    return radii, pmed
    
cutprofileHamu = Hamu.Algorithm(cutprofile)

def central_density(snap,dummy=None,rinner=None):
    '''
    Find central density
    NOTE - rinner, rcut are dummy parameters
    '''
    global rextent
    r,p = medianprofileHamu(snap,"rho")
    # Rescale from pc to cm
    pctocm = snap.RawData().info["unit_length"].express(C.cm) / \
        snap.RawData().info["unit_length"].express(C.pc)
    r *= pctocm
    # Sample inner cells
    boxlen = snap.RawData().info["boxlen"]
    if rextent is None:
        extent = 10
    else:
        extent = int(len(p)*rextent/boxlen)
    print("FINDING CENTRAL DENSITY USING ",extent, "CELLS")
    return np.nanmean(p[0:extent])

def mean_density(snap,centre,dummy=None,rinner=None):
    '''
    Find mean cloud density
    NOTE - rinner, rcut are dummy parameters
    '''
    global rextent
    r,p = medianprofileHamu(snap,"rho",centre)
    # Rescale from pc to cm
    pctocm = snap.RawData().info["unit_length"].express(C.cm) / \
        snap.RawData().info["unit_length"].express(C.pc)
    r *= pctocm
    # Sample inner cells
    boxlen = snap.RawData().info["boxlen"]
    if rextent is None:
        extent = 10
    else:
        extent = int(len(p)*rextent/boxlen)
    p = p[0:extent]
    r = r[0:extent]
    dr = np.diff(r)
    p = 0.5 * (p[1:] + p[:-1]) # mean at midpoint of dr
    r = 0.5 * (r[1:] + r[:-1])
    dm = p*4.0*np.pi*r**2*dr
    mass = np.sum(dm)
    dens = mass / (4.0/3.0 * np.pi * r[-1]**3)
    #print "FINDING CENTRAL DENSITY USING ",extent, "CELLS"
    return dens

def cloudEdge(snap,nHcut=3.0):
    r,p = makeprofsHamu(snap.hamusnap,"rho")
    radii = r[0:nray]
    profs = np.reshape(p,(nprofs, nray)).T # No idea
    # TODO: Replace for loop when less tired
    edges = np.zeros((nprofs))
    badcount = 0
    for iprof in range(0,nprofs):
        prof = profs[:,iprof]
        try:
            edges[iprof] = r[prof < nHcut].min()
        except ValueError:
            badcount += 1
            #print prof.min(), prof.max()
            edges[iprof] = 0.0
    print("Number of edges not found", badcount, "of", nprofs)
    edges = np.sort(edges)
    return edges

cloudEdgeHamu = Hamu.Algorithm(cloudEdge)    

def plot(simname,hydro,centre):
    global nray, nprofs,rays,rads
    sim = Hamu.Simulation(simname)
    for snap in sim.Snapshots():
        # Setup
        mkdir("../plots/rays")
        mkdir("../plots/rays/"+simname)
        path = "../plots/rays/"+simname+"/"+hydro
        mkdir(path)
        # Process
        r,p = makeprofsHamu(snap,hydro,centre)
        radii = r[0:nray]
        profs = np.reshape(p,(nprofs, nray)).T # No idea
        percs = np.percentile(profs,[0,25,50,75,100],axis=1).T
        pmin = percs[0:nprofs,0]
        p25p = percs[0:nprofs,1]
        pmed = percs[0:nprofs,2]
        p75p = percs[0:nprofs,3]
        pmax = percs[0:nprofs,4]
        # Plot
        plt.clf()
        # Absolute range and IQR shaded areas
        print(radii.shape)
        print(pmax.shape)
        taba = np.concatenate((radii,radii[::-1])) # there and back again
        minmax = np.concatenate((pmin,pmax[::-1]))
        iqr = np.concatenate((p25p,p75p[::-1]))
        plt.fill(taba,minmax,alpha=0.33,edgecolor='none',facecolor='r')
        plt.fill(taba,iqr,   alpha=0.33,edgecolor='none',facecolor='r')
        plt.plot(radii,pmed,"r")
        # Sample random profiles to overlay
        for i in range(0,3):
            rs, ps = sampleray(r,p)
            plt.plot(rs,ps,"k")
        plt.yscale(hydrofuncs.yscale(hydro))
        plt.xlabel("Radius / pc")
        plt.ylabel(hydrofuncs.hydro_label(hydro))
        plt.savefig(path+"/profile_"+str(snap.OutputNumber()).zfill(5)+".png")

def plotgradient(simname,hydro,time,centre,label=None,xlims=None,powfile=None,suffix=None):
    '''
    Like prof but with a higher gradient
    '''
    global nray, nprofs,rays,rads
    sim = Hamu.Simulation(simname)
    # HACK
    print("Running for sim", simname, "hydro", hydro)
    snaps = sim.Snapshots()
    tcode = time * Myrins / unit_t
    if not time is None:
        snaps = [sim.FindAtTime(tcode)]
    for snap in snaps: 
        # Setup
        mkdir("../plots/gradrays")
        mkdir("../plots/gradrays/"+simname)
        path = "../plots/gradrays/"+simname+"/"+hydro
        mkdir(path)
        # Process
        r,p = makeprofsHamu(snap,hydro,centre)
        radii = r[0:nray]
        profs = np.reshape(p,(nprofs, nray)).T # No idea
        grad = np.arange(0,51)*2.0
        percs = np.percentile(profs,grad,axis=1).T
        #pmin = percs[0:nprofs,0]
        #p25p = percs[0:nprofs,1]
        #pmed = percs[0:nprofs,2]
        #p75p = percs[0:nprofs,3]
        #pmax = percs[0:nprofs,4]
        # Plot
        plt.clf()
        # Shade areas
        taba = np.concatenate((radii,radii[::-1])) # there and back again
        for i in np.arange(0,25):
            there = percs[0:nprofs,i]
            back = percs[0:nprofs,50-i]
            pcontour = np.concatenate((there,back[::-1]))
            plt.fill(taba,pcontour,alpha=0.04,edgecolor='none',facecolor='r')
        # Make median line
        plt.plot(radii,percs[0:nprofs,25],"k")
        # Plot the mean line 
        nfit = 1e2
        rsph,psph = profilesphere.profileHamu(snap,hydro,10000000,centre,rcut=0.1)
        power, factor = profilesphere.powerlaw(snap,hydro,10000000,centre,rcut=0.1,nfit=nfit)
        print("POWER LAW:", power, factor)
        outnum = snap.OutputNumber()
        if powfile is not None:
            powfile.write(simname+" "+str(outnum)+" w, r0, nfit:"+str(power)+" "+str(factor)+" "+str(nfit)+"\n ")
            powfile.flush()
        # Mask out parts of the mean line with no values
        mask = psph != 0
        rsph = rsph[mask]
        psph = psph[mask]
        plt.plot(rsph,psph,"k--")
        # Make min/max lines
        plt.plot(radii,percs[0:nprofs,0], "k",alpha=0.3)
        plt.plot(radii,percs[0:nprofs,50],"k",alpha=0.3)
        # Plot the 25 and 75% lines
        plt.plot(radii,percs[0:nprofs,25//2], "k:")
        plt.plot(radii,percs[0:nprofs,75//2],"k:")
        # Sample random profiles to overlay
        #for i in range(0,3):
        #    rs, ps = sampleray(r,p)
        #    plt.plot(rs,ps,"k")
        plt.yscale(hydrofuncs.yscale(hydro))
        if hydro == "rho":
            plt.ylim([3e-1,1e7])
        if hydro == "vrad":
            plt.ylim([-15,25])
        plt.xlabel("Radius / pc")
        plt.ylabel(hydrofuncs.hydro_label(hydro))
        if xlims is not None:
            plt.xlim(xlims)
        plt.xscale("log")
        # Add text label if one exists:
        if not label is None:
            xlims = plt.gca().get_xlim()
            ylims = plt.gca().get_ylim()
            xr = xlims[1] - xlims[0]
            yr = ylims[1] - ylims[0]
            plt.text(xlims[1]-0.1*xr,
                     ylims[1]*0.5,
                     label,
                     horizontalalignment='right',
                     verticalalignment='top')
        if suffix is None:
            suffix = ""
        plt.savefig(path+"/profile_"+suffix+str(snap.OutputNumber()).zfill(5)+".pdf",rasterized=True,dpi=200)

def FindStarPositions(snap):
    stellar = stellars.FindStellar(snap.hamusnap)
    positions = []
    for istellar in range(0,len(stellar.mass)):
        sinkid = stellar.sinkid[istellar]
        sink = sinks.FindSinks(snap)
        sinkloc = np.where(sink.id == sinkid)
        boxlen = snap.info["boxlen"]
        pos = np.array([sink.x[sinkloc],sink.y[sinkloc],sink.z[sinkloc]])/boxlen
        positions.append(pos)
    return positions, stellar.mass

def _FindColumnsInSnap(snap):
    positions, masses = FindStarPositions(snap)
    amr = hydrofuncs.amr_source(snap,"rho")
    columns = [makecolumn(snap,centre,amr) for centre in positions]
    return columns, masses
FindColumnsByStarInSnap = Hamu.Algorithm(_FindColumnsInSnap)

# This is *much* more efficient where nsink << nstellar
def _FindColumnsInSnapBySink(snap):
    sink = sinks.FindSinks(snap)
    stellar = stellars.FindStellar(snap)
    amr = hydrofuncs.amr_source(snap,"rho")
    boxlen = snap.info["boxlen"]
    posns = []
    sinkids = []
    for isink in range(0,len(sink.id)):
        sinkid = sink.id[isink]
        pos = np.array([sink.x[isink],sink.y[isink],sink.z[isink]])/boxlen
        posns.append(pos)
        sinkids.append(sinkid)
    columns = []
    if len(posns) > 0:        
        columns = makecolumn(snap,posns,amr)
    return columns, sinkids
FindColumnsBySinkInSnap = Hamu.Algorithm(_FindColumnsInSnapBySink)

def FindColumnsInSim(simname):
    DEBUG = False
    sim = hamusims[simname]
    for snap in sim.Snapshots():
        columns, masses = FindColumnsBySinkInSnap(snap)
        if DEBUG and len(masses) > 0:
            print(columns, masses)
            return

# TODO: Work out how to plot this for initial tests
# TODO: Instead pick one LOS and draw rays from each star out to that

if __name__ == "__main__":
    # Make figure
    #powfile = open("../plots/powerlaws.txt","w")
    
    for simname in [allsims[1:]]:
        FindColumnsInSim(simname)
        #time = 3.38 # Time of first sink formation
        #starpos = findstarpos(simname,time)
        #plotgradient(simname,"rho",time,starpos,"Star's Position",
        #             xlims=[0.03,25.0],powfile=powfile,suffix="starpos")
        #plotgradient(simname,"rho",time,np.array([0.5,0.5,0.5]),"Geometric Centre",
        #             xlims=[0.03,25.0],powfile=powfile,suffix="centre")
    #powfile.close()
