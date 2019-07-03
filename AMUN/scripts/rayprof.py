'''
Make profiles in rays positioned uniformly on a sphere
Sam Geen, November 2014
'''

from startup import *

import isomesh
from pymses.utils import constants as C
import profilesphere

import scipy.optimize

rays = None
rads = None
nprofs = 10000
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

def makeprofs(snap, hydro,centre):
    '''
    Make profiles
    '''
    global rays, rads, nprofs, nray
    #if rays is None or rads is None:
    makerays(centre)
    amr = hydrofuncs.amr_source(snap,hydro)
    #import pdb; pdb.set_trace()
    dset = pymses.analysis.sample_points(amr,rays)
    #if dset.reorder_indices is not None:
    #    p = dset[hydro][dset.reorder_indices]
    #else:
    #    p = dset[hydro]
    p = hydrofuncs.scale_by_units(snap, hydro)(dset)
    r = rads*snap.info["unit_length"].express(C.pc)
    return r, p

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

def medianpowerlaw(snap,hydro,npoint=1000000,centre=[0.5,0.5,0.5],router=None,rinner=0.0,rfit=1.0):
    '''
    Calculates a power law fit r = factor * prof ** power 
    '''
    def func_powerlaw(x, m, c):
        return x**m * c
    r, prof = medianprofileHamu(snap,hydro,centre)
    if rinner > 0.0:
        prof = prof[r > rinner]
        r = r[r > rinner]
    if router is not None:
        prof = prof[r < router]
        r = r[r < router]
    notnan = np.logical_not(np.isnan(r*prof))
    r = r[notnan]
    prof = prof[notnan]
    prof[prof == 0.0] = prof[prof > 0.0].min() # Safety net
    #result = scipy.optimize.curve_fit(func_powerlaw,r,prof)
    #import pdb; pdb.set_trace()
    #power, factor = result[0]
    #import pdb; pdb.set_trace()
    fit = np.polyfit(np.log10(r),np.log10(prof),deg=1)
    power,factor = fit
    # n = n0 r0^w r^-w
    # --> log n = log n0 + w log r0 - w log r
    # --> m = -w, c = log n0 + w logr0
    # --> n0 = 10.0**(c - w log r0)
    #return power, ((10.0**factor)/nfit)**(-1/power)
    return power, 10.0**(factor + power * np.log10(rfit))

# Compatability hack for other scripts
profileHamu = Hamu.Algorithm(medianprofile)
    
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
    r,p = makeprofsHamu(snap.hamusnap,hydro,centre)
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
    print "FINDING CENTRAL DENSITY USING ",extent, "CELLS"
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
    print "Number of edges not found", badcount, "of", nprofs
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
        print radii.shape
        print pmax.shape
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
    print "Running for sim", simname, "hydro", hydro
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
        rsph,psph = profilesphere.profileHamu(snap,hydro,10000000,centre,rcut=0.1)
        #power, factor = medianpowerlaw(snap,hydro,10000000,centre,rinner=0.1,router=10.0,nfit=nfit)
        #print "POWER LAW:", power, factor
        outnum = snap.OutputNumber()
        if powfile is not None:
            rfit = 1.0
            power, factor = medianpowerlaw(snap,hydro,10000000,centre,rinner=0.1,router=1.0,rfit=rfit)
            powfile.write(simname+" "+str(outnum)+"INNER w, n0, rfit:"+str(power)+" "+str(factor)+" "+str(rfit)+"\n ")
            rfit = 10.0
            power, factor = medianpowerlaw(snap,hydro,10000000,centre,rinner=1.0,router=10.0,rfit=rfit)
            powfile.write(simname+" "+str(outnum)+"OUTER w, n0, rfit:"+str(power)+" "+str(factor)+" "+str(rfit)+"\n ")
            rfit = 10.0
            power, factor = medianpowerlaw(snap,hydro,10000000,centre,rinner=0.1,router=10.0,rfit=rfit)
            powfile.write(simname+" "+str(outnum)+"ALL w, n0, rfit:"+str(power)+" "+str(factor)+" "+str(rfit)+"\n ")
            powfile.flush()

            rfit = 10.0
            power, factor = medianpowerlaw(snap,hydro,10000000,centre,rinner=5.0,router=20.0,rfit=rfit)
            powfile.write(simname+" "+str(outnum)+"M5OUTER w, n0, rfit:"+str(power)+" "+str(factor)+" "+str(rfit)+"\n ")
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
        if hydro == "nH":
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

def findstarpos(simname,time):
    sim = hamusims[simname]
    tcode = time * Myrins / unit_t
    snap = sim.FindAtTime(tcode)
    stellar = stellars.FindStellar(snap)
    imax = np.where(stellar.mass == np.max(stellar.mass))[0][0]
    sinkid = stellar.sinkid[imax]-1
    sink = sinks.FindSinks(snap)
    boxlen = snap.RawData().info["boxlen"]
    pos = np.array([sink.x[sinkid],sink.y[sinkid],sink.z[sinkid]])/boxlen
    #import pdb; pdb.set_trace()
    return pos

if __name__ == "__main__":
    # Make figure
    powfile = open("../plots/powerlaws.txt","w")
    sims = ["IMF1_01","MASS_01"]
    labels = {"IMF1_01":"$10^4$ M$_{\odot}$ cloud",
              "MASS_01":"$10^5$ M$_{\odot}$ cloud"}
    #sims = ["MASS_04"]
    for simname in sims:
        if "IMF1" in simname:
            time = 3.38 # Time of first sink formation
        else:
            time = 1.2
        starpos = findstarpos(simname,time)
        #import pdb; pdb.set_trace()
        plotgradient(simname,"nH",time,starpos,labels[simname],
                     xlims=[0.03,25.0],powfile=powfile,suffix="starpos")
        #plotgradient(simname,"nH",time,np.array([0.5,0.5,0.5]),"Geometric Centre",
        #             xlims=[0.03,25.0],powfile=powfile,suffix="centre")
    powfile.close()
