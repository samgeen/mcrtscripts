'''
Make profiles in rays positioned uniformly on a sphere
Sam Geen, November 2014
'''

import customplot
import isomesh, hydrofuncs
import pymses, Hamu, os
import numpy as np
import matplotlib.pyplot as plt
from pymses.utils import constants as C

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

def makerays(centre=[0.5,0.5,0.5],length=0.5):
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
    # Make rays from (0,0,0)
    arrscale = np.arange(0,1.0,dr)*length
    arrscale = np.tile(arrscale, nprofs)
    for i in range(0,3):
        rays[0:nprofs*nray,i] *= arrscale
    # Find radial positions
    rads = np.sqrt(np.sum(rays**2,1))
    # Add centre positions to ray
    for i in range(0,3):
        rays[0:nprofs*nray,i] += centre[i]

def makeprofs(snap, hydro,centre,length):
    '''
    Make profiles
    '''
    global rays, rads, nprofs, nray
    #if rays is None or rads is None:
    # Always run because centre or length can change
    # TODO: make one set of rays and just scale it every time
    makerays(centre,length)
    amr = hydrofuncs.amr_source(snap,hydro)
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

def medianprofile(snap,hydro,centre,length):
    global nray, nprofs,rays,rads
    #sim = Hamu.Simulation(simname)
    r,p = makeprofsHamu(snap.hamusnap,hydro,centre,length)
    radii = r[0:nray]
    profs = np.reshape(p,(nprofs, nray)).T # No idea
    percs = np.percentile(profs,[0,25,50,75,100],axis=1).T
    pmed = percs[0:nprofs,2]
    return radii, pmed

#medianprofileHamu = Hamu.Algorithm(medianprofile)
# Compatability hack for other scripts
#profileHamu = Hamu.Algorithm(medianprofile)
    
def medianradius(snap,centre,length,xlim=0.1):
    '''
    Find the median ionisation front radius
    '''
    r,p = makeprofsHamu(snap,"xHII",centre,length)
    profs = np.reshape(p,(nprofs, nray)).T # No idea
    rlims = np.zeros((nprofs))
    for iprof in range(0,nprofs):
        prof = profs[iprof,:]
        if prof.max() > xlim:
            rlims[iprof] = r[profs < xlim]
    return np.percentile(rlims,[50])[0]

def cutprofile(snap,hydro,centre,length,hcut):
    '''
    Cut any rays where the profile is above a given value
    '''
    global nray, nprofs,rays,rads
    #sim = Hamu.Simulation(simname)
    r,p = makeprofsHamu(snap.hamusnap,hydro,centre,length)
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

def plot(snap,hydro,centre=[0.5,0.5,0.5],length=0.5,suffix="",simname=""):
    global nray, nprofs,rays,rads
    # Setup
    mkdir("../plots/rays")
    mkdir("../plots/rays/"+simname)
    path = "../plots/rays/"+simname+"/"+hydro
    mkdir(path)
    # Process
    r,p = makeprofsHamu(snap,hydro,centre,length)
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
    plt.savefig(path+"/profile_"+suffix+"_"+str(snap.OutputNumber()).zfill(5)+".png")

def plotgradient(simname,hydro,centre,length,outnum=None,label=None):
    '''
    Like prof but with a higher gradient
    '''
    global nray, nprofs,rays,rads
    sim = Hamu.Simulation(simname)
    # HACK
    print "Running for sim", simname, "hydro", hydro
    snaps = sim.Snapshots()
    if not outnum is None:
        snaps = [sim.Snapshots()[outnum-1]]
    for snap in snaps: 
        # Setup
        mkdir("../plots/gradrays")
        mkdir("../plots/gradrays/"+simname)
        path = "../plots/gradrays/"+simname+"/"+hydro
        mkdir(path)
        # Process
        r,p = makeprofsHamu(snap,hydro,centre,length)
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
        #rsph,psph = profilesphere.profileHamu(snap,hydro,1e6)
        #plt.plot(rsph,psph,"k--")
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
            plt.ylim([1e-2,1e8])
        if hydro == "vrad":
            plt.ylim([-15,25])
        plt.xlabel("Radius / pc")
        plt.ylabel(hydrofuncs.hydro_label(hydro))
        # Add text label if one exists:
        if not label is None:
            xlims = plt.gca().get_xlim()
            ylims = plt.gca().get_ylim()
            xr = xlims[1] - xlims[0]
            yr = ylims[1] - ylims[0]
            plt.text(xlims[1]-0.05*xr,
                     ylims[1]*0.1,
                     label,
                     horizontalalignment='right',
                     verticalalignment='top')
        plt.savefig(path+"/profile_"+str(snap.OutputNumber()).zfill(5)+".png")

if __name__ == "__main__":
    # Make figure
    sims = ["N48_M4_B00","N48_M4_B02"]
    plotgradient("N48_M4_B00","vrad",5)
    plotgradient("N00_M4_B02","vrad",5)
    plotgradient("N48_M4_B00","rho",5,"No B-field")
    plotgradient("N00_M4_B02","rho",5,"B-field included")
    '''
    #sims = ["N00_M4_B02","N47_M4_B02","N48_M4_B00",
    #        "N48_M4_B02","N48_M4_B02_HR","N49_M4_B02"]
    #sims = ["N48_M4_B02_C","N48_M4_B02_C2","N48_M4_B02_F2","N48_M4_B02F3"]
    sims = ["N00_M4_B02",
            "N48_M4_B02",
            "N48_M4_B02_C",
            "N48_M4_B02_C2",
            "N48_M4_B02_F2",
            "N48_M4_B02_F3"]
    #hydros = ["rho","P","T","xHII","vrad"]
    #hydros = ["Pram","vrad","P","rho","T","xHII"]
    sims = ["N00_M4_B02","N48_M4_B00"]
    hydros = ["rho","spd","vrad","Pram"]
    for sim in sims:
        for hydro in hydros:
            plotgradient(sim,hydro)
    '''
