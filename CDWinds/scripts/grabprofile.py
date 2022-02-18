'''
Get a specific profile from a simulation for more detailed investigations
Sam Geen, December 2021
'''

from startup import *

import rayprof

from collections import OrderedDict

nprofs = rayprof.nprofs
nray = rayprof.nray

def findstarpos(snap,useMyr=False):
    #sim = hamusims[simname]
    #if useMyr:
    #    tcode = time * Myrins / unit_t
    #else:
    #    tcode = time
    #snap = sim.FindAtTime(tcode)
    stellar = stellars.FindStellar(snap)
    if len(stellar.mass) == 0:
        print("No stellar objects in snapshot")
        return None
    imax = np.where(stellar.mass == np.max(stellar.mass))[0][0]
    sinkid = stellar.sinkid[imax]-1
    sink = sinks.FindSinks(snap)
    boxlen = snap.RawData().info["boxlen"]
    pos = np.array([sink.x[sinkid],sink.y[sinkid],sink.z[sinkid]])/boxlen
    return pos


def getgradients(simname,time=None,label=None,xlims=None,powfile=None,maxradius=None):
    '''
    Like prof but with a higher gradient
    '''
    global nray, nprofs,rays,rads
    sim = Hamu.Simulation(simname)
    # HACK
    print("Running for sim", simname)
    snaps = sim.Snapshots()
    if time is not None:
        tcode = time * Myrins / unit_t
        snaps = [sim.FindAtTime(tcode)]
    # Find ray where shock is most extended
    for snap in [snaps[-1]]:
        starpos = findstarpos(snap)
        if starpos is None:
            continue # skip this iteration
        centre = starpos
        r,nHs   = rayprof.makeprofsHamu(snap,"nH",centre)
        r,Ts    = rayprof.makeprofsHamu(snap,"T",centre)
        r,xHIIs = rayprof.makeprofsHamu(snap,"xHII",centre)
        radii = r[0:nray]
        nHs = np.reshape(nHs,(nprofs, nray)).T
        Ts = np.reshape(Ts,(nprofs, nray)).T
        xHIIs = np.reshape(xHIIs,(nprofs, nray)).T
        iraymax = None
        maxr = 0.0
        Tthresh = 1e6
        for iray in range(0,nray):
            rshock = radii[np.where(Ts[:,iray] > Tthresh)].max()
            if rshock > maxr:
                maxr = rshock
                iraymax = iray
    # Plot ray where shock is most extended at the earliest time there's a star
    print(maxr, iraymax)
    plt.clf()
    linestyles = ["-","--",":"]
    iline = -1
    nHlabel = "nH"
    Tlabel = "T"
    Plabel = "P"
    lnHs = []
    timelabels = []
    times = [0.1,0.2,0.3]
    snaps = [timefuncs.findsnapshot(sim,(t,"MyrFirstStar")) for t in times]
    for snap in snaps[::-1]:
        iline += 1
        line = linestyles[iline]
        starpos = findstarpos(snap)
        if starpos is None:
            continue # skip this iteration
        centre = starpos
        r,nHs   = rayprof.makeprofsHamu(snap,"nH",centre)
        r,Ts    = rayprof.makeprofsHamu(snap,"T",centre)
        r,xHIIs = rayprof.makeprofsHamu(snap,"xHII",centre)
        r,Ps = rayprof.makeprofsHamu(snap,"P",centre)
        radii = r[0:nray]
        nHs = np.reshape(nHs,(nprofs, nray)).T
        Ts = np.reshape(Ts,(nprofs, nray)).T
        xHIIs = np.reshape(xHIIs,(nprofs, nray)).T
        Ps = np.reshape(Ps,(nprofs, nray)).T
        lnH, = plt.plot(radii,nHs[:,iraymax],label=nHlabel,linestyle=line,color="k")
        plt.plot(radii,Ts[:,iraymax],label=Tlabel,linestyle=line,color="r")
        plt.plot(radii,Ps[:,iraymax],label=Plabel,linestyle=line,color="b")
        lnHs.append(lnH)
        tMyr = snaptime.Myr(snap) - timefuncs.FindTcreatedFirstStar(sim)
        timelabels.append("{:.1f}".format(tMyr)+" Myr")
        nHlabel = None
        Tlabel = None
        Plabel = None
        #plt.plot(radii,xHIIs[:,iraymax],label="xHII")
    legend1 = plt.legend(fontsize="small",ncol=3,frameon=False,loc="upper right")
    legend2 = plt.legend(lnHs, timelabels,fontsize="small",frameon=False,loc=(0.1,0.25))
    plt.gca().add_artist(legend1)
    if maxradius is not None:
        plt.xlim([0.0,maxradius])
    plt.xlabel("Radius / pc")
    plt.ylabel("Values at ray which has largest extent of wind")
    suffix = ""
    plt.xscale("linear")
    plt.yscale("log")
    #plt.ylim([1e-14,1e10])
    filename = "../plots/gradrays/"+simname+"/maxrprofiles.pdf"
    print("Writing "+filename)
    plt.savefig(filename,rasterized=True,dpi=200)

if __name__=="__main__":
    labels = OrderedDict()
    labels["SEED1_35MSUN_CDMASK_WINDUV"] = "Seed1, CDMask, Refine"
    labels["SEED0_35MSUN_CDMASK_WINDUV"] = "Seed0, CDMask, Refine"
    labels["SEED1_35MSUN_CDMASK_WINDUV"] = "Seed1, CDMask, Refine"
    labels["SEED2_35MSUN_CDMASK_WINDUV"] = "Seed2, CDMask, Refine"
    labels["SEED3_35MSUN_CDMASK_WINDUV"] = "Seed3, CDMask, Refine"
    labels["SEED4_35MSUN_CDMASK_WINDUV"] = "Seed4, CDMask, Refine"
    sims = labels.keys()
    #sims = ["MASS_04"]
    maxradius = 15.0 # in pc
    for simname in sims:
        getgradients(simname,maxradius=maxradius)
