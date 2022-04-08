'''
Plot the solid angle subtended by the wind bubble at each radius
Sam Geen, April 2022
'''

from startup import *

import rayprof

from collections import OrderedDict

nprofs = rayprof.nprofs
nray = rayprof.nray

def findstarpos(snap,useMyr=False):
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


def plotsolidangles(simname,times,label=None,xlims=None,powfile=None,maxradius=None):
    '''
    Like prof but with a higher gradient
    '''
    global nray, nprofs,rays,rads
    #sim = Hamu.Simulation(simname)
    sim = hamusims[simname]
    # HACK
    print("Running for sim", simname)
    snaps = []
    for timetuple in times:
        tcode = time * Myrins / unit_t
        snaps.append(timefuncs.findsnapshot(sim,timetuple)
    plt.clf()
    lines = ["-","--",":"]
    iline = -1
    # Find ray where shock is most extended
    for timetuple, snap in zip(times,snaps):
        iline += 1
        line = lines[iline]
        time, tunit = timetuple
        str = str(time)
        if "Myr" in tunit:
            tstr = str(time)+" Myr"
        starpos = findstarpos(snap)
        if starpos is None:
            continue # skip this iteration
        centre = starpos
        #r,nHs   = rayprof.makeprofsHamu(snap,"nH",centre)
        r,Ts    = rayprof.makeprofsHamu(snap,"T",centre)
        #r,xHIIs = rayprof.makeprofsHamu(snap,"xHII",centre)
        radii = r[0:nray]
        #nHs = np.reshape(nHs,(nprofs, nray)).T
        Ts = np.reshape(Ts,(nprofs, nray)).T
        #xHIIs = np.reshape(xHIIs,(nprofs, nray)).T
        iraymax = None
        maxr = 0.0
        Tthresh = 1e6
        windinray = radii*0.0
        for iray in range(0,nray):
            windinray[np.where(Ts[:,iray] > Tthresh)] += 1.0
        windinray *= 4.0 * np.pi / nray
        mask = np.where(windinray > 0.0)
        radii = radii[mask]
        windinray = windinray[mask]
        plt.plot(radii,windinray,linestyle=line,label=tstr)
    plt.xlabel("Radius / pc")
    plt.ylabel("$\Omega_{\mathrm{wind}}(R)$")
    suffix = ""
    plt.xscale("linear")
    plt.yscale("log")
    #plt.ylim([1e-14,1e10])
    filename = "../plots/windsolidangle_"+simname+".pdf"
    print("Writing "+filename)
    plt.savefig(filename)
    

if __name__=="__main__":
    timetuples = [(t,"MyrFirstStar") for t in (0.1,0.2,0.3)]
    simname = "SEED1_35MSUN_CDMASK_WINDUV"
    for simname in sims:
        plotsolidangles(simname)
