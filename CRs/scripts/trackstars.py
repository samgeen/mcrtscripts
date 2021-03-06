'''
Tracks the stars in a simulation via their sinks over time
Sam Geen, February 2018
'''

from startup import *

import glob

def _sinksinsnap(snap):
    stellar = stellars.FindStellar(snap)
    return stellar.sinkid
sinksinsnap = Hamu.Algorithm(_sinksinsnap)

def allsinks(sim):
    sinkids = set()
    # Find unique set of sinks with stars in all simulations
    for snap in sim.Snapshots():
        sinkidlist = sinksinsnap(snap)
        for sinkid in sinkidlist:
            sinkids.add(sinkid)
    return list(sinkids)

def maketrack_intrasnap(sim,sinkids):
    # Make empty list
    sinkpostime = {sinkid: list() for sinkid in sinkids}
    sinkfiles = []
    for snap in sim.Snapshots():
        ro = snap.RawData()
        folder = ro.output_repos+"output_"+str(ro.iout).zfill(5)+"/"
        sinkfiles += glob.glob(folder+"sink_?????.dat?????")
    for sinkfile in sinkfiles:
        if os.stat(sinkfile).st_size > 0:
            good = True
            try:
                dat = np.loadtxt(sinkfile,converters = {0: lambda s: float(s)},
                                 delimiter=",")
                good = False
            except:
                print "Error in file", sinkfile, ", skipping..."
            if good:
                ids = dat[:,0]
                ts = dat[:,1]
                xs = dat[:,3]
                ys = dat[:,4]
                zs = dat[:,5]
                for i,t,x,y,z in zip(ids,ts,xs,ys,zs):
                    sinkpostime[int(i)] = np.array([x,y,z,t])
            
def _maketrack_snapsonly(snap,sinkid):
    sink = sinks.FindSinks(snap)
    pos = None
    if len(sink.id) > 0:
        idloc = np.where(sink.id == sinkid)[0]
        if len(idloc) > 0:
            idloc = idloc[0]
            pos = (sink.x[idloc],sink.y[idloc],sink.z[idloc])
    return pos
maketrack_snapsonly = Hamu.Algorithm(_maketrack_snapsonly)

def runforsim(simname):
    sim = hamusims[simname]
    # Find sink IDS to track
    sinkids = allsinks(sim)
    # Track them
    track = {}
    for sinkid in sinkids:
        track[sinkid] = []
    for snap in sim.Snapshots():
        for sinkid in sinkids:
            sinkpos = maketrack_snapsonly(snap,sinkid)
            if sinkpos is not None:
                track[sinkid].append(sinkpos)
    for sinkid in sinkids:
        track[sinkid] = np.array(track[sinkid])
    return track
    
if __name__=="__main__":
    runforsim("IMF01")
