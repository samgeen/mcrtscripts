"""
Make extinction plots
Sam Geen, May 2020
"""

from startup import *
import rayprof


def extinctionbymass(simname):
    sim = hamusims[simname]
    masstracks = {}
    # Hack to get this to run on already processed data
    for snap in sim.Snapshots()[0:17]:
        c,m = rayprof.FindColumnsInSnap(snap)
        for columns, mass in zip(c,m):
            if mass not in masstracks:
                masstracks[mass] = []
            masstracks[mass].append(columns)
    for mass, track in masstracks.iteritems():
        masstracks[mass] = np.array(track)
        print masstracks[mass].shape

def extinctionforeachsnap(simname):
    sim = hamusims[simname]
    for snap in sim.Snapshots()[27:38]:
        extinctioninsnap(snap, simname)
        
def extinctioninsnap(snap, simname):
    stellar = stellars.FindStellar(snap)
    if len(stellar.mass) == 0:
        return
    plt.clf()
    # Load columns of each sink
    cols,isinks = rayprof.FindColumnsBySinkInSnap(snap)
    if len(isinks) < 0:
        return
    sinktocolumns = {}
    print isinks
    for col in cols:
        print "MIN COLS", col[0:10]
    for columns, isink in zip(cols,isinks):
        sinktocolumns[isink] = columns
    # Plot stellar objects
    for istellar in range(0,len(stellar.mass)):
        sinkid = stellar.sinkid[istellar]
        mass = stellar.mass[istellar]
        columns = sinktocolumns[sinkid]
        cols = np.sort(columns)
        exts = NHtoAv(cols)
        marr = exts*0.0 + mass
        # Make weights that peak at 1
        weights = np.concatenate([np.linspace(0.1,1,len(exts)//2),
                                  np.linspace(1,0.1,len(exts)//2)])
        plt.scatter(marr,exts,c=weights,alpha=0.1,edgecolors='none')
    plt.xlabel("Stellar Mass / Msun")
    plt.ylabel("V-band extinction")
    plt.yscale("log")
    figname = "../plots/extinctions/extinctioninsnap_"+simname+"_"+str(snap.OutputNumber()).zfill(5)+".png"
    print "Plotting", figname
    plt.savefig(figname)
                
        
if __name__=="__main__":
    extinctionforeachsnap("18_LEGO")
