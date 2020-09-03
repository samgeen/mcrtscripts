"""
Make extinction plots
Sam Geen, September 2020
"""

from startup import *
import extinctionplots, snaptime

from pymses.utils import constants as C  

from scipy import interpolate

AVtext = "A$_{\mathrm{V}}$"

def plotagevisiblevsstellarmass(simname,extinctionlimit=1):
    '''
    Plot the stellar age the star is visible versus its mass
    '''
    sim = hamusims[simname]
    simsnaps = sim.Snapshots()
    simtimes = [snaptime.Myr(snap) for snap in simsnaps]
    # Get the stellar object tracks
    tracks, isnaps = extinctionbyparticle(simname,wholesink=False)
    masses = tracks.keys()
    # Get ages each stellar object is visible at
    visibleages = []
    # TODO:
    # - Get birth time of each star
    # - Turn each track into an interp1d function of stellar age vs extinction
    # - Find age where extinction == extinctionlimit
    # - Plot ages vs masses
    for mass in masses:
        extinctions = tracks[mass]
        snapnums = isnaps[mass]
        firstsnap = simsnaps[snapnums[0]]
        stellar = stellars.FindStellar(firstsnaps)
        whichstar = np.where(stellar.mass == mass)
        tcreated = stellar.tcreated[whichstar]
        lifetime = stellar.lifetime[whichstar]
        ages = simtimes[snapnums] - tcreated
        # Interpolate extinction track
        agefunc = interpolate.interp1d(extinctions, ages)
        visibleage = agefunc[extinctionlimit]
        visibleages.append(visibleage)
    
    plt.clf()
    plt.scatter(masses, visibleages)
    limtxt = "extinctionlimit_"+str(extinctionlimit)
    plt.savefig("../plots/visibleages_"+limtxt+"_"+simname+".pdf")

if __name__=="__main__":
   plotagevisiblevsstellarmass("128_LEGO") 

