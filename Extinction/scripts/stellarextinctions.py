"""
Make extinction plots
Sam Geen, September 2020
"""

from startup import *
import extinctionplots, snaptime

from pymses.utils import constants as C  

from scipy import interpolate

AVtext = "A$_{\mathrm{V}}$"

def plotagevisiblevsstellarmass(simnames,extinctionlimit=1,scalewithlifetime=False,luminositylimit=None):
    '''
    Plot the stellar age the star is visible versus its mass
    '''
    plt.clf()
    for simname in simnames:
        sim = hamusims[simname]
        simsnaps = {snap.OutputNumber():snap for snap in sim.Snapshots()}
        simtimes = {snap.OutputNumber():snaptime.Myr(snap) for snap in simsnaps.values()}
        # Get the stellar object tracks
        tracks, isnaps = extinctionplots.extinctionbyparticle(simname,wholesink=False)
        masses = tracks.keys()
        # Get ages each stellar object is visible at
        visibleages = []
        visiblemasses = []
        # TODO:
        # - Get birth time of each star
        # - Turn each track into an interp1d function of stellar age vs extinction
        # - Find age where extinction == extinctionlimit
        # - Plot ages vs masses
        for mass in masses:
            # TODO: turn into probability plot based on los
            extinctions = np.array(tracks[mass][:,50])
            extinctions = NHtoAv(extinctions)
            snapnums = np.array(isnaps[mass])+1
            firstsnap = simsnaps[snapnums[0]]
            stellar = stellars.FindStellar(firstsnap)
            whichstar = np.where(stellar.mass == mass)
            tcreated = stellar.tcreated[whichstar][0]
            lifetime = stellar.lifetime[whichstar][0]
            ages = np.array([simtimes[num] for num in snapnums]) - tcreated
            # Get luminosities in the V band
            LVbands = np.array([singlestar.star_bandenergies(mass, ageins, 1.0) for ageins in ages*Myrins])
            # Get extincted LV values
            LVextincted = np.array([L * 10.0**(-extinction/2.5) for L, extinction in zip(LVbands, extinctions)])
            # Interpolate extinction track
            if len(ages) > 1:
                #agefunc = interpolate.interp1d(extinctions, ages)
                if extinctionlimit is not None:
                    if extinctionlimit > extinctions.max():
                        visibleage = 0.0
                    elif extinctionlimit < extinctions.min():
                        visibleage = lifetime
                    else:
                        visibleage = ages[np.where(extinctions < extinctionlimit)].min()
                if luminositylimit is not None:
                    if luminositylimit > LVextincted.max():
                        visibleage = lifetime
                    elif luminositylimit < LVextincted.min():
                        visibleage = 0.0
                    else:
                        visibleage = ages[np.where(LVextincted > luminositylimit)]
                if scalewithlifetime:
                    visibleage /= lifetime
                visibleages.append(visibleage)
                visiblemasses.append(mass)
        plt.scatter(visiblemasses, visibleages)
    plt.xlabel("Stellar Mass / Msun")
    if scalewithlifetime:
        plt.ylabel("Fraction of star's age until it is visible")
    else:
        plt.ylabel("Age star becomes visible / Myr")
    plt.xscale("log")
    plt.yscale("log")
    visibleages = np.array(visibleages)
    if scalewithlifetime:
        plt.ylim([visibleages[visibleages > 0.0].min(),1.0])
    else:
        plt.ylim([visibleages[visibleages > 0.0].min(),visibleages.max()])
        #if len(visibleages) > 0:
    #    plt.ylim([np.array(visibleages).min(),1.0])
    limtxt = "extinctionlimit_"+str(extinctionlimit)
    lifetxt = ""
    if scalewithlifetime:
        lifetxt = "_scaledwithlifetime"
    plt.savefig("../plots/visibleages_"+limtxt+lifetxt+"_"+simnames[0]+"_allsims.pdf",bbox_inches='tight')

if __name__=="__main__":
    luminositylimit = 10**4.5 * 2e33 # (limit / Lsolar) = 4.5 guestimate from Schootemeijer et al 2020
    for simnames in [imfsims, icsims]:
        for scalewithlifetime in [True, False]:
        plotagevisiblevsstellarmass(simnames,luminositylimit=luminositylimit,scalewithlifetime=scalewithlifetime)
        #plotagevisiblevsstellarmass(simnames,1,False) 
        #plotagevisiblevsstellarmass(simnames,1,True) 

