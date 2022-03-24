'''
Plot the ram pressure in different simulations
Sam Geen, March 2022
'''

from startup import *

import rayprof

from collections import OrderedDict

nprofs = rayprof.nprofs
nray = rayprof.nray

def plotallprofiles(simnames,hydro,centre,label=None,xlims=None,powfile=None,suffix=None,plotmean=False,xscale="log"):
    '''
    Plot the average (median or mean) of a given profile
    '''
    global nray, nprofs,rays,rads
    for simname in simnames:
        # Setup
        sim = Hamu.Simulation(simname)
        print("Running for sim", simname, "hydro", hydro)
        mkdir("../plots/gradrays")
        mkdir("../plots/gradrays/"+simname)
        path = "../plots/gradrays/"
        # Find snapshot to use
        snap = None
        snaps = sim.Snapshots()
        for currsnap in snaps: 
            mkdir(path)
            # Process
            if "starpos" in suffix:
                starpos = rayprof.findstarpos(simname,currsnap.Time())
                if starpos is None:
                    continue # skip this iteration
                hydrofuncs.allhydros.AddGlobals({"starpos":starpos})
                centre = starpos
                snap = currsnap
        r,p = makeprofsHamu(snap,hydro,centre)
        radii = r[0:nray]
        profs = np.reshape(p,(nprofs, nray)).T # No idea
        grad = np.arange(0,51)*2.0
        percs = np.percentile(profs,grad,axis=1).T
        # Plot
        plt.clf()
        # Make median line
        if not plotmean:
            plt.plot(radii,percs[0:nprofs,25],"k")
        # Plot the mean line 
        else:
            rsph,psph = profilesphere.profileHamu(snap,hydro,10000000,centre,rcut=0.1)
            # Mask out parts of the mean line with no values
            mask = psph != 0
            rsph = rsph[mask]
            psph = psph[mask]
            plt.plot(rsph,psph,"k--")
        # Output diagnostic information if requested
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
        plt.yscale(hydrofuncs.yscale(hydro))
        if hydro == "nH":
            plt.ylim([3e-2,1e7])
        if hydro == "vrad":
            plt.ylim([-15,25])
        if hydro == "P":
            plt.ylim([1e-17,1e-6])
        plt.xlabel("Radius / pc")
        plt.ylabel(hydrofuncs.hydro_label(hydro))
        if xlims is not None:
            plt.xlim(xlims)
        plt.xscale(xscale)
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
        plt.grid()
        #plt.gca().set_rasterized(True)
    plt.savefig(path+"/simsetprofiles_"+hydro+"_"+suffix+"_x"+xscale+".png",dpi=200)

if __name__ == "__main__":
    # Use no feedback runs
    labels["SEED0_35MSUN_CDMASK_NOFB"] = "Seed0, No Feedback"
    labels["SEED1_35MSUN_CDMASK_NOFB"] = "Seed1, No Feedback"
    labels["SEED2_35MSUN_CDMASK_NOFB"] = "Seed2, No Feedback"
    labels["SEED3_35MSUN_CDMASK_NOFB"] = "Seed3, No Feedback"
    labels["SEED4_35MSUN_CDMASK_NOFB"] = "Seed4, No Feedback"

    simnames = labels.keys()
    time = None
    starpos = None
    plotallprofiles(simnames,"nH",starpos,labels[simname],
                    xlims=xlims=[0.03,25.0],powfile=None,suffix="starpos_nofb",plotmean=False,xscale="log")
    plotallprofiles(simnames,"Pram",starpos,labels[simname],
                    xlims=xlims=[0.03,25.0],powfile=None,suffix="starpos_nofb",plotmean=False,xscale="log")
