'''
Plot the SFE in each simulation using simulation and "observation" estimates
Sam Geen, June 2016
Sam Geen, January 2020
'''

from startup import *

import sfeobs, sfesim, linestyles, rdmfile

from pymses.utils import constants as C

black = r"#000000"
red   = r"#ff0000"

TWOMASSFWHM = 1.0/6.0 # Degrees

Avhigh = 5.0

def _drawlines(sim,colour=black,multiangle=True,rdm=None):
    esims = []
    times = []
    utime = None
    effs = np.zeros((6))
    lsnap = len(sim.Snapshots())
    eobs = np.zeros((lsnap,6))
    itime = 0
    for snap in sim.Snapshots():
        # Fill the arrays
        if utime is None:
            utime = snap.RawData().info["unit_time"].express(C.Myr)
        esims.append(sfesim.SFEsimHamu(snap)*100.)
        # Do the observational estimates and sort in size
        effs[0] = sfeobs.SFEobsHamu(snap,'x',
                                    extincthigh=True,extinctlimit=Avhigh)
        if multiangle:
            effs[1] = sfeobs.SFEobsHamu(snap,'x',frontToBack=False,
                                        extincthigh=True,extinctlimit=Avhigh)
            effs[2] = sfeobs.SFEobsHamu(snap,'y',
                                        extincthigh=True,extinctlimit=Avhigh)
            effs[3] = sfeobs.SFEobsHamu(snap,'y',frontToBack=False,
                                        extincthigh=True,extinctlimit=Avhigh)
            effs[4] = sfeobs.SFEobsHamu(snap,'z',
                                        extincthigh=True,extinctlimit=Avhigh)
            effs[5] = sfeobs.SFEobsHamu(snap,'z',frontToBack=False,
                                        extincthigh=True,extinctlimit=Avhigh)
            effs = np.sort(effs)
        eobs[itime,:] = effs*100.0 # as a %
        times.append(snap.Time()*utime)
        itime += 1
    times = np.array(times)
    if rdm is not None:
        rdm.AddPoints(times,esims,label=sim.Name()+" "+"sim")
    plt.plot(times,esims,color=colour,linestyle="--")
    if not multiangle:
        plt.plot(times,eobs[:,0],color="r",linestyle="-")
        if rdm is not None:
            rdm.AddPoints(times,eobs[:,0],label=sim.Name()+" "+"obs")
    else:
        taba = np.concatenate((times,times[::-1]))
        outer  = np.concatenate((eobs[:,0],eobs[:,5][::-1]))
        middle = np.concatenate((eobs[:,1],eobs[:,4][::-1]))
        inner  = np.concatenate((eobs[:,2],eobs[:,3][::-1]))
        #import pdb; pdb.set_trace()
        plt.fill(taba,inner ,alpha=0.33,edgecolor='none',facecolor=colour)
        plt.fill(taba,middle,alpha=0.33,edgecolor='none',facecolor=colour)
        plt.fill(taba,outer ,alpha=0.33,edgecolor='none',facecolor=colour)

def makeplot(sims,dolegend=True,plotlabel=""):
    # Get no-RT sim
    #nortname = sim.Name().replace("UV","UVWIND")
    #nortsim = Hamu.Simulation(nortname)
    # Plot
    plt.clf()
    blank = "55"
    rdm = rdmfile.RDMFile(__file__)
    for sim in sims:
        _drawlines(sim,black,multiangle=False,rdm=rdm)
    #_drawlines(sim,red)
    # Make cloud size text
    ax = plt.gca()
    #labelsize = "10$^4$ M$_{\odot}$"
    #ax.text(0.05, 0.95, labelsize, 
    #        transform=ax.transAxes, fontsize="medium",
    #        verticalalignment='top')
    if dolegend:
        # Make estimate type legend
        lines1 = [mlines.Line2D([],[],color='k',linestyle='--'),
                  mlines.Line2D([],[],color='r',linestyle='-')]
        labels1 = [r"$\epsilon_{sim}$",
                   r"$\epsilon_{obs}$",]
        leg1 = plt.legend(lines1,labels1,fontsize="small",
                          frameon=False,loc="upper left")
        # Make RT/NRT legend
        #lines2 = [mlines.Line2D([],[],color=red,linestyle='-'),
        #          mlines.Line2D([],[],color=black,linestyle='-')]
        #labels2 = ["With UV radiation + Winds",
        #           "With UV radiation"]
        #leg2 = plt.legend(lines2,labels2,fontsize="small",
        #                  frameon=False,loc="center left")
        plt.gca().add_artist(leg1)
    plt.xlabel("Time / Myr")
    plt.ylabel("\% SFE")
    folder = "../plots/sfe/"
    limtxt = "_Avlim"+str(int(Avhigh))
    figname = "sfe_"+plotlabel+limtxt+".pdf"
    MakeDirs(folder)
    plt.savefig(folder+figname)
    rdm.Write(folder+figname)
    

if __name__=="__main__":
    makeplot([hamusims[sim] for sim in icsims],plotlabel="IC")
    makeplot([hamusims[sim] for sim in imfsims],plotlabel="IMF")
