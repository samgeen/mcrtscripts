'''
Plot the distribution of SFE
Sam Geen, March 2018
'''

from startup import *

import plotproperties

def makeplot(ax,simnames):
    sfes = []
    for simname in simnames:
        sim = hamusims[simname]
        snap = sim.Snapshots()[-1]
        sfe = plotproperties.tsfeinsnap(snap)
        sfes.append(sfe)
    sfes.sort()
    sfes = np.array(sfes)
    numsims = len(sfes)
    N = np.arange(0,numsims)+1
    ax.plot(sfes,N,"k")
    # Plot IQR,median
    n25 = 4#0.25*numsims+0.5
    n50 = 7#0.50*numsims+0.5
    n75 = 10#0.75*numsims+0.5
    xlims = np.array([0.05,0.24])
    ax.set_xlim(xlims)
    ax.plot(xlims,xlims*0.0+n25,"k:",alpha=0.7)
    ax.plot(xlims,xlims*0.0+n50,"k--",alpha=0.7)
    ax.plot(xlims,xlims*0.0+n75,"k:",alpha=0.7)
    return sfes[n25-1], sfes[n50-1], sfes[n75-1]
    
def run(simnamesets,plotlabels):
    plt.clf()
    fig, axes = plt.subplots(1,2,sharex=True,sharey=True)
    first = True
    iqrs = {}
    for ax, simnames, plotlabel in zip(axes,simnamesets,plotlabels):
        iqrs[plotlabel] = makeplot(ax,simnames)
        ax.set_xlabel("SFE (final)")
        if first:
            ax.set_ylabel("N($<$ SFE)")
            first = False
        ax.text(0.95,0.95,plotlabel,ha="right",va="top",transform=ax.transAxes)
        #ax.legend(frameon=False,fontsize="x-small",loc="upperleft")
    fig.subplots_adjust(wspace=0)
    fig.set_size_inches(14,6)
    fig.savefig(plotfolder+"sfedistribution_both.pdf", dpi=80)  
    print "IQRs", iqrs

if __name__=="__main__":
    run((imfsims,icsims),(linestyles.starsc,linestyles.turbsc))
