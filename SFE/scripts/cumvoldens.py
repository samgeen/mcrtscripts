'''
Cumulative volume density
Sam Geen, June 2016
'''

from startup import *
from pymses.filters import CellsToPoints
from pymses.utils import constants as C
import freefalltime, linestyles

def rebinarrays(x,y,nbin):
    bins = np.linspace(np.log10(0.999*x.min()),np.log10(1.001*x.max()),nbin)
    binr = bins.max() - bins.min()
    bini = ((np.log10(x)-bins.min())/binr*nbin).astype(int)
    counts = bins*0.0
    hist = bins*0.0
    counts[bini] += 1.0
    hist[bini] += y
    mask = counts > 0
    counts[counts==0] = 1.0
    hist /= counts
    bins = 10.0**bins
    hist = hist[mask]
    bins = bins[mask]
    return bins, hist

def _CumMassvsVolDens(snap,nbin=0):
    print "Finding mass vs density in snap", snap.iout
    # Extract raw data
    amr = snap.amr_source(["rho"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    vols = (cells.get_sizes())**3.0
    dens = cells["rho"]
    mass = dens*vols*snap.info["unit_mass"].express(C.Msun)
    # Sort arrays
    # (in reverse to cumulatively sum from smallest to largest density)
    mask = np.argsort(dens)[::-1]
    dens = dens[mask]
    mass = mass[mask]
    # Make cumulative mass
    mass = np.cumsum(mass)
    dens = dens[::-1]
    mass = mass[::-1]
    # If nbin > 0, rebin arrays
    if nbin > 0:
        dens, mass = rebinarrays(dens,mass,nbin)
    return dens, mass
CumMassvsVolDens = Hamu.Algorithm(_CumMassvsVolDens)

def PlotForSim(ax,sim,tfffact=1.0):
    print "Running for", sim.Name()
    # Load at tff
    tff = freefalltime.Tff(sim)
    # Find data to plot
    snap = sim.FindAtTime(tff*tfffact)
    # Trick Hamu into redoing stuff by slightly changing nbin
    dens, cmass = CumMassvsVolDens(snap,nbin=201)
    # Make plot lines
    col = linestyles.colour(sim.Name())
    line = '-'
    if tfffact > 1.0:
        line = '--'
    if tfffact > 2.0:
        line = '-.'
    ax.plot(dens,cmass,color=col,linestyle=line)

def PlotSims(ax,tffs):
    simnames = [s+"-RT" for s in ["L","M","S","XS"]]
    # 5 arcmin (PLANCKRES) for d = 150pc = 0.22pc
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        for tff in tffs:
            PlotForSim(ax,sim,tfffact=float(tff))
    txt = "Full Res"
    ax.text(0.05, 0.95,txt,
            horizontalalignment='left',
            verticalalignment='top',
            transform = ax.transAxes)

def MakePlot(ax=None,tffs=[1,2]):
    saveplot = True
    fig = None
    if ax is None:
        fig, ax = plt.subplots(1)
    else:
        saveplot = False
    PlotSims(ax,tffs=tffs)
    if saveplot:
        ax.set_ylabel("Cumulative Mass / "+Msolar)
    ax.set_xlabel("$n_{H}$ / cm$^{-3}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    # Do legend - simulations
    lines, labels = linestyles.sizelegend()
    leg1 = ax.legend(lines,labels,fontsize="x-small",frameon=False,loc="lower left")
    # Do legend - times
    tff = r"$t_{\mathrm{ff}}$"
    lines = [mlines.Line2D([], [], color="k",linestyle="-", label=tff),
             mlines.Line2D([], [], color="k",linestyle="--",label="2 "+tff)]
    if 4 in tffs:
        lines.append(mlines.Line2D([], [], color="k",
                                   linestyle="-.",label="4 "+tff))
    labels = [line.get_label() for line in lines]
    leg2 = ax.legend(lines,labels,fontsize="small",frameon=False,loc="upper right")
    ax.add_artist(leg1)
    ax.add_artist(leg2)
    folder = "../plots/planck/"
    MakeDirs(folder)
    outfile = "massvsvoldens.pdf"
    ax.set_xlim([3.0,2e9])
    ax.set_ylim([1,3e5])
    if saveplot:
        print "Saving figure to",folder+outfile
        fig.savefig(folder+outfile)

if __name__=="__main__":
    Hamu.Workspace("HIISFE")
    MakePlot()
