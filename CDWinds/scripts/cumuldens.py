'''
Cumulative density plot
Sam Geen, November 2018
'''

from startup import *
import columndensity, linestyles
import astropy.convolution as C
import astropy.modeling.models as M
import scipy.ndimage

documulative=True

# Get Sigma in Msun / pc^2 from NH in cm^-2
def NHtoSigma(NH):
    return NH * mH/X / (Msuning / pcincm**2)

def cumtohist(x,y,nbin=20):
    bins = np.linspace(np.log10(x.min()),np.log10(x.max()),nbin)
    binr = bins.max() - bins.min()
    newy = -np.diff(y) / np.diff(x)
    newx = np.log10(0.5*(x[1:]+x[:-1]))
    bini = ((newx-bins.min())/binr*nbin).astype(int)
    counts = bins*0.0
    hist = bins*0.0
    counts[bini] += 1.0
    hist[bini] += newy
    counts[counts==0] = 1.0 # prevent nans
    hist /= counts
    bins = 10.0**bins
    return bins, hist
            

def resampleplot(ax,x,y,color="k",linestyle="-",nresample=1000):
    '''
    Resample a line to plot to save plot size
    '''
    # HACK - eh, figure this out later
    '''
    if len(x) != len(y):
        print "Error: arrays different sizes in resampleplot"
        raise OverflowError
    if len(x) > nresample:
        downfact = int(len(x) / nresample)
        x = scipy.ndimage.zoom(x,1.0/downfact)
        y = scipy.ndimage.zoom(y,1.0/downfact)
    '''
    # Plot diff instead?
    mask = np.isfinite(x)*np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) == 0:
        return
    if not documulative:
        x, y = cumtohist(x,y)
    y2 = y[x > 100][-1]
    y3 = y[x > 1000][-1]
    print "Value at Sigma=100  Msun/pc^2:", y2
    print "Value at Sigma=1000 Msun/pc^2:", y3
    #import pdb; pdb.set_trace()
    ax.plot(x,y,color=color,linestyle=linestyle,zorder=2)

def PlotForSim(ax,sim,timeInMyr):
    print "Running for", sim.Name()
    # Load at time in Myr
    snap = sim.FindAtTime(timeInMyr * Myrins / unit_t)
    # Make plot lines
    col = linestyles.Colour(sim.Name())
    Avlow = 0.1
    def PlotForLOS(snap,los,line):
        NH, cmass = columndensity.CumulativeMassHamu(snap,los,"NH")
        Sigma = NHtoSigma(NH)
        resampleplot(ax,Sigma,cmass,color=col,linestyle = line)
        print "SIMNAME", sim.Name()
        if "MASS" in sim.Name():
            mlog = 5
        else:
            mlog = 4
        simtxt = r"$10^"+str(mlog)+"$ "+Msolar+" cloud"
        ax.text(Sigma[0], cmass[0], simtxt, fontsize="x-small",
                horizontalalignment='center',
                verticalalignment='top',
                backgroundcolor="white",
                zorder=1)
    line = '-'
    PlotForLOS(snap,'x',line)

def PlotSims(ax):
    simnames = ["IMF1_01","MASS_01"]
    times = [3.38,1.2]
    for simname, time in zip(simnames, times):
        sim = hamusims[simname]
        PlotForSim(ax,sim,time)

def MakePlot(ax=None):
    saveplot = True
    fig = None
    if ax is None:
        fig, ax = plt.subplots(1)
    else:
        saveplot = False
    PlotSims(ax)
    NHtxt = "$N_{\mathrm{H}}$"
    Sigmatxt = "$\Sigma$"
    if saveplot:
        #ax.set_ylabel("M($<$ "+NHtxt+" / "+Msolar)
        ax.set_ylabel("$M$($<$ "+Sigmatxt+") / "+Msolar)
    #ax.set_xlabel(NHtxt+" / cm$^{-2}$")
    ax.set_xlabel(Sigmatxt+" / "+Msolar+" pc$^{-2}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    # Do legend - simulations
    #if saveplot:
    #    lines, labels = linestyles.sizelegend()
    #    leg1 = ax.legend(lines,labels,fontsize="small",frameon=False,loc="lower left")
    # Do legend - times
    #tff = r"$t_{\mathrm{ff}}$"
    #lines = [mlines.Line2D([], [], color="k",linestyle="-", label=tff),
    #         mlines.Line2D([], [], color="k",linestyle="--",label="2 "+tff),
    #         mlines.Line2D([], [], color="k",linestyle=":", label="Gould Belt")]
    #if not PSFon:
    #    lines = lines[0:2]
    #if 4 in tffs:
    #    lines.append(mlines.Line2D([], [], color="k",
    #                               linestyle="-.",label="4 "+tff))
    #labels = [line.get_label() for line in lines]
    #leg2 = ax.legend(lines,labels,fontsize="small",frameon=False,loc="upper right")
    #if saveplot:
    #    ax.add_artist(leg1)
    #ax.add_artist(leg2)
    #folder = "../plots/planck/"
    #MakeDirs(folder)
    outfile = "../plots/cumuldens.pdf"
    xmin = NHtoSigma(7e20)
    xmax = NHtoSigma(5e24)
    ax.set_xlim(xmin, xmax)
    if documulative:
        # Do analytic gradient fit
        # n0 = 432 cm^-3, r0 = 1 pc
        # M = (2 * pi * rho0 * r0^2)^2 / Sigma
        anaconst = (2.0 * np.pi * 432.0*mH/X * (pcincm)**(3) / Msuning)**2.0
        xana = NHtoSigma(np.logspace(22,23.2))
        yana = anaconst/xana
        #yana *= 3e4*xana[0] # Normalise so the start is at 1e3
        ax.plot(xana, yana,"k--")
        anatxt1 = r"$M \propto \Sigma^{-1}$"
        anatxt2 = "Isothermal fit around first star"
        ax.text(xana[0], yana[0], anatxt1, fontsize="x-small",
                horizontalalignment='center',
                verticalalignment='bottom')
        ax.text(xana[-1], yana[-1], anatxt2, fontsize="x-small",
                horizontalalignment='center',
                verticalalignment='top')
        # Set y axis if cumulative and save
        ax.set_ylim([1,3e5])
    if saveplot:
        print "Saving figure to",outfile
        fig.savefig(outfile)

if __name__=="__main__":
    documulative=True # False is a mess so far
    MakePlot()
