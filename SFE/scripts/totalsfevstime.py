'''
Plot of total SFE (total stellar mass vs total mass)
'''

from startup import *

import columndensity, images, sfeobs, sinks, starsvsdensegas

from pymses.utils import constants as C

class SFEPlotter(starsvsdensegas.SFEPlotter):
    def Plot(self,xtff):
        # TODO: Better interface to base class member variables                 
        def yfunc(data):
            mtot = data.gas+data.stars
            if np.median(self._Aklow) == 0.0:
                mtot = 1e4
            return data.stars / mtot
        def xfunc(data):
            return data.times
        self.PlotError(xfunc,yfunc,xtff)

def MakePlot(sims,Aklow=[0.0],allstars=True,ysoage=[0.0],ax=None):
    print "PLOTTING TOTAL SFE VS TIME FOR AKLOW=",Aklow,"ALLSTARS=",allstars
    fig = None
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    # Simulations
    for sim in sims:
        plotter = SFEPlotter(sim,Aklow,allstars,ysoage,__name__,
                             axis = ax, figure = fig)
        plotter.Plot(xtff=True)
    # Lada+ 2010 fit (by eye, since they don't give an offset)
    #tlada = np.array([0.0,20.0]) 
    #sfelada = np.array([0.095,0.095])
    #plt.plot(tlada,sfelada,color=r"#888888",linestyle="--")
    # Do legends
    plotter.Legend(clouds="lower right",rt=False)
    # Set up figure and output
    mgas = "$M_{\mathrm{gas}}(>A_{k} = "+str(Aklow)+")$"
    mstar = "$M_{*}$"
    ysoagetxt = "2 Myr"
    if fig is not None:
        ax.set_xlabel("Time - $_t{ff}$ / Myr")
    ax.set_ylabel("MSTAR / (MSTAR + MGAS)")
    if Aklow[0] == 0.0:
        ax.set_ylabel("TSFE = "+mstar + "/ $10^4$ "+Msolar)
        ax.set_ylim([1e-3,1.0])
    if Aklow[0] == 1e-6:
        ax.set_ylabel(mstar + "/ "+mgas)
    #if allstars:
    #    plt.ylabel("$M_{*}$ (All stars) / "+Msolar)
    ax.set_yscale("log")
    ax.set_xlim([0,20])
    shade4Myr.run(ax)
    plotter.Save()
