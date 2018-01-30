'''
Plot of stellar mass / mass in gas A_k>0.8 / yso age vs time 
'''

from startup import *

import columndensity, images, sfeobs, sinks, starsvsdensegas
import freefalltime

from pymses.utils import constants as C

p=2

rfact = {"L":1.5**p,
         "M":1.0**p,
         "S":0.75**p,
         "X":0.5**p}

class SFEPlotter(starsvsdensegas.SFEPlotter):
    def Plot(self,usetff=False):
        # TODO!!!!!!!! CHECK WHY usetff WAS TRUE
        # TODO: Better interface to base class member variables
        def yfunc(data):
            tsf= 2.0 # Myr, from Lada et al
            #tff = freefalltime.Tff(self._sim,Myr=True)
            #if not usetff:
            tff = 1.0
            return data.stars / data.gas * tff # * rfact[self._sim.Name()[0]]# / tsf
        def xfunc(data):
            #tff = freefalltime.Tff(self._sim,Myr=True)
            #if not usetff:
            tff = 1.0
            return data.times / tff
        self.PlotError(xfunc,yfunc,usetff)

def MakePlot(sims,Aklow=[0.8],allstars=False,ysoage=[0.0],ax=None,
             usetff=True):
    print "PLOTTING SFR VS TIME FOR AKLOW=",Aklow,"ALLSTARS=",allstars
    fig = None
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    # Simulations
    for sim in sims:
        plotter = SFEPlotter(sim,Aklow,allstars,ysoage,__name__,
                             axis=ax,figure=fig)
        plotter.Plot(usetff=usetff)
    # Do legends
    if fig is not None:
        plotter.Legend()
    # Set up figure and output
    Aktxt = str(np.median(Aklow))
    mgas = "$M("+Aktxt+")$"
    mstar = "$M_{YSO}$"
    ysoagetxt = "2 Myr"
    tfftxt = "$t_{ff}$"
    ax.set_xlabel("Time / Myr")
    if usetff:
        ax.set_xlabel("Time / $t_{\mathrm{ff}}$")
    ylabel = "OSFE = "+mstar+" / "+mgas
    #if usetff:
    #    ylabel += r" $\times$ "+tfftxt
    ax.set_ylabel(ylabel)
    if allstars:
        ax.set_ylabel("$M_{*}$ (All stars) / "+Msolar)
    ax.set_yscale("log")
    if not usetff:
        ax.set_xlim([0,20])
        ax.set_ylim([1e-3,30.0])
    else:
        ax.set_xlim([-4.25,20-4.25])
        ax.set_ylim([1e-2,10.0])
    shade4Myr.run(ax)
    # Lada+ 2010 fit (by eye, since they don't give an offset)
    tlada = list(ax.get_xlim())
    sfelada = np.array([0.095,0.095])
    #if not usetff:
    ax.plot(tlada,sfelada,color=r"#888888",linestyle="--")

    plotter.Save()

if __name__=="__main__":
    simnames = ["L-RT","M-RT","S-RT","XS-RT",
            "L-NRT","M-NRT","S-NRT","XS-NRT"]
    sims = [Hamu.Simulation(s) for s in simnames]
    Aklows = np.array([0.7,0.72,0.74,0.76,0.78,0.8,
                       0.82,0.84,0.86,0.88,0.9])
    #Aklows /= 8.0
    ysoages = np.arange(-1.0,1.0001,0.2)+3.0
    MakePlot(sims,Aklow=Aklows,ysoage=ysoages,ax=None)
