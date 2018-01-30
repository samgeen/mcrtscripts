'''
Plot of stellar mass / mass in gas A_k>0.8 / yso age vs time 
'''

from startup import *

import columndensity, images, sfeobs, sinks, starsvsdensegas
import freefalltime

import ladadata

from pymses.utils import constants as C

p=2

rfact = {"L":1.5**p,
         "M":1.0**p,
         "S":0.75**p,
         "X":0.5**p}

class SFEPlotter(starsvsdensegas.SFEPlotter):
    def Plot(self,showerror):
        # TODO: Better interface to base class member variables
        def yfunc(data):
            return data.stars
        def xfunc(data):
            return data.gas
        #self.PlotError(xfunc,yfunc,linetype="scatter",dt=0.2,showerror=True)
        self.PlotError(xfunc,yfunc,xtff=False,
                       linetype="surfacesolid",dt=0.2,showerror=showerror)

def MakePlot(sims,Aklow=[0.8],allstars=False,ysoage=[0.0],ax=None):
    print "PLOTTING SFR VS TIME FOR AKLOW=",Aklow,"ALLSTARS=",allstars
    fig = None
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    # Simulations
    plotters = []
    for sim in sims:
        plotter = SFEPlotter(sim,Aklow,allstars,ysoage,__name__,
                             axis=ax,figure=fig)
        plotters.append(plotter)
    for plotter in plotters:
        plotter.Plot(showerror=True)
    # Lada+ 2010 data
    mgas = ladadata.lada["M0.8"]    
    myso = ladadata.lada["MYSO"]
    ax.scatter(mgas,myso,c="w",edgecolors="k",zorder=3,label="Lada et al (2010)")
    leg = ax.legend(frameon=False,fontsize="x-small",loc="lower left")
    # Lada+ 2010 fit (by eye, since they don't give an offset)
    tlada = np.array([0.0,20.0]) 
    sfelada = 0.095
    gaslada = np.array([1e1,2e4])
    starslada = gaslada * sfelada
    ax.plot(gaslada,starslada,color=r"#888888",linestyle="--")
    # Do legends
    plotter.Legend(rt=False,clouds="upper left",fontsize="x-small",
                   custom=leg)
    # Set up figure and output
    Aktxt = str(np.median(Aklow))
    mgas = "$M("+Aktxt+")$"
    mstar = "$M_{YSO}$"
    ax.set_xlabel(mgas+" / "+Msolar)
    ax.set_ylabel(mstar+" / "+Msolar)
    if allstars:
        ax.set_ylabel(mstar+" (All stars) / "+Msolar)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([40,2e4])
    ax.set_ylim([1,1e4])
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
