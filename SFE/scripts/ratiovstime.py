'''
Plot of stellar mass / mass in gas A_k>0.8 vs time 
'''

from startup import *

import columndensity, images, sfeobs, sinks, starsvsdensegas

from pymses.utils import constants as C

class SFEPlotter(starsvsdensegas.SFEPlotter):
    def Plot(self,los):
        data = self._MakeForLOS(los)
        x = data.times
        y = data.stars / (data.stars+data.gas)
        plt.plot(x,y,linestyle=self.Line(),color=self.Colour())
        plt.scatter(x[-1],y[-1],c=self.Colour(),edgecolors='none')

def MakePlot(sims,Aklow=0.8,allstars=False,ysoage=0.0):
    print "PLOTTING RATIO VS TIME FOR AKLOW=",Aklow,"ALLSTARS=",allstars
    plt.clf()
    # Simulations
    for sim in sims:
        plotter = SFEPlotter(sim,Aklow,allstars,ysoage,__name__)
        plotter.PlotAll()
    # Lada+ 2010 fit (by eye, since they don't give an offset)
    tlada = np.array([0.0,20.0]) 
    sfelada = np.array([0.1,0.1])
    plt.plot(tlada,sfelada,color=r"#888888",linestyle="--")
    # Do legends
    plotter.Legend()
    # Set up figure and output
    mgas = "$M_{\mathrm{gas}}(>A_{k} = "+str(Aklow)+")$"
    mstar = "$M_{*}$ / "
    plt.xlabel("Time / Myr")
    plt.ylabel(mstar+" / "+mgas)
    if allstars:
        plt.ylabel("$M_{*}$ (All stars) / "+Msolar)
    plt.yscale("log")
    plt.xlim([0,20])
    plt.ylim([1e-3,1.0])
    plotter.Save()
