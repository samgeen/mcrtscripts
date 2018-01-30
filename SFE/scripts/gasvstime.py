'''
Plot of mass in gas A_k>0.8 vs time 
'''

from startup import *

import columndensity, images, sfeobs, sinks, starsvsdensegas

from pymses.utils import constants as C

class SFEPlotter(starsvsdensegas.SFEPlotter):
    def __init__(self,sim,Aklow,allstars,ysoage,name,axis=None,figure=None):
        super(SFEPlotter, self).__init__(sim,Aklow,allstars,ysoage,name,
                                         axis=axis,figure=figure)

    def Plot(self,usetff=True):
        # TODO: Better interface to base class member variables                 
        def yfunc(data):
            return data.gas
        def xfunc(data):
            return data.times
        self.PlotError(xfunc,yfunc,usetff)

def MakePlot(sims,Aklow=0.8,allstars=False,ysoage=0.0,ax=None,usetff=True):
    print "PLOTTING GAS MASS VS TIME FOR AKLOW=",Aklow,"ALLSTARS=",allstars
    fig = None
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    # Simulations
    for sim in sims:
        plotter = SFEPlotter(sim,Aklow,allstars,ysoage,__name__,
                             axis=ax,figure=fig)
        plotter.Plot(usetff)
    # Do legends
    plotter.Legend(clouds="upper right",rt="lower right")
    # Set up figure and output                      
    Aktxt = str(np.median(Aklow))
    mgas = "$M_{"+Aktxt+"}$"
    mstar = "$M_{*}$"
    ysoagetxt = "(YSO Age / 1 Myr)"
    if fig is None:
        plt.xlabel("Time / Myr")
    ax.set_ylabel(mgas+" / "+Msolar)
    ax.set_yscale("log")
    ax.set_xlim([0-4.25,20-4.25])
    ax.set_ylim([1e2,1e4])
    shade4Myr.run(ax)
    # Lada+ 2010 fit (by eye, since they don't give an offset)
    tlada = list(ax.get_xlim())
    sfelada = np.array([1e3,1e3])
    ax.plot(tlada,sfelada,color=r"#888888",linestyle="--")

    plotter.Save()
