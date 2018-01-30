'''
Plot composite plots for the paper
Sam Geen, November 2016
'''

from startup import *
import totalsfevstime, eachsinkvstime
import sfrvstime, gasvstime
import planckcompare, cumvoldens

def Mstartot():
    '''
    Plot total Mstar, eachsink
    '''
    # SET UP
    simnames = ["L-RT","M-RT","S-RT","XS-RT",
            "L-NRT","M-NRT","S-NRT","XS-NRT"]
    sims = [Hamu.Simulation(s) for s in simnames]
    fig, axes = plt.subplots(2, sharex=True)
    # TOTAL SFE
    totalsfevstime.MakePlot(sims,Aklow=[0.0],ax=axes[0])
    # EACH SINK VS TIME
    sims = [Hamu.Simulation(s) for s in ["L-RT","L-NRT"]]
    eachsinkvstime.MakeOnePlot(sims,ax=axes[1])
    # SAVE
    axes[0].set_aspect('auto')
    axes[1].set_aspect('auto')
    axes[1].set_xlabel("Time - $t_{ff}$ / Myr")
    fig.set_size_inches(8.0, 12.0)
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig("../plots/mstartots.pdf")

def SFRPlot():
    '''
    Plot SFR, gas mass
    '''
    # SET UP
    simnames = ["L-RT","M-RT","S-RT","XS-RT",
            "L-NRT","M-NRT","S-NRT","XS-NRT"]
    sims = [Hamu.Simulation(s) for s in simnames]
    fig, axes = plt.subplots(2, sharex=True)
    Aklows = np.array([0.7,0.72,0.74,0.76,0.78,0.8,
                       0.82,0.84,0.86,0.88,0.9])
    ysoages = np.arange(-1.0,1.0001,0.2)+3.0
    # SFR
    sfrvstime.MakePlot(sims,Aklow=Aklows,ysoage=ysoages,ax=axes[0],
                       usetff=True)
    # GAS MASS
    gasvstime.MakePlot(sims,Aklow=Aklows,ax=axes[1],
                       usetff=True)
    # SAVE
    axes[0].set_aspect('auto')
    axes[1].set_aspect('auto')
    axes[1].set_xlabel("Time / Myr")
    fig.set_size_inches(8.0, 12.0)
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig("../plots/sfrplots.pdf")

def TffPlot():
    '''
    Plot SFE*tff vs time for 
    '''
    # SET UP
    simnames = ["L-RT","M-RT","S-RT","XS-RT",
            "L-NRT","M-NRT","S-NRT","XS-NRT"]
    sims = [Hamu.Simulation(s) for s in simnames]
    fig, axes = plt.subplots(1, sharex=True)
    Aklows = np.array([0.05,0.06,0.07,0.08,0.09,0.1,
                       0.11,0.12,0.13,0.14,0.15])
    ysoages = np.arange(-1.0,1.0001,0.2)+3.0
    # SFR
    sfrvstime.MakePlot(sims,Aklow=Aklows,ysoage=ysoages,ax=axes,usetff=True)
    # SAVE
    axes.set_aspect('auto')
    #fig.set_size_inches(8.0, 12.0)
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig("../plots/sfrtff.pdf")

def CumulDens():
    '''
    Plot cumulative densities
    '''
    # TODO: unfuck axis sharing
    fig, axes = plt.subplots(2,3, sharey=True)
    # Planck compare
    axes = axes.flatten()
    planckcompare.MakePlot(axes[0],PSFon=True,distance="near")
    planckcompare.MakePlot(axes[1],PSFon=True,distance="mid")
    planckcompare.MakePlot(axes[2],PSFon=True,distance="far")
    planckcompare.MakePlot(axes[3],PSFon=True,distance="herschel",
                           herschel=True)
    planckcompare.MakePlot(axes[4],PSFon=False,distance="NONE",tffs=[1,2,4])
    # Cumulative volume density
    cumvoldens.MakePlot(axes[5],tffs=[1,2,4])
    # SAVE
    axes[0].set_ylabel("Cumulative Mass / "+Msolar)
    axes[3].set_ylabel("Cumulative Mass / "+Msolar)
    axes[0].set_aspect('auto')
    axes[1].set_aspect('auto')
    axes[2].set_aspect('auto')
    fig.set_size_inches(17.0, 10.0)
    fig.subplots_adjust(wspace=0)
    fig.savefig("../plots/cumuldens.pdf")

if __name__=="__main__":
    # Test
    SFRPlot()
    #Mstartot()
    #CumulDens()
