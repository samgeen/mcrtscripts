'''
Run all functions for the project
Sam Geen, June 2016
'''

from startup import *
import compositeplots, columndensity, massvsav, plotsfe, radiusvsmass
import eachsinkvstime, images, gasvstime, mstarmgas2, mstarvstime
import ratiovstime, scattermstarmgas, sfrvstime, starsvsdensegas
import stellarmasscompare, totalsfevstime, obsvstotalmdense
import obsvstotalnyso, obsvstotalvoldens

Hamu.Workspace("HIISFE")

def coldensity():
    for simname in allsims:
        sim = Hamu.Simulation(simname)
        for snap in sim.Snapshots():
            # Pre-run all useful functions
            radius = columndensity.CloudRadius(snap) # los = z
            radius = columndensity.CloudRadius(snap,los='x')
            radius = columndensity.CloudRadius(snap,los='y')

def sfe(extlim):
    # TEST
    plotsfe.Akhigh = extlim
    simnames = [s+"-RT" for s in ["L","M","S","XS"]]
    dolegend = True
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        plotsfe.makeplot(sim,dolegend)
        dolegend = False

def cummassvsav():
    simnames = ["L-NRT","M-NRT","S-NRT","XS-NRT"]
    massvsav.Run(simnames,Akcutoff=5.0)

def radiusmass():
    simnames = ["L-NRT","M-NRT","S-NRT","XS-NRT"]
    radiusvsmass.Run(simnames)
    radiusvsmass.Run(simnames,Akcutoff=1.0)
    radiusvsmass.Run(simnames,Akcutoff=1.5)
    radiusvsmass.Run(simnames,Akcutoff=3.0)

def StarsvsDenseGas():
    simnames = ["L-NRT","M-NRT","S-NRT","XS-NRT",
                "L-RT","M-RT","S-RT","XS-RT"]
    sims = [Hamu.Simulation(s) for s in simnames]
    ysoages = np.arange(-1.0,1.0001,0.2)+3.0
    Aklows = np.array([0.05,0.06,0.07,0.08,0.09,0.1,
                       0.11,0.12,0.13,0.14,0.15])
    # Myso vs density ^ 2
    #mstarmgas2.MakePlot(sims,Aklow=Aklows,ysoage=ysoages,
    #                    powerindex=2.0)
    # Scatter plot of Myso vs Mgas
    simnames = ["L-RT","M-RT","S-RT","XS-RT"]
    Aklows = np.array([0.7,0.72,0.74,0.76,0.78,0.8,
                       0.82,0.84,0.86,0.88,0.9])

    sims = [Hamu.Simulation(s) for s in simnames]
    scattermstarmgas.MakePlot(sims,Aklow=Aklows,ysoage=ysoages)
    return
    mstarmgas2.MakePlot(sims,Aklow=Aklows,ysoage=ysoages)
    mstarmgas2.MakePlot(sims,ysoage=[3.0],powerindex=2.0,Aklow=[0.8])
    mstarmgas2.MakePlot(sims,ysoage=[3.0],powerindex=1.4,Aklow=[0.8])
    return
    eachsinkvstime.MakePlot(sims)
    totalsfevstime.MakePlot(sims,Aklow=1e-6)
    totalsfevstime.MakePlot(sims,Aklow=0.0)
    return
    mstarvstime.MakePlot(sims,ysoage=3.0)
    gasvstime.MakePlot(sims,Aklow=0.8)
    sfrvstime.MakePlot(sims,Aklow=0.8,ysoage=3.0)
    return
    ratiovstime.MakePlot(sims,ysoage=3.0,Aklow=0.1)
    ratiovstime.MakePlot(sims,ysoage=3.0,Aklow=0.8)
    dum = [mstarmgas2.MakePlot(sims,ysoage=a,powerindex=1.4) \
           for a in [0.5,1.0,2.0,3.0,4.0]]
    dum = [mstarmgas2.MakePlot(sims,ysoage=a,powerindex=2.0) \
           for a in [0.5,1.0,2.0,3.0,4.0]]
    gasvstime.MakePlot(sims,Aklow=0.1,ysoage=3.0)
    for mod in [ratiovstime,sfrvstime,starsvsdensegas]:#[gasvstime]:
        #mod.MakePlot(sims)
        #mod.MakePlot(sims,Aklow=0.0)
        #mod.MakePlot(sims,Aklow=0.1)
        #mod.MakePlot(sims,allstars=True)
        #mod.MakePlot(sims,Aklow=0.0,allstars=True)
        #mod.MakePlot(sims,Aklow=0.1,allstars=True)
        mod.MakePlot(sims,ysoage=0.5)
        mod.MakePlot(sims,ysoage=1.0)
        mod.MakePlot(sims,ysoage=2.0)
        mod.MakePlot(sims,ysoage=3.0)
        mod.MakePlot(sims,ysoage=4.0)

def Images():
    ysoage = 3.0
    times = [1.5,2.2,2.8,3.3]
    images.MakeFigure(["M-RT"],times,name="M-RT",los='z',ysoage=ysoage)
    times = [5.0,7.5,10.0,12.5]#,15.0]
    images.MakeFigure(["L-RT"],times,name="L-RT",los='z',ysoage=ysoage)
    images.MakeFigure(["L-NRT"],times,name="L-NRT",los='z',ysoage=ysoage)

def Production():
    # Comparing obs to total
    obsvstotalnyso.Run()
    obsvstotalmdense.Run()
    obsvstotalvoldens.Run()
    # Isolated SFR plots
    StarsvsDenseGas()
    # IMAGES
    ysoage = 3.0
    times = [4.2,8.4,12.6,16.0]
    images.MakeFigure("L-RT",times,name="L-RT",los='z',ysoage=ysoage,
                      nonamelabel=True,shape=(2,2))
    # List of things for final draft
    # TOTAL SINK MASSES
    #simnames = ["L-NRT","M-NRT","S-NRT","XS-NRT",
    #            "L-RT","M-RT","S-RT","XS-RT"]
    #simnames = ["L-NRT","M-NRT","L-RT","M-RT"]
    #sims = [Hamu.Simulation(s) for s in simnames]
    #eachsinkvstime.MakePlot(sims)
    # COMPOSITE FIGURES
    compositeplots.Mstartot()
    compositeplots.SFRPlot()

if __name__=="__main__":
    #sfe(5.0)
    #sfe(10.0)
    #sfe(20.0)
    #sfe(30.0)
    #cummassvsav()
    #radiusmass()
    #StarsvsDenseGas()
    #Images()
    Production()
