'''
Plot observed vs total SFE for the same clouds
Sam Geen, February 2017
'''

from startup import *
import totalsfevstime, sfrvstime

def Run():
    # Set up plot
    simnames = ["L-RT","M-RT","S-RT","XS-RT",
                "L-NRT","M-NRT","S-NRT","XS-NRT"]
    sims = [Hamu.Simulation(s) for s in simnames]
    fig, axes = plt.subplots(1,2,sharey=True)

    # Functions to use
    def tfunc(data):
        return data.stars / 1e4
    def ofunc(data):
        return data.stars / data.gas
    # Make lines
    for sim in sims:
        x = []
        yl = []
        ym = []
        yh = []
        # Total SFE
        tsfe = totalsfevstime.SFEPlotter(sim,Aklow=[0.0],allstars=True,
                                         ysoage=[0.0],name="")
        t, xmin, xmed, xmax = tsfe.MakeErrors(tfunc)
        # Observed SFE
        Aklows = np.array([0.7,0.72,0.74,0.76,0.78,0.8,
                           0.82,0.84,0.86,0.88,0.9])
        ysoages = np.arange(-1.0,1.0001,0.2)+3.0
        osfe = sfrvstime.SFEPlotter(sim,Aklow=Aklows,allstars=False,
                                    ysoage=ysoages,name="")
        t,ymin,ymed,ymax = osfe.MakeErrors(ofunc)
        # Observed SFE (Ak=0.1)
        Aklows = np.array([0.05,0.06,0.07,0.08,0.09,0.1,
                           0.11,0.12,0.13,0.14,0.15])
        ysoages = np.arange(-1.0,1.0001,0.2)+3.0
        osfe2 = sfrvstime.SFEPlotter(sim,Aklow=Aklows,allstars=False,
                                     ysoage=ysoages,name="")
        t,ymin2,ymed2,ymax2 = osfe2.MakeErrors(ofunc)
        # Plot OSFE / TSFE vs time
        ax = axes[0]
        x = t
        y = ymed / xmed
        yl = ymin / xmax
        yh = ymax / xmin
        line = linestyles.line(sim.Name())
        colour = linestyles.colour(sim.Name())
        ax.plot(x,y,linestyle=line,color=colour)
        ax.scatter(x[-1],y[-1],c=colour,edgecolors='none')
        # Plot fill around min/max                                          
        taba = np.concatenate((x,x[::-1])) # there and back again       
        pcontour = np.concatenate((yl,yh[::-1]))
        ax.fill(taba,pcontour,alpha=0.25,
                edgecolor='none',facecolor=colour)
        # TSFE vs OSFE (Ak=0.1)
        ax = axes[1]
        y = ymed2 / xmed
        yl = ymin2 / xmax
        yh = ymax2 / xmin
        line = linestyles.line(sim.Name())
        colour = linestyles.colour(sim.Name())
        ax.plot(x,y,linestyle=line,color=colour)
        ax.scatter(x[-1],y[-1],c=colour,edgecolors='none')
        # Plot fill around min/max                                          
        taba = np.concatenate((x,x[::-1])) # there and back again       
        pcontour = np.concatenate((yl,yh[::-1]))
        ax.fill(taba,pcontour,alpha=0.25,
                edgecolor='none',facecolor=colour)
    # OSFE vs TSFE
    ax = axes[0]
    # Lines showing 1:1 and 1:10
    #ax.plot([1e-5,10],[1e-5,10],linestyle=":",color="k")
    #ax.plot([1e-6,10],[1e-5,100],linestyle=":",color="k")
    # Legends
    lines, labels = linestyles.sizelegend()
    leg1 = axes[1].legend(lines,labels,fontsize="small",
                          frameon=False,loc="upper left")
    lines, labels = linestyles.rtlegend()
    leg2 = axes[0].legend(lines,labels,fontsize="small",
                          frameon=False,loc="lower center")
    #leg2.get_frame().set_linewidth(0.0)
    def MakeAktxt(ax,Ak):
        txt = "(A$_{\mathrm{k}}="+str(Ak)+")$"
        ax.text(0.95, 0.95,txt,
                horizontalalignment='right',
                verticalalignment='top',
                transform = ax.transAxes)
    ax.set_ylabel("OSFE / TSFE")
    ax.set_xlabel("Time / Myr")
    #ax.set_xlabel("TSFE")
    #ax.set_ylabel("OSFE")
    #ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([0,17])
    ax.set_ylim([5e-2,60])
    MakeAktxt(ax,0.8)
    ax.plot([0,25],[1,1],linestyle=":",color="k")
    # OSFE vs OSFE (Ak=0.1)
    ax = axes[1]
    #ax.set_ylabel("OSFE"+Aktxt+" / TSFE")
    ax.set_xlabel("Time / Myr")
    #ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([0,17])
    ax.set_ylim([2e-2,60])
    MakeAktxt(ax,0.1)
    ax.plot([0,25],[1,1],linestyle=":",color="k")
    # Output figure
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(15.0,5.0)
    fig.savefig("../plots/obsvstotalsfe.pdf")
        
        
if __name__=="__main__":
    Run()    
