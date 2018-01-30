'''
Plot observed vs total SFE for the same clouds
Sam Geen, February 2017
'''

from startup import *
import starsvsdensegas,totalsfevstime, sfrvstime

class Mdense(starsvsdensegas.SFEPlotter):
    def Plot(self):
        def yfunc(data):
            return data.gas
        def xfunc(data):
            return data.times
        self.PlotError(xfunc,yfunc)

def Run():
    # Set up plot
    #simnames = ["L-RT","M-RT","S-RT","XS-RT",
    #            "L-NRT","M-NRT","S-NRT","XS-NRT"]
    simnames = ["L-RT","M-RT","S-RT","XS-RT"]
    sims = [Hamu.Simulation(s) for s in simnames]
    fig, axes = plt.subplots(1,1)

    # Functions to use
    def ofunc(data):
        return data.gas
    # Make lines
    Aks = np.arange(0.1,1.0,0.05)
    ysoages = np.arange(-1.0,1.0001,0.2)+3.0
    for sim in sims:
        x = []
        yl = []
        ym = []
        yh = []
        # Total mass
        #tsfe = totalsfevstime.SFEPlotter(sim,Aklow=[0.0],allstars=True,
        #                                 ysoage=[0.0],name="")
        #t, xmin, xmed, xmax = tsfe.MakeErrors(tfunc)
        mtot = 1e4
        for Ak in Aks:
            print "Running for A_k=",Ak
            # Observed SFE
            Aklows = [Ak]
            mdense = Mdense(sim,Aklow=Aklows,allstars=False,
                            ysoage=ysoages,name="")
            t,ymin,ymed,ymax = mdense.MakeErrors(ofunc)
            x.append(Ak)
            cl = ymin/mtot
            cl = cl[np.isfinite(cl)]
            cl = cl[cl > 0]
            ch = ymax/mtot
            ch = ch[np.isfinite(ch)]
            yl.append(cl.min())
            ym.append((ymed/mtot)[len(ymed)//2])
            yh.append(ch.max())
        x = np.array(x)
        yl = np.array(yl)
        ym = np.array(ym)
        yh = np.array(yh)
        # Plot OSFE / TSFE vs time
        ax = axes
        line = linestyles.line(sim.Name())
        colour = linestyles.colour(sim.Name())
        ax.plot(x,ym,linestyle=line,color=colour)
        ax.scatter(x[-1],ym[-1],c=colour,edgecolors='none')
        # Plot fill around min/max                                          
        taba = np.concatenate((x,x[::-1])) # there and back again       
        pcontour = np.concatenate((yl,yh[::-1]))
        ax.fill(taba,pcontour,alpha=0.25,
                edgecolor='none',facecolor=colour)
    # OSFE vs TSFE
    ax = axes
    # Legends
    lines, labels = linestyles.sizelegend()
    leg1 = axes.legend(lines,labels,fontsize="small",
                       frameon=False,loc="upper left")
    #lines, labels = linestyles.rtlegend()
    #leg2 = axes.legend(lines,labels,fontsize="small",
    #                   frameon=False,loc="lower center")
    #ax.add_artist(leg1)
    ax.set_ylabel("$M_{A_k} / 10^4~$"+Msolar)
    ax.set_xlabel("$A_{\mathrm{k}}$")
    #ax.set_xlabel("TSFE")
    #ax.set_ylabel("OSFE")
    #ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([0,1])
    ax.set_ylim([2e-3,20])
    ax.plot([0,1],[1,1],linestyle=":",color="k")
    # Output figure
    #fig.subplots_adjust(wspace=0, hspace=0)
    #fig.set_size_inches(15.0,5.0)
    fig.savefig("../plots/akvary_mdensevsmsim.pdf")
        
        
if __name__=="__main__":
    Run()    