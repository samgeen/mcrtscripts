'''
Plot observed vs total NYSO for the same clouds
Sam Geen, February 2017
'''

from startup import *
import freefalltime, totalsfevstime, sfrvstime, starsvsdensegas

class NYSO(starsvsdensegas.SFEPlotter):
    def Plot(self):
        def yfunc(data):
            return data.stars
        def xfunc(data):
            return data.times
        self.PlotError(xfunc,yfunc)

class Mdense(starsvsdensegas.SFEPlotter):
    def Plot(self):
        def yfunc(data):
            return data.gas
        def xfunc(data):
            return data.times
        self.PlotError(xfunc,yfunc)

def Run():
    # Set up plot
    simnames = ["L-RT","M-RT","S-RT","XS-RT",
                "L-NRT","M-NRT","S-NRT","XS-NRT"]
    sims = [Hamu.Simulation(s) for s in simnames]
    fig, axes = plt.subplots(2,1,sharex=True)

    # Functions to use
    def sfunc(data):
        return data.stars
    def gfunc(data):
        return data.gas
    # Make Ak label
    def MakeAktxt(ax,Ak):
        txt = "(A$_{\mathrm{k}}="+str(Ak)+")$"
        ax.text(0.05, 0.05,txt,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform = ax.transAxes)
    # Make lines
    for iax in [0,1]:
        for sim in sims:
            ax = axes[iax]
            if iax == 0:
                Akmed = 0.8
                Aklows = np.array([0.7,0.72,0.74,0.76,0.78,0.8,
                                   0.82,0.84,0.86,0.88,0.9])
            else:
                Akmed = 0.1
                Aklows = np.array([0.05,0.06,0.07,0.08,0.09,0.1,
                                   0.11,0.12,0.13,0.14,0.15])
            x = []
            yl = []
            ym = []
            yh = []
            # Observed YSO count
            ysoages = np.array([3.0]) # Doesn't matter for gas mass
            ogas = Mdense(sim,Aklow=Aklows,allstars=False,
                          ysoage=ysoages,name="")
            t,ymin,ymed,ymax = ogas.MakeErrors(gfunc)
            # Initial gas mass
            xmin = ymin*0+1e4
            xmed = ymed*0+1e4
            xmax = ymax*0+1e4
            # Plot ratio vs time
            tff = freefalltime.Tff(sim,Myr=True)
            x = t-tff
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
            # Lines showing 1:1 and 1:10
            #ax.plot([1e-5,10],[1e-5,10],linestyle=":",color="k")
            #ax.plot([1e-6,10],[1e-5,100],linestyle=":",color="k")
        # Legends
        if iax == 0:
            lines, labels = linestyles.sizelegend()
            leg1 = ax.legend(lines,labels,fontsize="x-small",
                             frameon=False,loc="upper right")
        else:
            lines, labels = linestyles.rtlegend()
            leg2 = ax.legend(lines,labels,fontsize="x-small",
                             frameon=False,loc="lower right")
        MakeAktxt(ax,Akmed)
        ax.set_ylabel("$M_{"+str(Akmed)+"}/10^{4}~$"+Msolar)
        ax.set_xlabel("Time-$t_{ff}$ / Myr")
        ax.set_yscale("log")
        xl = [-4.25,17-4.25]
        yl = [5e-3,2]
        ax.plot(xl,[1,1],linestyle=":",color="k")
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        xl[0] = 4.0
        ax.fill([xl[0],xl[1],xl[1],xl[0]], 
                [yl[0],yl[0],yl[1],yl[1]], 
                color=r"#888888",hatch='\\',alpha=0.2)
    # Output figure
    fig.set_size_inches(8.0,12.0)
    fig.subplots_adjust(wspace=0, hspace=0)
    #fig.set_size_inches(15.0,5.0)
    fig.savefig("../plots/obsvstotalmdense.pdf")

if __name__=="__main__":
    Run()    