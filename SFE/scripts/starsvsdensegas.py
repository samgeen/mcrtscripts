'''
Plot of stellar mass vs mass in gas A_k>0.8 as in Fig 4 of Lada+ 2010
'''

from startup import *
import matplotlib

import scipy.spatial, concave

from pymses.filters import CellsToPoints
from pymses.utils import constants as C

from abc import ABCMeta, abstractmethod
import columndensity, freefalltime, images, sfeobs, sinks, ysos

def FindYSOMass(times,stars,ysoage=3.0):
    '''
    Find the mass in YSOs
    times in Myr, stars as cumulative mass in Msun, agemax in Myr
    '''
    ysomasses = []
    for t,s in zip(times,stars):
        ysomass = s
        if t > ysoage:
            ysomass = s - np.interp(t-ysoage,times,stars)
        ysomasses.append(ysomass)
    return np.array(ysomasses)

def islist(var):
    # Is var a list / other iterable with a zeroth element?
    try:
        a = var[0]
    except:
        return False
    return True

def _StarsInMap(snap,los,allstars=False):
    # Get dense gas mass
    dmap = columndensity.DensityMap(snap,los,NHlow=0.0)
    xl = columndensity.IMSIZE
    yl = columndensity.IMSIZE
    boxlen = snap.info["boxlen"]
    # Get sink positions
    sinkx, sinky, sinkm = images.ProjectSinks(snap.hamusnap,los)
    # Get all stars?
    if allstars:
        return np.sum(sinkm)
    def NHAtPos(x,y):
        # Find the column density in im at x,y
        # x,y position is in coordinates (0,boxlen)
        px = (xl*x/boxlen).astype(int)
        py = (yl*y/boxlen).astype(int)
        return dmap.NH()[px,py]
    # Get list of column densities
    NHs = NHAtPos(sinkx, sinky)
    # Find mass of sinks above NH threshold and below extinction limit
    NHthresh = AktoNH(0.1)
    ysomass = np.sum(sinkm[NHs >= NHthresh])
    # Return
    return ysomass
StarsInMap = Hamu.Algorithm(_StarsInMap)

def _GasInMap(snap,los,Aklows=[0.8]):
    dmap = columndensity.DensityMap(snap,los)
    masses = []
    for Ak in Aklows:
        dmap.NHlow(AktoNH(Ak))
        densemass = np.sum(dmap.Mass())
        masses.append(densemass)
    masses = np.array(masses)
    return masses
GasInMap = Hamu.Algorithm(_GasInMap)

def _GasInVolDens(snap,voldens):
    amr = snap.amr_source(["rho"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    vols = (cells.get_sizes())**3.0
    dens = cells["rho"]
    mass = dens*vols*snap.info["unit_mass"].express(C.Msun)
    if voldens > 0.0:
        mass = mass[dens > voldens]
    return np.sum(mass)
GasInVolDens = Hamu.Algorithm(_GasInVolDens)

class SimResults(object):
    times = None
    stars = None
    gas = None

class SFEPlotter(object):
    __metaclass__ = ABCMeta
    def __init__(self,sim,Aklow,allstars,ysoage,name,axis=None,figure=None,
                 voldens=None):
        self._sim = sim
        self._Aklow = np.atleast_1d(Aklow)
        self._allstars = allstars
        self._ysoage = np.atleast_1d(ysoage)
        self._voldens = voldens
        if self._ysoage.max() > 0.0:
            self._allstars = True
        self._colour = linestyles.colour(self._sim.Name())
        self._line = linestyles.line(self._sim.Name())
        self._name = name
        self._fig = figure
        if axis is None:
            self._fig = plt.figure()
            axis = self._fig.add_subplot(111)
        self._ax = axis
    
    def _MakeForLOS(self,los):
        '''
        Make values for a given LOS
        Aklow = optional value for varying errors
        ysoage = optional value for varying errors
        '''
        stars = []
        gas = []
        voldens = []
        times = []
        Aklow = self._Aklow
        ysoage = self._ysoage
        print "Plotting", self._sim.Name(), "for LOS", los
        for snap in self._sim.Snapshots():
            t = snaptime.Myr(snap)
            g = GasInMap(snap,los,self._Aklow)
            v = 0.0
            if self._voldens is not None:
                v = GasInVolDens(snap,self._voldens)
            v = g*0.0 + v # Make same size as gas
            if ysoage.max() > 0.0:
                ysos.sim = self._sim
                ds = ysos.MassInSnap(snap,los,ysoage)
                try:
                    s = ds.values()
                except:
                    s = ysoage*0.0
            else:
                s = np.array([StarsInMap(snap,los,allstars=self._allstars)])
            stars.append(s)
            gas.append(g)
            voldens.append(v)
            times.append(t)
        # Convert to numpy arrays, shape(len(times),len(ysoage or Aklow))
        stars = np.array(stars)
        gas = np.array(gas)
        voldens = np.array(voldens)
        times = np.array(times)
        # Make output structure
        results = SimResults()
        results.times = times
        results.stars = stars
        results.gas = gas
        results.voldens = voldens
        return results

    def Colour(self):
        return self._colour

    def Line(self):
        return self._line
    
    def Figure(self):
        return self._fig

    def Axis(self):
        return self._ax

    def Legend(self,clouds="lower right",rt="center right",
               custom=False,customtxt="",fontsize="small"):
        # Do legends
        ax = self._ax
        if clouds:
            lines, labels = linestyles.sizelegend()
            leg1 = ax.legend(lines,labels,fontsize=fontsize,
                             frameon=False,loc=clouds)
        if rt:
            lines, labels = linestyles.rtlegend()
            leg2 = ax.legend(lines,labels,fontsize=fontsize,
                             frameon=False,loc=rt)
        if custom:
            if type(custom) == type("string"):
                # Is a string for a location?
                leg3 = ax.legend([],[],fontsize=fontsize,
                                 title=customtxt,
                                 frameon=False,loc=custom)
            else:
                # Is a legend itself?
                leg3 = custom
        if clouds:
            ax.add_artist(leg1)
        if rt:
            ax.add_artist(leg2)
        if custom:
            ax.add_artist(leg3)

        
    def Save(self):
        if self._fig is not None:
            name = self.Filename()
            print "Saving",name,"..."
            self._fig.savefig(name)
        else:
            print "Figure not controlled by this class, not saving here..."

    def Filename(self):
        ysoage = self._ysoage
        Aklow = self._Aklow
        # Pick median value to avoid whole array in filename
        if len(ysoage) > 1:
            ysoage = ysoage[len(ysoage)//2]
        if len(Aklow) > 1:
            Aklow = Aklow[len(Aklow)//2]
        Aktxt = "_Ak"+str(Aklow).replace(".","p")
        allstarstxt = ""
        if self._allstars:
            allstarstxt = "_allstars"
        if self._ysoage.max() > 0.0:
            allstarstxt = "_ysos"+str(ysoage).replace(".","p")+"Myr"
        name = "../plots/"+self._name+Aktxt+allstarstxt+".pdf"
        return name
                    
    def MakeErrors(self,func):
        '''
        Make min, median, max of func vs time
        func: function with argument SimResults
        '''
        outputs = {}
        test = SimResults()
        # Make values for the given function
        for los in ["x","y","z"]:
            results = self._MakeForLOS(los)
            for itime in range(0,len(results.times)):
                gas = results.gas[itime,:]
                voldens = results.voldens[itime,:]
                stars = results.stars[itime,:]
                lg = len(gas)
                ls = len(stars)
                # Make a grid of all possible gas masses, yso masses
                test.stars = np.tile(stars,lg)
                test.gas   = np.repeat(gas,ls)
                test.voldens = np.repeat(voldens,ls)
                test.times = np.zeros((ls*lg))+results.times[itime]
                if not itime in outputs:
                    outputs[itime] = np.array([])
                outputs[itime] = np.concatenate((outputs[itime].flatten(),
                                                func(test).flatten()))
        # Now find min/median/max
        mins = []
        meds = []
        maxs = []
        for itime in range(0,len(results.times)):
            mins.append(outputs[itime].min())
            meds.append(np.median(outputs[itime]))
            maxs.append(outputs[itime].max())
        times = results.times
        mins = np.array(mins)
        meds = np.array(meds)
        maxs = np.array(maxs)
        return times, mins, meds, maxs

    def PlotError(self,xfunc,yfunc,xtff,
                  linetype="line",showerror=True,dt=None):
        # Plot median
        tff = freefalltime.Tff(self._sim,Myr=True)
        if linetype == "scatter" or linetype == "surface" or linetype == "surfacesolid":
            t, xlow, xmed, xhigh = self.MakeErrors(xfunc)
            t, ylow, ymed, yhigh = self.MakeErrors(yfunc)
            # Sample values over time
            if dt is not None:
                # Find new t spaced evenly until the last time value
                newt = np.arange(0.0,t[-1],dt)
                xhigh = np.interp(newt,t,xhigh)
                xmed = np.interp(newt,t,xmed)
                xlow = np.interp(newt,t,xlow)
                yhigh = np.interp(newt,t,yhigh)
                ymed  = np.interp(newt,t,ymed)
                ylow  = np.interp(newt,t,ylow)
            # TODO: OPTIMISE ERROR SPACE SEARCHING
            if not showerror:
                self._ax.scatter(xmed,ymed,c=self.Colour(),edgecolors='none')
            else:
                if linetype == "scatter":
                    self._ax.errorbar(xmed, ymed, 
                                      xerr=[xlow, xhigh], yerr=[ylow, yhigh], 
                                      fmt='o',
                                      color=self.Colour(),
                                      alpha=0.25,
                                      elinewidth=1)
                elif linetype == "surfacesolid":
                    self._PlotSurfaceErrorsSolid(xlow,xmed,xhigh,
                                                 ylow,ymed,yhigh,
                                                 self.Colour())
                else:
                    self._PlotSurfaceErrors(xlow,xmed,xhigh,
                                            ylow,ymed,yhigh,
                                            self.Colour())
        else:
            # Default: plot line
            # Only do y axis errors
            # d = dummy
            d1, d2, x, d3 = self.MakeErrors(xfunc)
            if xtff:
                x -= tff
            d4, low, med, high = self.MakeErrors(yfunc)
            self._ax.plot(x,med,linestyle=self.Line(),color=self.Colour())
            self._ax.scatter(x[-1],med[-1],c=self.Colour(),edgecolors='none')
            # Plot fill around min/max
            if showerror:
                taba = np.concatenate((x,x[::-1])) # there and back again
                pcontour = np.concatenate((low,high[::-1]))
                self._ax.fill(taba,pcontour,alpha=0.25,
                              edgecolor='none',facecolor=self.Colour())
                
    def _PlotSurfaceErrors(self,xlow,xmed,xhigh,ylow,ymed,yhigh,colour):
        '''
        Plot errors as a surface summing over ranges found
        Use a 2D Gaussian to populate the surface
        '''
        # The Object Orientation gods will smite me for this one for sure
        xmin = np.log10(xlow[xlow > 0.0].min())
        xmax = np.log10(xhigh.max())
        ymin = np.log10(ylow[ylow > 0.0].min())
        ymax = np.log10(yhigh.max())
        xr = xmax - xmin
        yr = ymax - ymin
        ngrid = 100
        numpts = 1000
        grid = np.zeros((ngrid,ngrid))
        xs = np.logspace(xmin,xmax,ngrid)
        ys = np.logspace(ymin,ymax,ngrid)
        xgrid, ygrid = np.meshgrid(xs,ys)
        def ctog(xin,yin):
            # Coordinate in xy space to grid space
            mask = (xin > 0.0) * (yin > 0.0)
            x = xin[mask]
            y = yin[mask]
            xind = ((np.log10(x) - xmin) / xr * ngrid).astype(int)
            yind = ((np.log10(y) - ymin) / yr * ngrid).astype(int)
            mask = (xind >= 0) * (yind >= 0) * (xind < ngrid) * (yind < ngrid)
            xind = xind[mask]
            yind = yind[mask]
            return (xind,yind)
        for xl,xm,xh,yl,ym,yh in zip(xlow,xmed,xhigh,ylow,ymed,yhigh):
            # Error is the average half-distance to min/max from median 
            xsig = 0.25 * (xh-xl)
            ysig = 0.25 * (yh-yl)
            # Make point cloud
            if xm > 0.0 and xsig > 0.0 and ym > 0.0 and ysig > 0.0:
                x = np.random.normal(xm,xsig,numpts)
                y = np.random.normal(ym,ysig,numpts)
                # HACK - just show median
                x = np.array([xm])
                y = np.array([ym])
                # Add to surface
                grid[ctog(x,y)] += 1.0
        # Plot the surface
        # Set up colour map
        mapname = "Custom"+colour
        r = int(colour[1:3],16)/255.0
        g = int(colour[3:5],16)/255.0
        b = int(colour[5:7],16)/255.0
        colours = [(1.0,1.0,1.0),(r,g,b)]
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(mapname,colours)
        # TODO: FIGURE OUT NORMALISATION ACROSS SIMULATIONS
        # DO IT ONCE AND FIGURE OUT LATER?
        self._ax.pcolormesh(xgrid,ygrid,grid.T,
                            alpha=0.25, cmap=cmap)


    def _PlotSurfaceErrorsSolid(self,xlow,xmed,xhigh,ylow,ymed,yhigh,colour):
        '''
        Plot errors as a surface summing over ranges found
        Use a 2D Gaussian to populate the surface
        '''
        dc = 0.001*np.pi
        c = np.arange(0.0,2.0*np.pi+dc+0.0001,dc)
        cx = np.sin(c)
        cy = np.cos(c)
        xr = xhigh - xlow
        yr = yhigh - ylow
        numpts = len(c)
        numvals = len(xlow)
        points = np.zeros((numpts*numvals,2))
        iv = 0
        # HACK
        ylim = [1,1e4]
        # Find all possible edge points around points
        def lighten(col):
            def lighthex(h):
                i = int("0x"+h,0)
                i = 255//2 + i//2
                return str(hex(i))[2:]
            r = lighthex(col[1:3])
            g = lighthex(col[3:5])
            b = lighthex(col[5:7])
            return r"#"+r+g+b
        def drawpoly(px,py,mx,my):
            nox = np.abs(np.mean(px) - mx) < 1e-4
            noy = np.abs(np.mean(py) - my) < 1e-4
            col = lighten(self.Colour())
            lw = 3
            if nox and noy:
                # No x or y variance, don't draw
                return
            elif nox:
                # Draw line in y
                self._ax.plot([mx,mx],[py.max(),py.min()],color=col,
                              linewidth=lw,zorder=1)
            elif noy:
                # Draw line in x
                self._ax.plot([px.max(),px.min()],[my,my],color=col,
                              linewidth=lw,zorder=1)
            else:
                ax = np.mean(px)
                ay = np.mean(py)
                px = np.concatenate((px,np.atleast_1d(mx-0.01*ax)))
                py = np.concatenate((py,np.atleast_1d(my-0.01*ay)))
                self._ax.fill(px,py,
                              edgecolor='none',
                              facecolor=col,
                              zorder=1)
        # Run through points
        for ip in range(0,numvals):
            xl = xlow[ip]
            yl = ylow[ip]
            xh = xhigh[ip]
            yh = yhigh[ip]
            xm = xmed[ip]
            ym = ymed[ip]
            if ym > ylim[0]:
                if yl <= 0.0:
                    yl = 1e-2
                # Quadrants on clock
                # 12:00 to 3:00
                drawpoly(cx[0:numpts//4+2]*(xh-xm)+xm,
                         cy[0:numpts//4+2]*(yh-ym)+ym,xm,ym)
                # 3 to 6
                drawpoly(cx[numpts//4:numpts//2+2]*(xh-xm)+xm,
                         cy[numpts//4:numpts//2+2]*(ym-yl)+ym,xm,ym)
                # 6 to 9
                drawpoly(cx[numpts//2:3*numpts//4+2]*(xm-xl)+xm,
                         cy[numpts//2:3*numpts//4+2]*(ym-yl)+ym,xm,ym)
                # 9 to 12
                drawpoly(cx[3*numpts//4:numpts]*(xm-xl)+xm,
                         cy[3*numpts//4:numpts]*(yh-ym)+ym,xm,ym)

        # Plot scatter for individual points
        self._ax.scatter(xmed,ymed,c=self.Colour(),edgecolors='none',
                         zorder=2)


    @abstractmethod
    def Plot(self):
        pass

class SFEPlotterStarsvsDenseGas(SFEPlotter):
    '''
    Plot stars vs dense gas
    '''
    def Plot(self,los):
        data = self._MakeForLOS(los)
        self._ax.plot(data.gas,data.stars,
                 linestyle=self.Line(),color=self.Colour())
        self._ax.scatter(data.gas[-1],data.stars[-1],
                    c=self.Colour(),edgecolors='none')

def MakePlot(sims,Aklow=0.8,allstars=False,ysoage=0.0):
    print "PLOTTING FOR AKLOW=",Aklow,"ALLSTARS=",allstars
    plt.clf()
    # Simulations
    for sim in sims:
        plotter = SFEPlotterStarsvsDenseGas(sim,Aklow,allstars,ysoage,__name__)
        plotter.PlotAll()
    # Lada+ 2010 fit (by eye, since they don't give an offset)
    mcloud = np.array([1e1,1e5]) 
    myso = 0.5 * np.array([1.9,1.9e4]) # Mass = 0.5 NYSOs
    plt.plot(mcloud, myso,color=r"#888888",linestyle="--")
    # Do legend
    plotter.Legend()
    # Set up figure and output
    plt.xlabel("$M_{\mathrm{gas}}(>A_{k} = "+str(Aklow)+")$ / "+Msolar)
    plt.ylabel("$M_{*}$ / "+Msolar)
    if allstars:
        plt.ylabel("$M_{*}$ (All stars) / "+Msolar)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim([10,1e5])
    plt.ylim([1,3e4])
    plotter.Save()
    
