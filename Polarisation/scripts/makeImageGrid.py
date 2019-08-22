'''
Make images of the cloud/sinks
Rebekka Bieri and Sam Geen, August 2018
'''

from startup import *

from pymses.utils import constants as C

import columndensity, sinks

from matplotlib import rc

rc('axes', labelsize=8)
rc('axes', labelsize=8,linewidth=0.5)
rc('xtick', labelsize=8)
rc('ytick', labelsize=8)

columndensity.IMSIZE = 512
IMSIZE = columndensity.IMSIZE

def FindLscale(boxlen):
    # Find the best match for the length scale bar
    lscales = np.array([1,2,5,10,20,50,100,200,500])
    # Ideal length is around 20% of the box length
    ideal = 0.2*boxlen
    # Find closest value to ideal
    diffs = np.abs(lscales - ideal)
    return int(lscales[diffs == diffs.min()][0])

def ProjectSinks(snap,los):
    up = columndensity.ups[los]
    across = columndensity.acrosses[los]
    sink = sinks.FindSinks(snap)
    sinkx = sink[across]
    sinky = sink[up]
    sinkm = np.atleast_1d(sink.mass)
    return sinkx, sinky, sinkm

def _createColDensMap(snap,los=None,wsink=False,imsize=IMSIZE):
    dmap = columndensity.DensityMap(snap,los,imsize=imsize)
    im = dmap.NH()
    boxlen = snap.info["boxlen"]
    if wsink:
        # Plot with sinks
        sinkx, sinky, sinkm = ProjectSinks(snap,los)
        return [im, sinkx, sinky, sinkm, boxlen]
    else:
        return [im, boxlen]

def MakeImageHamu(data,wsink,ax,dolengthscale,cmap,label=None):

    if wsink:
        im, sinkx, sinky, sinkm, boxlen = data
    else:
        im, boxlen = data
    
    # Plot image map
    hydro = "NH"
    yscale = hydrofuncs.yscale(hydro)
    if yscale == "log":
        im = np.log10(im)
        vmin, vmax = hydrofuncs.hydro_range(hydro)
    # Make image mesh
    xl, yl = im.shape
    xarr = np.arange(0,boxlen*1.0000001,boxlen/(xl-1.0))
    yarr = np.arange(0,boxlen*1.0000001,boxlen/(yl-1.0))
    
    # Plot the mesh
    #cmap = "YlGnBu_r"
    #cmap = "RdPu_r"
    cax = ax.pcolormesh(xarr,yarr,im,vmin=vmin,vmax=vmax,cmap=cmap)
    cax.set_rasterized(True)
    ax.set_xlim(0,boxlen)
    ax.set_ylim(0,boxlen)
   
    if wsink: 
        # Plot sinks
        try:
            lenm = len(np.atleast_1d(sinkm))
        except:
            lenm = 0
        if lenm > 0:
            area = np.pi * sinkm / 10.0
            ax.scatter(sinky,sinkx,s=area,c="w",alpha=0.5,edgecolors='w')

    # Add scale axis
    scalecol = "w"
    if dolengthscale:
        # length scale in pc (hopefully this is what units boxlen is in)
        lscale = FindLscale(boxlen)
        x1 = 0.1 * boxlen
        x2 = x1 + lscale
        y1 = 0.9 * boxlen
        y2 = y1
        ax.plot([x1,x2],[y1,y2],scalecol)
        ax.text(x2,y2, "  "+str(lscale)+" pc",color=scalecol,
                verticalalignment="center",fontsize="large")
    # Add label
    if label is not None:
        xt = 0.9 * boxlen
        yt = 0.1 * boxlen
        ax.text(xt,yt,label,
                horizontalalignment="right",
                verticalalignment="center",
                color=scalecol,fontsize="medium")
    return cax

def MakeFigure(simnames,time,name,los=None,wsink=False,nonamelabel=False,shape=None,dpi=100.0):
    nsims = len(simnames)
    ncols = 3
    nrows = 1
    while nrows*ncols < nsims:
        nrows += 1
        
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, 
                             sharex=True, sharey=True,
                             frameon=False)
    finches = IMSIZE/dpi

    createColDensMap = Hamu.Algorithm(_createColDensMap)
 
    # Run for all sims
    dolengthscale = True
    isim = -1
    for ax in axes.flatten():
        ax.set_axis_off()
        ax.set_aspect("equal","datalim")
        ax.axis("off")
    for simname in simnames:
        print "Running for sim", simname
        isim += 1
        cmap  = linestyles.ColourMap(simname)
        sim   = Hamu.Simulation(simname)
        snap  = sim.FindAtTime(time)
        ax    = axes.flatten()[isim]
       
        hydro = 'NH'
        data  = createColDensMap(snap,los,wsink,IMSIZE)       
        im    = MakeImageHamu(data,wsink,ax,dolengthscale,cmap,
                       label = linestyles.Label(simname))
        dolengthscale = False

    # Add colour bar for hydro variable
    print "Making colour bar..."
    # ---- This cax works well with nrow = 1, may need adjusting, could be improved 
    cax  = fig.add_axes([0.68, 0.1, 0.29, 0.03])
    cbar = fig.colorbar(im,cax=cax,orientation="horizontal")
    label = hydrofuncs.hydro_label((hydro))
    if not "xH" in hydro:
        label = "log("+label+")"
    cbar.set_label(label,fontsize="medium",color="k")
    cbar.ax.tick_params(labelsize="medium",labelcolor="k")
    cbar.solids.set_edgecolor("face")
    # Save figure
    suffix = ""
    if len(name) > 0:
        suffix = "_"+name
    #wsname = Hamu.SimData.Settings.Settings()["CurrentWorkspace"]
    folder = "../plots/vis/multiray/"
    MakeDirs(folder)
    figname = folder+"multiray_"+hydro+suffix+".png"
    print "Saving figure "+figname+"..."
    fig.subplots_adjust(hspace=0.0, wspace=0.0, 
                        left=0.0,right=1.0,
                        bottom=0.00,top=1.0)
    fig.set_size_inches(finches*ncols,finches*nrows)
    fig.savefig(figname,
                pad_inches=0,
                dpi=dpi)
    # Crop out borders
    #os.system("pdfcrop "+figname+" "+figname)
    print "Done!"

if __name__=="__main__":
    times = np.linspace(1,10,50) # in Myr
    simnames1 = ["IMF1_"+str(i).zfill(2) for i in [1,2,3,4,5]]
    simnames2 = ["IMF2_"+str(i).zfill(2) for i in [1,2,3,4,5]]
    simnamesm = ["MASS_"+str(i).zfill(2) for i in [1,2,3,4,5]]
    sims1 = [hamusims[simname] for simname in simnames1]
    sims2 = [hamusims[simname] for simname in simnames2]
    simsm = [hamusims[simname] for simname in simnamesm]
    wsink = True
    iout = 1
    for time in times:
        tcode = time/(unit_t/Myrins) 
        istr = str(iout).zfill(3)
        tstr = str(time)
        #MakeFigure(simnames1,tcode,name="rumgrid_1_"+tstr,los='z',wsink=wsink)
        #MakeFigure(simnames2,tcode,name="rumgrid_2_"+tstr,los='z',wsink=wsink)
        MakeFigure(simnamesm,tcode,name="rumgrid_m_"+tstr,los='z',wsink=wsink)
        iout += 1
    #MakeFigure(imf1sims,time,name="imf1",los='z',wsink=True)
    #MakeFigure(imf2sims,time,name="imf2",los='z',wsink=True)
