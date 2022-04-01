'''
Make images of the cloud/sinks
Rebekka Bieri and Sam Geen, August 2018
'''

from startup import *

from pymses.utils import constants as C

import columndensity, rayMap, sliceMap, sinks, ysos
import starrelations, listfigures, findproperties
import makefits

from matplotlib import rc
import matplotlib.cm as mplcm
import matplotlib.patheffects as PathEffects

import rdmfile

rc('axes', labelsize=8)
rc('axes', labelsize=8,linewidth=0.5)
rc('xtick', labelsize=8)
rc('ytick', labelsize=8)

IMSIZE = columndensity.IMSIZE
OUTLINEWIDTH = 1

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

def _createColDensMap(snap,los=None,hydro='rho',wsink=False,zoom=1.0):
    dmap = columndensity.DensityMap(snap,los,zoom=zoom)
    im = dmap.NH()
    boxlen = snap.info["boxlen"]
    return [im, boxlen*zoom]

def _createColDensMap_sink(snap,los=None,hydro='rho',wsink=True,zoom=1.0):
    ro = snap.RawData()
    dmap = columndensity.DensityMap(snap,los,zoom=zoom)
    cx, cy, im = dmap.NH()
    boxlen = ro.info["boxlen"]
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(ro,los)
    # Shift sink position if zooming in
    print("TODO: ADD CENTRED COLUMN DENSITY MAPS")
    #sinkx -= (cx - 0.5*zoom)*boxlen
    #sinky -= (cy - 0.5*zoom)*boxlen
    sinkx -= 0.5*(1.0-zoom)*boxlen
    sinky -= 0.5*(1.0-zoom)*boxlen
    return [im, sinkx, sinky, sinkm, boxlen*zoom]

def _createRayMap(snap,los=None,hydro='rho',wsink=False,zoom=1.0,starC=False):
    ro = snap.RawData()
    dmap = rayMap.RayTraceMap(snap,hydro,los,zoom=zoom,starC=starC)
    im   = dmap.getRaytraceMap()
    boxlen = ro.info["boxlen"]
    return [im, boxlen*zoom]

def _createRayMap_sink(snap,los=None,hydro='rho',wsink=True,zoom=1.0,starC=False):
    ro = snap.RawData()
    dmap = rayMap.RayTraceMap(snap,hydro,los,zoom=zoom,starC=starC)
    cx, cy,im   = dmap.getRaytraceMap()
    boxlen = ro.info["boxlen"]
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(ro,los)
    # Shift sink position if zooming in
    sinkx -= (cx - 0.5*zoom)*boxlen
    sinky -= (cy - 0.5*zoom)*boxlen
    return [im, sinkx, sinky, sinkm, boxlen*zoom]

def _createSliceMap(snap,hydro='rho',los=None,zoom=1.0,starC=False):
    ro = snap.RawData()
    dmap = sliceMap.SliceMap(snap,hydro,los,zoom=zoom,starC=starC)
    cx, cy, im   = dmap.getSliceMap()
    boxlen = ro.info["boxlen"] 
    return [im, boxlen*zoom]

def _createSliceMap_sink(snap,hydro='rho',los=None,zoom=1.0,starC=False):
    ro = snap.RawData()
    dmap = sliceMap.SliceMap(snap,hydro,los,zoom=zoom,starC=starC)
    cx, cy, im   = dmap.getSliceMap()
    boxlen = ro.info["boxlen"] 
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(ro,los)
    # Shift sink position if zooming in
    sinkx -= (cx - 0.5*zoom)*boxlen
    sinky -= (cy - 0.5*zoom)*boxlen
    #sinkx -= 0.5*(1.0-zoom)*boxlen
    #sinky -= 0.5*(1.0-zoom)*boxlen
    return [im, sinkx, sinky, sinkm, boxlen*zoom]

def rgb(r,g,b):
    return (float(r) / 255.0, float(g) / 255.0, float(b)/255.0)

def MakeImage(datas,hydros,snap,wsink,ax,dolengthscale,cmap,plottime=False,timeL=None,label=None,
                  starsink=None,rdm=None,contours=[],zoombox=-1,starC=False,zoom=1.0,Slice=False):

    # Get sink info
    stellar = stellars.FindStellar(snap)
    if len(stellar.mass) == 0:
        return 0.0
    imax = np.where(stellar.mass == np.max(stellar.mass))[0][0]
    sinkid = stellar.sinkid[imax]-1
    sink = sinks.FindSinks(snap)
    boxlen = snap.RawData().info["boxlen"]
    levelmax = snap.RawData().info["levelmax"]
    mincellsize = boxlen / 2.0**levelmax
    smass = stellar.mass
    starsinkid = stellar.sinkid[np.where(smass == smass.max())]
    starsink = np.where(sink.id == starsinkid)[0]
    starpos = np.array([sink.x[starsinkid],sink.y[starsinkid],sink.z[starsinkid]])/boxlen

    ims = []
    for data, hydro in zip(datas, hydros):
        # NOTE: for different hydro variables, only im should be different
        if wsink:
            im, sinkx, sinky, sinkm, zoomedboxlen = data
        else:
            im, zoomedboxlen = data

        # Specific hack for Lcool (ignore heating, cooling = +ve)
        if hydro == "Lcool" or "emission" in hydro:
            #if im.max() != im.min():
            try:
                minpositive = im[im > 0.0].min()*0.1
            except:
                minpositive = -1337.0
            im[im <= 0.0] = minpositive
        ims.append(im)

    finalim = None
    ihydro = -1
    # Make a colourblind-safer RGB map (3-class Dark2 in ColorBrewer2.org)
    threecolour = [rgb(27,158,119),
                    rgb(217,95,2),
                   rgb(117,112,179)][::-1]
    #threecolour = [rgb(27,158,119),
    #                rgb(217,95,2),
    #               rgb(231,41,138)][::-1]
    for im, hydro in zip(ims, hydros):
        ihydro += 1
        yscale = hydrofuncs.yscale(hydro)
        vmin, vmax = hydrofuncs.hydro_range(hydro)
        if yscale == "log":
            if im.max() != im.min():
                im = np.log10(im) # NOTE: make sure you don't have -ve or zero values here!
        # Make image mesh
        xl, yl = im.shape
        xarr = np.arange(0,zoomedboxlen*1.0000001,zoomedboxlen/(xl-1.0))
        yarr = np.arange(0,zoomedboxlen*1.0000001,zoomedboxlen/(yl-1.0))
        # Plot image map
        if len(ims) == 1:
            finalim = im
            print(hydro)
            print ("ORIGINAL",finalim.min(),finalim.max())
            finalim = hydrofuncs.axis_rescale_func(hydro)(finalim)
            print ("RESCALED",finalim.min(),finalim.max())
        else:
            col = threecolour[ihydro]
            if finalim is None:
                finalim = np.zeros((xl,yl,3))
            if im.min() < im.max():
                imscaled = (im - im.min()) / (im.max() - im.min()) # Normalise values
                finalim[:,:,0] += imscaled*col[0]
                finalim[:,:,1] += imscaled*col[1]
                finalim[:,:,2] += imscaled*col[2]
                finalim[finalim > 1.0] = 1.0

    flipsinks = False
    if len(ims) == 1:
        # imshow and pcolormesh display y axis inverted from each other
        flipsinks = True
        cax = ax.pcolormesh(xarr,yarr,np.flipud(finalim),vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        # Use the RGB map above
        flipsinks = True
        cax = ax.imshow(finalim,vmin=vmin,vmax=vmax,extent=(xarr.min(),xarr.max(),yarr.min(),yarr.max()))
    rdm.AddArray(finalim,label=label+" IMAGE")
        
    # Overlay any contours requested
    for contour in contours:
        if contour == "Wind":
            dmap = rayMap.RayTraceMap(snap,"xrayemission2",los,zoom=zoom) # Returns hot emission
            dum1, dum2, contourim   = dmap.getRaytraceMap()
            contourlims = [1e-37]
            contourcolour = "c"
        if contour == "Ionised":
            dmap = rayMap.RayTraceMap(snap,"xHIImax",los,zoom=zoom) # Returns max xHII
            dum1, dum2, contourim   = dmap.getRaytraceMap()
            contourlims = [1e-2] # Find any xHII above a small value
            contourcolour = "r"
        if contour == "WindSlice":
            dmap = sliceMap.SliceMap(snap,"xrayemission2",los=los,zoom=zoom,starC=starC)
            dum1, dum2, contourim   = dmap.getSliceMap()
            contourlims = [1e-29]
            contourcolour = "c"
        if contour == "IonisedSlice":
            dmap = sliceMap.SliceMap(snap,"xHII",los=los,zoom=zoom,starC=starC)
            dum1, dum2, contourim   = dmap.getSliceMap()
            contourlims = [1e-2]
            contourcolour = "r"
        if contour == "FreeStreamSlice":
            dmap = sliceMap.SliceMap(snap,"spd",los=los,zoom=zoom,starC=starC)
            dum1, dum2, contourim   = dmap.getSliceMap()
            contourlims = [1000]
            contourcolour = "m"
        if contour == "Damkoehler4Slice":
            lbubble = mincellsize*pcincm
            hydrofuncs.allhydros.AddGlobals({"Ldamkoehler":lbubble})
            hydrofuncs.allhydros.AddGlobals({"starpos":starpos})
            dmap = sliceMap.SliceMap(snap,"Damkoehler4",los=los,zoom=zoom,starC=starC)
            dum1, dum2, contourim   = dmap.getSliceMap()
            contourlims = [1]
            contourcolour = "b"
        if contour == "MaskCDSlice":
            dmap = sliceMap.SliceMap(snap,"MaskCD",los=los,zoom=zoom,starC=starC)
            dum1, dum2, contourim   = dmap.getSliceMap()
            contourlims = [1]
            contourcolour = "b"
        contourim = np.flipud(contourim)
        ax.contour(xarr,yarr,contourim,contourlims,colors=contourcolour,alpha=0.75,linewidths=1.5)
        rdm.AddArray(finalim,label=label+" contour"+contour)

    # Draw a box around a region we want to zoom in on
    if zoombox > 0.0:
        #import pdb; pdb.set_trace()
        zbstart = (1.0 - zoombox)/2.0 * zoomedboxlen
        zblength = zoombox * zoomedboxlen
        zbcolour = mplcm.get_cmap(cmap)(0.0)
        p = plt.Rectangle((zbstart, zbstart), zblength, zblength, fill=False,color=zbcolour)
        #p.set_transform(ax.transAxes)
        #p.set_clip_on(False)
        ax.add_patch(p)
        
    # Plot the mesh
    #cmap = "YlGnBu_r"
    #cmap = "RdPu_r"
    cax.set_rasterized(True)
    ax.set_xlim(0,zoomedboxlen)
    ax.set_ylim(0,zoomedboxlen)
   
    if wsink: 
        # Plot sinks
        try:
            lenm = len(np.atleast_1d(sinkm))
        except:
            lenm = 0
        if lenm > 0:
            area = np.pi * sinkm / 50.0
            if flipsinks:
                sinkx = zoomedboxlen - sinkx
            # Draw all sinks
            ax.scatter(sinky,sinkx,s=area,c="w",alpha=0.5,edgecolors='w')
            rdm.AddPoints(sinky,sinkx,label=label+" SINKS")
            # Draw star over main source
            if starsink is not None:
                starsinky = np.atleast_1d(np.array([sinky[starsink]]))
                starsinkx = np.atleast_1d(np.array([sinkx[starsink]]))
                stararea = np.atleast_1d(np.array([area[starsink]]))
                ax.scatter(starsinky,starsinkx,
                           marker="*",s=7*stararea,c="r",alpha=0.5,edgecolors="w")
            rdm.AddArray(np.array([starsinky, starsinkx, starsinky, starsinkx]).reshape(2,2),label=label+"STAR SINK")
            
    # Add scale axis
    scalecol = "w"
    if dolengthscale:
        # length scale in pc (hopefully this is what units boxlen is in)
        lscale = FindLscale(zoomedboxlen)
        x1 = 0.98 * zoomedboxlen
        x2 = x1 - lscale
        textx = x2 - 0.02 * zoomedboxlen
        textalign = "right"
        # HACK
        if not label or Slice:
            x1 = 0.02 * zoomedboxlen
            x2 = x1 + lscale
            textx = x2 + 0.02 * zoomedboxlen
            textalign = "left"
        if not Slice:
            y1 = 0.90 * zoomedboxlen
            y2 = y1
            verticalalignment="center"
        else:
            y1 = 0.04 * zoomedboxlen
            y2 = y1
            verticalalignment ="center"
        line = ax.plot([x1,x2],[y1,y2],"w",path_effects=[PathEffects.withStroke(linewidth=OUTLINEWIDTH+3, foreground='k')])
        txt = ax.text(textx,y2, "  "+str(lscale)+" pc",color=scalecol,
                horizontalalignment=textalign,
                verticalalignment=verticalalignment,fontsize="large")
        txt.set_path_effects([PathEffects.withStroke(linewidth=OUTLINEWIDTH, foreground='k')])
    # Add label
    if label:
        xt = 0.98 * zoomedboxlen
        yt = 0.02 * zoomedboxlen
        txt = ax.text(xt,yt,label,
                horizontalalignment="right",
                verticalalignment="bottom",
                color=scalecol,fontsize="large")
        txt.set_path_effects([PathEffects.withStroke(linewidth=OUTLINEWIDTH, foreground='k')])
    # Add time 
    if plottime:
        xt = 0.02 * zoomedboxlen
        yt = 0.02 * zoomedboxlen
        txt = ax.text(xt,yt,timeL,
                horizontalalignment="left",
                verticalalignment="bottom",
                color=scalecol,fontsize="large")
        txt.set_path_effects([PathEffects.withStroke(linewidth=OUTLINEWIDTH, foreground='k')])
    return cax

def MakeFigure(simnames,times,name,los=None,hydro="rho",Slice=False,wsink=False,starC=False,
               nonamelabel=False,timeL=None,shape=None,dpi=200.0,zoom=1.0,contours=[],forcerun=False,zoombox=-1,
               plotcolorbar=True,doplottime=False,velocitybins=False,altname=None):
    ncols = len(simnames)
    nrows = len(times)

    rdm = rdmfile.RDMFile(__file__)

    if type(hydro) == type("thisisastring"):
        hydros = [hydro]
    else:
        hydros = hydro
        
    # Make figure name
    suffix = ""
    if len(name) > 0:
        suffix = "_"+name
    #wsname = Hamu.SimData.Settings.Settings()["CurrentWorkspace"]
    if len(hydros) > 1:
        hname = ""
        for h in hydros:
            hname += h+"_"
    else:
        hname = hydros[0]    
            
    if Slice:
        folder = "../plots/vis/slice/"
        MakeDirs(folder)
        figname = "sliceTime_"+hname+suffix+".pdf"
    elif velocitybins:
        folder = "../plots/vis/velocitybins/"
        MakeDirs(folder)
        figname = "velocitybinsTime_"+hname+suffix+".pdf"
    else:
        folder = "../plots/vis/multiray/"
        MakeDirs(folder)
        figname = "multirayTime_"+hname+suffix+".pdf"
    # Check if figure needs to be made from figure list?
    figurelist = listfigures.makelist()
    if len(figurelist) > 0:
        if figname not in figurelist and not forcerun:
            print("Figure",figname,"not used by paper; returning without making figure")
            return
    figname = folder+figname
    # Set up figure
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, 
                             sharex=True, sharey=True,
                             frameon=False)
    finches = IMSIZE/dpi
    
    axes = np.atleast_1d(axes)
    
    if len(axes) <= 1:
        finches *= 1.5 # heuristic
    
    # Originally this was run through Hamu, but we don't really need another level of Hamu for otherwise quick functions
    if wsink:
        createColDensMap_sink = _createColDensMap_sink
        createRayMap_sink     = _createRayMap_sink
        createSliceMap_sink   = _createSliceMap_sink
    else:
        createColDensMap   = _createColDensMap
        createRayMap       = _createRayMap
        createSliceMap     = _createSliceMap
        createSliceMap_sink = _createSliceMap_sink

    # Run for all sims
    dolengthscale = False
    isim = 0
    for ax in axes.flatten():
        ax.set_axis_off()
        ax.set_aspect("equal", "datalim")
        ax.axis("off")

    for ii, timetuple in enumerate(times):
        for simname in simnames:
            sim   = Hamu.Simulation(simname)
            # Time stuff
            try:
                snap = sim.Snapshots()[0]
            except:
                print("ERROR FINDING ANY SNAPSHOTS IN", simname)
                print(sim.Folder())
                raise ValueError
            myr   = snap.RawData().info["unit_time"].express(C.Myr)
            time, timeunits = timetuple
            tcreated, sfe = starrelations.runforsim(simname,"firsttime")
            tcreatedcode = tcreated/myr
            if timeunits == "Myr":
                time /= myr
            if timeunits == "MyrFirstStar":
                time += tcreated # start from time first star created
                time /= myr
            if timeunits == "code":
                pass # already ok
            if timeunits == "codeFirstStar":
                time += tcreated
            if timeunits == "outputNumber":
                outsnaps = {snap.OutputNumber():snap for snap in sim.Snapshots()}
                time = outsnaps[int(time)].Time()
            # Simulation stuff
            snap  = sim.FindAtTime(time)
            ax    = axes.flatten()[isim]
            isim += 1
            # Get sink info
            stellar = stellars.FindStellar(snap)
            if len(stellar.mass) == 0:
                return 0.0
            imax = np.where(stellar.mass == np.max(stellar.mass))[0][0]
            sinkid = stellar.sinkid[imax]-1
            sink = sinks.FindSinks(snap)
            boxlen = snap.RawData().info["boxlen"]
            levelmax = snap.RawData().info["levelmax"]
            mincellsize = boxlen / 2.0**levelmax
            smass = stellar.mass
            starsinkid = stellar.sinkid[np.where(smass == smass.max())]
            starsink = np.where(sink.id == starsinkid)[0]
            starpos = np.array([sink.x[starsinkid],sink.y[starsinkid],sink.z[starsinkid]])/boxlen
            # Get global variables if needed
            if "Damkoehler" in hydro:
                #lbubblefunc = Hamu.Algorithm(findproperties.maxradiusatstarpos)
                #lbubble = lbubblefunc(snap)*pcincm
                # INSTEAD use cell size for outer scale (since we already include bubble-scale turbulence)
                lbubble = mincellsize*pcincm
                hydrofuncs.allhydros.AddGlobals({"Ldamkoehler":lbubble})
                hydrofuncs.allhydros.AddGlobals({"starpos":starpos})
            # One hydro variable?
            if type(hydro) == type("rho"):
                dohydrolist = False
                cmap  = linestyles.ColourMap(simname, hydro)
            # A list of them?
            else:
                dohydrolist = True
                cmap = None
                
            def MakeData(hydro):
                # Make sure we don't try to load the CD mask if it's not there
                if not "_CDMASK" in simname:
                    if hydro == "xrayemission3":
                        # Use function that doesn't include MaskCD
                        hydro = "xrayemission2"
                if Slice:
                    if wsink:
                        data  = createSliceMap_sink(snap,hydro,los,zoom,starC)      
                    else:
                        data  = createSliceMap(snap,hydro,los,zoom,starC)     
                else: 
                    #if hydro == 'NH':
                    #    if wsink: 
                    #        data  = createColDensMap_sink(snap,los,hydro,wsink,zoom,starC)      
                    #    else: 
                    #        data  = createColDensMap(snap,los,hydro,wsink,zoom,starC)       
                    #else:
                    if wsink: 
                        data = createRayMap_sink(snap,los,hydro,wsink,zoom,starC)    
                    else: 
                        data = createRayMap(snap,los,hydro,wsink,zoom,starC)   
                return data
            datas = [MakeData(hydro) for hydro in hydros]
                            
            if (simname == simnames[-1] and ii == len(times)-1):
                dolengthscale = True
            dolengthscale = True # For referee
            # Utter hack
            if simname == "UVWINDPRESS_30" and ii == 0 and hydro == "NH":
                dolengthscale = False
            if (simname == simnames[0]) and len(axes) > 1:
                plottime = True
            if not doplottime:
                plottime = False
            label = ""
            # Only plot name label for the first image in a simulation
            if not nonamelabel and ii == 0:
                label =  linestyles.Label(simname)
            #if len(axes) == 1:
            #    label = None
            # Make the pyplot image axis object
            im    = MakeImage(datas,hydros,snap,wsink,ax,dolengthscale,cmap,
                              plottime, timeL[ii],label = label,starsink=starsink,rdm=rdm,
                              contours=contours,zoombox=zoombox/zoom,starC=starC,zoom=zoom,Slice=Slice)
            plottime      = False
            dolengthscale = False

    # Add colour bar for hydro variable
    print("Making colour bar...")
    # ---- This cax works well with nrow = 1, may need adjusting, could be improved 
    # Colorbar at the bottom of the plots
    #cax  = fig.add_axes([0.2, -0.022, 0.4, 0.02])
    # Colorbar at the top of all the plots
    if len(hydros) == 1 and plotcolorbar:
        hydro = hydros[0]
        #if label:
        cax  = fig.add_axes([0.24, 0.98, 0.72, 0.02])
        #else:
        #    cax  = fig.add_axes([0.24, 0.1, 0.72, 0.02])
        cbar = fig.colorbar(im,cax=cax,orientation="horizontal")
        label = hydrofuncs.hydro_label((hydro))
        yscale = hydrofuncs.yscale(hydro)
        if yscale == "log":
            label = "log("+label+")"
        #cbar.set_label(label,fontsize="medium",color="k")
        #cbar.ax.tick_params(labelsize="medium",labelcolor="k")
        #cbar.solids.set_edgecolor("face")
        cbar.set_label(label,fontsize="large",color="w",
                       path_effects=[PathEffects.withStroke(linewidth=OUTLINEWIDTH, foreground='k')])
        cbar.ax.tick_params(labelsize="large",labelcolor="w")
        for t in cbar.ax.get_xticklabels()+cbar.ax.get_yticklabels():
            t.set_path_effects([PathEffects.withStroke(linewidth=OUTLINEWIDTH, foreground='k')])
        cbar.solids.set_edgecolor("face")
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')

        
    # Save figure
    print("Saving figure "+figname+"...")
    fig.subplots_adjust(hspace=0.0, wspace=0.0, 
                        left=0.20,right=1.0,
                        bottom=0.00,top=1.0)
    #fig.set_size_inches(finches*ncols/0.8,finches*(nrows+1))
    fig.set_size_inches(finches*ncols/0.8,finches*(nrows))
    fig.savefig(figname,
                pad_inches=0,
                dpi=dpi)
    if altname is not None:
        fig.savefig(altname,
                pad_inches=0,
                dpi=dpi)
    rdm.Write(figname)
    # Crop out borders
    os.system("pdfcrop "+figname+" "+figname)
    print("Done!")

def hotchampagneplot():
    # Slices
    linestyles.CURRSIMSET = "hotchampagne"
    simname = "SEED2_35MSUN_CDMASK_WINDUV"
    sim = hamusims[simname]
    # Choose outputs to plot
    outnums = [52,53,54,55]
    timetuples = [(o,"outputNumber") for o in outnums]
    # Get times in Myr
    outsnaps = {snap.OutputNumber():snap for snap in sim.Snapshots()}
    myr   = outsnaps[outnums[0]].RawData().info["unit_time"].express(C.Myr)
    tcreated = timefuncs.FindTcreatedFirstStar(sim)
    timesMyr = [outsnaps[o].Time() * myr - tcreated for o in outnums]
    timeL = ["%.2f" % x + r' Myr' for x in timesMyr]
    # Set up rest of plot
    los = "x"
    zoom = 0.25
    newsetname = "hotchampagne_zoom"+str(zoom)+"_"
    newsetname = newsetname.replace(".","p") # the extra dot confuses latex
    figname = newsetname+"_"+los
    # Make slices
    for hydro in ["vorticity2px_timescale","T","rho","xHII","P"]:
        MakeFigure([simname],timetuples,name=figname,los=los,hydro=hydro,
                    Slice=True,wsink=True,starC=True,doplottime=True,
                    timeL=timeL,zoom=zoom,forcerun=True)

if __name__=="__main__":

    # Should we force some figures to run?
    forcerun=True
    
    # Yes I did intend for this to sound weird
    hotchampagneplot()
    
    for setname, simset in simsets.items():
        linestyles.CURRSIMSET = setname
        #simset = ["NOFB","UV_"+smass,"UVWINDPRESS_"+smass]
        #setname = "windset_"+smass+"Msun"
        #simwindname = "UVWIND_"+smass
        #times = np.array([0.5, 0.75, 1.])
#        times = np.array([0.2]) # np.array([0.9]) # [0.9] # 3.5 Myr = tstarformed + 0.2 Myr
        times = np.array([0.3]) # np.array([0.9]) # [0.9] # 3.5 Myr = tstarformed + 0.2 Myr 
        zoom = 0.5
        #if dense:
        #    zoom = 1.0
        newsetname = setname+str(times[-1])+"Myr_"+"zoom"+str(zoom)+"_"
        newsetname = newsetname.replace(".","p") # the extra dot confuses latex
        #timeL = [str(x)+r' t$_{ff}$' for x in times]
        #timesin = [(time*tffcloud_code,"code") for time in times]
        timeL = [str(x)+r' Myr' for x in times]
        timesin = [(time,"MyrFirstStar") for time in times]
        for los in "yxz":
            figname = newsetname+"_"+los
            zoom2 = 0.25
            figname2 = figname.replace("zoom"+str(zoom).replace(".","p"),
                                        "zoom"+str(zoom2).replace(".","p"),)
            # Run for movie
            imovie = 0
            tmovies = []#np.linspace(0.1,1.0,50)
            for tmovie in tmovies:
                imovie += 1
                timein = (tmovie,"MyrFirstStar")
                tmovieL = str(tmovie)+r' Myr'
                coolhydros = ["coolemission","ionemission4","xrayemission2"]
                movienum = str(imovie).zfill(3)
                MakeFigure(simset,[timein],name=figname+"movie"+movienum,los=los,
                            hydro=coolhydros,Slice=False,wsink=True,
                            timeL=[tmovieL],zoom=zoom,forcerun=forcerun)
                #for hydro in ["rho","T","Lcool"]:
                #    MakeFigure([simset[-1]],[timein],name=figname+"movieslice",los=los,hydro=hydro,
                #               Slice=True,wsink=True,starC=True,
                #               timeL=[tmovieL],zoom=zoom,forcerun=True)
            # Merged emission map - just wind
            coolhydros = ["coolemission","ionemission4","xrayemission3"]
            timesmerged = [0.1,0.2,0.3]
            timesmergedIn = [(time,"MyrFirstStar") for time in timesmerged]
            timesmergedL = [str(x)+r' Myr' for x in timesmerged]
            # Doesn't really work, just stitch together sequences with a script
            #if mass == 30:
            #    windsimnames = ["UVWIND_"+x for x in ["30","60","120","120_DENSE"]]
            #    MakeFigure(windsimnames,timesmergedIn,name=allfigname+"windonly_sequence",los=los,hydro=coolhydros,Slice=False,wsink=True,
            #            timeL=timesmergedL,zoom=zoom,forcerun=True,doplottime=True)
            #    presssimnames = ["UVWINDPRESS_"+x for x in ["30","60","120","120_DENSE"]]
            #    MakeFigure(presssimnames,timesmergedIn,name=allfigname+"windpressonly_sequence",los=los,hydro=coolhydros,Slice=False,wsink=True,
            #            timeL=timesmergedL,zoom=zoom,forcerun=True,doplottime=True)

            # Forcing this figure to run
            #MakeFigure([simset[-1]],[timesin[-1]],name=figname+"ERF",los=los,hydro=coolhydros,
            #           Slice=False,wsink=True,starC=True,nonamelabel=True,
            #           timeL=[timeL[-1]],zoom=1.0,forcerun=True)
            #MakeFigure([simset[-1]],[timesin[-1]],name=figname+"ERF",los=los,hydro="T",
            #           Slice=True,wsink=True,starC=True,nonamelabel=True,
            #           timeL=[timeL[-1]],zoom=1.0,forcerun=True)
            #MakeFigure([simset[-1]],[timesin[-1]],name=figname+"ERF",los=los,hydro="NH",
            #           Slice=False,wsink=True,starC=True,nonamelabel=True,
            #           timeL=[timeL[-1]],zoom=1.0,forcerun=True,contours=["Wind"])

            DEBUG = False


            # Damkoehler comparison plot
            '''
            if setname == "single":
                for hydro in ["T","rho"]:
                    #contourslist = [[],["Damkoehler4Slice"]]
                    contourslist = [["MaskCDSlice"],[]]
                    Damzoom = 0.4
                    for contours in contourslist:
                        contoursname = ""
                        for contour in contours:
                            contoursname += contour
                        MakeFigure(simset,timesmergedIn,
                                   name=figname.replace(str(zoom).replace(".","p"),
                                                        str(Damzoom).replace(".","p")+"singleslice"+contoursname),
                                   los=los,hydro=hydro,Slice=True,wsink=True,starC=True,
                                   timeL=timesmergedL,zoom=Damzoom,forcerun=True,contours=contours)
            '''

        

            # Emission and NH maps
            for hydro in [coolhydros,"NH"]:
                contourslist = [[],["Wind","Ionised"]]
                for contours in contourslist:
                    contxt = ""
                    if len(contours) > 0:
                        contxt = "_contours"
                        MakeFigure(simset,timesmergedIn,name=figname+"windpressonly_sequence"+contxt,
                               los=los,hydro=hydro,Slice=False,wsink=True,
                               timeL=timesmergedL,zoom=zoom,forcerun=True,
                               doplottime=True,contours=contours,
                               plotcolorbar=True)

            # Single slices
            #for hydro in ["vorticity1px_timescale","vorticity1px_speedcompare"]:
            #    MakeFigure([simset[0]],[timesin[-1]],name=figname+"singleslice",los=los,hydro=hydro,
            #                Slice=True,wsink=True,starC=True,
            #                timeL=[timeL[-1]],zoom=zoom,forcerun=True)
            
            if setname == "single":
                for z in [zoom,zoom2]:
                    for hydro in ["T","rho","spd",
                                  "vdispersion1px","vdispersion1px_speedcompare",
                                  "vdispersion2px","vdispersion2px_speedcompare",
                                  "vorticity1px_timescale",
                                  "vorticity2px_timescale",
                                  "vorticity4px_timescale",
                                  "vorticity1px_speedcompare","vorticity2px_speedcompare","vorticity4px_speedcompare",
                                  "vorticity1px","vorticity2px","vorticity4px"]:
                        #"Lcool","T","rho","xHII","xHeII","xHeIII",
                        #"P","vradfrac3","vrad","vx","vy","vz"]:
                        MakeFigure([simset[0]],[timesin[-1]],name=figname+"singleslice",los=los,hydro=hydro,
                                   Slice=True,wsink=True,starC=True,
                                   timeL=[timeL[-1]],zoom=z,forcerun=True)

            # Slices
            for hydro in ["vorticity2px_timescale","Lcool","T","rho","xHII","xHeII","xHeIII","P"]:
                MakeFigure(simset,[timesin[-1]],name=figname,los=los,hydro=hydro,
                           Slice=True,wsink=True,starC=True,
                           timeL=[timeL[-1]],zoom=zoom,forcerun=True)
                MakeFigure(simset,timesmergedIn,name=figname,los=los,hydro=hydro,
                           Slice=True,wsink=True,starC=True,
                           timeL=timesmergedL,zoom=zoom,forcerun=True)
                                

            # Separate emission maps and related images
            for hydro in ["ionemission4","xrayemission2","coolemission","NH","xHIImax","fastmass6"][::-1]:
                MakeFigure(simset,[timesin[-1]],name=figname,los=los,hydro=hydro,Slice=False,wsink=True,
                            timeL=[timeL[-1]],zoom=zoom)

            for hydro in ["EkinperEtherm"]:
                #["Ekin","Etherm","EkinperEtherm","xrayemission2"]:
                MakeFigure(simset,[timesin[-1]],name=figname,los=los,hydro=hydro,
                            Slice=True,wsink=True,starC=True,
                            timeL=[timeL[-1]],zoom=zoom,forcerun=forcerun,
                            contours=["WindSlice","IonisedSlice","FreeStreamSlice"])


            #if DEBUG:
            #    MakeFigure([simwindname],timesmergedIn,name=figname+"windonly_sequence",los=los,hydro=coolhydros,Slice=False,wsink=True,
            #           timeL=timesmergedL,zoom=zoom,forcerun=forcerun,doplottime=True)
            #    MakeFigure([simset[-1]],[timesin[-1]],name=figname+"windonly",los=los,hydro=coolhydros,Slice=False,wsink=True,
            #                timeL=[timeL[-1]],zoom=zoom,forcerun=forcerun)
            #                     - All physics
            #if DEBUG:
            #    MakeFigure(simset,[timesin[-1]],name=figname+"allphysics",los=los,hydro=coolhydros,Slice=False,wsink=True,
            #                timeL=[timeL[-1]],zoom=zoom)

                    #MakeFigure([simset[-1]],[timesin[-1]],name=figname2+"windonly",los=los,hydro=hydro,Slice=True,wsink=True,starC=True,timeL=[timeL[-1]],zoom=zoom2)

            # Plot for Raphael
            #for z, zb in ([0.25,0.025],[0.025,0.004]):
            #    MakeFigure([simset[0]],[(0.1,"MyrFirstStar")],
            #               name=figname+"_hpcproposal_nofb"+str(z).replace(".","p").replace("zoom"+str(zoom).replace(".","p"),"")#,
#                               los=los,hydro='rho',Slice=False,wsink=False,nonamelabel=True,
#                               timeL=['0 Myr'],zoom=z,starC=True,forcerun=True,zoombox=zb

            #if DEBUG:
            #    # Column density
            #    MakeFigure([simset[-1]],[timesin[-1]],name=figname+"windonly",
            #                los=los,hydro='NH',Slice=False,wsink=True,
            #                timeL=[timeL[-1]],zoom=zoom,contours=["Wind","Ionised"],
            #                plotcolorbar=True)
            #    # Column density (all sims)
            #    MakeFigure(simset,timesin,name=figname,
            #                los=los,hydro='NH',Slice=False,wsink=True,
            #                timeL=timeL,zoom=zoom,contours=["Wind","Ionised"],
            #                plotcolorbar=(mass==30))
            #MakeFigure(simset,timesin,name=figname,los=los,hydro='maxT',Slice=False,wsink=True,timeL=timeL,zoom=zoom)

            # Temperature slice (all sims)
            MakeFigure(simset,timesin,name=figname,los=los,hydro='T',Slice=True,wsink=True,timeL=timeL,zoom=zoom,starC=True,forcerun=forcerun)
        
