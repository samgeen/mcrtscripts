'''
Make images of the cloud/sinks
Rebekka Bieri and Sam Geen, August 2018
'''

from startup import *

from pymses.utils import constants as C

import columndensity, rayMap, sliceMap, sinks, ysos, starrelations, listfigures

from matplotlib import rc
import matplotlib.cm as mplcm

import rdmfile

rc('axes', labelsize=8)
rc('axes', labelsize=8,linewidth=0.5)
rc('xtick', labelsize=8)
rc('ytick', labelsize=8)

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

def _createColDensMap(snap,los=None,hydro='rho',wsink=False,zoom=1.0):
    dmap = columndensity.DensityMap(snap,los,zoom=zoom)
    im = dmap.NH()
    boxlen = snap.info["boxlen"]
    return [im, boxlen*zoom]

def _createColDensMap_sink(snap,los=None,hydro='rho',wsink=True,zoom=1.0):
    ro = snap.RawData()
    dmap = columndensity.DensityMap(snap,los,zoom=zoom)
    im = dmap.NH()
    boxlen = ro.info["boxlen"]
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(ro,los)
    # Shift sink position if zooming in
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
    im   = dmap.getRaytraceMap()
    boxlen = ro.info["boxlen"]
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(ro,los)
    # Shift sink position if zooming in
    sinkx -= 0.5*(1.0-zoom)*boxlen
    sinky -= 0.5*(1.0-zoom)*boxlen
    return [im, sinkx, sinky, sinkm, boxlen*zoom]

def _createSliceMap(snap,hydro='rho',los=None,zoom=1.0,starC=False):
    ro = snap.RawData()
    dmap = sliceMap.SliceMap(snap,hydro,los,zoom=zoom,starC=starC)
    im   = dmap.getSliceMap()
    boxlen = ro.info["boxlen"] 
    return [im, boxlen*zoom]

def _createSliceMap_sink(snap,hydro='rho',los=None,zoom=1.0,starC=False):
    ro = snap.RawData()
    dmap = sliceMap.SliceMap(snap,hydro,los,zoom=zoom,starC=starC)
    im   = dmap.getSliceMap()
    boxlen = ro.info["boxlen"] 
    # Plot with sinks
    sinkx, sinky, sinkm = ProjectSinks(ro,los)
    # Shift sink position if zooming in
    sinkx -= 0.5*(1.0-zoom)*boxlen
    sinky -= 0.5*(1.0-zoom)*boxlen
    return [im, sinkx, sinky, sinkm, boxlen*zoom]

def rgb(r,g,b):
    return (float(r) / 255.0, float(g) / 255.0, float(b)/255.0)

def MakeImage(datas,hydros,snap,wsink,ax,dolengthscale,cmap,plottime=False,timeL=None,label=None,
                  starsink=None,rdm=None,contours=[],zoombox=-1,starC=False,zoom=1.0):

    ims = []
    for data, hydro in zip(datas, hydros):
        # NOTE: for different hydro variables, only im should be different
        if wsink:
            im, sinkx, sinky, sinkm, boxlen = data
        else:
            im, boxlen = data

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
                   rgb(117,112,179)]
    for im, hydro in zip(ims, hydros):
        ihydro += 1
        yscale = hydrofuncs.yscale(hydro)
        vmin, vmax = hydrofuncs.hydro_range(hydro)
        if yscale == "log":
            if im.max() != im.min():
                im = np.log10(im) # NOTE: make sure you don't have -ve or zero values here!
        # Make image mesh
        xl, yl = im.shape
        xarr = np.arange(0,boxlen*1.0000001,boxlen/(xl-1.0))
        yarr = np.arange(0,boxlen*1.0000001,boxlen/(yl-1.0))
        # Plot image map
        if len(ims) == 1:
            finalim = im
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
            contourim   = dmap.getRaytraceMap()
            contourlims = [1e-37]
            contourcolour = "r"
        if contour == "Ionised":
            dmap = rayMap.RayTraceMap(snap,"xHIImax",los,zoom=zoom) # Returns max xHII
            contourim   = dmap.getRaytraceMap()
            contourlims = [1e-2] # Find any xHII above a small value
            contourcolour = "c"
        if contour == "WindSlice":
            dmap = sliceMap.SliceMap(snap,"xrayemission2",los=los,zoom=zoom,starC=starC)
            contourim   = dmap.getSliceMap()
            contourlims = [1e-29]
            contourcolour = "r"
        if contour == "IonisedSlice":
            dmap = sliceMap.SliceMap(snap,"xHII",los=los,zoom=zoom,starC=starC)
            contourim   = dmap.getSliceMap()
            contourlims = [1e-2]
            contourcolour = "c"
        if contour == "FreeStreamSlice":
            dmap = sliceMap.SliceMap(snap,"spd",los=los,zoom=zoom,starC=starC)
            contourim   = dmap.getSliceMap()
            contourlims = [1000]
            contourcolour = "m"
        contourim = np.flipud(contourim)
        ax.contour(xarr,yarr,contourim,contourlims,colors=contourcolour,alpha=0.75,linewidths=1.5)
        rdm.AddArray(finalim,label=label+" contour"+contour)

    # Draw a box around a region we want to zoom in on
    if zoombox > 0.0:
        #import pdb; pdb.set_trace()
        zbstart = (1.0 - zoombox)/2.0 * boxlen
        zblength = zoombox * boxlen
        zbcolour = mplcm.get_cmap(cmap)(0.0)
        p = plt.Rectangle((zbstart, zbstart), zblength, zblength, fill=False,color=zbcolour)
        #p.set_transform(ax.transAxes)
        #p.set_clip_on(False)
        ax.add_patch(p)
        
    # Plot the mesh
    #cmap = "YlGnBu_r"
    #cmap = "RdPu_r"
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
            area = np.pi * sinkm / 50.0
            if flipsinks:
                sinkx = boxlen - sinkx
            # Draw all sinks
            ax.scatter(sinky,sinkx,s=area,c="w",alpha=0.5,edgecolors='w')
            rdm.AddPoints(sinky,sinkx,label=label+" SINKS")
            # Draw star over main source
            if starsink is not None:
                starsinky = [sinky[starsink]]
                starsinkx = [sinkx[starsink]]
                ax.scatter(starsinky,starsinkx,
                           marker="*",s=7*area,c="r",alpha=0.5,edgecolors="w")
            rdm.AddPoints(starsinky, starsinkx,label=label+"STAR SINK")
            
    # Add scale axis
    scalecol = "w"
    if dolengthscale:
        # length scale in pc (hopefully this is what units boxlen is in)
        lscale = FindLscale(boxlen)
        x1 = 0.1 * boxlen 
        x2 = x1 + lscale
        y1 = 0.1 * boxlen
        y2 = y1
        ax.plot([x1,x2],[y1,y2],scalecol)
        ax.text(x2,y2, "  "+str(lscale)+" pc",color=scalecol,
                verticalalignment="center",fontsize="large")
    # Add label
    if label:
        xt = 0.9 * boxlen
        yt = 0.1 * boxlen
        ax.text(xt,yt,label,
                horizontalalignment="right",
                verticalalignment="center",
                color=scalecol,fontsize="large")
    # Add time 
    if plottime:
        xt = 0.2 * boxlen
        yt = 0.9 * boxlen
        ax.text(xt,yt,timeL,
                horizontalalignment="right",
                verticalalignment="center",
                color=scalecol,fontsize="large")
    return cax

def MakeFigure(simnames,times,name,los=None,hydro="rho",Slice=False,wsink=False,starC=False,
               nonamelabel=False,timeL=None,shape=None,dpi=200.0,zoom=1.0,contours=[],forcerun=False,zoombox=-1,
               plotcolorbar=True):
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
    else:
        folder = "../plots/vis/multiray/"
        MakeDirs(folder)
        figname = "multirayTime_"+hname+suffix+".pdf"
    # Check if figure needs to be made from figure list?
    figurelist = listfigures.makelist()
    if len(figurelist) > 0:
        if figname not in figurelist and not forcerun:
            print "Figure",figname,"not used by paper; returning without making figure"
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
    doplottime      = False
    isim = -1
    for ax in axes.flatten():
        ax.set_axis_off()
        ax.set_aspect("equal", "datalim")
        ax.axis("off")

    for ii, timetuple in enumerate(times):
        for simname in simnames:
            sim   = Hamu.Simulation(simname)
            # Time stuff
            snap = sim.Snapshots()[0]
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
            # Simulation stuff
            snap  = sim.FindAtTime(time)
            ax    = axes.flatten()[isim]
            isim += 1
            # One hydro variable?
            if type(hydro) == type("rho"):
                dohydrolist = False
                cmap  = linestyles.ColourMap(simname, hydro)
            # A list of them?
            else:
                dohydrolist = True
                cmap = None
                
            def MakeData(hydro):
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
                            
            if (simname == simnames[-1] and ii == 0):
                dolengthscale = True 
            if (simname == simnames[0]) and len(axes) > 1:
                plottime = True
            if not doplottime:
                plottime = False
            label = ""
            if not nonamelabel:
                label =  linestyles.Label(simname)
            #if len(axes) == 1:
            #    label = None
            # Make the pyplot image axis object
            stellar = stellars.FindStellar(snap)
            smass = stellar.mass
            starsinkid = stellar.sinkid[np.where(smass == smass.max())]
            sink = sinks.FindSinks(snap)
            starsink = np.where(sink.id == starsinkid)[0]
            im    = MakeImage(datas,hydros,snap,wsink,ax,dolengthscale,cmap,
                              plottime, timeL[ii],label = label,starsink=starsink,rdm=rdm,
                              contours=contours,zoombox=zoombox/zoom,starC=starC,zoom=zoom)
            plottime      = False
            dolengthscale = False

    # Add colour bar for hydro variable
    print "Making colour bar..."
    # ---- This cax works well with nrow = 1, may need adjusting, could be improved 
    # Colorbar at the bottom of the plots
    #cax  = fig.add_axes([0.2, -0.022, 0.4, 0.02])
    # Colorbar at the top of all the plots
    if len(hydros) == 1 and plotcolorbar:
        hydro = hydros[0]
        cax  = fig.add_axes([0.4, 0.98, 0.4, 0.02])
        cbar = fig.colorbar(im,cax=cax,orientation="horizontal")
        label = hydrofuncs.hydro_label((hydro))
        if not "xH" in hydro:
            label = "log("+label+")"
        #cbar.set_label(label,fontsize="medium",color="k")
        #cbar.ax.tick_params(labelsize="medium",labelcolor="k")
        #cbar.solids.set_edgecolor("face")
        cbar.set_label(label,fontsize="large",color="w")
        cbar.ax.tick_params(labelsize="large",labelcolor="w")
        cbar.solids.set_edgecolor("face")
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')

    # Save figure
    print "Saving figure "+figname+"..."
    fig.subplots_adjust(hspace=0.0, wspace=0.0, 
                        left=0.20,right=1.0,
                        bottom=0.00,top=1.0)
    #fig.set_size_inches(finches*ncols/0.8,finches*(nrows+1))
    fig.set_size_inches(finches*ncols/0.8,finches*(nrows))
    fig.savefig(figname,
                pad_inches=0,
                dpi=dpi)
    rdm.Write(figname)
    # Crop out borders
    os.system("pdfcrop "+figname+" "+figname)
    print "Done!"

if __name__=="__main__":

    for dense in [False]:
        for mass in [120]:
            smass = str(mass)
            if not dense or mass == 120:
                simset = ["NOFB","UV_"+smass,"UVWIND_"+smass]
                if not dense and mass == 120:
                    simset += ["UVWINDCR25_120"]
                    simset += ["UVWINDCR26_120"]
                setname = "windset_"+smass+"Msun"
                if dense:
                    simset = [x+"_DENSE" for x in simset]
                    setname += "_dense"
                #times = np.array([0.5, 0.75, 1.])
                times = np.array([0.2]) # np.array([0.9]) # [0.9] # 3.5 Myr = tstarformed + 0.2 Myr 
                zoom = 0.25
                if dense:
                    zoom = 1.0
                setname = setname+str(times[-1])+"Myr_"+"zoom"+str(zoom)+"_"
                setname = setname.replace(".","p") # the extra dot confuses latex
                #timeL = [str(x)+r' t$_{ff}$' for x in times]
                #timesin = [(time*tffcloud_code,"code") for time in times]
                timeL = [str(x)+r' Myr' for x in times]
                timesin = [(time,"MyrFirstStar") for time in times]
                for los in "xyz":
                    figname = setname+"_"+los
                    zoom2 = 0.5
                    if dense:
                        zoom2 = 0.5
                    figname2 = figname.replace("zoom"+str(zoom).replace(".","p"),
                                               "zoom"+str(zoom2).replace(".","p"),)
                    # Run for movie
                    imovie = 0
                    tmovies = np.linspace(0.1,1.0,50)
                    for tmovie in tmovies:
                        imovie += 1
                        timein = (tmovie,"MyrFirstStar")
                        tmovieL = str(tmovie)+r' Myr'
                        coolhydros = ["coolemission","ionemission","xrayemission2"]
                        movienum = str(imovie).zfill(3)
                        MakeFigure([simset[-1]],[timein],name=figname+"movie"+movienum,los=los,
                                   hydro=coolhydros,Slice=False,wsink=True,
                                   timeL=[tmovieL],zoom=zoom,forcerun=True)
                        #for hydro in ["rho","T","Lcool"]:
                        #    MakeFigure([simset[-1]],[timein],name=figname+"movieslice",los=los,hydro=hydro,
                        #               Slice=True,wsink=True,starC=True,
                        #               timeL=[tmovieL],zoom=zoom,forcerun=True)
                    # Slices
                    for hydro in ["rho","Ekin","Etherm","EkinperEtherm","xrayemission2"]:
                        MakeFigure([simset[-1]],[timesin[-1]],name=figname2+"windonly",los=los,hydro=hydro,
                                   Slice=True,wsink=True,starC=True,
                                   timeL=[timeL[-1]],zoom=zoom2,
                                   contours=["WindSlice","IonisedSlice","FreeStreamSlice"])
                    for hydro in ["Lcool","T","rho","xHII","P"]:
                        MakeFigure([simset[-1]],[timesin[-1]],name=figname2+"windonly",los=los,hydro=hydro,
                                   Slice=True,wsink=True,starC=True,
                                   timeL=[timeL[-1]],zoom=zoom2)
                        MakeFigure([simset[-1]],[timesin[-1]],name=figname+"windonly",los=los,hydro=hydro,
                                   Slice=True,wsink=True,starC=True,
                                   timeL=[timeL[-1]],zoom=zoom)
                        MakeFigure(simset,[timesin[-1]],name=figname+"allphysics",los=los,hydro=hydro,
                                   Slice=True,wsink=True,starC=True,
                                   timeL=[timeL[-1]],zoom=zoom)
                        #MakeFigure([simset[-1]],[timesin[-1]],name=figname2+"windonly",los=los,hydro=hydro,Slice=True,wsink=True,starC=True,timeL=[timeL[-1]],zoom=zoom2)

                    # Plot for Raphael
                    #for z, zb in ([0.25,0.025],[0.025,0.004]):
                    #    MakeFigure([simset[0]],[(0.1,"MyrFirstStar")],
                    #               name=figname+"_hpcproposal_nofb"+str(z).replace(".","p").replace("zoom"+str(zoom).replace(".","p"),"")#,
    #                               los=los,hydro='rho',Slice=False,wsink=False,nonamelabel=True,
    #                               timeL=['0 Myr'],zoom=z,starC=True,forcerun=True,zoombox=zb)
                    # Column density
                    MakeFigure([simset[-1]],[timesin[-1]],name=figname+"windonly",
                               los=los,hydro='NH',Slice=False,wsink=True,
                               timeL=[timeL[-1]],zoom=zoom,contours=["Wind","Ionised"],
                               plotcolorbar=True)
                    # Column density (all sims)
                    MakeFigure(simset,timesin,name=figname,
                               los=los,hydro='NH',Slice=False,wsink=True,
                               timeL=timeL,zoom=zoom,contours=["Wind","Ionised"],
                               plotcolorbar=(mass==30))
                    #MakeFigure(simset,timesin,name=figname,los=los,hydro='maxT',Slice=False,wsink=True,timeL=timeL,zoom=zoom)

                    # Merged emission map - just wind
                    coolhydros = ["coolemission","ionemission","xrayemission2"]
                    MakeFigure([simset[-1]],[timesin[-1]],name=figname+"windonly",los=los,hydro=coolhydros,Slice=False,wsink=True,
                               timeL=[timeL[-1]],zoom=zoom)
                    #                     - All physics
                    MakeFigure(simset,[timesin[-1]],name=figname+"allphysics",los=los,hydro=coolhydros,Slice=False,wsink=True,
                               timeL=[timeL[-1]],zoom=zoom)
                    # Separate emission maps
                    for hydro in ["ionemission","xrayemission2","coolemission"][::-1]:
                        MakeFigure([simset[-1]],[timesin[-1]],name=figname+"windonly",los=los,hydro=hydro,Slice=False,wsink=True,
                                   timeL=[timeL[-1]],zoom=zoom)
                    # Temperature slice (all sims)
                    MakeFigure(simset,timesin,name=figname,los=los,hydro='T',Slice=True,wsink=True,timeL=timeL,zoom=zoom,starC=True,forcerun=True)
        
