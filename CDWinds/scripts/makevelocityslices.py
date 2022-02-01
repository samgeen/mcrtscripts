
from startup import *

from pymses.utils import constants as C

import makeImageGridTime

if __name__=="__main__":

    # Should we force some figures to run?
    forcerun=True

    setname = ["seed1only"]
    simset = ["SEED1_35MSUN_CDMASK_WINDUV"]

    zoom = 0.5
    time = 0.2 # Myr after star born
    setname = setname+str(time)+"Myr_"+"zoom"+str(zoom)+"_"
    setname = setname.replace(".","p") # the extra dot confuses latex
    timeL = [str(x)+r' Myr' for x in times]
    timesin = [(time,"MyrFirstStar") for time in times]
    for los in "y":
        figname = setname+"_"+los
        # Run for each velocity slice
        vdiff = 0.3
        for vmin in np.arange(-100,100.1) * vdiff:
            vmax = vmin + vdiff
            # Sample velocities within the velociyt range required
            def _velocityfunc(ro):
                def maskfunc(dset):
                    NH = dset["rho"]*ro.info["unit_density"].express(C.H_cc)*ro.info["unit_length"].express(C.cm)
                    uvel = ro.info["unit_velocity"].express(C.km/C.s)
                    vels = dset["vel"] * uvel
                    velmaps = {}
                    velmaps["x"] = vels[:,:,0]
                    velmaps["y"] = vels[:,:,1]
                    velmaps["z"] = vels[:,:,2]
                    vellos = velmaps[los]
                    # Cut away things outside the velocity bins specified
                    NH[vellos < vmin] = NH[vellos < vmin]*0.0
                    NH[vellos > vmax] = NH[vellos > vmax]*0.0
                    return NH
                return maskfunc
            velocityfunc = hydrofuncs.Hydro("$N_H("+str(vmin)+"< v/\mathrm{km/s} <"+str(vmax)_")$",_velocityfunc,
                                            ["rho","vel"] ,
                                            "GnBu_r","log",(20,23.7),surfacequantity=True)
            hydro = "velocityslice"+str(vmin)+"_v"+los
            hydrofuncs.allhydros[hydro] = velocityfunc
            MakeFigure(simset,timesin,name=figname,los=los,hydro=hydro,Slice=False,wsink=True,
                       timeL=timeL,zoom=zoom,velocitybins=True)
