'''
Make images of the cloud/sinks
Rebekka Bieri and Sam Geen, August 2018
'''

from startup import *

import makeImageGridTime

def runslices():
    for sim in allsims:
        for outtime in [0.1,0.15,0.2,0.4]:
            #simset = allsims
            simset = ["SHELL_CDMASK3"] 
            setname = "SHELLTEST"
            times = np.array([outtime])
            zoom = 0.5
            setname = setname+str(times[-1])+"Myr_"+"zoom"+str(zoom)+"_"
            setname = setname.replace(".","p") # the extra dot confuses latex
            timeL = [str(x)+r' Myr' for x in times]
            timesin = [(time,"MyrFirstStar") for time in times]
            for los in "y":
                figname = setname+"_"+los                    
                # Merged emission map - just wind
                coolhydros = ["coolemission","ionemission","xrayemission2"]
                #timesmerged = [0.1,0.2,0.3,0.4]
                #timesmergedIn = [(time,"MyrFirstStar") for time in timesmerged]
                #timesmergedL = [str(x)+r' Myr' for x in timesmerged]

                # Emission maps
                IMSIZE = 2048
                for hydro in [coolhydros]:
                    for sim in simset:
                        makeImageGridTime.MakeFigure([simset[-1]],[timesin[-1]],name=figname+"emission",
                                                     los=los,hydro=hydro,Slice=False,wsink=False,
                                                     timeL=[timeL[-1]],zoom=zoom,forcerun=True,
                                                     doplottime=False,contours=[],starC=True,
                                                     plotcolorbar=False)
                
                # Slices
                RUNSLICES = False
                if RUNSLICES:
                    for hydro in ["T","rho","xHII","P","xHeII","xHeIII","Bx","By","Bz","vx","vy","vz",
                                  "EXTRA1","EXTRA2","EXTRA3"][::-1]: # + Lcool
                        for sim in simset:
                            if not "EXTRA" in hydro or "CDMASK3" in sim:
                                makeImageGridTime.MakeFigure([sim],[timesin[-1]],
                                                             name=figname+sim,los=los,hydro=hydro,
                                                             Slice=True,wsink=True,starC=True,
                                                             timeL=[timeL[-1]],zoom=zoom,forcerun=True)

if __name__=="__main__":
    runslices()
