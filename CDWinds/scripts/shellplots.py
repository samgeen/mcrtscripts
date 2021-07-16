'''
Make images of the cloud/sinks for Marco Padovani project
Rebekka Bieri and Sam Geen, August 2018
'''

from startup import *

import makeImageGridTime

if __name__=="__main__":

    for dense in [False]:
        for mass in [30]: # [30,60,120][::-1]:
            for outtime in [0.1,0.15,0.2]:
                smass = str(mass)
                if not dense or mass == 120:
                    simset = ["UVWINDPRESS_30_ZOOM_NOGRAD","UVWINDPRESS_30_ZOOM_11GRAD","UVWINDPRESS_30_ZOOM_11GRAD_SHELL"] 
                    setname = "crset_"+smass+"Msun"
                    simwindname = "UVCR_"+smass
                    #if dense:
                    #    simset = [x+"_DENSE" for x in simset]
                    #    simwindname = "UVWIND_"+smass+"_DENSE"
                    #    setname += "_dense"
                    #times = np.array([0.5, 0.75, 1.])
                    times = np.array([outtime]) # np.array([0.9]) # [0.9] # 3.5 Myr = tstarformed + 0.2 Myr 
                    zoom = 0.25
                    setname = setname+str(times[-1])+"Myr_"+"zoom"+str(zoom)+"_"
                    setname = setname.replace(".","p") # the extra dot confuses latex
                    #timeL = [str(x)+r' t$_{ff}$' for x in times]
                    #timesin = [(time*tffcloud_code,"code") for time in times]
                    timeL = [str(x)+r' Myr' for x in times]
                    timesin = [(time,"MyrFirstStar") for time in times]
                    for los in "yxz":
                        figname = setname+"_"+los
                        allfigname = figname.replace(smass+"Msun","allstars")

                    
                        # Merged emission map - just wind
                        coolhydros = ["coolemission","ionemission","xrayemission2"]
                        timesmerged = [0.1,0.2,0.3,0.4]
                        timesmergedIn = [(time,"MyrFirstStar") for time in timesmerged]
                        timesmergedL = [str(x)+r' Myr' for x in timesmerged]

                        # Slices
                        for hydro in ["Lcool","T","rho","xHII","P","xHeII","xHeIII","Bx","By","Bz","vx","vy","vz"][::-1]:
                            for sim in simset:
                                makeImageGridTime.MakeFigure([sim],[timesin[-1]],
                                                             name=figname+"windshell"+sim,los=los,hydro=hydro,
                                                             Slice=True,wsink=True,starC=True,
                                                             timeL=[timeL[-1]],zoom=zoom,forcerun=True)

                        
