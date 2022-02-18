'''
Plots for the test of the uniform wind bubble with the CD mask
Sam Geen, February 2022
'''

from startup import *

import plotproperties, makeImageGridTime

uniformtestfolder = mainsimfolder+"91_maskcdtest/"
uniformtestturbfolder  = mainsimfolder+"92_maskcdtest_turb/"
uniformwindlum = 2e38 # ergs/s

def MakeThornton():
    global simfolders
    simfolders = {}
    simfolders["Th100_CD"] = uniformtestfolder+"06_thornton_cd/"
    simfolders["Th100_NoCD"] = uniformtestfolder+"07_thornton_nocd/"
    simfolders["Th100_Adiabatic"] = uniformtestfolder+"08_thornton_adiabatic/"
    _MakeSims()

def MakeWindConvergeSims():
    global simfolders
    simfolders = {}
    simfolders["LVL07_MASKCD"] = uniformtestfolder+"11_test_windonlyrefine_lvl7/"
    simfolders["LVL07_NOMASK"] = uniformtestfolder+"12_test_windonlyrefine_lvl7_nocd/"
    simfolders["LVL08_MASKCD"] = uniformtestfolder+"09_test_windonlyrefine_lvl8/"
    simfolders["LVL08_NOMASK"] = uniformtestfolder+"10_test_windonlyrefine_lvl8_nocd/"
    simfolders["LVL09_MASKCD"] = uniformtestfolder+"03_test_windonlyrefine/"
    simfolders["LVL09_NOMASK"] = uniformtestfolder+"04_test_windonlyrefine_nomaskcd/"
    simfolders["LVL10_MASKCD"] = uniformtestfolder+"13_test_windonlyrefine_lvl10/"
    simfolders["LVL10_NOMASK"] = uniformtestfolder+"14_test_windonlyrefine_lvl10_nocd/"
    _MakeSims()

def MakeTurbWindConvergeSims():
    global simfolders
    simfolders = {}
    simfolders["LVL07\_MASKCD"] = uniformtestturbfolder+"02_periodic_lvl7/"
    simfolders["LVL07\_NOMASK"] = uniformtestturbfolder+"06_lvl7_nocd/"
    simfolders["LVL08\_MASKCD"] = uniformtestturbfolder+"03_lvl8/"
    simfolders["LVL08\_NOMASK"] = uniformtestturbfolder+"07_lvl8_nocd/"
    simfolders["LVL09\_MASKCD"] = uniformtestturbfolder+"04_lvl9/"
    simfolders["LVL09\_NOMASK"] = uniformtestturbfolder+"08_lvl9_nocd/"
    simfolders["LVL10\_MASKCD"] = uniformtestturbfolder+"05_lvl10/"
    simfolders["LVL10\_NOMASK"] = uniformtestturbfolder+"09_lvl10_nocd/"
    _MakeSims()
    
def _MakeSims():
    global simfolders, hamusims
    hamusims = {}
    for simname, folder in simfolders.items():
        simexists = True
        try:
            sim = Hamu.Simulation(simname,silent=True)
        except KeyError:
            simexists = False
        if not simexists:
            label = simname
            sim = Hamu.MakeSimulation(simname,folder,label)
        hamusims[simname] = sim


def FindThorntont0(sim):
    # Find the maximum cooling time of the hot bubble
    t,Lwc = timefuncs.timefunc(sim,plotproperties.windLcoolinsnap,
                               processes=plotproperties.nprocs)
    tmaxcool = t[np.where(Lwc == np.max(Lwc))][0]
    return tmaxcool

def ERFinalThorntonForSim(sim):
    t0 = FindThorntont0(sim)
    tf = 13.0 * t0
    snap = timefuncs.findsnapshot(sim,(tf,"Myr"))
    ER = plotproperties.windenergyinsnap(snap)
    return ER

def ERFinalThornton(sims):
    for sim in sims:
        print(sim.Name(),ERFinalThorntonForSim(sim))

def WindConvergeEnergy(setname):
    plt.clf()
    nocdsims = [v for k, v in hamusims.items() if "NOMASK" in k]
    cdsims = [v for k, v in hamusims.items() if "MASKCD" in k]
    for line, sims in zip(["--","-"],[nocdsims, cdsims]):
        for col, sim in zip(["k","r","g","b","m"],sims):
            t, we = plotproperties.windenergy(sim)
            print("we.shape", we.shape)
            print("we.min(), we.max()", we.min(), we.max())
            wemitted = t * Myrins * uniformwindlum
            weretained = we[:,0] / wemitted
            plt.xlabel("Time / Myr")
            plt.ylabel("Fraction of Wind Energy Retained")
            plt.plot(t,weretained,linestyle=line,label=sim.Name(),color=col)
    plt.legend(frameon=False,fontsize="x-small")
    plt.ylim([0.0,1.0])
    plt.savefig("../plots/convergetest_windenergy_"+setname+".pdf")

def WindConvergeRadius(setname):
    plt.clf()
    nocdsims = [v for k, v in hamusims.items() if "NOMASK" in k]
    cdsims = [v for k, v in hamusims.items() if "MASKCD" in k]
    for line, sims in zip(["--","-"],[nocdsims, cdsims]):
        for col, sim in zip(["k","r","g","b","m"],sims):
            #t, we = plotproperties.windenergy(sim)
            #print("we.shape", we.shape)
            #print("we.min(), we.max()", we.min(), we.max())
            #wemitted = t * Myrins * uniformwindlum
            #weretained = we[:,0] / wemitted
            t, radius = plotproperties.radius(sim)
            plt.xlabel("Time / Myr")
            plt.ylabel("Radius / pc")
            plt.plot(t,radius,linestyle=line,label=sim.Name(),color=col)
    plt.legend(frameon=False,fontsize="x-small")
    plt.savefig("../plots/convergetest_windradius_"+setname+".pdf")

        
if __name__=="__main__":
    #MakeThornton()
    #ERFinalThornton(hamusims.values())
    #MakeWindConvergeSims()
    MakeTurbWindConvergeSims()
    #WindConvergeEnergy("uniformturb")
    WindConvergeRadius("uniformturb")
    MakeWindConvergeSims()
    #WindConvergeEnergy("uniformturb")
    WindConvergeRadius("uniform")
    
