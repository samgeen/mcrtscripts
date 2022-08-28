'''
Code to run on script startup for this project
Includes imports that can be used by other scripts
Sam Geen, June 2016
'''

import sys, os, glob
import numpy as np
# Import from the main scripts folder
sys.path.append("../../scripts")
sys.path.append("../../Cooling/f2py")
import customplot

sys.path.append("/home/stgeen0/Programming/")

import HamuLite as Hamu
Hamu.UNTARPATH = "../cache/"

print("Import 1...")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import hydrofuncs
import snaptime

import linestyles

print("Import 2...")
import pymses


print("Import 3...")
import sinks
import stellars
import timefuncs
import singlestar
import rtcooling

from pymses.utils import constants as C 

print("Import done")

#Hamu.Workspace("HIISFE")

# Physical conversions
X = 0.76
mH = 1.6735326381e-24
kB = 1.38062e-16
G = 6.674e-8
gamma = 1.4 # RAMSES hard-coded
pcincm = 3.086e18
Msuning = 1.9891e33
mHing = 1.66e-24
Myrins = 3.1556926e13


# Simulation units (to prevent need to dig into code)
unit_l      =  0.308000000000000E+19
unit_d      =  0.232474340000000E-23 
unit_t      =  0.253950790000000E+16

# Simulation global properties
tffcloud_Myr = 4.21875 # tff of cloud as a whole in Myr
tffcloud_code = tffcloud_Myr / (unit_t/Myrins)

# Some misc definitions
Msolar = "M$_{\odot}$"

# Plot folder
plotfolder = "../plots/"

# Simulation names
mainsimfolder = "/home/stgeen0/MCRT/runs/" # on ISMSIM1
cartsimfolder = "/home/samgeen/amun/runs/"
if os.path.exists(cartsimfolder):
    mainsimfolder = cartsimfolder

windshellfolder = mainsimfolder+"84_windshell/"
simfolders = {}
UseTestRuns = False
if UseTestRuns:
    simfolders["NOSHELL_NOMASK"] = windshellfolder+"11_uvwindpress_11grad_shellrttest/" # No wind shell imposed at t=0
    simfolders["SHELL_NOMASK"] = windshellfolder+"12_innershelltest/" # Wind shell imposed at t=0
    simfolders["SHELL_CDMASK"] = windshellfolder+"13_maskcdtest/" # As above but with contact discontinuity masked
    simfolders["SHELL_CDMASK2"] = windshellfolder+"14_maskcdtest2/" # Debugging
    simfolders["SHELL_CDMASK3"] = windshellfolder+"15_maskcdtest3/" # Debugging
    simfolders["SHELL_CDMASK4"] = windshellfolder+"16_maskcdtestnew/" # Debugging
else:
    # Different seeds
    simfolders["SEED0_35MSUN_CDMASK_WINDUV"] = windshellfolder+"20_maskcd_35Msun/"
    simfolders["SEED1_35MSUN_CDMASK_WINDUV"] = windshellfolder+"23_newseed_35Msun/"
    simfolders["SEED2_35MSUN_CDMASK_WINDUV"] = windshellfolder+"30_newseed2_35Msun/"
    simfolders["SEED3_35MSUN_CDMASK_WINDUV"] = windshellfolder+"31_newseed3_35Msun/"
    #simfolders["SEED4_35MSUN_CDMASK_WINDUV"] = windshellfolder+"32_newseed4_35Msun/"
    #simfolders["SEED5_35MSUN_CDMASK_WINDUV"] = windshellfolder+"33_newseed5_35Msun/"
    # NOTE: Run 32 is accidentally identical to run 31
    simfolders["SEED4_35MSUN_CDMASK_WINDUV"] = windshellfolder+"33_newseed5_35Msun/"
    # Big clouds
    simfolders["SEED0_35MSUN_CDMASK_WINDUV_BIGCLOUD"] = windshellfolder+"21_biggercloud_35Msun/"
    simfolders["SEED1_35MSUN_CDMASK_WINDUV_BIGCLOUD"] = windshellfolder+"22_biggercloud_newseed_35Msun/" 
    # Turning off physics tests
    simfolders["SEED1_35MSUN_NOCDMASK_WINDUV"] = windshellfolder+"24_newseed_35Msun_nocd/"
    simfolders["SEED1_35MSUN_CDMASK_WINDUV_NOB"] = windshellfolder+"25_newseed_35Msun_noB/"
    simfolders["SEED1_35MSUN_CDMASK_WINDUV_NOREFINE"] = windshellfolder+"26_newseed_35Msun_norefine/"
    simfolders["SEED1_35MSUN_NOCDMASK_WINDUV_NOREFINE"] = windshellfolder+"27_nocdnorefine/"
    simfolders["SEED1_35MSUN_WINDUV_WEAKB"] = windshellfolder+"34_newseed_35Msun_vweakB/"
    # Turn off feedback sources
    simfolders["SEED1_35MSUN_CDMASK_WIND"] = windshellfolder+"28_newseed_35Msun_justwind/"
    simfolders["SEED1_35MSUN_NOCDMASK_WIND"] = windshellfolder+"37_newseed_35Msun_justwind_nocd/"
    #simfolders["SEED1_35MSUN_CDMASK_UV"] = windshellfolder+"29_newseed_35Msun_justuv/"
    simfolders["SEED1_35MSUN_CDMASK_UV"] = windshellfolder+"Z29_newseed_35Msun_justuv_norefine/"
    simfolders["SEED1_35MSUN_CDMASK_NOFB"] = windshellfolder+"N23_nofb/"
    #Extra no feedback stuff
    simfolders["SEED0_35MSUN_CDMASK_NOFB"] = windshellfolder+"N20_nofb/"
    simfolders["SEED2_35MSUN_CDMASK_NOFB"] = windshellfolder+"N30_nofb/"
    simfolders["SEED3_35MSUN_CDMASK_NOFB"] = windshellfolder+"N31_nofb/"
    simfolders["SEED4_35MSUN_CDMASK_NOFB"] = windshellfolder+"N33_nofb/"
allsims = list(simfolders.keys())


# Simulation sets
seedlabels = {"SEED0_35MSUN_CDMASK_WINDUV":"Seed Ampere",
              "SEED1_35MSUN_CDMASK_WINDUV":"Seed Bellecour",
              "SEED2_35MSUN_CDMASK_WINDUV":"Seed Carnot",
              "SEED3_35MSUN_CDMASK_WINDUV":"Seed Duhamel",
              "SEED4_35MSUN_CDMASK_WINDUV":"Seed d'Enghien"}
seedset = list(seedlabels.keys())


physicslabels = {"SEED1_35MSUN_CDMASK_WINDUV": "Mask \& Refine",
                 "SEED1_35MSUN_CDMASK_WINDUV_NOREFINE": "Mask, No Refine",
                 "SEED1_35MSUN_NOCDMASK_WINDUV": "No Mask, Refine",
                 "SEED1_35MSUN_NOCDMASK_WINDUV_NOREFINE":"No Mask, No Refine"}
physicsset = list(physicslabels.keys())

fblabels = {"SEED1_35MSUN_CDMASK_WINDUV":"Wind \& UV", 
         "SEED1_35MSUN_CDMASK_WIND":"Wind Only", 
         "SEED1_35MSUN_CDMASK_UV":"UV Only",
         "SEED1_35MSUN_CDMASK_NOFB":"No Feedback"}
fbset = list(fblabels.keys())
# "SEED1_35MSUN_CDMASK_WINDUV_NOB" - no star yet

singlelabels = {"SEED1_35MSUN_CDMASK_WINDUV":"Fiducial Run"}
singleset = list(singlelabels.keys())


windonlylabels = {"SEED1_35MSUN_CDMASK_WIND":"Wind Only, Mask On",
               "SEED1_35MSUN_NOCDMASK_WIND":"Wind Only, Mask Off"}
windonlyset = list(windonlylabels.keys())

hotchampagnelabels = {"SEED2_35MSUN_CDMASK_WINDUV":""} # Only one simulation, no label needed
hotchampagneset = list(hotchampagnelabels.keys())

simsets = {"seeds":seedset,
           "windonly":windonlyset,
           "single":singleset,
           "fb":fbset, "physics":physicsset,
           "hotchampagne":hotchampagneset}
simlabels = linestyles.simlabels
simlabels.update({"single":singlelabels, "windonly":windonlylabels, "fb":fblabels, "physics":physicslabels, "seeds":seedlabels,"hotchampagne":hotchampagnelabels})
    
# Populate list of Hamu simulations
# TODO - make these on demand rather than on loadup?
#        (Might save time if we use 
hamusims = {}
def _MakeSims():
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
_MakeSims()

# Set up single star module
print("Importing singlestar...")

startableloc = "/home/stgeen0/StellarSources/Outputs/singlestar_z0.014"
singlestar.star_setup(startableloc)    

print("Done")

# Useful functions
# N_H to/from A_k (from Lombardi+ 2010)
akconst = 1.75e22
# N_H to/from A_v (from Savage & Mathis 1979)
avconst = 1.87e21
def AktoNH(Ak):
    return akconst*Ak

def NHtoAk(NH):
    return NH/akconst

def AvtoNH(Av):
    return avconst*Av

def NHtoAv(NH):
    return NH/avconst

def ColDenstoNH(dens):
    return dens*X/mHing

def NHtoColDens(NH):
    return NH*mHing/X

def ColDenstoAk(dens):
    return NHtoAk(ColDenstoNH(dens))

def AktoColDens(Ak):
    return NHtoColDens(AktoNH(Ak))

# Conversion constant from Kirk et al 2013
muH2 = 2.83
def NH2toColDens(NH2):
    return NH2*mHing*muH2

def ColDenstoNH2(dens):
    return dens/(mHing*muH2)

def NH2toNH(NH2):
    return NH2*X*muH2

def NHtoNH2(NH):
    return NH/(X*muH2)

def MakeDirs(folder):
    try:
        print("Making directory", folder)
        os.makedirs(folder)
    except:
        print("Not making directory", folder, "- it already exists!")
        pass

def alpha_B_HII(T):
    # HII recombination rate
    # input  : T in K 
    # output : HII recombination rate (in cm3 / s)
    l = 315614./T
    a = 2.753e-14 * l**1.5 / (1. + (l/2.74)**0.407)**2.242
    return a   

# Add bespoke derived hydro variables
def _ionemissfunc4(ro):
    def ionfuncinner(dset):
        xHII = dset["xHII"]
        emiss = (dset["rho"] * xHII)**2
        emiss[xHII < 1e-4] = 1e-20
        return emiss
    return ionfuncinner
ionemission4 = hydrofuncs.Hydro("Warm Emission",_ionemissfunc4,
                               ["rho","xHII"],"PuRd_r","log",(None,None))
hydrofuncs.allhydros["ionemission4"] = ionemission4
# ---
def _xrayfunc(ro):
    mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                          0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
    def xrayfunc(dset):
        temp = dset["P"]/dset["rho"]*ro.info["unit_temperature"].express(C.K)*mufunc(dset)
        emiss = temp*0.0
        emiss[(temp > 1e6)] = 1.0
        emiss *= dset["rho"]
        return emiss
    return xrayfunc
xrayemission = hydrofuncs.Hydro("X-Ray Emission Measure",_xrayfunc,
                                ["rho","P","xHII","xHeII","xHeIII"],"YlOrRd_r","log",(None,None))
hydrofuncs.allhydros["xrayemission"] = xrayemission
# ---
def _xrayfunc2(ro):
    mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                          0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
    coolfunc = rtcooling.dedtOnCells(ro)
    def xrayfunc(dset):
        shape = dset["rho"].shape
        temp = dset["P"].flatten()/dset["rho"].flatten()*ro.info["unit_temperature"].express(C.K)#*mufunc(dset)
        temp = np.reshape(temp,shape)
        cool = coolfunc(dset)
        cool[cool < 0.0] = 0.0
        emiss = cool+0.0
        emiss[temp < 1e6] = 0.0 # cool[cool > 0.0].min()*0.1
        return emiss
    return xrayfunc
xrayemission2 = hydrofuncs.Hydro("Hot Emission",_xrayfunc2,
                                ["rho","P","vel","xHII","xHeII","xHeIII","NpHII","NpHeII","NpHeIII"] ,
                                 "YlOrRd_r","log",(None,None))
hydrofuncs.allhydros["xrayemission2"] = xrayemission2

# ---
# Include MaskCD if possible
def _xrayfunc3(ro):
    mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                          0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
    coolfunc = rtcooling.dedtOnCells(ro)
    def xrayfunc(dset):
        shape = dset["rho"].shape
        temp = dset["P"].flatten()/dset["rho"].flatten()*ro.info["unit_temperature"].express(C.K)#*mufunc(dset)
        temp = np.reshape(temp,shape)
        cool = coolfunc(dset)
        cool[cool < 0.0] = 0.0
        if "MaskCD" in dset:
            maskcd = dset["MaskCD"] > 1.0
            cool[maskcd] = 0.0
        emiss = cool+0.0
        emiss[temp < 1e6] = 0.0 # cool[cool > 0.0].min()*0.1
        return emiss
    return xrayfunc
xrayemission3 = hydrofuncs.Hydro("Hot Emission",_xrayfunc3,
                                ["rho","P","vel","xHII","xHeII","xHeIII","NpHII","NpHeII","NpHeIII","MaskCD"] ,
                                 "YlOrRd_r","log",(None,None))
hydrofuncs.allhydros["xrayemission3"] = xrayemission2
# ---
def _lowtempfunc(ro):
    mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"]) + \
                          0.25*0.24*(1.+dset["xHeII"]+2.*dset["xHeIII"]))
    coolfunc = rtcooling.dedtOnCells(ro)
    def emissfunc(dset):
        shape = dset["rho"].shape
        temp = dset["P"].flatten()/dset["rho"].flatten()*ro.info["unit_temperature"].express(C.K)#*mufunc(dset)
        temp = np.reshape(temp,shape)
        cool = coolfunc(dset)
        cool[cool < 0.0] = 0.0
        emiss = cool+0.0
        emiss[temp > 1e3] = 0.0 # cool[cool > 0.0].min()*0.1
        return emiss
    return emissfunc
coolemission = hydrofuncs.Hydro("Cool Emission",_lowtempfunc,
                                ["rho","P","vel","xHII","xHeII","xHeIII","NpHII","NpHeII","NpHeIII"] ,
                                 "GnBu_r","log",(None,None))
hydrofuncs.allhydros["coolemission"] = coolemission
def _dustforming(ro):
    def maskfunc(dset):
        rho = dset["rho"]
        temp = dset["P"].flatten()/dset["rho"].flatten()*ro.info["unit_temperature"].express(C.K)
        mask = np.where(temp<1200.0)*np.where(rho<4.5e10)
        rho[not mask] = 0.0
        return rho
    return maskfunc
dustforming = hydrofuncs.Hydro("Dust Forming Gas",_dustforming,
                                ["rho","P","xHII","xHeII","xHeIII"] ,
                                 "GnBu_r","log",(None,None))
hydrofuncs.allhydros["dustforming"] = dustforming
# Fast-moving material, e.g.shells
def _fastmass6(ro):
    def maskfunc(dset):
        NH = dset["rho"]*ro.info["unit_density"].express(C.H_cc)*ro.info["unit_length"].express(C.cm)
        vels = dset["vel"]
        uvel = ro.info["unit_velocity"].express(C.cm/C.s)
        #print("VELS SHAPE", vels.shape)
        spds = np.sqrt(np.sum(vels**2.0,2))*uvel
        #print("SPD MINMAX", spds.min(), spds.max())
        # Check over 10 km/s
        #print("FASTMASS:")
        #print (NH.shape)
        #print("BEFORE", NH)
        #print("SPDS:", spds)
        #print(mask)
        spdlimit = 5 * 1e5 # 5 km/s
        NH[spds < spdlimit] = NH[spds < spdlimit]*1e-20
        #print("AFTER", NH)
        #print("---")
        return NH
    return maskfunc
fastmass6 = hydrofuncs.Hydro("$N_H(>5~$km/s$)$",_fastmass6,
                                ["rho","vel"] ,
                                 "GnBu_r","log",(20,23.7),surfacequantity=True)
hydrofuncs.allhydros["fastmass6"] = fastmass6
# Set new colours for fields
hydrofuncs.allhydros["P"].ColourMap("Reds")
hydrofuncs.allhydros["Lcool"].ColourMap("BuGn")
# ---
# Done!
print("Imported various project-global modules")
