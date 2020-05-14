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

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import hydrofuncs
import snaptime

import linestyles

import pymses

import sinks
import stellars
import timefuncs
import singlestar
import rtcooling

from pymses.utils import constants as C 

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
amunold = False
if amunold:
    numsims = 4
    imf1sims = ["IMF1_"+str(i).zfill(2) for i in range(1,numsims+1)]
    imf2sims = ["IMF2_"+str(i).zfill(2) for i in range(1,numsims+1)]
    massivesims = ["MASS_"+str(i).zfill(2) for i in range(1,numsims+1)]
    allsims = imf1sims + imf2sims + massivesims

    # Simulation locations
    #mainsimfolder = "/home/hd/hd_hd/hd_mp149/MCRT/runs/"
    mainsimfolder = "/home/stgeen0/MCRT/runs/"
    imf1simfolder = mainsimfolder+"/58_AMUN/"
    imf2simfolder = mainsimfolder+"/60_AMUN_IMF2/"
    mainsimfolder2 = "/home/stgeen0/sds/sd16j009/REBEKKA_BIERI/MCRT/runs/"
    imf1simfolder2 = mainsimfolder2+"/58_AMUN_IMF1/"
    imf2simfolder2 = mainsimfolder2+"/60_AMUN_IMF2/"
    masssimfolder = mainsimfolder+"64_AMUN_IMF2_M1e5/"

    simfolders = {}
    # IMF1 folders
    # Do this by hand since we're split between different user folders
    simfolders["IMF1_01"] = imf1simfolder+"01_nofbk/"
    simfolders["IMF1_02"] = imf1simfolder+"Y2_justUVphoto_fastlight"
    simfolders["IMF1_03"] = imf1simfolder+"Y3_justWind_fastlight"
    simfolders["IMF1_04"] = imf1simfolder+"Y4_UVphoto+Wind_fastlight"
    if numsims > 4:
        simfolders["IMF1_05"] = imf1simfolder2+"05_UVphoto+press+Wind"
        simfolders["IMF1_06"] = imf1simfolder2+"06_UVphoto+press+IRpress+Wind"

    # IMF2 folders
    # NOTE: no IMF2_01 because it's the same as IMF1_01
    simfolders["IMF2_01"] = imf1simfolder+"01_nofbk/"
    simfolders["IMF2_02"] = imf2simfolder+"Y2_justUVphoto_fastlight"
    simfolders["IMF2_03"] = imf2simfolder+"Y3_justwind_fastlight"
    simfolders["IMF2_04"] = imf2simfolder+"Y4_UVphoto+Wind_fastlight"
    if numsims > 4:
        simfolders["IMF2_05"] = imf2simfolder2+"05_UVphoto+press+Wind"
        simfolders["IMF2_06"] = imf2simfolder2+"06_UVphoto+press+IRpress+Wind"
        
    # Massive folders
    simfolders["MASS_01"] = masssimfolder+"01_nofbk"
    simfolders["MASS_02"] = masssimfolder+"Y2_justUVphoto_fastlight"
    simfolders["MASS_03"] = masssimfolder+"Y3_justwind_fastlight"
    simfolders["MASS_04"] = masssimfolder+"Y4_UVphoto+Wind"

if not amunold:
    mainsimfolder = "/home/stgeen0/MCRT/runs/"
    amunfolder = mainsimfolder+"69_AMUN_onestar/"
    simfolders = {}
    simfolders["NOFB"] = amunfolder+"01_nostar"
    simfolders["UV_120"] = amunfolder+"A2_uvonly120_correctposition"
    simfolders["UVWIND_120"] = amunfolder+"A3_uv+winds_120_correctposition"
    simfolders["NOFB_DENSE"] = amunfolder+"04_nostar_dense"
    simfolders["UV_120_DENSE"] = amunfolder+"05_uvonly_120_dense"
    simfolders["UVWIND_120_DENSE"] = amunfolder+"06_uv+winds_120_dense"
    simfolders["UV_60"] = amunfolder+"A7_uvonly_correctposition"
#    simfolders["UV_30"] = amunfolder+"08_uvonly_30"
    simfolders["UV_30"] = amunfolder+"A8_uvonly_correctposition"
    simfolders["UV_15"] = amunfolder+"09_uvonly_15"
    simfolders["SN_120"] = amunfolder+"10_snonly_120"
    simfolders["UVWIND_60"] = amunfolder+"11_uv+winds_60"
    simfolders["UVPRESS_120"] = amunfolder+"12_uv+press_120_dense"
    simfolders["UVWIND_30"] = amunfolder+"13_uv+winds_30"
    simfolders["UV_120_NOTURB"] = amunfolder+"X3_uvonly_120_spherical"
    simfolders["UVWINDPRESS_30"] = amunfolder+"14_uv+winds+press_30"
    simfolders["UVWINDPRESS_60"] = amunfolder+"16_uv+winds+press_60"
    simfolders["UVWINDPRESS_120"] = amunfolder+"15_uv+winds+press_120"
    allsims = simfolders.keys()
    
# Populate list of Hamu simulations
# TODO - make these on demand rather than on loadup?
#        (Might save time if we use 
hamusims = {}
def _MakeSims():
    for simname, folder in simfolders.iteritems():
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
startableloc = "/home/stgeen0/StellarSources/Compressed/singlestar_z0.014"
singlestar.star_setup(startableloc)

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
        print "Making directory", folder
        os.makedirs(folder)
    except:
        print "Not making directory", folder, "- it already exists!"
        pass

def alpha_B_HII(T):
    # HII recombination rate
    # input  : T in K 
    # output : HII recombination rate (in cm3 / s)
    l = 315614./T
    a = 2.753e-14 * l**1.5 / (1. + (l/2.74)**0.407)**2.242
    return a   

# Add bespoke derived hydro variables
_ionemissfunc = lambda ro: lambda dset: dset["rho"]**2 * dset["xHII"]
ionemission = hydrofuncs.Hydro("Warm Emission",_ionemissfunc,
                               ["rho","xHII"],"PuRd_r","log",(None,None))
hydrofuncs.allhydros["ionemission"] = ionemission
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
hydrofuncs.allhydros["xrayemission2"] = xrayemission2# ---
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
# Set new colours for fields
hydrofuncs.allhydros["P"].ColourMap("Reds")
hydrofuncs.allhydros["Lcool"].ColourMap("BuGn")
# ---
# Done!
print "Imported various project-global modules"
