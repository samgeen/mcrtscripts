
import os, sys
sys.path.append("/home/samgeen/Programming/Astro/Weltgeist")
sys.path.append("/home/samgeen/Programming/Astro/WindInUV")

import numpy as np

import weltgeist
import weltgeist.units as wunits # make this easier to type
import stars

def runforstar(mass,metal,rotating):
    star = stars.Star(mass, metal, rotating=rotating)
    times = star.spectimes
    Lionising = star.LIonising(times)
    Lnonionising = star.LNonIonising(times)           
    Teff = star.Teff(times) 
    out = np.array([Lionising,Lnonionising,Teff]).T
    print(out.shape)
    if rotating:
        rottxt = "rotating"
    else:
        rottxt = "nonrotating"
    filename = "../outputs/luminosities/luminosities_M"+str(mass)+"Msun_Z"+str(metal)+"_"+rottxt+".txt"
    np.savetxt(filename, out)


def runforall():
    masses = stars.starmasses
    metals = stars.starmetals
    for rotating in [True, False]:
        for mass in masses:
            for metal in metals:
                runforstar(mass,metal,rotating=rotating)


if __name__=="__main__":
    runforall()