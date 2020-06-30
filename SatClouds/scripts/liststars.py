'''
List all stars in the simulations
Sam Geen, March 2018
'''

from startup import *

def runforsim(simname,outfile):
    nl = "\n"
    outstr = "\\textbf{"+linestyles.Label(simname)+"} "
    masses = list()
    sim = hamusims[simname]
    FindStellarHamu = Hamu.Algorithm(stellars.FindStellar)
    for snap in sim.Snapshots():
        stellar = FindStellarHamu(snap)
        for mass in stellar.mass:
            if not mass in masses:
                masses.append(mass)
    for mass in masses:
        outstr += "{:.1f}".format(mass)+", "
    outstr = outstr[:-2]
    outstr += nl+nl
    print outstr
    outfile.write(outstr)

def run(simnames,plotname):
    outfile = open("../plots/starlist"+plotname+".tex","w")
    for simname in simnames:
        runforsim(simname,outfile)
    outfile.close()

if __name__=="__main__":
    run(imfsims,"imf")
    run(icsims,"ic")
