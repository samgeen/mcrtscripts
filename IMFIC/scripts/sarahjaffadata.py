'''
Get stellar object data for Sarah Jaffa
Sam Geen, September 2020
'''

from startup import *

import stellars, snaptime

FindStellarHamu = Hamu.Algorithm(stellars.FindStellar)
FindSinkHamu    = Hamu.Algorithm(sinks.FindSinks) 

def runforsim(simname):
    print "SIM", simname
    sim = hamusims[simname]
    for snap in sim.Snapshots():
        time = snaptime.Myr(snap)
        try:
            os.mkdir ("../tables/stars/"+simname)
        except:
            pass
        outfile = open("../tables/stars/"+simname+"/output"+str(snap.OutputNumber()).zfill(5)+".dat","w")
        def write(vals):
            nval = len(vals[0])
            for i in range(0,nval):
                for val in vals:
                    outfile.write(str(val[i])+" ")
                outfile.write("\n")
        stellar = FindStellarHamu(snap)
        sink = FindSinkHamu(snap)
        stellar.FixPosition(sink)
        outfile.write("x y z mass age lifetime \n")
        write([stellar.x,stellar.y,stellar.z,stellar.mass,time - stellar.tcreated,stellar.lifetime])
        outfile.close()

def run(simnames, suitename):
    for simname in simnames:
        runforsim(simname)

if __name__=="__main__":
    run(imfsims,"imf")
    run(icsims,"ic")
