'''
Make raw tables of TSFE versus other parameters
Sam Geen, February 2018
'''

from startup import *

import correlatestructure, starrelations

def run(simnames,relation,params,rhos=[None],times=[None],label=""):
    nsims = len(simnames)
    nparams = len(params)
    nrhos = len(rhos)
    ntimes = len(times)
    arr = np.zeros((nparams+1,nsims))
    if relation == "correlatestructure":
        arr = np.zeros((nparams*nrhos*ntimes+1,nsims))
    isim = -1
    titles = ["sfe"]
    for simname in simnames:
        isim += 1
        iparam = -1
        for param in params:
            for rho in rhos:
                for time in times:
                    iparam += 1
                    if relation == "starrelations":
                        vparam, sfe = starrelations.runforsim(simname,param,"sfe")
                        arr[0,isim] = np.log10(sfe)
                        arr[iparam+1,isim] = np.log10(vparam)
                        if isim == 0:
                            titles += [param]
                    elif relation == "correlatestructure":
                        correlatestructure.parameter = param
                        correlatestructure.rhofloor = float(rho)
                        correlatestructure.tff_fact = time
                        if isim == 0:
                            titles += [param+"_"+str(rho)+"_"+str(time)+"tff"]
                        sfe, vparam, dum = correlatestructure.runforsim(simname)
                        arr[0,isim] = np.log10(sfe*100.0)
                        arr[iparam+1,isim] = np.log10(vparam)
                    else:
                        print "No relation called", relation
                        raise ValueError

    np.savetxt("../plots/rawtsfetable_"+relation+"_"+label+".dat",arr.T,delimiter="\t",header="\t".join(titles))


if __name__=="__main__":
    run(imfsims,"starrelations",["alltimemax","compactness","firstmass","nphotons","nphotonstot","nphotonstff","tracklength"])
    run(icsims,"correlatestructure",["T","L","M","S","C","V"],[10],[0.5,1,1.5,2.0],label="n10")
    run(icsims,"correlatestructure",["T","L","M","S","C"],[10,100,1000],[0.5],label="0p5tff")
