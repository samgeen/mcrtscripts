'''
Plot the cloud lifetime against t_ff in the cloud at the start
Sam Geen, March 2018
'''

from startup import *

import tsfe, findproperties
import scipy.interpolate

meandensinsnap = Hamu.Algorithm(findproperties.meandensinsnap)

def tfffromnH(nH):
    tff = ((3 * np.pi) / (32 * G * nH * mH / X))**0.5
    return tff/Myrins

def tsfephases(sim,tstart):
    t,sfe = timefuncs.timefunc(sim,tsfe.tsfeinsnap)
    sfemax = np.max(sfe)
    sfe /= sfemax
    fsfe = scipy.interpolate.interp1d(sfe,t)
    t1 = fsfe(0.001)
    t2 = fsfe(120.0/1e4/sfemax)
    t3 = fsfe(240.0/1e4/sfemax)
    t4 = fsfe(0.999)
    return t1, t2, t3, t4

def runforsim(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tstart = 0.5*tffcloud_code
    tstartMyr = 0.5*tffcloud_Myr
    snap = sim.FindAtTime(tstart)
    #mdens = meandensinsnap(snap,nHthresh=50.0)
    t1,t2,t3,t4 = tsfephases(sim,tstart)
    #tff = tfffromnH(mdens)
    #print "TFF", tff
    #t1 -= tstartMyr
    #t2 -= tstartMyr
    #t3 -= tstartMyr
    f1 = (t2-t1)
    f2 = (t3-t2)
    f3 = (t4-t3)
    return f1, f2, f3

def run(simnames,plotname):
    plt.clf()
    isim = 0
    labels = []
    f1s = []
    f2s = []
    f3s = []
    simnums = range(1,len(simnames)+1)
    for simname in simnames:
        isim += 1
        f1, f2, f3 = runforsim(simname)
        f1s.append(f1)
        f2s.append(f2)
        f3s.append(f3)
        col = linestyles.Colour(simname)
        label = linestyles.Label(simname)
        labels.append(label)
        for f in [f1,f2,f3]:
            plt.scatter([isim],[f],s=80,marker="o",c=col,edgecolors="k",zorder=2)
    fmeantot = 0.0
    for txt,line,style in zip(["Type I", "Type II", "Type III"],[f1s,f2s,f3s],["-","--",":"]):
        fmean = np.sum(line)/float(len(line))
        fmeantot += fmean
        print txt," mean time", fmean," Myr"
        plt.plot(simnums,line,"k"+style,zorder=1,alpha=0.5)
        plt.plot(simnums,np.array(line)*0.0+fmean,"k"+style,zorder=1)
        plt.gca().annotate(txt,(simnums[-1]*1.02,fmean),fontsize="x-small",clip_on=True)
    print "Total time:", fmeantot, "Myr"
    plt.xticks(simnums,labels,fontsize="x-small",rotation="vertical")
    plt.xlabel("Simulation")
    plt.ylabel("Time of Phase / Myr")
    plt.xlim(0,len(simnames)+3)
    plt.savefig("../plots/cloudlifetime_"+plotname+".pdf")

if __name__=="__main__":
    run(icsims,"ic")
    run(imfsims,"imf")
