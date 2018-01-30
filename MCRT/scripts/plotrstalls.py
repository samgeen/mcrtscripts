'''
Find rstall in the HIISN model
Sam Geen, July 2015
'''

import customplot
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rc('text', usetex=True)

import numpy as np

import Hamu
Hamu.Workspace("MCRT")

import outflowmodel, stallradius, ragacompare, profilesphere, testprofilefit
import edgeradius

mH = ragacompare.mH
X = ragacompare.X

def stallrad(n0,r0,w,flux=1e48):
    r0 = r0*outflowmodel.pcincm
    r0w = r0**w
    rho0 = n0 * mH / X
    alphaB = ragacompare.beta2
    G = ragacompare.G
    ci = ragacompare.cs
    rstall = ((3.0*flux / (4*np.pi * alphaB))**0.5 * \
                  (1.0/r0w**2) * (3 - w) * ci**2 *mH/X / \
                  (8 * np.pi * G * rho0**2))**(2.0 / (7.0 - 4*w))
    return rstall/outflowmodel.pcincm

def fts(num):
    # Float to string
    return "{:.1f}".format(num) 

def makelabel(r0, n0):
    label = "$r_{0}$ = "+fts(r0)+" pc , $n_{0}$ = "+fts(n0)+" cm$^{-3}$"
    return label

def dotext(r0,n0,w,flux,colour="k"):
    rfact = n0 * r0**w
    rstall = stallrad(n0,r0,w,flux)
    plt.text(rfact,rstall,
             makelabel(r0,n0),
             fontsize="x-small",color=colour)

def domarker(r0,n0,w,flux,xlim,ylim,marker,colour):
    up = r"\uparrow"
    down = r"\downarrow"
    left = r"\leftarrow"
    right = r"\rightarrow"
    nw = r"\nwarrow"
    ne = r"\nearrow"
    sw = r"\swarrow"
    se = r"\searrow"
    arrows = [[sw,down,se],[left,"",right],[nw,up,ne]]
    rfact = n0 * r0**w
    rstall = stallrad(n0,r0,w,flux)
    text = r"$"+marker+r"$"
    xs = [xlim[0],rfact,xlim[1]]
    ys = [ylim[0],rstall,ylim[1]]
    ix = 1
    ix = 0 if rfact < xlim[0] else ix
    ix = 2 if rfact > xlim[1] else ix
    iy = 1
    iy = 0 if rstall < ylim[0] else iy
    iy = 2 if rstall > ylim[1] else iy
    ar = arrows[iy][ix]
    if ar:
        text = r"$"+marker +" "+ ar+r"$"
        print ix, iy, text
        #text = r"$"+ar+r"$"
    plt.text(xs[ix],ys[iy],
             text,
             fontsize="small",color=colour,
             horizontalalignment='center',
             verticalalignment='center')

def run(flux):
    r0 = 10.0**np.arange(-1,2.0001,0.02)
    n0 = 10.0**np.arange(2.0,5.0001,0.02)
    ws = np.array([0.0,0.5,1.0,1.5,2.0])
    colours = ["#fcc5c0",
               "#fa9fb5",
               "#f768a1",
               "#c51b8a",
               "#7a0177"]
    icol = 0
    flims = [10,1e6]
    rlims = [1e-2,1e3]
    plt.plot([10],[1e4],color="w",label="$w$")
    table = r"\begin{tabular}{ c c } "
    firsttable = True
    dotable = True
    rnms = [(2,100,r"\diamond"),
            (2,1e3,r"\circ"),
            (2,1e4,r"\triangle")]
    for w in ws:
        print "Calculating for w =",w
        rfact = n0 * r0**w
        rstalls = stallrad(n0,r0,w,flux)
        valid = (rstalls > rlims[0]) * (rstalls < rlims[1])
        rfact = rfact[valid]
        rstalls = rstalls[valid]
        r0c = r0[valid]
        n0c = n0[valid]
        plt.plot(rfact,rstalls,label=str(w),color=colours[icol])
        # Do labels
        #for il in [0,-1]:
        #    plt.text(rfact[il],rstalls[il],
        #             makelabel(r0c[il],n0c[il]),
        #             fontsize="x-small")
        for rnm in rnms:
            r, n, marker = rnm
            domarker(r,n,w,flux,
                     flims,rlims,
                     marker,colours[icol])
            if dotable:
                if firsttable:
                    firsttable = False
                else:
                    table += r" \\ "
                table += r"$"+marker+r"$ & "+makelabel(r,n)
        dotable = False
        #dotext(10,1e3,w,flux,colours[icol])
        icol += 1
    table += r"\end{tabular}"
    print "Making plot..."
    plt.xlabel("$(n_{0} $ / cm$^{-3}) (r_{0}$ / pc $) ^ {w}$")
    plt.ylabel("$r_{stall}$ / pc")
    plt.xlim(flims)
    plt.ylim(rlims)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize="small",frameon=False)
    plt.text(flims[0]*1.5,rlims[0]*1.5,table,
             horizontalalignment='left',
             verticalalignment='bottom',size="x-small")
    plt.savefig("../plots/rstalls/rstalls.pdf")
    print "Done!"
    
    

if __name__=="__main__":
    flux = 1e48
    run(flux)
