'''
Create diagram for paper explanation
Sam Geen, October 2020
'''

import customplot

import numpy as np

import matplotlib.pyplot as plt


def run():
    # Set up line data
    nrad = 1000
    radii = np.linspace(0.0,1.0,nrad)
    densities = np.zeros((nrad,))
    temperatures = np.zeros((nrad,))
    irw = nrad//10
    iri = 5*nrad//10
    irs = 9*nrad//10
    # Wind bubble
    nwind = 1e-1
    Twind = 1e7
    densities[:irw] = nwind
    temperatures[:irw] = Twind
    # Ionised shell
    iradii = radii[irw:iri]
    niw = 1e1
    nii = 1e2
    densities[irw:iri] = niw / (1.0 - 0.18*niw*iradii)
    temperatures[irw:iri] = 1e4
    # Neutral shell
    nshell = 1e3
    densities[iri:irs] = nshell
    temperatures[iri:irs] = 1e2
    # External medium
    densities[irs:] = 1e2
    temperatures[irs:] = 1e2
    # Plot
    plt.clf()
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # Do neutral line cut-off
    ax1.plot([radii[iri]]*2,[nwind,nshell],"k--",alpha=0.7)
    # Do temperature line
    ax2.plot(radii,temperatures,"r")
    ax2.set_yscale("log")
    ax2.set_ylabel("T / K",color="r")
    ax2.tick_params('y', colors="r")
    # Do density line
    ax1.plot(radii,densities,"b")
    ax1.set_yscale("log")
    ax1.set_ylabel("$n_{\mathrm{H}}$ / cm$^{-3}$",color="b")
    ax1.tick_params('y', colors='b')

    # Do text labels
    ax1.text(radii[(irw+iri)//2],1.0,"Ionised",horizontalalignment='center', verticalalignment='center')
    ax1.text(radii[(iri+irs)//2],1.0,"Neutral",horizontalalignment='center', verticalalignment='center')
    # Do other plotting stuff
    plt.xticks(radii[[irw,iri,irs]],["$r_w$","$r_i$","$r_s$"])
    ax1.set_xlabel("Radius")
    plt.savefig("../plots/diagram.pdf")


if __name__=="__main__":
    run()
    
