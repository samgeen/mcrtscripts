from startup import *

import matplotlib.patheffects as pe 

from matplotlib.lines import Line2D

import rdmfile

def plotformass(mass,vartoplot,rdm):
    print "Plotting",vartoplot," for mass",mass,"Msun"
    agesins = np.arange(0,10.0,0.02)*Myrins
    dtins = 1e-6*Myrins
    energies = []
    masslosses = []
    HIIs = []
    HeIIs = []
    HeIIIs = []
    mstr = str(int(mass))
    colour = linestyles.Colour(mstr)
    for ageins in agesins:
        if "wind" in vartoplot:
            if not "energy" in vartoplot:
                energymass = singlestar.star_winds(mass,ageins,dtins)
                energies.append(energymass[0]/dtins)
                masslosses.append(energymass[1]/dtins)
            else:
                energymass = singlestar.star_winds(mass,0,ageins)
                energies.append(energymass[0])
                masslosses.append(energymass[1])
        else:
            nphotons = singlestar.star_radiation(mass,ageins,dtins) 
            HIIs.append(nphotons[2]/dtins)
            HeIIs.append(nphotons[3]/dtins)
            HeIIIs.append(nphotons[4]/dtins)
    vars = {"windlum": energies,
            "windmassloss": masslosses,
            "windenergy": energies,
            "photoHII": HIIs,
            "photoHeII": HeIIs,
            "photoHeIII": HeIIIs}
    label = mstr+" "+Msolar
    if vartoplot != "photoall":
        x = agesins/Myrins
        y = vars[vartoplot]
        plt.plot(x,y,color=colour,label=label,alpha=0.9,
                 path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
        rdm.AddPoints(x,y,label=label)
        return None
    else:
        linenames = ["HII","HeII","HeIII"]
        lines = ["-","--",":"]
        for v, line in zip(linenames,lines):
            y = vars["photo"+v]
            plt.plot(agesins/Myrins,y,color=colour,linestyle=line,label=label,alpha=0.9,
                     path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
            label = None
        # Make other legend for groups
        legelements = [Line2D([0],[0],color="k",
                              linestyle=l,label=n,
                              alpha=0.9,
                              path_effects=[pe.Stroke(linewidth=5, foreground='k'),
                                            pe.Normal()]) for n, l in zip(linenames, lines)]
        legend2 = plt.gca().legend(handles=legelements,framealpha=0.0)
        return legend2
        
def makeplot(masses,vartoplot,ylabel):
    plt.clf()
    rdm = rdmfile.RDMFile(__file__) 
    for mass in masses[::-1]:
        leg2 = plotformass(mass,vartoplot,rdm)
    plt.legend(frameon=False,fontsize="small")
    if leg2 is not None:
        plt.gca().add_artist(leg2)
    plt.yscale("log")
    plt.xlabel("Age / Myr")
    plt.xlim([0.0,8.0])
    plt.ylabel(ylabel)
    figname = "../plots/stellarevolution_"+vartoplot+".pdf"
    plt.savefig(figname)
    rdm.Write(figname)
    
def run():
    masses = [30.0,60.0,120.0]
    makeplot(masses,"photoHII","Photon Emission Rate (HII) /s")
    makeplot(masses,"photoHeII","Photon Emission Rate (HeII) /s")
    makeplot(masses,"photoHeIII","Photon Emission Rate (HeIII) /s")
    makeplot(masses,"windlum","Wind Luminosity / erg/s")
    makeplot(masses,"windenergy","Cumulative Wind Energy / erg")
    makeplot(masses,"windmassloss","Wind Mass Loss / g/s")
    print "Done!"
        
if __name__=="__main__":
    run()
    
