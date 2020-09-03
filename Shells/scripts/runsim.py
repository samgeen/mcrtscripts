'''
Run the simulation
Sam Geen, July 2020
'''

import os, sys
sys.path.append("/home/samgeen/Programming/Astro/Weltgeist")

import numpy as np

#import derivedgeneva, readstar
import weltgeist
import weltgeist.units as wunits # make this easier to type
from scipy.interpolate import interp1d

# Some useful constants
lsolar = 3.828e33 # Solar luminosity / ergs/s
rsolar = 6.955e10 # Solar radius in cm
msolar = 1.98855e33 # Solar mass in g

G = 6.67e-8
sbconstant = 5.6704e-5 # Stefan-Boltzmann constant / erg cm^-2 s^-1 K^-4
c = 2.99792458e10

yrins = 3.154e+7

def run():
    # Make an integrator
    integrator = weltgeist.integrator.Integrator()
    # And the setup
    ncells = 512
    nanalytic = np.zeros((ncells))
    n0 = 100.0 # cm^-3
    T0 = 10.0 # K
    integrator.Setup(ncells = ncells, # 256 cells
            rmax = 10.0*wunits.pc, # 20 pc box
            n0 = n0, # atoms / cm^-3
            T0 = T0, # K
            gamma = 5.0/3.0) # monatomic gas (close enough...)
    hydro = integrator.hydro

    # # Use an isothermal profile
    isothermal = True
    if isothermal:
        r0 = 1.0 * wunits.pc
        hydro.nH[0:ncells] = n0 * (hydro.x[0:ncells] / r0)**(-2.0)
        # Get rid of the singularity at r=0
        hydro.nH[0] = hydro.nH[1]
        # Reset temperature
        hydro.T[0:ncells] = T0

    # Now let's add a star. 
    # First, and this is very important, set the table location
    # This tells the code where the stellar evolution tracks are
    # This will be different for your computer!
    # If you don't have the tables, email me
    # These tables are for Solar metallicity (Z = 0.014)
    # There are also tables for sub-Solar (Z = 0.002)
    # Depends if you want to model our Galaxy or the LMC/SMC
    weltgeist.sources.singlestarLocation = \
        "../../../StellarSources/Compressed/singlestar_z0.014"

    star = weltgeist.sources.TableSource(30.0,radiation=True,wind=True)
    weltgeist.sources.Sources().AddSource(star)

    # Turn cooling on
    weltgeist.cooling.cooling_on = True

    # Now let's render the temperature in red
    temperatureLine = weltgeist.graphics.Line(weltgeist.graphics.red,width=3.0)
    # And a blue line for an analytic fit
    densityLine = weltgeist.graphics.Line(weltgeist.graphics.blue,width=3.0)
    # And a black line for photoionisation
    ionisedLine = weltgeist.graphics.Line(weltgeist.graphics.black,width=3.0)
    # And make a rendering object again as in Exercise 3
    renderer = weltgeist.graphics.Renderer([temperatureLine,densityLine,ionisedLine])

    # Because of the way the rendering module pyglet works, it has to
    #  control the stepping. So let's make a function to give it
    def MyStep(dtRender):
        """
        Function for pyglet to run every timestep

        Parameters
        ----------

        dtRender : float
            wall clock timestep each frame of the rendering
        """
        # Get the hydrogen number density and radius and update the line
        # Note that you can also add points to the line over time
        # This is useful for plotting time-dependent functions like example 2
        x = hydro.x[0:ncells]/wunits.pc
        nH = hydro.nH[0:ncells]
        T = hydro.T[0:ncells]
        xhii = hydro.xhii[0:ncells]
        densityLine.Update(x,np.log10(nH))
        temperatureLine.Update(x,np.log10(T))
        ionisedLine.Update(x,xhii)
        steptime = integrator.time/yrins/1e3
        renderer.Text(str.format('{0:.3f}',steptime)+" kyr")
        # Step the integrator
        integrator.Step()
    
    # Now run the renderer and the simulation will evolve!
    # Press Escape or click the cross to exit
    # Note that this doesn't have axis labels, it's super simple so far
    renderer.Start(MyStep)

    # A few things to notice
    # 1) Winds are a giant pain. They're super fast (up to a few 
    #  1000 km/s) and they make super hot gas (up to 10^8-10^9 K) 
    # This is an issue for the Courant condition. A 3D simulation that 
    # takes 2 days without winds takes weeks with winds included
    # 2) What structures do you see evolving? How do they compare with
    #  https://ui.adsabs.harvard.edu/abs/1977ApJ...218..377W/abstract
    # Try turning radiation on to see what happens. Why?
    # Try plotting xHII instead of log(T) (or add a new line for it)
    # Now try setting the medium to uniform and see what happens

    # In the next example we explore saving the simulation state to file


if __name__=="__main__":
    run()