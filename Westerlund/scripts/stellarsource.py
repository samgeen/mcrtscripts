'''
Deal with stellar sources
Sam Geen, September 2020
'''


import os, sys
sys.path.append("/home/samgeen/Programming/Astro/Weltgeist")
sys.path.append("/home/samgeen/Programming/Astro/WindInUV")

import numpy as np

import weltgeist
import weltgeist.units as wunits # make this easier to type
import stars

class StellarSource(weltgeist.sources.AbstractSource):
    """
    Source of energy & photons based on a lookup table
    """
    def __init__(self,mass,metal,tbirth=0.0,radiation=True,wind=True,rotating=True):
        """
        Constructor
    
        Parameters
        ----------

        mass : float
            Mass of star in solar masses
        tbirth : float
            Birth time of the star in seconds
        radiation : bool
            Turn radiation on?
        wind : bool
            Turn winds on?
        rotating : bool
            Use the Geneva rotating tracks?
        """
        self._tbirth = tbirth
        self._mass = mass
        self._metal = metal
        self._radiation = radiation
        self._wind = wind
        self._star = stars.Star(mass, metal, rotating=rotating)

    def Inject(self,injector):
        """
        Injects feedback from single stars over the timestep
        
        Parameters
        ----------

        injector : _Injector object
            object that accepts values to inject to grid
        """

        # Calculate the star's current age
        integrator = weltgeist.integrator.Integrator()
        t = integrator.time
        dt = integrator.dt
        age = t-self._tbirth
        Teff = 0.0
        star = self._star

        # If it's been born, start the tracks
        if age > 0.0:
            if self._wind:
                # Get mass loss and energy from winds
                ml0 = star.WindMassLoss(age)
                lum0 = star.WindLuminosity(age)
                massloss = 0.5*(ml0+star.WindMassLoss(age+dt))*dt
                energy =   0.5*(lum0+star.WindLuminosity(age+dt))*dt
                # Add mass FIRST since KE needs to be added elastically
                injector.AddMass(massloss)
                # Add energy to grid as kinetic energy
                injector.AddKE(energy)
                # Add some thermal energy to account for star's temperature
                Teff = 0.5*(star.Teff(age)+star.Teff(age+dt))
                TE = 1.5 * wunits.kB * massloss/(wunits.mH/wunits.X)*Teff
                injector.AddTE(TE)
                # Set the Courant condition
                vwind = np.sqrt(2.0*lum0/ml0)
                integrator.CourantLimiter(vwind)
            if self._radiation:
                # Get radiation properties for ionising and non-ionising UV
                Lionising = 0.5*(star.LIonising(age)+star.LIonising(age+dt))
                Lnonionising = 0.5*(star.LNonIonising(age)+star.LNonIonising(age+dt))
                Eionising = 0.5*(star.EPhoton(age)+star.EPhoton(age+dt))               
                Teff = 0.5*(star.Teff(age)+star.Teff(age+dt)) 
                Tion = weltgeist.radiation.IonisedGasTemperature(Teff, self._metal)
                injector.AddPhotons(Lionising, Lnonionising, Eionising, Tion)