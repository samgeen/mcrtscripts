'''
Take a quantity from a list of stellar masses and produce an average for a pop
Sam Geen, November 2016
'''

import numpy as np
import scipy.interpolate

import chabrier, readstarburst, HamuFunc

masses = np.array(readstarburst.allmasses,dtype="float32")[::-1]

def Run(massesin,valuesin):
    '''
    Find the average value per solar mass in a population
    masses - array of sample masses
    values - array of sample values for each mass
    '''
    # Find IMF
    # mcen - mass in centre of each bin
    # imfcen - dN/dM for each bin
    # N - number-averaged fraction of mass in each bin
    # M - mass-averaged fraction in each bin
    mcen, imfcen, N, M = chabrier.makeimf()
    newvalues = np.interp(mcen,massesin,valuesin)
    popvalue = np.sum(newvalues * N)
    return popvalue

def _MakeGrid(time=1e6):
    # Make grid of spectra
    allspec = None
    imass = -1
    for mass in masses:
        print mass
        imass += 1
        spec = readstarburst.Track(mass).SpectrumAtTime(time)
        fluxes = spec.fluxes
        wls = spec.wls
        if allspec is None:
            allspec = np.zeros((fluxes.shape[0],masses.shape[0]))
        allspec[:,imass] = fluxes
    return allspec
MakeGrid = HamuFunc.Algorithm(_MakeGrid)

def _FindSpectrum(mass,time):
    return readstarburst.Track(mass).SpectrumAtTime(time)
FindSpectrum = HamuFunc.Algorithm(_FindSpectrum)

def _MakeSpectrum(time=1e6):
    '''
    Make a spectrum from Starburst 99
    time: time in years - NOTE - CHOOSES THE NEAREST ONE TODO CHANGE THIS
    '''
    # Find IMF
    # mcen - mass in centre of each bin
    # imfcen - dN/dM for each bin
    # N - number-averaged fraction of mass in each bin
    # M - mass-averaged fraction in each bin
    mcen, imfcen, N, M = chabrier.makeimf()
    # Make a grid of values
    grid = MakeGrid(time=time)
    # Interpolate grid to find the masses in the IMF
    wls = FindSpectrum(masses[0],time).wls
    f = scipy.interpolate.interp2d(wls,masses,grid.T,kind="linear")
    # Make a spectrum for each mcen, normalised by N
    spectrum = np.sum(f(wls,mcen)*N[:,None],axis=0)
    print spectrum.shape
    # Dun
    return wls, spectrum
    
MakeSpectrum = HamuFunc.Algorithm(_MakeSpectrum)

if __name__=="__main__":
    fluxes, spec =  MakeSpectrum(time=1e6)
    import matplotlib.pyplot as plt
    plt.clf()
    plt.plot(fluxes,spec)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
