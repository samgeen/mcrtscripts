"""
Read Teff tables
Sam Geen, April 2020
"""

import numpy as np
from scipy import interpolate

tables = {}

def moreatstartandend(x):
    # Gives more values at start and end
    # x must be 0 to 1
    return 0.5 - 0.5*np.cos(np.pi*x)

def Teff(age,starlifetime,starmass,starmetal,starrotating):
    """
    Get the effective temperature of a star at a given time
    NOTE: age and starmass must either be length 1, or equal length

    Parameters
    ----------

    age : float, numpy array
        stellar age

    starlifetime : float
        total stellar lifetime, same units as age
        (doesn't matter which units, we just divide one by the other)
    
    starmass : float, numpy array
        stellar mass in Msun

    starmetal : float
        metallicity of the star in absolute units (Sun = 0.014)

    starrotating : bool
        is the star rotating? True = 0.4x critical, False = 0.0x critical

    Return
    ------

    Teff : float, numpy array
        Effective temperature of the star
        Length is either 1, or same as age or starmetal if either or both is an array
    """
    global tables
    # Hard-coded properties used to make the tables
    masses = starmasses = np.arange(5,121,5) # every 5 Msun from 5 to 120 Msun
    nstars = len(masses)
    ntimes = 100000
    table = np.zeros((ntimes,nstars))
    times = np.linspace(0.0,1.0,ntimes) # range from 0.0 to 1.0
    times = moreatstartandend(times)
    # Read the table if necessary
    if not starmetal in tables:
        metalstr = str(starmetal)
        suffix = metalstr
        if not starrotating:
            suffix += "_nonrotating"
        filename = "Tefftable_Z"+suffix+".npy"
        try:
            # This reads the table from file
            tableIn = np.load(filename)
        except:
            # If we got here, we probably couldn't read the table
            print("Can't read file",filename,"does it exist?")
            raise FileNotFoundError
        # This makes an object that can interpolate across the table
        tables[starmetal] = interpolate.RectBivariateSpline(times,masses,tableIn)
    # Get the relevant tables and interpolate
    table = tables[starmetal]
    Teff = table(age / starlifetime, starmass)
    if len(Teff) == 1:
        Teff = Teff[0][0]
    return Teff


if __name__=="__main__":
    # Simple plotting routine to test the outputs
    import matplotlib.pyplot as plt
    import matplotlib
    test = Teff(0.5,1.0,120,0.014,True)
    times = np.linspace(0.0,1.0,200)
    Msolar = "M$_{\odot}$"
    c = 0.0
    cmap = matplotlib.cm.get_cmap('RdPu')
    for rotating in [True, False]:
        for m in [15,30,60,120]:
            c += 0.25
            rotline = {True: "-", False: "--"}
            rottxt = {True: ", rotating", False: ", non-rotating"}
            plt.plot(times,Teff(times,1.0,m,0.014,rotating), color=cmap(c),
                linestyle=rotline[rotating],
                label=str(m)+" "+Msolar+rottxt[rotating])
    plt.xlabel("Age / Stellar Lifetime")
    plt.ylabel("Effective temperature of star / K")
    plt.legend(frameon=False,fontsize="x-small")
    plt.savefig("testteff.pdf")