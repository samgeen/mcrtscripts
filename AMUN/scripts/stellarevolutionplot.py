from startup import *

def plotformass(mass):
    agesins = np.arange(0,10.0,0.02)*1e6*yrins
    dtins = 1.0*yrins
    energymass = singlestar.star_winds(mass,agesins,dtins)
    print energymass.shape

def run():
    masses = [30.0,60.0,120.0]
    for mass in masses:
        plotformass(mass)

if __name__=="__main__":
    run()
    