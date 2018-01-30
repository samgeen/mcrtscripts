'''
Plot the solution for outflow radius
Sam Geen, May 2015
'''

import timeplot, outflowmodel

def run():
    # Set up timeplot
    timeplot.sims["N00_M4_B02"] = "k"
    timeplot.sims["N47_M4_B02"] = "b"
    timeplot.sims["N48_M4_B02"] = "r"
    timeplot.sims["N49_M4_B02"] = "c"
    timeplot.sims["N48_M4_B00"] = "m"
    timeplot.starts["N00_M4_B02"   ] = 1.25
    timeplot.starts["N47_M4_B02"   ] = 1.25
    timeplot.starts["N48_M4_B02"   ] = 1.25
    timeplot.starts["N48_M4_B00"   ] = 1.25
    timeplot.starts["N49_M4_B02"   ] = 1.25
    for simname in timeplot.sims.iterkeys():
        timeplot.sims[simname] = linestyles.col(simname)



if __name__=="__main__":
    run()
