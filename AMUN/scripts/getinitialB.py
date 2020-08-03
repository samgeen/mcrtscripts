'''
Get the initial peak B field
Sam Geen, July 2020
'''

from startup import *
import plotproperties

def run(simname):
    t, B = plotproperties.maxBfield(simname)
    print "Initial max B in simname", simname, "Bmaxini =", B[0], "microGauss"

if __name__=="__main__":
    run("NOFB")
    run("NOFB_DENSE")
