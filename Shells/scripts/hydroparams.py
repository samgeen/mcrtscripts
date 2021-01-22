'''
Find parameters in a hydro snapshot in Weltgeist
Sam Geen, January 2021
'''

import os, sys

sys.path.append("../../scripts") 
import HamuLite as Hamu

sys.path.append("/home/samgeen/Programming/Astro/Weltgeist")

import weltgeist

def findshellradius():
    """
    Find the radius of the shell
    """
    integrator = weltgeist.integrator.Integrator()
    integrator.hydro
    radius = hydro.x[np.where(hydro.vel > 0)]
    return radius


def runanalysis(filename=None):
    """
    Run parameter finder

    Parameters
    ----------

    filename: string
        Filename to read
        If None, read the current state
    """

    # Get integrator object
    integrator = weltgeist.integrator.Integrator()
    # Read save file if we have to
    if filename is not None:
        integrator.Load(filename)
    # Do analysis tasks



sim = Hamu.MakeSimulation("TEST","../outputs/test/")
for snap in sim.Snapshots():
    print(snap.Time())
