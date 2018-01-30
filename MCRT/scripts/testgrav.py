'''
Test reading gravity
Sam Geen, June 2015
'''


import pymses
from pymses.sources.ramses.output import *

import Hamu

sim = Hamu.Simulation("N48_M4_B02")
snap = sim.Snapshots()[5]

print "Reading data..."
ro = snap.RawData()


ro.amr_field_descrs_by_file = \
    {
    "3D": {"hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2, 3]), 
                       Vector("B-left", [4, 5, 6]), 
                       Vector("B-right", [7, 8, 9]), 
                       Scalar("P", 10),
                       Scalar("xHII",11), 
                       Scalar("xHeII",12), Scalar("xHeIII",13)],
           "grav" : [Scalar("potential",0),Vector("force",[1,2,3])]}
    }

region = pymses.utils.regions.Sphere([0.5,0.5,0.5],0.1)

amr = ro.amr_source(["rho","potential"])
filt_amr = pymses.filters.RegionFilter(region,amr)
cells = pymses.filters.CellsToPoints(filt_amr).flatten()
print "MINMAX POTENTIAL", cells["potential"].min(),cells["potential"].max()

print "Done!"
