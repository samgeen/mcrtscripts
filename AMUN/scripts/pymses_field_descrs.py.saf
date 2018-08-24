# File descriptor in simulation
import pymses
from pymses.sources.ramses.output import *

myrtvars = []
ip = 0
ir = False

radbins = ["HII","HeII","HeIII"]
if ir:
    radbins = ["IR","Opt","HII","HeII","HeIII"]
for p in radbins:
    myrtvars.append(Scalar("Np"+p,ip))
    myrtvars.append(Vector("Fp"+p,[ip+1,ip+2,ip+3]))
    ip += 4

self.amr_field_descrs_by_file = \
    {
    "3D": {"hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2, 3]), 
                       Vector("B-left", [4, 5, 6]), 
                       Vector("B-right", [7, 8, 9]), 
                       Scalar("Pnontherm",10),
                       Scalar("P", 11), Scalar("Z",12),
                       Scalar("xHII",13), Scalar("xHeII",14), Scalar("xHeIII",15)],
           "rt"    : myrtvars
           }
    }
print "Read user field descriptors by Sam Geen"
#fields = pymses.sources.ramses.output.RamsesOutput.amr_field_descrs_by_file
#hs = fields["3D"]["hydro"]
#for h in hs:
#    print h.name
