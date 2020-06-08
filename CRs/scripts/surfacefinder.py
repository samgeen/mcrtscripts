'''
Find a surface fulfilling a certain criterion in a snapshot
Sam Geen, June 2018
'''

from startup import *

from pymses.filters import CellsToPoints
from pymses.utils import constants as C
from pymses.analysis import sample_points
from scipy.spatial import ConvexHull

def findsurface(snap,criterion="wind"):
    # Functions to determine which cells are in the surface
    def findwind(dset):
        temp = dset["P"]/dset["rho"]*snap.info["unit_temperature"].express(C.K) 
        return temp > 1e5
    def findHII(dset):
        return dset["xHII"] > 0.1
    # Get the cells
    amr = snap.amr_source(["rho","P","xHII","vel"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    mask = []
    if criterion == "wind":
        mask = findwind(cells)
    # Check if there are actually cells
    if len(mask) == 0:
        return None
    # Find surface
    points = cells.points
    return cells
    hull = ConvexHull(points[mask,:])
    # This is a list of indices of points on the surface
    surface = points[hull.vertices,:]
    # Can't be bothered to reverse mask to find the original list, so do this the longish way
    surfcells = sample_points(amr, surface)
    return surfcells


if __name__=="__main__":
    sim = hamusims["IMF2_04"]
    snap = sim.Snapshots()[39] # Get a random output with an evolved bubble
    surf = findsurface(snap.RawData(),"wind")
    import pdb; pdb.set_trace()
