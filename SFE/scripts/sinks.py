'''
Find the sinks in a snapshot
Sam Geen, June 2016
'''

from startup import *
from pymses.utils import constants as C

class SinkSnap(object):
    # A container for sink properties
    def __init__(self,i=[],x=[],y=[],z=[],mass=[],time=[]):
        self.id = np.atleast_1d(i)
        self.id = self.id.astype(np.int64)
        self.x = np.atleast_1d(x)
        self.y = np.atleast_1d(y)
        self.z = np.atleast_1d(z)
        self.mass = np.atleast_1d(mass)
        self.time = np.atleast_1d(time)
        self._dict = None
        self._Setup()

    def _Setup(self):
        self._dict = {"id":self.id,
                      "x": self.x,
                      "y": self.y,
                      "z": self.z,
                      "mass": self.mass,
                      "time": self.time}

    def __getitem__(self,key):
        return self._dict[key]

    def IsEmpty(self):
        return self.Length() == 0

    def Length(self):
        return len(np.atleast_1d(self.mass))
        

def FindSinks(snap):
    # Don't really need to cache this with Hamu as it's fast
    ro = snap.RawData()
    onum = str(ro.iout).zfill(5)
    sinkfile = ro.output_repos+"/output_"+onum+"/sink_"+onum+".csv"
    time = snap.Time()*ro.info["unit_time"].express(C.Myr)
    # Particle mass units don't include boxlen (thanks RAMSES)
    unitm = ro.info["unit_mass"].express(C.g)/ro.info["boxlen"]**3.0
    unitm /= Msuning
    praw = []
    if os.path.exists(sinkfile):
        if os.stat(sinkfile).st_size > 0:
            praw = np.genfromtxt(sinkfile, delimiter=',')
    pid = []
    pmass = []
    px = []
    py = []
    pz = []
    if len(praw) > 0:
        if len(praw.shape) == 1:
            pid = np.array(praw[0])
            pmass = np.array(praw[1])
            px = np.array([praw[2]])
            py = np.array([praw[3]])
            pz = np.array([praw[4]])
        else:
            pid = praw[:,0]
            pmass = praw[:,1]
            px = praw[:,2]
            py = praw[:,3]
            pz = praw[:,4]
        sinks = SinkSnap(pid,px,py,pz,pmass*unitm,time)
    else:
        sinks = SinkSnap() # Make an empty object
    return sinks
                
