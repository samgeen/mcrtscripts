'''
Created on 9 Nov 2014

@author: samgeen
'''
'''
Makes a mesh on a surface using a spiral
NOTE THIS ONE WORKS THE OTHERS DO NOT
Created on 8 Nov 2014
@author: samgeen
'''
import numpy as np
from scipy.spatial import ConvexHull

def Magnitude(points):
    '''
    Pop pop
    '''
    dist = np.sqrt(np.sum(points**2,1))
    for j in range(0,3):
        points[:,j] /= dist

class Mesh(object):
    
    def __init__(self, npoints):
        self._npoints = int(npoints)
        self._points = None
        self._nbors = None
        self._Setup()
        
    def Vertices(self):
        return self._points
    
    def Num(self):
        return self._npoints
    
    def Neighbours(self):
        return self._nbors
        
    def _Setup(self):
        num = self._npoints
        # Algorithm shamelessly nabbed from http://people.sc.fsu.edu/~jburkardt/f_src/sphere_grid/sphere_grid.html
        inds = np.arange(0,self._npoints)+1.0
        cosphi = ((2.0*inds-num-1.0))/(num-1.0)
        sinphi = np.sqrt(1.0 - cosphi*cosphi)
        theta = 3.6 / (sinphi*np.sqrt(num))
        theta[sinphi == 0.0] = 0
        theta = np.cumsum(theta)
        theta = np.mod(theta, 2.0*np.pi)
        theta[0] = 0.0
        theta[num-1] = 0.0
        self._points = np.array([sinphi*np.cos(theta),
                                 sinphi*np.sin(theta),
                                 cosphi]).T
        # Get neighbours
        hull = ConvexHull(self._points)
        nbors = hull.simplices
        #print nbors
        flat = np.sort(nbors.flatten())
        #self._nbors = hull
        print "Running mesh spacer"
        print nbors.max()
        #print self._npoints, nbors.shape
        self._nbors = list()
        for ivert in range(0,self._npoints):
            #print ivert
            #print nbors[np.where(nbors == ivert)[0],:]
            currnbors = np.unique(np.sort(nbors[np.where(nbors == ivert)[0],:].flatten()))
            currnbors = currnbors[currnbors != ivert]
            self._nbors.append(currnbors)
            
    def MakeHeightmap(self):
        smooth = 0.1
        print "Making heights"
        self._heights = np.random.rand(self._npoints)
        for iloop in range(0,10):
            for ivert in range(0,self._npoints):
                nbors = self._nbors[ivert]
                self._heights[ivert] = self._heights[ivert]*(1-smooth) + \
                    smooth*np.sum(self._heights[nbors])/(1.0+len(nbors))
                #self._heights[ivert] /= (1.0+len(nbors))
            print self._heights.max()
        pts = (self._points.T*(4.0*self._heights+0.0)).T
        #print pts
        return pts
        
        
        
def Draw():
    print "Drawing..."
    mesh = Mesh(300)  
    heights = mesh.MakeHeightmap()
    # Import these locally in case user doesn't have matplotlib installed
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.colors import colorConverter
    import matplotlib.pyplot as plt
    
    toplot = []
    verts = heights +  0.5 # Centre on matplotlib view
    z = np.zeros((mesh.Num()))
    '''
    for i in range(0,mesh.Num()):
        curr = np.zeros((3,3))
        centre = verts[i,:]
        curr[0,:] = verts[nbors[i,0]]
        curr[1,:] = verts[nbors[i,1]]
        curr[2,:] = verts[nbors[i,2]]
        #for j in range(0,3):
        #    curr[j,:] = 0.5*(verts[nbors[i,j]] + centre)
        #print centre
        #print verts[nbors[i,:]]
        #print np.sum(verts[nbors[i,:]]**2,1)
        #print "---"
        toplot.append(curr)
    '''
    #coll = Poly3DCollection(toplot,array=z)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    scat = ax.scatter(verts[:,0],verts[:,1],verts[:,2])
    #ax.add_collection3d(coll,zs=z)
    #ax.plot_trisurf(x,y,z,color=z)
    ax.autoscale_view()
    #fig.colorbar(coll,ax=ax)
    plt.show()
                
            
            
            
        
if __name__=="__main__":
    Draw()
