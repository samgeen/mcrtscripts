'''
Make a sphere map of whatever properties viewed from a central point
Sam Geen, June 2014
'''

import sys
import numpy as np
import matplotlib.pyplot as plt

# pymses gumf
from pymses.analysis import sample_points

totaljobs = 0.0
jobsdone = 0.0

def MollweideProj(x,y):
    '''
    Take x,y coordinates from 0->1 and return a lat and long in radians
    Return: Latitude, longitude (radians), 
       mask in lat, lon (true/false for valid pixels)
    c.f. http://en.wikipedia.org/wiki/Mollweide_projection
    '''
    theta2 = 2.0*np.arcsin(y/2**0.5)
    lat = np.arcsin((theta2 + np.sin(theta2))/np.pi)
    lon = np.pi * x / (2.0**1.5 * np.cos(theta2/2.0))
    print lat.min(), lat.max()
    print lon.min(), lon.max()
    pb2 = 0.5*np.pi
    latinrange = (lat >= -pb2) * (lat <= pb2)
    loninrange = (lon >=  -np.pi  ) * (lon <= np.pi)
    inan = np.logical_not(np.isfinite(lon)*np.isfinite(lat))
    mask  = (np.isfinite(lat)) * latinrange
    mask *= (np.isfinite(lon)) * loninrange
    mask = np.logical_not(latinrange*loninrange)
    lat[mask] = np.nan
    lon[mask] = np.nan
    return lat, lon

def OneRay(lray,radius):
    dr = radius*lray**(-1.0)
    oneray = np.arange(dr,radius+dr,dr)
    return oneray

def MakeRays(lat,lon,nray,radius):
    # Filter and repeat angle arrays
    good = np.isfinite(lat)
    lag = np.repeat(lat[good],nray)
    log = np.repeat(lon[good],nray)
    nangs = len(log)//nray
    # Set up ray array
    oneray = OneRay(nray,radius)
    oneray = np.tile(oneray,nangs)
    allrays = np.zeros((nangs*nray,3))
    # Now fill it
    allrays[:,0] = oneray*np.cos(log)*np.sin(lag)
    allrays[:,1] = oneray*np.sin(log)*np.sin(lag)
    allrays[:,2] = oneray*np.cos(lag)
    return allrays+0.5

def FuncOnRays(allrays, amr, func, nray):
    '''
    Run function on rays
    func should have args (sample_array, start, end)
    '''
    global jobsdone
    # Length of a ray
    numrays = allrays.shape[0]//nray
    samples = sample_points(amr,allrays)
    output = np.zeros(numrays)
    for i in np.arange(0,numrays*nray,nray):
        output[i//nray] = func(samples,i,i+nray)
        jobsdone += 1
        pc = jobsdone / totaljobs * 100.0
        sys.stdout.write("\r%d%%" %pc)
    return output 

def FuncOnLatLon(lat, lon, amr, func, nray, radius):
    '''
    Run function on rays
    func should have args (sample_array)
    '''
    global jobsdone
    # Length of a ray
    oneray = OneRay(nray,radius)
    ray = np.zeros((nray,3))
    # Now fill it
    ray[:,0 ] = oneray*np.cos(lon)*np.sin(lat)
    ray[:,1] = oneray*np.sin(lon)*np.sin(lat)
    ray[:,2] = oneray*np.cos(lat)
    samples = sample_points(amr,ray)
    output = func(samples)
    jobsdone += 1
    i = jobsdone / totaljobs * 100.0
    sys.stdout.write("\r%d%%" %i)
    return output 

def SphereMap(amr,func,nx=512,ny=256, nray=2e2,plotlog=False,radius=0.5):
    global totaljobs
    global jobsdone
    jobsdone = 0
    # Yes ny,nx (I'm about to start breaking everything in the room, pyplot)
    im = np.zeros((ny,nx))
    # Make an x,y value for each pixel
    x = (np.arange(0,nx)+0.5)/float(nx) - 0.5
    y = (np.arange(0,ny)+0.5)/float(ny) - 0.5
    x *= 2.0*np.pi
    y *= np.pi
    hemispheres = False
    if hemispheres:
        xin = np.repeat(x,ny)
        yin = np.tile(y,nx)
    else:
        xin = np.tile(x,ny)
        yin = np.repeat(y,nx)
    # Get latitude and longitude for each pixel in radians
    lat, lon = MollweideProj(xin,yin)
    # (test plot showing latitude)
    #lat /= np.pi; lat += 0.5
    #lon /= 2.0*np.pi;  lon += 0.5
    #for i in range(0,2):
    #    im = lat.reshape(ny,nx)
    # 
    # Make rays along each lat, lon
    # TODO: MOVE CENTRE, SET RADIUS LIMIT
    rays = MakeRays(lat,lon,nray,radius)
    # Make image array and normalise for plotting
    im = lat+0.0 # cheap copy
    good = np.isfinite(lat)
    lag = lat[good]
    log = lon[good]
    totaljobs = len(log)
    #im[good] = [FuncOnLatLon(a, o, amr, func, nray,radius) 
    #            for (a,o) in zip(lag,log)]
    im[good] = FuncOnRays(rays,amr,func,nray)
    im = im.reshape(ny,nx)
    if plotlog:
        im = np.log10(im)
    finim = im[np.isfinite(im)]
    try:
        imax, imin = finim.max(), finim.min()
    except:
        imax, imin = (0.0,0.0)
    #im = (im - imin) / (imax - imin)
    print imin, imax
    #import pdb; pdb.set_trace()
    # plot
    image = plt.imshow(im,origin="lower")
    plt.axis('off')
    plt.set_cmap('hot')
    return image, imin, imax

def MakeCBar(image, imrange,text):
    plt.colorbar(image, orientation='horizontal',label=text)

if __name__=="__main__":
    # NOTE - ONLY WORKS FOR CENTRE AT (0.5,0.5,0.5) SO FAR!
    import pymses
    ro = pymses.RamsesOutput("../runs/08_1000atcc_14lvl",19)
    amr = ro.amr_source(["rho"])
    def maxrho(samples,start,end):
        return samples["rho"][start:end].max()
    def rshell(samples,start,end):
        rhomax =  samples["rho"][start:end].max()
        imax = np.where(samples["rho"][start:end] == rhomax)
        radii = np.sum(samples.points[start:end,:]**2.0,axis=1)
        return np.sqrt(radii[imax])
    def cdens(samples,start,end):
        rhos =  samples["rho"][start:end]
        radii = np.sum(samples.points[start:end,:]**2.0,axis=1)
        lr = len(rhos)
        rhoc = 0.5*(rhos[0:lr-1]+rhos[1:lr])
        dr = np.abs(np.diff(radii))
        return np.sum((dr * rhoc)[rhoc > 1.5e3])
    # Max rho
        '''
    image, imax, imin = SphereMap(amr,maxrho,plotlog=True,radius=0.05)
    MakeCBar(image,(imin,imax),"log(Maximum Density)")
    plt.savefig("plots/maxrho.png",bbox_inches="tight")
    plt.clf()
    # Shell radius
    image, imax, imin = SphereMap(amr,rshell,plotlog=True,radius=0.05)
    MakeCBar(image,(imin,imax),"log(Shell radius)")
    plt.savefig("plots/shellradius.png",bbox_inches="tight")
    '''
    # Column density
    image, imax, imin = SphereMap(amr,cdens,plotlog=False,radius=0.05)
    MakeCBar(image,(imin,imax),"log(Column density)")
    plt.savefig("plots/columndensity.png",bbox_inches="tight")

