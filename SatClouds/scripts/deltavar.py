'''
Find the slope of the delta variance of the column density field
Sam Geen, January 2018
'''

from startup import *

from turbustat.statistics import DeltaVariance
from astropy.io import fits
from astropy import units

import columndensity, images

def _deltavarslope(snap,los='z',plotname=""):
    # Make column density map
    snap = snap.RawData()
    dmap = columndensity.DensityMap(snap,los)
    im = dmap.NH()
    # Make dummy FITS header for turbustat
    hdu = fits.PrimaryHDU(im)
    # Make lags by hand since something breaks if turbustat does it?
    # Possibly a floating point error where 512.0 > 1024/2.0 or something weird, anyway
    min_size = 3.0
    nlags = 25
    lags = np.logspace(np.log10(min_size),
                       np.log10(min(im.shape) / 2.), nlags) * 0.99999 * units.pix
    # Find the delta variance (see Ossenkopf et al. 2008)
    deltavar = DeltaVariance(im,header=hdu.header,lags=lags)
    deltavar.run(verbose=True,save_name="../plots/deltavar"+plotname+".pdf")
    import pdb; pdb.set_trace()

def runforsim(simname):
    print "Running for simulation", simname     
    sim = hamusims[simname]
    #tcode = tffcloud_code/2.0
    #tMyr = tffcloud_Myr/2.0
    #suffix = "0p50tff"
    tcode = tffcloud_code
    tMyr = tffcloud_Myr
    suffix = "1p00tff"
    snap = sim.FindAtTime(tcode)
    plt.clf()
    images.MakeFigure(simname,times=[tMyr],name="deltavarcompare"+suffix,los="z",shape=(1,1),nonamelabel=True)
    plt.clf()
    return _deltavarslope(snap,plotname=suffix)

def run(simnames,plotname):
    plt.clf()
    for simname in simnames:
        t, slope = runforsim(simname)


if __name__=="__main__":
    run(imfsims,"imf")
    run(icsims,"ic")
