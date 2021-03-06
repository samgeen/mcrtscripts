'''
Stitch a movie together from images of each simulation
Sam Geen, August 2017
'''

import glob

from PIL import Image, ImageChops

from collections import OrderedDict

from startup import *

bbox = None

def trim(im):
    global bbox
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    if bbox is None:
        diff = ImageChops.difference(im, bg)
        diff = ImageChops.add(diff, diff, 2.0, -100)
        bbox = diff.getbbox()
    return im.crop(bbox)

def listimages(folder):
    imnames = glob.glob(folder+"/*.png")
    imnames.sort()
    return imnames

def makeimage(imnames,imnum,hydro):
    # Make list of images to plot
    # Assume that they're ordered 1,2,3,etc
    imlist = []
    for ims in imnames.values():
        ims.sort()
        if len(ims) <= imnum:
            currim = ims[-1]
        else:
            currim = ims[imnum-1] # Starts from 1
        imlist.append(currim)
    # Make a square image
    nx = int(np.sqrt(len(imlist)))+1
    ny = nx
    # Open the images
    bigim = None
    ix = 0
    iy = 0
    w = 0
    h = 0
    for imname in imlist:
        # Load and trim
        im = Image.open(imname)
        im.load()
        im = trim(im)
        if bigim is None:
            w,h = im.size
            bw = w*nx
            bh = h*ny
            bigim = Image.new('RGB',(bw,bh))
        bigim.paste(im,(ix,iy))
        # Cycle position
        ix += w
        if ix >= bw:
            ix = 0
            iy += h
    print "Saving image", imnum,"snip size",w,h
    numtxt = str(imnum).zfill(5)
    bigim.save("../plots/vis/bigim_"+hydro+"_"+numtxt+".png")

def run(hydro):
    simnames = allsims
    # Get list of image names
    imnames = OrderedDict()
    for simname in simnames:
        folder = glob.glob("../plots/"+simname+"/projections/"+hydro)[0]
        imnames[simname] = listimages(folder)

    nimages = [len(y) for y in imnames.values()]
    ntot = max(nimages)
    print ntot
    for iim in range(1,ntot+1):
        makeimage(imnames,iim,hydro)
    
if __name__=="__main__":
    run("rho")
