'''
Create archive .fits file for tracking Research Data Management of figures
Sam Geen, December 2019
'''

import numpy as np

from astropy.io import fits

class RDMFile(object):
    '''
    This is a file object used to create .fits archive files for Research Data Management
    Useage:
    rdm = RDMFile()
    rdm.AddPoints(np.random(20),np.random(20),label="XY Random Noise")
    rdm.AddArray(np.zeros((300,200)),label="Empty image")
    rdm.Write("../plots/mycoolrdmfile.fits")
    Note that you can run this in parallel with your plotting to make archive files
    We could also skim pyplot objects for points, etc, but that sounds like a lot of work
    '''

    def __init__(self):
        '''
        Constructor
        '''
        # List of items to put in the fits file
        self._itemlist = []
        # 
        self._entries = {}
        self._first = True

    def AddPoints(self,x,y,z=None,label=None):
        '''
        Add points series from line or scatter plots
        '''
        array = np.array([x,y,z])
        self._addHDU(array,label)

    def AddArray(self,array,label=None):
        '''
        Add a single array (image data or a cube)
        '''
        self._AddHDU(array,label)

    def Write(self,filename):
        '''
        Write the fits archive file
        If filename contains ".pdf" or ".eps", it will be replaced with ".fits"
        This allows plotting tools to just pass the figure name for parallel creation of archive files
        '''
        if ".pdf" in filename or ".eps" in filename:
            filename = filename.replace(".pdf",".fits")
        hdul = fits.HDUList(self._itemlist)
        hdul.writeto(filename)
    
    def _AddHDU(self,array,label=None):
        '''
        Private method that creates 
        '''
        if label is None:
            label = ""
        if self._first:
            hdu = fits.PrimaryHDU(array)
        else:
            hdu = fits.ImageHDU(array)
        hdr = hdu.header
        hdr["LABEL"] = label
        self._itemlist.append(hdu)