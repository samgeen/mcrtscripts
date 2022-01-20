'''
MODIFIED FROM Algorithm.py IN HAMU TO WORK INDEPENDENTLY OF THE HAMU FRAMEWORK
Created on 18 Feb 2013

@author: samgeen
'''

# Cache file imports
import hashlib

import sys
if sys.version_info[0] < 3:
    import cPickle as pik
else:
    import pickle as pik

import os, glob

try:
    import pymses
except:
    print("HamuLite: pymses not found, not importing")

try:
    sys.path.append("/home/samgeen/Programming/Astro/Weltgeist")
    import weltgeist
except:
    print("HamuLite: Weltgeist not found, not importing")

import numpy as np

from collections import OrderedDict

# Used for cases where folders have been tarred
UNTARPATH = None

# Algorithm output cache location
CACHEPATH = None

# Force replacement of the cache regardless of whether one exists or not
GLOBALFORCEREPLACECACHE = False

# Error code to place in ignored error outputs
ERRORCODE = None

# Ignore errors from outputs not found?
IGNOREERRORS = False

# Class thats wrap a function around a possible error 
class ErrFunc(object):
    def __init__(self,func):
        self._func = func

    def __call__(self,snap,*args,**kwargs):
        # Call the function and check it works ok
        try:
            ret = func(snap,*args,**kwargs)
        except ValueError:
            print("ERROR FOUND, IGNORING OUTPUT...")
            return ERRORCODE

class CacheFile(object):
    def __init__(self, snapshot, algorithm):
        self._snapshot = snapshot
        self._algorithm = algorithm
        self._folder = snapshot.CachePath()
    
    def Save(self, data):
        # Save algorithm settings to a text file (for reference)
        pikfile = open(self._Filename("info"),"wb")
        pik.dump(str(self._algorithm),pikfile)
        pikfile.close()
        # Save the output data to a binary file
        pikfile = open(self._Filename("data"),"wb")
        pik.dump(data,pikfile)
        pikfile.close()
        
    def Load(self):
        # Load the (binary) data file
        if self.Exists():
            pikfile = open(self._Filename("data"),"rb")
            print("Loading from cache: snap", self._snapshot.OutputNumber(),\
                self._algorithm.FunctionName(), "...")
            try:
                output = pik.load(pikfile)
            except:
                print("Hamu: Error reading cache file", pikfile)
                raise ValueError
            pikfile.close()
            return output
        else:
            return None
        
    def Exists(self):
        '''
        Does the cache file exist?
        '''
        return os.path.exists(self._Filename())
        
    def _Filename(self,ext="data"):
        '''
        Cache file's filename
        '''
        return self._folder+self._algorithm.CacheFilename(ext)

class Algorithm(object):
    '''
    A wrapper class around a function that enables us to run functions and store their outputs for later use
    '''

    def __init__(self, function, *args, **kwargs):
        '''
        Constructor
        function: A Python function object to call that accepts snapshot.RawData() as its first argument, and arg/kwarg after
        arg, kwarg: arguments accepted by function (see, e.g., http://www.saltycrane.com/blog/2008/01/how-to-use-args-and-kwargs-in-python/)
        ''' 
        self._function = function
        self._args = args
        self._kwargs = kwargs
        self._force_replace_cache = False

    def ForceReplaceCache(self,on):
        self._force_replace_cache = on

    def FunctionName(self):
        '''
        Return the function's name
        '''
        return self._function.__name__
        
    def Run(self, snap, *args, **kwargs):
        '''
        Run for a single snapshot
        '''
        global GLOBALFORCEREPLACECACHE, IGNOREERRORS
        self._args = args
        self._kwargs = kwargs
        # First, get the cache filename and compare against existing files
        cache = CacheFile(snap, self)
        # Check to see if a cached dataset exists
        if cache.Exists() and not self._force_replace_cache and not GLOBALFORCEREPLACECACHE:
            try:
                output = cache.Load()
            except EOFError:
                # If the cache is corrupted, rerun anyway
                output = self._RunAlgorithm(snap)
                cache.Save(output)
        else:
            if not IGNOREERRORS:
                output = self._RunAlgorithm(snap)
                cache.Save(output)
            else:
                try:
                    output = self._RunAlgorithm(snap)
                    cache.Save(output)
                except ValueError:
                    print("Error loading file", snap.Folder(), "ignoring...")
                return ERRORCODE
        return output

    def Cached(self, snap, *args, **kwargs):
        '''
        Check that the cache for this data exists
        '''
        self._args = args
        self._kwargs = kwargs
        # First, get the cache filename and compare against existing files
        cache = CacheFile(snap, self)
        # Check to see if a cached dataset exists
        return cache.Exists()
        
    def _RunAlgorithm(self, snapshot):
        '''
        Unpack the algorithm and call the native python code
        '''
        raw = snapshot.RawData()
        output = self._function(raw, *self._args, **self._kwargs)
        return output
    
    def __str__(self):
        '''
        Parse this algorithm as a string (for writing to text file)
        '''
        out = ""
        out += "Function: "+self._function.__name__
        out += "\n"
        out += "Arguments: "+str(self._args)
        out += "\n"
        out += "Keyword Arguments: "+str(self._kwargs)
        out += "\n"
        return out
    
    def __call__(self, snapshot, *args, **kwargs):
        '''
        Allows the algorithm to be called like a function
        '''
        return self.Run(snapshot, *args, **kwargs)
    
    def CacheFilename(self, ext="data"):
        '''
        Return the name of this algorithm's cache file
        Uses a hash function to hash the arguments of the function
        ext - File extension (used for writing multiple cache files; default is ".data")
        '''
        objName = self._function.__name__
        hash = hashlib.sha1(str(self._args).encode('utf-8')).hexdigest()
        hash += hashlib.sha1(str(self._kwargs).encode('utf-8')).hexdigest()
        filepre = objName+hash
        return filepre+"."+ext

def _PymsesCacheTimeHamu(snap):
    return snap.info["time"]
CacheTime = Algorithm(_PymsesCacheTimeHamu)

class PymsesSnapshot(object):
    def __init__(self,folder,outnum,name=None):
        self._snap = None
        self._folder = folder
        self._outnum = outnum
        self._name = name

    def _Setup(self):
        if self._snap is None:
            self._snap = pymses.RamsesOutput(self._folder,self._outnum)
            # Hack to allow nested Algorithm calls
            self._snap.hamusnap = self 

    def RawData(self):
        self._Setup()
        return self._snap

    def Name(self):
        return self._name

    def CachePath(self):
        #if self._name is None:
        #    pathname = self._snap.output_repos
        #else:
        #    pathname = self._name
        #pathname = pathname.replace("/","__").replace("~","TILDE")
        #path = "./cache/"+pathname+"/output_"+str(self.OutputNumber()).zfill(5)+"/"
        if CACHEPATH is not None:
            pathname = CACHEPATH+"/"+self._folder
        else:
            pathname = self._folder
        # Put the cache files in the simulation folder itself
        path = pathname+"/cache/output_"+str(self.OutputNumber()).zfill(5)+"/"
        if not os.path.exists(path):
            try:
                os.makedirs(path)
            except:
                pass # Probably fine. Probably.
        return path

    def OutputNumber(self):
        return self._outnum

    def Time(self):
        '''
        Return the output time (for comparing outputs)
        TODO: Make this concept more concrete (i.e. make sure units/measurement methods match)
        '''
        return CacheTime(self) #self._snapshot.info["time"] 

    def Folder(self):
        return self._folder
                          
class WeltgeistSnapshot(object):
    def __init__(self,folder,filename,name=None):
        self._snap = None
        self._filename = filename
        self._folder = folder
        self._outnum = int(filename[-10:-5])
        self._name = name
        self._time = None

    def _Setup(self):
        if self._time is None:
            integrator = weltgeist.integrator.Integrator()
            integrator.Load(self._filename)
            self._time = integrator.time

    def RawData(self):
        integrator = weltgeist.integrator.Integrator()
        integrator.Load(self._filename)
        return integrator

    def Name(self):
        return self._name

    def CachePath(self):
        #if self._name is None:
        #    pathname = self._snap.output_repos
        #else:
        #    pathname = self._name
        #pathname = pathname.replace("/","__").replace("~","TILDE")
        #path = "./cache/"+pathname+"/output_"+str(self.OutputNumber()).zfill(5)+"/"
        if CACHEPATH is not None:
            pathname = CACHEPATH+"/"+self._folder
        else:
            pathname = self._folder
        # Put the cache files in the simulation folder itself
        path = pathname+"/cache/output_"+str(self.OutputNumber()).zfill(5)+"/"
        if not os.path.exists(path):
            try:
                os.makedirs(path)
            except:
                pass # Probably fine. Probably.
        return path

    def OutputNumber(self):
        return self._outnum

    def Time(self):
        '''
        Return the output time (for comparing outputs)
        TODO: Make this concept more concrete (i.e. make sure units/measurement methods match)
        '''
        self._Setup()
        return self._time

    def Folder(self):
        '''
        Return the simulation folder
        '''
        return self._folder


class _CodeFactory(object):
    """
    Code-specific factory class for making outputs
    """
    def __init__(self,folder,name):
        self._folder = folder
        self._name = name

    def Outputs(self):
        global UNTARPATH
        outputs = OrderedDict()
        pymseslist = glob.glob(self._folder+"/output_?????")
        pymseslist2 = glob.glob(self._folder+"/output_?????.tar")
        weltlist = glob.glob(self._folder+"/*.hdf5")
        # Process pymses (RAMSES) outputs
        if len(pymseslist) or len(pymseslist2) > 0:
            print("Making Pymses simulation...")
            mergedlist = pymseslist+pymseslist2
            mergedlist.sort()
            # Loop through list of folders and tar files
            for out in mergedlist:
                simfolder = self._folder
                outfolder = out
                if out[-4:] == ".tar":
                    outfolder = out[:-4]
                    # Check for existing mount point and mount if needed
                    if not os.path.exists(outfolder):
                        # Check if we need to untar to a specific place
                        # e.g. sometimes drives don't like being fusermounted to
                        print(UNTARPATH)
                        if UNTARPATH is not None:
                            outfolder = UNTARPATH+"/"+outfolder
                        # Check again for the new untar folder in case it's already been untarred
                        if not os.path.exists(outfolder):
                            try:
                                os.makedirs(outfolder)
                            except:
                                pass # Probably fine. Probably.
                            syscommand = "ratarmount "+out+" "+outfolder
                            print("Mounting tar file as: "+syscommand)
                        #import pdb; pdb.set_trace()
                        #os.system(syscommand)
                        simfolder = outfolder[:-len("output_?????")] 
                outnum = int(outfolder[-5:])
                outputs[outnum] = PymsesSnapshot(simfolder,outnum,self._name)
        # Process weltgeist outputs
        elif len(weltlist) > 0:
            print("Making Weltgeist simulation...")
            weltlist.sort()
            for out in weltlist:
                outnum = int(out[-10:-5])
                outputs[outnum] = WeltgeistSnapshot(self._folder,out,self._name)
        return outputs

# Singleton pattern to store simulations
folders = {}
labels  = {}

class Simulation(object):
    def __init__(self,name,silent=False):
        global folders, labels 
        self._name = name
        if not name in folders:
            if not silent:
                print("No simulation with this name stored in memory:", name)
                print("Currently available simulations:", folders.keys())
            raise KeyError
        self._folder = folders[name]
        self._label  = labels[name]
        self._snaps   = _CodeFactory(self._folder,name).Outputs()

    def Snapshots(self):
        return list(self._snaps.values())
                    
    def Name(self):
        return self._name
        
    def Label(self,label=None):
        global labels
        if label is not None:
            labels[self._name] = label
            self._label = label
        if self._label is None:
            return self._folder
        return self._label

    def Times(self):
        return np.array([snap.Time() for snap in self.Snapshots()])
    
    def Outputs(self,outnum=None):
        return self._snaps.keys()

    def FindAtTime(self,time,forceBefore=False):
        times = self.Times()
        diff = np.abs(times-time)
        best = np.where(diff == np.min(diff))
        best = best[0][0]
        # Force the chosen output to be before "time"?
        if times[best]-time > 0.0 and best != 0 and forceBefore:
            best -= 1
        pdiff = (times[best]-time) / time * 100.0 # Is a percentage, so multiply by 100
        snap = self.Snapshots()[best]
        print("Found match with output",snap.OutputNumber(),", %diff: ",pdiff, "at time ",snap.Time(),"in",self.Name())
        return snap
    
    def FindOutput(self,output):
        if output in self._snaps:
            print("Found output number", output)
            return self._snaps[output]
        else:
            print("Found no output number", output)
            return None

    def Folder(self):
        return self._folder
                                                                                            
def MakeSimulation(name,folder,label=None):
    global folders, labels
    folders[name] = folder
    labels[name] = label
    return Simulation(name)
