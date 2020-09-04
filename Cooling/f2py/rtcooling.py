# Calculates cooling rates from RAMSES RT's cooling module
# Sam Geen, August 2019

import numpy as np
import sys

from pymses.utils import constants as C 

import ramses
first = True

# Make sure this is the same as the Makefile and ramses.f90!
nvector = 500
ndim = 3
nIons = 3
nGroups = 5
T2buff = np.zeros(nvector)# + 1e7
xionbuff = np.zeros((nIons,nvector))#+0.999
Npbuff = np.zeros((nGroups,nvector))#+1e-30
Fpbuff = np.zeros((ndim,nGroups,nvector))
p_gasbuff = np.zeros((ndim,nvector))
nHbuff = np.zeros(nvector)# + 10000.0
Zsolarbuff = np.zeros(nvector)# + 1.0

def FinddTdt(T2,nH,xion,Zsolar,Np=None,Fp=None,p_gas=None,a_exp=np.array([1.0])):
    global first, nvector, ndim, nIons, nGroups, T2buff, xionbuff,Npbuff,Fpbuff,p_gasbuff,nHbuff,Zsolarbuff
    '''
    Rate of change in temperature in units of K/s
    # RT Cooling Fortran cheatsheet
    #   integer,parameter::nvectorIn=500
    #   integer,parameter::ndimIn=3
    #   integer,parameter::nIonsIn=3
    #   integer,parameter::nGroupsIn=5
    #   real(kind=8),intent(in),dimension(nvectorIn):: T2in
    #   real(kind=8),dimension(1:nIonsIn,1:nvectorIn),intent(in):: xionin
    #   real(kind=8),dimension(1:nGroupsIn,1:nvectorIn),intent(in):: Npin
    #   real(kind=8),dimension(1:ndimIn,1:nGroupsIn,1:nvectorIn),intent(in):: Fpin
    #   real(kind=8),dimension(1:ndimIn,1:nvectorIn),intent(in):: p_gasin
    #   real(kind=8),dimension(1:nvectorIn),intent(in):: nHin, Zsolarin
    '''
    # Only run startup once
    if first:
        ramses.data.startup()
        first = False
    # Number of buffer runs left
    nleft = T2.shape[0]
    # Initialise Np, Fp if empty
    if Np is None:
        Np = np.zeros((nGroups,nleft))
    if Fp is None:
        Fp = np.zeros((ndim,nGroups,nleft))
    if p_gas is None:
        p_gas = np.zeros((ndim,nleft))
    # 
    dt = np.array([3.1415e7]) # 1 year
    sys.stdout.flush()
    ndone = 0
    dTdt = T2*0.0
    # Loop on items of size nvector
    while nleft > 0:
        # Bookkeeping on how many values we're passing this round
        if nleft > nvector:
            ncell = nvector
        else:
            ncell = nleft
        nleft -= nvector
        sbuff = slice(0,ncell)
        sinput = slice(ndone,ndone+ncell)
        # Set up input buffers
        T2buff[sbuff] = T2[sinput]
        xionbuff[:,sbuff] = xion[:,sinput]
        Npbuff[:,sbuff] = Np[:,sinput]
        Fpbuff[:,:,sbuff] = Fp[:,:,sinput]
        p_gasbuff[:,sbuff] = p_gas[:,sinput]
        nHbuff[sbuff] = nH[sinput]
        Zsolarbuff[sbuff] = Zsolar[sinput]
        # Run the cooling!
        T2out, dNpdtout, dFpdtout = ramses.data.rt_cooling(T2buff, xionbuff, Npbuff, Fpbuff,
                                                     p_gasbuff,nHbuff, Zsolarbuff, dt, a_exp,ncell)
        ndone += ncell
        # Calculate temperature change rate
        dTdt[sinput] = (T2[sinput] - T2out[sbuff]) / dt
    return dTdt

def Finddedt(T2,nH,xion,Zsolar,Np=None,Fp=None,p_gas=None,a_exp=np.array([1.0])):
    '''
    Cooling "luminosity" density, i.e. rate of change in thermal energy per unit volume
    '''
    dTdt = FinddTdt(T2,nH,xion,Zsolar,Np=None,Fp=None,p_gas=None,a_exp=np.array([1.0]))
    
    gamma = 1.4 # RAMSES hard-coded
    X = 0.76
    mH = 1.6735326381e-24
    kB = 1.38062e-16
    # Calculate d(energy density)/dt based on change in temperature
    dedt = 1.0/(gamma - 1.0) * nH * kB * dTdt
    return dedt

def dedtOnCells(snap):
    '''
    Returns a function that acts on a list of cells (to fit the format of visualisation in Pymses)
    '''
    #print "Finding cooling luminosity in snap", snap.iout
    #import rtcooling
    #amr = snap.amr_source(["rho","P","vel","xHII","xHeII","xHeIII","NpHII","NpHeII","NpHeIII"])
    #cell_source = CellsToPoints(amr)
    #cells = cell_source.flatten()
    def coolfunc(cells):
        umass = snap.info["unit_mass"].express(C.g) 
        uvel = snap.info["unit_velocity"].express(C.cm/C.s)
        ue = umass*uvel**2
        nHs = cells["rho"]*snap.info["unit_density"].express(C.H_cc)
        toshape = False
        if len(nHs.shape) > 1:
            toshape = True
            newshape = nHs.shape
            nHs = nHs.flatten()
        vels = cells["vel"]
        #vols = (cells.get_sizes()*snap.info["unit_length"].express(C.cm))**3.0
        spds = np.sqrt(np.sum(vels**2.0,1))
        mufunc = lambda dset: 1./(0.76*(1.+dset["xHII"].flatten()) + \
                                    0.25*0.24*(1.+dset["xHeII"].flatten()+2.*dset["xHeIII"].flatten()))  
        temp = (cells["P"].flatten())/(cells["rho"].flatten())*snap.info["unit_temperature"].express(C.K)*mufunc(cells)
        xion = np.array([cells["xHII"].flatten(),cells["xHeII"].flatten(),cells["xHeIII"].flatten()])
        Zsolar = nHs*0.0 + 1.0
        Nps = np.array([Zsolar*0.0,Zsolar*0.0,cells["NpHII"].flatten(),cells["NpHeII"].flatten(),cells["NpHeIII"].flatten()])
        dedtpervol = Finddedt(temp,nHs,xion,Zsolar,Np=Nps,Fp=None,p_gas=None,a_exp=np.array([1.0]))
        Lcool = dedtpervol
        if toshape:
            Lcool = np.reshape(Lcool,newshape)
        return Lcool
    return coolfunc

if __name__=="__main__":
    ntest = 950
    T2 = np.zeros(ntest)+1e3
    nH = np.zeros(ntest)+1e3
    xion = np.zeros((3,ntest))+0.999
    Zsolar = np.zeros(ntest)+1.0
    print("Change in temperature in K/s:", run(T2,nH,xion,Zsolar))
