# Calculates cooling rates from RAMSES RT's cooling module
# Sam Geen, August 2019

import numpy as np
import sys

import ramses
first = True

def run():
    global first
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
    if first:
        ramses.data.startup()
        first = False
    nvector = 500
    ndim = 3
    nIons = 3
    nGroups = 5
    T2 = np.zeros(nvector) + 1e7
    xion = np.zeros((nIons,nvector))+0.999
    Np = np.zeros((nGroups,nvector))+1e-30
    Fp = np.zeros((ndim,nGroups,nvector))
    p_gas = np.zeros((ndim,nvector))
    nH = np.zeros(nvector) + 10000.0
    Zsolar = np.zeros(nvector) + 1.0
    dt = np.array([3.1415e7]) # 1 year
    a_exp = np.array([1.0])
    ncell = 500
    print "Running RT cooling"
    sys.stdout.flush()
    T2out, dNpdtout, dFpdtout = ramses.data.rt_cooling(T2, xion, Np, Fp, p_gas,nH, Zsolar, dt, a_exp,ncell)
    print T2out, dNpdtout, dFpdtout

if __name__=="__main__":
    run()