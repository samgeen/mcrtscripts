"""
Solve the hii region profiles like Draine 2011 or MartinezGonzalez+ 2014
Sam Geen, May 2021
"""

import os, sys, time

import numpy as np

import stars
import units

import matplotlib.pyplot as plt

import equilibriumtemperature
import ionisedtemperatures

import rdmfile

import analytics

gamma = equilibriumtemperature.gamma

Msun = "M$_{\odot}$"
cm3 = "cm$^{-3}$"

def integrateTrapezium(x,y):
    # Trapezium integration of function y in x
    a = 0.5*np.diff(x)*(y[1:]+y[:-1])
    return np.sum(np.concatenate(([0],a)))

def findnmean(ys,us,n0):
    # From Draine+ 2011
    # Modified to include ymin != 0
    ymin = ys.min()
    ymax = ys.max()
    vol = (ymax**3 - ymin**3) / 3.0 # factor 4 pi cancels below
    mean = n0 * integrateTrapezium(ys,ys*ys/us) / vol
    return mean

def findnrms(ys,us,n0):
    # From Draine+ 2011
    # Modified to include ymin != 0
    ymin = ys.min()
    ymax = ys.max()
    vol = (ymax**3 - ymin**3) / 3.0 # factor 4 pi cancels below
    rms = n0 * np.sqrt(integrateTrapezium(ys,ys*ys/(us*us)) / vol)
    return rms

def findshellmass(rs,ns,Omega):
    nsleft = ns[:-1]
    nsright = ns[1:]
    ravs = 0.5*(rs[:-1]+rs[1:])
    # Trapezium rule
    # = 0.5 * (nsleft + nsright) * dr
    f = Omega * units.mH/units.X * 0.5*(nsleft+nsright) * np.diff(rs) * ravs**2.0
    return np.sum(f)


def solvescalefree(uin,yin,gam,beta):
    '''
    We will do this in Draine 2011 scale-free units
    uin = value of u at inner edge
    yin = value of y at inner edge
    gam = scaling factor in Draine 2011
    beta = scaling factor in Draine 2011
    NOTES: 
    To solve, we need:
    phi(y=yin) = 1
    tau(y=yin) = 0
    phi(y=ymax) = 0
    u(y=yin) = defined by Eqn 4 plus the wind balance at r_wind
    '''
    phi = 1.0
    tau = 0.0
    u = uin+0.0
    y = yin+0.0
    phis = [phi]
    taus = [tau]
    us = [uin]
    ys = [yin]
    keepLooping = True
    # Loop through the steps
    dys = np.zeros(3)
    while keepLooping:
        # TODO: Update to RK4
        # Perform functions to calculate change in each function
        dudy = -1.0 - gam * (beta * np.exp(-tau) + phi) * u / (y*y)
        dphidy = -y*y/(u*u) - gam * phi/u
        dtaudy = gam/u
        # Pick a step that increases each by a small enough fraction
        fdy = 0.01
        dys[0] = np.abs(u/dudy)
        dys[1] = np.abs(phi/dphidy)
        dys[2] = np.abs(tau/dtaudy)
        dy = fdy * dys[dys > 0.0].min()
        #print(u, phi, tau, dudy, dphidy, dtaudy, dy)
        #import pdb; pdb.set_trace()
        # Trap the last step if we've nearly lost all our photons
        if phi + dphidy * dy < 1e-10:
            dy = -phi /dphidy
            keepLooping = False
        # Step
        u += dudy * dy
        phi += dphidy * dy
        tau += dtaudy * dy
        y += dy
        phis.append(phi)
        taus.append(tau)
        us.append(u)
        ys.append(y)
    # Return structure arrays
    phis = np.array(phis)
    taus = np.array(taus)
    us = np.array(us)
    ys = np.array(ys)
    return phis, taus, us, ys

def findprofile(star,cloud,rw,t):
    toplot = False
    # Collect values needed for equations
    alphaB = analytics.AlphaB(star,t)
    Tion = analytics.Tion(star,t)
    ci = analytics.SoundSpeed(star,t)
    Pw = analytics.WindPressure(star,cloud,rw,t)
    Lionising = star.LIonising(t)
    Lnonionising = star.LNonIonising(t)
    QH = star.NIonising(t)
    ephoton = star.EPhoton(t)
    # Get scale factors to produce scale-free values for calculation
    Qfac = QH / (4.0 * np.pi * alphaB)
    Tfac = 2.0 * units.c * units.kB * Tion / (alphaB * ephoton)
    n0 = Tfac*Tfac*Tfac / Qfac
    lambda0 = Qfac/(Tfac*Tfac)
    beta = Lnonionising/Lionising
    gam = Tfac * cloud.sigmaDust
    #beta = 0.0
    #gam = 0.0
    # Get values for n, r at inner edge from wind pressure balance
    nin = Pw * units.X / (units.mH * ci**2)
    # Convert inputs to scale-free units as in Draine 2011
    uin = n0 / nin
    yin = rw / lambda0
    # Find profile in scale-free units
    #import pdb; pdb.set_trace()
    phis, taus, us, ys = solvescalefree(uin,yin,gam,beta)
    # Rescale back to physical units
    ns = n0 / us
    rs = ys * lambda0
    # Get mean and rms density
    nmean = findnmean(ys,us,n0)
    nrms = findnrms(ys,us,n0)
    # Get shell mass
    mshell = findshellmass(rs,ns,cloud.Omega)
    #import pdb; pdb.set_trace()
    if toplot:
        print("PLOTTING")
        plt.clf()
        rs2 = rs+0
        rs2 -= rw
        rs2 /= rw
        plt.plot(rs2, ns,label="$n$")
        plt.xlabel("$(r - r_w) / r_w$")
        plt.ylabel("$n / $"+cm3)
        #plt.plot(rs, phis,label="$\\phi$")
        #plt.plot(rs, taus,label="$\\tau$")
        plt.scatter([rs[0]],[nin])
        plt.legend(fontsize="small",frameon=False)
        plt.yscale("log")
        plt.xscale("log")
        plt.savefig("../plots/testhiiprofile.pdf")
        print("DONE")
    return rs,ns,nmean,nrms,mshell

if __name__=="__main__":
    t = 1e5 * units.year
    starmass = 120.0
    star = stars.Star(starmass, 0.014, rotating=True)
    cloud = analytics.Cloud(3000,1.0*units.pc,2,1e-21)
    vw = analytics.dRdTforw2(star,cloud)
    rw = vw*t
    findprofile(star,cloud,rw,t)
    

