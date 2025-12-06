#!/usr/bin/env python

###########################
##       flux.py         ##
###########################

import numpy as np
from . import axes as axes
from . import fdrule as fd

def get_sampling_pes(time,signal):

    dt = time[1] - time[0]
    ww = np.fft.fftfreq(len(time),dt) * 2.0 * np.pi
    indx = np.where(ww>=0.0)
    ww = ww[indx]
    cpes = np.fft.ifft(signal) * dt / np.sqrt(2.0 * np.pi)
    pes = np.abs(cpes[indx])**2 * np.sqrt(ww)
    return ww, pes
    
def get_surface(psi,ax,xb,fdrule=2):
    indx = np.argmin(np.abs(ax.ax<xb))
    surf = psi[indx]

    rulepts = 2 * fdrule + 1
    maxfd = fdrule * ax.dx
    fdax = np.linspace(-maxfd,maxfd,rulepts)
    fdcoeff = fd.calculate_fdweights(0.0,fdax,2)
    
    ii = indx - fdrule
    ij = indx + fdrule + 1
    deriv_surf = sum(fdcoeff[:,2] * psi[ii:ij])
    return surf, deriv_surf

def get_tsurff(surf,deriv_surf,ax,xb,time,field,fdrule=2):
    dt = time[1]- time[0]
    ntime = len(time)
    kax = axes.kaxes(ax)
    bk = np.zeros(kax.maxptsx,dtype=complex)
    vphase = np.ones(kax.maxptsx,dtype=complex)
    
    for itime in range(ntime):
        intflux = np.zeros(kax.maxptsx,dtype=complex)
        intflux = calculate_tsurff_integrand(surf[itime],
                                             deriv_surf[itime],
                                             kax,xb,field[itime,0])
        vphase = get_volkov_phase(vphase,kax,field[itime,0],dt)
        
        intflux *= vphase
        bk += intflux 
    
    bk *= 1j * dt / np.sqrt(2.0*np.pi)
    return bk
    
def calculate_tsurff_integrand(surf,deriv_surf,
                               kax,xb,afield):
    integrand = -0.5 * 1j * kax.ax * surf
    integrand -= deriv_surf
    integrand -= 1j * afield * surf
    integrand *= np.exp(-1j*kax.ax*xb)
    return integrand
    

def get_volkov_phase(vphase,kax,afield,dt):
    term1 = (kax.ax - afield )**2
    newphase = np.exp(1j * 0.5 * term1 * dt)
    newphase *= vphase
    return newphase
    
