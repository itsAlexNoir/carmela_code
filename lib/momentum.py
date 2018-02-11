#!/usr/bin/env python
# -*- coding: utf-8 -*-


########################################
############### momentum ###############
########################################

import numpy as np
from . import axes as axes
from . import wavefunction as wave

def create_mask(axis,inner_rad,outer_rad):
    rad = np.sqrt(axis.ax**2)
    sigma = (outer_rad - inner_rad) / np.sqrt(-np.log(1e-8))
    
    mask = np.zeros(axis.maxptsx)
    mask[np.where(rad<inner_rad)] = 0.0
    mask[np.where(rad>=inner_rad)] = \
             1.0 - np.exp(-((rad-inner_rad)/sigma)**2)
    mask[np.where(rad>=outer_rad)] = 1.0
    return mask

def get_momentum_amplitude(axis,wavefunction,mask_radius):
    mask = create_mask(axis,mask_radius[0],mask_radius[1])
    kax = axes.kaxes(axis)
    psik = wave.wavefunction(kax)
    
    psifft = wavefunction.psi * mask
    psik.psi = np.fft.fft(psifft) * axis.dx / np.sqrt(2.0*np.pi)
    psik.psi = np.fft.fftshift(psik.psi)
    return kax,psik
