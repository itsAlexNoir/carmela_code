#!/usr/bin/env python

import numpy as np
from . import constants as const

#############################
##        laser.py
############################

class laser():
    def __init__(self,dt,intensity,wavelength,no_cycles,afterpulse=0,cep=0.0):

        self.dt        = dt
        self.intensity = intensity
        self.wavelength = wavelength
        self.no_cycles  = no_cycles
        self.cep       = cep
        self.afterpulse = afterpulse
        
        self.wavelength0 = self.wavelength / const.aulength_nm
        self.w0 = 2.0 * np.pi * const.speed_light_au / self.wavelength0
        self.period = 2.0 * np.pi / self.w0
        self.pulse_duration = self.period * self.no_cycles
        if(self.pulse_duration!=0):
            self.pulse_bandwidth = 4.0 * np.pi / self.pulse_duration
        else:
            self.pulse_bandwidth = 0.0
        self.pulse_bandwidtheV = self.pulse_bandwidth * const.energy_au_ev
        self.quiver = np.sqrt(self.intensity / const.intensity_au) / self.w0 / self.w0
        self.Up = (self.intensity / const.intensity_au) / 4. / self.w0 / self.w0
        self.E0 = np.sqrt(self.intensity / const.intensity_au)
        self.A0 = self.E0 / self.w0 

        self.afterpulsetime = self.afterpulse * self.period
        self.interactiontime = self.pulse_duration + self.afterpulsetime
        self.numtimesteps  = int(self.interactiontime / self.dt)

    def make_pulse(self,time):

        end_pulse = self.pulse_duration
        field = np.zeros(2)
        if(time<=0.0 or time >= end_pulse):
            envelope = 0.0
            envelopederiv = 0.0
        else:
            envelope = np.sin(np.pi * time / end_pulse)**2
            envelopederiv = 2.0 * np.pi / end_pulse * \
            np.sin(np.pi * time / end_pulse) * \
            np.cos(np.pi * time / end_pulse)

        field[1] = envelope * self.A0 * np.cos(self.w0 * time + self.cep)
        field[0] = envelope * self.E0 * \
                        np.sin(self.w0 * time + self.cep) - \
                        envelopederiv * self.A0 * np.cos(self.w0 * time + self.cep)
        return field
