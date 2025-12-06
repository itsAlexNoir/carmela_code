#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import fileinput as fileinp
import sys as sys
from . import constants as const

########################################
##           input_reader.py
########################################

def read_input(parameters,filepath):
    #for line in fileinp.input(sys.argv[1:]):
    for line in fileinp.input(filepath):
        if(line[0]!='#' and line[0]!='\n'):
            key, value = line.split("=")
            parameters[key.strip()] = value.strip()
            
    ########################################
    
    
class parameters:
    def __init__(self,inparams):
        self.dt = float(inparams['dt'])
        self.dx = float(inparams['dx'])
        self.maxptsx  = int(inparams['maxptsx'])
        self.maxptsk = int(inparams['maxptsk'])
            
        self.dke = float(inparams['dke'])
        self.dE = float(inparams['dE'])
        
        self.maxke = float(inparams['maxke'])
        self.maxE = float(inparams['maxE'])
        self.maxptske = int(self.maxke / self.dke)
        self.maxptsE = int(self.maxE / self.dE)
        
         # Let's convert the frequency to atomic units
        self.intensity = float(inparams['intensity'])         
        self.wavelength = float(inparams['wavelength'])
        self.afterpulse = int(inparams['after_pulse'])
        self.wavelength0 = self.wavelength / const.aulength_nm
        self.w0 = 2.0 * np.pi * const.speed_light_au / self.wavelength0
        self.period = 2.0 * np.pi / self.w0
        self.no_cycles = float(inparams['no_cycles'])
        self.pulse_duration = self.period * self.no_cycles
        if(self.pulse_duration!=0):
            self.pulse_bandwidth = 4.0 * np.pi / self.pulse_duration
        else:
            self.pulse_bandwidth = 0.0
        self.quiver = np.sqrt(self.intensity / const.intensity_au) / self.w0 / self.w0
        self.Up = (self.intensity / const.intensity_au) / 4. / self.w0 / self.w0
        self.E0 = np.sqrt(self.intensity / const.intensity_au)
        self.A0 = self.E0 / self.w0
        self.Ee_ejec_au = float(inparams['Ee_ejec_eV']) / const.energy_au_ev # in au
        self.ke_ejec = np.sqrt(2.0 * self.Ee_ejec_au)
        self.Eg = float(inparams['Eg'])
        self.Ip = self.Eg - 0.5

        self.surf_rad = float(inparams['surf_radius'])
        self.xsplit = float(inparams['xsplit'])
        self.Medge  = float(inparams['Medge'])

        self.makemovie = int(inparams['makemovie'])
        self.makeframe = int(inparams['makeframe'])
        self.frametime = int(inparams['frametime'])
        self.drawover = int(inparams['drawover'])

        self.calc_k = int(inparams['calc_k'])
        self.fftshift_on = int(inparams['fftshift_on'])
        self.inner_mask = float(inparams['inner_mask'])
        self.outer_mask = float(inparams['outer_mask'])
        if(int(inparams['fourier_on'])==1):
           self.fourier_on = True
        else:
            self.fourier_on = False
            
        if(int(inparams['tsurff_on'])==1):
           self.tsurff_on = True
        else:
            self.tsurff_on = False
            
        if(int(inparams['sampling_point_on'])==1):
            self.sampling_point_on = True
        else:
            self.sampling_point_on = False

        self.draw_field_wave = int(inparams['draw_field_wave'])
        self.draw_x  = int(inparams['draw_x'])
        self.draw_k  = int(inparams['draw_k'])
        self.draw_ke = int(inparams['draw_ke'])
        self.draw_E  = int(inparams['draw_E'])
        self.draw_total_cross = int(inparams['draw_total_cross'])
