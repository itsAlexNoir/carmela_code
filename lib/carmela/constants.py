"""
This module defines a set of physical constants commonly used in atomic units (au) 
and other unit systems. These constants are useful for scientific computations 
involving atomic and molecular physics.

Constants:
----------
- M1, M2: Masses of particles in atomic mass units (amu).
- femto: Conversion factor for femtoseconds (1e-15 seconds).
- autime_s: Atomic unit of time in seconds.
- autime_fs: Atomic unit of time in femtoseconds.
- aulength_m: Atomic unit of length in meters.
- aulength_nm: Atomic unit of length in nanometers.
- intensity_au: Intensity in atomic units (W/cm²).
- energy_au_ev: Energy in atomic units (eV).
- speed_light_au: Speed of light in atomic units.
- fine_struc_au: Fine structure constant in atomic units.
- waveconvert: Conversion factor to convert wavelengths from nanometers (nm) to atomic units (au).
"""
########################################
############## Constants ###############
########################################

M1 = 1836.0
M2 = 1836.0
femto = 1e-15
autime_s = 2.418884326505e-17 # in s
autime_fs = 2.418884326505e-2 # in fs
aulength_m = 5.2917721092e-11 # m
aulength_nm = 5.2917721092e-2 # nm
intensity_au = 3.5e16 # W/cm2
energy_au_ev = 27.211 # in eV
speed_light_au = 137.0
fine_struc_au = 1.0 / speed_light_au
waveconvert = 45.7720588235 # To convert from nm to au
