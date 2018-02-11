#!/usr/bin/env python

########################
##   wavefunction.py
########################
## A little class for
## wavefunctions
########################
import numpy as np

class wavefunction():
    def __init__(self,axis):

        self.axis = axis
        self.numpts = axis.maxptsx
        self.dx = axis.dx
        self.psi = np.zeros((self.numpts),dtype=complex)
        self.norm = 0.0
        
    def get_norm(self):
        norm = np.real(np.dot(np.conj(self.psi),self.psi))
        return norm * self.dx

    def normalize(self):
        norm = self.get_norm()
        self.psi /= np.sqrt(norm)

    def get_population(self,ion_radius):
        pop = np.real(np.dot(np.conj(
            self.psi[np.where(abs(self.axis.ax)>ion_radius)]),
                             self.psi[np.where(abs(self.axis.ax)>ion_radius)]))
        return pop

    def get_acceleration(self, pot_deriv, efield):
        acc = np.real(np.dot(np.conj(self.psi), (pot_deriv - efield)*self.psi))
        return acc
