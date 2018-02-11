#!/usr/bin/env python

import numpy as np
from . import axes as ax
from . import fdrule as fd

class hamiltonian():
    def __init__(self,axis,potential_type,fdrule=2,
                 gauge_type='velocity',omega=1.0):

        self.gauge_type = gauge_type
        if(gauge_type != 'velocity' and gauge_type != 'length'):
            print('We have no such gauge for you.')
        self.omega = omega
        self.axis = axis
        self.npts = axis.ax.shape[0]
        self.potential = np.zeros(self.npts)
        self.create_potential(potential_type)
        self.fdrule = fdrule
        self.rulepts = 2 * fdrule + 1
        maxfd = fdrule * self.axis.dx
        fdax = np.linspace(-maxfd,maxfd,self.rulepts)
        self.fdcoeff = fd.calculate_fdweights(0.0,fdax,2)

    def create_potential(self,potential_type):
        if(potential_type == 'free'):
            self.potential = np.zeros(self.npts)
            self.potnetial_deriv = np.zeros(self.npts)
        elif(potential_type == 'hydrogen'):
            self.potential = -1.0/ np.sqrt(self.axis.ax**2 + self.omega)
            exponent = 3.0/2.0
            self.potential_deriv = self.axis.ax / np.power(self.axis.ax**2
                                                           + self.omega, exponent)
        elif(potential_type == 'harmonic'):
            self.potential = 0.5 * self.omega * self.axis.ax**2
            self.potential_deriv = self.omega * self.axis.ax
        else:
            print('There is no such potential!')
            
    def apply_hamiltonian(self,psi,field):
        hampsi = self.apply_laplacian(psi)
        hampsi += self.potential * psi
        if(self.gauge_type == 'velocity'):
            auxpsi = np.zeros(len(psi)+2*self.fdrule,dtype=psi.dtype)
            auxpsi[self.fdrule:-self.fdrule] = psi
            for ik in range(len(psi)):
                ii = ik + self.rulepts
                hampsi[ik] -= 1.j* field[1] * \
                sum(self.fdcoeff[:,1] * auxpsi[ik:ii])
        elif(self.gauge_type == 'length'):
            hampsi += field[0] * self.axis.ax * psi
        return hampsi
         
    def apply_laplacian(self,psi):
        lappsi = np.zeros(len(psi),dtype=psi.dtype)
        auxpsi = np.zeros(len(psi)+2*self.fdrule,dtype=psi.dtype)
        auxpsi[self.fdrule:-self.fdrule] = psi
        
        for ik in range(len(psi)):
            ii = ik + self.rulepts
            lappsi[ik] = sum(self.fdcoeff[:,2] * auxpsi[ik:ii])
        lappsi *= -0.5
        return lappsi

    def diagonalize_hamiltonian(self):
        ham = np.zeros((self.npts,self.npts),dtype=complex)
        for j in range(self.npts):
            unit_vec = np.zeros(self.npts,dtype=complex)
            unit_vec[j] = 1.0
            ham[:,j] = self.apply_hamiltonian(unit_vec,
                                              ([0.0,0.0]))
        eigval,eigvec = np.linalg.eigh(ham)
  
        return eigval,eigvec

    def get_ground_state(self,noenergies=1):
        eigval, eigvec = self.diagonalize_hamiltonian()
        print('Eigenenergies:        \n')
        print('------------------------')
        for i in range(noenergies):
            print('Energy no'+str(i)+': '+str(eigval[i]))
        print('------------------------\n')
        return eigvec[:,0]


    def time_propagate(self,psi,dt,order,field):
        hamdtpsi = np.zeros(len(psi),dtype='complex128')
        dtpsi    = np.zeros(len(psi),dtype='complex128')
        dtpsi[:] = psi
        taylorfactor = 1.0
        # For each order of the taylor series apply the hamiltonian.
        for itaylor in range(1,order):
            hamdtpsi = self.apply_hamiltonian(dtpsi,field)
            # Increment taylor series. Add next higher order derivative 
            # of psi (hamdtpsi) to psi.
            taylorfactor *= -1j * dt / float(itaylor)
            psi += taylorfactor * hamdtpsi
            # Use the vector from applying the power order of the Hamiltonian
            dtpsi[:] = hamdtpsi
        return psi
