"""Carmela code's example

   Let's run the one-dimensional
   Schroedinger equation.
"""
#!/usr/bin/env python

import numpy as np
from carmela_code import constants as const
from carmela_code import input_reader as inp
from carmela_code import axes as axes
from carmela_code import wavefunction as wavef
from carmela_code import hamiltonian as ham
from carmela_code import momentum as mom
from carmela_code import laser
from carmela_code import absorber
from carmela_code import flux


print('Starting!\n')

# Fetch parameters form input file
inp_params = {}
inp.read_input(inp_params,'./input.inp')
params = inp.parameters(inp_params)

# Create mesh
print('Creating mesh')
ax = axes.axes(params.maxptsx,params.dx)

# Create wavefunction
print('Creating wavefunction')
psi = wavef.wavefunction(ax)

# Create hamiltonian
print('Creating hamiltonian')
hydro_ham = ham.hamiltonian(ax,'hydrogen',1,
                    gauge_type='velocity',omega=2.0)

# Create absorber
if(params.fourier_on==False):
    print('Creating absorber')
    absor = absorber.absorber(ax,params.xsplit,params.Medge)
    
# Diagonalize hamiltonian
print('Getting ground state of the system')
psi.psi = hydro_ham.get_ground_state(5)
## Normalize wavefunction
psi.normalize()

laser = laser.laser(params.dt,params.intensity,params.wavelength,
                    params.no_cycles,params.afterpulse)

# Prepare arrays for time evolution
poptotal = np.zeros(laser.numtimesteps)
popion   = np.zeros(laser.numtimesteps)
acc      = np.zeros(laser.numtimesteps)
time     = np.zeros(laser.numtimesteps)
field    = np.zeros((laser.numtimesteps,2))

sampling_indx = np.argmin(abs(ax.ax-params.surf_rad))
xb = ax.ax[sampling_indx]
if(params.sampling_point_on):
    detector   = np.zeros(laser.numtimesteps,dtype=complex)
if(params.tsurff_on):
    surf       = np.zeros(laser.numtimesteps,dtype=complex)
    deriv_surf = np.zeros(laser.numtimesteps,dtype=complex)

#####################
# Time-propagate!!
#####################
print('------------------------------')
print('Lets begin the time evolution!')
print('------------------------------')

### Set first the observables for the initial state.
poptotal[0]= psi.get_norm()
popion[0]  = psi.get_population(xb)
acc[0]     = psi.get_acceleration(hydro_ham.potential_deriv,0.0)

print('\n')
print('Time:                        ',str(time[0]*const.autime_fs))
print('Total population in the box: ',str(poptotal[0]))
print('Ionised population:          ',str(popion[0]))
print('Acceleration:                ',str(acc[0]))

for itime in range(1,laser.numtimesteps):
    time[itime] = float(itime) * laser.dt
    field[itime,:] = laser.make_pulse(time[itime])

    psi.psi = hydro_ham.time_propagate(psi.psi,params.dt,6,field[itime,:])
    
    if(params.fourier_on==False):
        psi.psi = absor.make_split(psi.psi)
        
    # Get surface values for momentum calculation
    if(params.tsurff_on):
        surf[itime], deriv_surf[itime] = flux.get_surface(psi.psi,ax,xb)
    if(params.sampling_point_on):
        detector[itime] = psi.psi[sampling_indx]
    poptotal[itime]= psi.get_norm()
    popion[itime]  = psi.get_population(xb)
    acc[itime]     = psi.get_acceleration(hydro_ham.potential_deriv, field[itime,0])

    print('\n')
    print('Time:                        ',str(time[itime]*const.autime_fs))
    print('Total population in the box: ',str(poptotal[itime]))
    print('Ionised population:          ',str(popion[itime]))
    print('Acceleration:                ',str(acc[itime]))



np.savetxt('field.dat',np.transpose((time,field[:,0],field[:,1])))
np.savetxt('populations.dat',np.transpose((time,poptotal,popion)))
np.savetxt('acceleration.dat',np.transpose((time,acc)))
np.savetxt('prob.dat',np.transpose((ax.ax,np.abs(psi.psi)**2)))


# Get momentum amplitude
if(params.fourier_on):
    mask_radius = np.array([params.inner_mask,params.outer_mask])
    kax, psik = mom.get_momentum_amplitude(ax,psi,mask_radius)
    fourier_energy = kax.ax**2*0.5
    fourier_pes = np.abs(psik.psi/kax.ax)**2
    np.savetxt('fourier_pes.dat',np.transpose((fourier_energy,fourier_pes)))

# Get pes from sampling point method
if(params.sampling_point_on):
    energy, pes = flux.get_sampling_pes(time,detector)
    np.savetxt('sp_pes.dat',np.transpose((energy,pes)))

if(params.tsurff_on):
    bk  = flux.get_tsurff(surf,deriv_surf,ax,xb,time,field)
    kax = axes.kaxes(ax)
    ts_ener = kax.ax**2*0.5
    tspes = np.abs(bk/kax.ax)**2
    tspes[np.where(kax.ax==0.0)] = 0.0
    np.savetxt('tsurff_pes.dat',np.transpose((ts_ener,tspes)))
    
if(params.draw_field_wave):
    timefs = time * const.autime_fs
    fig  = plt.figure(figsize=(14,10),facecolor='white')
    ax1 = plt.subplot()
    ax1.plot(timefs,field[:,1],linewidth=2)
    ax2 = ax1.twinx()
    ax2.plot(timefs,abs(surf)**2,'r',linewidth=2)
    ax1.axis('tight')
    ax1.set_xticks(np.arange(0,max(timefs),5.0))
    #ax1.yticks(np.arange(-20.,0.,1.))
    #ax1.xlim(0,4)
    #ax1.ylim(-9,0)
    ax1.tick_params(labelsize=20)
    ax2.tick_params(labelsize=20)
    ax1.set_xlabel('Time (fs)',fontsize=30)
    ax1.set_ylabel('E field (au)',fontsize=30)
    ax2.set_ylabel(r'$|\Psi(t)|^2$ (arb. units)',fontsize=30)
    #plt.title('Efield_wavefunction',fontsize=20)
    
    fig.tight_layout()
    
    if(params.makeframe):
        fig.savefig('field_wavefunction.pdf',format='pdf')


print('\nFin!')
