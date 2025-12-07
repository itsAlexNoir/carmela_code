import sys as sys
import numpy as np
import yaml
from pydantic import BaseModel, Field

from . import constants as const

########################################
##           config.py
########################################
class GridConfig(BaseModel):
    dt: float = Field(default=0.02, description="Time step size")
    dx: float = Field(default=0.2, description="Spatial step size")
    maxptsx: int = Field(default=2048, description="Number of spatial grid points")
    maxptsk: int = Field(default=2048, description="Number of momentum grid points")
    dke: float = Field(default=0.01, description="Momentum step size")
    dE: float = Field(default=0.01, description="Energy step size")
    maxke: float = Field(default=10.0, description="Maximum momentum value")
    maxE: float = Field(default=10.0, description="Maximum energy value")

class LaserConfig(BaseModel):
    wavelength: float = Field(default=20.0, description="Laser wavelength in nm")
    intensity: float = Field(default=1.0e12, description="Laser intensity in W/cm^2")
    no_cycles: float = Field(default=3.0, description="Number of cycles in the laser pulse")
    after_pulse: int = Field(default=0, description="Flag for after pulse")
    Eg: float = Field(default=0.5, description="Ground state energy")
    Ip: float = Field(default=1.100000922215923, description="Ionization potential")
    Ee_ejec_eV: float = Field(default=10.0, description="Ejected electron energy in eV")

class AbsorberConfig(BaseModel):
    surf_radius: float = Field(default=20.0, description="Surface radius for absorber")
    xsplit: float = Field(default=20.0, description="Split point for absorber")
    Medge: float = Field(default=0.2, description="Minimum edge value for absorber")
    inner_mask: float = Field(default=15.0, description="Inner mask radius")
    outer_mask: float = Field(default=30.0, description="Outer mask radius")

class MomentumConfig(BaseModel):
    fourier_on: bool = Field(default=False, description="Enable Fourier transform method")
    sampling_point_on: bool = Field(default=False, description="Enable sampling point method")
    tsurff_on: bool = Field(default=False, description="Enable TSURFF method")

class OutputConfig(BaseModel):
    draw_field_wave: bool = Field(default=False, description="Flag to draw field and wavefunction")
    draw_x: bool = Field(default=False, description="Flag to draw spatial representation")
    draw_k: bool = Field(default=False, description="Flag to draw momentum representation")
    draw_ke: bool = Field(default=False, description="Flag to draw kinetic energy representation")
    draw_E: bool = Field(default=False, description="Flag to draw energy representation")
    draw_total_cross: bool = Field(default=False, description="Flag to draw total cross section")
    makemovie: bool = Field(default=False, description="Flag to make movie")
    makeframe: bool = Field(default=True, description="Flag to make frame")
    frametime: int = Field(default=0, description="Frame time index")
    drawover: bool = Field(default=False, description="Flag to draw overlay")
    calc_k: bool = Field(default=False, description="Flag to calculate momentum space")
    fftshift_on: bool = Field(default=False, description="Enable FFT shift")

class Config(BaseModel):
    grid: GridConfig = Field(default_factory=GridConfig)
    laser: LaserConfig = Field(default_factory=LaserConfig)
    absorber: AbsorberConfig = Field(default_factory=AbsorberConfig)
    momentum: MomentumConfig = Field(default_factory=MomentumConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)

    calc_k: bool = Field(default=False, description="Flag to calculate momentum space")

    @classmethod
    def from_yaml(cls, yaml_path: str):
        with open(yaml_path, "r") as file:
            data = yaml.safe_load(file)

        grid = GridConfig(**data.get("grid", {}))
        laser = LaserConfig(**data.get("laser", {}))
        absorber = AbsorberConfig(**data.get("absorber", {}))
        momentum = MomentumConfig(**data.get("momentum", {}))
        output = OutputConfig(**data.get("output", {}))
        calc_k = data.get("calc_k", False)

        return cls(
            grid=grid,
            laser=laser,
            absorber=absorber,
            momentum=momentum,
            output=output,
            calc_k=calc_k,
        )


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
