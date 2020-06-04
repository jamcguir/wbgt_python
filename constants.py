#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 16:06:04 2020

@author: gray
"""

#%% Imports
import numpy as np

#%% Physics
# Constant:     Stefan-Boltzmann constant
# Units:        W/(m^2) * K^(-4)
stefan_boltzmann = 5.6704*(10**(-8))

#%% Solar
# Constant: 
# Units: 

# Constant:     Latitude of Tropic of Cancer
# Units:        Degrees 
toc_latitude_d = 23.45

# Constant:     Latitude of Tropic of Cancer
# Units:        Radians 
toc_latitude_r = toc_latitude_d*(np.pi/180)

# Constant:     Solar constant (average solar irradiance)
# Units:        W/(m^2) 
solar_irradiance = 1370

# Constant:     Atmospheric transfer coefficient
# Units:        Unitless
atc = 0.75

#%% Time
# Constant:     Average days per year
# Units:        Unitless
days_per_year = 365.25

# Constant:     Day of summer solstice 
# Units:        Unitless
summer_solstice = 173

#%% Gas Chemistry
# Constant:     Molecular weight of water vapor
# Units:        kg/kmol
mw_water_vapor = 18.015

# Constant:     Molecular weight of dry air
# Units:        kg/kmol
mw_dry_air = 28.97

# Constant:     Specific heat capacity of dry air  (0 C, 1 atm)
# Units:        J/(kg*K)
c_p = 1003.5

# Constant:     Gas constant (universal)
# Units:        J/(K*kmol)
R = 8314.34

# Constant:     Gas constant (dry air)
# Units:        J/(K*kmol)
R_dry_air = R/mw_dry_air

# Constant:     Prandtl number (dry air)
# Units:        Unitless
P_dry_air = c_p/(c_p + 1.25*R_dry_air)

#%% Wick Properties
# Constant:     Wick emissivity 
# Units:        Unitless
wick_emissivity = 0.95

# Constant:     Wick albedo 
# Units:        Unitless
wick_albedo = 0.4

# Constant:     Wick length 
# Units:        m
wick_length = 0.0254 

# Constant:     Wick diameter 
# Units:        m
wick_diameter = 0.007

#%% Globe Properties
# Constant:     Globe emissivity 
# Units:        Unitless
globe_emissivity = 0.95

# Constant:     Globe albedo 
# Units:        Unitless
globe_albedo = 0.05

# Constant:     Globe diameter 
# Units:        m
globe_diameter = 0.0508

#%% Surface Properties
# Constant:     Surface emissivity 
# Units:        Unitless
surface_emissivity = 0.999

# Constant:     Surface albedo 
# Units:        Unitless
surface_albedo = 0.45
