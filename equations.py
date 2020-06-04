#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Gray Martin
"""

#%% Imports
import numpy as np
from constants import *

#%% Solar
def solar_calc(lat_mask_RTMA, lon_mask_RTMA, jday_RTMA_mask, hour_RTMA,
               lat_mask_NDFD2, lon_mask_NDFD2, jday_NDFD2_mask, hour_NDFD2):    
    lat_rad_RTMA = (abs(lat_mask_RTMA))*(np.pi/180)
    lon_rad_RTMA =(abs(lon_mask_RTMA-360))*(np.pi/180)
    lat_rad_NDFD2 = (abs(lat_mask_NDFD2))*(np.pi/180)
    lon_rad_NDFD2 =(abs(lon_mask_NDFD2-360))*(np.pi/180)
    
    dec_RTMA = toc_latitude_r*np.cos((2*np.pi*(jday_RTMA_mask-summer_solstice))/days_per_year)    
    dec_NDFD2 = toc_latitude_r*np.cos((2*np.pi*(jday_NDFD2_mask-summer_solstice))/days_per_year)
    
    hra_RTMA = np.pi*hour_RTMA/12
    hra_NDFD2 = np.pi*hour_NDFD2/12
    
    elev_RTMA = np.arcsin((np.sin(dec_RTMA)*np.sin(lat_rad_RTMA)) 
                          - (np.cos(dec_RTMA)*np.cos(lat_rad_RTMA)*np.cos(hra_RTMA-lon_rad_RTMA)))
    elev_NDFD2 = np.arcsin((np.sin(dec_NDFD2)*np.sin(lat_rad_NDFD2)) 
                           - (np.cos(dec_NDFD2)*np.cos(lat_rad_NDFD2)*np.cos(hra_NDFD2-lon_rad_NDFD2)))
    
    zenith_RTMA = 90 - (elev_RTMA*(180/np.pi))
    zenith_NDFD2 = 90 - (elev_NDFD2*(180/np.pi))

    zenith_NDFD = np.concatenate((zenith_NDFD2[:,:,0:36], 
                                  zenith_NDFD2[:,:,38::3]),
                                 axis=2)
    return zenith_RTMA, zenith_NDFD, zenith_NDFD2
    
def rh_calc(temp, dewp):
    ea = (((7.5*dewp)/(237.3+dewp)))     
    es = (((7.5*temp)/(237.3+temp)))
    rat = 10**(ea-es)
    rh = (rat)*100.0
    return rh

def solar_rad(jday, hour, lat, lon, zenith, elev):    
    el = 90-zenith
    R = (1.0 + (0.0224*np.cos(2.0*np.pi*jday/365.0))) 
    sinal = np.sin(el*(np.pi/180))
    rm = ((((288.0 - 0.0065 * elev) / 288.0)**5.256) / (sinal + 0.15 * ((el + 3.885)**(-1.253))))
    toa = solar_irradiance * sinal / (R * R)
    solar_rad_max = np.where(sinal > 0, toa*(atc**rm), 0)
    return solar_rad_max

def solar_max(jday_mask, zenith):
    a = 149600000
    e = 0.017  
    
    theta = (jday_mask*360)/365.25
    theta_rad = theta*(np.pi/180)
    R = a*(1-e*e)/(1+e*np.cos(theta_rad))
    d = R/a
    Smax = np.where(zenith <= 89.5, 1367*np.cos(zenith*(np.pi/180))/(d*d), 
                    0.000000000001)
    return Smax

def direct_diffuse(zenith, sstar):
    fdb = np.where(zenith <= 89.5, np.exp(3 - 1.34*sstar - 1.65/sstar), 0)
    fdif = np.where(zenith <= 89.5, 1 - fdb, 0)
    
    return fdb, fdif

#%% Wind
def stability(night, speed, solar, dt=-0.1):
    stabt_day = np.zeros_like(solar)
    
    srad_1 = np.where(solar >= 925, True, False)
    srad_2 = np.where((solar < 925) & (solar >= 625), True, False)
    srad_3 = np.where((solar < 625) & (solar >= 175), True, False)
    srad_4 = np.where(solar <175, True, False)
    
    speed_2 = np.where(speed < 2, True, False)
    speed_3 = np.where((speed >= 2) & (speed < 3), True, False)
    speed_5 = np.where((speed >= 3) & (speed < 5), True, False)
    speed_6 = np.where((speed >= 5) & (speed < 6), True, False)
    speed_m = np.where(speed >= 6, True, False)
    
    stabt_day[(srad_1 & speed_2) | (srad_1 & speed_3) | (srad_2 & speed_2)] = 1
    stabt_day[(srad_1 & speed_5) | (srad_2 & speed_3) | (srad_2 & speed_5) | (srad_3 & speed_2)] = 2
    stabt_day[(srad_1 & speed_6) | (srad_1 & speed_m) | (srad_2 & speed_6) | (srad_3 & speed_3) | (srad_3 & speed_5)] = 3
    stabt_day[(srad_2 & speed_m) | (srad_3 & speed_6) | (srad_3 & speed_m) | (srad_4)] = 4

    stabt_day = stabt_day*night
    
    stabt_night = np.zeros_like(solar)
    
    speed_2 = np.where((speed < 2), True, False)
    speed_25 = np.where((speed >= 2) & (speed < 2.5), True, False)
    speed_m = np.where((speed >= 2.5), True, False)
    
    stabt_night[speed_2] = 5
    stabt_night[speed_25] = 4
    stabt_night[speed_m] = 4
    
    stabt_night = stabt_night*np.logical_not(night)
    stabt = stabt_day + stabt_night
    
    return stabt

def est_wind_speed(speed, stability, zspeed=10, urban=False):
    ref_height = 2
    if urban:
        expo = [0.15, 0.15, 0.2, 0.25, 0.3, 0.3]
    else:
        expo = [0.07, 0.07, 0.1, 0.15, 0.35, 0.55]
    
    expon = np.zeros_like(stability)
    
    for i in range(1, 7):
        expon[stability == i] = expo[i-1]
    
    est_speed = speed*(ref_height/zspeed)**expon
    est_speed = np.where(est_speed < 0.13, 0.13, est_speed) 
    
    return est_speed

#%% Natural Wet Bulb Temperature
def rh_fraction(rh_rtma, rh_mix):
    rh_rtma = rh_rtma/100.0
    rh_mix = rh_mix/100.0

def viscosity(temp):
    # The values of these constants weren't included in the original script
    sigma = 3.617
    eps_kappa = 97
    
    Tr = temp / eps_kappa
    omega = (Tr-2.9) / 0.4 * (-0.034) + 1.048
    return (2.6693*10**(-6) * np.sqrt(mw_dry_air*temp) / (sigma * sigma * omega))

def thermal_conductivity(temp):
    return ((c_p + 1.25 * R_dry_air) * viscosity(temp)) 

def h_cylinder_in_air(temp, speed, P=1000):
    # The values of these constants weren't included in the original script
    a = 0.56
    b = 0.281
    c = 0.4
    
    density = P*100 / (R_dry_air*temp)
    Re = speed * density * wick_diameter / viscosity(temp)      
    Nu = b*Re**(1-c)*P_dry_air**(1-a)                                   
    return(Nu * (thermal_conductivity(temp)/wick_diameter))

def esat(temp, phase=0, P=1000):    
    y = (temp - 273.15)/(temp - 32.18)
    es = 6.1121 * np.exp( 17.502 * y )
    
    #if (P > 800):
    #  es = 1.004 * es
    #else:
    #  es = 1.0034 * es
    return(es*1.004)

def diffusivity(temp, P=1000):
    Pcrit_air = 36.4
    Pcrit_h2o = 218
    Tcrit_air = 132
    Tcrit_h2o = 647.3
    a = 3.640*10**(-4)
    b = 2.2244
    
    Pcrit13  = (Pcrit_air * Pcrit_h2o)**(1/3)
    Tcrit512 = (Tcrit_air * Tcrit_h2o)**(5/12)
    Tcrit12  = np.sqrt(Tcrit_air * Tcrit_h2o)
    Mmix = np.sqrt((1/mw_dry_air) + (1/mw_water_vapor))
    Patm = P / 1013.2
    
    return(a * ((temp/Tcrit12)**b) * Pcrit13 * Tcrit512 * Mmix / Patm * (1*10**(-4)))

def evaporation(temp):
    return((313.15-temp)/30 * -71100 + 2.4073*(10**6))

def emis_atm(temp, RH):
    e = RH*esat(temp, 0)
    return (0.575 * (e**0.143))

def twb(temp, dewp, RH, wind_speed, srad, fdb, cza, rad, P=1000, a=0.56):
    Tsfc = temp
    sza = np.arccos(cza)
    eair = RH * esat(temp)
    twb_prev = dewp + 273.15
    converged = False
    max_iter = 25
    convergence = 0.02
    i = 0
    while i <= max_iter:
        i = i + 1

        Tref = 0.5*(twb_prev + temp)
        h = h_cylinder_in_air(temp, wind_speed)
        Fatm = stefan_boltzmann * wick_emissivity * (0.5*( emis_atm(temp, RH)*temp**4 + surface_emissivity*(Tsfc**4) - twb_prev**4)) + (1-wick_albedo) * srad * ((1-fdb)*(1+0.25*wick_diameter/wick_length) + fdb*((np.tan(sza)/np.pi)+0.25*wick_diameter/wick_length) + surface_albedo)
        ewick = esat(twb_prev)
        density = P * 100 / (R_dry_air * Tref)
        Sc = viscosity(Tref)/(density*diffusivity(Tref,P))
        Twb_new = temp - evaporation(Tref)/c_p*mw_dry_air/mw_water_vapor * (ewick-eair)/(P-ewick) * (P_dry_air/Sc)**a + (Fatm/h*rad)

        diff = np.abs(Twb_new-twb_prev)
        twb_prev = 0.9*twb_prev + 0.1*Twb_new
        if (diff < convergence).all():
            converged = True
            break
    return Twb_new-273.15

def h_sphere_in_air(temp, wind_speed, P=1000):
    density = P * 100 / (R_dry_air*temp)
    Re = wind_speed * density * globe_diameter / viscosity(temp)     
    Nu = 2 + 0.6*np.sqrt(Re) * (P_dry_air**0.3333)
    return(Nu * thermal_conductivity(temp) / globe_diameter)

def tglobe(temp, dewp, RH, wind_speed, srad, fdb, cza, P=1000, a=0.56):
    Tsfc = temp
    tglobe_prev = dewp + 273.15
    converged = False
    max_iter = 10
    convergence = 0.02
    i = 0
    while i <= max_iter:
        i = i + 1

        h = h_sphere_in_air(temp, wind_speed)
        tglobe_new = (0.5*(emis_atm(temp, RH)*temp**4 + surface_emissivity*(Tsfc**4)) - (h/(stefan_boltzmann*globe_emissivity)*(tglobe_prev-temp)) + (srad/(2*stefan_boltzmann*globe_emissivity)*(1-globe_albedo)*(fdb*(1/(2*cza)-1)+1 + surface_albedo)))**0.25

        diff = np.abs(tglobe_new-tglobe_prev)
        tglobe_prev = 0.9*tglobe_prev + 0.1*tglobe_new
        if (diff < convergence).all():
            converged = True
            break
    return tglobe_new-273.15