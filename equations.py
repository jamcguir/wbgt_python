#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Gray Martin
"""

#%% Imports
import numpy as np

#%% Solar
def solar_calc(lat_mask_RTMA, lon_mask_RTMA, jday_RTMA_mask, hour_RTMA,
               lat_mask_NDFD2, lon_mask_NDFD2, jday_NDFD2_mask, hour_NDFD2):
    phi_r = 23.45*(np.pi/180)
    solst = 173
    sconst = 1370
    num = 365.25
    
    lat_rad_RTMA = (abs(lat_mask_RTMA))*(np.pi/180)
    lon_rad_RTMA =(abs(lon_mask_RTMA-360))*(np.pi/180)
    lat_rad_NDFD2 = (abs(lat_mask_NDFD2))*(np.pi/180)
    lon_rad_NDFD2 =(abs(lon_mask_NDFD2-360))*(np.pi/180)
    
    dec_RTMA = phi_r*np.cos((2*np.pi*(jday_RTMA_mask-solst))/num)    
    dec_NDFD2 = phi_r*np.cos((2*np.pi*(jday_NDFD2_mask-solst))/num)
    
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
    ea = (6.11*10**((7.5*dewp)/(237.3+dewp)))     
    es = (6.11*10**((7.5*temp)/(237.3+temp)))
    rh = (ea/es)*100
    return rh

def solar_rad(jday, hour, lat, lon, zenith, elev):    
    el = 90-zenith  # solar elevation degrees from horizon
    
    R = (1.0 + (0.0224*np.cos(2.0*np.pi*jday/365.0))) # earth sun distance
    nrel = 1367.0 #solar constant
    
    sinal = np.sin(el*(np.pi/180))
    
    rm = ((((288.0 - 0.0065 * elev) / 288.0)**5.256) / (sinal + 0.15 * ((el + 3.885)**(-1.253))))
    toa = nrel * sinal / (R * R)
    
    atc = 0.75 # atmospheric transmission coefficient (0.7-0.91); Jordan found 0.75 gave better Tg and Twb
    
    solar_rad_max = np.where(sinal > 0, toa*(atc**rm), 0)
    
    return solar_rad_max

def solar_max(jday_mask, zenith):
    sig = 5.6704*10**(-8)
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