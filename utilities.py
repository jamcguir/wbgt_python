#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Gray Martin
"""

#%% Imports
import geopandas as gpd
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

import time
import datetime 
import logging
log = logging.getLogger(__name__)

from netCDF4 import Dataset

#%% Functions
#%%% Utility
def to_state(gdf, statelist=["North Carolina"]):
    """Returns an array of two-letter state abbreviations given arrays of 
    latitude and longitude."""
    
    log.debug("Processing {}".format(statelist))
    usa = gpd.read_file("resources/states/states.shp")
    states = usa.loc[usa["STATE_NAME"].isin(statelist)]
    points = gdf["geometry"]
    statelabels = [None]*len(states)
    statemasks = [None]*len(states)
    
    for row in range(0, len(states)):
        state = states.iloc[row]
        state_polygon = state["geometry"]
        statemasks[row] = points.within(state_polygon)
        statelabels[row] = state["STATE_ABBR"]
        log.debug("Finished processing {}".format(statelabels[row]))
    return statemasks, statelabels, states

def timing(z=6):
    log.debug("Processing timing for z={}".format(z))
    today_dt = datetime.datetime.utcnow()
    today_d = today_dt.date()
    
    yesterday_d = today_dt.date() - datetime.timedelta(days=1)
    fivedays_d = today_dt.date() + datetime.timedelta(days=5)
    
    z_duration = (datetime.datetime.min + datetime.timedelta(hours=z)).time()
    z1_duration = (datetime.datetime.min + datetime.timedelta(hours=z+1)).time()
    step = datetime.timedelta(hours = 1)
    
    # RTMA
    RTMA_end = datetime.datetime.combine(today_d, z_duration)
    RTMA_start = datetime.datetime.combine(yesterday_d, z_duration)
    
    RTMA_dates = [None]*25
    RTMA_dates_int = [None]*25
    RTMA_hours = [None]*25
    n = 0
    while RTMA_start + n*step <= RTMA_end:
        RTMA_date = RTMA_start + n*step
        RTMA_dates[n] = RTMA_date
        RTMA_dates_int[n] = RTMA_date.timetuple().tm_yday
        RTMA_hours[n] = RTMA_date.hour
        n += 1
    RTMA_dates = RTMA_dates[:n]
    RTMA_dates_int = RTMA_dates_int[:n]
    jday_RTMA = np.broadcast_to(RTMA_dates_int, (420, 370, 25))
    hour_RTMA = np.broadcast_to(RTMA_hours, (420, 370, 25))
    
    # NDFD 
    NDFD2_end = datetime.datetime.combine(fivedays_d, z_duration)
    NDFD2_start = datetime.datetime.combine(today_d, z1_duration)
    
    NDFD2_dates = [None]*120
    NDFD2_dates_int = [None]*120
    NDFD2_hours = [None]*120
    n = 0
    while NDFD2_start + n*step <= NDFD2_end:
        NDFD2_date = NDFD2_start + n*step
        NDFD2_dates[n] = NDFD2_date
        NDFD2_dates_int[n] = NDFD2_date.timetuple().tm_yday
        NDFD2_hours[n] = NDFD2_date.hour
        n += 1
    NDFD2_dates = NDFD2_dates[:n]
    NDFD2_dates_int = NDFD2_dates_int[:n]
    jday_NDFD2 = np.broadcast_to(NDFD2_dates_int, (420, 370, 120))
    hour_NDFD2 = np.broadcast_to(NDFD2_hours, (420, 370, 120))
    
    jday_NDFD = np.concatenate((jday_NDFD2[:,:,0:36], 
                                jday_NDFD2[:,:,38::3]),
                               axis=2)
    
    hour_NDFD = np.concatenate((hour_NDFD2[:,:,0:36], 
                                hour_NDFD2[:,:,38::3]),
                               axis=2)
    return jday_RTMA, hour_RTMA, jday_NDFD, hour_NDFD, jday_NDFD2, hour_NDFD2

#%%% Data Imports
def RTMA_import(filename):
    rootgrp = Dataset(filename, "a", format="NETCDF4")
    vs = [*rootgrp.variables.keys()]
    vars_d = dict.fromkeys(vs)
    data_d = dict.fromkeys(vs)
    unit_d = dict.fromkeys(vs)
    fill_d = dict.fromkeys(vs)
    for v in vs:
        vars_d[v] = rootgrp[v]
        data_d[v] = np.squeeze(rootgrp[v][:])
        data_d[v] = np.swapaxes(data_d[v], 0, 2)
        data_d[v].fill_value = np.nan
        unit_d[v] = rootgrp[v].getncattr("units")
        fill_d[v] = rootgrp[v].getncattr("_FillValue")
    return vars_d, data_d, unit_d, fill_d 

def NDFD2_import(filename):
    rootgrp = Dataset(filename, "a", format="NETCDF4")
    vs = [*rootgrp.variables.keys()]
    vars_d = dict.fromkeys(vs)
    data_d = dict.fromkeys(vs)
    unit_d = dict.fromkeys(vs)
    fill_d = dict.fromkeys(vs)
    for v in vs:
        vars_d[v] = rootgrp[v]
        data_d[v] = np.squeeze(rootgrp[v][:])
        data_d[v] = np.swapaxes(data_d[v], 0, 2)
        data_d[v].fill_value = np.nan
        data_d[v] = data_d[v][:,:,:46]
        unit_d[v] = rootgrp[v].getncattr("units")
        fill_d[v] = rootgrp[v].getncattr("_FillValue")
    return vars_d, data_d, unit_d, fill_d 

def small_import(filename):
    rootgrp = Dataset(filename, "a", format="NETCDF4")
    vs = [*rootgrp.variables.keys()]
    vars_d = dict.fromkeys(vs)
    data_d = dict.fromkeys(vs)
    for v in vs:
        vars_d[v] = rootgrp[v]
        data_d[v] = rootgrp[v][:]
    return vars_d, data_d

#%%% Bias Functions
def RTMA_bias(data_RTMA, z):
    if z == 12:
        bias_table = {
            "hour": np.array(list(range(12,23)) + list(range(0,12))),
            "temp_bias": np.array([1.86,1.74,1.26,0.74,0.28,-0.15,-0.5,-0.71,-0.81,-0.91,-1.11,-1.51,-1.99,-1.44,-1.07,-0.9,-0.79,-0.7,-0.67,-0.6,-0.57,-0.49,-0.2,1.14,1.86]) / 1.8,   #nums displayed are temp differences in F -> change to temp diff in K
            "dew_bias": np.array([-1.3,-1.85,-2.18,-2.35,-2.4,-2.4,-2.3,-2.19,-2.13,-1.96,-1.79,-1.56,-1.78,-1.76,-1.77,-1.72,-1.68,-1.69,-1.71,-1.68,-1.68,-1.67,-1.5,-0.93,-1.3]) / 1.8,
            "wind_bias": np.array([0.26,0.27,0.27,0.23,0.22,0.25,0.21,0.17,0.09,0.07,-0.01,-0.15,-0.02,0.09,0.1,0.11,0.1,0.12,0.12,0.14,0.14,0.15,0.2,0.26,0.26]),
            "srad_bias": np.array([82.56,35.49,-24.2,-87.89,-160.95,-229.6,-272.25,-287.4,-271.61,-221.73,-154.68,-68.81,-1.86,0.08,0,0,0,0,0,0,0,0.13,21.51,90.23,82.56])
            }
    elif z == 6:
        bias_table = {
            "hour": np.array(list(range(6,23)) + list(range(0,6))),
            "temp_bias": np.array([-0.67,-0.6,-0.57,-0.49,-0.2,1.14,1.86,1.74,1.26,0.74,0.28,-0.15,-0.5,-0.71,-0.81,-0.91,-1.11,-1.51,-1.99,-1.44,-1.07,-0.9,-0.79,-0.7,-0.67]) / 1.8,   #nums displayed are temp differences in F -> change to temp diff in K
            "dew_bias": np.array([-1.71,-1.68,-1.68,-1.67,-1.5,-0.93,-1.3,-1.85,-2.18,-2.35,-2.4,-2.4,-2.3,-2.19,-2.13,-1.96,-1.79,-1.56,-1.78,-1.76,-1.77,-1.72,-1.68,-1.69,-1.71]) / 1.8,
            "wind_bias": np.array([0.12,0.14,0.14,0.15,0.2,0.26,0.26,0.27,0.27,0.23,0.22,0.25,0.21,0.17,0.09,0.07,-0.01,-0.15,-0.02,0.09,0.1,0.11,0.1,0.12,0.12]),
            "srad_bias": np.array([0,0,0,0.13,21.51,90.23,82.56,35.49,-24.2,-87.89,-160.95,-229.6,-272.25,-287.4,-271.61,-221.73,-154.68,-68.81,-1.86,0.08,0,0,0,0,0])
            }
    for i in range(0, 25):  
        data_RTMA["TMP_P0_L103_GLC0"][:, :, i] = data_RTMA["TMP_P0_L103_GLC0"][:, :, i] + bias_table["temp_bias"][i]
        data_RTMA["DPT_P0_L103_GLC0"][:, :, i] = data_RTMA["DPT_P0_L103_GLC0"][:, :, i] + bias_table["dew_bias"][i]
        data_RTMA["WIND_P0_L103_GLC0"][:, :, i] = data_RTMA["WIND_P0_L103_GLC0"][:, :, i] + bias_table["wind_bias"][i]
    return data_RTMA

def srad_bias(data_srad, z):
    if z == 12:
        bias_table = {
            "hour": np.array(list(range(12,23)) + list(range(0,12))),
            "temp_bias": np.array([1.86,1.74,1.26,0.74,0.28,-0.15,-0.5,-0.71,-0.81,-0.91,-1.11,-1.51,-1.99,-1.44,-1.07,-0.9,-0.79,-0.7,-0.67,-0.6,-0.57,-0.49,-0.2,1.14,1.86]) / 1.8,   #nums displayed are temp differences in F -> change to temp diff in K
            "dew_bias": np.array([-1.3,-1.85,-2.18,-2.35,-2.4,-2.4,-2.3,-2.19,-2.13,-1.96,-1.79,-1.56,-1.78,-1.76,-1.77,-1.72,-1.68,-1.69,-1.71,-1.68,-1.68,-1.67,-1.5,-0.93,-1.3]) / 1.8,
            "wind_bias": np.array([0.26,0.27,0.27,0.23,0.22,0.25,0.21,0.17,0.09,0.07,-0.01,-0.15,-0.02,0.09,0.1,0.11,0.1,0.12,0.12,0.14,0.14,0.15,0.2,0.26,0.26]),
            "srad_bias": np.array([82.56,35.49,-24.2,-87.89,-160.95,-229.6,-272.25,-287.4,-271.61,-221.73,-154.68,-68.81,-1.86,0.08,0,0,0,0,0,0,0,0.13,21.51,90.23,82.56])
            }
    elif z == 6:
        bias_table = {
            "hour": np.array(list(range(6,23)) + list(range(0,6))),
            "temp_bias": np.array([-0.67,-0.6,-0.57,-0.49,-0.2,1.14,1.86,1.74,1.26,0.74,0.28,-0.15,-0.5,-0.71,-0.81,-0.91,-1.11,-1.51,-1.99,-1.44,-1.07,-0.9,-0.79,-0.7,-0.67]) / 1.8,   #nums displayed are temp differences in F -> change to temp diff in K
            "dew_bias": np.array([-1.71,-1.68,-1.68,-1.67,-1.5,-0.93,-1.3,-1.85,-2.18,-2.35,-2.4,-2.4,-2.3,-2.19,-2.13,-1.96,-1.79,-1.56,-1.78,-1.76,-1.77,-1.72,-1.68,-1.69,-1.71]) / 1.8,
            "wind_bias": np.array([0.12,0.14,0.14,0.15,0.2,0.26,0.26,0.27,0.27,0.23,0.22,0.25,0.21,0.17,0.09,0.07,-0.01,-0.15,-0.02,0.09,0.1,0.11,0.1,0.12,0.12]),
            "srad_bias": np.array([0,0,0,0.13,21.51,90.23,82.56,35.49,-24.2,-87.89,-160.95,-229.6,-272.25,-287.4,-271.61,-221.73,-154.68,-68.81,-1.86,0.08,0,0,0,0,0])
            }
    for i in range(0, 25):  
        data_srad[:, :, i] = data_srad[:, :, i] + bias_table["srad_bias"][i]
    return data_srad

#%%% Testing
def data_gen(datatype):
    if datatype == "RTMA":
        keys = ["TMP_P0_L103_GLC0",
                "DPT_P0_L103_GLC0",
                "WIND_P0_L103_GLC0",
                "TCDC_P0_L200_GLC0"]
        data_out = dict.fromkeys(keys)
        
        data_shape = (420, 370, 25)
        data_inside = np.full(data_shape, 100, dtype="float64")
        for key in keys:
            data_out[key] = data_inside
            
    elif datatype == "NDFD2":
        keys = ["TMP_P0_L103_GLC0",
                "DPT_P0_L103_GLC0",
                "WIND_P0_L103_GLC0",
                "TCDC_P0_L1_GLC0"]
        data_out = dict.fromkeys(keys)
        
        data_shape = (420, 370, 46)
        data_inside = np.full(data_shape, 100, dtype="float64")
        for key in keys:
            data_out[key] = data_inside
            
    elif datatype == "NBM": 
        keys = ["TMP_P0_L103_GLC0",
                "DPT_P0_L103_GLC0",
                "WIND_P0_L103_GLC0",
                "TCDC_P0_L1_GLC0"]
        data_out = dict.fromkeys(keys)
        
        data_shape = (420, 370, 18)
        data_inside = np.full(data_shape, 100, dtype="float64")
        for key in keys:
            data_out[key] = data_inside
            
    return data_out