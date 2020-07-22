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
import os
import subprocess
import logging
from datetime import timedelta
log = logging.getLogger(__name__)

from netCDF4 import Dataset

# Function to build data
# written by JAM July 2020
def build_data(date):
  print("Stuff goes here")

def build_data_RTMA(date_list,vdate):
  outfile = "rtma.nc4"
  outpath = "input/"+outfile
  catpath = outfile.replace(outfile,"rtma_cat.grib2")
  smallpath = outfile.replace(outfile,"rtma_small.grib2")
  file_list = []
  for tdate in date_list:
	# Filepath
        file_ym =  tdate.strftime("%Y%m")
        file_ymd =  tdate.strftime("%Y%m%d")
        file_h =  tdate.strftime("%H")
        file_path = "/data/ncep/rtma/"+file_ym+"/"+file_ymd+"/rtma2p5.t"+file_h+"z.2dvaranl_ndfd.grb2"
        if(os.path.isfile(file_path) == True):
          # Add to the file list
          file_list.append(file_path)
  cmd_list = [
	"/usr/bin/cat "+" ".join(file_list)+" > "+catpath+"",
	"/usr/local/bin/wgrib2 "+catpath+" -small_grib 275.0609:285.8117 32.98101:40.3277 "+smallpath+"",
	"/usr/local/bin/wgrib2 "+smallpath+" -s | /usr/bin/egrep '(:TMP:2|:DPT:2|:TCDC:|:WIND:10)' | /usr/local/bin/wgrib2 -i "+smallpath+" -netcdf "+outfile+"",
	"rm -rf "+catpath+" "+smallpath+"",
  ]
  for cmd in cmd_list:
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)

  return outfile

def build_data_NBM(date_list,vdate):
  outfile = "nbm.nc4"
  outpath = "input/"+outfile
  catpath = outfile.replace(outfile,"nbm_cat.grib2")
  smallpath = outfile.replace(outfile,"nbm_small.grib2")
  file_list = []
  vdatetime = datetime.datetime(vdate.year, vdate.month, vdate.day, 6)
  print(vdatetime)
  file_ym =  vdate.strftime("%Y%m")
  file_ymd =  vdate.strftime("%Y%m%d")
  file_h =  vdate.strftime("%H")
  for tdate in date_list:
        # Filepath
        vdatediff = tdate - vdatetime
        vdatediff_hr = int(divmod(vdatediff.total_seconds(),3600)[0])
        vdatediff_fhr = str(vdatediff_hr).zfill(3)
        file_ymdh =  tdate.strftime("%Y%m%d%H")
        print(vdatediff_fhr)
        file_path = "/data/nws/nbm/"+file_ym+"/"+file_ymd+"/06/blend.t06z.master.f"+vdatediff_fhr+".v"+file_ymdh+".co.grib2"
        if(os.path.isfile(file_path) == True):
          # Add to the file list
          file_list.append(file_path)

  cmd_list = [
	"/usr/bin/cat "+" ".join(file_list)+" > "+catpath+"",
	"/usr/local/bin/wgrib2 "+catpath+" -small_grib 275.0609:285.8117 32.98101:40.3277 "+smallpath+"",
	"/usr/local/bin/wgrib2 "+smallpath+" -s | /usr/bin/egrep '(:TMP:2|:DPT:2|:TCDC:|:WIND:10)' | /usr/local/bin/wgrib2 -i "+smallpath+" -netcdf "+outfile+"",
	"rm -rf "+catpath+" "+smallpath+"",
  ]
  for cmd in cmd_list:
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)

  return outfile


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

def file_list(z=6):
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

    files = []
    rtma_files = []
    rtma_dates = []
    n = 0
    while RTMA_start + n*step <= RTMA_end:
        RTMA_date = RTMA_start + n*step
        rtma_dates.append(RTMA_date)
        n += 1

    # Generate
    #build_data_RTMA(rtma_dates, today_d)

    # NBM
    NBM_end = datetime.datetime.combine(fivedays_d, z_duration)
    NBM_start = datetime.datetime.combine(today_d, z1_duration)

    nbm_dates = []
    n = 0
    while NBM_start + n*step <= NBM_end:
        NBM_date = NBM_start + n*step
        nbm_dates.append(NBM_date)
        n += 1

    # Generate
    build_data_NBM(nbm_dates, today_d)

    # NDFD - Get todays
    NDFD2_end = datetime.datetime.combine(fivedays_d, z_duration)
    NDFD2_start = datetime.datetime.combine(today_d, z1_duration)

    NDFD2_dates = [None]*120
    NDFD2_dates_int = [None]*120
    NDFD2_hours = [None]*120
    NDFD2_files = []
    n = 0
    while NDFD2_start + n*step <= NDFD2_end:
        NDFD2_date = NDFD2_start + n*step
        NDFD2_dates[n] = NDFD2_date
        # Filepath
        file_ym =  NDFD2_date.strftime("%Y%m")
        file_ymd =  NDFD2_date.strftime("%Y%m%d")
        file_h =  NDFD2_date.strftime("%H")
        file_path = "/data/ncep/rtma/"+file_ym+"/"+file_ymd+"/rtma2p5.t"+file_h+"z.2dvaranl_ndfd.grb2"
        # Add to the file list
        files.append({
          'source':'NDFD',
          'datetime':NDFD2_date,
          'file':file_path,
        })
        n += 1

    return files


def file_timing(z=6):
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

    return RTMA_dates,RTMA_dates_int,RTMA_hours,NDFD2_dates,NDFD2_dates_int,NDFD2_hours


def timing(z=6,nx=370,ny=420,nt=1):
    RTMA_dates,RTMA_dates_int,RTMA_hours,NDFD2_dates,NDFD2_dates_int,NDFD2_hours =  file_timing(z)
    
    # Convert to 3d grids
    # RTMA
    nt_rtma = len(RTMA_dates_int)
    jday_RTMA = np.broadcast_to(RTMA_dates_int, (nx,ny,nt_rtma))
    hour_RTMA = np.broadcast_to(RTMA_hours, (nx,ny,nt_rtma))
    jday_RTMA = np.swapaxes(jday_RTMA, 0, 2)
    hour_RTMA = np.swapaxes(hour_RTMA, 0, 2)

    # NDFD
    nt_ndfd = len(NDFD2_dates_int)
    jday_NDFD2 = np.broadcast_to(NDFD2_dates_int, (nx,ny,nt_ndfd))
    hour_NDFD2 = np.broadcast_to(NDFD2_hours, (nx,ny,nt_ndfd))
    jday_NDFD2 = np.swapaxes(jday_NDFD2, 0, 2)
    hour_NDFD2 = np.swapaxes(hour_NDFD2, 0, 2)

    jday_NDFD = np.concatenate((jday_NDFD2[0:36,:,:],
                                jday_NDFD2[38::3,:,:]),
                               axis=0)

    hour_NDFD = np.concatenate((hour_NDFD2[0:36,:,:],
                                hour_NDFD2[38::3,:,:]),
                               axis=0)

    return jday_RTMA, hour_RTMA, jday_NDFD, hour_NDFD, jday_NDFD2, hour_NDFD2

#%%% Dimension Imports. Written by JAM Jul 2020
def import_dims(filename):
    infile = Dataset(filename, "a", format="NETCDF4")
    nx = infile.dimensions['x'].size
    ny = infile.dimensions['y'].size
    nt = infile.dimensions['time'].size
    infile.close()
    return nx, ny, nt

def import_times(filename):
    infile = Dataset(filename, "a", format="NETCDF4")
    vtime = infile['time'].getncattr("reference_date")
    times = infile['time'][:]
    infile.close()
    return vtime, times


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
        #if(data_d[v].ndim == 3):
	#  # JUST DO 2 HOURS FOR NOW FOR TESTING
        #  data_d[v] = data_d[v][:2,:,:] 
        data_d[v].fill_value = np.nan
        unit_d[v] = rootgrp[v].getncattr("units")
        if '_FillValue' in rootgrp[v].ncattrs():
          fill_d[v] = rootgrp[v].getncattr("_FillValue")
        else:
          fill_d[v] = np.nan
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
        #data_d[v] = rootgrp[v][:]
        #if(data_d[v].ndim == 3):
        #  data_d[v] = data_d[v][:46,:,:] 
	#  ## JUST DO 2 HOURS FOR NOW FOR TESTING
        #  #data_d[v] = data_d[v][:2,:,:] 
        data_d[v].fill_value = np.nan
        if '_FillValue' in rootgrp[v].ncattrs():
          fill_d[v] = rootgrp[v].getncattr("_FillValue")
        else:
          fill_d[v] = np.nan
        if 'units' in rootgrp[v].ncattrs():
          fill_d[v] = rootgrp[v].getncattr("units")
        else:
          unit_d[v] = np.nan
    return vars_d, data_d, unit_d, fill_d 

def small_import(filename):
    rootgrp = Dataset(filename, "a", format="NETCDF4")
    vs = [*rootgrp.variables.keys()]
    vars_d = dict.fromkeys(vs)
    data_d = dict.fromkeys(vs)
    for v in vs:
        vars_d[v] = rootgrp[v]
        data_d[v] = rootgrp[v][:]
        #if(data_d[v].ndim == 3):
	#  # JUST DO 2 HOURS FOR NOW FOR TESTING
        #  data_d[v] = data_d[v][:2,:,:] 
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
            "dew_bias": np.array([-1.71,-1.68,-1.68,-1.67,-1.5,-0.93,-1.3,-1.85,-2.18,-2.35,-2.4,-2.4,-2.3,-2.19,-2.13,-1.96,-1.79,-1.56,-1.78,-1.76,-1.77,-1.72,-1.68,-1.69,-1.71]) / 1.8, #nums displayed are temp DeltaF, change to DeltaC/DeltaK
            "wind_bias": np.array([0.12,0.14,0.14,0.15,0.2,0.26,0.26,0.27,0.27,0.23,0.22,0.25,0.21,0.17,0.09,0.07,-0.01,-0.15,-0.02,0.09,0.1,0.11,0.1,0.12,0.12]),
            "srad_bias": np.array([0,0,0,0.13,21.51,90.23,82.56,35.49,-24.2,-87.89,-160.95,-229.6,-272.25,-287.4,-271.61,-221.73,-154.68,-68.81,-1.86,0.08,0,0,0,0,0])
            }
    ntime = data_RTMA["TMP_2maboveground"].shape[0]
    print("ntime for RTMA = ",ntime)
    #print(data_RTMA["TMP_2maboveground"][0,100,100])
    for i in range(0, ntime):  
        data_RTMA["TMP_2maboveground"][i,:,:] = data_RTMA["TMP_2maboveground"][i,:,:] + bias_table["temp_bias"][i]
        data_RTMA["DPT_2maboveground"][i,:,:] = data_RTMA["DPT_2maboveground"][i,:,:] + bias_table["dew_bias"][i]
        data_RTMA["WIND_10maboveground"][i,:,:] = data_RTMA["WIND_10maboveground"][i,:,:] + bias_table["wind_bias"][i]
    #print(data_RTMA["TMP_2maboveground"][0,100,100])
    return data_RTMA

def NDFD_bias(data_NDFD, z):
    bias_table = {
        "hour": np.array(list(range(0,23))),
        "temp_bias": np.array([1.5,1.3,1,0.8,0.8,0.7,0.8,0.6,0.5,0.5,0.4,0,-0.7,-1,-0.8,-0.5,-0.3,0,0.3,0.5,0.7,0.8,0.9,1.3]),   #nums displayed are temp differences in F -> change to temp diff in K
        "dew_bias": np.array([0.7,0.8,0.8,0.8,0.8,0.9,0.9,0.9,0.9,0.9,1,0.8,0.4,0.3,0.6,0.9,1,1,1,1,1,1,0.8,0.6]),
        "wind_bias": np.array([0.5,0.7,0.6,0.6,0.6,0.5,0.4,0.5,0.4,0.4,0.4,0.4,0,0.2,0.3,0.3,0.5,0.5,0.2,0.6,0.7,0.6,0.8,0.9]),
        "srad_bias": np.array([0.2,-0.2,0,0,0,0,0,0,0,0,-8,-64.1,-103.1,-78.9,-34.1,22.9,76.8,133,182.1,189.5,179.9,139.5,91.5,29.5])
        }
    ntime = data_NDFD["TMP_2maboveground"].shape[0]
    print("ntime for NDFD = ",ntime)
    print(data_NDFD["TMP_2maboveground"][0,100,100])
    for i in range(0, ntime):  
        data_NDFD["TMP_2maboveground"][i,:,:] = data_NDFD["TMP_2maboveground"][i,:,:] + bias_table["temp_bias"][i]
        data_NDFD["DPT_2maboveground"][i,:,:] = data_NDFD["DPT_2maboveground"][i,:,:] + bias_table["dew_bias"][i]
        data_NDFD["WIND_10maboveground"][i,:,:] = data_NDFD["WIND_10maboveground"][i,:,:] + bias_table["wind_bias"][i]
    print(data_NDFD["TMP_2maboveground"][0,100,100])
    return data_NDFD

def NBM_bias(data_NBM, z):
    bias_table = {
        "hour": np.array(list(range(0,23))),
        "temp_bias": np.array([1,0.7,0.5,0.4,0.4,0.4,0.3,0.3,0.3,0.2,0.1,-0.4,-1,-1,-0.8,-0.5,-0.3,-0.1,0,0.3,0.3,0.4,0.6,1.2]),   #nums displayed are temp differences in F -> change to temp diff in K
        "dew_bias": np.array([0.7,0.9,0.9,0.9,1,0.9,0.9,0.9,0.9,0.9,0.8,0.6,0.4,0.8,1.1,1.2,1.3,1.3,1.1,1.2,1.1,1,0.8,0.8]),
        "wind_bias": np.array([0,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.2,-0.3,-0.3,-0.2,-0.2,-0.2,-0.2,-0.3,-0.2,-0.1,-0.1,0,0.1]),
        "srad_bias": np.array([0.2,-0.2,0,0,0,0,0,0,0,0,-8,-64.1,-101.7,-74.1,-23.1,39.7,99,155.5,197.2,208,193.6,150.5,95.8,30.5])
        }
    ntime = data_NBM["TMP_2maboveground"].shape[0]
    print("ntime for NBM = ",ntime)
    print(data_NBM["TMP_2maboveground"][:,0,0])
    times = data_NBM["TMP_2maboveground"][:,0,0]
    #for i,t in enumerate(times):  
    for i in range(0, ntime):  
        data_NBM["TMP_2maboveground"][i,:,:] = data_NBM["TMP_2maboveground"][i,:,:] + bias_table["temp_bias"][i]
        data_NBM["DPT_2maboveground"][i,:,:] = data_NBM["DPT_2maboveground"][i,:,:] + bias_table["dew_bias"][i]
        data_NBM["WIND_10maboveground"][i,:,:] = data_NBM["WIND_10maboveground"][i,:,:] + bias_table["wind_bias"][i]
    return data_NBM

def srad_bias(data_srad, z=12, dataset="srad_bias_RTMA"):
    if z == 12:
        bias_table = {
            "srad_bias_RTMA": np.array([82.56,35.49,-24.2,-87.89,-160.95,-229.6,-272.25,-287.4,-271.61,-221.73,-154.68,-68.81,-1.86,0.08,0,0,0,0,0,0,0,0.13,21.51,90.23,82.56]),
            "srad_bias_NDFD": np.array([0.2,-0.2,0,0,0,0,0,0,0,0,-8,-64.1,-103.1,-78.9,-34.1,22.9,76.8,133,182.1,189.5,179.9,139.5,91.5,29.5]),
            "srad_bias_NBM": np.array([0.2,-0.2,0,0,0,0,0,0,0,0,-8,-64.1,-101.7,-74.1,-23.1,39.7,99,155.5,197.2,208,193.6,150.5,95.8,30.5])
            }
    elif z == 6:
        bias_table = {
            "srad_bias_RTMA": np.array([0,0,0,0.13,21.51,90.23,82.56,35.49,-24.2,-87.89,-160.95,-229.6,-272.25,-287.4,-271.61,-221.73,-154.68,-68.81,-1.86,0.08,0,0,0,0,0])
            }

    ntime = data_srad.shape[0]
    for i in range(0, ntime):  
        data_srad[i,:,:] = data_srad[i,:,:] + bias_table[dataset][i]
    return data_srad

#%%% Testing
def data_gen(datatype):
    if datatype == "RTMA":
        keys = ["TMP_2maboveground",
                "DPT_2maboveground",
                "WIND_10maboveground",
                "TCDC_entireatmosphere_consideredasasinglelayer_"]
        data_out = dict.fromkeys(keys)
        
        data_shape = (25,370,420)
        data_inside = np.full(data_shape, 100, dtype="float64")
        for key in keys:
            data_out[key] = data_inside
            
    elif datatype == "NDFD2":
        keys = ["TMP_2maboveground",
                "DPT_2maboveground",
                "WIND_10maboveground",
                "TCDC_P0_L1_GLC0"]
        data_out = dict.fromkeys(keys)
        
        data_shape = (46,370,420)
        data_inside = np.full(data_shape, 100, dtype="float64")
        for key in keys:
            data_out[key] = data_inside
            
    elif datatype == "NBM": 
        keys = ["TMP_2maboveground",
                "DPT_2maboveground",
                "WIND_10maboveground",
                "TCDC_P0_L1_GLC0"]
        data_out = dict.fromkeys(keys)
        
        data_shape = (18,370,420)
        data_inside = np.full(data_shape, 100, dtype="float64")
        for key in keys:
            data_out[key] = data_inside
            
    return data_out
