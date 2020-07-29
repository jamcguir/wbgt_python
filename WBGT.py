 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Gray Martin, July 2020 updates by John McGuire
"""

#%% Imports
import logging
import utilities
import equations
import datetime
import numpy as np
import pandas as pd
import geopandas as gpd
from netCDF4 import Dataset

#%% Setup
#%%% Logging
log = logging.getLogger(__name__)
log_format = "%(levelname)-7.7s - %(message)-60s\t - [%(lineno)d] %(module)s.%(funcName)s"
logging.basicConfig(level=logging.INFO, format=log_format)

#%%% Development
#TESTMODE = True
TESTMODE = False
np.warnings.filterwarnings('ignore')

#%%% Constants
Z = 6
statelist = ["North Carolina", "Virginia"]

# Get the file list
vdate = datetime.datetime.utcnow()
if(vdate.hour <= 6):
  vdate = datetime.datetime.utcnow() - datetime.timedelta(days=1)
print(vdate)
vdate_ymd =  vdate.strftime("%Y%m%d")
log.info("Generate files for"+vdate_ymd)#files = utilities.build_input_data(vdate, Z)
files = utilities.build_input_data(vdate, Z)
log.info("Done generate files")#files = utilities.build_input_data(vdate, Z)

# File input
nx,ny,nt = utilities.import_dims("input/rtma.nc")
vtime_rtma,times_rtma = utilities.import_times("input/rtma.nc")
vtime_ndfd,times_ndfd = utilities.import_times("input/ndfd.nc")
vtime_nbm,times_nbm = utilities.import_times("input/nbm.nc")

lat_array,lon_array =  utilities.import_latlon("input/rtma.nc")

# The Mixture Configuration is as follows:
#
# "We use the data [NDFD] until it reaches the 6-hourly mark, 
# then transition to NBM data to keep data at hourly then 3-hourly intervals.""
#
# According to this input data, it appears that it should be mixed as following
# 1-46, NDFD
# 47-64, NBM
vtime_mix=vtime_ndfd
times_mix=np.concatenate([times_ndfd[:46],times_nbm[46:]],axis=0)

RTMA_dates,RTMA_dates_int,RTMA_hours,mix_dates,mix_dates_int,mix_hours =  utilities.file_timing(Z)
#print(RTMA_dates)
#print(mix_dates)

source_rtma = ["rtma" for i in range(0,len(times_rtma))]
source_ndfd = ["ndfd" for i in range(0,len(times_ndfd))]
source_nbm = ["nbm" for i in range(0,len(times_nbm))]
source_mix = np.hstack([source_ndfd[:46],source_nbm[46:]]) 

# The final times should in theory be
# 1-25: RTMA
# 26-71: NDFD
# 72-89: NBM
#times = np.hstack((times_rtma, times_nbm))
#times_source = np.hstack((source_rtma, source_nbm))
times = np.hstack((times_rtma, times_mix))
times_source = np.hstack((source_rtma, source_mix))
#print(times)
#print(times_source)

#%%% Longitude and Latitude 
#log.info("Longitude and Latitude")
#with open("resources/lonlat.csv") as read_file:
#    lonlat = pd.read_csv(read_file)
##print(nx,ny,nt)
##quit()
#
## Set up arrays for both longitude and latitude
#lon = lonlat.lon.to_numpy()
#lon_array = lon.reshape(1,ny,nx)
#
#lat = lonlat.lat.to_numpy()
#lat_array = lat.reshape(1,ny,nx)

# Create "stacked" arrays spanning across time scales
lon_RTMA = np.repeat(lon_array, len(times_rtma), axis=0)
lat_RTMA = np.repeat(lat_array, len(times_rtma), axis=0)

lon_NDFD = np.repeat(lon_array, len(times_ndfd), axis=0)
lat_NDFD = np.repeat(lat_array, len(times_ndfd), axis=0)

lon_mix = np.repeat(lon_array, len(times_mix), axis=0)
lat_mix = np.repeat(lat_array, len(times_mix), axis=0)

## Create GeoPandas DataFrame containing latitude and longitude
#gdf = gpd.GeoDataFrame(lonlat, 
#                       geometry=gpd.points_from_xy(lonlat.lon-360, lonlat.lat))
#
##%%% Region Masking
#log.info("Region Masking")
#statemasks, statelabels, states = utilities.to_state(gdf, statelist=statelist)
#
#dimmasks = [None]*len(statemasks)
#combined_mask = np.full((1,ny,nx), False)
#for row in range(0, len(statemasks)):
#    statemask = statemasks[row]
#    mask = statemask.to_numpy()
#    dimmasks[row] = mask.reshape(1,ny,nx)
#    combined_mask = np.logical_or(combined_mask, dimmasks[row])
#
## Keeping mask because elevation data does not encompass entire model domain
log.info("DO NOT Mask Application")
combined_mask = np.full((1,ny,nx), True)

#%%% Mask Application
log.info("Mask Application")
# Apply masks to longitude and latitude arrays
lon_array_masked = np.where(combined_mask, lon_array, np.nan)
lat_array_masked = np.where(combined_mask, lat_array, np.nan)

# Create "stacked" arrays spanning across time scales
lon_mask_RTMA = np.repeat(lon_array_masked, len(times_rtma), axis=0)
lat_mask_RTMA = np.repeat(lat_array_masked, len(times_rtma), axis=0)

lon_mask_NDFD = np.repeat(lon_array_masked, len(times_ndfd), axis=0)
lat_mask_NDFD = np.repeat(lat_array_masked, len(times_ndfd), axis=0)

lon_mask_mix = np.repeat(lon_array_masked, len(times_mix), axis=0)
lat_mask_mix = np.repeat(lat_array_masked, len(times_mix), axis=0)
    
#%% Solar Calculations
#%%% Time
log.info("Timing")
jday_RTMA, hour_RTMA, jday_NDFD, hour_NDFD, jday_mix, hour_mix = utilities.timing(z=Z,nx=nx,ny=ny,nt=nt)
jday_RTMA_mask = np.where(combined_mask, jday_RTMA, np.nan)
jday_NDFD_mask = np.where(combined_mask, jday_NDFD, np.nan)
jday_mix_mask = np.where(combined_mask, jday_mix, np.nan)

zenith_RTMA, zenith_NDFD, zenith_mix = equations.solar_calc(lat_mask_RTMA, 
                                                              lon_mask_RTMA, 
                                                              jday_RTMA_mask, 
                                                              hour_RTMA, 
                                                              lat_mask_mix, 
                                                              lon_mask_mix, 
                                                              jday_mix_mask, 
                                                              hour_mix)

# Restrict data to 2 times
#log.info("Restricting items for two times FOR TESTING")
#jday_RTMA = jday_RTMA[:2,:,:]
#jday_RTMA_mask = jday_RTMA_mask[:2,:,:]
#hour_RTMA = hour_RTMA[:2,:,:]
#zenith_RTMA = zenith_RTMA[:2,:,:]
#jday_NDFD = jday_NDFD[:2,:,:]
#jday_NDFD_mask = jday_NDFD_mask[:2,:,:]
#hour_NDFD = hour_NDFD[:2,:,:]
#zenith_NDFD = zenith_NDFD[:2,:,:]
#jday_mix = jday_mix[:2,:,:]
#jday_mix_mask = jday_mix_mask[:2,:,:]
#hour_mix = hour_mix[:2,:,:]
#zenith_mix = zenith_mix[:2,:,:]

#%% Data Imports
log.info("Data Loading")
if not TESTMODE:
    #v1
    #vars_RTMA, data_RTMA, unit_RTMA, fill_RTMA = utilities.RTMA_import("resources/WBGT_RTMA.nc")
    #vars_NDFD, data_NDFD, unit_NDFD, fill_NDFD = utilities.NDFD_import("resources/WBGT_NDFD.nc")
    #vars_NBM, data_NBM = utilities.small_import("resources/WBGT_NBM.nc4")
    #v2 - Convert to NC4 before
    log.info("Loading Real Data")
    vars_RTMA, data_RTMA, unit_RTMA, fill_RTMA = utilities.RTMA_import("input/rtma.nc")
    vars_NDFD, data_NDFD, unit_NDFD, fill_NDFD = utilities.NDFD_import("input/ndfd.nc") #CHECK
    #vars_NBM, data_NBM = utilities.small_import("input/nbm.nc")
    vars_NBM, data_NBM,unit_NBM,fill_NBM = utilities.NBM_import("input/nbm.nc")
    #v3 - Native dataset format
    #vars_RTMA, data_RTMA, unit_RTMA, fill_RTMA = utilities.RTMA_import_grib("input/rtma_combo.grb2")
    #print(data_RTMA)

elif TESTMODE:
    log.info("Loading Test Data")
    data_RTMA = utilities.data_gen("RTMA")
    data_NDFD = utilities.data_gen("NDFD")
    data_NBM = utilities.data_gen("NBM")

# get the elevation data
#vars_elev, data_elev = utilities.small_import("resources/elevation_regrid_NCVA.nc")
vars_elev, data_elev = utilities.small_import("resources/elevation_sercc.nc")

##%%% RTMA Bias Correction
log.info("SERCC's RTMA Bias Correction")
data_RTMA = utilities.RTMA_bias(data_RTMA, z=Z)
#log.info("SERCC's NDFD Bias Correction")
#data_NDFD = utilities.NDFD_bias(data_NDFD, z=Z)
#log.info("SERCC's NBM Bias Correction")
#data_NBM = utilities.NBM_bias(data_NBM, z=Z)

#%%% Wind Speed Correction
log.info("Wind Speed Correction")
data_RTMA["WIND_10maboveground"] = np.where(data_RTMA["WIND_10maboveground"] < 0.5, 
                                          0.5, data_RTMA["WIND_10maboveground"])
data_NDFD["WIND_10maboveground"] = np.where(data_NDFD["WIND_10maboveground"] < 0.5, 
                                           0.5, data_NDFD["WIND_10maboveground"])
#data_NDFD["WIND_10maboveground"] = np.where(data_NDFD["WIND_10maboveground"] < 0.5, 
#                                           0.5, data_NDFD["WIND_10maboveground"])
data_NBM["WIND_10maboveground"] = np.where(data_NBM["WIND_10maboveground"] < 0.5, 
                                         0.5, data_NBM["WIND_10maboveground"])

#%%% Name Variables
log.info("Variable Refactoring")
lat_RTMA = data_RTMA["latitude"]
lon_RTMA = data_RTMA["longitude"]
temp_RTMA = data_RTMA["TMP_2maboveground"]
dewp_RTMA = data_RTMA["DPT_2maboveground"]
wind_RTMA = data_RTMA["WIND_10maboveground"]
#cldc_RTMA = data_RTMA["TCDC_surface"]
cldc_RTMA = data_RTMA["TCDC_entireatmosphere_consideredasasinglelayer_"]

lat_NDFD = data_NDFD["latitude"]
lon_NDFD = data_NDFD["longitude"]
temp_NDFD = data_NDFD["TMP_2maboveground"]
dewp_NDFD = data_NDFD["DPT_2maboveground"]
wind_NDFD = data_NDFD["WIND_10maboveground"]
cldc_NDFD = data_NDFD["TCDC_surface"]

lat_NBM = data_NBM["latitude"]
lon_NBM = data_NBM["longitude"]
temp_NBM = data_NBM["TMP_2maboveground"]
dewp_NBM = data_NBM["DPT_2maboveground"]
wind_NBM = data_NBM["WIND_10maboveground"]
cldc_NBM = data_NBM["TCDC_surface"]

# Increase the speed here
elev = np.atleast_3d(data_elev["var"])
# Get the elevation
elev = np.swapaxes(elev, 0, 1)
elev = np.swapaxes(elev, 0, 2)

##%%% Combine NDFD and NBM
log.info("Combine NDFD and NBM datasets")
temp_mix = np.concatenate((temp_NDFD, temp_NBM[46:,:,:]), axis=0)
dewp_mix = np.concatenate((dewp_NDFD, dewp_NBM[46:,:,:]), axis=0)
wind_mix = np.concatenate((wind_NDFD, wind_NBM[46:,:,:]), axis=0)
cldc_mix = np.concatenate((cldc_NDFD, cldc_NBM[46:,:,:]), axis=0)

# From this point forward, all calculations should be just 
# 1) RTMA
#    OR
# 2) mix and NOT NDFD/NBM

## DON'T COMBINE FOR NOW
#log.info("DO NOT Combine NDFD and NBM datasets. Using NBM instead")
#temp_mix =  temp_NBM
#dewp_mix =  dewp_NBM
#wind_mix =  wind_NBM
#cldc_mix =  cldc_NBM

#%%% Mask Arrays
log.info("SKIP Mask Imported Arrays")
temp_RTMA = np.where(combined_mask, temp_RTMA, np.nan)
dewp_RTMA = np.where(combined_mask, dewp_RTMA, np.nan)
wind_RTMA = np.where(combined_mask, wind_RTMA, np.nan)
cldc_RTMA = np.where(combined_mask, cldc_RTMA, np.nan)
elev_RTMA = np.where(combined_mask, elev, np.nan)
#
temp_mix = np.where(combined_mask, temp_mix, np.nan)
dewp_mix = np.where(combined_mask, dewp_mix, np.nan)
wind_mix = np.where(combined_mask, wind_mix, np.nan)
cldc_mix = np.where(combined_mask, cldc_mix, np.nan)
elev_mix = np.where(combined_mask, elev, np.nan)

#%%% Unit Conversion
log.info("Convert Temperature Units")
t_conv = 273.15

temp_RTMA -= t_conv
dewp_RTMA -= t_conv
cldc_RTMA = cldc_RTMA/100.0

#cldc_NDFD = cldc_NDFD/100.0

temp_mix -= t_conv
dewp_mix -= t_conv
cldc_mix = cldc_mix/100.0

#%% Solar Intensity
#%%% Calculate Relative Humidity
log.info("Calculate Relative Humidity")
rh_RTMA = equations.rh_calc(temp_RTMA, dewp_RTMA)
rh_mix = equations.rh_calc(temp_mix, dewp_mix)

#%%% Solar Radiation
log.info("Calculate Solar Radiation")
nght_RTMA = np.where((hour_RTMA <= 10) & (hour_RTMA >= 0), 0, 1)
nght_mix = np.where((hour_mix <= 10) & (hour_mix >= 0), 0, 1)

#print("zenith_rtma=",zenith_RTMA.shape)
#print("elev_rtma=",elev_RTMA.shape)
sr_RTMA = equations.solar_rad(jday_RTMA, hour_RTMA, 
                              lat_RTMA, lon_RTMA, 
                              zenith_RTMA, elev_RTMA)

sr_RTMA = np.where(combined_mask, sr_RTMA, np.nan)*nght_RTMA

sun_RTMA = utilities.srad_bias(sr_RTMA, z=Z)
shd_RTMA = utilities.srad_bias(sr_RTMA*(1-0.75*(1**3.4)), z=Z)
act_RTMA = utilities.srad_bias(sr_RTMA*(1-0.75*(np.power(cldc_RTMA, 3.4))), z=Z)


#print("zenith_ndfd=",zenith_mix.shape)
#print("elev_ndfd=",elev_mix.shape)
sr_mix = equations.solar_rad(jday_mix, hour_mix, 
                              lat_mix, lon_mix, 
                              zenith_mix, elev_mix)

sr_mix = np.where(combined_mask, sr_mix, np.nan)*nght_mix

sun_mix = sr_mix
shd_mix = sr_mix*(1-0.75*(1**3.4))
#print("sr_mix = ",sr_mix.shape)
#print("cldc_mix = ",cldc_mix.shape)
act_mix = sr_mix*(1-0.75*(np.power(cldc_mix, 3.4)))

#%%% Morning Shade
log.info("Calculate Morning Shade")
mshd_RTMA = np.where((hour_RTMA >= 10) & (hour_RTMA <= 14), True, False)
mshd_mix = np.where((hour_mix >= 10) & (hour_mix <= 14), True, False)

sun_RTMA = np.where(mshd_RTMA, shd_RTMA, sun_RTMA)
act_RTMA = np.where(mshd_RTMA, shd_RTMA, act_RTMA)

sun_mix = np.where(mshd_mix, shd_mix, sun_mix)
act_mix = np.where(mshd_mix, shd_mix, act_mix)

#%%% Theoretical Maximum Solar Radiation
log.info("Calculate Maximum Solar Radiation")
smax_RTMA = equations.solar_max(jday_RTMA_mask, zenith_RTMA)
smax_mix = equations.solar_max(jday_mix_mask, zenith_mix)

sun_RTMA = np.where(sun_RTMA >= smax_RTMA, smax_RTMA, sun_RTMA)
shd_RTMA = np.where(shd_RTMA >= smax_RTMA, smax_RTMA, shd_RTMA)
act_RTMA = np.where(act_RTMA >= smax_RTMA, smax_RTMA, act_RTMA)

sun_mix = np.where(sun_mix >= smax_mix, smax_mix, sun_mix)
shd_mix = np.where(shd_mix >= smax_mix, smax_mix, shd_mix)
act_mix = np.where(act_mix >= smax_mix, smax_mix, act_mix)

starsun_RTMA = sun_RTMA/smax_RTMA
starshd_RTMA = shd_RTMA/smax_RTMA
staract_RTMA = act_RTMA/smax_RTMA

starsun_mix = sun_mix/smax_mix
starshd_mix = shd_mix/smax_mix
staract_mix = act_mix/smax_mix

#%%% Diffuse and Direct Solar Radiation
log.info("Calculate Diffuse and Direct Solar Radiation")
fdb_sun_RTMA, fdif_sun_RTMA = equations.direct_diffuse(zenith_RTMA, starsun_RTMA)
fdb_shd_RTMA, fdif_shd_RTMA = equations.direct_diffuse(zenith_RTMA, starshd_RTMA)
fdb_act_RTMA, fdif_act_RTMA = equations.direct_diffuse(zenith_RTMA, staract_RTMA)

fdb_sun_mix, fdif_sun_mix = equations.direct_diffuse(zenith_mix, starsun_mix)
fdb_shd_mix, fdif_shd_mix = equations.direct_diffuse(zenith_mix, starshd_mix)
fdb_act_mix, fdif_act_mix = equations.direct_diffuse(zenith_mix, staract_mix) 

#%% Wind Speed   
#%%% Estimate Stability Class
log.info("Estimate Wind Stability Class")
stabt_sun_RTMA = equations.stability(nght_RTMA, wind_RTMA, sun_RTMA)
stabt_act_RTMA = equations.stability(nght_RTMA, wind_RTMA, act_RTMA)
stabt_shd_RTMA = equations.stability(nght_RTMA, wind_RTMA, shd_RTMA)

stabt_sun_mix = equations.stability(nght_mix, wind_mix, sun_mix)
stabt_act_mix = equations.stability(nght_mix, wind_mix, act_mix)
stabt_shd_mix = equations.stability(nght_mix, wind_mix, shd_mix)

#%%% Estimate Wind Speed
log.info("Estimate Wind Speed")
est_speed_sun_RTMA = equations.est_wind_speed(wind_RTMA, stabt_sun_RTMA)
est_speed_act_RTMA = equations.est_wind_speed(wind_RTMA, stabt_act_RTMA)
est_speed_shd_RTMA = equations.est_wind_speed(wind_RTMA, stabt_shd_RTMA)

est_speed_sun_mix = equations.est_wind_speed(wind_mix, stabt_sun_mix)
est_speed_act_mix = equations.est_wind_speed(wind_mix, stabt_act_mix)
est_speed_shd_mix = equations.est_wind_speed(wind_mix, stabt_shd_mix)

#%% Natural Wet Bulb Temp
#%%% Radiative Heating Switch
log.info("Radiative Heating")
rad_RTMA = np.where(nght_RTMA == 0, 1, 0)
rad_mix = np.where(nght_mix == 0, 1, 0)

#%%% Calculate Wet Bulb Temperature
log.info("Calculate Wet Bulb Temperature for RTMA Dataset")
twb_sun_RTMA = equations.twb(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_sun_RTMA, sun_RTMA, fdb_sun_RTMA, np.cos(zenith_RTMA*np.pi/180), rad_RTMA)
#print(twb_sun_RTMA[0,180:250,100])
twb_shade_RTMA = equations.twb(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_shd_RTMA, shd_RTMA, fdb_shd_RTMA, np.cos(zenith_RTMA*np.pi/180), rad_RTMA)
twb_actual_RTMA = equations.twb(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_act_RTMA, act_RTMA, fdb_act_RTMA, np.cos(zenith_RTMA*np.pi/180), rad_RTMA)

log.info("Calculate Wet Bulb Temperature for mix Dataset")
twb_sun_mix = equations.twb(temp_mix, dewp_mix, rh_mix, est_speed_sun_mix, sun_mix, fdb_sun_mix, np.cos(zenith_mix*np.pi/180), rad_mix)
#print(twb_sun_mix[0,180:250,100])
twb_shade_mix = equations.twb(temp_mix, dewp_mix, rh_mix, est_speed_shd_mix, shd_mix, fdb_shd_mix, np.cos(zenith_mix*np.pi/180), rad_mix)
twb_actual_mix = equations.twb(temp_mix, dewp_mix, rh_mix, est_speed_act_mix, act_mix, fdb_act_mix, np.cos(zenith_mix*np.pi/180), rad_mix)

#%%% Calculate Wet Globe Temperature
log.info("Calculate Wet Globe Temperature for RTMA Dataset")
tglobe_sun_RTMA = equations.tglobe(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_sun_RTMA, sun_RTMA, fdb_sun_RTMA, np.cos(zenith_RTMA*np.pi/180))
#print(tglobe_sun_RTMA[0,180:250,100])
tglobe_shade_RTMA = equations.tglobe(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_shd_RTMA, shd_RTMA, fdb_shd_RTMA, zenith_RTMA*np.pi/180)
tglobe_actual_RTMA = equations.tglobe(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_act_RTMA, act_RTMA, fdb_act_RTMA, zenith_RTMA*np.pi/180)

log.info("Calculate Wet Globe Temperature for mix Dataset")
tglobe_sun_mix = equations.tglobe(temp_mix, dewp_mix, rh_mix, est_speed_sun_mix, sun_mix, fdb_sun_mix, zenith_mix*np.pi/180)
#print(tglobe_sun_mix[0,180:250,100])
tglobe_shade_mix = equations.tglobe(temp_mix, dewp_mix, rh_mix, est_speed_shd_mix, shd_mix, fdb_shd_mix, zenith_mix*np.pi/180)
tglobe_actual_mix = equations.tglobe(temp_mix, dewp_mix, rh_mix, est_speed_act_mix, act_mix, fdb_act_mix, zenith_mix*np.pi/180)

#%% Wet Bulb Globe Temperature
#%%% Calculate Wet Bulb Globe Temperature
log.info("Combine Wet Bulb and Wet Globe Temperature for RTMA Dataset")
WBGT_sun_RTMA = 0.7*twb_sun_RTMA + 0.2*tglobe_sun_RTMA + 0.1*temp_RTMA
#print(WBGT_sun_RTMA[0,180:250,100])
WBGT_shade_RTMA = 0.7*twb_shade_RTMA + 0.2*tglobe_shade_RTMA + 0.1*temp_RTMA
WBGT_actual_RTMA = 0.7*twb_actual_RTMA + 0.2*tglobe_actual_RTMA + 0.1*temp_RTMA

log.info("Combine Wet Bulb and Wet Globe Temperature for mix Dataset")
WBGT_sun_mix = 0.7*twb_sun_mix + 0.2*tglobe_sun_mix + 0.1*temp_mix
#print(WBGT_sun_mix[0,180:250,100])
WBGT_shade_mix = 0.7*twb_shade_mix + 0.2*tglobe_shade_mix + 0.1*temp_mix
WBGT_actual_mix = 0.7*twb_actual_mix + 0.2*tglobe_actual_mix + 0.1*temp_mix

log.info("Combine RTMA and mix (NDFD/NBM) Datasets")
#Convert from C to F
WBGT_sun = np.concatenate((WBGT_sun_RTMA*(9/5)+32, WBGT_sun_mix*(9/5)+32), axis=0)
WBGT_shade = np.concatenate((WBGT_sun_RTMA*(9/5)+32, WBGT_sun_mix*(9/5)+32), axis=0)
WBGT_actual = np.concatenate((WBGT_sun_RTMA*(9/5)+32, WBGT_sun_mix*(9/5)+32), axis=0)

#%%% Export Wet Bulb Globe Temperature
log.info("Export NetCDF Version 4")
outfilename = "wbgt_"+vdate_ymd+".nc4"
outfile = Dataset(outfilename, "w", format="NETCDF4")

#outfile.title = 'Wet Bulb Globe Temperature (WBGT) forecast using RTMA, NDFD, and NBM. Written for SERCC by NC SCO'
outfile.institution = "Southeast Regional Climate Center"
outfile.source = "TBD"
outfile.Conventions = 'CF-1.5'
outfile.references = "TBD"
outfile.validtime = vtime_rtma

# Create the dimensions
#lon = outfile.createDimension("lat", nx)
#lat = outfile.createDimension("lon", ny)
#time = outfile.createDimension("time", None)
y = outfile.createDimension("y", ny)
x = outfile.createDimension("x", nx)
t = outfile.createDimension("t", None)

### create time axis
out_time = outfile.createVariable('time', 'f8', ('t'),zlib=True)
out_time.setncatts({
                    'standard_name': u"time",\
                    'long_name': u"time",\
                    #'units':u"Hours Since "+vtime+"",\
                    'units':u"Seconds since 1970-01-01 00:00:00.0 0:00",\
                    'coordinates':u'time',\
                    '_CoordinateAxisType':U'Time',\
                    'calendar':u'gregorian',\
                    'reference_date':u""+vtime_rtma+"",\
})
out_time[:]=times[:]

# create latitude axis
out_lat = outfile.createVariable('latitude', 'f8', ('y','x'),zlib=True)
out_lat.setncatts({
                    'standard_name': u"latitude",\
                    'long_name': u"latitude",\
                    'units':u"degrees_north",\
                    '_CoordinateAxisType':u"Lat",\
#                    'coordinates':u'latitude',\
})
out_lat[:,:]=lat_RTMA[:,:]

# create longitude axis
out_lon = outfile.createVariable('longitude', 'f8', ('y','x'),zlib=True)
out_lon.setncatts({
                    'standard_name': u"longitude",\
                    'long_name': u"longitude",\
                    'units':u"degrees_east",\
                    '_CoordinateAxisType':u"Lon",\
#                    'coordinates':u'longitude',\
})
out_lon[:,:]=lon_RTMA[:,:]-360

### Add the data source
out_source = outfile.createVariable('source', 'S10', ('t'),zlib=True)
out_source.setncatts({
                    'standard_name': u"Data Source Used",\
                    'long_name': u"Data Source Used",\
})
out_source[:]=times_source[:]

# create latitude axis

# Add the WBGT_SUN
log.info("Writing out WBGT Sun to file")
WBGT_sun_var = outfile.createVariable("wbgt_sun","f8",('t','y','x'),zlib=True)
#WBGT_sun_var = outfile.createVariable("wbgt_sun","f8",('y','x','t'),zlib=True)
WBGT_sun_var.setncatts({
                    'long_name': u"Sun WBGT",\
                    'units': u"degF", 'level_desc': u'Surface',\
                    'var_desc': u"Sun WBGT",\
                    'coordinates':u'latitude longitude',\
                    'level_desc':u"Surface",\
                    'min': 0,\
})
#print("var WBGT_sun=",WBGT_sun.shape)
#print("nc WBGT_sun=",WBGT_sun_var.shape)
WBGT_sun_var[:,:,:] = WBGT_sun[:,:,:]
#outfile.variables['wbgt_sun'][:] = WBGT_sun

# Add the WBGT_SHADE
#WBGT_shade_var = outfile.createVariable("wbgt_shade","f8",('lat','lon','time'),zlib=True)
#WBGT_shade_var = outfile.createVariable("wbgt_shade","f8",('y','x','t'),zlib=True)
log.info("Writing out WBGT Shade to file")
WBGT_shade_var = outfile.createVariable("wbgt_shade","f8",('t','y','x'),zlib=True)
WBGT_shade_var.setncatts({
                    'long_name': u"Shade WBGT",\
                    'units': u"degF", 'level_desc': u'Surface',\
                    'var_desc': u"Shade WBGT",\
                    'coordinates':u'latitude longitude',\
                    'level_desc':u"Surface",\
                    'min': 0,\
})
WBGT_shade_var[:,:,:] = WBGT_shade[:,:,:]

# Add the WBGT_ACTUAL
#WBGT_actual_var = outfile.createVariable("wbgt_actual","f8",('y','x','t'),zlib=True)
log.info("Writing out WBGT Actual to file")
WBGT_actual_var = outfile.createVariable("wbgt_actual","f8",('t','y','x'),zlib=True)
WBGT_actual_var.setncatts({
                    'long_name': u"Actual WBGT",\
                    'units': u"degF", 'level_desc': u'Surface',\
                    'var_desc': u"Actual WBGT",\
                    'coordinates':u'latitude longitude',\
                    'level_desc':u"Surface",\
                    'min': 0,\
})
WBGT_actual_var[:,:,:] = WBGT_actual[:,:,:]
outfile.close()
