 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Gray Martin
"""

#%% Imports
import logging
import utilities
import equations
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

#%%% Longitude and Latitude 
log.info("Longitude and Latitude")
with open("resources/lonlat.csv") as read_file:
    lonlat = pd.read_csv(read_file)

# Set up arrays for both longitude and latitude
lon = lonlat.lon.to_numpy()
lon_array = lon.reshape(420,370,1)

lat = lonlat.lat.to_numpy()
lat_array = lat.reshape(420,370,1)

# Create "stacked" arrays spanning across time scales
lon_RTMA = np.repeat(lon_array, 25, axis=2)
lat_RTMA = np.repeat(lat_array, 25, axis=2)

lon_NDFD = np.repeat(lon_array, 64, axis=2)
lat_NDFD = np.repeat(lat_array, 64, axis=2)

lon_NDFD2 = np.repeat(lon_array, 120, axis=2)
lat_NDFD2 = np.repeat(lat_array, 120, axis=2)

# Create GeoPandas DataFrame containing latitude and longitude
gdf = gpd.GeoDataFrame(lonlat, 
                       geometry=gpd.points_from_xy(lonlat.lon-360, lonlat.lat))

#%%% Region Masking
log.info("Region Masking")
statemasks, statelabels, states = utilities.to_state(gdf, statelist=statelist)

dimmasks = [None]*len(statemasks)
combined_mask = np.full((420,370,1), False)
for row in range(0, len(statemasks)):
    statemask = statemasks[row]
    mask = statemask.to_numpy()
    dimmasks[row] = mask.reshape(420,370,1)
    combined_mask = np.logical_or(combined_mask, dimmasks[row])

#%%% Mask Application
log.info("Mask Application")
# Apply masks to longitude and latitude arrays
lon_array_masked = np.where(combined_mask, lon_array, np.nan)
lat_array_masked = np.where(combined_mask, lat_array, np.nan)

# Create "stacked" arrays spanning across time scales
lon_mask_RTMA = np.repeat(lon_array_masked, 25, axis=2)
lat_mask_RTMA = np.repeat(lat_array_masked, 25, axis=2)

lon_mask_NDFD = np.repeat(lon_array_masked, 64, axis=2)
lat_mask_NDFD = np.repeat(lat_array_masked, 64, axis=2)

lon_mask_NDFD2 = np.repeat(lon_array_masked, 120, axis=2)
lat_mask_NDFD2 = np.repeat(lat_array_masked, 120, axis=2)
    
#%% Solar Calculations
#%%% Time
log.info("Timing")
jday_RTMA, hour_RTMA, jday_NDFD, hour_NDFD, jday_NDFD2, hour_NDFD2 = utilities.timing(z=Z)
jday_RTMA_mask = np.where(combined_mask, jday_RTMA, np.nan)
jday_NDFD_mask = np.where(combined_mask, jday_NDFD, np.nan)
jday_NDFD2_mask = np.where(combined_mask, jday_NDFD2, np.nan)

zenith_RTMA, zenith_NDFD, zenith_NDFD2 = equations.solar_calc(lat_mask_RTMA, 
                                                              lon_mask_RTMA, 
                                                              jday_RTMA_mask, 
                                                              hour_RTMA, 
                                                              lat_mask_NDFD2, 
                                                              lon_mask_NDFD2, 
                                                              jday_NDFD2_mask, 
                                                              hour_NDFD2)

# Restrict data to 2 times
log.info("Restricting items for two times FOR TESTING")
jday_RTMA = jday_RTMA[:,:,:2]
jday_RTMA_mask = jday_RTMA_mask[:,:,:2]
hour_RTMA = hour_RTMA[:,:,:2]
zenith_RTMA = zenith_RTMA[:,:,:2]
jday_NDFD = jday_NDFD[:,:,:2]
jday_NDFD_mask = jday_NDFD_mask[:,:,:2]
hour_NDFD = hour_NDFD[:,:,:2]
zenith_NDFD = zenith_NDFD[:,:,:2]
jday_NDFD2 = jday_NDFD2[:,:,:2]
jday_NDFD2_mask = jday_NDFD2_mask[:,:,:2]
hour_NDFD2 = hour_NDFD2[:,:,:2]
zenith_NDFD2 = zenith_NDFD2[:,:,:2]

#%% Data Imports
log.info("Data Loading")
if not TESTMODE:
    #v1
    #vars_RTMA, data_RTMA, unit_RTMA, fill_RTMA = utilities.RTMA_import("resources/WBGT_RTMA.nc")
    #vars_NDFD2, data_NDFD2, unit_NDFD2, fill_NDFD2 = utilities.NDFD2_import("resources/WBGT_NDFD.nc")
    #vars_NBM, data_NBM = utilities.small_import("resources/WBGT_NBM.nc4")
    #v2 - Convert to NC4 before
    log.info("Loading Real Data")
    vars_RTMA, data_RTMA, unit_RTMA, fill_RTMA = utilities.RTMA_import("resources/rtma.nc")
    vars_NDFD2, data_NDFD2, unit_NDFD2, fill_NDFD2 = utilities.NDFD2_import("resources/ndfd.nc")
    #vars_NBM, data_NBM = utilities.small_import("resources/nbm.nc")
    vars_NBM, data_NBM,unit_NBM,fill_NBM = utilities.NDFD2_import("resources/nbm.nc")
    #v3 - Native dataset format
    #vars_RTMA, data_RTMA, unit_RTMA, fill_RTMA = utilities.RTMA_import_grib("resources/rtma_combo.grb2")
    #print(data_RTMA)

elif TESTMODE:
    log.info("Loading Test Data")
    data_RTMA = utilities.data_gen("RTMA")
    data_NDFD2 = utilities.data_gen("NDFD2")
    data_NBM = utilities.data_gen("NBM")

vars_elev, data_elev = utilities.small_import("resources/elevation_regrid_NCVA.nc")

##%%% RTMA Bias Correction
log.info("Skipping RTMA Bias Correction since it doesn't work")
#log.info("RTMA Bias Correction")
#data_RTMA = utilities.RTMA_bias(data_RTMA, z=Z)
##data_RTMA = utilities.RTMA_grib_bias(data_RTMA, z=Z)

#%%% Wind Speed Correction
log.info("Wind Speed Correction")
data_RTMA["WIND_10maboveground"] = np.where(data_RTMA["WIND_10maboveground"] < 0.5, 
                                          0.5, data_RTMA["WIND_10maboveground"])
data_NDFD2["WIND_10maboveground"] = np.where(data_NDFD2["WIND_10maboveground"] < 0.5, 
                                           0.5, data_NDFD2["WIND_10maboveground"])
#data_NDFD2["WIND_10maboveground"] = np.where(data_NDFD2["WIND_10maboveground"] < 0.5, 
#                                           0.5, data_NDFD2["WIND_10maboveground"])
data_NBM["WIND_10maboveground"] = np.where(data_NBM["WIND_10maboveground"] < 0.5, 
                                         0.5, data_NBM["WIND_10maboveground"])

#%%% Name Variables
log.info("Variable Refactoring")
temp_RTMA = data_RTMA["TMP_2maboveground"]
dewp_RTMA = data_RTMA["DPT_2maboveground"]
wind_RTMA = data_RTMA["WIND_10maboveground"]
#cldc_RTMA = data_RTMA["TCDC_surface"]
cldc_RTMA = data_RTMA["TCDC_entireatmosphere_consideredasasinglelayer_"]

temp_NDFD2 = data_NDFD2["TMP_2maboveground"]
dewp_NDFD2 = data_NDFD2["DPT_2maboveground"]
wind_NDFD2 = data_NDFD2["WIND_10maboveground"]
cldc_NDFD2 = data_NDFD2["TCDC_surface"]

temp_NBM = data_NBM["TMP_2maboveground"]
dewp_NBM = data_NBM["DPT_2maboveground"]
wind_NBM = data_NBM["WIND_10maboveground"]
cldc_NBM = data_NBM["TCDC_surface"]

# Increase the speed here
elev = np.atleast_3d(data_elev["var"])
elev = np.swapaxes(elev, 0, 1)

##%%% Combine NDFD and NBM
#log.info("Combine NDFD and NBM datasets")
#temp_mix = np.concatenate((temp_NDFD2, temp_NBM), axis=2)
#dewp_mix = np.concatenate((dewp_NDFD2, dewp_NBM), axis=2)
#wind_mix = np.concatenate((wind_NDFD2, wind_NBM), axis=2)
#cldc_mix = np.concatenate((cldc_NDFD2, cldc_NBM), axis=2)
## DON'T COMBINE FOR NOW
log.info("DO NOTCombine NDFD and NBM datasets. Using NBM instead")
temp_mix =  temp_NBM
dewp_mix =  dewp_NBM
wind_mix =  wind_NBM
cldc_mix =  cldc_NBM

#%%% Mask Arrays
log.info("Mask Imported Arrays")
temp_RTMA = np.where(combined_mask, temp_RTMA, np.nan)
dewp_RTMA = np.where(combined_mask, dewp_RTMA, np.nan)
wind_RTMA = np.where(combined_mask, wind_RTMA, np.nan)
cldc_RTMA = np.where(combined_mask, cldc_RTMA, np.nan)
elev_RTMA = np.where(combined_mask, elev, np.nan)

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
nght_mix = np.where((hour_NDFD <= 10) & (hour_NDFD >= 0), 0, 1)

sr_RTMA = equations.solar_rad(jday_RTMA, hour_RTMA, 
                              lat_RTMA, lon_RTMA, 
                              zenith_RTMA, elev_RTMA)

sr_RTMA = np.where(combined_mask, sr_RTMA, np.nan)*nght_RTMA

sun_RTMA = utilities.srad_bias(sr_RTMA, z=Z)
shd_RTMA = utilities.srad_bias(sr_RTMA*(1-0.75*(1**3.4)), z=Z)
act_RTMA = utilities.srad_bias(sr_RTMA*(1-0.75*(np.power(cldc_RTMA, 3.4))), z=Z)

sr_NDFD = equations.solar_rad(jday_NDFD, hour_NDFD, 
                              lat_NDFD, lon_NDFD, 
                              zenith_NDFD, elev_mix)

sr_NDFD = np.where(combined_mask, sr_NDFD, np.nan)*nght_mix

sun_NDFD = sr_NDFD
shd_NDFD = sr_NDFD*(1-0.75*(1**3.4))
act_NDFD = sr_NDFD*(1-0.75*(np.power(cldc_mix, 3.4)))

#%%% Morning Shade
log.info("Calculate Morning Shade")
mshd_RTMA = np.where((hour_RTMA >= 10) & (hour_RTMA <= 14), True, False)
mshd_mix = np.where((hour_NDFD >= 10) & (hour_NDFD <= 14), True, False)

sun_RTMA = np.where(mshd_RTMA, shd_RTMA, sun_RTMA)
act_RTMA = np.where(mshd_RTMA, shd_RTMA, act_RTMA)

sun_NDFD = np.where(mshd_mix, shd_NDFD, sun_NDFD)
act_NDFD = np.where(mshd_mix, shd_NDFD, act_NDFD)

#%%% Theoretical Maximum Solar Radiation
log.info("Calculate Maximum Solar Radiation")
smax_RTMA = equations.solar_max(jday_RTMA_mask, zenith_RTMA)
smax_NDFD = equations.solar_max(jday_NDFD_mask, zenith_NDFD)

sun_RTMA = np.where(sun_RTMA >= smax_RTMA, smax_RTMA, sun_RTMA)
shd_RTMA = np.where(shd_RTMA >= smax_RTMA, smax_RTMA, shd_RTMA)
act_RTMA = np.where(act_RTMA >= smax_RTMA, smax_RTMA, act_RTMA)

sun_NDFD = np.where(sun_NDFD >= smax_NDFD, smax_NDFD, sun_NDFD)
shd_NDFD = np.where(shd_NDFD >= smax_NDFD, smax_NDFD, shd_NDFD)
act_NDFD = np.where(act_NDFD >= smax_NDFD, smax_NDFD, act_NDFD)

starsun_RTMA = sun_RTMA/smax_RTMA
starshd_RTMA = shd_RTMA/smax_RTMA
staract_RTMA = act_RTMA/smax_RTMA

starsun_NDFD = sun_NDFD/smax_NDFD
starshd_NDFD = shd_NDFD/smax_NDFD
staract_NDFD = act_NDFD/smax_NDFD

#%%% Diffuse and Direct Solar Radiation
log.info("Calculate Diffuse and Direct Solar Radiation")
fdb_sun_RTMA, fdif_sun_RTMA = equations.direct_diffuse(zenith_RTMA, starsun_RTMA)
fdb_shd_RTMA, fdif_shd_RTMA = equations.direct_diffuse(zenith_RTMA, starshd_RTMA)
fdb_act_RTMA, fdif_act_RTMA = equations.direct_diffuse(zenith_RTMA, staract_RTMA)

fdb_sun_NDFD, fdif_sun_NDFD = equations.direct_diffuse(zenith_NDFD, starsun_NDFD)
fdb_shd_NDFD, fdif_shd_NDFD = equations.direct_diffuse(zenith_NDFD, starshd_NDFD)
fdb_act_NDFD, fdif_act_NDFD = equations.direct_diffuse(zenith_NDFD, staract_NDFD) 

#%% Wind Speed   
#%%% Estimate Stability Class
log.info("Estimate Wind Stability Class")
stabt_sun_RTMA = equations.stability(nght_RTMA, wind_RTMA, sun_RTMA)
stabt_act_RTMA = equations.stability(nght_RTMA, wind_RTMA, act_RTMA)
stabt_shd_RTMA = equations.stability(nght_RTMA, wind_RTMA, shd_RTMA)

stabt_sun_NDFD = equations.stability(nght_mix, wind_mix, sun_NDFD)
stabt_act_NDFD = equations.stability(nght_mix, wind_mix, act_NDFD)
stabt_shd_NDFD = equations.stability(nght_mix, wind_mix, shd_NDFD)

#%%% Estimate Wind Speed
log.info("Estimate Wind Speed")
est_speed_sun_RTMA = equations.est_wind_speed(wind_RTMA, stabt_sun_RTMA)
est_speed_act_RTMA = equations.est_wind_speed(wind_RTMA, stabt_act_RTMA)
est_speed_shd_RTMA = equations.est_wind_speed(wind_RTMA, stabt_shd_RTMA)

est_speed_sun_NDFD = equations.est_wind_speed(wind_mix, stabt_sun_NDFD)
est_speed_act_NDFD = equations.est_wind_speed(wind_mix, stabt_act_NDFD)
est_speed_shd_NDFD = equations.est_wind_speed(wind_mix, stabt_shd_NDFD)

#%% Natural Wet Bulb Temp
#%%% Radiative Heating Switch
log.info("Radiative Heating")
rad_RTMA = np.where(nght_RTMA == 0, 1, 0)
rad_mix = np.where(nght_mix == 0, 1, 0)

#%%% Calculate Wet Bulb Temperature
log.info("Calculate Wet Bulb Temperature for RTMA Dataset")
twb_sun_RTMA = equations.twb(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_sun_RTMA, sun_RTMA, fdb_sun_RTMA, zenith_RTMA*np.pi/180, sr_RTMA)
#print(twb_sun_RTMA)
twb_shade_RTMA = equations.twb(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_shd_RTMA, shd_RTMA, fdb_shd_RTMA, zenith_RTMA*np.pi/180, sr_RTMA)
twb_actual_RTMA = equations.twb(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_act_RTMA, act_RTMA, fdb_act_RTMA, zenith_RTMA*np.pi/180, sr_RTMA)

log.info("Calculate Wet Bulb Temperature for NDFD Dataset")
twb_sun_NDFD = equations.twb(temp_mix, dewp_mix, rh_mix, est_speed_sun_NDFD, sun_NDFD, fdb_sun_NDFD, zenith_NDFD*np.pi/180, sr_NDFD)
twb_shade_NDFD = equations.twb(temp_mix, dewp_mix, rh_mix, est_speed_shd_NDFD, shd_NDFD, fdb_shd_NDFD, zenith_NDFD*np.pi/180, sr_NDFD)
twb_actual_NDFD = equations.twb(temp_mix, dewp_mix, rh_mix, est_speed_act_NDFD, act_NDFD, fdb_act_NDFD, zenith_NDFD*np.pi/180, sr_NDFD)

#%%% Calculate Wet Globe Temperature
log.info("Calculate Wet Globe Temperature")
tglobe_sun_RTMA = equations.tglobe(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_sun_RTMA, sun_RTMA, fdb_sun_RTMA, zenith_RTMA*np.pi/180)
tglobe_shade_RTMA = equations.tglobe(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_shd_RTMA, shd_RTMA, fdb_shd_RTMA, zenith_RTMA*np.pi/180)
tglobe_actual_RTMA = equations.tglobe(temp_RTMA, dewp_RTMA, rh_RTMA, est_speed_act_RTMA, act_RTMA, fdb_act_RTMA, zenith_RTMA*np.pi/180)

tglobe_sun_NDFD = equations.tglobe(temp_mix, dewp_mix, rh_mix, est_speed_sun_NDFD, sun_NDFD, fdb_sun_NDFD, zenith_NDFD*np.pi/180)
tglobe_shade_NDFD = equations.tglobe(temp_mix, dewp_mix, rh_mix, est_speed_shd_NDFD, shd_NDFD, fdb_shd_NDFD, zenith_NDFD*np.pi/180)
tglobe_actual_NDFD = equations.tglobe(temp_mix, dewp_mix, rh_mix, est_speed_act_NDFD, act_NDFD, fdb_act_NDFD, zenith_NDFD*np.pi/180)

#%% Wet Bulb Globe Temperature
#%%% Calculate Wet Bulb Globe Temperature
log.info("Combine Wet Bulb and Wet Globe Temperature")
WBGT_sun_RTMA = 0.7*twb_sun_RTMA + 0.2*tglobe_sun_RTMA + 0.1*temp_RTMA
WBGT_shade_RTMA = 0.7*twb_shade_RTMA + 0.2*tglobe_shade_RTMA + 0.1*temp_RTMA
WBGT_actual_RTMA = 0.7*twb_actual_RTMA + 0.2*tglobe_actual_RTMA + 0.1*temp_RTMA

WBGT_sun_NDFD = 0.7*twb_sun_NDFD + 0.2*tglobe_sun_NDFD + 0.1*temp_mix
WBGT_shade_NDFD = 0.7*twb_shade_NDFD + 0.2*tglobe_shade_NDFD + 0.1*temp_mix
WBGT_actual_NDFD = 0.7*twb_actual_NDFD + 0.2*tglobe_actual_NDFD + 0.1*temp_mix

log.info("Combine RTMA and NDFD Datasets")
WBGT_sun = np.concatenate((WBGT_sun_RTMA*(9/5)+32, WBGT_sun_NDFD*(9/5)+32), axis=2)
WBGT_shade = np.concatenate((WBGT_sun_RTMA*(9/5)+32, WBGT_sun_NDFD*(9/5)+32), axis=2)
WBGT_actual = np.concatenate((WBGT_sun_RTMA*(9/5)+32, WBGT_sun_NDFD*(9/5)+32), axis=2)

#%%% Export Wet Bulb Globe Temperature
log.info("Export NetCDF Version 4")
outfile = Dataset("test.nc", "w", format="NETCDF4")
time = outfile.createDimension("time", None)
lat = outfile.createDimension("lat", None)
lon = outfile.createDimension("lon", None)

# Add the WBGT_SUN
WBGT_sun_var = outfile.createVariable("wbgt_sun","f8",('time','lat','lon'),zlib=True)
WBGT_sun_var[:,:,:] = WBGT_sun[:,:,:]
# Add the WBGT_SHADE
WBGT_shade_var = outfile.createVariable("wbgt_shade","f8",('time','lat','lon'),zlib=True)
WBGT_shade_var[:,:,:] = WBGT_shade[:,:,:]
# Add the WBGT_ACTUAL
WBGT_actual_var = outfile.createVariable("wbgt_actual","f8",('time','lat','lon'),zlib=True)
WBGT_actual_var[:,:,:] = WBGT_actual[:,:,:]

