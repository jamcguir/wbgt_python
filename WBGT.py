 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Gray Martin
"""

#%% Imports
import utilities
import equations
import numpy as np
import pandas as pd
import geopandas as gpd

#%% Setup
TESTMODE = True
LOGGING = True

Z = 6
STATELIST = ["North Carolina", "Virginia"]

np.warnings.filterwarnings('ignore')

#%%% Longitude and Latitude Data
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
statelist = STATELIST
statemasks, statelabels, states = utilities.to_state(gdf, statelist=statelist)

dimmasks = [None]*len(statemasks)
combined_mask = np.full((420,370,1), False)
for row in range(0, len(statemasks)):
    statemask = statemasks[row]
    mask = statemask.to_numpy()
    dimmasks[row] = mask.reshape(420,370,1)
    combined_mask = np.logical_or(combined_mask, dimmasks[row])

#%%% Mask Application
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


#%% Data Imports
if not TESTMODE:
    vars_RTMA, data_RTMA, unit_RTMA, fill_RTMA = utilities.RTMA_import("resources/WBGT_RTMA.nc")
    vars_NDFD2, data_NDFD2, unit_NDFD2, fill_NDFD2 = utilities.NDFD2_import("resources/WBGT_NDFD.nc")
    vars_NBM, data_NBM = utilities.small_import("resources/WBGT_NBM.nc4")
elif TESTMODE:
    data_RTMA = utilities.data_gen("RTMA")
    data_NDFD2 = utilities.data_gen("NDFD2")
    data_NBM = utilities.data_gen("NBM")

vars_elev, data_elev = utilities.small_import("resources/elevation_regrid_NCVA.nc")

#%%% RTMA Bias Correction
data_RTMA = utilities.RTMA_bias(data_RTMA, z=Z)

#%%% Wind Speed Correction
data_RTMA["WIND_P0_L103_GLC0"] = np.where(data_RTMA["WIND_P0_L103_GLC0"] < 0.5, 
                                          0.5, data_RTMA["WIND_P0_L103_GLC0"])
data_NDFD2["WIND_P0_L103_GLC0"] = np.where(data_NDFD2["WIND_P0_L103_GLC0"] < 0.5, 
                                           0.5, data_NDFD2["WIND_P0_L103_GLC0"])
data_NBM["WIND_P0_L103_GLC0"] = np.where(data_NBM["WIND_P0_L103_GLC0"] < 0.5, 
                                         0.5, data_NBM["WIND_P0_L103_GLC0"])

#%%% Name Variables
temp_RTMA = data_RTMA["TMP_P0_L103_GLC0"]
dewp_RTMA = data_RTMA["DPT_P0_L103_GLC0"]
wind_RTMA = data_RTMA["WIND_P0_L103_GLC0"]
cldc_RTMA = data_RTMA["TCDC_P0_L200_GLC0"]

temp_NDFD2 = data_NDFD2["TMP_P0_L103_GLC0"]
dewp_NDFD2 = data_NDFD2["DPT_P0_L103_GLC0"]
wind_NDFD2 = data_NDFD2["WIND_P0_L103_GLC0"]
cldc_NDFD2 = data_NDFD2["TCDC_P0_L1_GLC0"]

temp_NBM = data_NBM["TMP_P0_L103_GLC0"]
dewp_NBM = data_NBM["DPT_P0_L103_GLC0"]
wind_NBM = data_NBM["WIND_P0_L103_GLC0"]
cldc_NBM = data_NBM["TCDC_P0_L1_GLC0"]

# Increase the speed here
elev = np.atleast_3d(data_elev["var"])
elev = np.swapaxes(elev, 0, 1)

#%%% Combine NDFD and NBM
temp_mix = np.concatenate((temp_NDFD2, temp_NBM), axis=2)
dewp_mix = np.concatenate((dewp_NDFD2, dewp_NBM), axis=2)
wind_mix = np.concatenate((wind_NDFD2, wind_NBM), axis=2)
cldc_mix = np.concatenate((cldc_NDFD2, cldc_NBM), axis=2)

#%%% Mask Arrays
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
t_conv = 273.15

temp_RTMA -= t_conv
dewp_RTMA -= t_conv
cldc_RTMA = cldc_RTMA/100.0

temp_mix -= t_conv
dewp_mix -= t_conv
cldc_mix = cldc_mix/100.0

#%% Solar Intensity
#%%% Calculate RH
rh_RTMA = equations.rh_calc(temp_RTMA, dewp_RTMA)
rh_mix = equations.rh_calc(temp_mix, dewp_mix)

#%%% Solar Radiation
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
mshd_RTMA = np.where((hour_RTMA >= 10) & (hour_RTMA <= 14), True, False)
mshd_mix = np.where((hour_NDFD >= 10) & (hour_NDFD <= 14), True, False)

sun_RTMA = np.where(mshd_RTMA, shd_RTMA, sun_RTMA)
act_RTMA = np.where(mshd_RTMA, shd_RTMA, act_RTMA)

sun_NDFD = np.where(mshd_mix, shd_NDFD, sun_NDFD)
act_NDFD = np.where(mshd_mix, shd_NDFD, act_NDFD)

#%%% Theoretical Maximum Solar Radiation
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
fdb_sun_RTMA, fdif_sun_RTMA = equations.direct_diffuse(zenith_RTMA, starsun_RTMA)
fdb_shd_RTMA, fdif_shd_RTMA = equations.direct_diffuse(zenith_RTMA, starshd_RTMA)
fdb_act_RTMA, fdif_act_RTMA = equations.direct_diffuse(zenith_RTMA, staract_RTMA)

fdb_sun_NDFD, fdif_sun_NDFD = equations.direct_diffuse(zenith_NDFD, starsun_NDFD)
fdb_shd_NDFD, fdif_shd_NDFD = equations.direct_diffuse(zenith_NDFD, starshd_NDFD)
fdb_act_NDFD, fdif_act_NDFD = equations.direct_diffuse(zenith_NDFD, staract_NDFD) 

#%% Wind Speed   
#%%% Estimate Stability Class
stabt_sun_RTMA = equations.stability(nght_RTMA, wind_RTMA, sun_RTMA)
stabt_act_RTMA = equations.stability(nght_RTMA, wind_RTMA, act_RTMA)
stabt_shd_RTMA = equations.stability(nght_RTMA, wind_RTMA, shd_RTMA)

stabt_sun_NDFD = equations.stability(nght_mix, wind_mix, sun_NDFD)
stabt_act_NDFD = equations.stability(nght_mix, wind_mix, act_NDFD)
stabt_shd_NDFD = equations.stability(nght_mix, wind_mix, shd_NDFD)

#%%% Estimate Wind Speed
est_speed_sun_RTMA = equations.est_wind_speed(wind_RTMA, stabt_sun_RTMA)
est_speed_act_RTMA = equations.est_wind_speed(wind_RTMA, stabt_act_RTMA)
est_speed_shd_RTMA = equations.est_wind_speed(wind_RTMA, stabt_shd_RTMA)