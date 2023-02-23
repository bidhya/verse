#!/usr/bin/env python
# coding: utf-8

# # Process SEUP NETCDF data [Run this only on  NASA Discover HPC]
# 1. Generate Daily Data by converting 3-hourly SEUP data
# 2. Concatenate monthly files to get water year file
# 3. Combine individual WY nc files files into one file with all varabiles
# 
# # TODO
# - check why this is not exactly 0.05 but rather DX=0.05000000074505806; same issue for DY  
# -- for now fix is to hard-code both DX and DY = 0.05 explicitly
# 
# 
# Sep 30, 2022: Updated to change mm to m for snowfall and swe  
# Feb 21, 2023
# 
# Info: ~ 2.2 million (2,198,176) pixels for North America at 0.05 degree lat/lon
# Time resolution : each file every 3 hours [00, 03, 06, 09, 12, 15, 18, 21 hours]  
# in_vars = ['MODSCAG', 'WRFG', 'WRFP', 'WRFSWE', 'WRFT']  
# 
# ## List of needed variables
# ### Based on list here.
# - Snowfall: Snowf : for now
# - Rainfall: Rainf  (for future)
# - SWE: SWE
# - Air temperature: Tair_f
# - Ground heat flux: Qg
# 
# 
# ## Issues
# - longitude and latitude are both 2-D arrays and embedded as variable. They should each be 1-D coords (and dims)  
#     - probably easier to construct lat/lon manually rather than extracting from 2-D array (but do verify at some point the values are correct)   
# - Extract datetime from filename itself
# - there is difference between groupby and resample operations; though data remains same (for sum) or similar (for mean)  
# ## Outstanding
# - background/no-data/masking area outside continent
# - cf-compliance? (nice to have)  
# 
# ## Conceptual
# ~~- Missing data inland (perhaps over wate bodies).~~  
# ~~- How to distinguish between inland missing and area outside continent (does not matter perhaps)~~
# 
# ## Outputs
# 1. Monthly nc files
# 2. Water Year nc files
# 3. Water Year all variables combined into one file

# ### 0: Snowf_tavg  
# units : kg m-2 s-1  
# standard_name : snowfall_rate  
# long_name : snowfall rate  
# 
# ~~### 1 : Rainf_tavg  
# units : kg m-2 s-1  
# standard_name : precipitation_rate  
# long_name : precipitation rate~~  
# 
# ### 2: SWE_tavg  :: for daily, statistics = max
# units : kg m-2  
# standard_name : liquid_water_content_of_surface_snow  
# long_name : snowfall rate  
# 
# ### 3: Tair_f_tavg  : for daily, statistics = mean  
# units :  K  
# standard_name : air_temperature  
# long_name : air temperature
# 
# ### 4: Qg_tavg  : for daily, statistics = sum?  should be mean  [july 29th]  
# units : W m-2  
# standard_name :  downward_heat_flux_in_soil  
# long_name :  soil heat flux  
# 

# In[1]:


# This cell is parameters passed from outside
cores = 10  # hardcoded to -1 below, so this has no effect currently


# In[2]:


import matplotlib.pyplot as plt

import platform
import os
import time
from datetime import timedelta

import numpy as np
import pandas as pd
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
import xarray as xr
import netCDF4

import hvplot.pandas
import hvplot.xarray

import logging
import argparse
from joblib import Parallel, delayed

# node is called host_machine in Julia
# ss, node, = platform.uname()[:2]
ss = platform.system()
node = platform.node()


# In[3]:


if 'STAFF-BY-M01' in node and "Windows" in ss:
    root_dir = "C:"
    base_folder = f"{root_dir}/Github/coressd/Blender/NoahMP"
    in_folder = f"{base_folder}/3hourly"
    out_folder = f"{base_folder}/temp"
elif 'STAFF-BY-M01' in node and "Linux" in ss:
    root_dir = "/mnt/c"  #"/home/yadav.111"
    base_folder = f"{root_dir}/Github/coressd/Blender/NoahMP"
    in_folder = f"{base_folder}/3hourly"
    out_folder = f"{base_folder}/temp"
elif "borg" in node and "Linux" in ss:
    # nasa-discover nodes
    base_folder = "/discover/nobackup/projects/SEUP/mwrzesie/KEEP_fromRhaeSung/SEUP/NoahMP/MERRA2/OL/OUTPUT/SURFACEMODEL"
    in_folder = base_folder
    out_folder = "/discover/nobackup/projects/coressd/Blender/NoahMP"
else:
    print("Unknow system. in_folder and out_folder cannot be set")
    assert(False)
print(f"in_folder  : {in_folder}")
print(f"out_folder : {out_folder}")


# # 1. Generate Daily Data
# By concatenating 3-hourly datasets  

# In[4]:


# yr_mon_list = os.listdir(f"{base_folder}/3hourly")
# or, Explicitly specify year month list to process
yr_mon_list = ['201510', '201511', '201512', '201601', '201602', '201603', '201604', '201605', '201606', '201607', '201608', '201609']
variables = ["Snowf_tavg", "SWE_tavg", "Tair_f_tavg", "Qg_tavg"]  #"Rainf_tavg", 


# In[5]:


def combine_seup(in_folder, month_subfolder, out_folder):
    """ Combine 3-hourly SEUP data to daily
    Need to pass these variables
    combined_folder = f'{out_folder}/combined/{var}'
    combined_nc_file = f'{combined_folder}/{month_subfolder}.nc' #appended monthly to indicate this is monthly mean data
    nc_files
    f"{folder}/{nc_file}"

    """
    folder = f"{in_folder}/{month_subfolder}"
    nc_files = [f for f in os.listdir(folder) if f.endswith(".nc")]
    # Sort in ascending order of data
    nc_files.sort(key = lambda x: pd.to_datetime(x.split("_")[2].split(".")[0]))
    for var in variables:
        print(f"\tVariable: {var}")
        dt_list = []
        # We will save the combined netcdf file here
        combined_folder = f'{out_folder}/combined/{var}'
        # if not os.path.exists(combined_folder):  # creates problem with joblib parallization, perhaps due to two threads tyring to create same directory
        os.makedirs(combined_folder, exist_ok=True)  # TODO: better to create this folders outside of this function to avoid joblib error
        combined_nc_file = f'{combined_folder}/{month_subfolder}.nc' #appended monthly to indicate this is monthly mean data
        if not os.path.exists(combined_nc_file):
            # Combine 3-hourly data for one varialbe into 1 file
            count = 0
            for nc_file in nc_files:
                # print(f"nc_file : {nc_file}")
                dt_list.append(pd.to_datetime(nc_file.split("_")[2].split(".")[0]))
                if count == 0:
                    ds = xr.open_dataset(f"{folder}/{nc_file}", engine='netcdf4')
                    da = ds[var]
                else:
                    ds = xr.open_dataset(f"{folder}/{nc_file}", engine='netcdf4')
                    da_temp = ds[var]
                    da = xr.concat([da, da_temp], dim='time')  #, dim='time'
                count += 1

            # Get (Assign) lat/lon
            # ds["X"] = np.int64(ds.X)
            # missing_value = ds.missing_value
            ll_lon = ds.attrs["SOUTH_WEST_CORNER_LON"]
            ll_lat = ds.attrs["SOUTH_WEST_CORNER_LAT"]
            # Hardcode DX and DY because extracting this from attrs adds some decimals (probably float related or real not sure)
            DX = 0.05 # ds.attrs["DX"]  #grid spacing in x-direction; TODO: check why this is not exactly 0.05 but rather DX=0.05000000074505806; same issue for DY
            DY = 0.05 #ds.attrs["DY"]  #grid spacing in y-direction; 
            print(f"ll_lon: {ll_lon} , ll_lat: {ll_lat}, DX: {DX}, DY: {DY}")
            # Construct lon and lat coordiantes
            lon = ds.east_west*DX + ll_lon
            lat = ds.north_south*DY + ll_lat
            # Coordinates can also be set or removed by using the dictionary like syntax:
            da["x"] = lon
            da["y"] = lat
            da["time"] = dt_list ## better?
            da = da.rename({"east_west":"x", "north_south":"y"})  #dims and coordinates name should match; new warning as of 10/25/2022
            print("\nResampling to Daily\n")
            if var in ["Snowf_tavg", "Rainf_tavg"]:
                print(var)
                # Converr the accumalation rate per second to per hour; then how much it accumulated in 3 hours
                da = da*3600*3  #units : kg m-2 s-1   to kg m-2 hr-1 convert rate to accumulation per second in 3 hours;ie, 3600 seconds x 3 hours
                daily = da.resample(time="1D").sum(skipna=False)  # or 24H; retains actual datetime
                daily.data = daily.data/1000.0  # convert from mm to meter
            elif var == "SWE_tavg":
                print(var)
                daily = da.resample(time="1D").max(skipna=False)  # or 24H; retains actual datetime 
                daily.data = daily.data/1000.0  # convert from mm to meter
            elif var == "Tair_f_tavg":
                print(var)
                daily = da.resample(time="1D").mean(skipna=False)

            elif var == "Qg_tavg":
                print(var)
                daily = da.resample(time="1D").mean(skipna=False)
            else:
                print("Error: Unknown Varialbe. Stopping processing")
                assert(False)
            # print(np.array_equal(daily.data, daily_g.data))  # at least data are equal when using resample and groupby
            # # Conver to Dataset, so this is proper nc file
            daily = xr.Dataset({var:daily})
            daily.to_netcdf(combined_nc_file)  #260 MB size


# In[6]:


# To procecess sequentially (in series)
# for month_subfolder in yr_mon_list[1:3]:
#     combine_seup(in_folder, month_subfolder, out_folder)
# Or in parallel 
# (replace -1 by cores, but -1 seems better option here as we don't have to worry about explicity passing cores; will be decided on by nunber of cores asked in slurm scheduler
Parallel(n_jobs=-1)(delayed(combine_seup) (in_folder, month_subfolder, out_folder) for month_subfolder in yr_mon_list[5:7])

# release memory; should work whether these variables exist or not
ds = None
da = None
da_temp = None
daily = None  


# # 2. Concatenate monthly files to get water year file
# Independent of above code
# Here, each variable is in sepate folder/file

# In[7]:


# Find way to manually generate this year_month file pre-fixes and in ascending order; 
# that way we do not need to depend on the way the monthly files are staged
# os.listdir(f"{out_folder}/combined/Snowf_tavg")

# Using list comprehension takes way too long time perhaps due to memory issues
ds_wy = xr.concat([xr.open_dataset(f"{base_folder}/combined/{var}/{f}", chunks=True) for f in monthly_nc_files], dim="time")  
# In[ ]:


year = "2016"  # water year 2016 is from Oct 2015 to Sept 2016 
for var in variables:
    print(var)
    # Get monthly input files
    monthly_nc_files = os.listdir(f"{out_folder}/combined/{var}")  # probably need to generate the year month string from 201510 to 201609
    wy_folder = f'{out_folder}/WY/{var}'
    # if not os.path.exists(wy_folder):
    os.makedirs(wy_folder, exist_ok=True)
    combined_nc_file = f'{wy_folder}/{year}.nc' #appended monthly to indicate this is monthly mean data
    # if not os.path.exists(combined_nc_file):  # need to create this regardless if file in combined folder is changed
    f = monthly_nc_files[0]
    wy_ds = xr.open_dataset(f"{out_folder}/combined/{var}/{f}", engine='netcdf4')
    for f in monthly_nc_files[1:]:
        ds_temp = xr.open_dataset(f"{out_folder}/combined/{var}/{f}", engine='netcdf4')  #, chunks={'time': 1}
        # da_temp = ds[var]
        wy_ds = xr.concat([wy_ds, ds_temp], dim='time')
    wy_ds.to_netcdf(combined_nc_file)  # ~ 3 GB size
wy_ds = None  # To prevent: Permission denied: b'C:\\Github\\coressd\\Blender\\NoahMP\\temp\\WY\\Snowf_tavg\\2016.nc'


# # 3. Combine individual WY nc files files into one file with all varabiles
# This is even bigger ~12 GB file
# Combine variables from different folder (STEP 2) into single file with all four variables for that water year  
# But may not be necessary ...  
# variables = ['Snowf_tavg', 'SWE_tavg', 'Tair_f_tavg', 'Qg_tavg']

# In[ ]:


year = "2016"  # water year 2016 is from Oct 2015 to Sept 2016 
# ds_merged = xr.merge([xr.open_dataset(f"{base_folder}/combined/{var}/{f}", chunks=True) for f in monthly_nc_files])  #chunks={'time': 1} # use only for different variables
wy_ds_merged = xr.merge([xr.open_dataset(f"{out_folder}/WY/{v}/{year}.nc") for v in variables])
# Here merged implies merged on variables (in Xarray paralay)
merged_wy_folder = f'{out_folder}/WY_merged'
# if not os.path.exists(merged_wy_folder):
os.makedirs(merged_wy_folder, exist_ok=True)
merged_nc_file = f'{merged_wy_folder}/{year}.nc' #appended monthly to indicate this is monthly mean data
# if not os.path.exists(merged_nc_file):
wy_ds_merged.to_netcdf(merged_nc_file)  #12.6 GB size
    


# # END
# End of processing SEUP data for whole of north America for one Water Year  
# 
# ---  
