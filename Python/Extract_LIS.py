""" Processes Land Information System (LIS) data.

USAGE: 
    >>> python /discover/nobackup/projects/coressd/Github/verse/Python/Extract_LIS.py 2016
    >>> python /discover/nobackup/projects/coressd/Github/verse/Python/Extract_LIS.py 2016 --RES 010

Currently it only works for LIS data with 0.01 degree resolution but the script can
be extended to other resoulutions as required in future.

This function takes in raw LIS data, performs necessary cleaning and transformations,
and returns a processed dataset ready for analysis. The processing steps may include
handling missing values, normalizing data, and extracting relevant features.
Outputs of this script will be used to run Blender Julia optimization routine.

Args:
    water_year, [RES].

Inputs:
    LIS model run outputs.
Returns: Saves netcdf files
    - Snowf_tavg
    - SWE_tavg 
    - Tair_f_tavg 
    - Qg_tavg
These 4 netcdf files along with scf (MODIS SCF) will be used for Julia Blender runs.

Notes:
    Hardcoded values for grid spacing (DX, DY).
    DX = 0.01  # grid spacing in x-direction
    DY = 0.01  # grid spacing in y-direction  
    chunk = {"x": 2**8, "y": 2**8}
    Scaling applied to Snowf_tavg, SWE_tavg, Tair_f_tavg variables.  

Updates:
    Nov 28, 2024 : Exported from Jupyter/Papermill workflow. Moving forward, this script will be used to extract LIS data for Blender optimization routine.

"""

# import sys
import platform
import os
import numpy as np
import pandas as pd
import xarray as xr
import logging
import argparse
# from dask.distributed import Client, LocalCluster
# from joblib import Parallel, delayed
# sys.path.append("/discover/nobackup/projects/coressd/Github/verse/Python")
from tictoc import tic, toc

# node is called host_machine in Julia
# ss, node, = platform.uname()[:2]
ss = platform.system()
node = platform.node()

parser = argparse.ArgumentParser(description='Extract LIS data.')
parser.add_argument('water_year', help='Water Year', type=str)
parser.add_argument('--RES', help='Grid size (Resolution) of SEUP Inputs', type=str, default="010")  # SEUP was 050 but may not work here.
# parser.add_argument('--cores', help='Number of cores for multiprocessing', type=int, default=12)  # don' use -1 on DISCOVER  
args = parser.parse_args()
water_year = args.water_year
RES = args.RES  # "050"  # Grid Resolution: 050 = 0.050 degree (5km); 025 = 0.025 degree; 010 = 0.01 degree; 001 = 0.001 degree (100 meters) etc
logging.basicConfig(filename=f"out/extract_lis_{water_year}_{RES}.log", level=logging.INFO, format='%(asctime)s : %(message)s')

# Hardcoded: Path to LIS data and Output folder. Will only work on Discover.
lis_folder = "/discover/nobackup/projects/coressd/mwrzesie/runOL/out.OL/SURFACEMODEL"  # new 1 km data
out_folder = f"/discover/nobackup/projects/coressd/Blender/Inputs_{RES}"  # output of this script will be saved in this parent folder.

variables = ["Snowf_tavg", "SWE_tavg", "Tair_f_tavg", "Qg_tavg"]  # LIS variables to extract for Blender optimization routine
# Generate list of water year and months. May seem a bit convoluted because of the way water year is defined.
months = ["10", "11", "12", "01", "02", "03", "04", "05", "06", "07", "08", "09"]
yr_mon_list1 = [f"{int(water_year) - 1}{mon}" for mon in months[:3]]  # include months "10", "11", "12"
yr_mon_list2 = [f"{water_year}{mon}" for mon in months[3:]]  # remaining months "01", "02", "03", "04", "05", "06", "07", "08", "09"
yr_mon_list = yr_mon_list1 + yr_mon_list2

# Take a sample netcdf file to extract a-priori data and attrs
# Assumption: All files follow the same structure. Final output may be all crap if this assumption does not hold
month_subfolder = yr_mon_list[6]
print(month_subfolder)
folder = f"{lis_folder}/{month_subfolder}"
nc_files = [f for f in os.listdir(folder) if f.endswith(".nc")]
nc_files = [f for f in nc_files if f.startswith("LIS_HIST_")]  # remove this extra LIS output generated Melissa.  
# # Sort in ascending order of data
# nc_files.sort(key = lambda x: pd.to_datetime(x.split("_")[2].split(".")[0]))
nc_file = nc_files[0]
# Read one sample file to extract attrs; cooridinates; and do manipulation as required
ds = xr.open_dataset(f"{folder}/{nc_file}", engine='netcdf4')

# Remove the Extra variables from LIS Run that we do not need for Blender optimization routine
# ds[variables]  # not good because it seems to drop the time dimension
# extra_variables = [v for v in list(ds.data_vars) if not v in variables]  # this also remove lon, lat variables which need for now. 
extra_variables = ['AvgSurfT_tavg', 'Albedo_tavg', 'SnowDepth_tavg', 'Snowcover_tavg', 'LAI_tavg', 'Greenness_inst', 'TotalPrecip_tavg']
print(f"Drop these variables from the dataset: \n {extra_variables}")
ds = ds.drop_vars(extra_variables)

# Manually genenate lan/lon values. Use this to replace unevenly spaced lat/lon from LIS runs. These likely persisited dut to floating point conversion.
ll_lon = ds.attrs["SOUTH_WEST_CORNER_LON"]
ll_lat = ds.attrs["SOUTH_WEST_CORNER_LAT"]
# Hardcode DX and DY because extracting this from attrs adds some decimals (probably float related or real not sure)
DX = 0.01  # ds.attrs["DX"]  #0.05 # grid spacing in x-direction; TODO: check why this is not exactly 0.05 but rather DX=0.05000000074505806; same issue for DY
DY = 0.01  # ds.attrs["DY"]  #0.05 # grid spacing in y-direction; 
print(f"ll_lon: {ll_lon} , ll_lat: {ll_lat}, DX: {DX}, DY: {DY}")

lon = ds.east_west.data * DX + ll_lon
lat = ds.north_south.data * DY + ll_lat
lon = np.around(lon, 3)
lat = np.around(lat, 3)

chunk = {"x": 2**8, "y": 2**8}
logging.info(f"Chunks: {chunk}")

dt_list = []
count = 0
tic()
for month_subfolder in yr_mon_list:
    print(month_subfolder)
    folder = f"{lis_folder}/{month_subfolder}"
    nc_files = [f for f in os.listdir(folder) if f.endswith(".nc")]
    nc_files = [f for f in nc_files if f.startswith("LIS_HIST_")]  # new files by Melissa has one new file. remove it.
    # Sort in ascending order of data
    nc_files.sort(key=lambda x: pd.to_datetime(x.split("_")[2].split(".")[0]))
    for nc_file in nc_files:
        # print(f"nc_file : {nc_file}")
        dt_list.append(pd.to_datetime(nc_file.split("_")[2].split(".")[0]))
        if count == 0:
            ds = xr.open_dataset(f"{folder}/{nc_file}", engine='netcdf4', chunks={"lon": 2**8, "lat": 2**8})
            attrs = ds.attrs  # we'll attach it later to the final processed data
            ds = ds.drop_vars(extra_variables)
            ds = ds.set_coords(["lon", "lat"])
        else:
            ds_temp = xr.open_dataset(f"{folder}/{nc_file}", engine='netcdf4', chunks={"lon": 2**8, "lat": 2**8})
            ds_temp = ds_temp.drop_vars(extra_variables)
            ds_temp = ds_temp.set_coords(["lon", "lat"])
            ds = xr.concat([ds, ds_temp], dim='time')  #, dim='time'
        count += 1
    toc()
tic()
# Replace the old data with new ones I created
ds["lon"].data = np.around(lon, 3)
ds["lat"].data = np.around(lat, 3)

# First rename lon lat to x, y respectively
ds = ds.rename({"lon": "x", "lat": "y"})
# Swap dims to x and y
# ds.swap_dims({"lon":"x", "lat":"y"})
# ds = ds.swap_dims({"east_west":"lon", "north_south":"lat"})
ds = ds.swap_dims({"east_west": "x", "north_south": "y"})
logging.info(f"File size before slice: {ds.nbytes / 1e9} GB")
# Truncate lower latitude pixels and match SEUP extent
# ds.sel(lat=slice(24.875, 72))  #match SEUP extent
ds = ds.sel(y=slice(30, 72))  # max value of latitude
# Update global attrs for changed latitude
ds.attrs["SOUTH_WEST_CORNER_LAT"] = ds.y.min().data.item(0)  # new we changes lon/lat at x/y
logging.info(f"File size after latitude slice (30, 72): {ds.nbytes / 1e9} GB")
logging.info("====================================================")

combined_folder = f"{out_folder}/lis/WY{water_year}"
os.makedirs(combined_folder, exist_ok=True)

int_variables = ["Snowf_tavg", "SWE_tavg", "Tair_f_tavg"]  # , "Qg_tavg". Use only the variables that will be saved to uint16 
# 1. Define nodata fill value that will be used for all 3 variables to be saved in uint32. 
d_type = np.uint16
nodata_value = np.iinfo(d_type).max  # 65535 for uint16
for var in int_variables:
    logging.info(f"Saving {var}")
    if var == "Snowf_tavg":
        print(f"Convert units for {var}")
        # 1. Convert unit of snowfall
        # kg m-2 s-1 to kg m-2 hr-1 convert rate to accumulation per second in 24 hours;ie, 3600 seconds x 24 hours
        ds1 = ds[var]  # .to_dataset() # keep as dataarray for calculations
        ds1.data = np.round(ds1.data * 3600 * 24)  # only for Snowf_tavg, convert to mm/day. divide by 1000 to get meters in Blender run               
    elif var == "SWE_tavg":
        # ds1 = ds[var].to_dataset() 
        # # ds1[var].data = ds[var].data/1000.0  # convert from mm to meter; prior error without chunks: # MemoryError: Unable to allocate 13.7 GiB for an array with shape (67, 4700, 11700) and data type float32
        ds1 = ds[var]  # 1kg/m2 = 1mm. So, in Blender run, divide by 1000 to get meters
        ds1.data = np.round(ds1.data)
    elif var == "Tair_f_tavg":
        ds1 = ds[var]  # da.to_dataset(name='values') will give name called value;  or xr.Dataset({var:da2})        
        ds1.data = np.round(ds1.data * 100)  # for do calculation; for Air temperate. divide by 100 to get kelvin with decimals in Blender run        
    # Re-define fillvalues and save dataset
    ds1 = ds1.fillna(nodata_value)
    ds1.data = ds1.data.astype(d_type)
    ds1 = xr.Dataset({var: ds1})
    ds1[var].attrs["_FillValue"] = nodata_value  # update fill value from nan to 255; do after creating dataset so it is nested with variable.              
    ds1 = ds1.chunk(chunks=chunk)
    ds1 = ds1.compute()  # ~10 minutes (~210 GB memory). compute loads in memory (Numpy array); load (akin to loading); persist: still in dask array. But xarray will load whenever required.  
    logging.info(f"{var}: File size (uncompressed): {ds1.nbytes / 1e9} GB")
    combined_file = f'{combined_folder}/{var}'
    # To save netcdf file. But takes forever to save!
    # encoding = {var: {'zlib': True} for var in ds.data_vars}
    encoding = {var: {'zlib': True}}  # 'complevel':5 
    ds1.to_netcdf(f"{combined_file}.nc", encoding=encoding)  # ~20 minutes with original floats.
    # ds1.to_netcdf(f"{combined_file}_uncompress.nc")
    logging.info(f"\tDone Saving: {var}")
    logging.info("---------------------------------------------------")

var = "Qg_tavg"  # Qg_tavg is the only variable that will be saved in float. This is used as template, mask, missing_val etc. in all of Julia workflow.
# ds1 = ds[var] #.to_dataset()
ds1 = ds[var].to_dataset()  # da.to_dataset(name='values') will give name called value;  or xr.Dataset({var:da2})        
ds1 = ds1.chunk(chunks=chunk)
ds1 = ds1.compute()  # ~10 minutes (~210 GB memory). compute loads in memory (Numpy array); load (akin to loading); persist: still in dask array. But xarray will load whenever required.  
logging.info(f"{var}: File size (uncompressed): {ds1.nbytes / 1e9} GB")
combined_file = f'{combined_folder}/{var}'
# To save netcdf file. But takes forever to save!
# encoding = {var: {'zlib': True} for var in ds.data_vars}
encoding = {var: {'zlib': True}}  # 'complevel':5 
ds1.to_netcdf(f"{combined_file}.nc", encoding=encoding)  # ~20 minutes.
# ds1.to_netcdf(f"{combined_file}_uncompress.nc")
logging.info(f"\tDone Saving: {var} \n")
logging.info(f"Output saved here: {combined_file}")
logging.info(f"Finished:  {toc()}")
