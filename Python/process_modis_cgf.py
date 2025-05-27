#!usr/bin/env python

""" Extract MODIS CGF over North America matching the LIS/SEUP Data 

    Usage
    ====== 
    python ../verse/Python/process_modis_cgf.py 2016 --cores 8  
    sbatch b_process_modis_cgf.sh $water_year (Slurm script for Discover)

    Description
    ===========
    This script extracts MODIS CGF data for North America matching the LIS data resolution (1 km).
    It uses MOD44B tree cover fraction data to correct for forest cover fraction.
    
    STEPS
    -----
    1. Use sample templates for "SWE_tavg" for matching MOD10A1F to LIS resolution
    2. Generate year_doy_list in MODIS naming convention
    3. Read MOD44B and divide by 100 to convert 0 to 1 range
    4. Merge daily adjascent MOD10A1F tiles (40 for North America)
        - convert to float32
        - replace 211 flag (Nights) with nans
        - convert from NDSI to Snow cover fraction
        - Correct for cover fraction
    5. Reproject and match to LIS template
    6. Use nan mask from LIS and replace in MODIS pixels  
    7. Combine daily data to produce annual data
    8. Save the final output "SCF.nc" inside the same folder with LIS input files  
    After this, data in ready for Blender run

    Compute Resources
    =================
    Milan node with 32 cores and 300 GB RAM. 
    More than this causes NodeFailError, likely due to memory limit.   
    Runtime ~1 hour for 1 year with 32 cores with parallelization.

    Input
    =====
    Raw downloaded files in original sinusoidal projection
    modis_download_folder   = /discover/nobackup/projects/coressd/OSU/MOD10A1F.061/MODIS_Proc/download_snow     # /2016001/001
    modis_mosaic_template   = /discover/nobackup/projects/coressd/Blender/Template/MOD10A1F_NA_mosaic.tif  # used only to get bounds for mosaicking 40 tiles
    lis_template            = /discover/nobackup/projects/coressd/Blender/Inputs/WY2016/SWE_tavg.nc
    TreeCoverFraction       = /discover/nobackup/projects/coressd/Blender/Modis/MOD44B/Reproj/MOD44B.A{water_year}065.061.nc

    Outputs 
    =======
    Intermediate Folders: 
    - /discover/nobackup/projects/coressd/Blender/Modis/MOD10A1F  
        - mosaic2016: temporary mosaics only and consider deleting (or not saving)
            - MOD10A1F.061_2010261.tif # notice tif file retained instead of nc. Is a daily North America mosaic in original Sinusoidal projection
        - clipped2016_010: clipped files
            - MOD10A1F.061_2010265.nc # clipped and reprojected to LIS resolution and extent. 15 GB.
    Final Output: 
    - /discover/nobackup/projects/coressd/Blender/Inputs/WY2016/SCF.nc (~12 GB)  
    - Final water year concatenated MODIS CGF data for Blender run.

    TODO
    =====
    - Create a different layer that holds info on QA/QC related to the data and processing
    - Can use bitmask layer to hold this information
    - Do AND operation on the bitmask layer for all the pixels that are covered by LIS data
    - Idea is to use this layer in Blender script to make necessary corrections
    - For example, if SCF is zero, and bitmask 3 is set, then we give some reasonable value to SCF
    - Might be much harder to implement on a short timeline  
    - There are systematic issues with MOD10A1F data 
    - Thus these actions points can also be ignored and just convert all flags to zero except polar nights (as is already implemented)

    History
    =======
    June 01, 2024:
    - Tested saving the zarr files. Though zarr is faster, it is not currently usable we NetCDF files for Blender (Julia) run. 
    - Switching back to saving NetCDF files. 
    June 12, 2024: 
    - Saving final ouput as uint8 (0 - 100 and _FillValue=255)and renamed as SCF.nc => updated call_Blender_v18 with this name as well
    - Will require settig 255 to nodata/missing and dividing by 100 to get fraction between 0 and 1..
    - saving intermediate clipped modis files as np.float32. Reduced folder size from 14 GB to 6.2 GB.
    June 21, 2024: 
    - Revert back to saving final ouput as float32 but with rounding and values between (0 - 100)
    - rouding saves space
    - divide by 100 to get fraction between 0 and 1 in Blender run.  
    - Error for unit8 was caused in Estimate_v59.jl, likely due to mix of missing due to actual nodata and polar-nights nodata.  
    May 27, 2025:
    - removed  and f.startswith("MOD10A1F") so that MYD can be used to replace missing MODIS data    
    - Correction of Tree cover fraction done at higer resolution (~500 m MOD10A1F) instead of older LIS resolution (1 km)
    - Replace default nearest-neighbor resampling with average resampling. This may help reduce isolated zero SCF pixels.
    - Added parallelization using joblib. Reduced runtime from ~6 hours to 1 to 2.5 hours without increasing compute cost.

    Author
    ======
    Bidhya N Yadav
    The Ohio State University
    2022-2025

"""
import os
import datetime
import time
import pandas as pd
from tictoc import tic, toc
import numpy as np
import xarray as xr
import rioxarray
from rioxarray.merge import merge_arrays
from rasterio.enums import Resampling  # for average, cubic, bilinear etc. resampling during reprojection
import logging
import argparse
from joblib import Parallel, delayed

parser = argparse.ArgumentParser(description='Extract MODIS CGF files for North America.')
parser.add_argument('water_year', help='Water Year for processing', type=str)
parser.add_argument('--cores', help='Number of cores for parallelization', type=int, default=2)  # don' use -1 on DISCOVER MILAN node  
args = parser.parse_args()
water_year = args.water_year
cores = args.cores
chunk = {"x": 2**8, "y": 2**8}

# Create "OUT" folder for Slurm output and Python logging
os.makedirs("OUT", exist_ok=True)  # must be created manually for Slurm output
logging.basicConfig(filename=f"OUT/process_modis_cgf_{water_year}.log", level=logging.INFO, format='%(asctime)s : %(message)s')
logging.info('  ')
logging.info('-------------------------START LOGGING--------------------------------------')
logging.info(f'water_year={water_year}     \ncores={cores}')
logging.info(f'chunks={chunk}')

# List of tiles that cover North America (our AOI)
modis_tile_list = [
    'h07v03', 'h07v05', 'h07v06', 'h08v03', 'h08v04', 'h08v05', 'h08v06', 'h09v02', 'h09v03', 'h09v04', 'h09v05', 'h09v06', 'h10v02', 'h10v03', 
    'h10v04', 'h10v05', 'h10v06', 'h11v02', 'h11v03', 'h11v04', 'h11v05', 'h11v06', 'h12v01', 'h12v02', 'h12v03', 'h12v04', 'h12v05', 'h13v01', 
    'h13v02', 'h13v03', 'h13v04', 'h14v01', 'h14v02', 'h14v03', 'h14v04', 'h15v01', 'h15v02', 'h15v03', 'h16v01', 'h16v02']
logging.info(f"Number modis tiles required: {len(modis_tile_list)}")


def doy_to_date(modis_date, reverse=False):
    """Convert modis year and doy (day or year) to standard python datetime
        date of year should be contexualized vis-a-vis year
        because day of year in all years will not be same (leap year for example)

        modis_date [string] : Modis native date that include year and year of day
        reverse: NOT CURRENTLY USED. Placeholder for converting regular date to modis date  
    """
    begin_year = int(modis_date[:4])
    begin_days = int(modis_date[4:])
    begin_ymd = datetime.datetime(begin_year, 1, 1) + datetime.timedelta(begin_days - 1)
    return begin_ymd


coressd_folder = "/discover/nobackup/projects/coressd"
modis_download_folder = f"{coressd_folder}/OSU/MOD10A1F.061/MODIS_Proc/download_snow"  # /2016001/001
# Manually generate a date range for water year and convert to DOY format of MODIS naming convention
start_date = datetime.datetime.strptime(f"{int(water_year) - 1}-10-01", "%Y-%m-%d")  # 
end_date = datetime.datetime.strptime(f"{water_year}-09-30", "%Y-%m-%d")  # 
date_generated = pd.date_range(start_date, end_date)
# convert to DOY format of MODIS naming convention
year_doy_list = list(date_generated.strftime("%Y%j"))  # Does this match Day of year for MODIS nameing convention?  
# year_doy_list = ["2016001", "2016002", "2016003"]
# year_doy_list = os.listdir(modis_download_folder)

# Global variablees
# Get a template raster for reprojection match. This is the SEUP data. Any nc file from ../coressd/Blender/Inputs_050/combined/SWE_tavg will work.
lis_folder = f"{coressd_folder}/Blender/Inputs"  # _{RES} NoahMP
lis_template = xr.open_dataset(f'{lis_folder}/WY{water_year}/SWE_tavg.nc', chunks="auto")  
lis_template = lis_template.isel(time=0)  # Just one day of data sufficient. Use this template for reprojection match of MODIS data
lis_template.rio.write_crs("EPSG:4326", inplace=True)  # TODO : check if this should Plate Caree  projection  
nan_mask = np.isnan(lis_template["SWE_tavg"].data)  # we'll use this nan mask to replace MODIS data so they match exactly with noahmp
DATAFIELD_NAME = "CGF_NDSI_Snow_Cover"  # we need this variable from MOD10AF

# Create folder to save intermediate modis files. 
modis_folder = f"{coressd_folder}/Blender/Modis"  
# mosaic_folder = f"{modis_folder}/MOD10A1F/mosaic{water_year}"  # uncomment if Saving mosaics for QA/QC and testing only. Not required for final run.  
# os.makedirs(mosaic_folder, exist_ok=True)
clip_folder = f"{modis_folder}/MOD10A1F/clipped{water_year}"  # intermediate daily Geotiff files saved here.
os.makedirs(clip_folder, exist_ok=True)

# Get Tree cover fraction data. This new MOD44B data is already matched to MOD10A1F Snow Cover data.
daF = rioxarray.open_rasterio(f'{coressd_folder}/Blender/Modis/MOD44B/Reproj/MOD44B.A{water_year}065.061.tif').squeeze()
logging.info(f"daF Shape = {daF.shape}")
daF.data = daF.data.astype(np.float32)  # default might be 64 bit, so being explicit here.
daF.data = daF.data / 100

logging.info("Start Reproj match of Modis (MOD10A1F).")
# New May 16, 2025: Use the a-priori saved template raster to get the bounds when mosaicking the 40 tiles.
MOD10A1F_mosaic_template = rioxarray.open_rasterio(f"{coressd_folder}/Blender/Template/MOD10A1F_NA_mosaic.tif")
bounds_tuple = MOD10A1F_mosaic_template.rio.bounds()
logging.info(f"MOD10A1F_mosaic_template Bounds tuple = {bounds_tuple}")
del MOD10A1F_mosaic_template


def extract_modis(download_folder, daF):
    """ Extract the MODIS CGF matching the SEUP data ready for Blender Julia run
    file_path
    =========
    Full path to netcdf file that was processed at NASA discover
    not tested: hdf4 file not tested but should work seamlessly
    """
    files = [f for f in os.listdir(download_folder) if f.endswith(".hdf")]  # removed  and f.startswith("MOD10A1F") so that MYD can be used to replace missing MODIS data
    # Get the subset of modis 10 degree files for this day
    subset = [f for f in files if f.split(".")[2] in modis_tile_list]
    hdf_filename = subset[0]  # TODO: can be problematic if first file is the replacement AQUA (MYD) file
    product, year_doy, tile, version, _, _ = hdf_filename.split(".")  # Create filenme to one matching the pattern processed by Sarith
    logging.info(f"\t For {year_doy}, total number of tiles = {len(files)} and subset file count = {len(subset)}")

    out_tif_name = f"{product}.{version}_{year_doy[1:]}"
    # if not os.path.exists(f"{clip_folder}/{out_tif_name}.nc"):
    da = merge_arrays([rioxarray.open_rasterio(f'{download_folder}/{tif}')[DATAFIELD_NAME] for tif in subset], bounds=bounds_tuple)  # chunks={"x": 2**8, "y": 2**8}
    da = da.squeeze()
    # da.rio.to_raster(f"{mosaic_folder}/{out_tif_name}.tif", driver="COG")  # saving temporarily only to check if everything is working.

    # Convert to float and replace various flags to appropriate values
    da.data = da.data.astype(np.float32)  # float
    # da.data[da.data >100] = np.nan  #this causes most of the pixels to be nan. Hence, cannot be used in Blender
    # da.data[da.data == da._FillValue] = np.nan
    da.data[da.data == 211] = np.nan  # Night 
    # da.data[da.data == 239] = np.nan  # Ocean . This nan seems to problem around edges because LIS will have valid data but MODIS may have some missing data.  
    da.data[da.data > 100] = 0  # Give rest of the flags value of zero snow! TODO: change to nan because blender can handle missing now  
    da.data = da.data / 100  # to get zero to one range. This is required for calculations in next steps.
    # a) May 16, 2022: Convert NDSI to snow cover fraction (FRA); FRA = 0.06 + 1.21 * NDSI. (Salomonson and Appel, 2004)
    da.data = 0.06 + 1.21 * da.data
    # Correct of biases introduced by above equation
    da.data[da.data <= 0.06] = 0
    da.data[da.data > 1] = 1
    # Apply correction for Forest Tree Cover Fraction. Memory error for 1 km run, hence, using Milan node on Discover.
    # da.data = da.data / (1 - (daF + np.finfo(np.float32).eps))  # 2.220446049250313e-16 for float64; 1.1920929e-07 for float32
    # da.data = da.data / (1 - (daF - np.finfo(np.float32).eps))
    da.data = da.data / (1 - daF)  # no need for eps because none of daF values is greater than 100 (checked manually) 
    da.data[da.data > 1] = 1  # to correct any bias introduced by above equation

    # Match to LIS template using reproj match   
    da_clipped = da.rio.reproject_match(lis_template, resampling=Resampling.average)  # can have memory error    
    da_clipped.data[da_clipped.data == da_clipped._FillValue] = np.nan  # because reproj match introduced fill value of 255
    da_clipped.data[da_clipped.data > 1] = 1
    # TODO: We will get some 255 values for nodata; check if other values also included (due to interpolation)
    da_clipped.data[nan_mask] = np.nan  # done again, because reproj match will fill new nodata with 255 (perhaps this value from attrs)
    # da_clipped.data = da_clipped.data / 100
    # OPTIONAL: Fix attrs, no-data etc here so that data is self-describing
    # Conver to Dataset, so this is proper nc file
    ds_clipped = xr.Dataset({DATAFIELD_NAME: da_clipped})
    encoding = {DATAFIELD_NAME: {'zlib': True}}
    ds_clipped.to_netcdf(f"{clip_folder}/{out_tif_name}.nc", encoding=encoding)


# Make of List of downloaded sub-folder. To make it easy to parallelize using Joblib
download_folder_list = []
for year_doy in year_doy_list:
    year = year_doy[:4]  # "2016"
    doy = year_doy[4:]  # "215"
    download_folder = f"{modis_download_folder}/{year}{doy}/{doy}"
    download_folder_list.append(download_folder)
    # extract_modis(download_folder)

# ------------------------------------------------------------------------------------------------
# # This is no more required as we are using the bounds from MOD10A1F_mosaic_template that was created as part of Tree Cover Fraction processing
# # Lower left (h07v06) and Upper Right (h16v01) MODIS Granule for North America
# # Use bounds from this two files in mosaicking. esle there is one extra pixel in x, y or both direction.
# # This should help constrain and we always get same size of output raster
# download_folder = download_folder_list[182]  # April 01, 2015
# # DATAFIELD_NAME = "CGF_NDSI_Snow_Cover"
# files = [f for f in os.listdir(download_folder) if f.endswith(".hdf") and f.startswith("MOD10A1F")]
# lower_left_granule = [f for f in files if "h07v06" in f][0]
# upper_right_granule = [f for f in files if "h16v01" in f][0]
# lower_left_da = rioxarray.open_rasterio(f"{download_folder}/{lower_left_granule}")[DATAFIELD_NAME]
# upper_right_da = rioxarray.open_rasterio(f"{download_folder}/{upper_right_granule}")[DATAFIELD_NAME]
# bounds_tuple = lower_left_da.rio.bounds()[:2] + upper_right_da.rio.bounds()[2:]  # left, bottom, right, top: float  
# del lower_left_da, upper_right_da, lower_left_granule, upper_right_granule, files, download_folder
# ------------------------------------------------------------------------------------------------

# for download_folder in download_folder_list:
#     extract_modis(download_folder)  # Serial processing
tic()
Parallel(n_jobs=cores)(delayed(extract_modis)(download_folder, daF) for download_folder in download_folder_list)
time.sleep(5)  # wait for all threads to finish
# with n_jobs=-1 give out of memory error on Discover, likely because it is using all cores even though only a few
# requested by Slurm. Hence, write a explicit number rather than -1. Keep cores below 30 to avoid NodeFailError  
# Runtime ~10 mins for 1 year with 26 cores.
logging.info(f"Finished Reproj match of Modis (MOD10A1F) to LIS Resolution. {toc()}")
# Cleanup
del lis_template
del nan_mask
del daF
# ========================================================================================================
# ========================================================================================================

logging.info("Start Concatenation Reproj matched Modis (MOD10A1F) along time dimension.")
# nc_files = [f for f in os.listdir(modis_folder) if f.endswith(".nc4") and f.startswith("MOD10A1F")]  # sarith has nc4 extension
nc_files = [f for f in os.listdir(clip_folder) if f.endswith(".nc") and f.startswith("M")]  # no more MOD10A1F because we may have MYD files
nc_files.sort(key=lambda x: int(x.split(".")[1].split("_")[-1]))  # Sort ascending order
logging.info(f"Count of netcdf files: {len(nc_files)}")


def concat_along_time(nc_file):
    """ Concatenate the MODIS CGF data along time dimension
        Required for parallelization using joblib
    """
    begin_dt = nc_file.split(".")[1].split("_")[-1]
    begin_year = int(begin_dt[:4])
    begin_days = int(begin_dt[4:])
    begin_ymd = datetime.datetime(begin_year, 1, 1) + datetime.timedelta(begin_days - 1)
    da_temp = xr.open_dataset(f"{clip_folder}/{nc_file}", decode_coords="all")["CGF_NDSI_Snow_Cover"].squeeze()
    da_temp["time"] = begin_ymd
    return da_temp


# Here da is a list of results. The results should include data arrays
da = Parallel(n_jobs=cores)(delayed(concat_along_time)(nc_file) for nc_file in nc_files)
da = xr.concat(da, dim='time')
logging.info(f"\tFinished concatenation of modis data. Shape = {da.shape}")

# da.data = np.round(da.data * 100)  # Just rounding also saves space for float32 compared to unrounded. Also used uint8.
# da = da.fillna(255)  # a)required uint8; else we only get 0 and 1 values upon conversion. 
# da.data = da.data.astype(np.uint8)  # b)required uint8. np.float32; np.(int8 or ubyte): 0-255 (unsigned); np.byte: -128 to 127)  data was previously float64; but even int8 should be enough (TODO)
da.data = da.data.astype(np.float32)
ds = xr.Dataset({"SCF": da})  # MODSCAG so we can save to netcdf or append to SEUP dataset
# ds["SCF"].attrs["_FillValue"] = 255  # c)required uint8. update fill value from nan to 255; do after creating dataset so it is nested with variable.
logging.info("Saving concatenated MODIS_CGF")
ds = ds.chunk(chunks=chunk)
ds = ds.compute()  # Load numpy array in memory before saving. Faster but can create memory error. If memory error, then comment this line
ds.to_netcdf(f"{lis_folder}/WY{water_year}/SCF.nc", encoding={"SCF": {'zlib': True}}, format="NETCDF4", engine="h5netcdf")  # lis/ keep filename same as variable name. required for new rasterstack in Julia.
# ds.to_zarr(f"{combined_modis_folder}/{water_year}_modis.zarr")  # time consuming; 1.5 hours
logging.info(f"\t Finished saving final Modis concatenated file: {lis_folder}/WY{water_year}/SCF.nc")
logging.info("Finished Job.")

if __name__ == "__main__":
    """ Call the main function to get MODIS_CGF north America matching SEUP resolution
    """
    # main()
