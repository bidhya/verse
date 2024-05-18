#!usr/bin/env python

""" My custorm script to extract MODIS CGF over North America matching the SEUP Data 

    Usage       : python ../process_modis_cgf.py 2016 --cores 8
    slurm_scipt : ../Github/Slurm_Blender/b_process_modis_cgf.sh [params = 2015, 22]

    Resources
    =========
    Need ~28 GB per core
    For 1 year with 32 GB/8cores runtime ~20 minutes, with output of ~6GB

    NB:
    ===
    New May 07, 2023: Save the original DNs. Float, convertion etc in subsequent script. This will also save space
    ie, no more reprojection match with SEUP in this script

    TODO: Decide on what to do with cores; pass it or hardcode to -1

    Input
    =====
    Raw downloaded files in original sinusoidal projection
    modis_download_folder   = /discover/nobackup/projects/coressd/OSU/MOD10A1F.061/MODIS_Proc/download_snow     # /2016001/001
    template_raster         = /discover/nobackup/projects/coressd/Blender/Inputs/combined/SWE_tavg/201509.nc    # HARDCODED and required      
    Tree cover fraction     = /discover/nobackup/projects/coressd/Blender/Modis/MOD44B/Percent_Tree_Cover/NA/MOD44B.A{water_year}065.061.nc

    Outputs 
    =======
    Folder: /discover/nobackup/projects/coressd/Blender/Modis/CGF_NDSI_Snow_Cover/
    NA2015_mosaic
    ------------- 
        this is temporary file only and consider deleting (or not saving)
        each file is a daily mosaic over North America in original Sinusoidal projection
        filename example: MOD10A1F.061_2015249.nc

    NA2015
    ------ 
        mosaics from previous step is resampled and matched to SEUP data
        filename example: MOD10A1F.061_2015249.nc

"""
import os
# from datetime import timedelta
import datetime
import pandas as pd

# from giuh_helpers import tic, toc
import numpy as np
import xarray as xr
import rioxarray
from rioxarray.merge import merge_arrays
# from rasterio.enums import Resampling  # for cubic, bilinear resampling etc during reprojection
import logging
import argparse
from joblib import Parallel, delayed

parser = argparse.ArgumentParser(description='Extract MODIS CGF files for North America.')
parser.add_argument('water_year', help='Water Year for processing', type=str)
parser.add_argument('--log_name', help='Name of Log file', type=str, default='process_modis_cgf')
parser.add_argument('--cores', help='Number of cores to use for multiprocessing', type=int, default=2)  # don' use -1 on DISCOVER  
parser.add_argument('--RES', help='Grid size (Resolution) of SEUP Inputs', type=str)
args = parser.parse_args()
water_year = args.water_year
cores = args.cores
log_name = args.log_name
log_name = f"{log_name}{water_year}.log"
RES = args.RES  # "050"  # Grid Resolution: 050 = 0.050 degree (5km); 025 = 0.025 degree; 010 = 0.01 degree; 001 = 0.001 degree (100 meters) etc

logging.basicConfig(filename=f'out/{log_name}', level=logging.INFO, format='%(asctime)s:%(message)s')
logging.info('  ')
logging.info('-------------------------START LOGGING--------------------------------------')
logging.info(f'water_year={water_year}     \ncores={cores}')

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
        NOT CURRENTLY USED  
    """
    begin_year = int(modis_date[:4])
    begin_days = int(modis_date[4:])
    begin_ymd = datetime.datetime(begin_year, 1, 1) + datetime.timedelta(begin_days - 1)
    return begin_ymd


base_folder = "/discover/nobackup/projects/coressd"
modis_download_folder = f"{base_folder}/OSU/MOD10A1F.061/MODIS_Proc/download_snow"  # /2016001/001
# Generate the date range for WY2016 and convert to DOY format of MODIS naming convention
start_date = datetime.datetime.strptime(f"{int(water_year) - 1}-10-01", "%Y-%m-%d")  # 
end_date = datetime.datetime.strptime(f"{water_year}-09-30", "%Y-%m-%d")  # 
date_generated = pd.date_range(start_date, end_date)
# print(date_generated.strftime("%d-%m-%Y"))
year_doy_list = list(date_generated.strftime("%Y%j"))  # Does this match Day of year for MODIS nameing convention?  
# year_doy_list = ["2016001", "2016002", "2016003"]
# year_doy_list = os.listdir(modis_download_folder)

# Global variablees
# Get a template raster for reprojection match. This is the SEUP data. Any nc file from ../coressd/Blender/Inputs_050/combined/SWE_tavg will work.
seup_folder = f"{base_folder}/Blender/Inputs_{RES}"  # NoahMP
template_raster = xr.open_dataset(f'{seup_folder}/combined/SWE_tavg/{water_year}04.nc')  # Use April month  #OLD Harcoded: 201509, 200109  
# Select just one day data
# noah_ds_clip = noah_ds_clip.sel(time = noah_ds_clip2.time)
template_raster = template_raster.isel(time=0)
template_raster.rio.write_crs("EPSG:4326", inplace=True)  # TODO : check if this should Plate Caree  projection  
nan_mask = np.isnan(template_raster["SWE_tavg"].data)  # we'll use this nan mask to replace MODIS data so they match exactly with noahmp
# variables = ["CGF_NDSI_Snow_Cover", "Cloud_Persistence", "Basic_QA", "Algorithm_Flags_QA", "MOD10A1_NDSI_Snow_Cover"]
# DATAFIELD_NAME = variables[0]  # for now just check one important varialbes; may need to check all

DATAFIELD_NAME = "CGF_NDSI_Snow_Cover"  # or just hardcode the variable of interest
# Create clip folder to save clipped modis files. Initially clipping from global mosaic, hence, named clipped. 
# No more strictly clipping in new workflow because I manually select subset of MODIS tiles  
mosaic_folder = f"{base_folder}/Blender/Modis/MOD10A1F/mosaic_{water_year}_{RES}"  # Feb 07, 2024: inserted temp for mosaics  
os.makedirs(mosaic_folder, exist_ok=True)
clip_folder = f"{base_folder}/Blender/Modis/MOD10A1F/NA{water_year}_{RES}"
os.makedirs(clip_folder, exist_ok=True)

# Get Tree cover fraction data
# daF = rioxarray.open_rasterio(f'{base_folder}/Blender/Modis/MOD44B/Percent_Tree_Cover/NA/MOD44B.A{water_year}065.061.nc').squeeze()
daF = rioxarray.open_rasterio(f'{base_folder}/Blender/Modis/MOD44B/Reproj/MOD44B.A{water_year}065.061.tif').squeeze()
logging.info(f"daF Shape = {daF.shape}")
daF.data = daF.data.astype(float)
daF.data[daF.data == daF._FillValue] = np.nan
# # da.data[da.data == 200] = np.nan  # water ; maybe it water implies zero forest cover?; so don't use nan here
daF.data[daF.data > 100] = 0  # Give rest of the flags value of zero tree cover percent
# daF.data = daF.data / 100
# daF = daF.rio.reproject_match(template_raster)  # match directly to seup resolution   or TODO better to match to MODIS SCF resolution

# # Replace MOD44B with VCF5KYR (0.05 degree 1 global file per year, ending in 2016)
# daF = rioxarray.open_rasterio(f'{base_folder}/Blender/Modis/VCF5KYR/downloads/VCF5KYR_{water_year}.tif')  # squeeze not required here.  
# daF = daF.sel(band=1)  # tree-cover fraction in band 1.  
# daF = daF.rio.reproject_match(template_raster)  # match directly to seup resolution   
# daF.data = daF.data / 100  # opposite of above. but works without having to convert to float.  todo: more verification.  

# TODO: First just merge and save data in original format (both cgf and tree cover), ie, as follows
""" Even same just one merged files at temporary template: required for daF reproj match
da = merge_arrays([rioxarray.open_rasterio(f'{download_folder}/{tif}')[DATAFIELD_NAME] for tif in subset])  # chunks={'x':2**13, 'y':2**13}
da = da.squeeze()
# New May 07, 2023: Save the original DNs. Float, convertion etc in subsequent script. This will also save space
# Conver to Dataset, so this is proper nc file
da = xr.Dataset({DATAFIELD_NAME: da})
da.to_netcdf(f"{clip_folder}/{out_tif_name}.nc")

"""
logging.info("Start Reproj match of MODIS.")


def extract_modis(download_folder, daF):
    """ Extract the MODIS CGF matching the SEUP data ready for Blender Julia run
    file_path
    =========
    Full path to netcdf file that was processed at NASA discover
    not tested: hdf4 file not tested but should work seamlessly
    """
    files = [f for f in os.listdir(download_folder) if f.endswith(".hdf") and f.startswith("MOD10A1F")]
    # Get the subset of modis 10 degree files for this day
    subset = [f for f in files if f.split(".")[2] in modis_tile_list]
    hdf_filename = subset[0]
    product, year_doy, tile, version, _, _ = hdf_filename.split(".")  # Create filenme to one matching the pattern processed by Sarith
    logging.info(f"\t For {year_doy}, total number of tiles = {len(files)} and subset file count = {len(subset)}")

    out_tif_name = f"{product}.{version}_{year_doy[1:]}"
    if not os.path.exists(f"{clip_folder}/{out_tif_name}.nc"):
        # da = merge_arrays([rioxarray.open_rasterio(f'{download_folder}/{tif}')[DATAFIELD_NAME] for tif in subset])  # chunks={'x':2**13, 'y':2**13}
        da = merge_arrays([rioxarray.open_rasterio(f'{download_folder}/{tif}')[DATAFIELD_NAME] for tif in subset], bounds=bounds_tuple)
        da = da.squeeze()
        da.rio.to_raster(f"{mosaic_folder}/{out_tif_name}.tif", driver="COG")
        # # New May 07, 2023: Save the original DNs. Float, convertion etc in subsequent script. This will also save space
        # # Conver to Dataset, so this is proper nc file
        # # TODO Feb 07, 2024: Comment next two lines. Was saved only to inspect data. Though 1 sample mosaic seems used for MOD44B data.  
        # ds_mosaic = xr.Dataset({DATAFIELD_NAME: da})
        # encoding = {DATAFIELD_NAME: {'zlib': True}}
        # ds_mosaic.to_netcdf(f"{mosaic_folder}/{out_tif_name}.nc", encoding=encoding)

        # da.rio.write_crs("epsg:4326", inplace=True)  # not used yet
        # mosaic_tif = f"{mosaic_tif_folder}/{out_tif_name}.tif"  # This only temporary to make sure everything checks out
        # da.rio.to_raster(mosaic_tif)
        # More processing before reprojection_match operation
        da.data = da.data.astype(float)
        # Replace with Nans: They both look to give same result
        # da.data[da.data >100] = np.nan  #this causes most of the pixels to be nan. Hence, cannot be used in Blender
        da.data[da.data == da._FillValue] = np.nan
        da.data[da.data == 211] = np.nan  # Night 
        da.data[da.data == 239] = np.nan  # Ocean 

        da.data[da.data > 100] = 0  # Give rest of the flags value of zero snow! TODO: change to nan because blender can handle missing now  
        da.data = da.data / 100  # May 17, 2023: moved this from below: da_clipped.data = da_clipped.data / 100
        # a) May 16, 2022: Convert NDSI to snow cover fraction (FRA); FRA = 0.06 + 1.21 * NDSI
        da.data = 0.06 + 1.21 * da.data
        # Correct of biases introduced by above equation
        da.data[da.data <= 0.06] = 0
        da.data[da.data > 1] = 1
        # New 04/20/2024: Reproj match Tree cover data here, so we get higher resolution tree cover fraction with MODIS CGF
        # daF = daF.rio.reproject_match(da)  # match directly to seup resolution   or TODO better to match to MODIS SCF resolution
        # daF.data[daF.data > 100] = 0  # Give rest of the flags value of zero tree cover percent
        daF.data = daF.data / 100

        da.data = da.data / (1 - (daF + np.finfo(float).eps))  # Correct for Forest Tree Cover Fraction. Memory error for 1 km run, hence, using Milan node on Discover.    

        # da.data = da.data/100  # mabye do at the end
        # Clip on original globa data was >2 minutes
        da_clipped = da.rio.reproject_match(template_raster)  # seems similar to ds ; error in rasterio 1.3.8 so keep 1.3.7; was just memory error    
        da_clipped.data[da_clipped.data == da_clipped._FillValue] = np.nan  # because reproj match introduced fill value of 255
        # # b) Correct for Forest Tree Cover Fraction
        # da_clipped.data = da_clipped.data / (1 - (daF + np.finfo(float).eps))  # ValueError: operands could not be broadcast together with shapes (14401,24000) (14401,24001)
        da_clipped.data[da_clipped.data > 1] = 1

        # TODO: We will get some 255 values for nodata; check if other values also included (due to interpolation)
        da_clipped.data[nan_mask] = np.nan  # done again, because reproj match will fill new nodata with 255 (perhaps this value from attrs)
        # da_clipped.data = da_clipped.data / 100
        # TODO: Fix attrs, no-data etc here so that data is self-describing
        # Conver to Dataset, so this is proper nc file
        ds_clipped = xr.Dataset({DATAFIELD_NAME: da_clipped})
        encoding = {DATAFIELD_NAME: {'zlib': True}}
        ds_clipped.to_netcdf(f"{clip_folder}/{out_tif_name}.nc", encoding=encoding)
        # noah_ds_clip[var] = da_clipped  # Append Modis data to clipped dataset


# Make of List of downloaded sub-folder. To make it easy to parallelize using Joblib
download_folder_list = []
for year_doy in year_doy_list:
    year = year_doy[:4]  # "2016"
    doy = year_doy[4:]  # "215"
    download_folder = f"{modis_download_folder}/{year}{doy}/{doy}"
    download_folder_list.append(download_folder)
    # extract_modis(download_folder)
# ------------------------------------------------------------------------------------------------
# Lower left (h07v06) and Upper Right (h16v01) MODIS Granule for North America
# Use bounds from this two files in mosaicking. esle there is one extra pixel in x, y or both direction.
# This should help constrain and we always get same size of output raster
download_folder = download_folder_list[182]  # April 01, 2015
# DATAFIELD_NAME = "CGF_NDSI_Snow_Cover"
files = [f for f in os.listdir(download_folder) if f.endswith(".hdf") and f.startswith("MOD10A1F")]
lower_left_granule = [f for f in files if "h07v06" in f][0]
upper_right_granule = [f for f in files if "h16v01" in f][0]
lower_left_da = rioxarray.open_rasterio(f"{download_folder}/{lower_left_granule}")[DATAFIELD_NAME]
upper_right_da = rioxarray.open_rasterio(f"{download_folder}/{upper_right_granule}")[DATAFIELD_NAME]
bounds_tuple = lower_left_da.rio.bounds()[:2] + upper_right_da.rio.bounds()[2:]  # left, bottom, right, top: float  
del lower_left_da, upper_right_da, lower_left_granule, upper_right_granule, files, download_folder
# ------------------------------------------------------------------------------------------------
# for download_folder in download_folder_list:
#     extract_modis(download_folder)  # Serial processing
Parallel(n_jobs=cores)(delayed(extract_modis)(download_folder, daF) for download_folder in download_folder_list)
# with n_jobs=-1 give out of memory error on Discover, likely because it is using all cores even though only a few
# requested by Slurm. Hence, write a explicit number rather than -1.  
logging.info("Finished Reproj match of MODIS to SEUP Resolution.")

# ========================================================================================================
# ========================================================================================================

logging.info("Start Concatenating MODIS CGF along time dimension.")
# May 08, 2024: Moving MODIS concatenation part from merge_modis_seup.py to here
# TODO: check variable folder names that do not crash existing names ...
# Part I: Read daily MODIS-CGF North America mosiaics and create WaterYear data [ie, concatenate by time]
# ========================================================================================================
# modis_folder = f"{base_folder}/MOD10A1F.061_clip/WY16"  # C:/Github/Blender/MOD10A1F/WY16"  # ones processed by sarith copied here for prototyping convenience
# modis_folder = f"{base_folder}/Blender/Modis/MOD10A1F/NA{water_year}_{RES}"  # Read daily North America mosaicked ModisCGF. [Older = CGF_NDSI_Snow_Cover]  
# modis_folder = clip_folder  #=  f"{base_folder}/Blender/Modis/MOD10A1F/NA{water_year}_{RES}"
# Following two folders will be created by this script as required (below).
# combined_folder = f"{base_folder}/Blender/Inputs_{RES}/WY"  # Modis/MOD10A1F/combined temp To save concatenated with origianl flags for QA/QC in future; but not really necessary for this or other workflow
# seup_folder =     f"{base_folder}/Blender/Inputs_{RES}"  # NoahMP

# hdf_files = [f for f in os.listdir(modis_folder) if f.endswith(".nc4") and f.startswith("MOD10A1F")]  # sarith has nc4 extension
hdf_files = [f for f in os.listdir(clip_folder) if f.endswith(".nc") and f.startswith("MOD10A1F")]  # I have nc extension
# Sort in ascending order of data
# nc_files.sort(key = lambda x: pd.to_datetime(x.split("_")[2].split(".")[0]))
# hdf_files.sort(key=lambda x: int(x.split(".")[1][1:]))
hdf_files.sort(key=lambda x: int(x.split(".")[1].split("_")[-1]))
# hdf_files = hdf_files[40:270]  # 60:240 for testing; remove later
logging.info(len(hdf_files))
# Idea : Just extract and merge; then later convert to float, nans etc
count = 0
# 1. First concat the files along time dimension
# da = None  # just declaring for red lines below
for hdf_file in hdf_files:
    # these are actually NetCDF files
    # begin_dt = hdf_file.split('.')[1][1:]
    begin_dt = hdf_file.split(".")[1].split("_")[-1]  # hdf_file.split('.')[1][1:]
    begin_year = int(begin_dt[:4])
    begin_days = int(begin_dt[4:])
    begin_ymd = datetime.datetime(begin_year, 1, 1) + datetime.timedelta(begin_days - 1)
    # dt_list.append(begin_ymd)
    if count == 0:
        # da = rioxarray.open_rasterio(f"{clip_folder}/{hdf_file}").squeeze()
        da = xr.open_dataset(f"{clip_folder}/{hdf_file}", decode_coords="all")["CGF_NDSI_Snow_Cover"].squeeze()
        # da = da["CGF_NDSI_Snow_Cover"].squeeze()
        da["time"] = begin_ymd
    else:
        # da_temp = rioxarray.open_rasterio(f"{clip_folder}/{hdf_file}").squeeze()
        da_temp = xr.open_dataset(f"{clip_folder}/{hdf_file}", decode_coords="all")["CGF_NDSI_Snow_Cover"].squeeze()
        # da_temp = ds["CGF_NDSI_Snow_Cover"].squeeze()
        da_temp["time"] = begin_ymd
        da = xr.concat([da, da_temp], dim='time')  # , dim='time'
    count += 1
logging.info(da.shape)
# da = da.drop(["band", "spatial_ref"])  # after shifting to xarray, these coords were introduced; so remove before combining with seup
da = da.drop_vars(["band"])  # ValueError: These variables cannot be found in this dataset: ['spatial_ref'] 
da = da.drop_duplicates(dim="time")  # required for Sarith script; may not be required for mine as there should be no duplicate dates 
ds = xr.Dataset({"MODSCAG": da})  # so we can save to netcdf or append to SEUP dataset
ds = ds.chunk(chunks={"x": 2**8, "y": 2**8})
ds.to_netcdf(f"{seup_folder}/WY/{water_year}_modis.nc", encoding={"MODSCAG": {'zlib': True}})
ds.to_zarr(f"{seup_folder}/WY/{water_year}_modis.zarr")

logging.info(f"Saved concatenated MODIS_CGF: {seup_folder}/WY/{water_year}_modis.nc")
logging.info("Finished Job.")

if __name__ == "__main__":
    """ Call the main function to get MODIS_CGF north America matching SEUP resolution
    """
    # main()
