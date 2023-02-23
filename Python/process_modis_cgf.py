#!usr/bin/env python

""" My custorm script to extract MODIS CGF over North America matching the SEUP Data
    Input
    =====
    Raw downloaded files in original sinusoidal projection

    Usage: python ../process_modis_cgf.py 2016 --cores 8
    Resources
    =========
    Need ~28 GB per core
    For 1 year with 32 GB/8cores runtime ~20 minutes, with output of ~6GB

    TODO: Decide on what to do with cores; pass it or hardcode to -1
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
from rasterio.enums import Resampling  # for cubic, bilinear resampling etc during reprojection
import logging
import argparse
from joblib import Parallel, delayed

parser = argparse.ArgumentParser(description='Extract MODIS CGF files for North America.')
parser.add_argument('year', help='Year of modis to process', type=str)
parser.add_argument('--log_name', help='Name of Log file', type=str, default='modis_cgf_clip.log')
parser.add_argument('--cores', help='Number of cores to use for multiprocessing', type=int, default=-1)

# parser.add_argument('--end_idx', help='PGC ID', type=int)
# # parser.add_argument('--logname', help='Number of cores', type=str, default='20')
args = parser.parse_args()
year = args.year
log_name = args.log_name
cores = args.cores
# cores = -1  # defualt should be -1 or just not pass
# log_name = 'modis_cgf_clip.log'

logging.basicConfig(filename=f'{log_name}', level=logging.INFO, format='%(asctime)s:%(message)s')
logging.info('  ')
logging.info('-------------------------START LOGGING--------------------------------------')

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
    """
    begin_year = int(modis_date[:4])
    begin_days = int(modis_date[4:])
    begin_ymd = datetime.datetime(begin_year, 1, 1) + datetime.timedelta(begin_days - 1)
    # print(modis_date, begin_ymd)
    return begin_ymd


def main():
    # root_dir = "C:"  # "/mnt/c"
    # base_folder = f"{root_dir}/Github/coressd/Blender"
    # # modis_folder = f"{base_folder}/MOD10A1F/WY16"  #Github/Blender
    # # year = "2016"
    # modis_folder = f"{base_folder}/MOD10A1F_061/{year}"  #Github/Blender
    # clip_folder = f"{base_folder}/MOD10A1F_061_clip/{year}"  # We will save the combined netcdf file here

    # For NCCS Discover
    # base_folder = "/discover/nobackup/projects/coressd/Blender"
    # modis_folder = f"/discover/nobackup/projects/coressd/OSU/MOD10A1F.061/{year}"
    # clip_folder = f"/discover/nobackup/projects/coressd/OSU/MOD10A1F.061_clip/{year}"
    # New folders
    # root_dir = "C:"  # "/mnt/c"
    # base_folder = f"{root_dir}/Github/coressd"
    base_folder = "/discover/nobackup/projects/coressd"
    modis_download_folder = f"{base_folder}/OSU/MOD10A1F.061/MODIS_Proc/download_snow"  # /2016001/001
    # Generate the date range for WY2016 and convert to DOY format of MODIS naming convention
    start_date = datetime.datetime.strptime("2015-10-01", "%Y-%m-%d")  # 2015-10-01
    end_date = datetime.datetime.strptime("2016-09-30", "%Y-%m-%d")  # 2016-09-30
    date_generated = pd.date_range(start_date, end_date)
    # print(date_generated.strftime("%d-%m-%Y"))
    year_doy_list = list(date_generated.strftime("%Y%j"))
    # year_doy_list = ["2016001", "2016002", "2016003"]
    # year_doy_list = os.listdir(modis_download_folder)

    # Global variablees
    # Get a template raster for reprojection match
    seup_folder = f"{base_folder}/Blender/NoahMP"
    template_raster = xr.open_dataset(f'{seup_folder}/combined/SWE_tavg/201609.nc')
    # Select just one day data
    # noah_ds_clip = noah_ds_clip.sel(time = noah_ds_clip2.time)
    template_raster = template_raster.isel(time=0)
    template_raster.rio.write_crs("EPSG:4326", inplace=True)
    nan_mask = np.isnan(template_raster["SWE_tavg"].data)  # we'll use this nan mask to replace MODIS data so they match exactly with noahmp
    # variables = ["CGF_NDSI_Snow_Cover", "Cloud_Persistence", "Basic_QA", "Algorithm_Flags_QA", "MOD10A1_NDSI_Snow_Cover"]
    # DATAFIELD_NAME = variables[0]  # for now just check one important varialbes; may need to check all
    DATAFIELD_NAME = "CGF_NDSI_Snow_Cover"  # or just hardcode the variable of interest
    # Create clip folder to save clipped modis files
    clip_folder = f"{base_folder}/Blender/{DATAFIELD_NAME}/NA"
    os.makedirs(clip_folder, exist_ok=True)
    # mosaic_tif_folder = f"{base_folder}/{DATAFIELD_NAME}/mosaic_tif"
    # os.makedirs(mosaic_tif_folder, exist_ok=True)

    def extract_modis(download_folder):
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
            da = merge_arrays([rioxarray.open_rasterio(f'{download_folder}/{tif}')[DATAFIELD_NAME] for tif in subset])  # chunks={'x':2**13, 'y':2**13}
            da = da.squeeze()
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
            # da.data = da.data/100  # mabye do at the end
            # Clip on original globa data was >2 minutes
            da_clipped = da.rio.reproject_match(template_raster)  # seems similar to ds 
            # TODO: We will get some 255 values for nodata; check if other values also included (due to interpolation)
            da_clipped.data[nan_mask] = np.nan  # done again, because reproj match will fill new nodata with 255 (perhaps this value from attrs)
            da_clipped.data = da_clipped.data / 100
            # TODO: Fix attrs, no-data etc here so that data is self-describing
            # Conver to Dataset, so this is proper nc file
            ds_clipped = xr.Dataset({DATAFIELD_NAME: da_clipped})
            ds_clipped.to_netcdf(f"{clip_folder}/{out_tif_name}.nc")
            # noah_ds_clip[var] = da_clipped  # Append Modis data to clipped dataset

    # Make of List of downloaded sub-folder. To make it easy to parallelize using Joblib
    download_folder_list = []
    for year_doy in year_doy_list:
        year = year_doy[:4]  # "2016"
        doy = year_doy[4:]  # "215"
        download_folder = f"{modis_download_folder}/{year}{doy}/{doy}"
        download_folder_list.append(download_folder)
        # extract_modis(download_folder)
    # for download_folder in download_folder_list:
    #     extract_modis(download_folder)  # Serial processing
    Parallel(n_jobs=cores)(delayed(extract_modis) (download_folder) for download_folder in download_folder_list)


if __name__ == "__main__":
    """ Call the main function to get MODIS_CGF north America matching SEUP resolution
    """
    main()