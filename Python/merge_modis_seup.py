#!usr/bin/env python

""" Clip SEUP and MODIS CGF over North America by watershed polygon to generate Analysis Ready Data for Blender
    Usage: python ../merge_modis_seup.py   # OLD: clip_by_watershed.py
    Nov 04, 2022: Populated missing days of data from previous day; this was causing error in Blender run
    Nov 06, 2022: Updated with ffill and bfill with limit=1 for missing days. 
        TODO: Check if better to interpolate between two values
    Nov 08, 2022: Updating for North America MODIS CGF that I newly generated using rioxarray array_merge
    Nov 19, 2022: Updating to generate Blender ready data for whole North America; ie, NoahMP + MODIS_CGF in on nc file
    Feb 25, 2023: Last successful run incorporating Polar Nights Fix
    Jun 29, 2023: using xarray to open and concatenate MODIS files instead of Rioxarray due to stripping on outputs  

    1. Concatenate NA-MODIS files along time dimension
    2. Back and Forward fill missing data
    3. Combine with SEUP data and save the file that is ready for Blender Run
    4. [optional] Clip with watershed is currently not active
    5. Runtime: 15 mins with 30 GB and 1 core on Discover

    Inputs
    =======
    Daily MODIS mosaics over NA processed by "process_modis_cgf.py" script 
    1. modis_folder = f"{base_folder}/Modis/CGF_NDSI_Snow_Cover/NA{water_year}"   


    Outputs
    =======
    1. ../Blender/Modis/CGF_NDSI_Snow_Cover/modis_scf2016_original.nc   # intermediate/temporary
    2. ../Blender/CGF_NDSI_Snow_Cover/modis_scf2016.nc                  # intermediate
    3. ../Blender/Inputs/WY_merged/2016_seup_modis.nc                   # Final. Use for Blender run.  


"""
import platform
import os
import sys
import time
import datetime
import numpy as np

# from shapely.geometry import mapping 
# import geopandas as gpd
# import rioxarray
import xarray as xr
import logging
import argparse

ss = platform.system()
node = platform.node()
logging.basicConfig(filename='OUT/merge_modis_seup.log', level=logging.INFO, format='%(asctime)s:%(message)s')

parser = argparse.ArgumentParser(description='Merge SEUP and MODIS data.')
parser.add_argument('water_year', help='Water Year for processing', type=str)
parser.add_argument('--RES', help='Grid size (Resolution) of SEUP Inputs', type=str)
args = parser.parse_args()
water_year = args.water_year
RES = args.RES  # "050"  # Grid Resolution: 050 = 0.050 degree (5km); 025 = 0.025 degree; 010 = 0.01 degree; 001 = 0.001 degree (100 meters) etc


if 'L-JY0R5X3' in node and "Windows" in ss:
    coressd_folder = "D:/coressd"
elif 'L-JY0R5X3' in node and "Linux" in ss:
    coressd_folder = "/mnt/d/coressd"
elif "borg" in node and "Linux" in ss:
    coressd_folder = "/discover/nobackup/projects/coressd"
elif "asc.ohio-state.edu" in node and "Linux" in ss:
    coressd_folder = "/fs/project/howat.4/yadav.111/coressd"
elif ".osc.edu" in node and "Linux" in ss:
    coressd_folder = "/fs/ess/PAS1785/coressd"
else:
    print("Unknow computer system. coressd folder NOT set")
    assert(False)
base_folder = f"{coressd_folder}/Blender"
logging.info(f"base_folder: {base_folder}")

# Following two folders will be created by this script as required (below).
temp_folder = f"{base_folder}/Modis/MOD10A1F/temp"  # To save concatenated with origianl flags for QA/QC in future; but not really necessary for this or other workflow
out_folder = f"{base_folder}/Inputs_{RES}/WY_merged"  # To save the final file that is ready for Blender Run
logging.info(f"out_folder: {out_folder}")

# Part I: Read daily MODIS-CGF North America mosiaics and create WaterYear data [ie, concatenate by time]
# ========================================================================================================
# modis_folder = f"{base_folder}/MOD10A1F.061_clip/WY16"  # C:/Github/Blender/MOD10A1F/WY16"  # ones processed by sarith copied here for prototyping convenience
modis_folder = f"{base_folder}/Modis/MOD10A1F/NA{water_year}_{RES}"  # Read daily North America mosaicked ModisCGF. [Older = CGF_NDSI_Snow_Cover]  
# hdf_files = [f for f in os.listdir(modis_folder) if f.endswith(".nc4") and f.startswith("MOD10A1F")]  # sarith has nc4 extension
hdf_files = [f for f in os.listdir(modis_folder) if f.endswith(".nc") and f.startswith("MOD10A1F")]  # I have nc extension
# Sort in ascending order of data
# nc_files.sort(key = lambda x: pd.to_datetime(x.split("_")[2].split(".")[0]))
# hdf_files.sort(key=lambda x: int(x.split(".")[1][1:]))
hdf_files.sort(key=lambda x: int(x.split(".")[1].split("_")[-1]))
hdf_files = hdf_files[40:270]  # 60:240 for testing; remove later
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
        # da = rioxarray.open_rasterio(f"{modis_folder}/{hdf_file}").squeeze()
        da = xr.open_dataset(f"{modis_folder}/{hdf_file}", decode_coords="all")["CGF_NDSI_Snow_Cover"].squeeze()
        # da = da["CGF_NDSI_Snow_Cover"].squeeze()
        da["time"] = begin_ymd
    else:
        # da_temp = rioxarray.open_rasterio(f"{modis_folder}/{hdf_file}").squeeze()
        da_temp = xr.open_dataset(f"{modis_folder}/{hdf_file}", decode_coords="all")["CGF_NDSI_Snow_Cover"].squeeze()
        # da_temp = ds["CGF_NDSI_Snow_Cover"].squeeze()
        da_temp["time"] = begin_ymd
        da = xr.concat([da, da_temp], dim='time')  # , dim='time'
    count += 1
logging.info(da.shape)
# da = da.drop(["band", "spatial_ref"])  # after shifting to xarray, these coords were introduced; so remove before combining with seup
da = da.drop_vars(["band"])  # ValueError: These variables cannot be found in this dataset: ['spatial_ref'] 
da = da.drop_duplicates(dim="time")  # required for Sarith script; may not be required for mine as there should be no duplicate dates 
os.makedirs(temp_folder, exist_ok=True)
# ds = xr.Dataset({"MODSCAG": da})  # so we can save to netcdf or append to SEUP dataset
# ds.to_netcdf(f"{temp_folder}/modis_scf{water_year}_original_{RES}.nc", encoding={"MODSCAG": {'zlib': True}})
# logging.info(f"Saved NA WaterYear MODIS_CGF with original flags: {temp_folder}/modis_scf{water_year}_original_{RES}.nc")
# del ds  # clear from memory
# Memeory error after this point somewhere for 1 km run.
# -----------------------------------------------------------------------------------------------------------------------------------

# # Fill in All_Nan days: due to all nans, this creates problem in Blender, so fill them with previous day of data
# da.sel(time="2015-10-24").data[:] = da.sel(time="2015-10-23").data
# da.sel(time="2015-12-09").data[:] = da.sel(time="2015-12-08").data
# da.sel(time="2016-05-02").data[:] = da.sel(time="2016-05-01").data
# da.sel(time="2016-06-24").data[:] = da.sel(time="2016-06-23").data
# logging.info(da.shape)
# logging.info(f"CRS for concatenated da: {da.rio.crs}")
# Replace fillvalues with Nans and various flags with zero [Replacements already done in prior script]
# da.data = da.data.astype(float)
# Replace with Nans: They both look to give same result
# da.data[da.data >100] = np.nan  #this causes most of the pixels to be nan. Hence, cannot be used in Blender
# da.data[da.data == da._FillValue] = np.nan  # [Alreay done in prior script] TODO : comment this line
# da.data[da.data > 100] = 0  # Give rest of the flags value of zero snow! [Alreay done in prior script]
# da.data[da>100] = np.nan  # this does not work for multidimensional data (ie, da concatenated with more than 1 time period)
# logging.info("2 converted nans ans zeros")
# # May 02, 2022: Convert NDSI to snow cover fraction (FRA); FRA = 0.06 + 1.21 * NDSI
# da.data = 0.06 + 1.21 * da.data
# # Correct of biases introduced by above equation
# da.data[da.data <= 0.06] = 0
# da.data[da.data > 1] = 1  # fix this in process_modis_cgf.py script

# # Now ffil and/or Bfill to missing/corrupt pixels  
# # ffill : Fill NaN values by propagating values forward; need bottleneck
# da = da.ffill(dim="time", limit=1)
# da = da.bfill(dim="time", limit=1)  # required if first day empty. If two consecutive days empty, the second empty will be filled
# logging.info(f"CRS for concatenated da: {da.rio.crs}")
# Nov 19, 2022: Save this one round before chaning pixel values; this may directly be used in final run after updating no-data pixel values
ds = xr.Dataset({"MODSCAG": da})  # so we can save to netcdf or append to SEUP dataset
ds.to_netcdf(f"{temp_folder}/modis_scf{water_year}_{RES}.nc", encoding={"MODSCAG": {'zlib': True}})
logging.info(f"Saved NA WaterYear MODIS_CGF with original flags: {temp_folder}/modis_scf{water_year}_{RES}.nc")
ds = ds.reset_coords(drop=True)  # required for combining with SEUP
time.sleep(5)  # to avoid memory error; not sure if this is needed
del da, da_temp  # clear from memory
time.sleep(5)  # to avoid memory error; not sure if this is needed

# Part II: Read SEUP Data and combine MODIS-CGF WY to create 1-WY-ARD for Blender Run
# Aside: Reproj match with SEUP was already performed in "process_modis_cgf.py" script  
# ========================================================================================================
seup_ds = xr.open_dataset(f"{base_folder}/Inputs_{RES}/WY/{water_year}_seup.nc", decode_coords="all", chunks="auto")  # added: chunks="auto" for memory error (May 05, 2024)
# Merge Seup and Modis Data
seup_ds = seup_ds.sel(time=ds.time)  # if same, this should also organize data in some order
# MODSCAG_clipped = MODSCAG_clipped.drop_duplicates(dim="time")  # we have one duplicate time value; remove else error in next step
# TODO: duplicate values were from mapping MOD10A1F.061_2015366.nc4 and MOD10A1F.061_2016001.nc4 to 2016/01/01. Serious issue. 
# Investigate further what is causing this: original MODIS data or the processed ones  
seup_ds["MODSCAG"] = ds["MODSCAG"]  # ValueError: cannot reindex or align along dimension 'time' because the (pandas) index has duplicate values
time.sleep(5)
del ds  # clear from memory
time.sleep(5)

# seup_ds_clipped = seup_ds_clipped.sel(time=noah_ds_clip2.time)  # TODO: or select time explicitly here
# seup_ds_clipped = seup_ds_clipped.sel(time=slice("2015-10-01", "2016-09-29"))  # KeyError: "cannot represent labeled-based slice indexer for coordinate 'time' with a slice over integer positions; the index is unsorted or non-unique"
# seup_ds = seup_ds.isel(time=slice(None, -2))  # perhaps only for Sarith's
logging.info("Creating Directory")

os.makedirs(out_folder, exist_ok=True)  # NEW July 3, 2023. this folder no more created in Extract_Seup notebook
# seup_ds.to_netcdf(f"{base_folder}/Inputs/WY_merged/2016_noahmp_cgf.nc")  # use this for Blender run for NoahMP
logging.info("Saving Output nc file")

seup_ds.to_netcdf(f"{out_folder}/{water_year}_seup_modis.nc")  # use this for Blender run for NoahMP
logging.info("Fished preparing ARD for Blender")

# UPTO THIS PART FOR CONCATENATING MODIS CGF FOR NORTH AMERICA
sys.exit(0)
# -----------------------------------------------------------------------------------------------------------------------------------

# # Get Watershed boundary. Make sure this is in geographic crs
# gdf = gpd.read_file(f"{base_folder}/coordinates/TuolumneWatershed/GlobalWatershed.shp")  # already projected: NAD_1983_Albers
# gdf.to_crs("EPSG:4326", inplace=True)  # unproject to lat/lon

# Clip Annual merged SEUP data


# seup_ds.rio.write_crs(gdf.crs, inplace=True)  # Must for rio.clip operation but not necessary after decoding_coords
# seup_ds.rio.write_crs("epsg:4326", inplace=True)  # Explicit crs is better  

# # Why this is needed here?
# seup_ds_clipped = seup_ds.rio.clip(gdf.geometry.apply(mapping), gdf.crs)  # A list of geojson geometry dicts
# seup_ds_clipped.to_netcdf(f'{base_folder}/Inputs/WY_merged/2016_clip.nc')

# ================================================================================================
# Clip MODIS CGF for watershed using reprojection match
# noah_ds_clip = rioxarray.open_rasterio(f'{noahmp_folder}/WY_merged/2016_clip.nc')
# noah_ds_clip = rioxarray.open_rasterio('C:/Github/coressd/Blender/Delete/2016_clip.nc')
# # Select one sample variable and just one time
# template_raster = noah_ds_clip.isel(time=0)["SWE_tavg"]
# ================================================================================================



modis_ws_folder = f"{base_folder}/CGF_NDSI_Snow_Cover/ws"  # Save MODIS data clipped to watershed; just intermediate step; not directly used by Blender
os.makedirs(modis_ws_folder, exist_ok=True)

# CONTINUE NEXT PORTION FOR CLIPPING BY WATERSHED
# Reprojection match with SEUP Data and Clip by Tuolumne Watershed (cannot yet work with whole North America as we do not have all the tiles merged to NA
# This should take care of size, shape, resolution, clips etc.
# # Get Clipped NoahMP file: we will match MODIS to this data (resolution, spatial ref etc)
# noahmp_folder = f"{base_folder}/Inputs"
# # clipped = clipped.sel(time = ds.time)  # match the time period from MODSCAG data; it has one less
# # Alternative: clipped.sel(time=slice("2015","2016-09-29"))
# # template_raster = clipped
# # noah_ds_clip = rioxarray.open_rasterio(f'{noahmp_folder}/WY_merged/2016_clip.nc')
seup_ds_clipped = rioxarray.open_rasterio(f'{base_folder}/Inputs/WY_merged/2016_clip.nc')  # open again if processed separate from NoahMP/SEUP
# Select one sample variable and just one time
template_raster = seup_ds_clipped.isel(time=0)["SWE_tavg"]
# template_raster.data==np.nan
nan_mask = np.isnan(template_raster.data)  # we will use this mask to mask out areas outside the watershed for MODIS data

# MODSCAG_clipped = ds.rio.reproject_match(template_raster)  # resampling=Resampling.bilinear; still need to convert to float etc for nan, so use da
MODSCAG_clipped = da.rio.reproject_match(template_raster)  # seems similar to ds
logging.info(f"After Reproj match, MODSCAG_clipped CRS = {MODSCAG_clipped.rio.crs}")

# Fill outside are with Nans extracted from NoahMP clipped data
# MODSCAG_clipped.data[nan_mask] = np.nan  # didn't work for mulitdimentional data
# MODSCAG_clipped1= xr.where(nan_mask, np.nan, MODSCAG_clipped, keep_attrs=True) #This also works
MODSCAG_clipped = xr.where(~nan_mask, MODSCAG_clipped, np.nan, keep_attrs=True)
# MODSCAG_clipped.data = MODSCAG_clipped.data/100  # Blender requires cloud fraction in zero to one but this step seems already done
MODSCAG_clipped = xr.Dataset({"MODSCAG":MODSCAG_clipped})  # so we can save to netcdf or append to SEUP dataset
MODSCAG_clipped.to_netcdf(f"{modis_ws_folder}/CGF_clipped.nc")  # 260 MB size
# MODSCAG_clipped.to_netcdf(f"{base_folder}/MOD10A1F.061_clip/CGF_clipped.nc")  # 260 MB size
# clip_folder = f"{base_folder}/{DATAFIELD_NAME}/NA"
logging.info("4 Saved CGF_clipped")


# Finally, append MODIS Snow cover data to NoahMP clips
# noahmp_folder = f"{base_folder}/Inputs"
# clipped = clipped.sel(time = ds.time)  # match the time period from MODSCAG data; it has one less
# Alternative: clipped.sel(time=slice("2015","2016-09-29"))
# template_raster = clipped
# noah_ds_clip = rioxarray.open_rasterio(f'{noahmp_folder}/WY_merged/2016_clip.nc')
# noah_ds_clip2 = xr.open_dataset(f'{noahmp_folder}/WY_merged/2016_clip_noahmp_modscag.nc')  # to get time from here [should not be necessary in future. For now due to incomplete old record (due to leap year)]
# noah_ds_clip =    xr.open_dataset(f'{noahmp_folder}/WY_merged/2016_clip.nc')  # temp only, to subset time

# seup_ds_clipped = xr.open_dataset(f'{base_folder}/Inputs/WY_merged/2016_clip.nc')  # open again if processed separate from NoahMP/SEUP
# MODSCAG_clipped = xr.open_dataset(f"{modis_ws_folder}/CGF_clipped.nc")  # new from my workflow
# replaced with rioxarray open_rasterio method because of discrepancy in y-axis coords resulting in some Nans in final product
seup_ds_clipped = rioxarray.open_rasterio(f'{base_folder}/Inputs/WY_merged/2016_clip.nc') 
MODSCAG_clipped = rioxarray.open_rasterio(f"{modis_ws_folder}/CGF_clipped.nc")

# MODSCAG_clipped = xr.open_dataset(f"{base_folder}/MOD10A1F.061_clip/CGF_clipped.nc")  # from Sarith's workflow
logging.info("5 opened clipped files")
# Drop band from coordinates
# MODSCAG_clipped.drop_vars("band") : this also works
MODSCAG_clipped = MODSCAG_clipped.reset_coords(drop=True)
seup_ds_clipped = seup_ds_clipped.sel(time=MODSCAG_clipped.time)  # if same, this should also align time in some order
# MODSCAG_clipped = MODSCAG_clipped.drop_duplicates(dim="time")  # we have one duplicate time value; remove else error in next step
# TODO: duplicate values were from mapping MOD10A1F.061_2015366.nc4 and MOD10A1F.061_2016001.nc4 to 2016/01/01. Serious issue. 
# Investigate further what is causing this: original MODIS data or the processed ones  
seup_ds_clipped["MODSCAG"] = MODSCAG_clipped["MODSCAG"]  # ValueError: cannot reindex or align along dimension 'time' because the (pandas) index has duplicate values
# seup_ds_clipped = seup_ds_clipped.sel(time=noah_ds_clip2.time)  # TODO: or select time explicitly here
# seup_ds_clipped = seup_ds_clipped.sel(time=slice("2015-10-01", "2016-09-29"))  # KeyError: "cannot represent labeled-based slice indexer for coordinate 'time' with a slice over integer positions; the index is unsorted or non-unique"
seup_ds_clipped = seup_ds_clipped.isel(time=slice(None, -2))  # perhaps only for Sarith's
seup_ds_clipped.to_netcdf(f"{base_folder}/Inputs/WY_merged/2016_clip_noahmp_cgf.nc")  # use this for Blender run for NoahMP
logging.info("6 Fished preparing ARD for Blender")


if __name__ == "__main__":
    """ Call the main function to get MODIS_CGF north America matching SEUP resolution
    """
    # main()
    print(f"Create Analysis Ready Data for Blender: {datetime.datetime.now()}")
