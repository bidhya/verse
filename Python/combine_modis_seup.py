#!usr/bin/env python

""" Combine with Merged LIS variables with MODIS SCF to produce a 
    final NetCDF file that is ready for Blender Run.

    Runtime: Exclusive Milan node ~ 480 GB ram for 1 km run

    Inputs
    =======
    Daily MODIS mosaics clipped to SEUP resolution over NA. Processed by "process_modis_cgf.py" script 
    modis = ../Blender/Inputs_010/WY/2010_modis.nc     
    seup  = ../Blender/Inputs_010/WY/2010_seup.nc

    Outputs
    =======
    ../Blender/Inputs_010/WY_merged/2016_seup_modis.nc                   
    Use this file for Blender run.  
"""
import platform
import os
# import time
# import datetime, tempfile, sys
# import numpy as np
import xarray as xr
import logging
import argparse
import shutil


ss = platform.system()
node = platform.node()

parser = argparse.ArgumentParser(description='Merge SEUP and MODIS data.')
parser.add_argument('water_year', help='Water Year for processing', type=str)
parser.add_argument('--RES', help='Grid size (Resolution) of SEUP Inputs', type=str)
parser.add_argument('--suffix', help='Temporary. Suffix to add to output filename', type=str, default="")
parser.add_argument('--filetype', help='nc or zarr', type=str, default="nc")

args = parser.parse_args()
water_year = args.water_year
RES = args.RES  # "050"  # Grid Resolution: 050 = 0.050 degree (5km); 025 = 0.025 degree; 010 = 0.01 degree; 001 = 0.001 degree (100 meters) etc
suffix = args.suffix
filetype = args.filetype  # "zarr"  # "nc" or "zarr"

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

# https://betterstack.com/community/questions/how-to-change-time-format-in-python-logging/
logging.basicConfig(filename=f"out/lis_modis_combine{water_year}_{RES}{suffix}.log", level=logging.INFO, format='%(asctime)s : %(message)s')
base_folder = f"{coressd_folder}/Blender"
# Following two folders will be created by this script as required (below).
out_folder = f"{base_folder}/Inputs_{RES}/WY_merged"  # To save the final file that is ready for Blender Run

logging.info(f"base_folder: {base_folder}")
logging.info(f"out_folder: {out_folder}")

tmpdir = os.environ.get("TMPDIR")
# tmpdir = os.environ.get("LOCAL_TMPDIR")    # for NASA Discover
tse_tmpdir = os.environ.get("TSE_TMPDIR")  # os.environ.get("TSE_TMPDIR")  # 250 GB. Fast scratch filesystem using NVMe (flash) 
logging.info(f"copying to HPC/local/scratch folder: {tse_tmpdir}")
# Move MODIS and SEUP nc files to HPC node (where this notebook is running)
# shutil.copy(f"{base_folder}/Modis/MOD10A1F/temp/modis_scf{water_year}_{RES}.nc", tmpdir)
if filetype == "nc":
    shutil.copy(f"{base_folder}/Inputs_{RES}/WY/{water_year}_modis.nc", tse_tmpdir)
    shutil.copy(f"{base_folder}/Inputs_{RES}/WY/{water_year}_seup.nc", tse_tmpdir)
    # Load nc files
    ds = xr.open_dataset(f"{tse_tmpdir}/{water_year}_modis.nc")
    ds = ds.reset_coords(drop=True)  # required for combining with SEUP. TODO check again if this is necessary
    # ds = ds.sel(time=slice(f"{water_year}-04-01", f"{water_year}-04-03"))  # for testing only
    seup_ds = xr.open_dataset(f"{tse_tmpdir}/{water_year}_seup.nc", decode_coords="all")  # , chunks={"time": 50}, chunks="auto"
elif filetype == "zarr":
    shutil.copytree(f"{base_folder}/Inputs_{RES}/WY/{water_year}_modis.zarr", f"{tse_tmpdir}/{water_year}_modis.zarr")
    shutil.copytree(f"{base_folder}/Inputs_{RES}/WY/{water_year}_seup.zarr", f"{tse_tmpdir}/{water_year}_seup.zarr")
    # Load zarr files
    ds = xr.open_zarr(f"{tse_tmpdir}/{water_year}_modis.zarr")
    ds = ds.reset_coords(drop=True)
    seup_ds = xr.open_zarr(f"{tse_tmpdir}/{water_year}_seup.zarr", decode_coords="all")
logging.info(f"Finished copying and loading {filetype} to HPC local folder: {tse_tmpdir}")

# Align SEUP and MODIS data temporally
seup_ds = seup_ds.sel(time=ds.time)  # if same, this should also organize data in some order
seup_ds["MODSCAG"] = ds["MODSCAG"]
seup_ds = seup_ds.chunk(chunks={"x": 2**8, "y": 2**8})
data_size = seup_ds.nbytes / 1e9
logging.info(f"Data size: {data_size:.2f} GB")
del ds
# seup_ds = seup_ds.sel(time=slice("2009-04-01", "2009-04-05"))
# seup_ds = seup_ds.persist()

# -------------------------------------------------------------------------------------------------------------------------------
# Save Combined Output file
os.makedirs(out_folder, exist_ok=True)
if filetype == "zarr":
    logging.info("Save Output ZARR file on HPC/Scratch folder")
    out_fname = f"{water_year}_seup_modis{suffix}.zarr"
    seup_ds.to_zarr(f"{tmpdir}/{out_fname}")
    logging.info("Move ZARR file from HPC/Scratch to project directory")
    shutil.copytree(f"{tmpdir}/{out_fname}", f"{out_folder}/{out_fname}")  # Copy back to project directory
else:
    logging.info("Save Output NetCDF file on HPC/Scratch folder")
    out_fname = f"{water_year}_seup_modis{suffix}.nc"
    encoding = {var: {'zlib': True} for var in seup_ds.data_vars}
    seup_ds.to_netcdf(f"{tmpdir}/{out_fname}", encoding=encoding)
    logging.info("Move NetCDF file from HPC/Scratch to project directory")
    shutil.copy(f"{tmpdir}/{out_fname}", f"{out_folder}/{out_fname}")  # Copy back to project directory
logging.info("Finished Run")
