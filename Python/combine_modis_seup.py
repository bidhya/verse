#!usr/bin/env python

""" Combine with Merged LIS variables with MODIS SCF to produce a final NetCDF file that is ready for Blender Run.
    
    May 2024: Even with zarr inputs, this script is not producing final NetCDF in 24 hours.
    TODO: Try compute to save numpy array on memory (if within 488 GB memory limit of Milan nodes).
    Currently failing in runtime and memory for 1km resolution.
    So plan moving forward is not to use this script. 
    And keep lis and modis variables separate.

    Runtime: Exclusive Milan node ~ 480 GB ram for 1 km run

    Inputs
    =======
    Daily MODIS mosaics clipped to SEUP resolution over NA. Processed by "process_modis_cgf.py" script 
    modis = ../Blender/Modis/WaterYear/2014_modis.zarr #Inputs_010/WY/2010_modis.nc     
    lis   = ../Blender/Inputs_010/lis/2014_lis.zarr  # WY/2010_seup.nc

    Outputs
    =======
    ../Blender/Modis/modis_wy_combined/  # Old for SEUP data: ../Inputs_010/WY_merged/2016_seup_modis.nc                   
    Was planned to sse this file for Blender run. But not working for 1 km resolution.  
"""
import platform
import os
import time
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
# parser.add_argument('--suffix', help='Temporary. Suffix to add to output filename', type=str, default="")
parser.add_argument('--filetype', help='zarr or nc', type=str, default="zarr")  
# process zarr files by default but maybe save the nc as final output for Julia Blender run

args = parser.parse_args()
water_year = args.water_year
RES = args.RES  # "050"  # Grid Resolution: 050 = 0.050 degree (5km); 025 = 0.025 degree; 010 = 0.01 degree; 001 = 0.001 degree (100 meters) etc
filetype = args.filetype  # "zarr"  # "nc" or "zarr"

test_run = True  # True False  # for testing only

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
    assert False


def main():
    # https://betterstack.com/community/questions/how-to-change-time-format-in-python-logging/
    logging.basicConfig(filename=f"out/lis_modis_combine_{water_year}_{RES}.log", level=logging.INFO, format='%(asctime)s : %(message)s')
    base_folder = f"{coressd_folder}/Blender"
    # Following two folders will be created by this script as required (below).
    out_folder = f"{base_folder}/Inputs_{RES}/lis_modis"  # WY_merged To save the final file that is ready for Blender Run

    # Inputs
    lis_folder = f"{base_folder}/Inputs_{RES}/lis"  # NoahMP
    modis_folder = f"{base_folder}/Modis/WaterYear"  # modis_wy_combined
    # chunk = {"x": 2**8, "y": 2**8}
    # chunk = {"y": 2**7}  # use this same chunk size everywhere
    chunk = "auto"

    logging.info(f"base_folder: {base_folder}")
    logging.info(f"out_folder: {out_folder}")

    # tmpdir = os.environ.get("TMPDIR")
    tmpdir = os.environ.get("LOCAL_TMPDIR")    # for NASA Discover
    # tse_tmpdir = os.environ.get("TSE_TMPDIR")  # 250 GB. Fast scratch filesystem using NVMe (flash) 
    logging.info(f"copying to HPC/local/scratch folder: {tmpdir}")
    # Move MODIS and SEUP nc files to HPC node (where this notebook is running)
    shutil.copytree(f"{modis_folder}/{water_year}_modis.zarr", f"{tmpdir}/{water_year}_modis.zarr")
    logging.info("Finished MODIS copying.")
    shutil.copytree(f"{lis_folder}/{water_year}_lis.zarr", f"{tmpdir}/{water_year}_lis.zarr") 
    logging.info("Finished LIS copying. Now opening MODIS file.")
    # Load zarr files
    modis_ds = xr.open_zarr(f"{tmpdir}/{water_year}_modis.zarr")
    # modis_ds = xr.open_zarr(f"{modis_folder}/{water_year}_modis.zarr")  # , chunks=chunk. don't use here becuase we did not chunk while saving it.
    modis_ds = modis_ds.reset_coords(drop=True)
    # modis_ds = modis_ds.sel(time=slice(f"{water_year}-04-02", f"{water_year}-04-05"))

    logging.info("Loading LIS file.")
    lis_ds = xr.open_zarr(f"{tmpdir}/{water_year}_lis.zarr", decode_coords="all", chunks=chunk)
    # lis_ds = xr.open_zarr(f"{lis_folder}/{water_year}_lis.zarr", decode_coords="all", chunks=chunk)
    logging.info("Finished opening LIS file")
    # elif filetype == "nc":
    # shutil.copy(f"{base_folder}/Inputs_{RES}/WY/{water_year}_modis.nc", tse_tmpdir)
    # shutil.copy(f"{base_folder}/Inputs_{RES}/WY/{water_year}_seup.nc", tse_tmpdir)
    # modis_ds = xr.open_dataset(f"{tse_tmpdir}/{water_year}_modis.nc")
    # modis_ds = modis_ds.reset_coords(drop=True)  # required for combining with SEUP. TODO check again if this is necessary
    # # modis_ds = modis_ds.sel(time=slice(f"{water_year}-04-01", f"{water_year}-04-03"))  # for testing only
    # lis_ds = xr.open_dataset(f"{tse_tmpdir}/{water_year}_seup.nc", decode_coords="all")  # , chunks={"time": 50}, chunks="auto"
    # Align SEUP and MODIS data temporally
    lis_ds = lis_ds.sel(time=modis_ds.time)  # if same, this should also organize data in some order
    lis_ds["MODSCAG"] = modis_ds["MODSCAG"]

    # chunk = {"x": 2**8, "y": 2**8}
    # lis_ds = lis_ds.chunk(chunks=chunk)
    data_size = lis_ds.nbytes / 1e9
    logging.info(f"Data size: {data_size:.2f} GB")
    # del modis_ds
    # if test_run:
    #     print("Test Run: Slicing data for testing only")
    # lis_ds = lis_ds.persist()  # takes long time (~5 hours). Test with client.persist(lis_ds)

    # -------------------------------------------------------------------------------------------------------------------------------
    # Save Combined Output file
    os.makedirs(out_folder, exist_ok=True)
    # if filetype == "zarr":
    # out_fname = f"{water_year}_lis_modis.zarr"
    # lis_ds.to_zarr(f"{tmpdir}/{out_fname}")
    # shutil.copytree(f"{tmpdir}/{out_fname}", f"{out_folder}/{out_fname}")  # Copy back to project directory
    # shutil.rmtree(f"{tmpdir}/{out_fname}")
    # else:
    logging.info(f"Saving NetCDF file project directory: {out_folder}")
    out_fname = f"{water_year}_lis_modis.nc"
    encoding = {var: {'zlib': True} for var in lis_ds.data_vars}
    lis_ds.to_netcdf(f"{out_folder}/{out_fname}", encoding=encoding)  # Save to project directory directly
    # lis_ds.to_netcdf(f"{tse_tmpdir}/{out_fname}", encoding=encoding)
    # shutil.copy(f"{tse_tmpdir}/{out_fname}", f"{out_folder}/{out_fname}")  # rather copy from bash itself

    # # Alternative: save as a mfdataset for faster writing
    # # Create a list with the single chunked dataset
    # datasets = [lis_ds]
    # # Create a list with the desired output file path
    # paths = [f"{out_folder}/{out_fname}"]
    # # paths = [f"{tse_tmpdir}/{out_fname}"]

    # # Save the dataset in parallel using save_mfdataset
    # xr.save_mfdataset(datasets, paths, engine="netcdf4", compute=False)
    # # Trigger the parallel write operation
    # delayed_obj = xr.save_mfdataset(datasets, paths, engine="netcdf4", compute=False, encoding=encoding)
    # delayed_obj.compute()
    logging.info("Finished Run")


if __name__ == "__main__":
    """ In some cases, the location of the code that creates the LocalCluster can cause issues. 
        It's generally recommended to create the LocalCluster and Client within a function 
        or the if __name__ == "__main__": block, rather than at the module level.
    """
    # dask cluster must be started in separate module, like main. Else, it will not work for Slurm job.
    from dask.distributed import Client, LocalCluster
    # cluster = LocalCluster()
    cluster = LocalCluster(n_workers=30, memory_limit="15GB")
    client = Client(cluster)
    time.sleep(10)
    # Call the main function.
    main()
