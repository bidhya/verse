#!usr/bin/env python

""" Diagnose error in input file for script running on nasa-discover

    Usage: python /home/yadav.111/Github/giuh/scripts/arctic_hydro_job.py 98 --cores 30 --memory 180gb --runtime 60:00:00 --suffix a
	   python /home/yadav.111/Github/giuh/scripts/arctic_hydro_job.py 95 --cores 10 --memory 80gb --runtime 20:00:00 --maskdir /fs/project/howat.4/moortgat/Masks/2022/Multi/Masks_gage95 --dlext _mask.tif
    
    Notebook to call this script for multiple gages can be found here: ../02_Durand/Unity/00_c_Prepare_gages_4_run.ipynb
"""
import os
# import shutil
import argparse
import logging
# from symbol import pass_stmt
from pyhdf.SD import SD, SDC

parser = argparse.ArgumentParser(description='Check corrupt MODIS CGF corrupt files.')
parser.add_argument('--start_idx', help='Starting Index of downloaed folder (sorted)', type=int)
parser.add_argument('--end_idx', help='Ending index', type=int)
# parser.add_argument('--logname', help='Number of cores', type=str, default='20')

args = parser.parse_args()
start_idx = args.start_idx
end_idx = args.end_idx

base_folder = "/discover/nobackup/projects/coressd/OSU/MOD10A1F.061/MODIS_Proc"
download_folder = f"{base_folder}/download_snow"  # /2016215/215   Github/Blender
logging.basicConfig(filename=f"/discover/nobackup/projects/coressd/Github/Slurm/Modis/out/check_corrupt_{start_idx}.log", level=logging.INFO, format="%(message)s")

variables = ["CGF_NDSI_Snow_Cover", "Cloud_Persistence", "Basic_QA", "Algorithm_Flags_QA", "MOD10A1_NDSI_Snow_Cover"]
DATAFIELD_NAME = variables[0]  # for now just check one important varialbes; may need to check all
# Total folders to process = 8269

# for folder in sorted(os.listdir(download_folder))[start_idx:end_idx]:
#     doy = folder[-3:]
#     modis_folder = f"{download_folder}/{folder}/{doy}"
#     files = [f for f in os.listdir(modis_folder) if f.endswith(".hdf") and f.startswith("MOD10A1F")]
#     logging.info(f"{folder}: {len(files)}")
#     corrupt_files = []
#     for hdf_file in files:
#         hdf = SD(f"{modis_folder}/{hdf_file}", SDC.READ)
#         try:
#             data2D = hdf.select(DATAFIELD_NAME)
#         except:
#             corrupt_files.append(hdf_file)
#     if len(corrupt_files) > 0:
#         with open(f"/home/byadav/Github/giuh/out/{folder}.txt", "a") as fi:
#             for line in corrupt_files:
#                 fi.write(line)


def check_hdf_files(folder):
    doy = folder[-3:]
    modis_folder = f"{download_folder}/{folder}/{doy}"
    files = [f for f in os.listdir(modis_folder) if f.endswith(".hdf") and f.startswith("MOD10A1F")]
    logging.info(f"{folder}: {len(files)}")
    corrupt_files = []
    for hdf_file in files:
        hdf = SD(f"{modis_folder}/{hdf_file}", SDC.READ)
        try:
            data2D = hdf.select(DATAFIELD_NAME)
        except:
            corrupt_files.append(hdf_file)
    if len(corrupt_files) > 0:
        # with open(f"/home/byadav/Github/giuh/out/{folder}.txt", "a") as fi:
        with open(f"/discover/nobackup/projects/coressd/Github/Slurm/Modis/out/{folder}.txt", "a") as fi:
            for line in corrupt_files:
                fi.write(line)


for folder in sorted(os.listdir(download_folder))[start_idx:end_idx]:
    check_hdf_files(folder)


def main():
    pass


if __name__ == "__main__":
    """ Call the main function to generate the slurm job script with input 
        files to run the matlab workflow ArcticDEM Hydrology
    """
    main()
