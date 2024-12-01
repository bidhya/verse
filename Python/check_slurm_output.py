""" Parse Slum output files to to check  for any errors.  
    Usage: python verse/Python/check_slurm_output.py 2016

    The script can be run directly on login node.
    NB: Recommended to run before combining the intermediate netcdf files using Julia.  
    That is the following scripts:
        d_combine_out_nc_files.sh that calls: julia verse/Julia/combine_nc_files.jl WY$water_year $RES
"""
import os
# from joblib import Parallel, delayed
import logging
# import shutil
import argparse
import platform
ss = platform.system()
node = platform.node()
# print(ss, node)

if 'L-JY0R5X3' in node and "Windows" in ss:
    root_dir = "D:/coressd"
elif 'L-JY0R5X3' in node and "Linux" in ss:
    root_dir = "/mnt/d/coressd"
elif "borg" in node and "Linux" in ss:
    root_dir = "/discover/nobackup/projects/coressd"
elif "discover" in node and "Linux" in ss:
    # when running the script directly on login node in discover
    root_dir = "/discover/nobackup/projects/coressd"
elif "asc.ohio-state.edu" in node and "Linux" in ss:
    root_dir = "/fs/project/howat.4/yadav.111/coressd"
elif ".osc.edu" in node and "Linux" in ss:
    # root_dir = "/fs/scratch/PAS1785/coressd"  # if processing on scratch directory  
    root_dir = "/fs/ess/PAS1785/coressd"
else:
    print("Unknow computer system. coressd folder NOT set")
    assert False

parser = argparse.ArgumentParser(description='Check Slurm job out files for Water Year.')
parser.add_argument('wy', help='Water Year for processing', type=str)
parser.add_argument('--RES', help='Grid size (Resolution) of SEUP Inputs', type=str, default="010")  # SEUP was 050 but may not work here.

args = parser.parse_args()
wy = args.wy
RES = args.RES

# ##############################################################################
log_name = "out_slurm_job.log"  # f"out_slurm_job_{wy}.log"
logging.basicConfig(filename=f'{log_name}', level=logging.INFO, format='%(asctime)s:%(message)s')
logging.info(f'============================START LOGGING {wy} ==================================')

# 1. First check the temp_outputs folder for the number of files
run_folder = f"/discover/nobackup/projects/coressd/Blender/Runs/{RES}/WY{wy}/temp_nc"
temp_nc_files = len(os.listdir(run_folder))
logging.info(f"Number of temp_nc_files: {temp_nc_files}")
print(f"Number of temp_nc_files: {temp_nc_files}")

# 2. Check the slurm jobs created for the water year
slurm_job_folder = f"{root_dir}/Github/slurm_jobs/{RES}/{wy}"
slurm_files = os.listdir(slurm_job_folder)
slurm_files = [f for f in slurm_files if f.endswith(".job")]
logging.info(f"Number of slurm files: {len(slurm_files)}")
print(f"Number of slurm files: {len(slurm_files)}")

# 3. Check Slurm job output files 
OUTDIR = f"{root_dir}/Github/slurm_jobs/{RES}/{wy}/.out"
files = os.listdir(OUTDIR)
logging.info(f"Number of files: {len(files)}")
print(f"Number of files: {len(files)}")
# Remove some of thes files that were renamed x1, x2_ etc during separate test or wshed runs
files = [f for f in files if not f.startswith("x")]  # Aside only when some files are renamed x1, x2_ etc
files = [f for f in files if not f.startswith("tes")]  # Aside only
files = [f for f in files if not f.startswith("ws")]  # Aside only
files = [f for f in files if not f.startswith("err")]  # Aside only
logging.info(f"Number of files after filtering files: {len(files)}")
print(f"Number of files after filtering files: {len(files)}")

# Check for errors in slurm job output file
err_file_list = []
for file in files:
    file_path = f"{OUTDIR}/{file}"
    with open(file_path, 'r') as file_obj:
        lines = file_obj.readlines()
        for line in lines:
            if "ERROR" in line or "error" in line:
                err_file_list.append(file)
err_file_list = list(set(err_file_list))
count_err_files = len(err_file_list)
print(f"Number Error files: {count_err_files}")
logging.info(f"Number Error files: {count_err_files}")
if count_err_files > 0:
    logging.info(f"Error files: {err_file_list}")
logging.info(f'============================END LOGGING {wy} ==================================\n')
