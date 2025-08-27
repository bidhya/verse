""" Delete temporary files and folders
    Manually change here before running
    B. Yadav
    May 18, 2024

    Usage: python delete_files.py 2013 010
    Very low-overhead. So call directly from bash with python command.
"""

import os
import shutil
import argparse
import platform
from joblib import Parallel, delayed
ss = platform.system()
node = platform.node()
print(ss, node)

parser = argparse.ArgumentParser(description='Delete temp outputs for water year.')
parser.add_argument('water_year', help='Water Year for processing', type=str)
parser.add_argument('--RES', help='Grid size (Resolution) of SEUP Inputs', type=str, default="010")
args = parser.parse_args()
water_year = args.water_year
RES = args.RES

if 'L-JY0R5X3' in node and "Windows" in ss:
    root_dir = "D:/coressd"
elif 'L-JY0R5X3' in node and "Linux" in ss:
    root_dir = "/mnt/d/coressd"
elif "borg" in node and "Linux" in ss:
    root_dir = "/discover/nobackup/projects/coressd"
elif "discover" in node and "Linux" in ss:
    # if running on login node
    root_dir = "/discover/nobackup/projects/coressd"
elif "asc.ohio-state.edu" in node and "Linux" in ss:
    root_dir = "/fs/project/howat.4/yadav.111/coressd"
elif ".osc.edu" in node and "Linux" in ss:
    root_dir = "/fs/ess/PAS1785/coressd"
else:
    print("Unknow computer system. coressd folder NOT set")
    assert False
base_folder = f"{root_dir}/Blender"
print(f"base_folder: {base_folder}")

months = ["10", "11", "12", "01", "02", "03", "04", "05", "06", "07", "08", "09"]
yr_mon_list1 = [f"{int(water_year) - 1}{mon}" for mon in months[:3]]
yr_mon_list2 = [f"{water_year}{mon}" for mon in months[3:]]
yr_mon_list = yr_mon_list1 + yr_mon_list2
print("yr_mon_list:")
print(yr_mon_list)

variables = ["Snowf_tavg", "SWE_tavg", "Tair_f_tavg", "Qg_tavg"]  #"Rainf_tavg", 
# 1. Delete files from combined folder
# for month_subfolder in yr_mon_list:
for var in variables:
    print(var)
    combined_folder = f"{base_folder}/Inputs_{RES}/combined/{var}"
    for year_month in yr_mon_list:
        nc_file = f"{combined_folder}/{year_month}.nc"
        if os.path.exists(nc_file):
            print(f"Deleting {nc_file}")
            os.remove(nc_file)
            # shutil.rmtree(f"{combined_folder}/{year_month}")

# OLD: No more needed
# outputs_txt_folder = f"{base_folder}/Runs/x_WY2013/outputs_txt"
# outputs_txt_folder = f"{base_folder}/Runs/WY{water_year}/outputs_txt"
# pixels = os.listdir(outputs_txt_folder)
# print(len(pixels))
# # for pix in pixels[:10000]:
# #     # if output_str in os.listdir(f"{outputs_txt_folder}/{pix}"):
# #     shutil.rmtree(f"{outputs_txt_folder}/{pix}")
# def delete_folder(pix):
#     """ To delete in parallel use function.
#         os.remove(f'{download_folder}/{file}') : to delete on file
#     """
#     shutil.rmtree(f"{outputs_txt_folder}/{pix}")
# Parallel(n_jobs=-1)(delayed(delete_folder)(pix) for pix in pixels)
# print("After Deleting Files")
# pixels = os.listdir(outputs_txt_folder)
# if len(pixels) == 0:
#     print(f"Removing Empty directory: {outputs_txt_folder}")
#     os.rmdir(outputs_txt_folder)
# print(len(pixels))
