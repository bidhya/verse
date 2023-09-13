""" Check the intermediate text files of Blender Runs
    Delete pixel folders where run did not complete (need uncommmenting)
"""
import os
# from joblib import Parallel, delayed
import logging
# import shutil
import argparse

parser = argparse.ArgumentParser(description='Count Blender created text files for Water Year.')
parser.add_argument('water_year', help='Water Year for processing', type=str)
args = parser.parse_args()
water_year = args.water_year

# ##############################################################################
log_name = 'file_check.log'
logging.basicConfig(filename=f'{log_name}', level=logging.INFO, format='%(asctime)s:%(message)s')
logging.info('  ')
logging.info('-------------------------START LOGGING--------------------------------------')

base_folder = "/discover/nobackup/projects/coressd"
tmp_txtDir = f"{base_folder}/Blender/Runs/WY{water_year}/outputs_txt"    # Outputs are saved here
pixels = os.listdir(tmp_txtDir)
pixels = [pix for pix in pixels if pix.startswith("Pix")]
print(len(pixels))
logging.info(f"Number pixel folder: {len(pixels)}")

# Uncomment this section to delete pixels folders where run did not complete
# But not required anymore because similar checks already happening in "call_Blender.jl" script
# for pix in pixels:
#     file_count = os.listdir(f"{tmp_txtDir}/{pix}")
#     if len(file_count) != 2:  # ie, has only one file. Most likely this file is 'log.txt'
#         logging.info(f"Deleting: {tmp_txtDir}/{pix}")
#         # TODO: uncomment next line to delete the folder
#         # shutil.rmtree(f"{tmp_txtDir}/{pix}")

missing_count = 0
processed_count = 0
for pix in pixels:
    # file_count = os.listdir(f"{tmp_txtDir}/{pix}")
    fname = f"{tmp_txtDir}/{pix}/out_vars.txt"
    if not os.path.exists(fname):
        logging.info(f"Check: {tmp_txtDir}/{pix}")
        missing_count += 1
    else:
        processed_count += 1
        # to check if some files are only partially processed (ie, runtime terminated pre-maturely)
        # This kind of error message from Julia Blender Run:  nested task error: DimensionMismatch: tried to assign 124-element array to 1×1×366 destination
        with open(fname) as fi:
            # Runtime ~6 hours with this condition 
            if len(fi.readlines()) < 364:
                logging.info(pix)
                # TODO: move the this pixel subfolder to temporary location for further examination
                # will also help in Blender run because the pixel will be processed and hence the text file can be combined to nc file  

logging.info(f'Processed_count = {processed_count}')
logging.info(f'Missing_count = {missing_count}')
