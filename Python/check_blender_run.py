""" Check the intermediate text files of Blender Runs
    Delete pixel folders where run did not complete (need uncommmenting)
"""
import os
# from joblib import Parallel, delayed
import logging
# import shutil

# ##############################################################################
log_name = 'file_check.log'
logging.basicConfig(filename=f'{log_name}', level=logging.INFO, format='%(asctime)s:%(message)s')
logging.info('  ')
logging.info('-------------------------START LOGGING--------------------------------------')

base_folder = "/discover/nobackup/projects/coressd"
tmp_txtDir = f"{base_folder}/Blender/Runs/WY2016/outputs_txt"    # Outputs are saved here
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
    if not os.path.exists(f"{tmp_txtDir}/{pix}/out_vars.txt"):
        logging.info(f"Check: {tmp_txtDir}/{pix}")
        missing_count += 1
    else:
        processed_count += 1

logging.info(f'Processed_count = {processed_count}')
logging.info(f'Missing_count = {missing_count}')
