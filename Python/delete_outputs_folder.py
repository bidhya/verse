""" Delete outputs subfolder
    B. Yadav
    Jun 12, 2022
"""


import os
import shutil
from joblib import Parallel, delayed

import platform
ss = platform.system()
node = platform.node()
print(ss, node)

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
    assert(False)
base_folder = f"{root_dir}/Blender"


outputs_txt_folder = f"{base_folder}/Runs/x_WY2013/outputs_txt"
pixels = os.listdir(outputs_txt_folder)
print(len(pixels))

# for pix in pixels[:10000]:
#     # if output_str in os.listdir(f"{outputs_txt_folder}/{pix}"):
#     shutil.rmtree(f"{outputs_txt_folder}/{pix}")


def delete_folder(pix):
    """ To delete in parallel use function.
        os.remove(f'{download_folder}/{file}') : to delete on file
    """
    shutil.rmtree(f"{outputs_txt_folder}/{pix}")


Parallel(n_jobs=-1)(delayed(delete_folder)(pix) for pix in pixels)

print("After Deleting Files")
pixels = os.listdir(outputs_txt_folder)
print(len(pixels))
