# To execute a notebook using paperpill

import os
import papermill as pm
import argparse
parser = argparse.ArgumentParser(description='Run Jupyter notebook using papermill.')  # noqa
parser.add_argument('source', help='Fullpath to jupyter notebook', type=str)
parser.add_argument('--cores', help='Number of cores to use for multiprocessing', type=int, default=1)
parser.add_argument('--log_name', help='Name of Log file', type=str, default='hydrogen.log')

# parser.add_argument('--start_index', help='Region name to start running code', type=int)
args = parser.parse_args()
source = args.source
# src="/discover/nobackup/byadav/Github/giuh/11_a_Extract_SEUP_monthly.ipynb"
cores = args.cores  # -1 number of cores for pararallel processing
# log_name = args.log_name

# Generate full path to destination notebook
# destination = source.replace('giuh', 'giuh/paper_outputs')
destination = source.replace('verse/Python', 'Slurm/paper_outputs')

print(source)
print(destination)
# destination = "/discover/nobackup/byadav/Github/giuh/paper_outputs"
# source='/home/yadav.111/Github/giuh/notebooks/osu/01_Howat/Hydrogen/02_c_Clip_mosaics_by_AOI-param.ipynb'
# destination='/home/yadav.111/Github/giuh/notebooks/osu/01_Howat/Hydrogen/02_c_Clip_mosaics_by_AOI_param-out.ipynb'
# pm.execute_notebook('path/to/input.ipynb', 'path/to/output.ipynb', parameters=dict(alpha=0.6, ratio=0.1))
# pm.execute_notebook('/home/yadav.111/Github/papermill_test.ipynb', '/home/yadav.111/Github/test/test_out.ipynb', parameters=dict(x=500))

# pm.execute_notebook(source, destination, parameters=dict(start=start, end=end, cores=cores, region=region))
# pm.execute_notebook(source, destination)
pm.execute_notebook(source, destination, parameters=dict(cores=cores))
