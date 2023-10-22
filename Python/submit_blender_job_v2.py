#!usr/bin/env python

""" Generate and submit Slurm job to call_Blender_v10.jl  
    - Must run from command line; this script will thus submit the slurm jobs
    - Only for HPC.
    - Manually update this script everytime and make sure to change/udpate everything in this script here

    Usage: python ~/coressd/Github/verse/Python/submit_blender_job.py 2016  

    Approach
    ======== 
    Process a subset of pixels using start and end index. 
    Requires ~48 GB ram for 46 cores. discover allocates all cores to for single user/job!
    Can use the final_check run with more more memobry (120GB) and JULIA_NUM_THREADS=10 for converting text to nc files in parallel. -- still problematic    


"""
import os
import logging
import time
import platform
import argparse

parser = argparse.ArgumentParser(description='Create and Submit Blender jobs for Water Year.')
parser.add_argument('water_year', help='Water Year for processing', type=str)
parser.add_argument('start', help='Staring Row', type=int, default=0)
parser.add_argument('end', help='Ending Row', type=int, default=941)
parser.add_argument('step', help='Number of pixels to process in in job', type=int)

args = parser.parse_args()
water_year = args.water_year
start = args.start
end = args.end + 1  # for python we need one more for last index
step = args.step


def mkdir_p(folder):
    '''make a (sub) directory (folder) if it doesn't exist'''
    if not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)


# Create (hidden) folder to save slurm outputs and julia logs
""" Will be created relative where this script is submitted from"""
mkdir_p(f"slurm_jobs/{water_year}/.out")
os.chdir(f"slurm_jobs/{water_year}")


def create_job(hpc, jobname='test', cores=15, memory='50gb', runtime='12:00:00', out_subfolder="NoahMP_CGF", region=None, start_idx=1, end_idx=100):
    """ Generate and submit slurm job"""
    logging.info(f'jobname = {jobname}    out_subfolder = {out_subfolder}  start_idx = {start_idx} end_idx = {end_idx}')

    # for lizard in lizards:
    # job_file = os.path.join(job_directory, f"{lizard}.job")
    # job_file = f"/home/yadav.111/slurm_jobs/{jobname}.sh"
    # Or create job where the job is submitted from, including ouputs/logs
    #  relative to this location
    job_file = os.path.join(os.getcwd(), f"{jobname}.job")
    print(f"Jobfile: {job_file}")

    # lizard_data = os.path.join(data_dir, lizard)
    # Create lizard directories
    # mkdir_p(lizard_data)
    with open(job_file, 'w') as fh:
        # TODO: Need if condition in many of this writelines code
        fh.writelines("#!/usr/bin/env bash\n\n")
        if hpc == "discover":
            fh.writelines("#SBATCH --account=s2701\n")  # for DISCOVER
            fh.writelines('#SBATCH --constraint="[cas|sky]"\n')
            # fh.writelines("#SBATCH --qos=long\n")  # for jobs 12-24 hours.
        elif hpc == "osc":
            fh.writelines("#SBATCH --account=PAS1785\n")
        else:
            pass
            # fh.writelines(f"#SBATCH --partition=howat-ice\n")  # for Unity  eg: howat-ice
        fh.writelines(f"#SBATCH --job-name={jobname}.job\n")
        fh.writelines(f"#SBATCH --output=.out/{jobname}.out\n")  # Directory must exist; created above 
        # fh.writelines(f"#SBATCH --output={jobname}.out\n")
        # fh.writelines(f"#SBATCH --error=.out/{jobname}.err\n")
        fh.writelines(f"#SBATCH --time={runtime}\n")        
        # fh.writelines(f"#SBATCH --cpus-per-task={cores}\n")
        # fh.writelines(f"#SBATCH --nodes=1 --ntasks-per-node={cores}\n")
        fh.writelines(f"#SBATCH --nodes=1 --ntasks={cores}\n")
        fh.writelines(f"#SBATCH --mem={memory}\n")
        # fh.writelines("#SBATCH --qos=normal\n")
        fh.writelines("#SBATCH --mail-type=ALL\n")
        # fh.writelines("#SBATCH --mail-user=$USER@stanford.edu\n")
        fh.writelines("#SBATCH --mail-user=yadav.111@osu.edu\n")
        fh.writelines("\n")

        # first thing we do when the job starts is to "change directory to the place where the job was submitted from".
        # fh.writelines("cd $SLURM_SUBMIT_DIR\n")
        fh.writelines("echo $SLURM_SUBMIT_DIR\n")
        if hpc == "discover":
            fh.writelines("cd $LOCAL_TMPDIR\n")
        else:
            fh.writelines("cd $TMPDIR\n")
        fh.writelines("date; hostname; pwd\n")         # Add host, time, and directory name for later troubleshooting
        fh.writelines("echo ========================================================\n")
        # upto this point, it can be common to other jobs as well

        # Julia script specific inputs and parameters
        fh.writelines(f"echo Blender run for {out_subfolder} start_idx = {start_idx} end_idx = {end_idx}\n\n")
        # uncomment next line on thread if using threads in final run for creating final nc files from txt files 
        # fh.writelines("export OPENBLAS_NUM_THREADS=1\n")
        # fh.writelines(f"export JULIA_NUM_THREADS={cores}\n")  # 10
        # fh.writelines("export JULIA_NUM_THREADS=$SLURM_NTASKS\n")  # or $SLURM_JOB_CPUS_PER_NODE
        fh.writelines("sleep 5\n")  # for slurm error when scheduling on multi nodes
        # fh.writelines(f"cores={cores} #we can only run 30 jobs concurrently on Unity\n")
        # fh.writelines(f"log_name=.out/{jobname}.log\n")
        fh.writelines("\n")
        # Call the main julia script to run by this slurm script
        if hpc == "discover":
            fh.writelines("cp -r /discover/nobackup/projects/coressd/Github/verse .\n")
            # fh.writelines(f"julia verse/Julia/call_Blender_v12.jl {out_subfolder} {start_idx} {end_idx}\n\n")
            fh.writelines(f"julia verse/Julia/call_Blender_v14.jl {out_subfolder} {start_idx} {end_idx}\n\n")
            # fh.writelines(f"julia /discover/nobackup/byadav/Github/giuh/scripts/verse/Julia/call_Blender_v8.jl NA_temp {start_idx} {end_idx}\n\n")
            # fh.writelines(f"julia /discover/nobackup/projects/coressd/Github/verse/Julia/call_Blender_v10.jl {out_subfolder} {start_idx} {end_idx}\n\n")
        else:
            fh.writelines("cp -r ~/Github/verse .\n")
            # fh.writelines(f"julia verse/Julia/call_Blender_v12.jl {out_subfolder} {start_idx} {end_idx}\n\n")
            fh.writelines(f"julia verse/Julia/call_Blender_v14.jl {out_subfolder} {start_idx} {end_idx}\n\n")
            # fh.writelines(f"julia --threads $SLURM_NTASKS verse/Julia/call_Blender_v12.jl {out_subfolder} {start_idx} {end_idx}\n\n")  # for threads but give error
        # fh.writelines("cp *.log $SLURM_SUBMIT_DIR \n")  # if log files were saved on compute node
        fh.writelines("squeue --job $SLURM_JOBID \n")
        fh.writelines("echo List of files on TMPDIR\n")
        fh.writelines("echo ========================================================\n")
        fh.writelines("ls -ltrh\n")
        fh.writelines("ls logs|wc -l\n")
        fh.writelines("echo Finished Slurm job \n")
    # submit the job
    os.system(f"sbatch {job_file}")


def main():
    """ Main script to parameterize, generate, and dispatch slurm jobs
        Make this suite self-sufficient so everything can be controled from 
        here rather than submitting from the command line
        ie, this script needs to be updated everything the job is run 
        --cores 1 --memory 16gb --runtime 36:00:00 --jobname s2_download

        NB: 4 nodes (48 cores/192GB (187 max allowed)); 2 nodes (40 cores/192GB); and 4 nodes (26 cores/96GB)
        In the new workflow, ~2.5GB per core seems enough; check further

    """
    # 1. MAIN BLENDER RUN JOB
    import math
    import numpy as np
    node = platform.node()
    print(f"Node Name: {node}")
    runtime = "08:00:00"  # 12 "24:00:00"  # default (max for Discover)
    memory = "180gb"
    if "discover" in node:
        # on login node node name is discover, not borg
        hpc_name = "discover"
        cores = "45"  # so one core can be used to monitor run using srun/htop.
        # memory = "180gb"
    elif "asc.ohio-state.edu" in node:
        hpc_name = "unity"
        cores = "39"  # 24 cores with 96GB memory for old node
        # memory = "180gb"
        # runtime = '12:00:00'
    elif ".osc.edu" in node:
        hpc_name = "osc"
        cores = "40"
        memory = "175gb"
        # runtime = "11:00:00"  # "23:00:00"
    else:
        print("Unknow computer system. coressd folder NOT set")
        assert(False)
    logging.basicConfig(filename='job_submission.log', level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')
    logging.info('\n')
    logging.info('--------------------------------------Job Creation/Submission info----------------------------------------------')   
    # =======================================================================================================================================
    # Run for one or multiple regions based on integer index; 
    # the indexing is based on ascending order sorted glacier area (done inside the main script)
    # start at 0 in the beginning. Later can be be changed to whatever value depending on 
    # avilability of cores, runtime limitation etc.
    # start = 850  # 910  0 625000  
    # # step = 85000  # number of pixels to process; for DISCOVER: 46 cores X 12 job = 552 (<=560 allowed for long qos)
    # # job_count = math.ceil((1011329 + 1) / step)  # number of slurm jobs to process all this pixels 
    # # job_count = math.ceil((941 + 1) / step)  #(2336, 941)
    # # end = start + job_count * step  # 15 11  1011329 + 1
    # end = 910 + 1 # 941 + 1
    job_count = math.ceil((end - start) / step)
    logging.info(f'step = {step}. job_count = {job_count}. end = {end}')
    for i in np.arange(start, end, step):
        start_idx = i + 1
        end_idx = i + step
        print(start_idx, end_idx)
        jobname = f"{start_idx}_{end_idx}"
        create_job(hpc=hpc_name, jobname=jobname, out_subfolder=f"V14_WY{water_year}", start_idx=start_idx, end_idx=end_idx, cores=cores, memory=memory, runtime=runtime)  # 120gb  96 64 48gb
        # for Discover, usable node: Haswell=28; Skylake=36; Cascade=46
        # logging.info(f"jobname={jobname}, start_idx={start_idx}, end_idx={end_idx}, cores=36, memory=144gb, runtime=12:00:00 ")
        time.sleep(1)
    # # Final Sanity Check to see if all pixels are processed. Also for creating nc files in new version using thread, thus use more memory here
    # create_job(hpc=hpc_name, jobname="final_check", out_subfolder="WY2016", start_idx=1, end_idx=1011329, cores=46, memory='120gb', runtime='12:00:00')
    logging.info("Job submission complete\n")


if __name__ == "__main__":
    """ Call the main function """
    main()
