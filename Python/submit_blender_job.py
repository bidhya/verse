#!usr/bin/env python

""" Generate and submit Slurm job to call_Blender_x.jl script
    - Only for HPC.
    - Manually update this script everytime and make sure to change/udpate everything in this script here

    Usage: python ~/Github/verse/submit_txt2nc_job.py
    Best usage: Run this python from command line; this script will thus submit the slurm jobs

    Approach
    ======== 
    Process a subset of pixels using start and end index. 

"""
import os
import logging
import time


def mkdir_p(folder):
    '''make a (sub) directory (folder) if it doesn't exist'''
    if not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)


# Create (hidden) folder to save slurm outputs and julia logs
mkdir_p('slurm_jobs/.out')
os.chdir("slurm_jobs")


def create_job(hpc, jobname='test', cores=15, memory='60gb', runtime='12:00:00', out_subfolder="NoahMP_CGF", region=None, start_idx=1, end_idx=100):
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
        elif hpc=="osc":
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
        fh.writelines(f"#SBATCH --nodes=1 --ntasks-per-node={cores}\n")
        fh.writelines(f"#SBATCH --mem={memory}\n")
        # fh.writelines("#SBATCH --qos=normal\n")
        fh.writelines("#SBATCH --mail-type=ALL\n")
        # fh.writelines("#SBATCH --mail-user=$USER@stanford.edu\n")
        fh.writelines("#SBATCH --mail-user=yadav.111@osu.edu\n")
        fh.writelines("\n")

        # first thing we do when the job starts is to "change directory to the place where the job was submitted from".
        fh.writelines("cd $SLURM_SUBMIT_DIR\n")
        fh.writelines("date; hostname; pwd\n")         # Add host, time, and directory name for later troubleshooting
        fh.writelines("echo ========================================================\n")
        # upto this point, it can be common to other jobs as well

        # Julia script specific inputs and parameters
        fh.writelines(f"echo Blender run for {out_subfolder} start_idx = {start_idx} end_idx = {end_idx}\n\n")
        # fh.writelines(f"cores={cores} #we can only run 30 jobs concurrently on Unity\n")
        # fh.writelines(f"log_name=.out/{jobname}.log\n")
        fh.writelines("\n")
        # Call the main julia script to run by this slurm script
        # fh.writelines(f"julia /discover/nobackup/byadav/Github/giuh/scripts/verse/Julia/call_Blender_v8.jl NA_temp {start_idx} {end_idx}\n\n")        
        # fh.writelines(f"julia /discover/nobackup/byadav/Github/verse/call_Blender_v9.jl NA_temp {start_idx} {end_idx}\n\n")        
        fh.writelines(f"julia ~/Github/verse/Julia/call_Blender_v9.jl {out_subfolder} {start_idx} {end_idx}\n\n")        
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
    import numpy as np
    logging.basicConfig(filename='job_submission.log', level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')
    logging.info('\n')
    logging.info('--------------------------------------Job Creation/Submission info----------------------------------------------')   
    # =======================================================================================================================================
    # Run for one or multiple regions based on integer index; 
    # the indexing is based on ascending order sorted glacier area (done inside the main script)

    # start at 0 in the beginning. Later can be be changed to whatever value depending on 
    # avilability of cores, runtime limitation etc.
    start = 920000  # 625000  # 375000  # 250000  # 0
    step = 30000  # 20000 number of pixels to process
    end = start + 4 * step  # 500000  # 1011329 + 1
    for i in np.arange(start, end, step):
        start_idx = i + 1
        end_idx = i + step
        print(start_idx, end_idx)
        jobname = f"{start_idx}_{end_idx}"
        # hpc = unity osc discover
        create_job(hpc="unity", jobname=jobname, out_subfolder="NA4", start_idx=start_idx, end_idx=end_idx, cores=40, memory='56gb', runtime='15:15:00')
        logging.info(f"jobname={jobname}, start_idx={start_idx}, end_idx={end_idx}, cores=40, memory=56gb, runtime=12:00:00 ")
        time.sleep(2)

    # Final Sanity Check to see if all pixels are processed
    # create_job(jobname="final_check.job", start_idx=1, end_idx=1011329, cores=24, memory='56gb', runtime='12:00:00')
    # create_job(jobname="1_600000", start_idx=1, end_idx=600000, cores=34, memory='56gb', runtime='12:00:00')
    logging.info("Job submission complete\n")


if __name__ == "__main__":
    """ Call the main function """
    main()
