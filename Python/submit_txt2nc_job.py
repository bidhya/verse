#!usr/bin/env python

""" Generate and submit Slurm job for combining temporary text outputs to NETCDF file
    - Only for HPC.
    - Manually update this script everytime and make sure to change/udpate everything in this script here

    Usage: python ~/Github/verse/submit_txt2nc_job.py
    Best usage: Run this python from command line; will submit all slurm jobs
"""
import os

def mkdir_p(folder):
    '''make a (sub) directory (folder) if it doesn't exist'''
    if not os.path.exists(folder):
        os.mkdir(folder)


# Create (hidden) folder to save slurm outputs and julia logs
mkdir_p('slurm_jobs/.out')


def create_job(jobname='txt2nc.job', cores=1, memory='16gb', runtime='06:00:00', out_subfolder="NoahMP_CGF", region=None, var_name=None, var_idx=None):
    """ Since this is only designed for clipping, the default value of clip_flag is 1, ie, do the clipping"""
    job_file = os.path.join(os.getcwd(), f"slurm_jobs/{jobname}.job")

    # lizard_data = os.path.join(data_dir, lizard)
    # Create lizard directories
    # mkdir_p(lizard_data)
    with open(job_file, 'w') as fh:
        # TODO: Need if condition in many of this writelines code
        fh.writelines("#!/usr/bin/env bash\n\n")
        fh.writelines("#SBATCH --account=s2701\n")        
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
        # fh.writelines("#SBATCH --partition=howat\n")

        # first thing we do when the job starts is to "change directory to the place where the job was submitted from".
        fh.writelines("cd $SLURM_SUBMIT_DIR\n")
        fh.writelines("date; hostname; pwd\n")         # Add host, time, and directory name for later troubleshooting
        fh.writelines("echo ========================================================\n")
        # upto this point, it can be common to other jobs as well

        # Julia script specific inputs and parameters
        fh.writelines(f"echo Blender run for {out_subfolder} var_name = {var_name} end_idx = {var_idx}\n\n")
        # fh.writelines(f"cores={cores} #we can only run 30 jobs concurrently on Unity\n")
        # fh.writelines(f"log_name=.out/{jobname}.log\n")
        fh.writelines("\n")
        # Call the main julia script to run by this slurm script
        # fh.writelines(f"julia /discover/nobackup/byadav/Github/giuh/scripts/verse/Julia/call_Blender_v8.jl NA2 {start_idx} {end_idx}\n\n")
        # fh.writelines(f"julia /discover/nobackup/byadav/Github/giuh/scripts/verse/Julia/combine_txt2csv.jl NA2 {var_name} {var_idx}\n\n")        
        fh.writelines(f"julia /discover/nobackup/byadav/Github/verse/Julia/combine_txt2csv.jl NA2 {var_name} {var_idx}\n\n")
        fh.writelines("echo Finished Slurm job \n")
    # submit the job
    os.system(f"sbatch {job_file}")


def main():
    """ Main script to parameterize, generate, and dispatch slurm jobs
        combine_txt2csv.jl to COMBINE TEMPORARY TEXT FILES TO FINAL NETCDF FILE
    """
    create_job(jobname="txt2csv1", var_name="SWE", var_idx=1, cores=1, memory='16gb', runtime='05:00:00')
    create_job(jobname="txt2csv2", var_name="Gmelt", var_idx=2, cores=1, memory='16gb', runtime='05:00:00')
    create_job(jobname="txt2csv3", var_name="G", var_idx=3, cores=1, memory='16gb', runtime='05:00:00')
    create_job(jobname="txt2csv4", var_name="Precip", var_idx=4, cores=1, memory='16gb', runtime='05:00:00')
    create_job(jobname="txt2csv5", var_name="Us", var_idx=5, cores=1, memory='16gb', runtime='05:00:00')
    create_job(jobname="txt2csv6", var_name="Gpv", var_idx=6, cores=1, memory='16gb', runtime='05:00:00')
    create_job(jobname="txt2csv7", var_name="Gmeltpv", var_idx=7, cores=1, memory='16gb', runtime='05:00:00')
    create_job(jobname="txt2csv8", var_name="Upv", var_idx=8, cores=1, memory='16gb', runtime='05:00:00')
    create_job(jobname="txt2csv9", var_name="SWEpv", var_idx=9, cores=1, memory='16gb', runtime='05:00:00')


if __name__ == "__main__":
    """ Call the main function """
    main()
