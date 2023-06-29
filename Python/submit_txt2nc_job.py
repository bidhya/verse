#!usr/bin/env python

""" Generate and submit Slurm job to call_Blender_x.jl script
    - Only for HPC.
    - Manually update this script everytime and make sure to change/udpate everything in this script here

    Usage: python ~/Github/verse/submit_txt2nc_job.py 2016
    Best usage: Run this python from command line; this script will thus submit the slurm jobs

    Approach
    ======== 
    Process a subset of pixels using start and end index. 

"""
import os
import logging
import platform
import argparse

parser = argparse.ArgumentParser(description='Create and Submit Blender jobs for Water Year.')
parser.add_argument('water_year', help='Water Year for processing', type=str)
args = parser.parse_args()
water_year = args.water_year


def mkdir_p(folder):
    '''make a (sub) directory (folder) if it doesn't exist'''
    if not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)


# # Create (hidden) folder to save slurm outputs and julia logs
# mkdir_p('slurm_jobs/.out')
# os.chdir("slurm_jobs")
mkdir_p(f"slurm_jobs/{water_year}/.out")
os.chdir(f"slurm_jobs/{water_year}")


def create_job(hpc, jobname='txt2nc', cores=1, memory='16gb', runtime='07:00:00', out_subfolder="WY2016", var_name=None, var_idx=None):
    """ Generate and submit slurm job"""
    logging.info(f'jobname = {jobname}    out_subfolder = {out_subfolder}  var_name = {var_name} var_idx = {var_idx}')

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
        elif hpc == "osc":
            fh.writelines("#SBATCH --account=PAS1785\n")
        else:
            # pass
            fh.writelines(f"#SBATCH --partition=howat-ice\n")  # for Unity  eg: howat-ice

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
        fh.writelines(f"echo Blender run for {out_subfolder} var_name = {var_name} end_idx = {var_idx}\n\n")
        # fh.writelines("sleep 60\n")  # for slurm error when scheduling on multi nodes
        # fh.writelines(f"cores={cores} #we can only run 30 jobs concurrently on Unity\n")
        # fh.writelines(f"log_name=.out/{jobname}.log\n")
        fh.writelines("\n")
        # Call the main julia script to run by this slurm script
        if hpc == "discover":
            fh.writelines(f"julia /discover/nobackup/projects/coressd/Github/verse/Julia/combine_txt2nc.jl {out_subfolder} {var_name} {var_idx}\n\n")
        else:
            fh.writelines(f"julia ~/Github/verse/Julia/combine_txt2nc.jl {out_subfolder} {var_name} {var_idx}\n\n")     
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
    node = platform.node()
    print(f"Node Name: {node}")
    if "discover" in node:
        # on login node node name is discover, not borg
        hpc_name = "discover"
    elif "asc.ohio-state.edu" in node:
        hpc_name = "unity"
    elif ".osc.edu" in node:
        hpc_name = "osc"
    else:
        print("Unknow computer system. coressd folder NOT set")
        assert(False)
    logging.basicConfig(filename='txt2nc.log', level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')
    logging.info('\n')
    logging.info('--------------------------------------Job Creation/Submission info----------------------------------------------')   
    out_subfolder = f"WY{water_year}"  # "NA_2016" 
    create_job(hpc=hpc_name, jobname="txt2nc1", var_name="SWE", var_idx=1, out_subfolder=out_subfolder, cores=1, memory='16gb', runtime='07:30:00')
    create_job(hpc=hpc_name, jobname="txt2nc2", var_name="Gmelt", var_idx=2, out_subfolder=out_subfolder, cores=1, memory='16gb', runtime='07:30:00')
    create_job(hpc=hpc_name, jobname="txt2nc3", var_name="G", var_idx=3, out_subfolder=out_subfolder, cores=1, memory='16gb', runtime='07:30:00')
    create_job(hpc=hpc_name, jobname="txt2nc4", var_name="Precip", var_idx=4, out_subfolder=out_subfolder, cores=1, memory='16gb', runtime='07:30:00')
    create_job(hpc=hpc_name, jobname="txt2nc5", var_name="Us", var_idx=5, out_subfolder=out_subfolder, cores=1, memory='16gb', runtime='07:30:00')
    create_job(hpc=hpc_name, jobname="txt2nc6", var_name="Gpv", var_idx=6, out_subfolder=out_subfolder, cores=1, memory='16gb', runtime='07:30:00')
    create_job(hpc=hpc_name, jobname="txt2nc7", var_name="Gmeltpv", var_idx=7, out_subfolder=out_subfolder, cores=1, memory='16gb', runtime='07:30:00')
    create_job(hpc=hpc_name, jobname="txt2nc8", var_name="Upv", var_idx=8, out_subfolder=out_subfolder, cores=1, memory='16gb', runtime='07:30:00')
    create_job(hpc=hpc_name, jobname="txt2nc9", var_name="SWEpv", var_idx=9, out_subfolder=out_subfolder, cores=1, memory='16gb', runtime='07:30:00')


if __name__ == "__main__":
    """ Call the main function """
    main()
