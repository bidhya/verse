# This is minimum set of commands to run Blender on Discover HPC  
cd /discover/nobackup/projects/coressd  
cd Github  
1. julia verse/Julia/submit_slurm.jl 2015 1  
    Wait for a few days to finish all the runs.  
cd ../Github/Slurm_Blender  
2. sbatch c_combine_out_nc_files.sh 2015  
    This will combine individual nc_files  into PAN American netcdf file  

[OPTIONAL]
# Following commands are for generating Input files. 
cd ../coressd/Github/Slurm_Blender
1. sbatch a_lis_process.sh 2015  
2. sbatch b_process_modis_cgf.sh 2015  
