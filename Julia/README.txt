Files here are listed in the order it should be run. Only the required files are described. Refer to the main readme on how this files should be run. Many of these scripts are called by slurm job.   

submit_slurm.jl 
===============
- create a slurm job (~ 500) and run on Discover.
- uses `call_Blender.jl` 
- call this directly on a bash terminal from Github folder on Discover.
- takes a few minutes to generate and also submit slurm jobs.
- check the created jobs in the following subfolder: ../Github/slurm_job/2015
- this subfolder will be created automatically by the script

call_Blender.jl
===================  
This reads the input files, selects a subset of pixel (by slice) can calls Julia optimization routine.  

Estimate.jl
===============  
- This is main JuMP/Ipopt optimzation code.  

combine_nc_files.jl  
===================
Run this script at the end of Blender run. This will combine the individual netcdf outputs for different slice into a PAN-America NetCDF file.  
