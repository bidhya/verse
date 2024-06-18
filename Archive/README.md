# Blender
Order of Script Execution for coressd project
=============================================   
## I: To Prepare data  
Basefolder for Data: ../coressd/Blender/  

1. paper_seup.py : Process SEUP data using papermill/jupyter notebook combination
    Use multiple cores
    Outputs:
        ../Blender/Inputs/combined/.. : intermediate  
        WY/[varname]/2016.nc etc  # varname = one of SEUP variables. eg SWE_tavg  
        WY/2016_seup.nc  All daily seup data combined to one file    

2. process_modis_cgf.py  
    No dependency with 1.  
    Extract NA (North America) scale daily Modis_CGF matching the SEUP rasters
    Outputs: ../coressd/Blender/Modis/CGF_NDSI_Snow_Cover/
        - temp/..    # temporary mosaic of modis tiles over NA. can be safelydeleted
        - NA2016/..  # clipped and matched to seup resolution  

3. merge_modis_seup.py  
    - not required anymore for 1 km run
    Concatenate MODIS along time dim then append to SEUP varaibles. # This file is a final input for Blender run.  
    Outputs: 
        ../Blender/Inputs/WY_merged/2016_seup_modis.nc  

## II: To Run Blender Code
==========================  
- Locally for small area: julia call_Blender_v15.jl WY2016 start_index end_index
- on HPC: the following julia script will generate slurm jobs for a water_year.  
    - julia submit_slurm.jl 2015 40  # to run for WY2015 with each job processing horizontal slice of 40 pixels.  

    Following are old python based job-submission scripts. Designed for text-based outputs. Currently not used.  
    - python submit_blender_job.py -> This will generate several slurm jobs that calls julia call_Blender_vx.jl output_folder start_index end_index  
		Modify and update this Python script before running.
        ~ 48 GB for 45 cores; 11 separate jobs with 100,000 pixels per job  
	python submit_txt2nc_job.py -> Combine temporary text files into a netcdf file.  Separately for each variable.   
