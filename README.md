# Blender
The notes here are now only specific to LIS related 1km runs.  
010 refers to 0.01 degree (1km) resolution. Future runs can thus be readily adapted to any higher/lower resolution as required.  

## Project folder: /discover/nobackup/projects/coressd  
- everything is relative to this main folder on Discover  

### Subfolder
OSU: Modis snow cover data downloaded here.  
mwrzesie: LIS data generated by Melisa  
Github: codes related to Input generation, Blender runs, and post-processing.  
1. verse : scripts (both Python and Julia)
2. Slurm_Blender : slrum jobs related to generation of input files and post-processing of Blender runs.
3. slurm_jobs : will be created by blender run automatically as required.
4. Blender_notebooks: Jupyter notebooks. Mostly for prototyping and analysis of data. But two imporant ones for Blender workflow:
    - LIS/Extract_LIS_1km.ipynb : used to process LIS data. Uses papermill to run as a slurm job. Plan is to convert it to python script when everything is tested and finalized.
    - Modis/Forest/05_MOD44B_on_Discover.ipynb : use directly from Jupyterhub on discover. used to process annual tree-cover data. Since, data already processed no need to run again.
    - paper_outputs : not required to check per se. Output notebooks created by papermill is saved here.
    - rest of notebooks were created during various phases of project and mostly used for protyping of ideas and analysis.  

Blender: Blender related data processed by various scripts.
1. Inputs_010/lis: LIS and MODIS files staged by water year. Input files for Blender run.
2. Modis: MODIS intermediate files only
    - MOD44B : modis tree-cover reprojected and matched to LIS data  
    - MOD10A1F : modis snow-cover data corrected for tree and clipped to match LIS data  
Runs: Blender runs saved relative to resolution. 010 => 1km  
    - 010: 1km (1 degree) runs using LIS data staged by WaterYear  
    - 010/WY2015/temp_nc: temporary netcdf files for each slice or run
    - 010/WY2015/outputs: Final PAN Americal Outputs created by merging files from "temp_nc" folder  

Order of Script Execution for coressd project
=============================================   
## I: To Prepare data  
1. a_lis_process.sh : Process LIS data using papermill/jupyter notebook combination  
    Use multiple cores  
    Outputs: All daily LIS netcdf files saved separately as follows
        ../Blender/Inputs_010/lis/WY2015
        - Snowf_tavg.nc
        - SWE_tavg.nc
        - Tair_f_tavg.nc
        - Qg_tavg.nc
        - SCF.nc (will be created by next slurm job (b_process_modis_cgf.sh))  

2. b_process_modis_cgf.sh: Extract NA (North America) scale daily Modis_CGF matching the LIS rasters
- uses process_modis_cgf.py  
- uses ../OSU/MOD10A1F.061/MODIS_Proc/download_snow/.. for snow-cover data
- uses MOD44B for tree-cover correction  
- Hard dependency to use "Qg_tavg.nc" as template and mask    
- intermediate outputs: ../coressd/Blender/Modis/MOD10A1F/clipped2015_010: clipped and matched to seup resolution
- ../Blender/Inputs_010/lis/WY2015/SCF.nc : ie, same location as LIS outputs as above.  

## II: To Run Blender (Julia) Code
================================  
Install the following Julia packages:  
- add JuMP Ipopt Rasters NCDatasets CSV LoggingExtras Distributions  

1. On HPC (Discover): the following julia script will generate slurm jobs for a water_year.  
- cd to Github folder  
- julia verse/Julia/submit_slurm.jl 2015 2 010  # example: for WY2015 with resolution=010 and adaptive slice = 2 rows.
    - generate ~600 slurm jobs that calls the following julia script:
    - julia verse/Julia/call_Blender_v18.jl WY2010 1 17 010  : example script for Blender run for WY2010 for slice 1 to 17
    - slurm job folder: 
        - ../slurm_jobs/010/2015/
            - .out/ : slurm outputs saved here  
Blender Outputs: ../Blender/Runs/010/WY2015/temp_nc/ : These netcdf files for each slice will be combined by next script, after which this folder can be deleted. 

2. Post-processing (Required for continental map): Combine individual nc_files into PAN American netcdf file  
- uses:  Github/Slurm_Blender/d_combine_out_nc_files.sh  
- calls: julia verse/Julia/combine_nc_files.jl WY$water_year $RES

### Edge Cases
For test run:
- Use "test" prefix to WY. Example: test_WY2015 while calling call_Blender julia script.   
- folder and sub-folders created within script to save outputs.  
    - usage: julia verse/Julia/call_Blender_v19.jl "test/test_WY2015" 3204 3204 010 2

For 1 pixel run:
- Use "pixel" prefix to WY. Example: pixel_WY2015 while calling call_Blender julia script.   
- folder and sub-folders created within script to save outputs.  
    - usage: julia verse/Julia/call_Blender_v19.jl "pixel/pixel_WY2015" 4200 4200 010 2

For watershed run:
- need a text file of watershed bounding box (../coressd/Blender/coordinates/wshed.csv)  
- select watershed by passing the index (1-based in Julia) of watershed  
- usage: julia verse/Julia/call_Blender_v19.jl "wshed/Tuolumne/wshed_WY2015" 100 100 010 2 1
