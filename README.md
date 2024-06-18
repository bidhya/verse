# Blender
Order of Script Execution for coressd project
=============================================   
## I: To Prepare data  
Basefolder for Data: ../coressd/Blender/  

1. paper_seup.py : Process SEUP data using papermill/jupyter notebook combination
    Use multiple cores
    Outputs: All daily seup data combined to one file
        ../Blender/Inputs_010/lis/WY2015
        - Snowf_tavg.nc
        - SWE_tavg.nc
        - Tair_f_tavg.nc
        - Qg_tavg.nc
        - SCF.nc

2. process_modis_cgf.py  
    Hard dependency with 1.  
    Extract NA (North America) scale daily Modis_CGF matching the LIS rasters
    Folders: ../coressd/Blender/Modis/
        - MOD44B : annual tree-fraction data processed by jupyter notebook
        - MOD10A1F:
            - clipped2015_010: clipped and matched to seup resolution

## II: To Run Blender Code
==========================  
- Locally for small area: 
    - julia verse/Julia/call_Blender_v18.jl test_WY2010 3005 3006 010  
- on HPC: the following julia script will generate slurm jobs for a water_year.  
    - julia submit_slurm.jl 2015 010 3  # example: for WY2015 with resolution=010 and adaptive slice = 3 rows.
To combine nc_files
- uses:  Github/Slurm_Blender/d_combine_out_nc_files.sh  
    - calls: julia verse/Julia/combine_nc_files.jl WY$water_year $RES
