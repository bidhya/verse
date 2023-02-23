# blender

Feb 23, 2023 
============ 
After using data with missing values on SWE due to polar nights, we have problem in 
Estimate_v54.jl line 93. Will have similar errors for other if condtion checks.  
    # 1.3 Match up SWE and MSCF
    for i=1:nt
        if MSCF[i]==0  # Nov 04, 2022: "ERROR: LoadError: TypeError: non-boolean (Missing) used in boolean context" because some days there was no MODIS data
            WRFSWE[i] = 0
        end
    end


Feb 20, 2023 
============ 
Refactoring and reorganization codes into Julia, Python subfolders
Main Julia Scripts staged inside the Julia subfolder. These include:  
Estimate_vxx.jl : Main Blender Script that operate on one pixel at a time
call_Blender_vx.jl : To call Estimate_vxx.jl for a given netcdf files. Will also outpuf the final netcdf file is processing is finished for all pixels.  
combine_txt2csv.jl : [Optional] Combine the temporary text Blender outputs into a NetCDF file.
	May be necessary for large area where processing is not finished in one call.
	Still takes long time (~6 hour for North America)
    Thus combining each variable as a separate job

Nov 11, 2022
============ 
Order of Script Execution for coressd project
I: To Prepare data  
==================  
1. paper_seup.py : Process SEUP data using papermill/jupyter notebook combination
    Use multiple cores
    Outputs:
        combined/
        WY/
        WY_merged/

2. process_modis_cgf.py
    Extract NA (North America) scale daily Modis_CGF matching the SEUP rasters

3. merge_modis_seup.py [old name: clip_by_watershed.py]
    Part I: Concatenate MODIS along time dim then append to SEUP varaibles, making it ARD for Blender run
        and save nc file: ../NoahMP/WY_merged/2016_clip_noahmp_cgf.nc
    Part II: [Optional] Clip by watershed and 
    
TODO: Generate Analysis ready data for North America by merging SEUP and MODIS_CGF
    This can be intermediate step of clip_by_watershed script
    So that extraction part can be separated  
    Or even better 
    ============== 
    Move ../coressd/CGF_NDSI_Snow_Cover/NA inside the ..coressd/Blender/NoahMP/combined/ folder
                                                    or,  ..coressd/Blender/NoahMP/WY/ folder : 1 file for one variable for 1 Water Year
                                                    or, to ../WY_merged : one giant file with all variable and should currently be usable in Blender
    and organize monthly tiles and name it something called Modis_cgf 
    That way we have monthly north America tiles in separate variables
    Modify Blender Scirp the run on this kind of data structure
    This is manageable beacuse each monthly NA file is ~300 MB, so a total of 300 * 5 variables = 1.5 GB per month
    or ~17 GB per year
    but we will be reading one nc files for each variable separately
    
II: To Run
===========  
- Locally for small area: julia call_Blender_vx.jl output_folder start_index end_index
- on HPC 				: python submit_blender_job.py -> This will generate several slurm jobs that calls julia call_Blender_vx.jl output_folder start_index end_index
							Modify and update this Python script before running.
						: python submit_txt2nc_job.py -> Combine temporary text files into a netcdf file. Separately for each variable.   