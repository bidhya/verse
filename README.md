# Blender
Order of Script Execution for coressd project
=============================================   
## I: To Prepare data  
1. paper_seup.py : Process SEUP data using papermill/jupyter notebook combination
    Use multiple cores
    Outputs:
        combined/../..  
        WY/[varname]/2016.nc etc  
        WY/2016_seup.nc  OLD: WY_merged/2016_seup.nc  

2. process_modis_cgf.py  
    Extract NA (North America) scale daily Modis_CGF matching the SEUP rasters
    Outputs: /coressd/Blender/Modis/CGF_NDSI_Snow_Cover/
        - NA2016_mosaic/ [temp]
        - NA2016/


3. merge_modis_seup.py [old name: clip_by_watershed.py]  
    Part I: Concatenate MODIS along time dim then append to SEUP varaibles, making it ARD for Blender run
        Outputs: 
            ../NoahMP/WY_merged/2016_seup_modis.nc  #2016_clip_noahmp_cgf.nc
    Part II: [Optional] Clip by watershed  

## II: To Run
===========  
- Locally for small area: julia call_Blender_v10.jl output_folder start_index end_index
- on HPC
    - python submit_blender_job.py -> This will generate several slurm jobs that calls julia call_Blender_vx.jl output_folder start_index end_index  
		Modify and update this Python script before running.
        ~ 48 GB for 45 cores; 11 separate jobs with 100,000 pixels per job  
	python submit_txt2nc_job.py -> Combine temporary text files into a netcdf file.  Separately for each variable.   


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

    Rename Folders:
        NoahMP --> seup  
        WY_merged --> Inputs; even move outside the NoahMP folder

## Nov 07, 2023
v14.jl : Final to save nc files without intermediate text files. log files left on node (deleted at end of job).
submit_slurm.jl : to submit v14 using Julia and considering the number of pixels for runtime.  
New workflow is memory intensive so use all available RAM on nodes (~180 GB)  
Use all cores and memory in exclusive mode for slurm jobs.  

## Oct 14, 2023
v13.jl with SharedArray. working fine but not fully tested   
v14.jy : Created to make updates to V13 for Blender run

## Sep 06, 2023
v12.jl and Estimate_v56.jl frozen and merged to main. Good for runs with text files
v13.jl to prototype saving NetCDF files directly   

## Sep 06, 2023
Major refactoring going-on inside v12.jl and Estimate_v56.jl  
Separate outputs_txt and logs folders  
Modifiy code log to save only Info and higer messages  
Incorporate test: use "test" in substring for output folder to run test quickly  

## Aug 30, 2023
Freeze v11.jl and Estimate_v55.jl with cleanup of some comments  
Created v12.jl to make modifications    
Create Estimate_v56.jl for v12.jl

## Aug 07, 2023
Adding lazy=true when reading rasterstack because of updates to rasters and dimensional data. example,  
    A = RasterStack("$tmpdir/$water_year" * "_seup_modis.nc", lazy=true)  ## https://github.com/rafaqz/Rasters.jl/issues/449
    lazy=true may fix problem with newer raster and Julia 1.9.2

### Aug 03, 2023
synchronizec _v10 and _v11. So run _v11 everywhere including discover  
moved old _v11 to archive  

### Aug 02, 2023
change TMPDIR to LOCAL_TMPDIR for Discover  
log file no more working: try copy it at end of script from LOCAL_TMPDIR to ../.out/  
running _v10 on Discover (WY2010 to 2016) and _v11 for OSC  


### July 15, 2023
Merged all branches and deleted old branches  
Create call_blender_v11.jl for new prototyping. But still using v10 for all runs.  


### June 16, 2023
Call txt2nc python script for submitting slurm jobs from withing call_blender.jl file
Filxed Stripping issue when combining SEUP with MODIS data
    was caused by using mixing Rioxarray and Xarray in "merge_modis_seup.py" script  
    now using only Xarray  
Adapted for OSC  
Using same input file for all Slurm Systems  
    
### June 16, 2023
Pass water_year directly from slurm job for SEUP processing  
Process save Modis data by water year  
Appending older Tree Cover data for higher latitudes  
Running WY2013  

### May 23, 2023  
Uses forest fraction cover for MODIS CGF data  


### May 02, 2023  
Refractor: update file names; used newer xarray (warning but data still same), pandas etc.    
Convert NDSI to Snow cover fraction    

### May 02, 2023  
Refractor: update file names; used newer xarray (warning but data still same), pandas etc.    
Convert NDSI to Snow cover fraction    


### Mar 01, 2023  
Arctic Polar Nights issue fixed and working  
Updated python job submisison script to set HPC system automatically  
Run code on more than one node (use sleep time for nodes to be ready )
TODO: run text2nc in parallel  

### Feb 23, 2023  
After using data with missing values on SWE due to polar nights, we have problem in 
Estimate_v54.jl line 93. Will have similar errors for other if condtion checks.  
    # 1.3 Match up SWE and MSCF
    for i=1:nt
        if MSCF[i]==0  # Nov 04, 2022: "ERROR: LoadError: TypeError: non-boolean (Missing) used in boolean context" because some days there was no MODIS data
            WRFSWE[i] = 0
        end
    end


### Feb 20, 2023  
Refactoring and reorganization codes into Julia, Python subfolders
Main Julia Scripts staged inside the Julia subfolder. These include:  
Estimate_vxx.jl : Main Blender Script that operate on one pixel at a time
call_Blender_vx.jl : To call Estimate_vxx.jl for a given netcdf files. Will also outpuf the final netcdf file is processing is finished for all pixels.  
combine_txt2csv.jl : [Optional] Combine the temporary text Blender outputs into a NetCDF file.
	May be necessary for large area where processing is not finished in one call.
	Still takes long time (~6 hour for North America)
    Thus combining each variable as a separate job

