# Updates to codebase in reverse chronological order
=============================================   
## Jun 22, 2024  
- Merged Jun branch to main.  
- Even bigger workflow overhaul for LIS 1km run.   
- Uses five input variables separately that is combined in Julia
- Qg_tavg replaces SWE_tavg as template raster for Blender run and combining
- MODSCAG renamed SCF with values between 0 and 1
- Adaptively generate Slurm job based on number of pixels.  


## Jun 02, 2024  
- Merged Apr07 branch to main.  
- Complete overhaul of processing workflow for using LIS 1km run.   
- Tried zarr. Was fast bud did not work for final NetCDF file.  

## Mar 20, 2024  
- Remove tar/copy from julia script and move to bash. This was crashing julia for files >3GB.   
- increase runtime in submit_slurm job    

## Feb 21, 2024  
- putting garbage collector inside Estimate_v59.jl, Estimate_v58.jl  
- and removing gc from call_Blender_v16.jl, call_Blender_v15.jl      
- to overcome many of the memory issues plaguing the run at large scale    

## Feb 20, 2024  
- creating call_Blender_v16.jl and Estimate_v59.jl    
- to get over bottleneck due to sharedarrays  
- will create text files inside compute node and convert these to NetCDF  
- this hybrid approach incorporates v12 and v15  
- that is, still run blender with horizontal slice and save intermediate nc files  

## Feb 11, 2024  
- using garbage collection 2% for for loop  
- this did reduce the memory footprint due to sharredarrays    
- current scripts are: call_blender_v15.jl with Estimate_v58.jl  

## Feb 07, 2024  
- updating documentation for generating various inputs  
- process_modis_cgf.py : moving mosaics intside temp subfolder.  
- compressing both mosaics and final output. 10x reduction in filesize.    


## Feb 02, 2024  
- update σWRFP using number of snowy days  
- move log tar can copy to end of file due to copying error, perhaps because larger file size in Milan Nodes  

## Jan 29, 2024  
- Diagnosed and fixed nan/inf error
- seems for only one pixel for WY2014 (row 888)
- caused by σWRFG becoming zero; which was caused by WRFG=0 (a rare but valid valid)  
- Fix: set minimum value of σWRFG = 25  

## Jan 18, 2024  
- merging Jan11_2024 branch to main  
- increased runtime for jobs  
- adding Milan (mil) to list of list of usable clusters on Discover    
- Blender error for one pixel in WY2014 from new updates by Jack remain  

## Jan 14, 2024  
- tar and gzip log folder on compute node  
- move the tar.gz file to ../Runs/WY20xx/logs/  
    - this is record of runs  
    - archive can be kept for long-term storage  
- tar from within julia code. But gzip will not work on windows. 
- there is also option to tar from slurm job but this is used anymore  
- seperate notebook (outside this verse git) to extract relevant optimization messages from tar and save to csv
- ../03_Blender/10_Analyze_output_log_files.ipynb  
    - run this ipynb with papermill from slurm job (../Github/Slurm_jobs/log_paper.job) passing WY as argument  

## Jan 11, 2024  
- Created Jan11_2023 branch  
- Changes to Incorporate
    - Manually create "log" folder either on compute node or on output_folder (for debugging purpose)  
    - call_Blender_v15.jl
        - control creation of log file from this script  
        - pass the location of this log folder from  
    - Estimate_v58.jl  
        - receive the "log" folder as function argument  

## Jan 10, 2024  
- Reverted to main branch from github with dec03, 2023 updates  
- due to many kinds of errors and long running times  
- working on discover. challenges on ther systems  

## Nov 21, 2023
Created JacksUpdateNov21 to incorporate changes from Jack.
New Files:
- Estimate_v58.jl     : changes from Jack
- call_Blender_v15.jl : calls Estimate_v58.jl

## Nov 20, 2023
Final commit to main before incorporating changes from Jack.
Final codes
    - submit_slurm.jl  
    - call_Blender_v14.jl  
    - combine_nc_files.jl  
Successfully ran using 35 to 40 rows and all available memory (>180 GB)  

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

