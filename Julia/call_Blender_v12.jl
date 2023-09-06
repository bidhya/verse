"""
USAGE: Hardcoded with requirments for output_folder_YEAR, start_idx, and end_idx for running on Discover
- julia call_Blender_v12.jl test2016 1 100

SWE estimation using Blender algorithm (Prof. Mike Durand).
This is a helper script to organize inputs for calling the Blender function (Estimate_v3.jl).
Therefore, manually update input/output directories here for seperate runs (say NoahMP vs WRF etc)
Currently, only output_subfolder can optionally be passed as an argument.
Input/Output are currently saved in same main folder, but can easily be saved at separate locations [dig down below].

The script is currently capable of running on my following platforms.
- Windows machine
- Ubuntu (WSL2-based)
- OSC HPC
- Ohio-State Unity HPC 
Julia Versions
- 1.8.x (1, 2, and 3)
- 1.9.3
Minimum extra Julia Packages required
- JuMP, Ipopt, Rasters, NCDatasets 

Approach
========
This script will call Estimate_v56.jl for each script at a time, and save temporary outputs to text file. Thus can readily be parallelzed.
Threads was stright forward but did not work here due to known limitation of Ipopt and threads module.
Working on alternative parallelization scheme using DISTRIBUTED module
Finally, all text output files are assembled into a nc file. The script can thus run on parts of pixels at different times, and finally combined into one nc file.
==============================================================================================
Dec 03, 2022 : Trying multinodes with ClusterManagers.jl
Sep 04, 2023 : text and log files saved in separate folders; created call_Blender_v12.jl uses Estimate_v56.jl
"""
arg_len = length(ARGS)
out_subfolder = ARGS[1]  # WY2016. output subfolder relative to input files; temp_text and nc_outputs saved here
water_year = out_subfolder[end-3:end] #last 4 chars are assumed year, else error. out_subfolder[3:end]
start_idx = ARGS[2]  # this is string
end_idx = ARGS[3]
start_idx = parse(Int64, start_idx)
end_idx = parse(Int64, end_idx)
if occursin("test", out_subfolder)
    # if test substring is part of output subfolder then do the test run on subset of pixels
    test_run = true
else
    test_run = false
end

using Logging, LoggingExtras
using Distributed  # otherwise everywhere macro won't work
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Setup inputs and output directories and file locations
# select the root directory (this will be different on windows, linux or Mac) 
host_machine = gethostname()
# if occursin("STAFF-BY-M", host_machine)
if occursin("L-JY0R5X3", host_machine)
    if Sys.iswindows()
        root_dir = "D:"  #for windows machine
    elseif Sys.islinux()
        root_dir = "/mnt/d"  #for Ubuntu (WSL2 on windows machine
    # else
    #     @info("Unknown system on windows, manually add root directory before proceeding. Exiting code")
    #     exit(1)    
    end
    base_folder = "$root_dir/coressd/Blender"
    # Get info on system machine. We use this info to farm number of processors. 
    system_machine = "Windows"  # a bit misnomer flag because this flag also works for WLS linux. Better flag could be local vs hpc/remote execution
elseif occursin("borg", host_machine)  # TODO: discover
    root_dir = "/discover/nobackup/projects/coressd"
    base_folder = "$root_dir/Blender"
    tmpdir =  ENV["LOCAL_TMPDIR"]  #tempdir()
    system_machine = "Slurm"
elseif occursin(".osc.edu", host_machine)
    root_dir = "/fs/ess/PAS1785/coressd"  # "/fs/scratch/PAS1785/coressd"
    base_folder = "$root_dir/Blender"
    tmpdir =  ENV["TMPDIR"]  #tempdir()
    system_machine = "Slurm"
elseif occursin("asc.ohio-state.edu", host_machine)  # .unity
    root_dir = "/fs/project/howat.4/yadav.111/coressd"  # homedir()  #  Unity
    # base_folder = "/home/yadav.111/Github/Blender"  # old
    base_folder = "$root_dir/Blender"  # "$root_dir/Github/coressd/Blender"
    tmpdir =  ENV["TMPDIR"]  #tempdir()
    system_machine = "Slurm"
    # @info("SLURM_SUBMIT_DIR: $(ENV["SLURM_SUBMIT_DIR"])")
else
    @info("Unknown computer, manually add root directory before proceeding. Exiting code")  # will output directly to console, ie like print statement
    exit(1)  # if error, comment this line, add root and base_folder below and run again
end
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Sep 05, 2023: Removing SlurmClusterManager prototyping 
if system_machine == "Windows" # addproc function works for both windows and WSL linux. Hence, this code is more for local vs remote machine
    log_filename = string(start_idx, "_", end_idx, ".log")  # on HPC created inside computer local node, so move to outside at end of job
    addprocs()  ## works on laptop, but on cluster shows all cores for node; not only the number of cpus requested
elseif system_machine == "Slurm"
    log_filename = string(ENV["SLURM_SUBMIT_DIR"], "/",start_idx, "_", end_idx, ".log")  # on HPC created inside computer local node, so move to outside at end of job
    # cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"])  # ERROR: LoadError: KeyError: key "SLURM_CPUS_PER_TASK" not found [from Princenton website]
    # ntasks = parse(Int, ENV["SLURM_NTASKS"])  # gave 1, when cpus-per-task was passed by slurm; not defined when just nodes=1 given.  
    # cores = parse(Int, ENV["SLURM_JOB_CPUS_PER_NODE"])  # if one node provided; else we have parse as follows
    subs = Dict("x"=>"*", "(" => "", ")" => "");
    cores = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))  # for 1 or more nodes (generic)
    # println("SLURM_NTASKS: $ntasks")
    # println("SLURM_JOB_CPUS_PER_NODE: $ntasks")
    @info("SLURM Cores: $cores")
    # addprocs()  # Not working. using all cores even though not allocated for use
    addprocs(cores)  # ; exeflags="--project"subtract one becuase master already has one; but seems to work with higher number as well
    # using SlurmClusterManager  # to pick all nodes and cores allocated in slurm job
    # # With slurmcluster manager, hopefully next few lines of selecting cores is not necessary  
    # sleep(10)  # 60+cores To prevent scheduling job on more than 1 node on Discover Slurm cluster  
    # addprocs(SlurmManager())  # to use all available nodes and cores automatically. comment this line and uncomment one above this to match _v8.jl
else
    @info("Must be windows or linux system, aborting")
    exit() # <- you can provide exit code if you want
end
# addprocs()
logger = FormatLogger(open(log_filename, "w"), 0) do io, args
    # Write the module, level and message only
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end
info_only_logger = MinLevelLogger(logger, Logging.Info);  # Logging.Error
global_logger(info_only_logger)  # Set the global logger to logger; else logs doing directly to console

@info("base_folder : $base_folder")
@info("Host computer machine: $host_machine")

# println("Number of procs: $(nprocs())")  
@info("Number of processes: $(nprocs())")
@info("Number of workers: $(nworkers())")
# Everywhere should come after cores are already added
@everywhere using Rasters, NCDatasets
@everywhere include("Estimate_v56.jl")  #https://docs.julialang.org/en/v1/manual/code-loading/; evaluated in global scope

# base_folder = "$root_dir/Github/Blender"
# DataDir= "$base_folder/nc_files"
DataDir = "$base_folder/NoahMP"  # must exist
# # use only for HPC to copy input file here (hopefully for faster i/o operation)
# if occursin("borg", host_machine)
#     tmpdir =  ENV["LOCAL_TMPDIR"]  #tempdir()  
# else
#     # TODO: on discover if node does not start with "borg" the code can still jump here!
#     tmpdir =  ENV["TMPDIR"]  #tempdir()

# Folder for saving outputs of run. out_subfolder can be passed as ARGS. Folder/subfolders will be created if non-existent
if occursin(".osc.edu", host_machine)
    # out_folder = "/fs/ess/PAS1785/coressd/Blender/Runs/$out_subfolder" #  "$base_folder/Runs/$out_subfolder"  # "$DataDir/Runs/$out_subfolder"
    out_folder = "/fs/scratch/PAS1785/coressd/Blender/Runs/$out_subfolder" #  "$base_folder/Runs/$out_subfolder"  # "$DataDir/Runs/$out_subfolder"
# elseif occursin("asc.ohio-state.edu", host_machine)  # .unity
#     out_folder = "$tmpdir/$out_subfolder"  # "$DataDir/Runs/$out_subfolder"
else
    out_folder = "$base_folder/Runs/$out_subfolder"  # "$DataDir/Runs/$out_subfolder"
end
@info("Output_folder : $out_folder")

# Make a folder insise HPC node because we want to copy existing files there
exp_dir = "$out_folder/outputs_txt"    # tmp_txtDir(old name) To save text outputs for each pixel
logDir = "$out_folder/logs"    # Save logs in separate folder
mkpath(exp_dir)  # mkdir
mkpath(logDir)  # mkdir
# cp("$base_folder/Runs/$out_subfolder/outputs_txt", "$tmp_txtDir/$water_year")  # copy to local machine; error if running the first time as this dir would not exist

nc_outDir = "$out_folder/outputs"         # To convert text outputs to netcdf file
# OutDir = "$DataDir/$out_subfolder/outputs_txt"  # To save text outputs for each pixel
# nc_exp_dir = "$DataDir/$out_subfolder/outputs_nc"  # To convert text outputs to netcdf file

# # 2. Read the Input netCDF file
# A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_modscag.nc")  #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
# Following check are for prototyping only when running code locally, because I do not yet have NorthAmerica netcdf file
if occursin("L-JY0R5X3", host_machine) || test_run # || testing = True; use this for testing as well
    @info("Test Run only.")
    # A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_cgf.nc", lazy=true)  # Error now (Aug 15, 2023). maybe it is old file/format
    A = RasterStack("$(DataDir)/WY_merged/$(water_year)_seup_modis.nc", lazy=true)
    subset = A[X(Between(-120.5, -120)), Y(Between(60, 60.5))]  # use small chip for prototyping
    subset_fname = "$(DataDir)/WY_merged/subset_$(water_year)_seup_modis.nc"
    isfile(subset_fname) || write(subset_fname, subset)  # save. we need for post-processing analysis of results
    A = RasterStack(subset_fname, lazy=true)  # read again (why). for consistency
# elseif occursin("borg", host_machine)  # TODO: discover
#     # A = RasterStack("$DataDir/WY_merged/2013_seup_modis.nc")  # 2016_noahmp_cgf 2016_clip_noahmp_cgf #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
#     A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
# elseif occursin(".osc.edu", host_machine)
#     # A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_cgf.nc")
#     A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
# elseif occursin("asc.ohio-state.edu", host_machine)
#     A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
else
    # A = RasterStack("$DataDir/WY_merged/2016_seup_modis.nc")
    # A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
    cp("$DataDir/WY_merged/$water_year" * "_seup_modis.nc", "$tmpdir/$water_year" * "_seup_modis.nc")  # copy to local machine
    A = RasterStack("$tmpdir/$water_year" * "_seup_modis.nc", lazy=true)  ## https://github.com/rafaqz/Rasters.jl/issues/449
    # exit(1)  # println("Exiting code. Manually set A (rasterstack) around line 188")
end

# # For Testing
# A = RasterStack("$(DataDir)/WY_merged/$(water_year)_seup_modis.nc", lazy=true)
# subset = A[X(Between(-120.5, -120)), Y(Between(60, 60.5))]  # use small chip for prototyping
# subset_fname = "$(DataDir)/WY_merged/subset_$(water_year)_seup_modis.nc"
# isfile(subset_fname) || write(subset_fname, subset)  # save. we need for post-processing analysis of results
# A = RasterStack(subset_fname, lazy=true)  # read again (why). for consistency


# A = RasterStack("$DataDir/WY_merged/2016_clip3.nc")  # for NoahMP with CGF MODIS
# Subset only the required variables because the nc file can have extraneous vars that cause problem with julia
A = A[(:Snowf_tavg, :SWE_tavg, :Tair_f_tavg, :Qg_tavg, :MODSCAG)];  # to remove spatial ref that cause problem during subsetting

# Select a subset for prototyping
# A[X(Between(2, 4)), Y(Between(50, 60)), Ti(5)]
# A = A[X(4:6), Y=4:6] #, Ti(1:10)  # For prototyping only, select a smaller subset of pixels
# A = A[X(38:40), Y=(38:40)]  # For prototyping only, select a smaller subset of pixels
# println("A: ", A)

# Now using cartesian index for iteration

start_time = time_ns()
# Extract cartesian index for non-missing (ie, all) data using one-day of data 
valid_pix_ind = findall(!ismissing, A["SWE_tavg"][Ti=1])
valid_pix_count = length(valid_pix_ind)
@info("Total valid pixel count = $(valid_pix_count)")
if start_idx > valid_pix_count
    # return nothing
    exit(0)
    @info("Exiting without Run because start_idx is greater than valid_pix_idx", start_idx, valid_pix_count)
end
if end_idx > valid_pix_count
    # if the passed index is more than valid number of pixels
    # in julia end index cannot be larger than this
    end_idx = valid_pix_count
end
# TODO: Select a-priori what pixels are already processed. We want this because script take hours to check files 
# that are already processed in the section below.

# ind = valid_pix_ind  # to process all pixels in one go
ind = valid_pix_ind[start_idx:end_idx]
# TODO (May 26, 2026): check the files that are already processed, so processing of those can be skipped

@info("Starting with loop.")
@info("==========================================================")
@info("Processing $(length(ind)) of $(length(valid_pix_ind))")
# @info("Processing ", length(ind), " out of ", length(valid_pix_ind))
# exit(0)
""" General info about pixel counts
xind: 2336 , yind: 941 tind:364
Non-missing pixel count = 1011329
"""
# @sync makes the code wait for all processes to finish their part of the computation before moving on from the loop
# ie, without @sync, the on of the processors may move to next while loop is still running, creating error and crash whole script prematurely 
# Threads.@threads for ij in ind
# @distributed for ij in ind
@sync @distributed for ij in ind
# for ij in ind
    i = ij[1]
    j = ij[2]
    # @info("ij = ", ij)
    # exp_dir = string("$tmp_txtDir/", "Pix_", i, "_", j)  # full path to folder for saving (temporary) text files   
    # if !isdir(exp_dir)  # process only if the pixel is not already processed
    # if !ispath(exp_dir * "/SWEpv.txt") # ispath is generic for both file and folder
    # if !isfile(exp_dir * "/SWEpv.txt")  # This is better check: process the pixel only if the last file exported by optimzer (SWEpv.txt) does not yet exist
    if !isfile("$(exp_dir)/Pix_$(i)_$(j).txt")  # may still be problem if scrip ended prematurely without writing the full file
        # extract each of the input variables separately (trying to match of processing was done in prior version with text inputs)
        # but this approach is also useful if we decide to save each input netcdf file separately
        WRFSWE = A["SWE_tavg"][X=i, Y=j].data  # Here "data" is an AbstractArray.
        WRFP = A["Snowf_tavg"][X=i, Y=j].data
        WRFG = A["Qg_tavg"][X=i, Y=j].data
        AirT = A["Tair_f_tavg"][X=i, Y=j].data
        MSCF = A["MODSCAG"][X=i, Y=j].data;
        # @info(exp_dir * "/SWEpv.txt")
        # Call blender for the pixel. This is the only required line
        # rest of the codes are for houskeeking, preprocssing, post-processing
        blender(out_folder, i, j, WRFSWE, WRFP, WRFG, MSCF, AirT)
        GC.gc()
        # cp(exp_dir, "$base_folder/Runs/$out_subfolder/outputs_txt", force=True)
        # cp(exp_dir, "$base_folder/Runs/$out_subfolder/outputs_txt", force=True)
        # count +=1
        # @info("After calling BLENDER function")
        # here exp_dir = export directory for holding outputs text and log files (same as older version of code)
    end
end
# without sync above one of the processor to the next step (combining text files to netcdf) which will cause error
end_time = time_ns()
running_time = (end_time - start_time)/1e9/3600
@info("Blender Running Time = $(round(running_time, digits=3)) hours")
# exit(0)  # exit here to avoid running 2nd part (text2nc)

# Step 4. PostProcessing: Combine text files into a grid and save as netcdf
function text2nc(var, idx, outRaster)
    """Convert a text file to netcdf
    var: the name of variable (example SWE, Gmelt, G etc.)
    idx : the index where the variable of exist in the file [1 through 9]
    TODO: move this fuction inside the Estimate_v53.jl script
    TODO: seems like nc file produced here cannot by plotted using xarray/hvplot
    """
    # outRaster = copy(A[:MODSCAG])  # TODO: check for NoahMP; use other variable because it creates MODIS attrs
    # outRaster = copy(A[:SWE_tavg])
    # Initial resulting array
    outRaster[:,:,:] .= missing
    # update the name of variable
    outRaster = rebuild(outRaster; name=var)  #:SWEhat
    # @sync @distributed for pix in pixels  # worked
    Threads.@threads for pix in pixels  # using this nested thread too seems to help processing faster   
    # for pix in pixels  # this worked
        pix_xy = split(pix, ".txt")[1]
        pix_xy = split(pix_xy, "_")
        x = parse(Int16, pix_xy[2])
        y = parse(Int16, pix_xy[3])
        out_vars = readdlm("$(exp_dir)/$(pix)") # readdlm("$(exp_dir)/Pix_$(x)_$(y).txt")
        # out_vars = readdlm("$exp_dir/$pix/out_vars.txt")  # $(exp_dir)/Pix_$(i)_$(j).txt   read whole combined file for that pixel
        # out_vars = readdlm("$tmpdir/outputs_txt/$pix/out_vars.txt")  # read whole combined file for that pixel
        # outRaster[X=x, Y=y] = readdlm("$exp_dir/$pix/$var.txt");  # Older workflow
        outRaster[X=x, Y=y] = out_vars[:,idx]
    end
    # Save to nc file; 
    mkpath(nc_outDir)
    write("$nc_outDir/$var.nc", outRaster)  # Aug 08, 2023: Error perhaps due to updates to raster/ncdataset/etc. (error tested on windows and wsl2)
end

# Count if all pixels are processed before calling the following section for converting text files to netcdf file
sleep(1)

pixels = readdir(exp_dir)  # Danger: Error if we have outputs from prior runs that do not match current A dimensions
pixels = [pix for pix in pixels if startswith(pix, "Pix")];
# if length(pixels) == valid_pix_count && system_machine == "Slurm"
#     run(`python $root_dir/Github/verse/Python/submit_txt2nc_job.py $water_year`)
#     # run(`python $tmpdir/verse/Python/submit_txt2nc_job.py $water_year`)
#     println("Submitted python script for converting text to nc file")
# end

if length(pixels) == valid_pix_count #&& system_machine == "Windows"
    if system_machine == "Slurm"
        # copy text files to node; but this seems to increase runtime on OSC; so comment next two lines to leave txt_files on original location  
        t1 = time_ns()
        cp(exp_dir, "$tmpdir/outputs_txt")  
        exp_dir = "$tmpdir/outputs_txt"
        t2 = time_ns()
        running_time = (t2 - t1)/1e9/3600
        @info("Total to copy text time files to tmpdir (hours) = $running_time")
    end
    @info("Creating OUTPUT NETCDF FILES")
    outRaster = copy(A[:SWE_tavg])
    var_idx_tuple = (("SWE", 1), ("Gmelt", 2), ("G", 3), ("Precip", 4), ("Us", 5), ("Gpv", 6), ("Gmeltpv", 7), ("Upv", 8), ("SWEpv", 9))
    for var_idx in var_idx_tuple  # Threads.@threads 
        var_name = var_idx[1]
        var_idx = var_idx[2]
        text2nc(var_name, var_idx, outRaster);
    end
    # # Call the function for creating netcdf for each of the text output files separately
    # # Use the name of text file that was saved; Change the filename in main script if so required
    # # TODO Parallelize this. Convenient when using multiple nodes to solve whole year in one go. 
    # text2nc("SWE", 1, outRaster)
    # text2nc("Gmelt", 2, outRaster)
    # text2nc("G", 3, outRaster)
    # text2nc("Precip", 4, outRaster)
    # text2nc("Us", 5, outRaster)
    # text2nc("Gpv", 6, outRaster)
    # text2nc("Gmeltpv", 7, outRaster)
    # text2nc("Upv", 8, outRaster)
    # text2nc("SWEpv", 9, outRaster);
else
    @info("All pixels not yet processed, so OUTPUT NETCDF FILES not yet created")
end
end_time = time_ns()
running_time = (end_time - start_time)/1e9/60  # minutes
if running_time < 60
    @info("Grant Total Running Time = $(round(running_time, digits=2)) minutes")
else
    running_time = running_time/60 # convert to hours
    @info("Grant Total Running Time = $(round(running_time, digits=2)) hours")
end
