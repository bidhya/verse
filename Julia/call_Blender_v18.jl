"""
USAGE: julia verse/Julia/call_Blender_v18.jl test_WY2010 3005 3006 010
    - start and end indices are for y-axis only. All x's selected by default  

SWE estimation using Blender algorithm.
This is a helper script to organize inputs for calling the Blender function (Estimate_v59.jl).

The script is currently capable of running on my following platforms.
- Windows machine
- Ubuntu (WSL2-based)
- OSC HPC
- Ohio-State Unity HPC 
Julia Versions >= 1.10.5 
Minimum extra Julia Packages required
- JuMP, Ipopt, Rasters, NCDatasets 

Approach
========
This script will call Estimate_v59.jl for each script at a time, and save temporary outputs to text file. Thus can readily be parallelzed.
Threads was stright forward but did not work here due to known limitation of Ipopt and threads module.
Working on alternative parallelization scheme using DISTRIBUTED module
Finally, all text output files are assembled into a nc file. The script can thus run on parts of pixels at different times, and finally combined into one nc file.
==============================================================================================
Sep 04, 2023 : text and log files saved in separate folders; created call_Blender_v12.jl uses Estimate_v56.jl
Oct 29, 2023 : Save nc files to temp_nc folder using call_Blender_v14.jl uses Estimate_v57.jl
Dec 03, 2023 : New call_Blender_v15.jl uses Estimate_v58.jl (modifications by Jack)
Jun 20, 2024 : call_Blender_v18.jl uses Estimate_v59.jl
    - read five input files separately rather than one merged file: required for 1 km run due to data volume
    - uses Qg_tavg as template replacing SWE_tavg.
    - renames MODSCAG to SCF
    - SCF is np.uint8, thus divide by 100 to get fraction
    - Snowf_tavg, SWE_tavg, Tair_f_tavg are np.unit16
        - Snowf_tavg and SWE_tavg divide by 1000 get floating point values in meters.
        - Tair_f_tavg divide by 100 to get floating point values in Kelvin.

"""
arg_len = length(ARGS)
out_subfolder = ARGS[1]  # WY2016. output subfolder relative to input files; temp_text and nc_outputs saved here
water_year = out_subfolder[end-3:end] #last 4 chars are assumed year, else error. out_subfolder[3:end]
start_idx = ARGS[2]  # this is string
end_idx = ARGS[3]
start_idx = parse(Int64, start_idx)
end_idx = parse(Int64, end_idx)
RES = ARGS[4] # "050"  # Grid Resolution: 050 = 0.050 degree (5km); 025 = 0.025 degree; 010 = 0.01 degree; 001 = 0.001 degree (100 meters) etc

if occursin("test", out_subfolder)  # TODO: either replace this with wshed or create different criteria for watershed run.  
    # if test substring is part of output subfolder then do the test run on subset of pixels
    test_run = true
    wshed_run = false
elseif occursin("wshed", out_subfolder)  # if wshed prefix substring is part of output subfolder
    wshed_run = true
    test_run = false
else
    # we define both these varaibles to avoid below where we check for both test or wshed run.
    test_run = false
    wshed_run = false
end
start_time = time_ns()

using Logging, LoggingExtras
using Tar
using DelimitedFiles
using Distributed  # otherwise everywhere macro won't work; Oct 8 2024: Error in Julia 1.11.0: Illegal instruction
# using SharedArrays  # move this line below addprocs(cores), else won't work in Julia 1.10.0
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Setup inputs and output directories and file locations
# select the root directory (this will be different on windows, linux or Mac) 
host_machine = gethostname()
# if occursin("STAFF-BY-M", host_machine)
if occursin("L-JY0R5X3", host_machine)
    system_machine = "Windows"  # a bit misnomer flag because this flag also works for WLS linux. Better flag could be local vs hpc/remote execution
    if Sys.iswindows()
        root_dir = "D:"  #for windows machine
    elseif Sys.islinux()
        root_dir = "/mnt/d"  #for Ubuntu (WSL2 on windows machine
    end
    base_folder = "$root_dir/coressd/Blender"
    DataDir = "$root_dir/coressd/Blender/Inputs_$(RES)"  # must exist  (Old = NoahMP)
    OUTDIR = "$base_folder/Runs/$(RES)"
    log_filename = string(start_idx, "_", end_idx, ".log")  # on HPC created inside computer local node, so move to outside at end of job
    addprocs()
else  
    # Prepare folder, environment and workers for parallel compute for SLURM
    system_machine = "Slurm"
    log_filename = string(ENV["SLURM_SUBMIT_DIR"], "/",start_idx, "_", end_idx, ".log")  # on HPC created inside computer local node, so move to outside at end of job
    # cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"])  # ERROR: LoadError: KeyError: key "SLURM_CPUS_PER_TASK" not found [when not supplied on slurm scipt]
    cores = parse(Int, ENV["SLURM_NTASKS"])  # pick ntasks from slurm job script. must be provided.    
    # # cores = parse(Int, ENV["SLURM_JOB_CPUS_PER_NODE"])  # if one node provided; else we have parse as follows
    # subs = Dict("x"=>"*", "(" => "", ")" => "");
    # cores = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))  # for 1 or more nodes (generic); Int64 
    @info("SLURM Cores: $cores")
    if occursin("borg", host_machine)  # TODO: discover
        root_dir = "/discover/nobackup"
        base_folder = "$root_dir/projects/coressd/Blender"
        # base_folder = "$root_dir/byadav/coressd/Blender"  # used when project directory was full.  
        tmpdir =  ENV["LOCAL_TMPDIR"]  #tempdir() to save tempoary text files on hpc node. 
        # DataDir = "$root_dir/projects/coressd/Blender/Inputs"  # INDIR. must exist  (Old = NoahMP)
        DataDir = "$root_dir/projects/coressd/Blender/Inputs_$(RES)"  # INDIR. must exist  (Old = NoahMP)
        OUTDIR = "$base_folder/Runs/$(RES)"  # Blender Run outputs saved here. (../temp_nc/, ../logs/, ../Outputs/ etc)
    elseif occursin(".osc.edu", host_machine)
        root_dir = "/fs/ess/PAS1785"  # "/fs/scratch/PAS1785/coressd"
        base_folder = "$root_dir/coressd/Blender"
        tmpdir =  ENV["TMPDIR"]  #tempdir()
        DataDir = "$root_dir/coressd/Blender/Inputs_$(RES)"  # Inputs
        OUTDIR = "$base_folder/Runs/$(RES)"
    elseif occursin("asc.ohio-state.edu", host_machine)  # .unity
        root_dir = "/fs/project/howat.4/yadav.111"  # homedir()  #  Unity
        # base_folder = "/home/yadav.111/Github/Blender"  # old
        base_folder = "$root_dir/coressd/Blender"  # "$root_dir/Github/coressd/Blender"
        tmpdir =  ENV["TMPDIR"]  #tempdir()
        DataDir = "$root_dir/coressd/Blender/Inputs_$(RES)"
        OUTDIR = "$base_folder/Runs/$(RES)"
        # @info("SLURM_SUBMIT_DIR: $(ENV["SLURM_SUBMIT_DIR"])")
    else
        @info("Unknown computer system. Aborting ...  ")
        exit() # <- you can provide exit code if you want
    end
    addprocs(cores)  # ; exeflags="--project"subtract one becuase master already has one; but seems to work with higher number as well
    # using SlurmClusterManager  # to pick all nodes and cores allocated in slurm job
    # # With slurmcluster manager, hopefully next few lines of selecting cores is not necessary  
    # sleep(10)  # 60+cores To prevent scheduling job on more than 1 node on Discover Slurm cluster  
    # addprocs(SlurmManager())  # to use all available nodes and cores automatically. comment this line and uncomment one above this to match _v8.jl
end
# addprocs()  # Works but on cluster will use all cores, even though we not asked for. 
# also with exclusive option will use all cores which may result in memory error.  
logger = FormatLogger(open(log_filename, "w"), 0) do io, args
    # Write the module, level and message only
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end
info_only_logger = MinLevelLogger(logger, Logging.Info);  # Logging.Error
global_logger(info_only_logger)  # Set the global logger to logger; else logs doing directly to console

@info("base_folder : $base_folder")
@info("Host computer machine: $host_machine")
@info("Number of processes (procs): $(nprocs())")
@info("Number of workers: $(nworkers())")

# Everywhere should come after cores are already added
@everywhere using Rasters, NCDatasets, Distributions  # Distributions for garbage collection using uniform distribution
@everywhere include("Estimate_v59.jl")  #https://docs.julialang.org/en/v1/manual/code-loading/; evaluated in global scope

# DataDir = "$base_folder/Inputs"  # must exist  (Old = NoahMP)
# DataDir = "$root_dir/projects/coressd/Blender/Inputs"  # must exist  (Old = NoahMP)

# Make a folder insise HPC node because we want to copy existing files there
exp_dir = "outputs_txt_$(start_idx)_$(end_idx)"  #"$out_folder/outputs_txt"  # tmp_txtDir(old name) To save text outputs for each pixel
mkpath(exp_dir)  # mkdir
logDir = "logs_$(start_idx)_$(end_idx)"  # "logs"  "$out_folder/logs"   # Save logs in separate folder
# logDir = "$base_folder/Runs/$out_subfolder/logs/logs_$(start_idx)_$(end_idx)"  # For Debug: save on same folder and outputs  
mkpath(logDir)  # mkdir; must create here, else error in the current setup  
# Error on Discover (Nov 08, 2023): rm: cannot remove '/lscratch/tdirs/batch/slurm.24967584.byadav/logs': Directory not empty; Solution: #SBATCH --no-requeue

# cp("$base_folder/Runs/$out_subfolder/outputs_txt", "$tmp_txtDir/$water_year")  # copy to local machine; error if running the first time as this dir would not exist
nc_outDir = "$OUTDIR/$out_subfolder/temp_nc/outputs_$(start_idx)_$(end_idx)"
# if !ispath(exp_dir * "/SWEpv.txt") # ispath is generic for both file and folder
# if !isfile(exp_dir * "/SWEpv.txt")  # This is better check: process the pixel only if the last file exported by optimzer (SWEpv.txt) does not yet exist
if isdir(nc_outDir)  # process only if the pixel is not already processed
    @info("Exiting because following folder already processed: $nc_outDir")
    exit(0)
end

# # # 2. Read the Input netCDF file
# # A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_modscag.nc")  #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
# # Following check are for prototyping only when running code locally, because I do not yet have NorthAmerica netcdf file
# if occursin("L-JY0R5X3", host_machine) #|| test_run  # use this for testing as well
#     @info("Test Run only.")
#     A = RasterStack("$(DataDir)/WY_merged/$(water_year)_seup_modis.nc", lazy=true)
#     # subset = A[X(Between(-120.5, -120)), Y(Between(60, 60.5))]  # 10 by 10 pixels; use small chip for prototyping
# else
#     # # copy to local Node (machine) on slurm
#     # # cp("$DataDir/WY_merged/$water_year" * "_seup_modis.nc", "$tmpdir/$water_year" * "_seup_modis.nc", force=true)  #force=true will first remove an existing dst.
#     # # cp("$DataDir/lis/WY$(water_year)/SWE_tavg.nc", "$tmpdir/WY$(water_year)/SWE_tavg.nc", force=true)  #force=true will first remove an existing dst.
#     # # true required for discover when running on same node again after node_failure.
#     # # A = RasterStack("$tmpdir/$water_year" * "_seup_modis.nc", lazy=true)  ## , lazy=true https://github.com/rafaqz/Rasters.jl/issues/449
#     # A = RasterStack("$tmpdir/WY$(water_year)/SWE_tavg.nc", lazy=true)

#     # June 04, 2024
#     # cp("$DataDir/lis/WY$(water_year)", "$tmpdir/WY$(water_year)", force=true)  #force=true will first remove an existing dst. New for v18
#     # Full path to Datadir
#     files = (
#         "$DataDir/lis/WY$(water_year)/SCF.nc",
#         "$DataDir/lis/WY$(water_year)/Snowf_tavg.nc",
#         "$DataDir/lis/WY$(water_year)/SWE_tavg.nc",
#         "$DataDir/lis/WY$(water_year)/Tair_f_tavg.nc",
#         "$DataDir/lis/WY$(water_year)/Qg_tavg.nc"
#         )
#     A = RasterStack(files; lazy=true)
# end

files = (
    "$DataDir/lis/WY$(water_year)/SCF.nc",
    "$DataDir/lis/WY$(water_year)/Snowf_tavg.nc",
    "$DataDir/lis/WY$(water_year)/SWE_tavg.nc",
    "$DataDir/lis/WY$(water_year)/Tair_f_tavg.nc",
    "$DataDir/lis/WY$(water_year)/Qg_tavg.nc"
    )
A = RasterStack(files; lazy=true)


# end_time = time_ns()
running_time = (time_ns() - start_time)/1e9/60
@info("Time until copying input netcdf to Node = $(round(running_time, digits=2)) minutes")

# Subset only the required variables because the nc file can have extraneous vars that cause problem with julia
A = A[(:Snowf_tavg, :SWE_tavg, :Tair_f_tavg, :Qg_tavg, :SCF)];  # to remove spatial ref that cause problem during subsetting
""" General info about pixel counts
    For 0.050 degree resolution, size(A) = (2336, 941); Non-missing pixel count = 1011329
    For 0.010 degree resolution, size(A) = (11700, 6500); Non-missing pixel count = 

"""
szY = size(A, 2)  # get size in Y-dimension; here dim = 2
if end_idx > szY
    end_idx = szY
end
if test_run
    @info("TEST RUN ONLY  ")
    # Aside: Get test set of data. This part needs re-writing before using for test [TODO].
    A = A[1101:1109, start_idx:end_idx, :]  # for test run
elseif wshed_run
    @info("Watershed Run  ")
    # For watershed: give lower left and upper right longitude/latitude as corners of the bounding box (or hardcode here).
    # TODO: read bounding box from a textfile or use shapefile.   
    if arg_len > 4
        ws_idx = ARGS[5]  # read the index from command line
        ws_idx = parse(Int64, ws_idx)
        data_cells, header_cells = readdlm("$base_folder/coordinates/wshed.csv", ',', header=true, skipblanks=true)
        wshed_name, x0, x1, y0, y1 = data_cells[ws_idx, :]  
        @info("Watershed and bounding coords: $wshed_name $x0 $x1 $y0 $y1")
    else
        # Hardcoded to Tuolumne watershed
        x0, x1 = -119.66, -119.19  # longitude
        y0, y1 = 37.73, 38.12      # latitude    
    end
    # We can use either of the following two methods to subset the data. But using view was much slower.
    # A = view(A, X(x0 .. x1), Y(y0 .. y1))  # using view to select a small chip around watershed, but seems slower
    A = A[X(Between(x0, x1)), Y(Between(y0, y1))]  # Another way of getting the same chip around watershed
    # To save the clipped input file for future use
    subset_folder = "$OUTDIR/$out_subfolder/Inputs"  # file name with fullpath for subset
    mkpath(subset_folder)
    write("$subset_folder/wshed.nc", A, force=true)
else
    @info("DEFAULT RUN  ")
    A = A[1:end, start_idx:end_idx, :]  # use for default runs
end
# Use cartesian index for iteration
# valid_pix_ind = findall(!ismissing, A["SWE_tavg"][Ti=1])  # Extract cartesian index for non-missing (ie, all) data using one-day of data 
valid_pix_ind = findall(!ismissing, A["Qg_tavg"][Ti=1])  # Extract cartesian index for non-missing (ie, all) data using one-day of data 
valid_pix_count = length(valid_pix_ind)
ind = valid_pix_ind
@info("Total valid pixel count = $(valid_pix_count)")
# println("Total valid pixel count = $(valid_pix_count)")  # debug

running_time = (time_ns() - start_time)/1e9/60
@info("Pre-generating output nc files = $(round(running_time, digits=2)) minutes")
@info("Starting with blender loop.")
@info("==========================================================")
@info("Processing $(length(ind)) of $(length(valid_pix_ind))")
# @sync makes the code wait for all processes to finish their part of the computation before moving on from the loop
# ie, without @sync, the on of the processors may move to next while loop is still running, creating error and crash whole script prematurely 
# Threads.@threads for ij in ind
# @distributed for ij in ind
@sync @distributed for ij in ind
# for ij in ind
    i = ij[1]
    j = ij[2]
    # TODO replace the following 5 lines with the direct access to file
    WRFSWE = A["SWE_tavg"][X=i, Y=j].data/1000  # Here "data" is an AbstractArray.
    WRFP = A["Snowf_tavg"][X=i, Y=j].data/1000
    WRFG = A["Qg_tavg"][X=i, Y=j].data
    AirT = A["Tair_f_tavg"][X=i, Y=j].data/100
    MSCF = A["SCF"][X=i, Y=j].data; #/100 Convert to fraction. multiplication and devision works without dot. but add/subtract will need dot.
    # blender(out_folder, i, j, WRFSWE, WRFP, WRFG, MSCF, AirT)  # Call blender for the pixel. OLD.
    blender(i, j, WRFSWE, WRFP, WRFG, MSCF, AirT, logDir, exp_dir)
end
# without sync above one of the processor to the next step (combining text files to netcdf) which will cause error
sleep(10)
@everywhere GC.gc()

end_time = time_ns()
running_time = (end_time - start_time)/1e9/3600
@info("Running Time (blender For Loop) = $(round(running_time, digits=2)) hours")
# exit(0)  # exit here to avoid running 2nd part (text2nc)
# ===================================================================================================================================

# Step 4. PostProcessing: Combine text files into a grid and save as netcdf
@info("Combining Text Ouputs to NetCDF FILES")
@everywhere function text2nc(var, idx, outRaster, exp_dir, nc_outDir)
    """Convert a text file to netcdf
    NB: to use with distributed, pass all vars that are used here.  
    var: the name of variable (example SWE, Gmelt, G etc.)
    idx : the index where the variable of exist in the file [1 through 9]
    TODO: move this fuction inside the Estimate_v53.jl script
    TODO: seems like nc file produced here cannot by plotted using xarray/hvplot
    """
    # Initial resulting array
    outRaster[:,:,:] .= missing
    # update the name of variable
    outRaster = rebuild(outRaster; name=var)  #:SWEhat
    pixels = readdir(exp_dir)  # read text-based blender output files if not passed a function parameter.  
    pixels = [pix for pix in pixels if startswith(pix, "Pix")];
    # @sync @distributed for pix in pixels  # worked but seems slow, even slower than vanilla for loop
    # Threads.@threads for pix in pixels  # can be much faster but causing memory related Segmentation Fault when combined with Distributed.    
    for pix in pixels  # this worked
        pix_xy = split(pix, ".txt")[1]
        pix_xy = split(pix_xy, "_")
        x = parse(Int16, pix_xy[2])
        y = parse(Int16, pix_xy[3])
        out_vars = readdlm("$(exp_dir)/$(pix)") # readdlm("$(exp_dir)/Pix_$(x)_$(y).txt")
        outRaster[X=x, Y=y] = out_vars[:,idx]
    end
    # Save to nc file; 
    # mkpath(nc_outDir)  # error in distributed
    write("$nc_outDir/$var.nc", outRaster)  # ERROR1 Aug 08, 2023: Error perhaps due to updates to raster/ncdataset/etc. (error tested on windows and wsl2)
    # June 05, 2024: NetCDF error: NetCDF: Not a valid data type or _FillValue type mismatch (NetCDF error code: -45)
end

sleep(1)
# Count if all pixels are processed before calling the following section for converting text files to netcdf file
pixels = readdir(exp_dir)
pixels = [pix for pix in pixels if startswith(pix, "Pix")];
processed_pix_count = length(pixels)
@info("Pixels processsed = $processed_pix_count out of $valid_pix_count")
if processed_pix_count == valid_pix_count #&& system_machine == "Windows"  # ERROR1
    outRaster = copy(A[:Qg_tavg])  # copy(A[:SWE_tavg])
    var_idx_tuple = (("SWE", 1), ("Gmelt", 2), ("G", 3), ("Precip", 4), ("Us", 5), ("Gpv", 6), ("Gmeltpv", 7), ("Upv", 8), ("SWEpv", 9))
    mkpath(nc_outDir)
    # Threads.@threads (Roughly same running for threaded and non-thread version on Discover. Likely due to I/O bottleneck)
    # Threads.@threads for var_idx in var_idx_tuple
    @sync @distributed for var_idx in var_idx_tuple
        var_name = var_idx[1]
        var_idx = var_idx[2]
        text2nc(var_name, var_idx, outRaster, exp_dir, nc_outDir);  # ERROR1
    end
else
    @info("All pixels not yet processed, so OUTPUT NETCDF FILES not yet created")
end
@info("Finished: Combining Text Ouputs to NetCDF FILES")
# end_time = time_ns()
# running_time = (end_time - start_time)/1e9/3600
# @info("Running Time (convert text to nc) = $(round(running_time, digits=3)) hours")

# # Create path for saving log files. tar and copy will be done in the bash script.  
# mkpath("$OUTDIR/$out_subfolder/logs")  # also works if path already exists
# @info("Logdir: $OUTDIR/$out_subfolder/logs")

end_time = time_ns()
running_time = (end_time - start_time)/1e9/60  # minutes
if running_time < 60
    @info("Grant Total Running Time = $(round(running_time, digits=2)) minutes")
else
    running_time = running_time/60 # convert to hours
    @info("Grant Total Running Time = $(round(running_time, digits=2)) hours")
end
