"""
USAGE: julia verse/Julia/call_Blender_v20.jl test_WY2010 3005 3006 2
    - start and end indices are for y-axis only. All x's selected by default  

SWE estimation using Blender algorithm.
This is a helper script to organize inputs for calling the Blender function (Estimate_v61.jl).

The script is currently capable of running on my following platforms.
- Windows machine
- Ubuntu (WSL2-based)
- OSC HPC
- Ohio-State Unity HPC 
Julia Versions >= 1.11.5 
Minimum extra Julia Packages to install  
- JuMP, Ipopt, Rasters, NCDatasets, LoggingExtras, Distributions 

Approach
========
This script will call Estimate_v61.jl for each script at a time, and save temporary outputs to text file. Thus can readily be parallelzed.
Threads was stright forward but did not work here due to known limitation of Ipopt and threads module.
Working on alternative parallelization scheme using DISTRIBUTED module
Finally, all text output files are assembled into a nc file. The script can thus run on parts of pixels at different times, and finally combined into one nc file.
==============================================================================================
Jun 20, 2024 : call_Blender_v18.jl uses Estimate_v59.jl
    - read five input files separately rather than one merged file: required for 1 km run due to data volume
    - uses Qg_tavg as template replacing SWE_tavg.
    - renames MODSCAG to SCF
    - SCF is np.uint8, thus divide by 100 to get fraction
    - Snowf_tavg, SWE_tavg, Tair_f_tavg are np.unit16
        - Snowf_tavg and SWE_tavg divide by 1000 get floating point values in meters.
        - Tair_f_tavg divide by 100 to get floating point values in Kelvin.
Nov 20, 2024 : call_Blender_v19.jl uses Estimate_v60.jl . Currently in prototyping phase on a GFix branch.  
Dec 17, 2024 : Adding option to run for just one pixel for testing. 
Mar 02, 2025 : Removing "RES" everywhere to simplify the script. 
Mar 20, 2025 : Removing test runs because both pixel and watershed runs are now available and can be concidered as test run.
Apr 20, 2025 : remove opt and modis fix flag from arguments list to simplify and prepare for new continental runs.
"""
arg_len = length(ARGS)
out_subfolder = ARGS[1]  # WY2016. output subfolder relative to input files; temp_text and nc_outputs saved here
water_year = out_subfolder[end-3:end] #last 4 chars are assumed year, else error. 
start_idx = ARGS[2]  # this is string
stop_idx = ARGS[3]
start_idx = parse(Int64, start_idx)
stop_idx = parse(Int64, stop_idx)
# Extra arguments for various test runs: pixel, watershed, modis_fix etc.
if arg_len > 3
    ws_pix_idx = parse(Int64, ARGS[4])  # watershed or pixel index. hardcoded below if run for wshed.
end
# fix_modis_flag = 0  # default=0 means don't apply fix to MODIS.
# if arg_len > 4
#     fix_modis_flag = parse(Int64, ARGS[5]) # 0=false, 1=true  # to fix MODIS using Jack's approach.
# end

# Get path to the directory of the script. This is only used to access the test data for pixel and watershed runs
verse_dir = joinpath(splitpath(@__DIR__)[1:end-1])  # get parent folder of current script; later used to retrieve pixel or wshed csv file

pixel_run = false
wshed_run = false
if occursin("pixel", out_subfolder)  
    # if pixel substring is part of output subfolder then run on one pixel
    pixel_run = true
elseif occursin("wshed", out_subfolder)  # if wshed prefix substring is part of output subfolder
    wshed_run = true
end
start_time = time_ns()

# randint = 2  # abs(rand(Int8(1)))  # random integer for log file name (wshed and pixel run)

using Logging, LoggingExtras
using Tar
using DelimitedFiles
using Distributed  # otherwise everywhere macro won't work; Oct 8 2024: Error in Julia 1.11.0: Illegal instruction
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1. Setup inputs and output directories and file locations
# select the root directory (this will be different on windows, linux or Mac) 
host_machine = gethostname()
if occursin("L-JY0R5X3", host_machine)
    system_machine = "Windows"  # my laptop (also works for WLS linux)
    if Sys.iswindows()
        root_dir = "D:"  # Windows 10 machine
    elseif Sys.islinux()
        root_dir = "/mnt/d"  # Ubuntu (WSL2 on windows machine)
    end
    base_folder = "$root_dir/coressd/Blender"
    DataDir = "$root_dir/coressd/Blender/Inputs"  # must exist
    OUTDIR = "$base_folder/Runs"  # will be created if missing
    log_filename = string("LOGS/", start_idx, "_", stop_idx, ".log")  # on HPC created inside computer local node, so move to outside at end of job
    addprocs()
else
    # Prepare folder, environment and workers for parallel compute for SLURM
    system_machine = "Slurm"
    log_filename = string(ENV["SLURM_SUBMIT_DIR"], "/LOGS/",start_idx, "_", stop_idx, ".log")  # on HPC created inside computer local node, so move to outside at end of job
    if wshed_run
        log_filename = string(ENV["SLURM_SUBMIT_DIR"], "/LOGS/wshed", water_year, ws_pix_idx, ".log")
    elseif pixel_run
        log_filename = string(ENV["SLURM_SUBMIT_DIR"], "/LOGS/pixel", water_year, ws_pix_idx, ".log")
    end
    # cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"])  # ERROR: LoadError: KeyError: key "SLURM_CPUS_PER_TASK" not found [when not supplied on slurm scipt]
    cores = parse(Int, ENV["SLURM_NTASKS"])  # pick ntasks from slurm job script. must be provided.    
    # # cores = parse(Int, ENV["SLURM_JOB_CPUS_PER_NODE"])  # if one node provided; else we have parse as follows
    # subs = Dict("x"=>"*", "(" => "", ")" => "");
    # cores = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))  # for 1 or more nodes (generic); Int64 
    if occursin("borg", host_machine)  # TODO: discover
        root_dir = "/discover/nobackup"
        base_folder = "$root_dir/projects/coressd/Blender"
        tmpdir =  ENV["LOCAL_TMPDIR"]  # tempdir() to save tempoary text files on hpc node. 
        DataDir = "$root_dir/projects/coressd/Blender/Inputs"  # must exist
        OUTDIR = "$base_folder/Runs"  # will be created if missing
    elseif occursin(".osc.edu", host_machine)
        root_dir = "/fs/ess/PAS1785"  # "/fs/scratch/PAS1785/coressd"
        base_folder = "$root_dir/coressd/Blender"
        tmpdir =  ENV["TMPDIR"]  #tempdir()
        DataDir = "$root_dir/coressd/Blender/Inputs"
        OUTDIR = "$base_folder/Runs"
    elseif occursin("asc.ohio-state.edu", host_machine)
        root_dir = "/fs/project/howat.4/yadav.111"  # homedir()
        base_folder = "$root_dir/coressd/Blender"
        tmpdir =  ENV["TMPDIR"]  # tempdir()
        DataDir = "$root_dir/coressd/Blender/Inputs"
        OUTDIR = "$base_folder/Runs"
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
# exclusive option will also use all cores which may result in memory error.  
mkpath(dirname(log_filename))  # dirname gets the directory part of name (ie, without the filename) ; mkpath("logs")
logger = FormatLogger(open("$log_filename", "a+"), 0) do io, args  # w=write (original); a(a+)=(read), write,create,append
    # Write the module, level and message only
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end
info_only_logger = MinLevelLogger(logger, Logging.Info);  # Logging.Error
global_logger(info_only_logger)  # Set the global logger to logger; else logs doing directly to console

using Dates
@info("================================== START ==================================")
@info("ARGS: $ARGS" )  # print the script name and arguments to the console
@info("PROGRAM_FILE: $(PROGRAM_FILE)")  # print the script name and arguments to the console
@info("@__DIR__  $(@__DIR__)")   # directory of this script
@info("@__FILE__ $(@__FILE__)" )  # name of this script
@info("verse_dir: $verse_dir" )
@info("Script Run on $(Dates.now())")
@info("Water Year : $water_year")
@info("base_folder : $base_folder")
@info("Host computer machine: $host_machine")
@info("Number of processes (procs): $(nprocs())")
@info("Number of workers: $(nworkers())")

# Everywhere should come after cores are already added
@everywhere using Rasters, NCDatasets, Distributions  # Distributions for garbage collection using uniform distribution
@everywhere include("Estimate_v61.jl")  #https://docs.julialang.org/en/v1/manual/code-loading/; evaluated in global scope

# Make a folder insise HPC node because we want to copy existing files there
exp_dir = "$out_subfolder/outputs_txt_$(start_idx)_$(stop_idx)"  # To save text outputs for each pixel
mkpath(exp_dir)  # mkdir
logDir = "$out_subfolder/logs_$(start_idx)_$(stop_idx)"  # Save logs created by JuMP. 
# logDir = "$base_folder/Runs/$out_subfolder/logs/logs_$(start_idx)_$(stop_idx)"  # For Debug: save on same folder and outputs  
mkpath(logDir)  # mkdir; must create here, else error in the current setup  

# cp("$base_folder/Runs/$out_subfolder/outputs_txt", "$tmp_txtDir/$water_year")  # copy to local machine; error if running the first time as this dir would not exist
nc_outDir = "$OUTDIR/$out_subfolder/temp_nc/outputs_$(start_idx)_$(stop_idx)"
# if !ispath(exp_dir * "/SWEpv.txt") # ispath is generic for both file and folder
# if !isfile(exp_dir * "/SWEpv.txt")  # This is better check: process the pixel only if the last file exported by optimzer (SWEpv.txt) does not yet exist
if isdir(nc_outDir)  # process only if the pixel is not already processed
    @info("Exiting because following folder already processed: $nc_outDir")
    exit(0)
end

# 2. Read the Input netCDF file
files = (
    "$DataDir/WY$(water_year)/SCF.nc",
    "$DataDir/WY$(water_year)/Snowf_tavg.nc",
    "$DataDir/WY$(water_year)/SWE_tavg.nc",
    "$DataDir/WY$(water_year)/Tair_f_tavg.nc",
    "$DataDir/WY$(water_year)/Qg_tavg.nc"
    )
A = RasterStack(files; lazy=true)  # lazy=true https://github.com/rafaqz/Rasters.jl/issues/449

running_time = (time_ns() - start_time)/1e9/60
@info("Time until creating Rasterstack from separate input netcdf files = $(round(running_time, digits=2)) minutes")

# Subset only the required variables because the nc file can have extraneous vars that cause problem with julia
A = A[(:Snowf_tavg, :SWE_tavg, :Tair_f_tavg, :Qg_tavg, :SCF)];  # to remove spatial ref that cause problem during subsetting
""" General info about pixel counts
    For 0.050 degree resolution, size(A) = (2336, 941); Non-missing pixel count = 1011329
    For 0.010 degree resolution, size(A) = (11700, 6500); Non-missing pixel count = 

"""
szY = size(A, 2)  # get size in Y-dimension; here dim = 2
if stop_idx > szY
    stop_idx = szY
end
if pixel_run
    @info("Pixel Run  ")
    # To select just one pixel. Index of coordinates must be passed as vector (ie [1101] not 1101), esle it will loose the dimension of X, Y and subsequent code with cartesian index will not work.
    # A = A[[1101], [3001], :]  # using index for X and Y respectively
    # A = A[X=Near([-100.2]), Y=Near([50.1])]  # Using lon/lat. Ti is optional here.  
    # A = A[X=Near([-119.348099]), Y=Near([37.876408])]  # Using lon/lat. Ti is optional here.  
    # Passing pixels using csv file
    # data_cells, header_cells = readdlm("$base_folder/coordinates/pixels.csv", ',', header=true, skipblanks=true)
    pixel_csv = joinpath(verse_dir, "data", "pixels.csv")  # append the rest of the path
    data_cells, header_cells = readdlm(pixel_csv, ',', header=true, skipblanks=true)

    id, pix_name, x0, y0, wshed, wsid = data_cells[ws_pix_idx, :]
    @info("Pixel Index, Name and bounding coords: $id $pix_name $x0  $y0 $wshed $wsid")
    A = A[X=Near([x0]), Y=Near([y0])]  # Using lon/lat. Ti is optional here.  
    # To save inputs for making plotting and analysis 
    subset_folder = "$OUTDIR/$out_subfolder/Inputs"  # file name with fullpath for subset
    mkpath(subset_folder)
    write("$subset_folder/pixel_$ws_pix_idx.nc", A, force=true)  # !isfile("$subset_folder/pixel.nc") && ## error when two processors try to write at the same time
elseif wshed_run
    @info("Watershed Run  ")
    # # A = A[X=Near([-100.2]), Y=Near([50.1])]  # Using lon/lat. Ti is optional here.  
    # For watershed: give lower left and upper right longitude/latitude as corners of the bounding box (or hardcode here).
    # data_cells, header_cells = readdlm("$base_folder/coordinates/wshed.csv", ',', header=true, skipblanks=true)
    wshed_csv = joinpath(verse_dir, "data", "wshed.csv")  # append the rest of the path
    data_cells, header_cells = readdlm(wshed_csv, ',', header=true, skipblanks=true)
    id, wshed_name, x0, x1, y0, y1 = data_cells[ws_pix_idx, :]
    @info("Watershed Index and bounding coords: $id $wshed_name $x0 $x1 $y0 $y1")
    # # Hardcoded to Tuolumne watershed
    # x0, x1 = -119.66, -119.19  # longitude
    # y0, y1 = 37.73, 38.12      # latitude    
    # We can use either of the following two methods to subset the data. But using view was much slower.
    # A = view(A, X(x0 .. x1), Y(y0 .. y1))  # using view to select a small chip around watershed, but seems slower
    A = A[X(Between(x0, x1)), Y(Between(y0, y1))]  # Another way of getting the same chip around watershed
    # To save the clipped input file for future use
    subset_folder = "$OUTDIR/$out_subfolder/Inputs"  # file name with fullpath for subset
    mkpath(subset_folder)
    write("$subset_folder/wshed_$ws_pix_idx.nc", A, force=true)  ## error when two processors try to write at the same time
else
    @info("DEFAULT RUN  ")
    # A = A[1:end, start_idx:stop_idx, :]  # use for default runs
    A = A[X(1:end), Y(start_idx:stop_idx)]
    # A = A[2109:5329, start_idx:stop_idx, :]  # Uncomment this for debugging and testing.
end
# Use cartesian index for iteration
# valid_pix_ind = findall(!ismissing, A["SWE_tavg"][Ti=1])  # Extract cartesian index for non-missing (ie, all) data using one-day of data 
valid_pix_ind = findall(!ismissing, A["Qg_tavg"][Ti=1])  # Extract cartesian index for non-missing (ie, all) data using one-day of data 
valid_pix_count = length(valid_pix_ind)
ind = valid_pix_ind
@info("Total valid pixel count = $(valid_pix_count)")
# println("Total valid pixel count = $(valid_pix_count)")  # debug
# y_dim = dims(A, Y)  # # Returns a Y-dimensional (latitude) lookup array
# @info("y_dim = $(collect(y_dim))")  # collect will remove extra metadata and just print values

running_time = (time_ns() - start_time)/1e9/60
@info("Elapsed Time = $(round(running_time, digits=2)) minutes")
@info("Starting with blender loop.")
@info("---------------------------")
@info("Processing $(length(ind)) of $(length(valid_pix_ind))")

# Next few lines are for testing twindow and σWRFGmin only.
# σWRFGmin_list = collect(0:5:200)  # as range 1:5:120
# # σWRFGmin_list[1] = 1  # to ensure that 1 is included in the list
# # σWRFGmin_list = [1, 5, 10, 15, 25, 50, 100, 200, 300, 400,500]
# σWRFGmin_list = [10]
# @info("Fix MODIS flag = $fix_modis_flag")
# twindow_list = collect(0:1:15)  # for smoothing SCF.
# @info(twindow_list)

# @sync makes the code wait for all processes to finish their part of the computation before moving on from the loop
# ie, without @sync, the on of the processors may move to next while loop is still running, creating error and crash whole script prematurely 
# Threads.@threads for ij in ind
# @distributed for ij in ind
@sync @distributed for ij in ind
# for ij in ind
    # @info("Processing pixel $(ij[1]) $(ij[2])")
    i = ij[1]
    j = ij[2]
    # lat = y_dim[j]  # pixel latitude. can be used to correct MODIS SCF based on latitude 
    # TODO replace the following 5 lines with the direct access to file
    WRFSWE = A["SWE_tavg"][X=i, Y=j].data/1000  # Here "data" is an AbstractArray.
    WRFP = A["Snowf_tavg"][X=i, Y=j].data/1000
    WRFG = A["Qg_tavg"][X=i, Y=j].data
    AirT = A["Tair_f_tavg"][X=i, Y=j].data/100
    MSCF = A["SCF"][X=i, Y=j].data; #/100 Convert to fraction. multiplication and devision works without dot. but add/subtract will need dot.
    # TODO: extract latitude and pass to blender function. can be used to spatially fix MODIS SCF
    # blender(out_folder, i, j, WRFSWE, WRFP, WRFG, MSCF, AirT)  # Call blender for the pixel. OLD.
    # @sync @distributed for twindow in twindow_list  # σWRFGmin in σWRFGmin_list
    #     # blender(i, j, WRFSWE, WRFP, WRFG, MSCF, AirT, logDir, exp_dir, opt, σWRFGmin, fix_modis_flag)
    #     blender(i, j, WRFSWE, WRFP, WRFG, MSCF, AirT, logDir, exp_dir, twindow)
    # end
    # println("MSCF: ", MSCF)
    # println("WRFG: ", WRFG)
    blender(i, j, WRFSWE, WRFP, WRFG, MSCF, AirT, logDir, exp_dir)
end
# without sync above one of the processor to the next step (combining text files to netcdf) which will cause error
sleep(10)
@everywhere GC.gc()

end_time = time_ns()
running_time = (end_time - start_time)/1e9/60  # minutes. or divide by 3600 for seconds
efficiency = (running_time * 60 * nworkers()) / valid_pix_count  # processing time (seconds) per pixel  
@info("Pixel processing Efficiency (seconds per pixel) = $(round(efficiency, digits=2)) \n")
if running_time < 60
    @info("Running Time (blender For Loop) = $(round(running_time, digits=2)) minutes")
else
    running_time = running_time/60 # convert to hours
    @info("Running Time (blender For Loop) = $(round(running_time, digits=2)) hours")
end
if pixel_run
    exit(0)  # exit here to avoid running 2nd part (text2nc)
end
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
    # var_idx_tuple = (("SWE", 1), ("Gmelt", 2), ("G", 3), ("Precip", 4), ("Us", 5), ("Gpv", 6), ("Gmeltpv", 7), ("Upv", 8), ("SWEpv", 9))
    # var_idx_tuple = (("SWE", 1), ("Gmelt", 2), ("G", 3), ("Precip", 4), ("Us", 5), ("Gpv", 6), ("Gmeltpv", 7), ("Upv", 8), ("SWEpv", 9), ("Gmelt_prior", 10), ("sigmaWRFG", 11)) 
    var_idx_tuple = (("SWE", 1), ("Gmelt", 2), ("G", 3), ("Precip", 4), ("Us", 5), ("Gpv", 6), ("Gmeltpv", 7), ("Upv", 8), ("SWEpv", 9), ("SCFobs", 10))
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
@info("================================== END ==================================\n")
