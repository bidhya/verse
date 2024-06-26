"""
USAGE: julia call_Blender_v5.jl [output_subfolder]

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
- 7.x (1, 2, and 3)
- 8.x (only ran on Ubuntu-WSL)
Minimum extra Julia Packages required
- JuMP, Ipopt, Rasters 

Approach
========
This script will call Estimate_v5x.jl for each script at a time, and save temporary outputs to text file. Thus can readily be parallelzed.
Threads was stright forward but did not work here due to known limitation of Ipopt and threads module.
Working on alternative parallelization scheme using DISTRIBUTED module
Finally, all text output files are assembled into a nc file. The script can thus run on parts of pixels at different times, and finally combined into one nc file.
TODO: put a flag for this second part of script that combines individual text files into a nc file. Currently, this part is executed for every run.

==============================================================================================

Call the blender fucntion :Estimate_v53.jl
Here, we subset what pixels to  process
Decide whether to process in Serial or Distributed way depending on the number of cores we have available

Sep 11, 2022 (Organizing input and output folders with NoahMP run)
Nov 11, 2022 : Last successful run for MODIS-CGF
Nov 11, 2022 : Updating next (_v6) to incorporate cartesian index for intering with just one loop
    helpful for parallel code
Nov 12, 2022 : Successfully incorporated cartesian index
Nov 19, 2022 : Starting with North America Run
Nov 20, 2022 : checking for existence of file rather than folder before calling Blender function. 
    Though not 100% full-proof yet, this will help process the pixel(s) at the point of failure due to early termination of runtime 
Nov 21, 2022 : incorporate logging 
Nov 21, 2022 : v6 is still good; but implementing logic here to process a subset of pixels due to limitation of 12 hours for Discover
Nov 22, 2022 : v8 is using Estimate_v54 where 9 text outputs are combined into 1 text file
Nov 22, 2022 : Updating text2nc function to create nc file from the combined text file
Nov 27, 2022 : Checking the output file before loading Rasters array; should save (signficantly) time for pixels already processed


#datadir is a moving directory containing an argument given in the job file
# in the shell script the arg is increased by a step counter
# same arg moves the pix number directory in each seq run

TODO 
====
Find better way of indexing/iterating through 2D-array rather than with two for loops with i and j.
"""
# DataDir= "/users/PAS0503/jldechow/foursix/AllTuolomne46/Pix3252"
# DataDir = ARGS[1]
# cores = ARGS[1]
# cores = parse(Int8, cores)  # maximum of Int8 is 127 only

# # This start, end index is no more useful because we are processing from netcdf rather than text files separately
# start_idx = parse(Int16, ARGS[1])
# end_idx = parse(Int16, ARGS[2])
# Base.ARGS :: Vector{String} >>>>> An array of the command line arguments passed to Julia, as strings.

arg_len = length(ARGS)
if arg_len == 0
    # no argument passed from command line
    out_subfolder = "NoahMP_CGF"
else
    out_subfolder = ARGS[1]     # output subfolder relative to input files; temp text and nc_outputs saved here
end
start_idx = ARGS[2]
end_idx = ARGS[3]
println(typeof(start_idx))
println(typeof(start_idx))

start_idx = parse(Int64, start_idx)
end_idx = parse(Int64, end_idx)
println(typeof(start_idx))
println(typeof(start_idx))

log_filename = string("log_", out_subfolder, "_", start_idx, "_", end_idx, ".txt")  #construct a log filename
# log_filename = string(".out/log_", out_subfolder, "_", start_idx, "_", end_idx, ".txt")  #construct a log filename
# Nov 20, 2022: Updated logger for finer control; still not working with distributed
using Logging, LoggingExtras
# io = open("log.txt", "a")
# logger = SimpleLogger(io)  # Create a simple logger
# logger = FileLogger("logfile.txt")
# logger = FileLogger("log.txt"; append = false)  # true
logger = FormatLogger(open(log_filename, "w"), 0) do io, args
    # Write the module, level and message only
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end
global_logger(logger)  # Set the global logger to logger; else logs doing directly to console


using Distributed  #otherwise everywhere macro won't work
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. Setup inputs and output directories and file locations
# select the root directory (this will be different on windows, linux or Mac) 
host_machine = gethostname()
println("Host computer machine: $host_machine")
if occursin("STAFF-BY-M", host_machine)
    if Sys.iswindows()
        root_dir = "C:"  #for windows machine
    elseif Sys.islinux()
        root_dir = "/mnt/c"  #for Ubuntu (WSL2 on windows machine
    else
        println("Unknown system on windows, manually add root directory before proceeding. Exiting code")
        exit(1)    
    end
    base_folder = "$root_dir/Github/coressd/Blender"
    # Get info on system machine. We use this info to farm number of processors. 
    system_machine = "Windows"  # a bit misnomer flag because this flag also works for WLS linux. Better flag could be local vs hpc/remote execution
elseif occursin("borg", host_machine)  # TODO: discover
    root_dir = "/discover/nobackup/projects/coressd"
    base_folder = "$root_dir/Blender"
    system_machine = "Slurm"
elseif occursin(".osc.edu", host_machine)
    root_dir = "/fs/ess/PAS1785/coressd" #homedir()  # OSC
    base_folder = "$root_dir/Blender"
    system_machine = "Slurm"

elseif occursin("asc.ohio-state.edu", host_machine)  # .unity
    root_dir = homedir()  #  Unity
    base_folder = "$root_dir/Github/Blender"
    base_folder = "$root_dir/Github/coressd/Blender"
    system_machine = "Slurm"
else
    println("Unknown computer, manually add root directory before proceeding. Exiting code")
    exit(1)  # if error, comment this line, add root and base_folder below and run again
end
println("base_folder : $base_folder")
# Supply custom root and base_folder for new machine
# root_dir = ""  
# base_folder = "$root_dir/..."
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if system_machine == "Windows" # addproc function works for both windows and WSL linux. Hence, this code is more for local vs remote machine
#   root_dir = "C:"
  addprocs()  ## works on laptop, but on cluster shows all cores for node; not only the number of cpus requested
elseif system_machine == "Slurm" #Sys.islinux()
  # For HPC
#   root_dir = homedir()
#   num_cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"])  # ERROR: LoadError: KeyError: key "SLURM_CPUS_PER_TASK" not found [from Princenton website]
#   println("No of cores = $num_cores")
  ntasks = parse(Int, ENV["SLURM_NTASKS"])  # gave 1, when cpus-per-task was passed by slurm
  # cores = parse(Int, ENV["SLURM_JOB_CPUS_PER_NODE"])  # if one node provided; else we have parse as follows
  subs = Dict("x"=>"*", "(" => "", ")" => "");
  cores = sum(eval(Meta.parse(replace(ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s]))))
  println("SLURM_NTASKS: $ntasks")
#   println("SLURM_JOB_CPUS_PER_NODE: $ntasks")
  println("SLURM Cores: $cores")
  addprocs(cores)  # ; exeflags="--project"subtract one becuase master already has one; but seems to work with higher number as well
#   addprocs()  # Dec 03, 2022 : testing new
else
    println("Must be windows or linux system, aborting")
    exit() # <- you can provide exit code if you want
end
# println("Number of procs: $(nprocs())")  
println("Number of processes: ", nprocs())
println("Number of workers: ", nworkers())
# Everywhere should come after cores are already added
@everywhere using Rasters
@everywhere include("Estimate_v54.jl")  #https://docs.julialang.org/en/v1/manual/code-loading/; evaluated in global scope


# base_folder = "$root_dir/Github/Blender"
# DataDir= "$base_folder/nc_files"
DataDir = "$base_folder/NoahMP"  # must exist

# Folder for saving outputs of run. out_subfolder can be passed as ARGS. Folder/subfolders will be created if non-existent
out_folder = "$DataDir/Runs/$out_subfolder"
println("Output_folder : $out_folder")

tmp_txtDir = "$out_folder/outputs_txt"    # To save text outputs for each pixel
nc_outDir = "$out_folder/outputs"         # To convert text outputs to netcdf file
# OutDir = "$DataDir/$out_subfolder/outputs_txt"  # To save text outputs for each pixel
# nc_exp_dir = "$DataDir/$out_subfolder/outputs_nc"  # To convert text outputs to netcdf file

# # 2. Read the Input netCDF file
# A = RasterStack("$base_folder/nc_files/inputs/merged.nc");  #merged_proj.nc
# For NoahMP
# A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_modscag.nc")  #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
# Following check are for prototyping only when running code locally, because I do not yet have NorthAmerica netcdf file
if occursin("STAFF-BY-M", host_machine)
    A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_cgf.nc")  #2016_clip_noahmp_cgf #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
elseif occursin("borg", host_machine)  # TODO: discover
    A = RasterStack("$DataDir/WY_merged/2016_noahmp_cgf.nc")  #2016_clip_noahmp_cgf #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
else
    A = RasterStack("$DataDir/WY_merged/2016_noahmp_cgf.nc")
    # println("Exiting code. Manually set A (rasterstack) around line 188")
    # exit(1)
end


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
# Extract cartesian index for non-missing data using one-day of data 
valid_pix_ind = findall(!ismissing, A["SWE_tavg"][Ti=1])
valid_pix_count = length(valid_pix_ind)
if start_idx > valid_pix_count
    # return nothing
    exit(0)
    println("Exiting without Run because start_idx is greater than valid_pix_idx", start_idx, valid_pix_count)
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

println("Non-missing pixel count = ", length(valid_pix_ind))
@info("Starting with loop: \n")
println("Processing ", length(ind), " out of ", length(valid_pix_ind))
# @info("Processing ", length(ind), " out of ", length(valid_pix_ind))
# exit(0)
""" General info about pixel counts
xind: 2336 , yind: 941 tind:364
Non-missing pixel count = 1011329
"""
# @sync makes the code wait for all processes to finish their part of the computation before moving on from the loop
# ie, without @sync, the on of the processors may move to next while loop is still running, creating error and crash whole script prematurely 
# Threads.@threads for pix in pixels
# @distributed for ij in ind
@sync @distributed for ij in ind
# for ij in ind
    i = ij[1]
    j = ij[2]
    # println("ij = ", ij)
    exp_dir = string("$tmp_txtDir/", "Pix_", i, "_", j)  # full path to folder for saving (temporary) text files   
    # if !isdir(exp_dir)  # process only if the pixel is not already processed
    # if !ispath(exp_dir * "/SWEpv.txt") # ispath is generic for both file and folder
    # if !isfile(exp_dir * "/SWEpv.txt") # This is better check: process the pixel only if the last file exported by optimzer (SWEpv.txt) does not yet exist
    if !isfile(exp_dir * "/out_vars.txt") # This is better check: process the pixel only if the last file exported by optimzer (SWEpv.txt) does not yet exist
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
        blender(exp_dir, WRFSWE, WRFP, WRFG, MSCF, AirT)
        # println("After calling BLENDER function")
        # here exp_dir = export directory for holding outputs text and log files (same as older version of code)
    end
end
end_time = time_ns()
running_time = (end_time - start_time)/1e9/60
println("Running Time (minutes) = $running_time")

#=
# Get index to iterate over
# Seems here x is rows and y are columns
# tind, yind, xind  = size(A)  # need a better way because order can be changed sometimes!; ie, use more explicit way of extracting x and y index
# This is more explicit
xind  = size(A, X)
yind  = size(A, Y)
tind  = size(A, Ti)
println("xind: ", xind, " , yind: ", yind, " tind:", tind)
# TODO: Here only outer loop is parallelzed, inner loop run sequentially, then wait for another outer loop to start
# https://github.com/JuliaLang/julia/issues/30343
start_time = time_ns()
# Threads.@threads for pix in pixels
@sync @distributed for i in 1:xind
# for i in 1:xind
    for j in 1:yind
        # println("i = ", i, " j = ", j)
        val = A[X=i, Y=j][:MODSCAG]  # pick any varaible: if one is missing, all will be missing
        # val = A[X=i, Y=j][:SWE_tavg]  # pick any varaible: if one is missing, all will be missing
        # println("val = ", val[1])
        if !ismissing(val[1])  # !isnan(val[1]) || 
            # print("Inside IF")
            # A_pt = A[X=i, Y=j]  # Get all variables for all time for one pixel
            # extract each of the input variables separately (trying to match of processing was done in prior version with text inputs)
            # but this approach is also useful if we decide to save each input netcdf file separately

            # WRFSWE = A[X=i, Y=j]["SWE_tavg"].data  # Here "data" is an AbstractArray.
            # WRFP = A[X=i, Y=j]["Snowf_tavg"].data
            # WRFG = A[X=i, Y=j]["Qg_tavg"].data
            # AirT = A[X=i, Y=j]["Tair_f_tavg"].data
            # MSCF = A[X=i, Y=j]["MODSCAG"].data;

            WRFSWE = A["SWE_tavg"][X=i, Y=j].data  # Here "data" is an AbstractArray.
            WRFP = A["Snowf_tavg"][X=i, Y=j].data
            WRFG = A["Qg_tavg"][X=i, Y=j].data
            AirT = A["Tair_f_tavg"][X=i, Y=j].data
            MSCF = A["MODSCAG"][X=i, Y=j].data;

            
            exp_dir = string("$tmp_txtDir/", "Pix_", i, "_", j)  # full path to folder for saving (temporary) text files
            # process only if the pixel is not already processed 
            if !isdir(exp_dir)
                # Call blender for the pixel. This is the only required line
                # rest of the codes are for houskeeking, preprocssing, post-processing
                blender(exp_dir, WRFSWE, WRFP, WRFG, MSCF, AirT)
                # println("After calling BLENDER function")
                # here exp_dir = export directory for holding outputs text and log files (same as older version of code)
            end
        end
    end
end
end_time = time_ns()
=#

# Step 4. PostProcessing: Combine text files into a grid and save as netcdf
function text2nc(var, idx)
    """Convert a text file to netcdf
    var: the name of variable (example SWE, Gmelt, G etc.)
    idx : the index where the variable of exist in the file [1 through 9]
    TODO: move this fuction inside the Estimate_v53.jl script
    TODO: seems like nc file produced here cannot by plotted using xarray/hvplot
    """
    # outRaster = copy(A[:MODSCAG])  # TODO: check for NoahMP; use other variable because it creates MODIS attrs
    outRaster = copy(A[:SWE_tavg])
    # Initial resulting array
    outRaster[:,:,:] .= missing
    # update the name of variable
    outRaster = rebuild(outRaster; name=var)  #:SWEhat
    for pix in pixels
        pix_xy = split(pix, "_")
        x = parse(Int16, pix_xy[2])
        y = parse(Int16, pix_xy[3])
        out_vars = readdlm("$tmp_txtDir/$pix/out_vars.txt")  # read whole combined file for that pixel
        # outRaster[X=x, Y=y] = readdlm("$tmp_txtDir/$pix/$var.txt");  # Older workflow
        outRaster[X=x, Y=y] = out_vars[:,idx]
    end
    # Save to nc file; 
    mkpath(nc_outDir)
    write("$nc_outDir/$var.nc", outRaster)
end

# without sync above one of the processor may move here causing error
pixels = readdir(tmp_txtDir)  # Danger: Error if we have outputs from prior runs that do not match current A dimensions
pixels = [pix for pix in pixels if startswith(pix, "Pix")];
# TODO (Nov 20, 2022): Count if all pixels are processed before invoking this following section below for converting text files to netcdf file
if length(pixels) == valid_pix_count  # length(valid_pix_ind)
    @info("Creating OUTPUT NETCDF FILES")
    # Call the function for creating netcdf for each of the text output files separately
    # Use the name of text file that was saved; Change the filename in main script if so required
    text2nc("SWE", 1)
    text2nc("Gmelt", 2)
    text2nc("G", 3)
    text2nc("Precip", 4)
    text2nc("Us", 5)
    text2nc("Gpv", 6)
    text2nc("Gmeltpv", 7)
    text2nc("Upv", 8)
    text2nc("SWEpv", 9);
else
    @info("All pixels not yet processed, so OUTPUT NETCDF FILES not yet created")
end
running_time = (end_time - start_time)/1e9/60
println("Total Running Time (minutes) = $running_time")
