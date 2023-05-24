"""
Create NetCDF files from the temporary text files (for HPC only)  
- Slurm job submitted from "submit_blender_job.py" python script
- Creating one nc file at a time (because the it it too time consuming)

USAGE:   julia combine_txt2csv.jl sub_folder var_name var_index
example: julia combine_txt2csv.jl new_test1 SWE 1

This function should only be required on HPC where large number of pixels
    are being processed. We will only have partial outputs on such situation. 
    This script should thus be called at the end to combine text files to csv file.
    Processing of whole NA is so time consuming (~6 hours), with this approach we can process
    each of the nine ouput variables (NetCDF files) as a separate job
Still Hardcoded in script
- root_dir
- base_folder
- DataDir = %base_folder/NoahMP
- out_folder = %DataDir/Runs/%out_subfolder
- tmp_txtDir = %out_folder/outputs_txt" # Fullpath to processed text files
- nc_outDir = %out_folder/outputs"      # Fullpath to save output NETCDF files; created here
- A = RasterStack(%DataDir/WY_merged/2016_noahmp_cgf.nc)
    
"""
# DataDir= "/users/PAS0503/jldechow/foursix/AllTuolomne46/Pix3252"
# DataDir = ARGS[1]
# cores = ARGS[1]
# cores = parse(Int8, cores)  # maximum of Int8 is 127 only

# # This start, end index is no more useful because we are processing from netcdf rather than text files separately
# start_idx = parse(Int16, ARGS[1])
# end_idx = parse(Int16, ARGS[2])
# Base.ARGS :: Vector{String} >>>>> An array of the command line arguments passed to Julia, as strings.

# arg_len = length(ARGS)
# if arg_len == 0
#     # no argument passed from command line
#     out_subfolder = "NA2"  # "NoahMP_CGF"
# else
#     out_subfolder = ARGS[1]     # output subfolder relative to input files; temp text and nc_outputs saved here
# end
out_subfolder = ARGS[1]         # output subfolder relative to input files; temp text and nc_outputs saved here
var_name = ARGS[2]              # name of the output variable (SWE for example)
var_idx = parse(Int8, ARGS[3])  # index of variable [1, 2, ..., 9]. Follow strictly as the file was created during Blender run

using DelimitedFiles
using Rasters

log_filename = string("txt2nc_", out_subfolder, "_", var_name, "_", var_idx, ".log")  #construct a log filename
# Updated logger for finer control; still not working with distributed
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


# 1. Setup inputs and output directories and file locations
# select the root directory (this will be different on windows, linux or Mac) 
host_machine = gethostname()
println("Host computer machine: $host_machine")
# if occursin("STAFF-BY-M", host_machine)
if occursin("L-JY0R5X3", host_machine)
    if Sys.iswindows()
        root_dir = "D:"  #for windows machine
    elseif Sys.islinux()
        root_dir = "/mnt/d"  #for Ubuntu (WSL2 on windows machine
    else
        println("Unknown system on windows, manually add root directory before proceeding. Exiting code")
        exit(1)    
    end
    base_folder = "$root_dir/coressd/Blender"
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
    root_dir = "/fs/project/howat.4/yadav.111/coressd"  # homedir()  #  Unity
    # base_folder = "/home/yadav.111/Github/Blender"  # old
    base_folder = "$root_dir/Blender"  # "$root_dir/Github/coressd/Blender"
    system_machine = "Slurm"
else
    println("Unknown computer, manually add root directory before proceeding. Exiting code")
    exit(1)  # if error, comment this line, add root and base_folder below and run again
end
println("base_folder : $base_folder")

# base_folder = "$root_dir/Github/Blender"
# DataDir= "$base_folder/nc_files"
DataDir = "$base_folder/NoahMP"  # must exist
# Folder for saving outputs of run. out_subfolder can be passed as ARGS. Folder/subfolders will be created if non-existent
out_folder = "$base_folder/Runs/$out_subfolder"  # "$DataDir/Runs/$out_subfolder"
println("Output_folder : $out_folder")

tmp_txtDir = "$out_folder/outputs_txt" # Fullpath to processed text files
nc_outDir = "$out_folder/outputs1"      # Fullpath to save output NETCDF files

# 2. Read the original Input netCDF file
# # 2. Read the Input netCDF file
# A = RasterStack("$base_folder/nc_files/inputs/merged.nc");  #merged_proj.nc
# For NoahMP
# A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_modscag.nc")  #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
# Following check are for prototyping only when running code locally, because I do not yet have NorthAmerica netcdf file
if occursin("L-JY0R5X3", host_machine)  # STAFF-BY-M
    A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_cgf.nc")  #2016_clip_noahmp_cgf #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
elseif occursin("borg", host_machine)  # TODO: discover
    A = RasterStack("$DataDir/WY_merged/2016_seup_modis.nc")  # 2016_noahmp_cgf 2016_clip_noahmp_cgf #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
elseif occursin(".osc.edu", host_machine)
    A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_cgf.nc")
else
    A = RasterStack("$DataDir/WY_merged/2016_seup_modis.nc")
    # A = RasterStack("$DataDir/WY_merged/2016_noahmp_cgf.nc")
    # A = RasterStack("$DataDir/WY_merged/ak_polar_fix.nc")
    # A = RasterStack("$DataDir/WY_merged/ak_polar_fix_no.nc")
    # println("Exiting code. Manually set A (rasterstack) around line 188")
    # exit(1)
end


# A = RasterStack("$DataDir/WY_merged/2016_clip3.nc")  # for NoahMP with CGF MODIS
# Subset only the required variables because the nc file can have extraneous vars that cause problem with julia
A = A[(:Snowf_tavg, :SWE_tavg, :Tair_f_tavg, :Qg_tavg, :MODSCAG)];  # to remove spatial ref that cause problem during subsetting

start_time = time_ns()
# Extract cartesian index for non-missing data using one-day of data 
valid_pix_ind = findall(!ismissing, A["SWE_tavg"][Ti=1])
valid_pix_count = length(valid_pix_ind)

# println("Non-missing pixel count = ", length(valid_pix_ind))
# @info("Starting with loop: \n")
# println("Processing ", length(ind), " out of ", length(valid_pix_ind))
# if !isfile(exp_dir * "/out_vars.txt") # This is better check: process the pixel only if the last file exported by optimzer (SWEpv.txt) does not yet exist

# Step 4. PostProcessing: Combine text files into a grid and save as netcdf
function text2nc(var, idx)
    """Convert a text file to netcdf
    var: the name of variable (example SWE, Gmelt, G etc.)
    idx : column index of variable [1 through 9]
    TODO: move this fuction inside the Estimate_v53.jl script
    TODO: seems like nc file produced here cannot by plotted using xarray/hvplot
    """
    if !isfile("$nc_outDir/$var.nc")
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
    else
        @info("$var already process, skipping ... ")
    end
end

pixels = readdir(tmp_txtDir)  # Danger: Error is we have outputs from prior runs that do not match current A dimensions
pixels = [pix for pix in pixels if startswith(pix, "Pix")];
# println("Processed pix count = ", length(pixels), "  Total pix count = ", valid_pix_count)
# TODO (Nov 20, 2022): Count if all pixels are processed before invoking this following section below for converting text files to netcdf file
if length(pixels) == valid_pix_count  # length(valid_pix_ind)
    @info("Creating OUTPUT NETCDF FILES: $var_name : $var_idx")
    # # Call the function for creating netcdf for each of the text output files separately
    # # Use the name of text file that was saved; Change the filename in main script if so required
    # text2nc("SWE", 1)
    # text2nc("Gmelt", 2)
    # text2nc("G", 3)
    # text2nc("Precip", 4)
    # text2nc("Us", 5)
    # text2nc("Gpv", 6)
    # text2nc("Gmeltpv", 7)
    # text2nc("Upv", 8)
    # text2nc("SWEpv", 9);
    text2nc(var_name, var_idx);

else
    # @info("OUTPUT NETCDF file not create because all the pixels not yet processed")
    @info("Creating NETCDF with incomplete run")
    text2nc(var_name, var_idx);
end
end_time = time_ns()
running_time = (end_time - start_time)/1e9/60
println("Total Running Time (minutes) = $running_time")
