"""
USAGE: Hardcoded with requirments for output_folder_YEAR, start_idx, and end_idx for running on Discover
- julia combine_nc_files.jl WY2015 010
Oct 29, 2023
Bidhya N Yadav

To combine small parts of nc files into a PAN-america-wide nc file for all the Blender output variables
Slurm hint: 20 tasks with ~90 GB Ram
parallelized using threads  
June 09, 2024: removed parallelization due to memory limitation for 1 km run. 
"""

arg_len = length(ARGS)
out_subfolder = ARGS[1]  # WY2016. output subfolder relative to input files; temp_text and nc_outputs saved here
water_year = out_subfolder[end-3:end] #last 4 chars are assumed year, else error. out_subfolder[3:end]
RES = ARGS[2] # "050"  # Grid Resolution: 010 = 0.01 degree; 050 = 0.050 degree (5km); 025 = 0.025 degree etc


using Logging, LoggingExtras, Dates
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

log_filename = string(ENV["SLURM_SUBMIT_DIR"], "/out/combine_ncfiles.log")
# log_filename = "out/combine_ncfiles.log"
logger = FormatLogger(open(log_filename, "w"), 0) do io, args
    # Write the module, level and message only
    println(io, args._module, " | ", "[", args.level, "] ", args.message, " | ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
end
info_only_logger = MinLevelLogger(logger, Logging.Info);  # Logging.Error
global_logger(info_only_logger)  # Set the global logger to logger; else logs doing directly to console

@info("base_folder : $base_folder")
@info("Host computer machine: $host_machine")

# Everywhere should come after cores are already added
using Rasters, NCDatasets #, ArchGDAL

# base_folder = "$root_dir/Github/Blender"
# DataDir= "$base_folder/nc_files"
DataDir = "$base_folder/Inputs_$(RES)"  # must exist (Old = NoahMP)
# DataDir = "$root_dir/projects/coressd/Blender/Inputs_$(RES)"  # INDIR. must exist  (Old = NoahMP)

start_time = time_ns()

# A = RasterStack("$(DataDir)/WY_merged/$(water_year)_seup_modis.nc", lazy=true)  # withoout lazy=true, NODE_FAIL error. Likely by loading all data in memory  
# A = A[(:Snowf_tavg, :SWE_tavg, :Tair_f_tavg, :Qg_tavg, :MODSCAG)];  # to remove spatial ref that cause problem during subsetting
files = (
    "$DataDir/lis/WY$(water_year)/SCF.nc",
    "$DataDir/lis/WY$(water_year)/Qg_tavg.nc"
    )
A = RasterStack(files; lazy=true)
# A = Raster("$(DataDir)/lis/WY$(water_year)/MODSCAG.nc", lazy=true)  # 64 bit float hence output was big.
# A = Raster("$(DataDir)/lis/WY$(water_year)/SWE_tavg.nc", lazy=true)  # this is 32 bit but may not work
# Create output rasters based in input raster template
outRasterTemplate = copy(A[:Qg_tavg])  # only for rasterstack
# outRasterTemplate = copy(A)  # for Raster. TODO: can use A directly without this copying part because now we have just one raster.
outRasterTemplate[:,:,:] .= missing 

A = nothing
GC.gc()

nc_outDir = "$base_folder/Runs/$(RES)/$out_subfolder/temp_nc"  # /outputs_$(start_idx)_$(end_idx)

function combine_nc(nc_outDir, var)
    # Combine nc files
    outRaster = deepcopy(outRasterTemplate)
    folders = readdir(nc_outDir)
    # TODO Sort these files in ascending ordering starting at index 1.
    sym_var = Symbol(var) # uppercase(var)
    # Threads.@threads for folder in folders  # TODO: June 10, 2024: Investigate error when using thread.
    for folder in folders
            B = RasterStack("$(nc_outDir)/$(folder)/$(var).nc", lazy=false)  # TODO: change to Raster
        outRaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[sym_var]  # B[:SWE]
        sleep(1)  # this may help with unknown random error
    end
    return outRaster
end

final_nc_outDir = "$base_folder/Runs/$(RES)/$out_subfolder/outputs"  # TODO local_temp_dir then tse_tempdir
mkpath(final_nc_outDir)
out_vars = ("SWE", "Gmelt", "G", "Precip", "Us", "Gpv", "Gmeltpv", "Upv", "SWEpv")
# out_vars = ("SWE", "Gmelt")  # for testing only

for var in out_vars
    # var = "SWE"
    if !isfile("$(final_nc_outDir)/$(var).nc")
        @info("Combine nc files: $var")
        outRaster = combine_nc(nc_outDir, var)
        sym_var = Symbol(var)
        @info("Rebuild raster.")
        outRaster = rebuild(outRaster; name=sym_var)
        @info("Save netcdf")
        write("$(final_nc_outDir)/$(var).nc", outRaster, force=true)
        @info("Finished saving netcdf file.")
        println(logger.stream, "")  # include blank line in the output log after logging this message
        # outRaster = Nothing
        # GC.gc()
    else
        @info("$var  NC file alreay processed before, skipping")
    end
end

end_time = time_ns()
# running_time = (end_time - start_time)/1e9/3600
running_time = (end_time - start_time)/1e9/60  # minutes
if running_time < 60
    @info("Grant Total Running Time = $(round(running_time, digits=2)) minutes")
else
    running_time = running_time/60 # convert to hours
    @info("Grant Total Running Time = $(round(running_time, digits=2)) hours")
end
