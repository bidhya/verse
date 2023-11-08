"""
USAGE: Hardcoded with requirments for output_folder_YEAR, start_idx, and end_idx for running on Discover
- julia combine_nc_files.jl V14_WY2015
Oct 29, 2023
Bidhya N Yadav

To combine small parts of nc files into a PAN-america-wide nc file for all the Blender output variables
Slurm hint: 20 tasks with ~90 GB Ram
parallelized using threads  
"""

arg_len = length(ARGS)
out_subfolder = ARGS[1]  # WY2016. output subfolder relative to input files; temp_text and nc_outputs saved here
water_year = out_subfolder[end-3:end] #last 4 chars are assumed year, else error. out_subfolder[3:end]

using Logging, LoggingExtras
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

log_filename = string(ENV["SLURM_SUBMIT_DIR"], "/combine_nc.log")
# log_filename = "combine_nc.log"
logger = FormatLogger(open(log_filename, "w"), 0) do io, args
    # Write the module, level and message only
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end
info_only_logger = MinLevelLogger(logger, Logging.Info);  # Logging.Error
global_logger(info_only_logger)  # Set the global logger to logger; else logs doing directly to console

@info("base_folder : $base_folder")
@info("Host computer machine: $host_machine")

# Everywhere should come after cores are already added
using Rasters, NCDatasets, ArchGDAL

# base_folder = "$root_dir/Github/Blender"
# DataDir= "$base_folder/nc_files"
DataDir = "$base_folder/NoahMP"  # must exist

# # # 2. Read the Input netCDF file
# # A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_modscag.nc")  #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
# # Following check are for prototyping only when running code locally, because I do not yet have NorthAmerica netcdf file
# if occursin("L-JY0R5X3", host_machine) #|| test_run  # use this for testing as well
#     @info("Test Run only.")
#     # A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_cgf.nc", lazy=true)  # Error now (Aug 15, 2023). maybe it is old file/format
#     A = RasterStack("$(DataDir)/WY_merged/$(water_year)_seup_modis.nc", lazy=true)
#     # subset = A[X(Between(-120.5, -120)), Y(Between(60, 60.5))]  # 10 by 10 pixels; use small chip for prototyping
#     # subset_fname = "$(DataDir)/WY_merged/subset_$(water_year)_seup_modis.nc"
#     # # isfile(subset_fname) || write(subset_fname, subset)  # save. we need for post-processing analysis of results
#     # write(subset_fname, subset)
#     # A = RasterStack(subset_fname, lazy=true)  # read again (why). for consistency
# # elseif occursin("borg", host_machine)  # TODO: discover
# #     # A = RasterStack("$DataDir/WY_merged/2013_seup_modis.nc")  # 2016_noahmp_cgf 2016_clip_noahmp_cgf #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
# #     A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
# # elseif occursin(".osc.edu", host_machine)
# #     # A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_cgf.nc")
# #     A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
# # elseif occursin("asc.ohio-state.edu", host_machine)
# #     A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
# else
#     # A = RasterStack("$DataDir/WY_merged/2016_seup_modis.nc")
#     # A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
#     cp("$DataDir/WY_merged/$water_year" * "_seup_modis.nc", "$tmpdir/$water_year" * "_seup_modis.nc")  # copy to local machine
#     A = RasterStack("$tmpdir/$water_year" * "_seup_modis.nc")  ## , lazy=true https://github.com/rafaqz/Rasters.jl/issues/449
#     # if test_run
#     #     # Aside: Get test set of data
#     #     A = A[X(Between(-130, -120)), Y(Between(60, 65))]
#     #     subset_fname = "$(DataDir)/WY_merged/subset_$(water_year)_seup_modis.nc"
#     #     write(subset_fname, A)
#     # end
# end

start_time = time_ns()

A = RasterStack("$(DataDir)/WY_merged/$(water_year)_seup_modis.nc")
A = A[(:Snowf_tavg, :SWE_tavg, :Tair_f_tavg, :Qg_tavg, :MODSCAG)];  # to remove spatial ref that cause problem during subsetting
# Create output rasters based in input raster template
outRasterTemplate = copy(A[:SWE_tavg])
outRasterTemplate[:,:,:] .= missing 

A = nothing
GC.gc()

# SWERaster = copy(A[:SWE_tavg])
# SWERaster[:,:,:] .= missing
# GmeltRaster = deepcopy(SWERaster)  # use deepcopy; simple copy will create general array, not a sharedarray, hence, won't work in parallel  
# GRaster = deepcopy(SWERaster)
# PrecipRaster = deepcopy(SWERaster)
# UsRaster = deepcopy(SWERaster)
# GpvRaster = deepcopy(SWERaster)
# GmeltpvRaster = deepcopy(SWERaster)
# UpvRaster = deepcopy(SWERaster)
# SWEpvRaster = deepcopy(SWERaster);


nc_outDir = "$base_folder/Runs/$out_subfolder/temp_nc"  # /outputs_$(start_idx)_$(end_idx)

function combine_nc(nc_outDir, var)
    # Combine nc files
    outRaster = deepcopy(outRasterTemplate)

    folders = readdir(nc_outDir)
    # TODO Sort these files in ascending ordering starting at index 1.
    sym_var = Symbol(var) # uppercase(var)
    # Threads.@threads for folder in folders
    for folder in folders
            B = RasterStack("$(nc_outDir)/$(folder)/$(var).nc", lazy=false)        
        outRaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[sym_var]  # B[:SWE]
        sleep(1)  # this may help with unknown random error
    end
    return outRaster
end

final_nc_outDir = "$base_folder/Runs/$out_subfolder/outputs"
mkpath(final_nc_outDir)
out_vars = ("SWE", "Gmelt", "G", "Precip", "Us", "Gpv", "Gmeltpv", "Upv", "SWEpv")
for var in out_vars
    # var = "SWE"
    if !isfile("$(final_nc_outDir)/$(var).nc")
        @info("Processing nc_file: $var")
        outRaster = combine_nc(nc_outDir, var)
        sym_var = Symbol(var)
        outRaster = rebuild(outRaster; name=sym_var)
        write("$(final_nc_outDir)/$(var).nc", outRaster, force=true)
        # outRaster = Nothing
        # GC.gc()
    else
        @info("$var  NC file alreay processed before, skipping")
    end
end
# exit(0)


# folders = readdir(nc_outDir)
# Threads.@threads for folder in folders
#     @info("Processing nc_file: $folder")

#     B = RasterStack("$(nc_outDir)/$(folder)/SWE.nc", lazy=true)
#     SWERaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[:SWE]

#     B = RasterStack("$(nc_outDir)/$(folder)/Gmelt.nc", lazy=true)
#     GmeltRaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[:Gmelt]

#     B = RasterStack("$(nc_outDir)/$(folder)/G.nc", lazy=true)
#     GRaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[:G]

#     B = RasterStack("$(nc_outDir)/$(folder)/Precip.nc", lazy=true)
#     PrecipRaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[:Precip]
    
#     B = RasterStack("$(nc_outDir)/$(folder)/Us.nc", lazy=true)
#     UsRaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[:Us]
    
#     B = RasterStack("$(nc_outDir)/$(folder)/Gpv.nc", lazy=true)
#     GpvRaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[:Gpv]
    
#     B = RasterStack("$(nc_outDir)/$(folder)/Gmeltpv.nc", lazy=true)
#     GmeltpvRaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[:Gmeltpv]
    
#     B = RasterStack("$(nc_outDir)/$(folder)/Upv.nc", lazy=true)
#     UpvRaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[:Upv]
    
#     B = RasterStack("$(nc_outDir)/$(folder)/SWEpv.nc", lazy=true)
#     SWEpvRaster[At(lookup(B, X)), At(lookup(B, Y)), :] = B[:SWEpv]    
# end

# nc_outDir = "$base_folder/Runs/$out_subfolder/outputs"
# mkpath(nc_outDir)
# @info("nc_outDir : $nc_outDir")

# SWERaster = rebuild(SWERaster; name=:SWE)
# write("$nc_outDir/SWE.nc", SWERaster)
# GmeltRaster = rebuild(GmeltRaster; name=:Gmelt)
# write("$nc_outDir/Gmelt.nc", GmeltRaster)
# GRaster = rebuild(GRaster; name=:G)
# write("$nc_outDir/G.nc", GRaster)
# PrecipRaster = rebuild(PrecipRaster; name=:Precip)
# write("$nc_outDir/Precip.nc", PrecipRaster)
# UsRaster = rebuild(UsRaster; name=:Us)
# write("$nc_outDir/Us.nc", UsRaster)
# GpvRaster = rebuild(GpvRaster; name=:Gpv)
# write("$nc_outDir/Gpv.nc", GpvRaster)
# GmeltpvRaster = rebuild(GmeltpvRaster; name=:Gmeltpv)
# write("$nc_outDir/Gmeltpv.nc", GmeltpvRaster)
# UpvRaster = rebuild(UpvRaster; name=:Upv)
# write("$nc_outDir/Upv.nc", UpvRaster)
# SWEpvRaster = rebuild(SWEpvRaster; name=:SWEpv)
# write("$nc_outDir/SWEpv.nc", SWEpvRaster)

end_time = time_ns()
# running_time = (end_time - start_time)/1e9/3600
running_time = (end_time - start_time)/1e9/60  # minutes
if running_time < 60
    @info("Grant Total Running Time = $(round(running_time, digits=2)) minutes")
else
    running_time = running_time/60 # convert to hours
    @info("Grant Total Running Time = $(round(running_time, digits=2)) hours")
end
