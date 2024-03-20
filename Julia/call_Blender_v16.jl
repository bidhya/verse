"""
USAGE: Hardcoded with requirments for output_folder_YEAR, start_idx, and end_idx for running on Discover
- start and end indices are for y-axis only. All x's selected by default  
- julia call_Blender_v14.jl test2016 1 10  

SWE estimation using Blender algorithm (Prof. Mike Durand).
This is a helper script to organize inputs for calling the Blender function (Estimate_v57.jl).
Therefore, manually update input/output directories here for seperate runs (say NoahMP vs WRF etc)
Currently, only output_subfolder can optionally be passed as an argument.
Input/Output are currently saved in same main folder, but can easily be saved at separate locations [dig down below].

The script is currently capable of running on my following platforms.
- Windows machine
- Ubuntu (WSL2-based)
- OSC HPC
- Ohio-State Unity HPC 
Julia Versions
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
Oct 29, 2023 : Save nc files to temp_nc folder using call_Blender_v14.jl uses Estimate_v57.jl
Dec 03, 2023 : New call_Blender_v15.jl uses Estimate_v58.jl (modifications by Jack)

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
start_time = time_ns()

using Logging, LoggingExtras
using Tar
using Distributed  # otherwise everywhere macro won't work
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
        root_dir = "/discover/nobackup/projects/coressd"
        base_folder = "$root_dir/Blender"
        tmpdir =  ENV["LOCAL_TMPDIR"]  #tempdir()
    elseif occursin(".osc.edu", host_machine)
        root_dir = "/fs/ess/PAS1785/coressd"  # "/fs/scratch/PAS1785/coressd"
        base_folder = "$root_dir/Blender"
        tmpdir =  ENV["TMPDIR"]  #tempdir()
    elseif occursin("asc.ohio-state.edu", host_machine)  # .unity
        root_dir = "/fs/project/howat.4/yadav.111/coressd"  # homedir()  #  Unity
        # base_folder = "/home/yadav.111/Github/Blender"  # old
        base_folder = "$root_dir/Blender"  # "$root_dir/Github/coressd/Blender"
        tmpdir =  ENV["TMPDIR"]  #tempdir()
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
# addprocs()  # # Works but on cluster will use all cores, even though we not asked for. also with exclusive option will use all cores which may result in memory error.  
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

# base_folder = "$root_dir/Github/Blender"
# DataDir= "$base_folder/nc_files"
DataDir = "$base_folder/Inputs"  # must exist  (Old = NoahMP)
# # use only for HPC to copy input file here (hopefully for faster i/o operation)
# if occursin("borg", host_machine)
#     tmpdir =  ENV["LOCAL_TMPDIR"]  #tempdir()  
# else
#     # TODO: on discover if node does not start with "borg" the code can still jump here!
#     tmpdir =  ENV["TMPDIR"]  #tempdir()

# # Folder for saving outputs of run. out_subfolder can be passed as ARGS. Folder/subfolders will be created if non-existent
# # Oct 20, 2023: This out_folder is no more used in V14. but verify
# if occursin(".osc.edu", host_machine)
#     # out_folder = "/fs/ess/PAS1785/coressd/Blender/Runs/$out_subfolder" #  "$base_folder/Runs/$out_subfolder"  # "$DataDir/Runs/$out_subfolder"
#     out_folder = "/fs/scratch/PAS1785/coressd/Blender/Runs/$out_subfolder" #  "$base_folder/Runs/$out_subfolder"  # "$DataDir/Runs/$out_subfolder"
# elseif occursin("asc.ohio-state.edu", host_machine)  # .unity
#     out_folder = "$tmpdir/$out_subfolder"  # "$DataDir/Runs/$out_subfolder"
# else
#     out_folder = "$base_folder/Runs/$out_subfolder"  # "$DataDir/Runs/$out_subfolder"
# end
# @info("Output_folder : $out_folder")

# Make a folder insise HPC node because we want to copy existing files there
exp_dir = "outputs_txt_$(start_idx)_$(end_idx)"  #"$out_folder/outputs_txt"  # tmp_txtDir(old name) To save text outputs for each pixel
mkpath(exp_dir)  # mkdir
logDir = "logs_$(start_idx)_$(end_idx)"  # "logs"  "$out_folder/logs"   # Save logs in separate folder
# logDir = "$base_folder/Runs/$out_subfolder/logs/logs_$(start_idx)_$(end_idx)"  # For Debug: save on same folder and outputs  
mkpath(logDir)  # mkdir; must create here, else error in the current setup  
# Error on Discover (Nov 08, 2023): rm: cannot remove '/lscratch/tdirs/batch/slurm.24967584.byadav/logs': Directory not empty; Solution: #SBATCH --no-requeue

# cp("$base_folder/Runs/$out_subfolder/outputs_txt", "$tmp_txtDir/$water_year")  # copy to local machine; error if running the first time as this dir would not exist
# nc_outDir = "$out_folder/outputs"         # To convert text outputs to netcdf file
# nc_outDir = "$base_folder/Runs/$out_subfolder/outputs"
nc_outDir = "$base_folder/Runs/$out_subfolder/temp_nc/outputs_$(start_idx)_$(end_idx)"
# OutDir = "$DataDir/$out_subfolder/outputs_txt"  # To save text outputs for each pixel
# nc_exp_dir = "$DataDir/$out_subfolder/outputs_nc"  # To convert text outputs to netcdf file
# Oct 20, 2023: exit here if folder alreay exists
# if !ispath(exp_dir * "/SWEpv.txt") # ispath is generic for both file and folder
# if !isfile(exp_dir * "/SWEpv.txt")  # This is better check: process the pixel only if the last file exported by optimzer (SWEpv.txt) does not yet exist
if isdir(nc_outDir)  # process only if the pixel is not already processed
    @info("Exiting because following folder already processed: $nc_outDir")
    exit(0)
end

# # 2. Read the Input netCDF file
# A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_modscag.nc")  #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
# Following check are for prototyping only when running code locally, because I do not yet have NorthAmerica netcdf file
if occursin("L-JY0R5X3", host_machine) #|| test_run  # use this for testing as well
    @info("Test Run only.")
    # A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_cgf.nc", lazy=true)  # Error now (Aug 15, 2023). maybe it is old file/format
    A = RasterStack("$(DataDir)/WY_merged/$(water_year)_seup_modis.nc", lazy=true)
    # subset = A[X(Between(-120.5, -120)), Y(Between(60, 60.5))]  # 10 by 10 pixels; use small chip for prototyping
    # subset_fname = "$(DataDir)/WY_merged/subset_$(water_year)_seup_modis.nc"
    # # isfile(subset_fname) || write(subset_fname, subset)  # save. we need for post-processing analysis of results
    # write(subset_fname, subset)
    # A = RasterStack(subset_fname, lazy=true)  # read again (why). for consistency
# elseif occursin("borg", host_machine)  # TODO: discover
#     # A = RasterStack("$DataDir/WY_merged/2013_seup_modis.nc")  # 2016_noahmp_cgf 2016_clip_noahmp_cgf #, mappedcrs=EPSG(4326); for NoahMP with MODSCAG mapped to NoahMP resolution
#     A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
# elseif occursin(".osc.edu", host_machine)
#     # A = RasterStack("$DataDir/WY_merged/2016_clip_noahmp_cgf.nc")
#     A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
# elseif occursin("asc.ohio-state.edu", host_machine)
#     A = RasterStack("$DataDir/WY_merged/$water_year" * "_seup_modis.nc", lazy=true)
else
    # A = RasterStack("$DataDir/WY_merged/2016_seup_modis.nc")
    # A = RasterStack("$DataDir/WY_merged/" * water_year * "_seup_modis.nc")
    # copy to local Node (machine) on slurm
    cp("$DataDir/WY_merged/$water_year" * "_seup_modis.nc", "$tmpdir/$water_year" * "_seup_modis.nc", force=true)  #force=true will first remove an existing dst.
    # true required for discover when running on same node again after node_failure.
    A = RasterStack("$tmpdir/$water_year" * "_seup_modis.nc", lazy=true)  ## , lazy=true https://github.com/rafaqz/Rasters.jl/issues/449
end
# end_time = time_ns()
running_time = (time_ns() - start_time)/1e9/60
@info("Time until copying input netcdf to Node = $(round(running_time, digits=2)) minutes")

# A = RasterStack("$DataDir/WY_merged/2016_clip3.nc")  # for NoahMP with CGF MODIS
# Subset only the required variables because the nc file can have extraneous vars that cause problem with julia
A = A[(:Snowf_tavg, :SWE_tavg, :Tair_f_tavg, :Qg_tavg, :MODSCAG)];  # to remove spatial ref that cause problem during subsetting
""" General info about pixel counts
    size(A) = (2336, 941)
    xind: 2336 , yind: 941 tind:364
    Non-missing pixel count = 1011329
"""
szY = size(A, 2)  # get size in Y-dimension; here dim = 2
if end_idx > szY
    end_idx = szY
end
if test_run
    @info("TEST RUN ONLY  ")
    # # Aside: Get test set of data. This part needs re-writing before using for test [TODO].
    # A = A[X(Between(-130, -120)), Y(Between(60, 65))]
    # subset_fname = "$(DataDir)/WY_merged/subset_$(water_year)_seup_modis.nc"
    # write(subset_fname, A)
    # A = A[1:end, 100:102, :]
    A = A[1101:1120, start_idx:end_idx, :]  # change this as required for selecting a particular region/watershed to test on
else
    @info("DEFAULT RUN  ")
    A = A[1:end, start_idx:end_idx, :]  # use for default runs
end

# Now using cartesian index for iteration
# Extract cartesian index for non-missing (ie, all) data using one-day of data 
valid_pix_ind = findall(!ismissing, A["SWE_tavg"][Ti=1])
valid_pix_count = length(valid_pix_ind)
ind = valid_pix_ind
@info("Total valid pixel count = $(valid_pix_count)")
println("Total valid pixel count = $(valid_pix_count)")  # debug

# @everywhere using SharedArrays  # must come after addprocs(cores) else won't work for julia 1.10.0 version.  
# sizeA = size(A)  # to create sharedarrays
# println("sizeA = $(sizeA)") # debug
# # SWERaster = SharedArray{Float32}(10, 10, 366)
# SWERaster = SharedArray{Float32}(sizeA)
# SWERaster[:,:,:] .= NaN
# GmeltRaster = deepcopy(SWERaster)  # use deepcopy; simple copy will create general array, not a sharedarray, hence, won't work in parallel  
# GRaster = deepcopy(SWERaster)
# PrecipRaster = deepcopy(SWERaster)
# UsRaster = deepcopy(SWERaster)
# GpvRaster = deepcopy(SWERaster)
# GmeltpvRaster = deepcopy(SWERaster)
# UpvRaster = deepcopy(SWERaster)
# SWEpvRaster = deepcopy(SWERaster);

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
    WRFSWE = A["SWE_tavg"][X=i, Y=j].data  # Here "data" is an AbstractArray.
    WRFP = A["Snowf_tavg"][X=i, Y=j].data
    WRFG = A["Qg_tavg"][X=i, Y=j].data
    AirT = A["Tair_f_tavg"][X=i, Y=j].data
    MSCF = A["MODSCAG"][X=i, Y=j].data;
    # blender(out_folder, i, j, WRFSWE, WRFP, WRFG, MSCF, AirT)  # Call blender for the pixel. OLD.
    blender(i, j, WRFSWE, WRFP, WRFG, MSCF, AirT, logDir, exp_dir)
    # SWEhat, GmeltHat, Ghat,  Phat, Ushat, G_pv, Gmelt_pv, U_pv, SWEpv = blender(i, j, WRFSWE, WRFP, WRFG, MSCF, AirT, logDir, exp_dir)
    # SWERaster[i, j, :] = SWEhat
    # GmeltRaster[i, j, :] = GmeltHat
    # GRaster[i, j, :] = Ghat
    # PrecipRaster[i, j, :] = Phat
    # UsRaster[i, j, :] = Ushat
    # GpvRaster[i, j, :] = G_pv
    # GmeltpvRaster[i, j, :] = Gmelt_pv
    # UpvRaster[i, j, :] = U_pv
    # SWEpvRaster[i, j, :] = SWEpv;    
    # GC.gc()  # With GC, memory footprint is smaller, even though runtime was same
    # rand(Uniform(0,1)) < 0.1 && GC.gc();  # with 0.02 garbage collection only 2% of time; need Distributions package
    # cp(exp_dir, "$base_folder/Runs/$out_subfolder/outputs_txt", force=True)
end
# without sync above one of the processor to the next step (combining text files to netcdf) which will cause error
sleep(10)
@everywhere GC.gc()

end_time = time_ns()
running_time = (end_time - start_time)/1e9/3600
@info("Running Time (blender For Loop) = $(round(running_time, digits=2)) hours")

# # Oct 08, 2023. To save output nc files. Create Raster matching the input and save to nc file 
# mkpath(nc_outDir)
# outRaster = copy(A[:SWE_tavg])  # output raster template
# outRaster[:, :, :] = SWERaster
# SWERaster = rebuild(outRaster; name=:SWE)
# write("$nc_outDir/SWE.nc", SWERaster)

# outRaster = copy(A[:SWE_tavg])
# outRaster[:, :, :] = GmeltRaster
# GmeltRaster = rebuild(outRaster; name=:Gmelt)
# write("$nc_outDir/Gmelt.nc", GmeltRaster)

# outRaster = copy(A[:SWE_tavg])
# outRaster[:, :, :] = GRaster
# GRaster = rebuild(outRaster; name=:G)
# write("$nc_outDir/G.nc", GRaster)

# outRaster = copy(A[:SWE_tavg])
# outRaster[:, :, :] = PrecipRaster
# PrecipRaster = rebuild(outRaster; name=:Precip)
# write("$nc_outDir/Precip.nc", PrecipRaster)
# GC.gc()

# outRaster = copy(A[:SWE_tavg])
# outRaster[:, :, :] = UsRaster
# UsRaster = rebuild(outRaster; name=:Us)
# write("$nc_outDir/Us.nc", UsRaster)

# outRaster = copy(A[:SWE_tavg])
# outRaster[:, :, :] = GpvRaster
# GpvRaster = rebuild(outRaster; name=:Gpv)
# write("$nc_outDir/Gpv.nc", GpvRaster)

# outRaster = copy(A[:SWE_tavg])
# outRaster[:, :, :] = GmeltpvRaster
# GmeltpvRaster = rebuild(outRaster; name=:Gmeltpv)
# write("$nc_outDir/Gmeltpv.nc", GmeltpvRaster)

# outRaster = copy(A[:SWE_tavg])
# outRaster[:, :, :] = UpvRaster
# UpvRaster = rebuild(outRaster; name=:Upv)
# write("$nc_outDir/Upv.nc", UpvRaster)

# outRaster = copy(A[:SWE_tavg])
# outRaster[:, :, :] = SWEpvRaster
# SWEpvRaster = rebuild(outRaster; name=:SWEpv)
# write("$nc_outDir/SWEpv.nc", SWEpvRaster)

# exit(0)  # exit here to avoid running 2nd part (text2nc)

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
    # outRaster = copy(A[:MODSCAG])  # TODO: check for NoahMP; use other variable because it creates MODIS attrs
    # outRaster = copy(A[:SWE_tavg])
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
        # out_vars = readdlm("$exp_dir/$pix/out_vars.txt")  # $(exp_dir)/Pix_$(i)_$(j).txt   read whole combined file for that pixel
        # out_vars = readdlm("$tmpdir/outputs_txt/$pix/out_vars.txt")  # read whole combined file for that pixel
        # outRaster[X=x, Y=y] = readdlm("$exp_dir/$pix/$var.txt");  # Older workflow
        outRaster[X=x, Y=y] = out_vars[:,idx]
    end
    # Save to nc file; 
    # mkpath(nc_outDir)  # error in distributed
    write("$nc_outDir/$var.nc", outRaster)  # Aug 08, 2023: Error perhaps due to updates to raster/ncdataset/etc. (error tested on windows and wsl2)
end

sleep(10)
# Count if all pixels are processed before calling the following section for converting text files to netcdf file
pixels = readdir(exp_dir)
pixels = [pix for pix in pixels if startswith(pix, "Pix")];
processed_pix_count = length(pixels)
@info("Pixels processsed = $processed_pix_count out of $valid_pix_count")
if processed_pix_count == valid_pix_count #&& system_machine == "Windows"
    outRaster = copy(A[:SWE_tavg])
    var_idx_tuple = (("SWE", 1), ("Gmelt", 2), ("G", 3), ("Precip", 4), ("Us", 5), ("Gpv", 6), ("Gmeltpv", 7), ("Upv", 8), ("SWEpv", 9))
    mkpath(nc_outDir)
    # Threads.@threads (Roughly same running for threaded and non-thread version on Discover. Likely due to I/O bottleneck)
    # Threads.@threads for var_idx in var_idx_tuple
    @sync @distributed for var_idx in var_idx_tuple
        var_name = var_idx[1]
        var_idx = var_idx[2]
        text2nc(var_name, var_idx, outRaster, exp_dir, nc_outDir);  # text2nc("SWE", 1, outRaster)
    end
else
    @info("All pixels not yet processed, so OUTPUT NETCDF FILES not yet created")
    # TODO list out what pixels were not processed (and why? at a later point)
    # if system_machine == "Windows"
    #     Tar.create(exp_dir, "$(exp_dir).tar.gz")
    # else
    #     Tar.create(exp_dir, pipeline(`gzip -9`, "$(exp_dir).tar.gz"))
    # end
    # mkpath("$base_folder/Runs/$out_subfolder/outputs_txt")
    # sleep(2)
    # cp("$(exp_dir).tar.gz", "$base_folder/Runs/$out_subfolder/outputs_txt/$(exp_dir).tar.gz", force=true)  #ERROR: LoadError: IOError: sendfile: Unknown system error -122 (Unknown system error -122)
    # @info("Moved outputs_txt: $base_folder/Runs/$out_subfolder/outputs_txt/$(exp_dir).tar.gz")
end
@info("Finished: Combining Text Ouputs to NetCDF FILES")
end_time = time_ns()
running_time = (end_time - start_time)/1e9/3600
@info("Running Time (convert text to nc) = $(round(running_time, digits=3)) hours")

mkpath("$base_folder/Runs/$out_subfolder/logs")  # also works if path already exists

# Altert : Creating Tar (large (>3GB)) file and/or copying them using Julia is unstable. Hence, this code is commeted out. Use bash/slurm script for this.  
# # Tar move log files from compute node to permanent location
# if system_machine == "Windows"
#     # TODO file not really .gz on windows but works. Might be a good idea to remore .gz suffix for windows files
#     @info("Create tar on Windows")
#     Tar.create(logDir, "$(logDir).tar.gz")  # no native compression for tar in julia and no pipeline in windows
# else
#     @info("Create tar on Slurm node")
#     Tar.create(logDir, pipeline(`gzip -9`, "$(logDir).tar.gz"))  # compress with gzip through pipeline in linux
#     # TODO ERROR: LoadError: unsupported file type: "logs_1_110/Pix_1112_37.txt"
#     # Undiagnosed error for now but keep an eye if it occurs again on Discover.  
# end
# sleep(2)
# cp("$(logDir).tar.gz", "$base_folder/Runs/$out_subfolder/logs/$(logDir).tar.gz", force=true)  #TODO: ERROR: LoadError: IOError: sendfile: Unknown system error -175626013 (Unknown system error -175626013)
# # Note: program crash at this point is OK. logDir will be copied so we don't loose anything.
# @info("Moved logfile to : $base_folder/Runs/$out_subfolder/logs/$(logDir).tar.gz")

end_time = time_ns()
running_time = (end_time - start_time)/1e9/60  # minutes
if running_time < 60
    @info("Grant Total Running Time = $(round(running_time, digits=2)) minutes")
else
    running_time = running_time/60 # convert to hours
    @info("Grant Total Running Time = $(round(running_time, digits=2)) hours")
end
