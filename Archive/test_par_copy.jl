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
# log_filename = string(start_idx, "_", end_idx, ".log")  #construct a log filename

logger = FormatLogger(open("/users/PZS0724/bidhya/Github/slurm_jobs/mylog.log", "w"), 0) do io, args
    # Write the module, level and message only
    println(io, args._module, " | ", "[", args.level, "] ", args.message)
end
info_only_logger = MinLevelLogger(logger, Logging.Info);  # Logging.Error
global_logger(info_only_logger)  # Set the global logger to logger; else logs doing directly to console

@info("base_folder : $base_folder")
@info("Host computer machine: $host_machine")

DataDir = "$base_folder/NoahMP"  # must exist
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
pixels = readdir(exp_dir)  # Danger: Error if we have outputs from prior runs that do not match current A dimensions
pixels = [pix for pix in pixels if startswith(pix, "Pix")];
pixels = pixels[1:2000]
println(length(pixels))
start_time = time_ns()

if system_machine == "Slurm"
    # copy text files to node; but this seems to increase runtime on OSC; so comment next two lines to leave txt_files on original location  
    t1 = time_ns()
    test_dir = "$out_folder/copy_txt/outputs_txt"    # tmp_txtDir(old name) To save text outputs for each pixel
    mkpath(test_dir)
    Threads.@threads for pix in pixels  # 
        cp("$exp_dir/$pix", "$test_dir/$(pix)", force=true)
        sleep(1)
        # exp_dir = "$tmpdir/outputs_txt"
    end
    t2 = time_ns()
    running_time = (t2 - t1)/1e9/3600  # 5 hours on Discover just to copy files. Thus, remove this copy part; But once files moved to tmpdir; creation of nc files is very fast (3 mins vs 30 mins (without copying))  
    @info("Total to copy text time files to tmpdir (hours) = $(round(running_time, digits=2))")
end
end_time = time_ns()
running_time = (end_time - start_time)/1e9/60  # minutes
if running_time < 60
    @info("Grant Total Running Time = $(round(running_time, digits=2)) minutes")
else
    running_time = running_time/60 # convert to hours
    @info("Grant Total Running Time = $(round(running_time, digits=2)) hours")
end
