"""
USAGE: pass Water Year and stepsize (number of rows)
- julia ../verse/Julia/submit_slurm.jl 2016 35

Create and submit Blender job using Julia

"""
arg_len = length(ARGS)
out_subfolder = ARGS[1]  # WY2016. output subfolder relative to input files; temp_nc saved here
water_year = out_subfolder[end-3:end] #last 4 chars are assumed year, else error
# step is the count of rows (y-direction) of netcdf file to process in one job/run 
step = ARGS[2] #|| 35 on discover
step = parse(Int32, step)
memory = "184gb"

# mkpath(logDir)
mkpath("slurm_jobs/$(water_year)/.out")
cd("slurm_jobs/$(water_year)")

# 1. Setup inputs and output directories and file locations
# select the root directory (this will be different on windows, linux or Mac) 
host_machine = gethostname()
if occursin("discover", host_machine) #|| occursin("borg", host_machine)
    root_dir = "/discover/nobackup/projects/coressd"
    base_folder = "$root_dir/Blender"
    hpc_name = "discover"
    cores = 46  # use 45 so one core can be used to monitor run using srun/htop.
    memory = "0" #"184gb"
elseif occursin(".osc.edu", host_machine)
    root_dir = "/fs/ess/PAS1785/coressd"  # "/fs/scratch/PAS1785/coressd"
    base_folder = "$root_dir/Blender"
    hpc_name = "osc"
    cores = 40
    memory = "0" #"175gb"
elseif occursin("asc.ohio-state.edu", host_machine)  # .unity
    root_dir = "/fs/project/howat.4/yadav.111/coressd"  # homedir()  #  Unity
    # base_folder = "/home/yadav.111/Github/Blender"  # old
    base_folder = "$root_dir/Blender"  # "$root_dir/Github/coressd/Blender"
    hpc_name = "unity"
    cores = 40  # 39 24 cores with 96GB memory for old node
    memory = "186gb"  # 180gb seems 186 max allowed
else
    @info("Unknown computer, manually add root directory before proceeding. Exiting code")  # will output directly to console, ie like print statement
    exit(1)  # if error, comment this line, add root and base_folder below and run again
end

function create_job(hpc, jobname, cores, memory, runtime, out_subfolder, start_idx, end_idx, valid_pix_count, begin_delay)
    # job_file = "$(pwd())/$(jobname).job"  # this also works same
    job_file = "$(jobname).job"  # this looks much cleaner because files will created relative to where the script is called from.
    open(job_file, "w") do f
        # write(f, "A, B, C, D\n")
        write(f, "#!/usr/bin/env bash\n\n")
        # fh.writelines("#SBATCH --account=s2701\n")  # for DISCOVER
        # fh.writelines('#SBATCH --constraint="[cas|sky]"\n')
        if hpc == "discover"
            write(f, "#SBATCH --account=s2701\n")
            write(f, """#SBATCH --constraint="[mil]"\n""")  # |cas|sky
            if runtime > 12
                write(f, "#SBATCH --qos=long\n")  # for jobs 12-24 hours.
                if runtime > 24
                    runtime = 24  # this is max allowed in Discover
                end
            end
            write(f, "#SBATCH --no-requeue\n")
        elseif hpc == "osc"
            write(f, "#SBATCH --account=PAS1785\n")
        elseif hpc == "unity"
            write(f, """#SBATCH --constraint="[cascade|skylake]"\n""")
        end
        write(f, "#SBATCH --job-name=$(jobname).job\n")
        write(f, "#SBATCH --output=.out/$(jobname).out\n")
        write(f, "#SBATCH --time=$(runtime):00:00\n")
        write(f, "#SBATCH --nodes=1 --ntasks=$(cores)\n")  # 
        write(f, "#SBATCH --exclusive\n")  #  use whole node with all cores without sharing 
        write(f, "#SBATCH --mem=$(memory)\n")  # zero implies use all available memory
        write(f, "#SBATCH --mail-type=ALL\n")
        write(f, "#SBATCH --mail-user=yadav.111@osu.edu\n")
        # Delay slurm job start by a few minutes so that the memory related problem of all nodes reading the same input nc file is resolved  
        write(f, "#SBATCH --begin=now+$(begin_delay)minutes\n")  ## minute, hour
        write(f, "\n")
        write(f, "echo ==============================================================================\n")
        write(f, "echo \$SLURM_SUBMIT_DIR\n")    
        if hpc == "discover"
            write(f, "cd \$LOCAL_TMPDIR\n")
        else
            write(f, "cd \$TMPDIR\n")
        end
        write(f, "date; hostname; pwd\n")
        write(f, "echo Blender run for $(out_subfolder) start_idx = $(start_idx) end_idx = $(end_idx) valid_pix_count = $(valid_pix_count). \n\n")
        # write(f, "export JULIA_NUM_THREADS=\$SLURM_NTASKS\n")
        write(f, "sleep 5\n")
        write(f, "julia --version\n")
        if hpc == "discover"
            write(f, "cp -r /discover/nobackup/projects/coressd/Github/verse .\n")
        else
            write(f, "cp -r ~/Github/verse .\n")
        end
        write(f, "julia verse/Julia/call_Blender_v15.jl $(out_subfolder) $(start_idx) $(end_idx)\n\n")
        write(f, "squeue --job \$SLURM_JOBID \n")
        write(f, "echo List of files on TMPDIR\n")
        write(f, "echo ---------------------------------------------------------------------\n")
        write(f, "ls -ltrh\n")
        write(f, "ls logs*|wc -l\n")  # list log files on node
        # write(f, "tar -czf logs_$(start_idx)_$(end_idx).tar.gz logs*\n")  # list log files on node
        # write(f, "cp logs_$(start_idx)_$(end_idx).tar.gz \$SLURM_SUBMIT_DIR\n")  # list log files on node
        write(f, "echo Finished Slurm job \n\n")
    end
    run(`sbatch $(job_file)`)  # submit the job
end

using Rasters, NCDatasets
DataDir = "$base_folder/Inputs"  # must exist. Old --> NoahMP

# 2. Read the Input netCDF file
A = RasterStack("$(DataDir)/WY_merged/$(water_year)_seup_modis.nc", lazy=true)
A = A[:SWE_tavg][Ti=1];
szY = size(A, 2)  # get size in Y-dimension; here dim = 2
# szY = 10 # only for debug and testing, use smaller number
job_count = 0
delay_multiplier = 2  # delay the consecutive slurm job by ~3 minutes
if occursin("asc.ohio-state.edu", host_machine)
    delay_multiplier = 7  # longer delay on unity becuase of slow speeds in moving data (network related)
end
total_runtime = 0
# for i in StepRange(501, step, 600)  # for testing
for i in StepRange(1, step, szY)
    start_idx = i
    end_idx = start_idx + step - 1
    if end_idx > szY
        end_idx = szY
    end
    # begin_delay = Int((job_count) * step/3)  
    begin_delay = Int(job_count * delay_multiplier)
    jobname = "$(start_idx)_$(end_idx)"
    B = A[1:end, start_idx:end_idx]  # ERROR: LoadError: NetCDF error: NetCDF: Start+count exceeds dimension bound (NetCDF error code: -57)
    """
    # TODO Save of copy of sub-array A to hard-drive somewhere. This will by copied by call_blender_v14 script. To prevent slurm jobs failing on Unity.
    # likely due to multiple jobs copying the same file to nodes 
    A = A[1:end, start_idx:end_idx, :]    
    """
    valid_pix_ind = findall(!ismissing, B)
    valid_pix_count = length(valid_pix_ind)
    # estimate runtime as function of the number of cores used. 210 seconds = 3.5 mins; ie 1 pixel processing time ~ 0.3 min.
    runtime = Int(ceil(valid_pix_count/(cores*180)))  # 150 120; 210; with exclusive, how many cores we get is not certain, but just an estimate
    global total_runtime += runtime
    # 240 seconds : timeout error on OSC, so redude time
    # runtime = "08:00:00"  # 12 "24:00:00"  # default (max for Discover)
    # runtime = "$(approx_time):00:00"
    println(start_idx, " ", end_idx, " ", valid_pix_count, " hours ", runtime)

    # # Use next 5 lines only if re-running few folders due to error/timeout. Still testing. 
    # nc_outDir = "$base_folder/Runs/WY$(water_year)/temp_nc/outputs_$(start_idx)_$(end_idx)"
    # if isdir(nc_outDir)  # process only if the pixel is not already processed
    #     println("Exiting because following folder already processed: $nc_outDir")
    #     continue
    # end
    
    # create_job(hpc_name, jobname, cores, memory, runtime, "V14x_WY$(water_year)", start_idx, end_idx, valid_pix_count)
    create_job(hpc_name, jobname, cores, memory, runtime, "WY$(water_year)", start_idx, end_idx, valid_pix_count, begin_delay)
    sleep(1)
    global job_count += 1
end
println("Total Slurm jobs submitted = $job_count")
println("Total Runtime = $total_runtime")
# =============================================================================================================================
# =============================================================================================================================
