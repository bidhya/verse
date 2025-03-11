""" Create and submit Julia Blender jobs to Slurm on Discover, OSC HPCs.

    USAGE: pass Water_Year and stepsize (number of rows)
    julia verse/Julia/submit_slurm.jl 2016 1
        good values for stepsize: 1, 2, ... upto 10. last tested with 3. Use 1 for more finer control.  

    Note: 1. This script will call the main script call_Blender_v18.jl
          2. The script is called with 3 arguments: water_year, resolution, stepsize
          3. The script reads the input netcdf file, and creates slurm jobs for each set of rows (stepsize) in the Y-direction.
          4. The script calculates the number of pixels in each set of rows, and based on a threshold, creates a slurm job.
          5. Calculates the runtime for each job based on the number of pixels but currently hard-coded to 12 or 24 (long QOS).
          7. Mechanism of a delay in the start of each job to prevent memory related errors on Discover.
    Changelog:  
    - increase number of cores to 125
    - replacing julia installed with conda with the one installed in ../coressd/installs/julia
    - remove RES from the arguments to simplify the script. RES is not used in the script.
"""
arg_len = length(ARGS)
out_subfolder = ARGS[1]  # WY2016. output subfolder relative to input files; temp_nc saved here
water_year = out_subfolder[end-3:end] #last 4 chars are assumed year, else error
# step is the count of rows (y-direction) of netcdf file to process in one job/run 
step = ARGS[2] #|| 35 on discover
step = parse(Int32, step)
# Grid Resolution: 050 = 0.050 degree (5km); 025 = 0.025 degree; 010 = 0.01 degree; 001 = 0.001 degree (100 meters) etc
# RES = ARGS[3] # 010
memory = "184gb"

# mkpath(logDir)
mkpath("slurm_jobs/$(water_year)/.out")  # /$(RES)
cd("slurm_jobs/$(water_year)")  # /$(RES)

# 1. Setup inputs and output directories and file locations
# select the root directory (this will be different on windows, linux or Mac) 
host_machine = gethostname()
if occursin("discover", host_machine) #|| occursin("borg", host_machine)
    root_dir = "/discover/nobackup"
    base_folder = "$root_dir/projects/coressd/Blender"
    # base_folder = "$root_dir/byadav/coressd/Blender"  # aside: temporary when main storage was full.
    DataDir = "$root_dir/projects/coressd/Blender/Inputs"  # _$(RES)  INDIR. must exist  (Old = NoahMP)
    hpc_name = "discover"
    cores = 125  # 110 for stepsize 125 for 5 km run.
    # if RES == "010"
    #     cores = 125  # 120, 95, 80. # maybe use less cores to prevent NODE_FAIL error for 1km run (TBD). 
    # end
    memory = "0" #"184gb"
elseif occursin(".osc.edu", host_machine)
    root_dir = "/fs/ess/PAS1785/coressd"  # "/fs/scratch/PAS1785/coressd"
    base_folder = "$root_dir/Blender"
    DataDir = "$base_folder/Inputs"  # _$(RES)
    hpc_name = "osc"
    cores = 48  # 40
    memory = "0" #"175gb"
elseif occursin("asc.ohio-state.edu", host_machine)  # .unity
    root_dir = "/fs/project/howat.4/yadav.111/coressd"  # homedir()  #  Unity
    # base_folder = "/home/yadav.111/Github/Blender"  # old
    base_folder = "$root_dir/Blender"  # "$root_dir/Github/coressd/Blender"
    DataDir = "$base_folder/Inputs"  # _$(RES)
    hpc_name = "unity"
    cores = 40  # 39 24 cores with 96GB memory for old node
    memory = "0"  # 186gb seems max allowed
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
        write(f, "#SBATCH --error=.out/$(jobname).err\n")
        write(f, "#SBATCH --time=$(runtime):00:00\n")
        write(f, "#SBATCH --nodes=1 --ntasks=$(cores)\n")  # 
        write(f, "#SBATCH --exclusive\n")  #  use whole node with all cores without sharing 
        write(f, "#SBATCH --mem=$(memory)\n")  # zero implies use all available memory
        write(f, "#SBATCH --mail-type=FAIL\n")  #   ALL; BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_90 (reached 90 percent of time limit)
        write(f, "#SBATCH --mail-user=yadav.111@osu.edu\n")
        # Delay slurm job start by a few minutes so that the memory related problem of all nodes reading the same input nc file is resolved  
        write(f, "#SBATCH --begin=now+$(begin_delay)minutes\n")  ## minute, hour
        write(f, "\n")
        write(f, "echo ==============================================================================\n")
        write(f, "echo \$SLURM_SUBMIT_DIR\n")    
        if hpc == "discover"
            write(f, "cd \$LOCAL_TMPDIR\n")
            # # New for using Julia installed using conda on Discover (June 20, 2024). Will not work any
            # write(f, "module load anaconda\n")
            # write(f, "conda activate julia\n")
        else
            write(f, "cd \$TMPDIR\n")
        end
        write(f, "date; hostname; pwd\n")
        write(f, "echo Blender run for $(out_subfolder) start_idx = $(start_idx) end_idx = $(end_idx) valid_pix_count = $(valid_pix_count). \n\n")
        # write(f, "export JULIA_NUM_THREADS=\$SLURM_NTASKS\n")  # mixing thread with distributed and GC.gc might be creating Segmeentation Fault problem.  
        write(f, "sleep 3\n")
        # write(f, "julia --version\n")
        write(f, "which julia\n")
        if hpc == "discover"
            write(f, "cp -r /discover/nobackup/projects/coressd/Github/verse .\n")
        else
            write(f, "cp -r ~/Github/verse .\n")
        end
        write(f, "julia verse/Julia/call_Blender_v18.jl $(out_subfolder) $(start_idx) $(end_idx)\n\n")  #  $(RES)
        write(f, "squeue --job \$SLURM_JOBID \n")
        write(f, "echo List of files on TMPDIR\n")
        write(f, "echo ---------------------------------------------------------------------\n")
        write(f, "ls -ltrh\n")
        write(f, "echo -n Count of outputs_txt files: ; ls outputs_txt_*|wc -l\n")  # more specific: outputs_txt_$(start_idx)_$(end_idx)
        # write(f, "echo -n Logs: ; ls logs_$(start_idx)_$(end_idx)|wc -l\n")  # logs_* can be  misleading as it matches other files with logs prefix, so be explicit.  
        
        # To move log files from node to main folder, create a tar file and copy to main folder.
        write(f, "ECODE=\$?\n")
        write(f, "if [ \$ECODE -eq 0 ]; then\n")
        write(f, "\ttar -czf logs_$(start_idx)_$(end_idx).tar.gz logs*\n")  # list log files on node
        write(f, "\tmkdir -p $base_folder/Runs/$out_subfolder/logs\n")  # /$(RES) will copy logs here.  
        write(f, "\tcp logs_$(start_idx)_$(end_idx).tar.gz $base_folder/Runs/$out_subfolder/logs\n")  # /$(RES) 
        # write(f, "cp logs_$(start_idx)_$(end_idx).tar.gz \$SLURM_SUBMIT_DIR\n")
        # TODO Write/Call a python script (or Jupyter notebook) to QA/QC the run: count, tables, figures.  
        write(f, "\techo Finished Slurm job successfully\n")
        # Cleanup the files on the node (not strictly necessary but good practice)
        write(f, "\techo Cleaning up files on node\n")
        write(f, "\trm -rf verse\n")
        write(f, "\trm -rf outputs_txt_*\n")
        write(f, "\trm -rf logs*\n")
        write(f, "else\n")
        write(f, "\techo Error in Blender run\n")
        write(f, "fi\n\n")

        # if hpc == "discover"
        #     write(f, "conda deactivate\n")  # New for Julia installed using conda
        # end
        write(f, "echo ==============================================================================\n")
    end
    run(`sbatch $(job_file)`)  # submit the job
end

using Rasters, NCDatasets
# DataDir = "$base_folder/Inputs"  # must exist. Old --> NoahMP

# 2. Read the Input netCDF file
# A = RasterStack("$(DataDir)/WY_merged/$(water_year)_seup_modis.nc", lazy=true) # OLD: 1 file with all variables
# A = Raster("$DataDir/lis/WY$(water_year)/Qg_tavg.nc", lazy=true) # ERROR: LoadError: ArgumentError: invalid index: :Qg_tavg of type Symbol
files = (
    "$DataDir/WY$(water_year)/SCF.nc",
    "$DataDir/WY$(water_year)/Qg_tavg.nc"
    )
A = RasterStack(files; lazy=true)

A = A[:Qg_tavg][Ti=1];
szY = size(A, 2)  # get size in Y-dimension; here dim = 2. 6500 for 1km run.
# szY = 10 # only for debug and testing, use smaller number
job_count = 0
delay_multiplier = 0.5 # delay the consecutive slurm job by ~1 minutes
if occursin("asc.ohio-state.edu", host_machine)
    delay_multiplier = 5  # longer delay on unity becuase of slow speeds in moving data (network related)
end
total_runtime = 0
cum_valid_pix_count = 0
row_count = 0
start_idx = 1  # Alert: if "i" changed in for loop below, update this value to the first value of i, else job will start at 1 causing node_fail error.
# start_idx should normally be 1, but for testing, it can be set to a higher value (matching the starting value in for loop below)
# test_start_idx = 1000
# start_idx = test_start_idx
# for i in StepRange(test_start_idx, step, 1003)  # for testing
# for i in StepRange(6001, step, szY)
for i in StepRange(1, step, szY)
    # start_idx = i
    end_idx = i + step #- 1
    if end_idx > szY
        end_idx = szY
    end
    B = A[1:end, i:end_idx]  # ERROR: LoadError: NetCDF error: NetCDF: Start+count exceeds dimension bound (NetCDF error code: -57)
    """
    # TODO Save of copy of sub-array A to hard-drive somewhere. This will by copied by call_blender_v14 script. To prevent slurm jobs failing on Unity.
    # likely due to multiple jobs copying the same file to nodes 
    A = A[1:end, start_idx:end_idx, :]    
    """
    valid_pix_ind = findall(!ismissing, B)
    valid_pix_count = length(valid_pix_ind)
    row_count += step
    cum_valid_pix_count += valid_pix_count
    # println("Processing rows: ", i, " to ", i+step, "and cum_valid_pix_count = ", cum_valid_pix_count)
    pix_count_threshold = 80000  # 80000 for 5km run; 20000 for 1km run
    # Update the pix_count_threshold based on row index (ie, latitude) to account for the processing time which is higher at lower latitude. This needs further verification.
    # dimensions(sizes): x(11700), y(4700), time(366)
    if start_idx < 100  # needs most memory. old=120.
        pix_count_threshold = 60000  # 59000
    elseif start_idx < 2000
        pix_count_threshold = 80000  # 79000
    elseif start_idx < 3000
        pix_count_threshold = 90000  # 89000
    else
        pix_count_threshold = 110000  # 109000
    end

    if cum_valid_pix_count > pix_count_threshold || end_idx == szY  # last job may not have enough pixels to pass the pix count threshold
        # prepare slurm job
        # begin_delay = Int((job_count) * step/3)  
        begin_delay = Int(round(job_count * delay_multiplier))
        jobname = "$(start_idx)_$(end_idx)"
        # Use next 5 lines only if re-running few folders due to error/timeout. Still testing. 
        nc_outDir = "$base_folder/Runs/WY$(water_year)/temp_nc/outputs_$(start_idx)_$(end_idx)"
        if isdir(nc_outDir)  # process only if the pixel is not already processed
            println("Job Not created because output_folder $nc_outDir already exist")
            continue
        end
        # estimate runtime as function of the number of cores used. 270 seconds = 4.5 mins; ie 1 pixel processing time ~ 0.3 min.
        runtime = Int(ceil(valid_pix_count/(cores*420)))  # hours. 270 for v_15=180; with exclusive, how many cores we get is not certain, but just an estimate
        runtime = max(runtime, 12)  # bump smaller runtime to 12 hours (max allowed by default on Discover). Due to timeout errors for short runs.   
        # if runtime < 10 && hpc_name=="discover"
        #     # keep getting runtime errors when less number of pixels are processed.
        #     # say in less than 4 hours, likely all cores are not being utilized efficiently
        #     runtime = 11
        # end
        global total_runtime += runtime
        # 240 seconds : timeout error on OSC, so redude time
        # runtime = "08:00:00"  # 12 "24:00:00"  # default (max for Discover)
        # runtime = "$(approx_time):00:00"
        println(start_idx, " ", end_idx, " rows:", row_count, " pixels:", cum_valid_pix_count)  # , " hours ", runtime
        
        # create_job(hpc_name, jobname, cores, memory, runtime, "V14x_WY$(water_year)", start_idx, end_idx, valid_pix_count)
        create_job(hpc_name, jobname, cores, memory, runtime, "WY$(water_year)", start_idx, end_idx, cum_valid_pix_count, begin_delay)
        # sleep(1)
        global job_count += 1
        # reset the counter
        global row_count = 0
        global cum_valid_pix_count = 0
        global start_idx = end_idx + 1
        if end_idx >= szY
            println("Exiting loop. End index already included in last job.")
            break
        end
    # else
    #     println("No slurm job created. Cumulative pixel count = ", cum_valid_pix_count, " less than threshold = ", pix_count_threshold)   
    end
end
println("Total Slurm jobs submitted = $job_count")
println("Total Runtime = $total_runtime")
# =============================================================================================================================
# =============================================================================================================================
