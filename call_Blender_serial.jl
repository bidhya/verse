"""
This is a master control file. All the inputs, outputs, pre-processing, post-processing will be controlled from here
User do not have to change anything inside the "Estimate_v53.jl" script unless updates to main code is required.

USAGE: julia call_Blender_serial.jl

The script can be run directly from a code editor, called from terminal, or through Slurm job scheduler.

B. Yadav: July 17, 2022

Inputs 
====== 
A netCDF file with all the variables
    NB: We could easily have separate netCDF file for each of the inputs

Outputs 
======= 
tmp_txtDir = intermediate outputs incluing log files
nc_outDir = Final outputs by convering the text files to netcdf
    G.nc  Gmelt.nc  Gmeltpv.nc  Gpv.nc  Precip.nc  SWE.nc  SWEpv.nc  Upv.nc  Us.nc

Steps 
=====
Step1: Set input and output directories
Step2: Using Step1 path setups read input netcdf file --> 3D Array [X, Y, and time]
Step3: Call the Blender function
    blender(exp_dir, WRFSWE, WRFP, WRFG, MSCF, AirT) :: This is the only line required for this whole script
        rest of the sections are just organzing inputs and outputs.
Step4: Convert the text files from exp_dir folder to netcdf file
"""

using Rasters
include("Estimate_v53.jl")  #https://docs.julialang.org/en/v1/manual/code-loading/; evaluated in global scope

# 1. Setup inputs and output directories and file locations
# select the root directory (this will be different on windows, linux or Mac) 
host_machine = gethostname()
if occursin("STAFF-BY-M", host_machine)
        root_dir = "C:"  #for windows machine
        base_folder = "$root_dir/Github/Blender"
elseif occursin("hpc.osc.edu", host_machine)
    root_dir = "/fs/ess/PAS1785/coressd" #homedir()  # OSC
    base_folder = "$root_dir/Blender"

elseif occursin("asc.ohio-state.edu", host_machine)
    root_dir = homedir()  #  Unity
    base_folder = "$root_dir/Github/Blender"
else
    println("Unknown computer, manually add root directory before proceeding")
end
# root_dir = ""  # supply custom root directory

DataDir= "$base_folder/nc_files"
# TODO: Give users to pass their own output directory through ARGS
tmp_txtDir = "$DataDir/outputs_txt"    # To save text outputs for each pixel
nc_outDir = "$DataDir/outputs"         # To convert text outputs to netcdf file

# 2. Read the Input netCDF file
A = RasterStack("$base_folder/nc_files/inputs/merged_proj.nc");
A = A[X(38:40), Y=(38:40)]  # For prototyping only, select a smaller subset of pixels

# Get index to iterate over
# need a better/explicit way of extracting x, y, and time indices
tind, yind, xind  = size(A)  
println(tind, yind, xind)

start_time = time_ns()
# 3. Iterate over each pixel using the for loop
# no-data values are called missing in Julia
# In each iteration, if a a missing value is found, we just skip the processing of the pixel
for i in 1:xind
    for j in 1:yind
        val = A[X=i, Y=j][:MODSCAG]  # pick any varaible: if one is missing, all will be missing
        if !ismissing(val[1])
            A_pt = A[X=i, Y=j]  # Get all variables for all time for one pixel
            # extract each of the input variables separately (trying to match of processing was done in prior version with text inputs)
            # but this approach is also useful if we decide to save each input netcdf file separately
            WRFSWE = A_pt["WRFSWE"].data  # Here "data" is an AbstractArray.
            WRFP = A_pt["WRFP"].data
            WRFG = A_pt["WRFG"].data
            MSCF = A_pt["MODSCAG"].data
            AirT = A_pt["WRFT"].data;
            exp_dir = string("$tmp_txtDir/", "Pix_", i, "_", j)  # full path to folder for saving (temporary) text files
            # process only if the pixel is not already processed 
            if !isdir(exp_dir)
                # Call blender for the pixel. This is the only required line
                # rest of the codes are for houskeeking, preprocssing, post-processing
                blender(exp_dir, WRFSWE, WRFP, WRFG, MSCF, AirT)
                # here exp_dir = export directory for holding outputs text and log files (same as older version of code)
            end
        end
    end
end
end_time = time_ns()
running_time = (end_time - start_time)/1e9/60
println("Running Time for Blender (minutes) = $running_time")

# Step 4. PostProcessing: Combine text files into a grid and save as netcdf
pixels = readdir(tmp_txtDir)
pixels = [pix for pix in pixels if startswith(pix, "Pix")];
function text2nc(var)
    """Convert a text file to netcdf
    var: the name of variable (example SWE, Gmelt, G etc.)
    TODO: move this fuction inside the Estimate_v53.jl script
    """
    outRaster = copy(A[:WRFSWE])
    # Initial resulting array
    outRaster[:,:,:] .= missing
    # update the name of variable
    outRaster = rebuild(outRaster; name=var)  #:SWEhat
    for pix in pixels
        pix_xy = split(pix, "_")
        x = parse(Int16, pix_xy[2])
        y = parse(Int16, pix_xy[3])
        outRaster[X=x, Y=y] = readdlm("$tmp_txtDir/$pix/$var.txt");
    end
    # Save to nc file; 
    mkpath(nc_outDir)
    write("$nc_outDir/$var.nc", outRaster)
end

# Call the function for creating netcdf for each of the text output files separately
# Use the name of text file that was saved; Change the filename in main script if so required
text2nc("SWE")
text2nc("Gmelt")
text2nc("G")
text2nc("Precip")
text2nc("Us")
text2nc("Gpv")
text2nc("Gmeltpv")
text2nc("Upv")
text2nc("SWEpv");
running_time = (end_time - start_time)/1e9/60
println("Total Running Time (minutes) = $running_time")
