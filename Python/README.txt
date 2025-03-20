Files here are listed in the order it should be run. Only the required files are described. Refer to the main readme on how this files should be run. Many of these scripts are called by slurm job.   
These python scripts are used generate the input file for Blender run.  

Extract_LIS.py 
===============
- extract LIS file in format required by Blender script
- save the files inside ../Blender/Inputs/ folder

process_modis_cgf.py
===================  
- process MODIS files
- apply tree-fraction correction  
- resample to LIS resolution
- save SCF.nc inside the same ../Blender/Inputs/ folder 

At the end of this scripts we will have five Input NetCDF files to start the Blender Runs.

