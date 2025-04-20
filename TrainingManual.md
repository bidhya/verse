Training and Handover of Blender Plan for Mike and Jack  

Week1  
======  
Julia Setup: MD and JD  
1-pixel test run: BY  

## Julia Setup: 
Use installation documentation: https://julialang.org/downloads/  
Use the latest available version (Julia 1.11.4 as of this writing)  
Steps  
```
cd /discover/nobackup/your-user-name
wget https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.5-linux-x86_64.tar.gz
tar zxvf julia-1.11.5-linux-x86_64.tar.gz
cd ~  
vim .bashrc  
Press i to enter Insert Mode and paste following 2 lines (towards the end)  
export PATH=$PATH:/discover/nobackup/your-user-name/julia-1.11.4/bin  
export JULIA_DEPOT_PATH=/discover/nobackup/your-user-name/.julia  
Hit Escape  
type :wq # write and quit  
Logout and Login. Good practice to logout about any changes to .bashrc  
```

### Adding required Julia packages  
```
julia
]
add JuMP Ipopt Rasters NCDatasets CSV LoggingExtras Distributions  
Press Backspace key to exit out  
```

Blender test run for pixel: copy paste the following    
```
cd /discover/nobackup/projects/coressd  
cd Github/wshed # we have slurm scripts for test runs (pixel.job and wshed.job)  
sbatch pixel.job 2016 1  
sbatch pixel.job 2016 2  
.  
.  
.  
```


Week2
======  
Watershed run  
Blender test run for a watershed: copy paste the following    
```
cd /discover/nobackup/projects/coressd  
cd Github/wshed # we have slurm scripts for test runs (pixel.job and wshed.job)  
sbatch wshed.job 2016 1  
sbatch wshed.job 2016 2  
.  
.  
.  
```

Week3
======  
Review and fix any outstanding issues from Week1 and 2  
Learn the workflow  
Check slrum output and logs files  
Check results  

Week4
======  
Continental Run  
Combine NetCDF slices  

Optional  
Python Env setup: MD and JD  
Create Inputs for Blender  
QA/QC of Blender runs  
Analysis  
Plotting  
Challenges/Unknowns/Pitfalls  
