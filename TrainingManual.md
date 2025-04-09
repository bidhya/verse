Training and Handover of Blender Plan for Mike and Jack  

Week1  
======  
Julia Setup: MD and JD  
1-pixel test run: BY  

Julia Setup: Use what I have already setup. 
Just update .bashrc file to use the same Julia setup that I have    
Copy/past the following  
```
cd ~  
vim .bashrc  
Press i to enter Insert Mode and paste following 2 lines (towards the end)  
export PATH=$PATH:/discover/nobackup/projects/coressd/installs/julia-1.11.4/bin  
export JULIA_DEPOT_PATH=/discover/nobackup/projects/coressd/installs/.julia  
Hit Escape  
type :wq # write and quit  
Logout and Login. Good practice to logout about any changes to .bashrc  
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
=========
Watershed run
Learn the workflow
Check slrum output and logs files
Check results
Week3

Continental Run
Combine NetCDF slices
Optional

Python Env setup: MD and JD
Create Inputs for Blender
QA/QC of Blender runs
Analysis
Plotting
Challenges/Unknowns/Pitfalls

Many changes and updates in last 2 months
Availability of Nodes