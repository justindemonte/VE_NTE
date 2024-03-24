This repository contains scripts for the simulation study and data analysis in the manuscript "Assessing COVID-19 vaccine effectiveness in observational studies via nested trial emulation" by DeMonte et al.

Simulations were conducted in R and Python on the longleaf High Performance Computing Cluster at UNC-Chapel Hill.  The shell scripts in this repository are used to specify global parameters and submit jobs to the Slurm Workload Manager.  To replicate the simulations, all files under 'sims', including shell scripts, R and Python scripts, should be placed in the same directory.  The following instructions assume that this folder is in the location ```[user directory]/sims/```.

Before beginning, users need to create a python virtual environment, ensure that the needed R packages are installed, and create a directory structure according to the following 4 steps. 

Step 1: create python virtual environment by running the following commands:

```
module add anaconda/2023.03
python3 -m venv [user directory]/virtualEnv/bin/activate
module add r/4.2.2
source [user directory]/virtualEnv/bin/activate

python3.10 -m pip install --upgrade pip
python3.10 -m pip install --upgrade 'delicatessen'
python3.10 -m pip install --upgrade 'pandas'
```

Step 2: ensure necessary R packages are installed.
This can be done via an interactive R prompt using ```install.packages()```.  The following packages are required:
```simEstVE_xxxx.R``` requires the packages ```rms, dplyr, tibble, tidyr, purrr,``` and ```reticulate```.
simGather.R requires the package ```tidyverse```.

Step 3: edit shell scripts (where applicable).
If the following line appears near the top of the shell script, change the directory to match the one specified in Step 1:
I.e., change
```source /nas/longleaf/home/demjus01/virtualEnv/bin/activate```
to 
```source [user directory]/virtualEnv/bin/activate```

Step 4: create directory structure for simulation results.
The R scripts create output for each replication of the simulation.  The output is stored in a series of subdirectories below the directory which contains all the scripts.  Please create the following (empty) directory structure:

```
[user directory]/sims/results/spline/scenario1/results
[user directory]/sims/results/spline/scenario2/results
[user directory]/sims/results/spline/scenario3/results

[user directory]/sims/results/poly/scenario1/proposed/results
[user directory]/sims/results/poly/scenario1/naive/results
[user directory]/sims/results/poly/scenario2/proposed/results
[user directory]/sims/results/poly/scenario2/naive/results
[user directory]/sims/results/poly/scenario3/proposed/results
[user directory]/sims/results/poly/scenario3/naive/results
```

There are three R scripts and one Python script for running the simulations. 

```00_balancing_numerical_solver_source_code.R``` was adapted from https://github.com/serobertson/BalancingInterceptSolver.  This script is used to calculate the "balancing intercept" for the true logit model for the outcome hazard.

```simGenerate.R``` generates a simulated dataset according to the data generating process described in the manuscript. 

```simEstVE_xxxx.R``` calculates point estimates for the simulated dataset created by simGenerate.R and generates output files.  There are two versions of this script: one for the simulations appearing in the main text of the manuscript ('poly') and one for the simulations appearing in the supplemental materials ('spline').

```simEstVar_xxxx.py``` calculates variance estimates.  Similarly to simEstVE_xxxx.R, there are two versions: 'poly' and 'spline'. 

To replicate the simulations, the user calls the script 'simRun.sh'.  Once all simulations are completed, running the script 'simGather.R' will generate LaTeX format tables of the results which appear in the manuscript (main text and supplement).
	
The script ```LexisDiagram/LexisDiagram.R``` generates Figure 1 in the manuscript.  This small script can be run on any R setup assuming that the necessary packages (```tidyverse``` and ```patchwork```) are installed.  Scripts used to analyze the Abruzzo data are contained in ```VEapplication```.  The analysis was carried out using R and Python on the longleaf High Performance Computing Cluster at UNC-Chapel Hill.  These scripts depend on the Abruzzo COVID-19 dataset, which is the property of Dr. Lamberto Manzoli and colleagues and was used under Dr. Manzoli's written permission.  

Users are encouraged to submit any issues.  For other questions, please email demjus01@live.unc.edu.

