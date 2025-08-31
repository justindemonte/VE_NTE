This repository contains scripts for the simulation study and data analysis in the manuscript "Assessing vaccine effectiveness in observational studies via nested trial emulation" by [DeMonte et al.](https://arxiv.org/abs/2403.18115)

Simulations were conducted in R and Python on on a SLURM-managed High Performance Computing (HPC) cluster, and the shell scripts used are specific to SLURM.  

The following steps can be used to reproduce the simulation results.  All files in the gitHub repository should be placed under the same user directory.  For example, if the user directory is titled 'userDirectory', then the top level of that directory must contain the following:

userDirectory/sims

userDirectory/virtualEnv

Users can create a python virtual environment by running the following commands:

```
module add anaconda/2023.03
python3 -m venv [user directory]/virtualEnv/bin/activate
module add r/4.2.2
source [user directory]/virtualEnv/bin/activate

python3.10 -m pip install --upgrade pip
python3.10 -m pip install --upgrade 'delicatessen'
python3.10 -m pip install --upgrade 'pandas'
```

Each .sh shell script performs a set of simulations under certain parameters.  The following naming convention is used for the shell scripts.  The first letter is 'm' for the "main" simulations (described in Section 3 of the main text and Section S2 of the supplementary materials) or 'e' for "extension" simulations (described in Section 4.2 of the main text and Section S3 of the supplementary materials); the second letter is 'n' for the "naive" or comparator analysis or 'p' for the "proposed" analysis; the third letter is 's' for "splines" or 'p' for "polynomials", specifying the functional form of time used in the analysis models; the number refers to the scenario.  For example, running the shell file sim_enp1.sh will reproduce results in the "Unstandardized" columns and "Scenario 1" rows of Table S3 in the supplemental materials (extension methods, naive analysis, using polynomial functions of time, scenario 1).  Scripts truM.sh and truE.sh are used to approximate the true parameter values for scenarios 1-3 of the "main" and "extension" simulations, respectively.  

The shell scripts for reproducing the simulations must be run from the working directory userDirectory/sims.  Note that 'userDirectory/sims' contains several empty subdirectories.  Running the shell scripts will populate these subdirectories with results files.  Once all .sh files have run successfully, the R program 204gather.R can be run (from the working directory userDirectory/sims) to produce LaTeX code which in turn can be run to reproduce results tables S1-S3 in the supplemental materials.  The directory 'virtualEnv' contains the python virtual environment containing the python version and all dependencies under which the results in the manuscript were generated.  The required R version and dependencies are as follows:

R version 4.4.0

tidyverse v2.0

matrixStats v1.5.0

reticulate v1.43.0

VGAM v1.1-13

Hmisc v5.2-3

xtable v1.8-4

The directory 'VEapplication' contains scripts used to produce the applied data analysis results in the manuscript.  These scripts depend on the Abruzzo COVID-19 dataset, which is the property of Dr. Lamberto Manzoli and colleagues and cannot be shared publicly.  These data were used under Dr. Manzoli's written permission.  

For any questions, please email demjus01@live.unc.edu.
