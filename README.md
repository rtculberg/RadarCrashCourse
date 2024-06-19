# RadarCrashCourse
Practical exercise notebooks for the CHIRP radar crash course.

## Installation Instructions

This notebook is designed to be run on G2 through tunneled session. 
[1] Log into G2 with tunneled session:`ssh netID@g2-login.coecis.cornell.edu -L PORT#:culberg-cpu-01:PORT#`. Replace `PORT#` with a value between 8000-10000.  
[2] In your home directory on G2, clone the training repository: `git clone git@github.coecis.cornell.edu:CHIRP-Lab/RadarCrashCourse.git`
[3a] If you already have a conda environment setup for using the RESAnalysisTools notebooks, activate that environment. Then execut `pip install mat73` to install an additional package needed for this notebook. 
[3b] If you do not have a conda environment for radar processing, install the complete environment from the yml file. Change directory into your new `RadarCrashCourse` directory. Then execute: `conda env create -f environment.yml`.
[4] Start an interactive session on our compute node: `srun -p culberg-interactive -n 1 --mem=32G --pty --nodelist=culberg-cpu-01 /bin/bash`
[5] Activate your conda environment: `conda activate notebooks` (or whatever you named your environment if you didn't use the yml file for install).
[6] Start a Jupyter notebook server: `XDG_RUNTIME_DIR=/tmp/netID jupyter-notebook --ip=0.0.0.0 --port=PORT#`. Don't forget to replace `netID` with your netID and `PORT#` with the port number you used when logging in to G2.
[7] Open the Jupyter server in your notebook using the URL starting with `https://127.0.0.0....`. Now you are ready to code!
