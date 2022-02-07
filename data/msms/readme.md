This directory estimates the MSM count matrices at a series of lag times.  

- `convert_traj.sh` subsamples the trajectories and converts to the relatively light-weight `xtc` format.
- `hyperparameters.sh` runs the command to generate a hdf5 file (`hyperparameters.h5`) full of hyperparameters. 
This should be the same for each protein. Uses the file `searchspace.yaml`.
- `msm.sh`  runs the command to estimate the count matrices.
