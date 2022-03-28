#!/usr/bin/bash
#
export PYEMMA_NJOBS=1
export OMP_NUM_THREADS=1

lag=41
n_cores=6
hps=(53 52 60 47 86 81)
for i in "${!hps[@]}";
do
  echo "${hps[$i]}"
  msmsense project-evs \
    -l ${lag} \
    -i hpsample.h5 \
    -d /home/rob/Data/DESRES \
    -t  /home/rob/Data/DESRES/DESRES-Trajectory_1FME-0-protein/1FME-0-protein/protein.pdb \
    -g '*1FME*/**/*.xtc' \
    -n ${n_cores} \
    -o 1fme/project_evs \
    -s 49587 \
    -h "${hps[$i]}"  \
    -p 1fme/evs \
    52 53 47 60 81 86
done

