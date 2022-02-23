#!/usr/bin/bash
#
export PYEMMA_NJOBS=1
export OMP_NUM_THREADS=1

#hps=(53 52 60 47 86 81)
#states=(2 2 2 2 2 2)
hps=(86)
states=(2)
lag=100
n_samples=100

for i in "${!states[@]}";
do
  msmsense sample-metastable \
    -h "${hps[$i]}" \
    -n "${states[$i]}" \
    -l ${lag} \
    -a ${n_samples} \
    -i hpsample.h5 \
    -d /home/rob/Data/DESRES \
    -t  /home/rob/Data/DESRES/DESRES-Trajectory_1FME-0-protein/1FME-0-protein/protein.pdb \
    -g '*1FME*/**/*.xtc' \
    -o 1fme \
    -s 49587
done
