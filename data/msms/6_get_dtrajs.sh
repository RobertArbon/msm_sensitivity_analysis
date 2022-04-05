#!/usr/bin/bash
#
export PYEMMA_NJOBS=1
export OMP_NUM_THREADS=1

msmsense dump-dtrajs \
  -i hpsample.h5 \
  -d /home/rob/Data/DESRES \
  -t  /home/rob/Data/DESRES/DESRES-Trajectory_1FME-0-protein/1FME-0-protein/protein.pdb \
  -g '*1FME*/**/*.xtc' \
  -o 1fme/dtrajs \
  -s 49587 \
  120 24
