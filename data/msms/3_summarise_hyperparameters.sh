#!/usr/bin/zsh

#  -h, --help            show this help message and exit
#  -i HP_SAMPLE, --hp-sample HP_SAMPLE
#                        Path to file that contains the hyperparameter samples
#  -d DATA_DIR, --data-dir DATA_DIR
#                        Base directory used to determine trajectory and topology paths
#  -t TOPOLOGY_PATH, --topology-path TOPOLOGY_PATH
#                        Topology path
#  -g TRAJECTORY_GLOB, --trajectory-glob TRAJECTORY_GLOB
#                        Trajectory glob string relative to --data-dir
#  -r NUM_REPEATS, --num-repeats NUM_REPEATS
#                        Number of bootstrap samples
#  -n NUM_CORES, --num-cores NUM_CORES
#                        Number of cpu cores to use.
#  -l LAGS, --lags LAGS  Lags as a Python range specification start:end:stride
#  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
#                        Path to output directory
#  -s SEED, --seed SEED  Random seed
#
export PYEMMA_NJOBS=1
export OMP_NUM_THREADS=1

#time msmsense score \
# -i hpsample.h5 \
# -d /home/rob/Data/DESRES \
# -t  /home/rob/Data/DESRES/DESRES-Trajectory_1FME-0-protein/1FME-0-protein/protein.pdb \
# -g '**/*.xtc' \
# -r 100 \
# -n 4 \
# -l 1:102:10 \
# -o 1fme \
# -s 49587 \
# {85..99}

if [ ! -f "1fme/summary.h5" ]; then
  sensetools summarise  1fme
fi

sensetools select lag 1fme/summary.h5 1fme -c 0.01
sensetools select dominant 1fme/summary.h5 1fme -c 5
sensetools plot vamps-ranked 1fme/summary.h5 hpsample.h5 1fme