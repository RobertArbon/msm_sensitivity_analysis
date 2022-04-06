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

out_dir=/home/rob/Research/msm_sensitivity_analysis/analysis/1fme
if [ ! -f "${out_dir}/summary.h5" ]; then
  sensetools summarise -r 1fme/scores ${out_dir}
fi
sensetools select lag ${out_dir}/summary.h5 ${out_dir} -c 0.01
#sensetools select dominant ${out_dir}/summary.h5 ${out_dir} -c 5
#sensetools select models ${out_dir}/summary.h5 hpsample.h5 -c 0.02
#sensetools plot vamps-ranked ${out_dir}/summary.h5 hpsample.h5 ${out_dir}
#sensetools plot vamps-vs-gap ${out_dir}/summary.h5 hpsample.h5 ${out_dir}
#sensetools plot its ${out_dir}/summary.h5 ${out_dir}
