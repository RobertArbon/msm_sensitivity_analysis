#!/usr/bin/bash
#
export PYEMMA_NJOBS=1
export OMP_NUM_THREADS=1

#hps=(53 52 60 47 86 81)
#@click.option('-l', '--lag', type=int, help='Lag of model')
#@click.option('-k', '--processes', type=int, help='Number of process to compare')
#@click.option('-q', '--num-cuts', type=int, help='Number of quantiles to cut the EV into')
#@click.option('-i', '--hp-sample', type=Path, help='Path to file that contains the hyperparameter samples')
#@click.option('-d', '--data-dir', type=Path, help='Base directory used to determine trajectory and topology paths')
#@click.option('-t', '--topology-path', type=Path, help='Topology path')
#@click.option('-g', '--trajectory-glob', type=str, help='Trajectory glob string relative to --data-dir')
#@click.option('-r', '--num-repeats', type=int, help='Number of bootstrap samples')
#@click.option('-n', '--num-cores', type=int, help='Number of cpu cores to use.', default=1)
#@click.option('-o', '--output-dir', type=Path, help='Path to output directory')
#@click.option('-s', '--seed', type=int, help='Random seed', default=None)
#@click.argument('hp-ixs', type=int, nargs=-1)
states=2
n_cuts=50
lag=41
n_samples=2
n_cores=2

msmsense sample-evs \
  -l ${lag} \
  -k ${states} \
  -q ${n_cuts} \
  -i hpsample.h5 \
  -d /home/rob/Data/DESRES \
  -t  /home/rob/Data/DESRES/DESRES-Trajectory_1FME-0-protein/1FME-0-protein/protein.pdb \
  -g '*1FME*/**/*.xtc' \
  -r ${n_samples} \
  -n ${n_cores} \
  -o 1fme/evs \
  -s 49587 \
  53

#  52 60 47 86 81
