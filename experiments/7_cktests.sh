#!/usr/bin/zsh

export PYEMMA_NJOBS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

#@click.option('-n', '--n-metastable', type=int, help='Number of metastable states')
#@click.option('-l', '--lag', type=int, help='Lag of model')
#@click.option('-a', '--num-repeats', type=int, help='Number of samples')
#@click.option('-n', '--num-cores', type=int, help='Number of cpu cores to use.', default=1)
#@click.option('-i', '--hp-sample', type=Path, help='Path to file that contains the hyperparameter samples')
#@click.option('-d', '--data-dir', type=Path, help='Base directory used to determine trajectory and topology paths')
#@click.option('-t', '--topology-path', type=Path, help='Topology path')
#@click.option('-g', '--trajectory-glob', type=str, help='Trajectory glob string relative to --data-dir')
#@click.option('-o', '--output-dir', type=Path, help='Path to output directory')
#@click.option('-s', '--seed', type=int, help='Random seed', default=None)
#@click.argument('hp-ixs', type=int, nargs=-1)

nohup msmsense score \
 -n 2 \
 -i hpsample.h5 \
 -d /home/rob/Data/DESRES \
 -t /home/rob/Data/DESRES/DESRES-Trajectory_1FME-0-protein/1FME-0-protein/protein.pdb \
 -g '*1FME*/**/*.xtc' \
 -r 100 \
 -n 1 \
 -l 41 \
 -o 1fme/cktest \
 -s 49587 \
 1 2
