#!/usr/bin/bash
#
export PYEMMA_NJOBS=1
export OMP_NUM_THREADS=1

#@cli.command()
#@click.option('-l', '--lag', type=int, help='Lag of model')
#@click.option('-i', '--hp-sample', type=Path, help='Path to file that contains the hyperparameter samples')
#@click.option('-d', '--data-dir', type=Path, help='Base directory used to determine trajectory and topology paths')
#@click.option('-t', '--topology-path', type=Path, help='Topology path')
#@click.option('-g', '--trajectory-glob', type=str, help='Trajectory glob string relative to --data-dir')
#@click.option('-n', '--num-cores', type=int, help='Number of cpu cores to use.', default=1)
#@click.option('-o', '--output-dir', type=Path, help='Path to output directory')
#@click.option('-s', '--seed', type=int, help='Random seed', default=None)
#@click.option('-h', '--hp-ix', type=int, help='HP index to model')
#@click.option('-p', '--project-dir', type=Path, help='Directory containing the ev projections')
#@click.argument('project-ixs', type=int, nargs=-1)

lag=41
n_cores=2

msmsense project-evs \
  -l ${lag} \
  -i hpsample.h5 \
  -d /home/rob/Data/DESRES \
  -t  /home/rob/Data/DESRES/DESRES-Trajectory_1FME-0-protein/1FME-0-protein/protein.pdb \
  -g '*1FME*/**/*.xtc' \
  -n ${n_cores} \
  -o 1fme/project_evs \
  -s 49587 \
  -h 53 \
  -p 1fme/evs \
  52 53 47 60 81 86


#  52 60 47 86 81
