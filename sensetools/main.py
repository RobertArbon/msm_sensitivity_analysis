from pathlib import Path

import click

from .summary import main as summary
from .selections import select_dominant, select_lag


@click.group()
def cli():
    pass


@cli.command()
@click.argument('directory',
                type=click.Path(exists=True, file_okay=False, resolve_path=True, path_type=Path))
def summarise(directory):
    summary(directory)

@cli.group()
def select():
    pass


@select.command()
@click.argument('summary',
                type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path))
@click.argument('output_directory',
                type = click.Path(exists=False, file_okay=False, resolve_path=True, path_type=Path))
@click.option('-c', '--cutoff', required=True, help='Cutoff gradient for convergence', default=0.01, type=float,
              show_default=True )
def lag(summary, output_directory, cutoff):
    output_directory.mkdir(exist_ok=True, parents=True)
    select_lag(summary, output_directory, cutoff)


@select.command()
@click.argument('summary',
                type=click.Path(exists=True, dir_okay=False, resolve_path=True, path_type=Path))
@click.argument('output_directory',
                type=click.Path(exists=False, file_okay=False, resolve_path=True, path_type=Path))
@click.option('-c', '--cutoff', required=True, help='Cutoff for number of timescales', default=1, type=int,
              show_default=True)
def dominant(summary, output_directory, cutoff):
    output_directory.mkdir(exist_ok=True, parents=True)
    select_dominant(summary, output_directory, cutoff)
