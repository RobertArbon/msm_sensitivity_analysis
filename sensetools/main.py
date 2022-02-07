from pathlib import Path

import click

from .summary import main as summary


@click.group()
def cli():
    pass


@cli.group()
@click.argument('directory', help='Directory to trials',
                type=click.Path(exists=False, file_okay=False, resolve_path=True, path_type=Path)
                )
def summarise(directory):
    summary(directory)



