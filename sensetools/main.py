from pathlib import Path

import click

from .summary import main as summary


@click.group()
def cli():
    pass


@cli.command()
@click.argument('directory',
                type=click.Path(exists=False, file_okay=False, resolve_path=True, path_type=Path))
def summarise(directory):
    summary(directory)



