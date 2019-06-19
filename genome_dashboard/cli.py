# -*- coding: utf-8 -*-

"""Console script for genome_dashboard."""

import sys
import click


@click.command()
def main(args=None):
    """Console script for genome_dashboard."""
    click.echo("Replace this message by putting your code into "
               "genome_dashboard.cli.main")
    click.echo("See click documentation at http://click.pocoo.org/")
    return 0


# if __name__ == "__main__":
#     sys.exit(main())  # pragma: no cover


@click.command()
@click.option('--count', default=1, help='Number of greetings.')
@click.option('--name', prompt='Your name',
              help='The person to greet.')
def hello(count, name):
    """Simple program that greets NAME for a total of COUNT times."""
    for x in range(count):
        click.echo('Hello %s!' % name)


if __name__ == '__main__':
    hello()
