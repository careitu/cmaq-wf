#! /usr/bin/env python
# -*- coding: utf-8 -*-
import click
from run import run_cmd

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(
    help="\n  Manage CMAQ Workflow (v0.0.1)\n",
    epilog="Try 'cmaq run mcip', 'cmaq settings --show', or 'cmaq run bcon'",
    context_settings=CONTEXT_SETTINGS)
def cli() -> None:
    pass


@cli.command("version", short_help="Show cmaq-wf version")
def version_cmd() -> None:
    print('v0.0.1')


@click.command()
def settings() -> None:
    click.echo('CMAQ global settings')


cli.add_command(version_cmd)
cli.add_command(run_cmd)
cli.add_command(settings)

if __name__ == '__main__':
    cli()
