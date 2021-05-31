#! /usr/bin/env python
# -*- coding: utf-8 -*-
import click
from cmaqwf.settings import setting as s

proj = s.get_active_project()

@click.group("run", short_help="Run CMAQ Workflow components")
def run_cmd() -> None:
    pass


@run_cmd.command("mcip", short_help="mcip component")
@click.option(
    "-dom",
    "--domain",
    help="Set the domain id(s)",
    multiple=True,
    type=int,
    default=None,
    show_default=True,
)
@click.option(
    "-y",
    "--year",
    help="Set the year(s) to process",
    multiple=True,
    type=int,
    default=None,
    show_default=True,
)
@click.option(
    "-m",
    "--month",
    help="Set the month(s) to process",
    multiple=True,
    type=int,
    default=None,
    show_default=True,
)
@click.option(
    "-d",
    "--day",
    help="Set the day(s) to process",
    multiple=True,
    type=int,
    default=None,
    show_default=True,
)
def mcip_cmd(domain, year, month, day) -> None:
    domain = list(domain)
    year = list(year)
    month = list(month)
    day = list(day)
    click.echo('mcip_run_script')

