import logging
import argparse
import sys
from .linear_to_graph_mapper import LinearToGraphMapper


def run_map_linear_to_graph(args):
    chromosomes = args.chromosomes.split(",")
    LinearToGraphMapper(args.fasta, args.reference, args.data_dir, chromosomes)


def main():
    run_argument_parser(sys.argv[1:])


def run_argument_parser(args):
    parser = argparse.ArgumentParser(
        description='Rough Graph Mapper',
        prog='rough_graph_mapper',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    subparsers = parser.add_subparsers(help='Subcommands')
    subparser = subparsers.add_parser(
        "map_linear_to_graph",
        help="...")

    subparser.add_argument("-r", "--reference", help="Linear reference fasta")
    subparser.add_argument("-f", "--fasta", help="Input fasta file")
    subparser.add_argument("-d", "--data-dir", help="")
    subparser.add_argument("-c", "--chromosomes", help="Comma-separated list of chromosomes")
    subparser.set_defaults(func=run_map_linear_to_graph)
    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)

    if hasattr(args, 'func'):
        args.func(args)

