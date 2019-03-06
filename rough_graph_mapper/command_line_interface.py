import logging
import shutil
import argparse
import sys
from .linear_to_graph_mapper import LinearToGraphMapper
from .filter_graphalignments import filter_graphalignments


def remove_reads_from_fasta(args):
    with open(args.alignments) as f:
        remove_ids = set((line.split()[0] for line in f))
    logging.info("%d ids found to remove" % len(remove_ids))

    with open(args.fasta) as f:
        skip = False
        for line in f:
            if line.startswith(">"):
                id = line.strip().replace(">", "")
                if id in remove_ids:
                    skip = True
                    continue
                skip = False

            if not skip:
                print(line.strip())


def run_filter(args):
    alignments = args.alignments
    min_mapq = args.min_mapq
    filter_graphalignments(alignments, min_mapq)

def run_map_linear_to_graph(args):
    chromosomes = args.chromosomes.split(",")
    LinearToGraphMapper(args.fasta, args.reference, args.data_dir, chromosomes, n_threads=args.n_threads,
                        write_final_alignments_to_file=args.output_file, skip_mapq_adjustment=args.skip_mapq_adjustment,
                        minimum_mapq_to_graphalign=args.min_mapq)

def main():
    run_argument_parser(sys.argv[1:])


def run_argument_parser(args):

    if shutil.which("bwa") is None:
        logging.error("BWA MEM cannot be found in path. Make sure BWA is installed.")
        sys.exit()
    if shutil.which("minimap2") is None:
        logging.error("minimap2 cannot be found in path. Make sure minimap2 is installed")
        sys.exit()

    parser = argparse.ArgumentParser(
        description='Rough Graph Mapper',
        prog='rough_graph_mapper',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    subparsers = parser.add_subparsers(help='Subcommands')
    subparser_linear = subparsers.add_parser("map_linear_to_graph", help="Graphmapping by first mapping to a "
                                    "linear reference genome and then fitting these alignments to the graph.")
    subparser_filter = subparsers.add_parser("filter", help="Filtering alignments.")

    for subparser in [subparser_linear]:
        subparser.add_argument("-r", "--reference", help="Linear reference fasta", required=True)
        subparser.add_argument("-f", "--fasta", help="Input fasta file", required=True)
        subparser.add_argument("-d", "--data-dir", help="", required=True)
        subparser.add_argument("-c", "--chromosomes", help="Comma-separated list of chromosomes", required=True)
        subparser.add_argument("-t", "--n_threads", help="Number of threads to use", type=int, default=8, required=False)
        subparser.add_argument("-m", "--min-mapq", help="Minimum mapq to keep", required=False, type=int, default=40)
        subparser.add_argument("-o", "--output-file", help="Write final alignments to this file. If not set, will write to stdout", default=None, required=False)
        subparser.add_argument("-q", "--skip-mapq-adjustment", default=False, required=False,
                help="Set to true to skip running mininap 2 adjust mapq's. Takes less time with worse performance.")


    subparser_remove_from_fasta = subparsers.add_parser("remove_reads_from_fasta", help="Outputs fasta entries with ID not in --alignments file")
    subparser_remove_from_fasta.add_argument("-f", "--fasta", required=True)
    subparser_remove_from_fasta.add_argument("-a", "--alignments", help="File with one alignment per line, id as first column.", required=True)

    subparser_filter.add_argument("-m", "--min-mapq", help="Minimum mapq to keep", required=False, type=int, default=50)
    subparser_filter.add_argument("-a", "--alignments", help="Graphalignments file to filter", required=True)

    subparser_linear.set_defaults(func=run_map_linear_to_graph)
    subparser_filter.set_defaults(func=run_filter)
    subparser_remove_from_fasta.set_defaults(func=remove_reads_from_fasta)


    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)

    if hasattr(args, 'func'):
        args.func(args)

