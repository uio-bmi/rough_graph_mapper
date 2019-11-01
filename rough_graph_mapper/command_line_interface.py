import logging
logging.basicConfig(level=logging.INFO)
import shutil
import argparse
import sys
from .linear_to_graph_mapper import LinearToGraphMapper
from .filter_graphalignments import filter_graphalignments
from offsetbasedgraph import Graph, NumpyIndexedInterval, SequenceGraph
from .mdz_aligner import mdz_align_bam_file
from .util import split_sam_by_chromosomes, select_lowest_mapq_from_two_sam_files, merge_single_line_sams, merge_sams2
from multiprocessing import Process


def merge_sams(args):
    if args.single_line:
        merge_single_line_sams(args.sam1, args.sam2)
    else:
        merge_sams2(args.sam1, args.sam2, scores_are_double=args.scores_are_double, only_score_lowering=args.only_score_lowering)
        #improve_mapping_with_two_sams(args.sam1, args.sam2)
    #logging.info("Merging ...")
    #select_lowest_mapq_from_two_sam_files(args.sam1, args.sam2)


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

def mdz_align_single_chromosome(graph_dir, sam_file, chromosome, output_file):
    out_file_name = output_file + "_" + chromosome + ".edgecounts"
    sam_file_name = sam_file  # sam_base_name + "_chr" + chromosome + ".sam"

    graph = Graph.from_file(graph_dir + "/" + chromosome + ".nobg")
    sequence_graph = SequenceGraph.from_file(graph_dir + "/" + chromosome + ".nobg.sequences")
    reference_path = NumpyIndexedInterval.from_file(graph_dir + "/" + chromosome + "_linear_pathv2.interval")
    linear_reference_nodes = reference_path.nodes_in_interval()

    mdz_align_bam_file(sam_file_name, graph, sequence_graph, reference_path, linear_reference_nodes,
                       out_file_name, limit_to_chromosome=chromosome)

def mdz_align_wrapper(args):
    graph_dir = args.data_dir
    bam_file = args.bam_file
    sam_base_name = ".".join(bam_file.split(".")[0:-1])

    #split_sam_by_chromosomes(bam_file, args.chromosomes.split(","))

    processes = []
    for chromosome in args.chromosomes.split(","):
        p = Process(target=mdz_align_single_chromosome, args=(graph_dir, bam_file, chromosome, args.output_file))
        p.start()

    for p in processes:
        p.join()



def run_filter(args):
    alignments = args.alignments
    min_mapq = args.min_mapq
    filter_graphalignments(alignments, min_mapq)

def run_map_linear_to_graph(args):
    chromosomes = args.chromosomes.split(",")
    LinearToGraphMapper(args.fasta, args.reference, args.data_dir, chromosomes, n_threads=args.n_threads,
                        write_final_alignments_to_file=args.output_file, skip_mapq_adjustment=args.skip_mapq_adjustment,
                        minimum_mapq_to_graphalign=args.min_mapq, only_run_minimap=args.only_run_minimap)

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
        subparser.add_argument("-M", "--only-run-minimap", help="Only run Minimap (no BWA-MEM)", default=False, required=False)
        subparser.add_argument("-q", "--skip-mapq-adjustment", default=False, required=False,
                help="Set to true to skip running mininap to adjust mapq's. Takes less time with worse performance.")


    subparser_remove_from_fasta = subparsers.add_parser("remove_reads_from_fasta", help="Outputs fasta entries with ID not in --alignments file")
    subparser_remove_from_fasta.add_argument("-f", "--fasta", required=True)
    subparser_remove_from_fasta.add_argument("-a", "--alignments", help="File with one alignment per line, id as first column.", required=True)

    subparser_filter.add_argument("-m", "--min-mapq", help="Minimum mapq to keep", required=False, type=int, default=50)
    subparser_filter.add_argument("-a", "--alignments", help="Graphalignments file to filter", required=True)

    subparser_linear.set_defaults(func=run_map_linear_to_graph)
    subparser_filter.set_defaults(func=run_filter)
    subparser_remove_from_fasta.set_defaults(func=remove_reads_from_fasta)

    # mdz align
    subparser_mdzalign = subparsers.add_parser("mdz_align_bam", help="Align a bam file to graph by using cigar/mdz align.")
    subparser_mdzalign.add_argument("-b", "--bam-file", help="Input bam file", required=True)
    subparser_mdzalign.add_argument("-d", "--data-dir", help="Directory containing graphs", required=True)
    subparser_mdzalign.add_argument("-o", "--output-file", help="Write edge counts to this file "
                                                                "(will be used as a base file name)", required=True)
    subparser_mdzalign.add_argument("-c", "--chromosomes", help="Comma-separated list of chromosomes", required=True)
    subparser_mdzalign.set_defaults(func=mdz_align_wrapper)

    # merge_sams
    subparser_merge_sams = subparsers.add_parser("merge_sams", help="Improve mapping performance by merging two SAMs.")
    subparser_merge_sams.add_argument("sam1", help="Sam file 1")
    subparser_merge_sams.add_argument("sam2", help="Sam file 2")
    subparser_merge_sams.add_argument("-t", "--scores-are-double", required=False, type=bool, default=False, help="Set to True if sam2 scores are double (default with minimap)")
    subparser_merge_sams.add_argument("-r", "--only-score-lowering", required=False, type=bool, default=False)
    subparser_merge_sams.add_argument("-s", "--single-line", required=False, type=bool, default=False,
                                      help="Set to True if both same files have only one line per read. Is faster")
    subparser_merge_sams.set_defaults(func=merge_sams)

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)

    if hasattr(args, 'func'):
        args.func(args)

