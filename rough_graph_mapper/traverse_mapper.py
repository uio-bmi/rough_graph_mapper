from multiprocessing import Process
from .linear_to_graph_mapper import LinearToGraphMapper
from .util import read_sam, read_fasta_to_numeric_sequences
import os
import mappy as mp
from tqdm import tqdm
import logging
from collections import defaultdict
from offsetbasedgraph import Interval, SequenceGraph, Graph, NumpyIndexedInterval
from .traverser import traverse

def process_linear_reads_fitted_to_graph(file_name):
    # Todo: Only add reads on this chromosome (with node ids within graph range)
    edge_counts = defaultdict(int)
    read_ids = set()

    f = open(file_name)
    for line in f:
        l = line.split("\t")
        read_ids.add(l[0])
        if l[1] == "None":
            continue
        alignment = Interval.from_file_line(l[1])
        rps = alignment.region_paths
        for j in range(0, len(rps) - 1):
            edge_counts["%d-%d" % (rps[j], rps[j+1])] += 1

    return edge_counts, read_ids


def map_single_chromosome(base_name, graph_dir, chromosome, linear_ref_bonus=1, seed_length=100):
    sequence_graph = SequenceGraph.from_file(graph_dir + chromosome + ".nobg.sequences")
    graph = Graph.from_file(graph_dir + chromosome + ".nobg")
    edge_counts, read_ids = process_linear_reads_fitted_to_graph(base_name + "_chr" + chromosome + ".sam.graphalignments")
    aligner = mp.Aligner(base_name + "_without_those_aligned_well_linear.fa", n_threads=2, k=19, w=1, best_n=300, min_cnt=2, preset="sr")
    reads = read_fasta_to_numeric_sequences(base_name + "_without_those_aligned_well_linear.fa", sequence_graph)
    linear_path = NumpyIndexedInterval.from_file(graph_dir + "/%s_linear_pathv2.interval" % chromosome)
    linear_path_nodes = linear_path.nodes_in_interval()

    traverse(chromosome, graph, aligner, linear_path_nodes, sequence_graph, reads, seed_length=seed_length,
                        linear_ref_bonus=linear_ref_bonus, edge_counts=edge_counts)


class TraverseMapper:
    def __init__(self, fasta_file_name, linear_reference_file_name, graph_dir, chromosomes,
                 minimum_mapq_to_graphalign=60, skip_run_linear_to_graph=False, write_final_alignments_to_file=None,
                 n_threads=8, skip_mapq_adjustment=False):
        self.fasta_file_name = fasta_file_name
        self.base_name = '.'.join(self.fasta_file_name.split(".")[:-1])
        self.chromosomes = chromosomes
        self.graph_dir = graph_dir
        self.write_final_alignments_to_file = write_final_alignments_to_file

        # First run linear to graph mapper
        if not skip_run_linear_to_graph:
            LinearToGraphMapper(fasta_file_name, linear_reference_file_name, graph_dir, chromosomes,
                                minimum_mapq_to_graphalign=minimum_mapq_to_graphalign,
                                write_final_alignments_to_file=self.base_name + "_from_linear.graphalignments",
                                n_threads=n_threads, skip_mapq_adjustment=skip_mapq_adjustment)
        else:
            logging.info("Will use %s_chr*.sam.graphalignments as linear rads." % self.base_name)
            for chromosome in self.chromosomes:
                assert os.path.isfile(self.base_name + "_chr" + chromosome + ".sam.graphalignments")

        #assert os.path.isfile(self.base_name + "_from_linear.graphalignments"), \
        #   "Graphalignments from linear alignments does not exist. Run without -s or run map_linear_to_graph first."

        self._create_fasta_without_those_aligned_well_linear()
        self.run_graphtraverse_mapping()

    def _create_fasta_without_those_aligned_well_linear(self):
        ignore_reads = set()
        logging.info("Removing reads that mapped well to linear ref")
        # We want to add all that had mapq 60 originally, because we have attempted to fit these to graph
        for record in read_sam(self.base_name + ".sam"):
            if record.mapq == 60 and record.score > 0:
                ignore_reads.add(record.name)

        for chromosome in tqdm(self.chromosomes):
            with open(self.base_name + "_chr" + str(chromosome) + ".sam.multimapping.txt") as f:
                for line in f:
                    ignore_reads.add(line.strip())

        logging.info("Will ignore %d reads when graphmapping" % len(ignore_reads))

        # Create new fasta without these reads
        skip = False
        fasta = open(self.fasta_file_name)
        new_fasta = open(self.base_name + "_without_those_aligned_well_linear.fa", "w")
        for line in fasta:
            if skip:
                skip = False
                continue

            if line.startswith(">"):
                if line.replace(">", "").strip() in ignore_reads:
                    skip = True
                    continue

            new_fasta.writelines([line])

        new_fasta.close()

    def run_graphtraverse_mapping(self):
        processes = []
        for chromosome in self.chromosomes:
            process = Process(target=map_single_chromosome,
                              args=(self.base_name, self.graph_dir, chromosome))
            process.start()

            processes.append(process)

        for p in processes:
            p.join()

        out_file = None
        if self.write_final_alignments_to_file is not None:
            out_file = open(self.write_final_alignments_to_file, "w")

        logging.info("Done with all chromosomes. Merging alignments with those from linear")
        for chromosome in self.chromosomes:
            with open(chromosome + ".graphalignments") as f:
                for line in f:
                    if self.write_final_alignments_to_file is None:
                        print(line.strip())
                    else:
                        out_file.writelines([line])

            with open(self.base_name + "_chr" + chromosome + ".sam.graphalignments") as f:
                for line in f:
                    if self.write_final_alignments_to_file is None:
                        print(line.strip())
                    else:
                        out_file.writelines([line])

        if self.write_final_alignments_to_file is not None:
            logging.info("Wrote all final graphalignments to %s" % self.write_final_alignments_to_file)
            out_file.close()
        logging.info("Done mapping reads")
