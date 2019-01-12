import numpy as np
from tqdm import tqdm
from .util import number_of_lines_in_file
import logging
from offsetbasedgraph import Graph, SequenceGraph, NumpyIndexedInterval
import mappy as mp
from .single_read_aligner import SingleSequenceAligner


class SamToGraphAligner:
    def __init__(self, graph_dir, chromosome, sam_file_name):
        self.graph_dir = graph_dir
        self.chromosome = chromosome
        self.sam_file_name = sam_file_name

        self.n_skipped_supplementary = 0
        self.n_skipped_low_mapq = 0
        self.n_aligned = 0
        self.n_did_not_align = 0
        self._read_graph_data()

    def _read_graph_data(self):
        chromosome = self.chromosome
        graph = Graph.from_file(self.graph_dir + chromosome + ".nobg")
        self.graph = graph
        self.sequence_graph = SequenceGraph.from_file(self.graph_dir + chromosome + ".nobg.sequences")
        linear_path = NumpyIndexedInterval.from_file(self.graph_dir + chromosome + "_linear_pathv2.interval")
        self.linear_path = linear_path

    def _align_sam_line(self, line):
        l = line.split()
        chrom = l[2]

        flag = int(l[1])
        if flag == 2048 or flag == 2064:
            self.n_skipped_supplementary += 1
            return

        mapq = int(l[4])
        read_name = l[0]

        if mapq < 60:
            self.n_skipped_low_mapq += 1
            return

        sequence = l[9]
        linear_pos = int(l[3]) - 1
        position = self.linear_path.position_at_offset(linear_pos)

        if l[1] == "16":
            # Is mapped to reverse strand
            sequence = mp.revcomp(sequence)
            #logging.warning("Found reverse alignment. Not implemented now, ignoring")

        sequence = self.sequence_graph._letter_sequence_to_numeric(np.array(list(sequence.lower())))
        # logging.info("Seq: %s" % sequence)
        aligner = SingleSequenceAligner(self.graph, self.sequence_graph, position.region_path_id,
                                        int(position.offset), sequence, n_mismatches_allowed=10, print_debug=False)
        aligner.align()
        a = aligner.get_alignment()
        if a:
            self.n_aligned += 1
            print("%s\t%s\t%d" % (read_name, a.to_file_line(), aligner.n_mismatches_so_far))
        else:
            self.n_did_not_align += 1

    def align_sam(self):
        f = open(self.sam_file_name)
        if self.chromosome == "X":
            progress_position = 23
        else:
            progress_position = int(self.chromosome)

        for i, line in enumerate(tqdm(f, total=number_of_lines_in_file(self.sam_file_name), position=progress_position)):
            if line.startswith("@"):
                continue
            self._align_sam_line(line)

            #if i % 1000 == 0:
            #    logging.info("%d sam records processed. N aligned: %d, N skipped low mapq: %d. N did not align: %d" %
            #    (i, self.n_aligned, self.n_skipped_low_mapq, self.n_did_not_align))


