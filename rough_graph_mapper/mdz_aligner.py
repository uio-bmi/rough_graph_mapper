import logging
import pysam
from tqdm import tqdm
from offsetbasedgraph import Graph, NumpyIndexedInterval, SequenceGraph
from collections import defaultdict
from .util import number_of_lines_in_file, read_sam
import pickle


def mdz_align_bam_file(bam_file_name, graph, sequence_graph, reference_path, linear_reference_nodes,
                       out_file_name, limit_to_chromosome=None):

    edge_counts = defaultdict(int)
    #samfile = pysam.AlignmentFile(bam_file_name, "r", check_sq=True)

    #for line in tqdm(samfile.fetch(until_eof=True), total=samfile.count()):
    #n_lines = number_of_lines_in_file(bam_file_name)

    if limit_to_chromosome is None:
        position = 1
    elif limit_to_chromosome == "X":
        position = 23
    else:
        position = int(limit_to_chromosome)

    n_reads_cover_variants = 0
    n_deletions = 0
    n_insertions = 0
    n_substitutions = 0
    i = 0
    for line in tqdm(read_sam(bam_file_name, return_pysam_objects=True), position=position,
                     desc="Chromosome " + str(limit_to_chromosome)):
        if limit_to_chromosome is not None and line.reference_name != limit_to_chromosome:
            continue


        aligner = MdzAligner(line, graph, sequence_graph, reference_path, linear_reference_nodes)
        edges = aligner.align_and_get_variant_nodes()
        if len(edges) > 0:
            n_reads_cover_variants += 1
        for edge in edges:
            edge_counts[edge] += 1

        n_substitutions += aligner.n_subsitutions
        n_deletions += aligner.n_deletions
        n_insertions += aligner.n_insertions

        #logging.info("%d reads processed, %d reads cover variants. Substitutions: %d, insertions: %d, deletions: %d" %
        #             (i, n_reads_cover_variants, n_substitutions, n_insertions, n_deletions))
        i += 1

    logging.info("Chr %s: Found %d variant edges" % (limit_to_chromosome, len(edge_counts.keys())))
    logging.info("Chr %s: Found %d substitutions, %d insertions and %d deletiosn" % (limit_to_chromosome, n_substitutions, n_insertions, n_deletions))
    logging.info("Writing edges to file")
    with open(out_file_name, "wb") as f:
        pickle.dump(edge_counts, f)
        logging.info("Wrote edge counts to file %s" % out_file_name)



# Aligns using the cigar and mdz tag of the sam (bam file). Mainly a wrapper around pysam's method for doing this
# Aims only at discovering nodes relevant to variants in the graph
class MdzAligner:
    def __init__(self, bam_entry, graph, sequence_graph, reference_path, linear_reference_nodes):
        self.bam_entry = bam_entry
        self.graph = graph
        self.sequence_graph = sequence_graph
        self.reference_path = reference_path
        self.linear_reference_nodes = linear_reference_nodes
        self._deletion_nodes_processed = set()
        self._variant_edges_detected = set()
        self.n_subsitutions = 0
        self.n_insertions = 0
        self.n_deletions = 0

    def align_and_get_variant_nodes(self):
        pairs = self.bam_entry.get_aligned_pairs(with_seq=True)
        substitutions = [pair for pair in pairs if pair[2] is not None and pair[2].islower()]
        # Need to get the base before the insertion in order to know where the insertion is relative to the reference
        # NB: This does not handle insertion at the first base in the read
        # NB2: requiring pairs[i][2] is not None ignores insertion right after an insertion, which is what we want
        # since we only care about the first base in the insertion
        insertions = [(pairs[i], pair) for i, pair in enumerate(pairs[1:])
                      if pair[1] is None and pairs[i][2] is not None]
        deletions = [pair for pair in pairs if pair[0] is None]

        for substitution in substitutions:
            self._process_substitution(substitution[0], substitution[1], substitution[2])


        for deletion in deletions:
            did_find = self._process_deletion(deletion[1], deletion[2])

        for insertion in insertions:
            #print("Insertion: " + str(insertion))
            self._process_insertion(insertion[0][1], insertion[1][0])

        return set(self._variant_edges_detected)

    def _process_substitution(self, read_offset, ref_offset, ref_base):
        node = self.reference_path.get_node_at_offset(ref_offset)
        node_offset = self.reference_path.get_node_offset_at_offset(ref_offset)
        if node_offset != 0:
            # We ignore. This is most likely a variant not represented in the graph (or sequencing error)
            return

        prev_node = self.reference_path.get_node_at_offset(ref_offset - 1)
        read_base = self.bam_entry.query_sequence[read_offset]
        #print("Node: %d, offset: %d, node size: %d. Prev node: %d. Read base: %s" %
        #(node, node_offset, graph.blocks[node].length(), prev_node, read_base))

        # Try to find next node that matches read base
        for potential_next in self.graph.adj_list[prev_node]:
            if potential_next == node:
                continue
            node_seq = self.sequence_graph.get_sequence(potential_next, 0, 1)
            #print("  Next node %d has seq %s" % (potential_next, node_seq))
            if node_seq.lower() == read_base.lower():
                # Found a match!
                self._variant_edges_detected.add((prev_node, potential_next))
                self.n_subsitutions += 1
                break

    def _process_deletion(self, ref_offset, ref_base):
        node = self.reference_path.get_node_at_offset(ref_offset)
        if node in self._deletion_nodes_processed:
            return False
        self._deletion_nodes_processed.add(node)
        node_offset = self.reference_path.get_node_offset_at_offset(ref_offset)

        #print("Processing deltion %s, %s, node offset %d" % (ref_offset, ref_base, node_offset))
        if node_offset != 0:
            # Ignore, deletion not in graph
            return

        prev_node = self.reference_path.get_node_at_offset(ref_offset - 1)

        # Using simple heuristic, take next edge that goes back to linear reference
        for potential_next in self.graph.adj_list[prev_node]:
            if potential_next in self.linear_reference_nodes:
                self._variant_edges_detected.add((prev_node, potential_next))
                self.n_deletions += 1
                #print("Found next node %d " % potential_next)
                break

        #print("Node: %d, node offset: %d, node size: %d" % (node, node_offset, graph.blocks[node].length()))
        return True

    def _process_insertion(self, ref_offset, read_offset):
        base = self.bam_entry.query_sequence[read_offset]
        node = self.reference_path.get_node_at_offset(ref_offset)
        node_offset = self.reference_path.get_node_offset_at_offset(ref_offset)
        node_size = self.graph.blocks[node].length()
        if node_offset != node_size - 1:
            # We are not at end of node, insertion is not represented in the graph, ignore
            return False

        # Find out which next node matches the insertion
        for potential_next in self.graph.adj_list[node]:
            if potential_next in self.linear_reference_nodes:
                continue  # Next should not be in linear ref

            #print("Processing insertion at ref offset %d with base %s" % (ref_offset, base))
            #print("  Node %d with offset %d. Node size: %d" % (node, node_offset, self.graph.blocks[node].length()))

            next_base = self.sequence_graph.get_sequence(potential_next, 0, 1).upper()
            if next_base == base.upper():
                self._variant_edges_detected.add((node, potential_next))
                self.n_insertions += 1
                #print("  Found next node %d with seq %s" % (potential_next, next_base))
                return




if __name__ == "__main__":
    graph = Graph.from_file("21.nobg")
    sequence_graph = SequenceGraph.from_file("21.nobg.sequences")
    reference_path = NumpyIndexedInterval.from_file("21_linear_pathv2.interval")
    linear_reference_nodes = reference_path.nodes_in_interval()
    samfile = pysam.AlignmentFile("test2.sorted.bam", "rb")

    for line in tqdm(samfile.fetch(), total=1000000):
        if line.reference_name != "21":
            continue


        aligner = MdzAligner(line, graph, sequence_graph, reference_path, linear_reference_nodes)
        aligner.align_and_get_variant_nodes()


