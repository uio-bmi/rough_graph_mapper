import logging
logging.basicConfig(level=logging.INFO)
import pysam
from tqdm import tqdm
from offsetbasedgraph import Graph, NumpyIndexedInterval, SequenceGraph


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

    def align_and_get_variant_nodes(self):
        pairs = self.bam_entry.get_aligned_pairs(with_seq=True)
        substitutions = [pair for pair in pairs if pair[2] is not None and pair[2].islower()]
        insertions = [pair for pair in pairs if pair[1] is None]
        deletions = [pair for pair in pairs if pair[0] is None]

        for substitution in substitutions:
            self._process_substitution(substitution[0], substitution[1], substitution[2])


        for deletion in deletions:
            did_find = self._process_deletion(deletion[1], deletion[2])
            if not did_find:
                break

        return self._variant_edges_detected

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
        for potential_next in graph.adj_list[prev_node]:
            if potential_next == node:
                continue
            node_seq = self.sequence_graph.get_sequence(potential_next, 0, 1)
            #print("  Next node %d has seq %s" % (potential_next, node_seq))
            if node_seq.lower() == read_base:
                # Found a match!
                self._variant_edges_detected.add((prev_node, potential_next))
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
                #print("Found next node %d " % potential_next)
                break

        #print("Node: %d, node offset: %d, node size: %d" % (node, node_offset, graph.blocks[node].length()))

        return True

    def _process_insertion(self, insertion):
        pass



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


