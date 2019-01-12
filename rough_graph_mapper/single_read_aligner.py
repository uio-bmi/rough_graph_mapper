import numpy as np
from offsetbasedgraph import Graph, SequenceGraph, Interval
import logging


class SingleSequenceAligner:
    def __init__(self, graph, sequence_graph, node, offset, sequence, n_mismatches_allowed=0, n_mismatches_init=0,
                 print_debug=False):
        assert isinstance(sequence, np.ndarray), "Sequence must be numpy array"
        self.graph = graph
        self.sequence_graph = sequence_graph
        self.node = node
        self.offset = offset
        self.sequence = sequence
        self.n_mismatches_allowed = n_mismatches_allowed

        self._current_node_offset = offset
        self.sequence_offset = 0
        self.nodes = []
        self.start_interval_offset = offset
        self.sequence_length = len(sequence)
        self._last_node = None
        self.print_debug = print_debug
        self.n_mismatches_so_far = n_mismatches_init

    @staticmethod
    def count_mismatches_between_sequences(seq1, seq2):
        return np.count_nonzero(seq1 != seq2)

    def get_alignment(self):
        if self._last_node is None:
            return False
        return Interval(self.start_interval_offset, self.graph.blocks[self._last_node].length(), self.nodes)

    def _node_matches(self, node):
        assert node > 0, "Only positive nodes supported in numeric"
        if self.print_debug:
            logging.debug(" Checking if node match against node %d" % node)
        node_size = self.graph.blocks[node].length()
        rest_of_node_sequence = self.sequence_graph.get_numeric_node_sequence(node)[
                                self._current_node_offset:node_size]

        subsequence = self.sequence[self.sequence_offset:self.sequence_offset + (node_size - self._current_node_offset)]
        rest_of_node_sequence = rest_of_node_sequence[0:len(subsequence)]
        if self.print_debug:
            logging.debug("     Rest of node sequence for node %d: %s" % (node, rest_of_node_sequence))
        if self.print_debug:
            logging.debug("     Subsequence:                       %s" % subsequence)
        n_mismatches = SingleSequenceAligner.count_mismatches_between_sequences(rest_of_node_sequence, subsequence)
        if self.print_debug:
            logging.debug("     N mismatches: %d" % n_mismatches)
        return n_mismatches

        # if rest_of_node_sequence[0:len(subsequence)] != subsequence:
        #    return False

        # return True

    def align(self):
        current_node = self.node

        start_interval_offset = self.offset
        n_mismatches = self._node_matches(current_node)
        if n_mismatches > self.n_mismatches_allowed:
            if self.print_debug:
                logging.debug("First node mismatch")
            return False

        self.n_mismatches_so_far += n_mismatches

        while True:
            # Find best next node
            self.nodes.append(current_node)
            self.sequence_offset += (self.graph.blocks[current_node].length() - self._current_node_offset)
            if self.print_debug:
                logging.debug("    sequence offset: %d" % self.sequence_offset)
            self._current_node_offset = 0

            if self.sequence_offset >= self.sequence_length:
                if self.print_debug:
                    logging.debug("Not searching more since sequence offset >= sequence length")
                break

            next_node = None
            adj_list = self.graph.adj_list
            if current_node < 0:
                adj_list = self.graph.reverse_adj_list

            if self.print_debug:
                logging.debug("Current node: %d" % current_node)
            fewest_mismatches = 1000
            best_next_node = None
            for potential_next in adj_list[current_node]:
                if self.print_debug:
                    logging.debug("Potential next node: %d" % potential_next)
                n_mismatches = self._node_matches(potential_next)
                if n_mismatches == 0:
                    # Choose this and don't care about rest
                    fewest_mismatches = 0
                    best_next_node = potential_next
                    if self.print_debug:
                        logging.debug("Chose next node as %d" % best_next_node)
                    break

                if self.print_debug:
                    logging.debug("    N mismatches: %d, fewest: %d, so far: %d, allowed: %d" % (
                    n_mismatches, fewest_mismatches, self.n_mismatches_so_far, self.n_mismatches_allowed))
                if n_mismatches < fewest_mismatches \
                        and n_mismatches + self.n_mismatches_so_far <= self.n_mismatches_allowed:
                    fewest_mismatches = n_mismatches
                    best_next_node = potential_next
                    if self.print_debug:
                        logging.debug(
                            "    Found best next node as %d with %d mismatches" % (best_next_node, n_mismatches))
                else:
                    if self.print_debug:
                        logging.debug("    Not accepting because too many mismatches")
                    pass

            self.n_mismatches_so_far += fewest_mismatches
            next_node = best_next_node
            if next_node is None:
                # print("Did not find next node")
                return False

            current_node = next_node

        self._last_node = current_node


