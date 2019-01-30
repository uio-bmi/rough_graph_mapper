from .single_read_aligner import SingleSequenceAligner
import numpy as np
import mappy as mp
import logging
from collections import defaultdict
from tqdm import tqdm

DEBUG_ID = "2178a3e2d74e352d"

def traverse(chrom, graph, aligner, linear_path_nodes, sequence_graph, reads, seed_length=100, linear_ref_bonus=2, n_mismatches_allowed=6, edge_counts=None):
    first_nodes = graph.get_first_blocks()
    assert len(first_nodes) == 1
    traversed_sequence = np.zeros(300000000, dtype=np.uint8)
    n_traversed = 150

    out_file = open(chrom + ".graphalignments", "w")

    node = first_nodes[0]
    #node = 170000
    n_alignments = 0
    n_skipped_q_end_short = 0
    n_skipped_start_mismatching = 0
    n_skipped_because_attempted_recently = 0
    n_hit_read_end = 0
    path = []
    i = 0
    reads_aligned = {}  # read id to node id
    best_alignments = {}  # read hit, number of mismatches
    previous_alignment_attempt = defaultdict(int)
    n_times_aligned = defaultdict(int)  # read id to n times aligned
    n_attempted_aligned = 0
    prev_offset_attempted = 0

    logging.info("Allowing %d mismatches " % n_mismatches_allowed)

    if chrom == "X":
        progress_position = 23
    else:
        progress_position = int(chrom)

    for i in tqdm(range(len(graph.blocks)), total=len(linear_path_nodes), position=progress_position, desc="Chromosome " + chrom):
        path.append(node)
        if i % 500000 == 0 and chrom == "20":
            n_nodes_not_in_linear = len(path) - len(set(path).intersection(linear_path_nodes))
            logging.info \
                ("Chrom %s: %d nodes processed. n aligned: %d (%d unique reads). N attempted aligned: %d. N nodes not in linear ref so far: %d. N skipped because q end is too short: %d. N hitting far end of read: %d. N skipped because start is mismatching: %d. N alignment attempts saved because tried recently to align and skipped: %d" %
                (chrom, i, n_alignments, len(reads_aligned), n_attempted_aligned, n_nodes_not_in_linear,
                n_skipped_q_end_short, n_hit_read_end, n_skipped_start_mismatching,
                n_skipped_because_attempted_recently))
        i += 1


        # Add this node's sequence
        assert node is not None
        node_sequence = sequence_graph.get_numeric_node_sequence(node)
        traversed_sequence[n_traversed:n_traversed + len(node_sequence)] = node_sequence
        n_traversed += len(node_sequence)

        next_nodes = graph.adj_list[node]

        if len(next_nodes) == 0:
            logging.info("Done traversing")
            break

        # Get a seed from this linear sequence
        traversed_sequence_numeric = traversed_sequence[n_traversed - seed_length:n_traversed]
        sub_sequence = ''.join(sequence_graph._letters[traversed_sequence_numeric])

        if node > 170242 and node < 170242 + 10:
            print("Node: %d, n_traversed: %d, sub_sequence: %s" % (node, n_traversed, sub_sequence.upper()))
            print("prev attempt: %d" % prev_offset_attempted)

        n_hits = 0

        # Ignore, we have aligned stuff recently
        if prev_offset_attempted < n_traversed - 5 and linear_ref_bonus < 100:  # Hack when testing, don't try to find next node if we want linear ref anyway
            prev_offset_attempted = n_traversed
                
            if node == 170249:
                logging.info("Trying to align %s" % sub_sequence)

            # Align this seed to the reads
            for read_hit in aligner.map(sub_sequence):

                if node == 170249:
                    print(read_hit)

                if read_hit.ctg == DEBUG_ID:
                    logging.info("Hit on node %d: %s " % (node, read_hit))

                if n_times_aligned[read_hit.ctg] > 50:
                    # Ignore when we already have more than 20 matches for read. This read is probably multimapping a lot, and we don't want to use it to predict the path
                    # This is only for speedup
                    continue

                if read_hit.ctg in reads_aligned:
                    if reads_aligned[read_hit.ctg] > n_traversed - 100:
                        continue  # Ignore if already aligned to a very recent node. This means we are probably on the same alignment

                if previous_alignment_attempt[read_hit.ctg] > 0 and previous_alignment_attempt[
                    read_hit.ctg] > n_traversed - 40:
                    n_skipped_because_attempted_recently += 1
                    continue
                n_hits += 1
                score = seed_length - read_hit.NM
                if read_hit.NM > n_mismatches_allowed:
                    if read_hit.ctg == DEBUG_ID:
                        logging.info("Skipped because too low score")
                    continue  # Skip, too low score

                # Ingore if short match. If short match, we know we will probably get this match on next node again, so we'll wait until we have more sequence matching in order to
                # let minimap to more of the job for us
                if read_hit.q_en - read_hit.q_st < seed_length - 32:
                    if read_hit.ctg == DEBUG_ID:
                        logging.info("Short match, ignore")
                    continue

                full_read_sequence, numeric_read_sequence = reads[read_hit.ctg]
                if read_hit.strand < 0:
                    full_read_sequence = mp.revcomp(full_read_sequence)
                    # TODO: THis can be done much more efficient than converting the text sequence again:
                    numeric_read_sequence = sequence_graph._letter_sequence_to_numeric(
                        np.array(list(full_read_sequence.lower())))
                    read_start = len(full_read_sequence) - read_hit.r_en
                    read_end = len(full_read_sequence) - read_hit.r_st
                else:
                    read_start = read_hit.r_st
                    read_end = read_hit.r_en

                # Ignore when query end is not seed length, then we don't have match towards end of seed
                if read_hit.q_en + (5 - read_hit.NM) < seed_length:
                    n_skipped_q_end_short += 1
                    continue

                # Ignore when query start and read start is high, means we have a big mismatch in beginning
                if read_hit.q_st > (n_mismatches_allowed - read_hit.NM) and read_start > (
                        n_mismatches_allowed - read_hit.NM):
                    continue

                # (Don't) Ignore if read end is as long as the read, this means we don't have any more read left to predict the path using this read
                if read_end == len(full_read_sequence):
                    if read_hit.ctg == DEBUG_ID:
                        logging.info("Ignoring, hitting read end")
                    n_hit_read_end += 1

                n_mismatches_so_far = read_hit.NM

                # If read start is higher than zero, first check that we have an okay match backwards with our traversed path
                if read_start > 0:
                    pre_sequence_read = numeric_read_sequence[0:read_start]
                    # pre_sequence_read = sequence_graph._letter_sequence_to_numeric(np.array(list(pre_sequence_read.lower())))
                    path_start = int(n_traversed - seed_length + read_hit.q_st - read_start)
                    path_end = int(path_start + read_start)
                    path_sequence = traversed_sequence[path_start:path_end]
                    if read_hit.ctg == DEBUG_ID:
                        logging.info("Trying to match starts:")
                        logging.info(full_read_sequence)
                        logging.info(numeric_read_sequence)
                        logging.info(pre_sequence_read)
                        logging.info(path_sequence)
                    n_mismatches = np.count_nonzero(path_sequence != pre_sequence_read)
                    # logging.info("N mismatches: %d" % n_mismatches)
                    if n_mismatches > (n_mismatches_allowed - read_hit.NM):
                        n_skipped_start_mismatching += 1
                        continue

                    n_mismatches_so_far += n_mismatches

                if read_hit.ctg == DEBUG_ID:
                    logging.info("Read start %d, read end %d, q start %d, q end %d" % (
                    read_hit.r_st, read_hit.r_en, read_hit.q_st, read_hit.q_en))
                    logging.info("Reversed seq:\n%s" % full_read_sequence.upper())
                    logging.info("Seed     seq:\n%s" % sub_sequence.upper())
                    logging.info("Read start: %d" % read_start)
                    logging.info("Read end: %d" % read_end)
                    logging.info("Seed start: %d" % read_hit.q_st)
                    logging.info("Seed end: %d" % read_hit.q_en)

                # print(read_hit.q_en - read_hit.q_st)
                # Fine-align this match, we want to start aligning one basepair before end of node we are at
                sub_sequence_start = read_hit.q_st
                dist_to_node_end = seed_length - sub_sequence_start

                read_start_align_pos = read_start + dist_to_node_end - 1
                # read_sequence_to_align = sequence_graph._letter_sequence_to_numeric(np.array(list(full_read_sequence[read_start_align_pos:].lower())))
                read_sequence_to_align = full_read_sequence[read_start_align_pos:]
                numeric_read_sequence_to_align = numeric_read_sequence[read_start_align_pos:]

                # print("Will try to align %s starting from end of node %d" % (read_sequence_to_align, node))
                print_debug = False
                if read_hit.ctg == DEBUG_ID:
                    print_debug = True
                    logging.info("Starting with %d mismatches" % n_mismatches_so_far)
                    logging.info("Previous nodes in path: %s" % str([("%d-%d" % (n1, n2), edge_counts["%d-%d" % (n1, n2)]) for n1, n2 in zip(path[-6:-1], path[-5:])]))

                simple_aligner = SingleSequenceAligner(graph, sequence_graph, node, graph.blocks[node].length() - 1,
                                                       numeric_read_sequence_to_align,
                                                       n_mismatches_allowed=n_mismatches_allowed,
                                                       n_mismatches_init=n_mismatches_so_far, print_debug=print_debug)
                simple_aligner.align()
                alignment = simple_aligner.get_alignment()
                n_attempted_aligned += 1
                previous_alignment_attempt[read_hit.ctg] = n_traversed
                if not alignment:
                    if read_hit.ctg == DEBUG_ID:
                        logging.info("Not aligned, mismatching")

                    continue
                    # print(" No valid alignment")
                else:
                    assert simple_aligner.n_mismatches_so_far <= simple_aligner.n_mismatches_allowed
                    if read_hit.ctg == DEBUG_ID:
                        logging.info("ALIGNED !!!! n mismatches total: %d" % simple_aligner.n_mismatches_so_far)
                    n_alignments += 1
                    reads_aligned[read_hit.ctg] = n_traversed
                    n_times_aligned[read_hit.ctg] += 1

                    out_file.writelines(["%s\t%s\t%d\n" % (read_hit.ctg, alignment.to_file_line(), simple_aligner.n_mismatches_so_far)])
                    # Add edge counts to all edges followed
                    # Only add edge counts if this is the best alignment for this read
                    if read_hit.ctg not in best_alignments or (best_alignments[read_hit.ctg] >= simple_aligner.n_mismatches_so_far):
                        best_alignments[read_hit.ctg] = simple_aligner.n_mismatches_so_far
                        rps = alignment.region_paths
                        for j in range(0, len(alignment.region_paths) - 1):
                            if rps[j] == 170242:
                                logging.info("        ADDING 1 to %d - %d for alignment %s for read %s" % (rps[j], rps[j+1], alignment, read_hit.ctg))
                            edge_counts["%d-%d" % (rps[j], rps[j + 1])] += 1

        # Choose next node based on edge counts
        highest_count = 0
        best_next = None
        if len(next_nodes) == 1:
            best_next = next_nodes[0]
        else:

            for potential_next in next_nodes:
                count = edge_counts["%d-%d" % (node, potential_next)]
                if potential_next in linear_path_nodes:
                    count += linear_ref_bonus

                if count > highest_count or count == highest_count and potential_next in linear_path_nodes:
                    highest_count = count
                    best_next = potential_next

        node = best_next

    logging.info("Done traversing graph")

    nodes_chosen = set(path)
    n_on_linear = len(nodes_chosen.intersection(linear_path_nodes))
    n_not_on_linear = len(nodes_chosen) - n_on_linear

    # logging.info("N ambigious choices: %d" % n_ambigious)
    logging.info("Total nodes in linear ref: %d" % len(linear_path_nodes))
    logging.info("N nodes chosen that are not in linear ref: %d " % n_not_on_linear)
    logging.info("N nodes chosen that are in linear ref: %d " % n_on_linear)

    # logging.info("N special case: %d" % n_special_case)

    logging.info("N nodes in path: %d" % len(path))

    return path

