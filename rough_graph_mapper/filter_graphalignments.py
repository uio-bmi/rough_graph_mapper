from collections import defaultdict
import logging


def filter_graphalignments(file_name, min_mapq=50):
    f = open(file_name)
    alignments = defaultdict(list)  # read id to list of alignments

    for i, line in enumerate(f):
        if i % 1000000 == 0:
            logging.info("Read %d lines" % i)
        l = line.split("\t")
        alignments[l[0]].append((l[1], int(l[2])))

    n_single = 0
    n_skipped_multiple_good = 0
    n_skipped_multiple_good_unique = 0
    n_kept_multiple_good = 0
    n_kept_multiple_good_unique = 0
    n_filtered_bad_score = 0
    n_filtered_linear_many_alternative = 0

    for read_name, hits in alignments.items():
        sorted_hits = sorted(hits, key=lambda x: x[1])
        if sorted_hits[0][1] >= 7:
            # Skip because too bad score
            n_filtered_bad_score += 1
            continue

        if len(sorted_hits) == 1:
            best_hit = sorted_hits[0]
            score = 60 - best_hit[1] * 2
            n_single += 1
        else:
            # Only keep if second best has much more mismatches
            n_mismatches = sorted_hits[0][1]
            if sorted_hits[1][1] <= n_mismatches + 0:
                # skip
                n_skipped_multiple_good_unique += 1
                n_skipped_multiple_good += len(sorted_hits)
                continue
            else:
                best_hit = sorted_hits[0]
                n_kept_multiple_good_unique += 1
                n_kept_multiple_good += len(sorted_hits)

            score = 5
            if len(hits) == 2:
                score = 58
            elif len(hits) == 3:
                score = 56
            elif len(hits) >= 4:
                score = 54

            score -= best_hit[1] * 1 + (sorted_hits[1][1] - sorted_hits[0][1])

        if score >= min_mapq:
            print("%s\t%s\t%s" % (read_name, best_hit[0], score))

    logging.info("%d reads filtered because multimapping on linear ref" % n_filtered_linear_many_alternative)
    logging.info("%d reads filtered because too low score" % n_filtered_bad_score)
    logging.info("%d reads have only a single good alignment" % n_single)
    logging.info("%d reads filtered out because they have multiple good mappings (%d in total)" %
                 (n_skipped_multiple_good_unique, n_skipped_multiple_good))
    logging.info("%d reads with multiple alignments kept, because second best is bad (%d in total)" %
                 (n_kept_multiple_good_unique, n_kept_multiple_good))


