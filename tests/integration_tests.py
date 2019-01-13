from offsetbasedgraph import Graph, NumpyIndexedInterval, SequenceGraph, Block, Interval
import subprocess
import logging
logging.basicConfig(level=logging.DEBUG)


def create_test_data():
    node1 = "ATTACTACTACTAGGGATCGGGACTATATAACCCATCATTTTACTTTACGGGACTACGAGCGGGAGCTTATCACTACTGACGGGACTACTACTTTACGGAGCTACGAGCATCTAGCGACGACATCAGAATGAGCTATCGTACGTATTATTCTACTATCATGCGTAGCATGTCTTTATCTATTACTATTACTATC"
    node2 = "TTTTTTTT"
    node3 = "GGGGGGGG"
    node4 = "TAAATATATCATACTGCAGCTACGAGCATAGCAGCAGCTACTATCTACGACGATCAGCAGCATCAGCATCTACGGTTGTGATGCTTGGTAGGGACTGGATCTATCTATCATCTATGATGCTGATGCATGCTGATGCTATTATTTTTAAAACGCGGCACTATCTGATTACTATGCGAC"
    linear_ref = node1 + node2 + node4

    graph = Graph({1: Block(len(node1)),
                   2: Block(len(node2)),
                   3: Block(len(node3)),
                   4: Block(len(node4))},
                  {1: [2, 3],
                   3: [4],
                   2: [4]})

    graph.convert_to_numpy_backend()
    graph.to_file("testgraphs/1.nobg")
    sequence_graph = SequenceGraph.create_empty_from_ob_graph(graph)
    sequence_graph.set_sequence(1, node1)
    sequence_graph.set_sequence(2, node2)
    sequence_graph.set_sequence(3, node3)
    sequence_graph.set_sequence(4, node4)
    sequence_graph.to_file("testgraphs/1.nobg.sequences")

    linear_ref_interval = Interval(0, len(node4), [1, 2, 4], graph)
    indexed_interval = linear_ref_interval.to_numpy_indexed_interval()
    indexed_interval.to_file("testgraphs/1_linear_pathv2.interval")

    testreads = [node1[10:160], node1[-100:] + node3 + node4[0:42], node1[-90:] + node3 + node4[0:52]]
    with open("testreads.fa", "w") as readsfile:
        for i, testread in enumerate(testreads):
            readsfile.writelines([">read" + str(i) + "\n", testread + "\n"])

    with open("linear_ref_test.fa", "w") as linear_ref_file:
        linear_ref_file.write(">1\n" + linear_ref + "\n")

    bwa_index_command = "bwa index linear_ref_test.fa"
    output = subprocess.check_output(bwa_index_command.split())
    print(output)


def test_mappers():
    from rough_graph_mapper import LinearToGraphMapper, TraverseMapper
    for mapper in [LinearToGraphMapper, TraverseMapper]:
        mapper("testreads.fa", "linear_ref_test.fa", "testgraphs/", ["1"], minimum_mapq_to_graphalign=60, write_final_alignments_to_file="test.graphalignments")
        with open("test.graphalignments") as results_file:
            alignment1_nodes = Interval.from_file_line(results_file.readline().split("\t")[1]).region_paths
            assert alignment1_nodes == [1]

            alignment2_nodes = Interval.from_file_line(results_file.readline().split("\t")[1]).region_paths
            assert alignment2_nodes == [1, 3, 4]

            alignment3_nodes = Interval.from_file_line(results_file.readline().split("\t")[1]).region_paths
            assert alignment3_nodes == [1, 3, 4]


if __name__ == "__main__":
    create_test_data()
    test_mappers()
