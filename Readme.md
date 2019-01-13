


# A simple proof-of-concept graph mapper
Rough Graph Mapper is a simple proof-of-concept graph-mapper implemented in Python. 

## Install
The only requirement is BWA. Makes sure that you have BWA installed.

You install by cloning and running `python setup.py`

```bash
git clone ...
cd rough_graph_mapper
python3 setup.py install
```

## How to run
You will need a directory with graphs and the linear reference of which the graphs has been built from.
### Linear to graph mapping
This is the simplest and fastest form of mapping. 
A combination of BWA MEM and Minimap is run against the linear reference genome 
first, and then these alignments are "fitted" to the graph.
```bash
rough_graph_mapper map_linear_to_graph -r linear_reference.fa -f reads.fa -d graphs_dir/ --chromosomes 1,2,3 > mapped.graphalignments
```

### Traversemapper
This method first runs linear to graph mapping, and then tries to map the rest of the reads 
(those that did not get a good mapping in the first step) by traversing the graph and fitting the graph to the reads.
This method is a lot slower and requires a lot of memory.
```bash
rough_graph_mapper traversemapper -r linear_reference.fa -f reads.fa -d graphs_dir/ --chromosomes 1,2,3 > mapped.graphalignments
```

### Filter the final alignments
This is a super simple filtering script, aimed at removing bad alignments and alignments from reads that multimap.
```bash
rough_graph_mapper filter -r mapped.graphalignments --min-mapq 50
```
