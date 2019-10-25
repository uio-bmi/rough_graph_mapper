


# A simple proof-of-concept graph mapper
Rough Graph Mapper is a simple proof-of-concept naive graph-mapper implemented in Python, that align reads to a graph
by first aligning them to a linear reference genome and then locally aligning using [GSSW](https://github.com/vgteam/gssw).

## Install
You will need BWA and Minimap2 installed first. The rest of the dependencies are 
handled when you run the python install script.

You install by using pip:

```bash
pip3 install rough_graph_mapper
```

Also, add the scripts dir to your path by adding the following line to your `~/.profile` file:
```
export PATH="$PATH:/path/to/rough_graph_mapper/scripts/"
```



## How to run
You will need a directory with graphs and the linear reference of which the graphs has been built from.

A combination of BWA MEM and Minimap is run against the linear reference genome 
first, and then these alignments are "fitted" to the graph.
```bash
rough_graph_mapper map_linear_to_graph -r linear_reference.fa -f reads.fa -d graphs_dir/ --chromosomes 1,2,3 > mapped.graphalignments
```


### Filter the final alignments
This is a super simple filtering script, aimed at removing bad alignments and alignments from reads that multimap.
```bash
rough_graph_mapper filter -r mapped.graphalignments --min-mapq 50
```
