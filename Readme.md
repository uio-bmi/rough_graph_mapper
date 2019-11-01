


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
The linear reference will have to be indexed with BWA-MEM version 2. If you want to use BWA-MEM version 1, just modify
scripts/map_linear file.

A combination of BWA MEM and Minimap is run against the linear reference genome 
first, and then these alignments are "fitted" to the graph.
```bash
map_linear number_of_threads bwa_index_base.fa reference_genome.fa reads.fa > mapped.sam

# We now have a sam file that we can fit to the graph
# $chromosomes should be a comma-separated list of chromosomes, and obg_graph_dir the directory containing the graphs
rough_graph_mapper mdz_align_bam -b mapped.sam -d $obg_graph_dir -o mdzaligned -c $chromosomes

```

