# Barcode error-correction

C++ code for sequencing error correction of fixed-length DNA or RNA sequences (barcodes).  
This is useful for creating a custom feature library for single-cell RNA-seq data analysis with feature barcoding.  See [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/feature-bc) for more information.

## Description

Extract groups of sequences sharing the same barcode using a graph approach.  
The assumption is that the reads sequenced from the same barcode B differ from each other by a Hamming distance <= *D*, and that, among those, the "true" sequence is the one occurring most often in the dataset.  
The method works as follows. First, the sequence graph is built, where nodes are sequences and edges connect sequences at distance <= *D*. Then, a greedy procedure identifies stars in the graph, as follows. Starting from the highest abundant sequence *S*, it creates a group consisting of *S* plus all its neighbours (sequences at distance <= *D*) whose abundance is at most *F* times the abundance of *S*. The "true" sequence B is the centre of the star (in this case, *S*). The procedure stops when >= 20% of the neighbours of the current sequence *S* are also neighbours of a previously considered sequence. 
For more details on the method, please refer to our [preprint](https://doi.org/10.1101/2023.06.28.546923). 

Below we assume that the sequences and their abundance have already been computed (using *e.g.* [seqkit](https://bioinf.shenwei.me/seqkit/)). See also a sample input file in ```example/in/```.

## Requirements

* a C++ compiler with C++11 support
* R (for plotting)

## Download and installation

Either clone the repository:

```
$ git clone https://github.com/fnadalin/barcode_groups.git
```

or download and decompress the tarball archive:

```
$ wget https://github.com/fnadalin/barcode_groups/-/archive/master/barcode_groups-master.tar.gz
$ tar -xzf barcode_groups-master.tar.gz barcode_groups
```

To compile the code:

```
$ cd barcode_groups
$ make
```

## Instructions

Folder ```example/``` contains a sample input and output.  
First define the variables (here, DIST is *D* and FRAC is *F*; *F* = 1 means that we do not put any constraint on the abundance of the neighbours):

```
$ DIST="1"
$ FRAC="1"
$ input="example/in/in_with_counts.GBC.txt"
$ outdir="example/out"
```

then run the example instance:

```
$ ./barcode-groups ${input} ${DIST} ${FRAC} ${outdir} > ${outdir}/STDOUT 2> ${outdir}/STDERR
```

and generate plots from the output files:

```
$ MIN_GROUP_COUNT=$(grep "Min group count" ${outdir}/STDOUT | sed "s/Min group count: //g")
$ Rscript barcode-groups-plots.R ${outdir} ${outdir}/plots "my sample" ${MIN_GROUP_COUNT} > ${outdir}/R.STDOUT 2> ${outdir}/R.STDERR
```

## Citing

If you find this software useful, please cite:

Nadalin *et al.* Multi-omic lineage tracing predicts the transcriptional, epigenetic and genetic determinants of cancer evolution. *biorxiv*. https://doi.org/10.1101/2023.06.28.546923
