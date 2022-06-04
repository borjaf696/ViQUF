# ViQUF

New algorithm for full viral haplotype reconstruction and abundance estimation. It is an alternative approach for the first developed approach viaDBG which was based entirely on **de Bruijn** graphs. ViQUF is based on flow networks which allows us to do a proper and successful estimation of the strains frequencies. 

The overall workflow is as follows:

* Building assembly graph (BCALM)
* Polishing the assembly graph by:
	* Classifying edges as weak and strong, and removing the weak ones called as filigree edges.
	* Removing isolated nodes.
	* Removing short tips.
* Paired-end association.
* Paired-end polishing removing as many wrong associations as possible.	
* Core algorithm:
	* For every pair of adjacent nodes A -> B, we built a DAG from their paired-end information.
	* DAG is translated into a flow network and a min-cost flow is solved.
	* The flow is translated into paths via a "greedy" path heuristic. 
* Final strains are build following two rules:
	* Standard contig traversion (deprecated)
	* Min-cost flow over the Approximate Paired de Bruijn Graph built from the core algorithm.

# Depedencies

* Python 3.* - we encourage you to build a conda environment and still all dependencies via conda: conda create -n ViQUF-env python=3.6
	* Biopython, altair, gurobi 
	* matplotlib, scipy, numpy
* C++17
* SDSL
* BCALM
* quast 4.3 or quast 5.0.1 to evaluate the results
* gatb-library:
	* cd lib/ && rm -r gatb-core
	* git clone https://github.com/GATB/gatb-core.git
	* follow the instructions in: https://github.com/GATB/gatb-core
* lemon 1.3.1:
	* cd ..
	* wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz
	* tar xvf lemon-1.3.1.tar.gz
	* follow instructions in https://lemon.cs.elte.hu/trac/lemon/wiki/InstallLinux (IMPORTANT: install lemon in /usr/local otherwise the Makefile will not get the lemon package)
* Compile: make clean && make


## Updates (20/10/2021)

* Amplicons - new approach to deal with amplicons (On progress).
* Third Generation Sequencing - new approaches to deal with this type of data. Right now, we are able to infer a valid flow from the data (all test over simlord high depth simulated data).
	* TODO: Apply Flow decomposition with subpath constraints.

## Dockerfile
Exists a Dockerfile which automatizes the installation procedure. To use it just run:

* sudo docker rm [your_decision_docker_name]
* sudo docker build -t [your_decision_docker_name] . --no-cache
* sudo docker run -d --name [your_decision_docker_name] [your_decision_docker_name]

## Command line standard:

The file **execution-script** contains an example about how to execute the code.

* python scripts/testBcalm.py $1 $2 10 ngs $3 $4 --no-meta
	* $1 - folder with NGS reads
	* $2 - kmer size
	* $3 - --correct/--no-correct to perform correction or not respectively
	* $4 - --joined, has pear been executed? If so --joined otherwise --no-join
* ./bin/output.out tmp $2 tmp/unitigs.graph tmp/unitigs.unitigs.fa tmp/unitigs-viadbg.fa tmp/Ownlatest/append.fasta $5 $6 --virus
	* $5 - complete set of reads (it is not mandatory but recommended)
	* $6 - --debug or not.
* python scripts/post-process.py
	* Linear programming algorithm to adjust contigs frequencies, there is not mandatory but suggested.


