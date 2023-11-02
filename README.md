# ViQUF

A new algorithm for full viral haplotype reconstruction and abundance estimation. It is an alternative approach for the first developed approach viaDBG which was based entirely on **de Bruijn** graphs. ViQUF is based on flow networks which allows us to do a proper and successful estimation of the strain frequencies. 

The overall workflow is as follows:

* Building assembly graph (BCALM)
* Polishing the assembly graph by:
	* Classifying edges as weak and strong, and removing the weak ones called filigree edges.
	* Removing isolated nodes.
	* Removing short tips.
* Paired-end association.
* Paired-end polishing removing as many wrong associations as possible.	
* Core algorithm:
	* For every pair of adjacent nodes A -> B, we built a DAG from their paired-end information.
	* DAG is translated into a flow network and a min-cost flow is solved.
	* The flow is translated into paths via a "greedy" path heuristic. 
* Final strains are built following two rules:
	* Standard contig traversing (deprecated)
	* Min-cost flow over the Approximate Paired de Bruijn Graph built from the core algorithm.

# Depedencies

* Python 3.* - we encourage you to build a conda environment and still all dependencies via conda: conda create -n ViQUF-env python=3.6
	* Biopython, altair, gurobi 
	* matplotlib, scipy, numpy
* C++17
* SDSL
* BCALM
* quast 4.3 or quast 5.0.1 to evaluate the results
* gatb-library: follow the instructions in: https://github.com/GATB/gatb-core
```bash
cd lib/ && rm -r gatb-core
git clone https://github.com/GATB/gatb-core.git
```
* lemon 1.3.1: follow instructions in https://lemon.cs.elte.hu/trac/lemon/wiki/InstallLinux (IMPORTANT: install lemon in `/usr/local` otherwise the Makefile will not get the lemon package)
```bash
cd ..
wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz
tar xvf lemon-1.3.1.tar.gz
```
* Compile:
```bash
make clean && make
```


## Updates (20/10/2021)

* Amplicons - a new approach to deal with amplicons (In progress).
* Third Generation Sequencing - new approaches to deal with this type of data. Right now, we are able to infer a valid flow from the data (all test over simlord high depth simulated data).
	* TODO: Apply Flow decomposition with subpath constraints.

## Dockerfile
Exists a Dockerfile which automatizes the installation procedure. To use it just run:
```bash
sudo docker rm [your_decision_docker_name]
sudo docker build -t [your_decision_docker_name] . --no-cache
sudo docker run -d --name [your_decision_docker_name] [your_decision_docker_name]
```

## Command line standard:

The file **execution-script** contains an example of how to execute the code.
```bash
python scripts/testBcalm.py [you folder name] [kmer size] ngs [--correct/--no-correct] [--join] --no-meta
./bin/output.out tmp [kmer size] tmp/unitigs.graph tmp/unitigs.unitigs.fa tmp/unitigs-viadbg.fa tmp/Ownlatest/append.fasta [complete set of reads (optional)] [--debug] --virus
python scripts/post-process.py
```

The last step runs a "linear programming algorithm" to adjust contigs frequencies, it is not mandatory but suggested.


