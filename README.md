# ViQUF

New algorithm for full viral haplotype reconstruction and abundance estimation. It is an alternative approach from the original methodology viaDBG:

* Building assembly graph
* Polishing the assembly graph by:
	* Classifying edges as weak and strong, and removing the weak ones called as filigree edges.
	* Removing isolated nodes
	* Removing short tips
* Paired-end association and polishing removing both:
	* Extremely low and high frequency pairs.	
* Core algorithm:
	* For every pair of adjacent nodes A -> B, we built a DAG from their paired-end information.
	* DAG is transform into a flow network and the min-cost flow problem is solved.
	* The flow is translated into paths via maximal allowed flow. 
* Final strains are build following two rules:
	* Standard contig traversion
	* Maximum flow (Edmund Karp) + Min-cost flow problem - it allows to estimate relative and absolute abundance for every contig reported.

## Command line standard:

screen -L -Logfile BCAMLTest.log python test/testBcalm.py ../Datasets/Helsinki2.0/10-strain-HCV-20000x 121 10 ngs
screen -L -Logfile execution3.log ./bin/output.out tmp/unitigs.FM placements tmp 121 tmp/unitigs.graph tmp/unitigs.unitigs.fa tmp/unitigs-viadbg.fa tmp/Ownlatest/append.fasta --debug


