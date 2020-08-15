# viaDBG-2.0

New methodology for viral quasispecies inference:
	
	* From de Bruijn Graphs to Assembly graphs
	* From naïve polishing to more sofisticated methodology (p.e. filigree edges)
	* Paired-end associated bottle neck corrected
	* From cliques graph to network flows
	* New debug mode based on logging the whole process

## Command line standard:

screen -L -Logfile BCAMLTest.log python test/testBcalm.py ../Datasets/Helsinki2.0/10-strain-HCV-20000x 121 10 ngs
screen -L -Logfile execution3.log ./bin/output.out tmp/unitigs.FM placements tmp 121 tmp/unitigs.graph tmp/unitigs.unitigs.fa tmp/unitigs-viadbg.fa tmp/Ownlatest/append.fasta --debug


