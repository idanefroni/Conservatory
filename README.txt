*************************************************************************************************************************
Conservatory v2.0.1

Last update: 4/9/2025
*************************************************************************************************************************

  Conservatory is a structural-variants and polyploidization-aware program for the identification of conserved non-coding
	sequences (CNS) in a large number of genomes. Conservatory identifies cis-elements and protein sequence microsynteny to derive
	multiple alignments of upstream and downstream regulatory sequences in closely related species. It then reconstructs
	the ancestral sequences for each CNS and searches for them in syntenic regions of distant genomes. 
	Conservatory utilizes multiple genomes per species to account for genome-specific variants or annotation problems.
	Conservatory uses a reference genome to define the CNS. After running conservation for each gene separately,
	 conservatory merges all the CNS for a given reference genome.

	Conservatory can also utilize multiple reference genomes. A separate script that merges the CNS from all references
         to generate a single consensus file of CNSs.

	Conservatory was developed as a collaboration between the Efroni (HUJI) and Lippman (CSHL) labs.

	This is version 2.0.1.

	Amundson K.R., Hendelman H., Ciren D., Yang H., de Neve1 A.E., Tal S., Sulema A., Jackson D., Bartlett M.E.,
		Lippman Z.B., Efroni I. (2025) A deep-time landscape of plant cis-regulatory sequence evolution.
		bioXriv doi:

	The processed CNS for 314 genomes detailed in the paper is available from www.conservatorycns.com.


	The previous version (1.5) was described in this publication:

	Hendelman A., Zebell S., Rodriguez-Leal D., Dukler N., Robitaille G., Wu X., Kostyun J., Tal L., Wang P.,
		Bartlett M.E., Eshed Y., Efroni I*, Lippman ZB.* (2021) Conserved pleiotropy of an ancient plant
			homeobox gene uncovered by cisregulatory dissection. Cell, 184:1724-1739.

	For questions\comments email Idan Efroni @ idan.efroni@mail.huji.ac.il \ Zachary Lippman lippman @ cshl.edu
	
**********************************************************************************************************************************************************

Setup V2.0.1

*** We use bioconda for all installation and environment setup. 

1. Download the program files from GitHub (https://github.com/idanefroni/Conservatory)

2. Make sure the following programs are installed and in the conda environment
	* lastz
	* mafft
	* phast
	* samtools
	* bedtools
	* wig2bed
	* bioperl
	* blastp
	* agat
	* seqkit
	* bgzip

	You can also build the conda environment using the ./buildCondaEnvironment script in the Conservatory directory.

3. Download the genomes to the genomes/family directory. For each genome, conservatory requires three files: the genome fasta file (bgziped),
   an annotation GFF3 file, and the protein sequence file. Files should be named. <genomename>.fasta.gz, <genomename>.genes.gff3 and <genomename>.proteins.fasta, respectively.
   
4. Once all genomes are downloaded, set up the genome_database.csv file. This is a comma-delimited text file with the following fields:
 Genome - a unique genome name (must match the name of the genome data files)
 Species - the names of the species. Multiple genomes can be assigned to one species. Conservatory will select the best genome for the species on a gene-by-gene basis
 Family - the name of the species family. A family (which does not have to correspond to a taxonomic family) is a group of genomes that can be compared to one reference
 	      genome by sequence alignment. For plants, this means usually less than 60 million years of divergence.
 ReferenceGenome - <Depreceated>
 UpstreamLength - length of gene upstream region to search for the specific genome. Leave this blank, and the processGenome script will fill this field automatically based on genome characteristics.
 DownstreamLength - length of gene downstream region to search for the specific genome. Leave this blank, and the processGenome script will fill this field automatically based on genome characteristics.
 GeneNameField - the field in the GFF3 file CDS line containing the gene name
 GeneProcessingRegEx - A regular expression that specifies what characters to delete from the gene name to generate a simple gene name which is used to match gene and protein names.
 Gene2SpeciesIdentifiedString - <Deprecated>
 ProteinProcessingRegEx - A regular expression that specifies what characters to delete from the protein name to generate a simple protein name. This is used to match gene and protein names.
 Classification - Six-level classification, separated by dashes, that specify the phylogenic classification of the genome. e.g. Embryophyta-Tracheophytes-Angiosperms-Dicots-Asterids-Solanaceae

 The genome_database.csv file provided in the github directory has the 314 genomes currently used.

5. For each family, generate a phylogenetic tree of all species in the family in Newick format and store it as "<family>.tree" in the genomes directory. 
   This is best done using an alignment of known ortholog protein sequences, but can also be constructed using prior knowledge. The tree does not need to be calibrated - this will be done by Conservatory. 
	The directory currently has the trees for all families used in the paper.

6. Build the genome database for a reference genome by running from the conservatory directory:

	./processGenomes --reference <referenceGenomeName> --threads <number of threads>
	
	This will process all the genomes in the genome_database.csv file for a specific reference. It will create orthologs files and regulatory sequence databases.
	 processGenomes will not overwrite existing genome databases and indexes. To do so, use the --force or --force-orthology flag. For more info, see processGenomes --help
		 
	*** This can be a very long process but only need to be run once. On a single node 16-core Xeon(R) CPU E5-2630 v3 @ 2.40GHz, each genome requires about an hour of computer time.

	If you want to use multiple reference genomes, this script has to be run for each one separately

7. Generate a neutral model for phyloP under the name <family>.mod and place it in the genomes directory.
   Conservatory can try a generate a neutral model automatically based on 4-fold codon variants. For that, you need the CDS sequences for each genome named <genomename>.cds.fasta in the family directory.
   If such a file does not exist, processGenome will attempt to generate it automatically during setup.
   
   * First, determine the orthology by CNS alignments by running 'buildConservation --just-family-alignment' for all genes in the genome.
	This can be done by looping over all genes in the reference and running
	./scripts/buildConservation --reference <referenceGenomeName> --just-family-alignment --locus <geneID> [--tmp-dir $TMPDIR]

	This computes the alignments for all genes in the genome. Depending on genome size and computer resources, this can take hours to days.
	This part is best performed in parallel. We use the processConservation job located in the jobs directory.

   * After alignments have been computed, you can build the model and calibrate the tree using
	./script/buildModel --reference <referenceGenomeName> --tree <treeFileName>

   The model will be placed in the genomes directory. The models in the github directory are those of the 10 reference genome families currently used by us.
   
7. Once the genome database has been generated, you can build the CNS for any given gene using:
	./scripts/buildConservation --reference <referenceGenome> --locus <geneName> --min-phylop-score [score]

	The min-phylop-score is the minimum -log-pvalue used to define conservation. This can vary greatly between different reference genomes and families, and depends on the phylogenetic coverage in your family
		and the required sensitivity.
	
	  For our dataset, the default was 1.75, but if multiple references are used, this parameter is best determined empirically. The cutoffs.csv file has the scores we used for the 10 reference genomes

	You can run this script by iterating over all genes in the references, or in parallel on a computing cluster.

	Output alignment files will be in alignments/<referenceName>/<geneName>
	Output CNS files will be in CNS/<referenceName>/<geneName>

	There are two comma-delimited files generated for each gene: a CNS file, containing the list of unique CNSs found, and a MAP file, listing all the instances of the CNS in all genomes in the family.

	The fields in the CNS file are:
		ReferenceGenomeName
		UniqueCNSID
		Locus
		Relative position in the reference genome
		CNS length
		CNS Age (deepest node in the tree)
		Number of species having the CNS

	The fields in the MAP file are:
		UniqueCNSID
		Species Name
		Locus
		Relative position to the gene
		Strand (+ or -) relative to the reference CNS
		Start of alignment position in the reference CNS
		CNS sequence in the current genome
		CNS sequence in reference genome
		Chromosome
		Absolute position
		Chromosome strand
		Unique location name (optional)
	
8. Conservatory can be run for individual genes, but is more powerful when the CNS data is merged into a unified dataset. To merge all CNS for a reference genome:

	Concat all CNS/<referenceName>/*.cns.csv directory to a single file <referenceName>.cns.csv
	Concat all CNS/<referenceName>/*.map.csv directory to a single file <referenceName>.map.csv

	then run
	./scripts/mergeCNS --in-cns <referenceName>.cns.csv --in-map <referenceName>.map.csv --out-cns <referenceName>.merged.cns.csv --out-map <referenceName>.merged.map.csv

	This will generate a single merged file with unique names for all CNS in the reference genome.

9. If you want to merge CNS dataset from multiple references:

	Concat all the <referenceName>.merged.cns.csv files to conservatory.cns.csv
	Concat all the <referenceName>.merged.map.csv files to conservatory.map.csv

	Then merge using:

	./scripts/mergeCNS --merge-deep --sort --in-cns conservatory.cns.csv --in-map conservatory.map.csv --out-cns conservatory.merged.cns.csv --out-map conservatory.merged.map.csv

	Running mergeCNS on a cluster:
	  * mergeCNS can be run as a single process, but it can take a long time for a large number of CNSs. Running time can be improved by running a helper parallel alignment process
		to do that user add the following: --align-CNS-process --tmp-dir  [temporary directory]
        execute the script. It will generate all the files requiring alignments and will let you know it is waiting for the aligner. After all files were generated,
		split the conservatory.cns.csv file into N separate files (N being the number of parallel jobs you plan to run).

		Then run in parallel:
		./scripts/reconstructCNSSequences --in-cns <PartialCNSFile> --in-map conservatory.map.csv --out-cns <outputPartialCNSFile> --out-map <outputPartialMapFile> --break --reconstruct

		Once all files were aligned by the help, mergeCNS will continue to process the data.
	* Note: mergeCNS loads all CNS data to memory. Merging a family can be done on a 16GB RAM machine. Merging the 314 genome dataset required 192GB

10. After merging, CNS data needs to be filtered. CNS can be renamed so they have more meaningful UniqID. To do that run
		./scripts/CNSUtils --in-cns conservatory.merged.cns.csv --in-map conservatory.merged.map.csv --out-cns conservatory.final.cns.csv --out-map conservatory.final.map.csv --filter-long --filter-deep-cns --filter-repetitive --filter-plastid --rename-cns
							--rename-positions --verbose --out-reject conservatoryV10.reject.cns.csv

   It is recommended to verify the integrity of the data after the merge procedure using:
		./scripts/CNSUtils --in-cns conservatory.final.cns.csv --in-map conservatory.final.map.csv --fix-positions --sort --out-cns conservatory.final.verified.cns.csv --out-map conservatory.final.verified.map.csv --out-reject conservatory.final.rejected.map.csv

11. Finally, you can convert the raw CNS and MAP files to GFF3 files that can be used as track on genome browser. To do so, run:

	./script/CNSUtils --in-cns conservatory.final.cns.csv --out-cns conservatory.final.map.csv --gff-reference <genomeToGenerateGFFTo> --make-gff <outputGFFFileName>

	
**********************************************************************************************************************************************************

** Update History

4/9/2025: Release 2.0.1
3/12/2022: Version 2.0.1 beta version

