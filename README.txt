*************************************************************************************************************************
Conservatory v2.0.1

Last update: 4/12/2023
*************************************************************************************************************************

  Conservatory is a structural-variants and polyploidization-aware program for the identification of conserved non-coding
	sequences (CNS) in a large number of genomes. Conservatory identifies cis-elements and protein sequence microsynteny to derive
	multiple alignments of upstream and downstream regulatory sequences in closely related species. It then reconstructs
	the ancestral sequences for each CNS and searches for that in syntenic regions of distant genomes. 
	Conservatory utilizes multiple genomes per species to account for genome-specific variants or annotation problems.
	Conservatory uses a reference genome to define the CNS. After running conservation for each gene separately,
	 conservatory merges all the CNS for a given reference genome.

	Conservatory can also utilize multiple reference genomes. A separate script than merges the CNS from all references
         to generate a single concensus file of CNSs.

	Conservatory was developed as a collaboration between the Efroni (HUJI) and Lippman (CSHL) labs.

	This is version 2.0.1.

	Previous version (1.5) was described in this publication:

	Hendelman A., Zebell S., Rodriguez-Leal D., Dukler N., Robitaille G., Wu X., Kostyun J., Tal L., Wang P.,
		Bartlett M.E., Eshed Y., Efroni I*, Lippman ZB.* (2021) Conserved pleiotropy of an ancient plant
			homeobox gene uncovered by cisregulatory dissection. Cell, 184:1724-1739.


	For questions\comments email Idan Efroni @ idan.efroni@mail.huji.ac.il \ Zachary Lippman lippman@cshl.edu
	
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

3. Download the genomes to the genomes/family directory. For each genome, conservatory requires three files: the genome fasta file (bgziped),
   an annotation GFF3 file and the protein sequence file. Files should be named. <genomename>.fasta.gz, <genomename>.genes.gff3 and <genomename>.proteins.fasta, respectively.
   
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

5. For each family, generate a phylogenetic tree in Newick format and store it as "<family>.tree" in the genomes directory.

6. Build the genome database for a reference genome by running from the conservatory directory:

	./processGenomes --reference <referenceGenomeName> --threads <number of threads>
	
	This will process all the genomes in the genome_database.csv file for a specific reference. It will create orthologs files and regulatory sequence databases.
	 processGenomes will not overwrite existing genome databases and indexes. To do so, use the --force or --force-orthology flag. For more info, see processGenomes --help
		 
	*** This can be a very long process but only need to be run once. On a single node 16-core Xeon(R) CPU E5-2630 v3 @ 2.40GHz, each genome requires about an hour of computer time.

	If you want to use multiple reference genomes, this script has to be run for each one separately

7. Generate a neutral model for phyloP under the name <family>.mod and place in the genomes directory.
   Conservatory can try a generate a neutral model automatically based on 4-fold codon variants. For that, you need the CDS sequences for each genome named <genomename>.cds.fasta in the family directory
   
   * First, determine the orthology by CNS alignments using:
	./processConservation --reference <referenceGenomeName> --just-family-alignment
	* this computes the alignments for all genes in the genome. Depending on genome size and computer resources, this can take hours to days.
   * Build the model using
	./script/buildModel --reference <referenceGenomeName> --tree <treeFileName>

   The model will be placed in the genome directory.

   The models in the github directory are those of the 10 reference genome families currently used by us.
   
7. Once the genome database has been generated, you can build the CNS for any given gene using:
	./scripts/buildConservation --reference <referenceGenome> --locus <geneName> --min-phylop-score [score]

	The min-phylop-score is the minimum -log-pvalue used to define conservation. This can vary greatly between different reference genomes and families.
	  this parameter has to be determined empirically and we detail a method for that in the paper. The git hub cutoffs.csv file has the scores we used for the 10 reference genomes

	Output alignment files will be in alignments\<familyName>\<geneName>
	Output CNS files will be in CNS\<familyName>\<geneName>

	rocess all genes in the genome. This is a long process.
	
	See --help for additional options.
	

Conservatory can be run for individual genes. Data can be merged to generate a complete list for a given reference genome or merged between different reference genome
  to generate a composite CNS database.


**********************************************************************************************************************************************************

** Update History

3/12/2022: Version 2.0.1 beta version

