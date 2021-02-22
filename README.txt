*********************************************************************************************************************************************************************
Conservatory v1.0

Last update: 2/22/2021
*********************************************************************************************************************************************************************

  Conservatory is structural variants and polyploidization-aware program for identification of conserved non-coding sequences.
	conservatory is based on identification of cis-elements and protein sequence microsynteny to assembled multiple alignments
	of upstream and downstream regulatory sequences. Conservatory then relies on phyloP to derive a model-based conservation
	scores and assembled short regions of conservation.
	Conservatory is also able to utilize multiple genome per species as pan genomes to account for genome-specific variants or
	annotation problems.
	
	The program is described in:
	Hendelman A., Zebell S., Rodriguez-Leal D., Dukler N., Robitaille G., Wu X., Kostyun J., Tal L., Wang P., Bartlett M.E., Eshed Y., Efroni I*, Lippman ZB. (2021) 
		Conserved pleiotropy of an ancient plant homeobox gene uncovered by cisregulatory dissection. Cell () 

	For questions\comments email Idan Efroni @ idan.efroni@mail.huji.ac.il
	
*********************************************************************************************************************************************************************

Setup

*** We use bioconda for all installation and environment setup. 

1. Download the program files from github (https://github.com/idanefroni/Conservatory)

2. Make sure the following programs are install and in the conda environment
	* lastz
	* multiz
	* phyloP
	* samtools
	
3. Install vista (http://genome.lbl.gov/vista/downloads.shtml) and MafFilter (https://jydu.github.io/maffilter/) to the program subdirectory.
    if these programs are installed at a different path, update the conservatory.parameters.txt file

4. Update regulatory sequence length paramters in conservatory.parameters.txt file (default is 50K upstream and 5K downstream)
    the program can easily process longer sequences, but it will require additional diskspace.
  
5. Download the genomes to the genome/family directory. For each genome, conservatory requires three files: the genome fasta file (bziped),
   an annotation GFF3 file and the protein sequence file. Files should be named. <genomename>.fasta.gz, <genomename>.genes.gff3 and <genomename>.proteins.fasta, respectively.
   
6. Once all genomes are downloaded, set up the genome_database.csv file. This is a comma delimited text file with the following fields:
 Genome - the genome name (must match the name of the genome data files)
 Species - the names of the species. Multiple genomes can be assigned to one species. Conservatory will select the best genome for the species on a gene-by-gene basis
 Family - the name of the species family. A family (which does not has to correspond to a taxonimic family) is a group of genomes which can be compared to one reference
 	      genome. Conservation is also computed within a specific family, although multiple family can exist in the same conservatory directory.
 ReferenceGenome - the genome name of the reference genome for the specific family (genome must exists in the file).
 MinQuality - minimum alignment quality to consider (default:1000)
 GeneNameField - the field in the GFF3 file CDS lines containing the gene name
 GeneProcessingRegEx - A regular expression that will be deleted from the gene name. This is used to match gene and protein names.
 Gene2SpeciesIdentifiedString - A regular expression the should appear at the begining of the gene name and which will uniqely identify the genome
 ProteinProcessingRegEx - optional. A regular expression that will be deleted from the protein name. This is used to match gene and protein names. 

The genome specification for all genomes used in Hendelman et al. are already in the file.

7. Build the genome database and index by running from the conservatory directory:

	./processGenomes --threads <thread-num> --verbose
	
	This will process all the genomes in the genome_database.csv file. The processing can be limited to specific families or genomes. processGenomes will not overwrite existing
	 genome databases and indexes. To do so, use the --force flag. For more info, see processGenomes --help
	 
	This process can take several hours per genome, depending on the number of CPU cores and genome size.
	
8. Build the phyloP model. For this you will require a phylogenetic tree of each family and an alignment of the an alignment of neutrally evolving sequence. For the paper,
    we used the fourfold codon position for a large number of orthologous genes. Other options exists.
    
	phyloFit --tree <family>.tre --nrates 4 -o genomes/<family> --subst-mod HKY85+Gap <neutralAlignment>.fasta

8. Once the genome database has been generated, you can build the CNS for any given gene using:
	./scripts/buildConservationFormFamily1 --family <familyName> --locus <genename>
	
	or generate CNS for all genes using the multithreaded script:
	./processConservation --family <family> --version 1
	
	Calculating conservation for entire families takes several hours, depending on number of CPU cores and genome size.
	
	See --help for additional options.
	
At the end, the alignment directory contains the multiple alignments for the cis elements of each gene (as MAF files) and the conserved non-coding regions as a BED file and a CSV table in the CNS directory.

*********************************************************************************************************************************************************************

Update 2/22/2021:

The ortholog selection and alignment algorithm has been modified from the one used in the paper. To use the modified algorithm, use '--version 2' in processConservation.

