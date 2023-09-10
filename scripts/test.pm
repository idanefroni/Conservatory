#!/usr/bin/perl

use POSIX;
use strict;

use Cwd;
use Cwd 'abs_path';

use lib './scripts';
use ConservatoryUtils;
use GenomeDatabase;
use MappingDatabase;
use CNSDatabase;
use ConservatoryTree;
use CREOrthologList;


$|=1;

###############################################################################
######### Setup

my $conservatoryDir=abs_path(".");


my $alignmentSeqOverlapToMerge=0;


my $CNSDB = new CNSDatabase();
my $mappingDB = new MappingDatabase();
my $genomeDB = new GenomeDatabase($conservatoryDir);
my $conservatoryTree = new ConservatoryTree($genomeDB);
my $tseq1 =   "CTCAATTAGTTATAATTCTCATGCATTATTGAGTCT-ATGATTTTCAACATTGAAAACTATTTGTACACAACCATAACGAACTTGTTTTTGCCAATTAGAATCATATTTCAACTTATAAAAACGATAGATGGACAATTTC-ATATGTAACAAGTATTTGGACCACCACGTTAGATAGTCAAAAATCACCCTACAACCAAGGACATCAATCAT";
my $refSeq1 = "CTTAGCTAGTGATAATTCTTATGCATTATTGAGACTTGTGATTTTCAACATTGAAAAGTATTTATACACAACCATAACGAACCTGTTTTTGCCAATTAGAATCATATTTCAACTTATAAAAACGATAGATGGATAATTTTGATATGTAACAAGTA-TTGGACCACCACGTTAGATAGTCAGAAAT--------AACCAAGGACATCAATCAT";


my $mapping = new Mapping("CNSID", "Aalpina", "Aa_G754810", 2661, "-", 2612, "ATCGTTTAGACCATAGTACATTTTAGCCTTCATAGAGATGCTAAAAAAGAAATTGTTAGTTGACATTTAATTACAAC--TTACAAGCATAT--GTACAG-AATTTATAGAGATGTAAG",
"ATCGTTAAGACCATATGCCTTTTTTGCCTTCATAGAGATGCT-AAAAAGAAATGGTTAGTTGACA--TAGTTACAACCATTA-AATCATATAAGAATATCGATTCATAGAGATGTGGG");



$mapping->print;

my $start = 2631;
my $end = 2631+68;
my $subMapping = $mapping->createSubSetMapping($start, $end);
$subMapping->print;


2661

2697

2689
#$mappingDB = new MappingDatabase("./CNS/Athaliana/AT2G17950.map.csv");
#$mappingDB->orderMappingByGenome();

#my $camMappings = new MappingDatabase();
#foreach my $curMapping (@{ $mappingDB->getMappingsByOrder() }) {
#	if($curMapping->getGenome() eq "Ptrichocarpa" $curMapping->getGenome() eq "Athaliana") {
#		$curMapping->fillAbsoluteCoordiantes($genomeDB);
#		$curMapping->print;
#		$camMappings->add($curMapping);
#	}
#}

#$camMappings->writeDatabase("tmp/z.map.csv");


