#!/usr/bin/perl
use strict;
use warnings;

if( @ARGV < 1) {
   print "collapseCDSGFF <GeneGffFile> [item]\n";
   exit();
}

open (GFF, "grep CDS " . $ARGV[0] . "|") || die "can't open file $ARGV[1]";

my $curgene = "";
my $start = -1;
my $end = -1;

my $cdsname = "CDS";
if(@ARGV>1) {
	$cdsname = $ARGV[1];
}

my @oldgenefields;
while(<GFF>)
{
  if(substr($_, 1,1) ne "#") {
  	chomp($_);
  	my @genefields = split(/\t/, $_);
  	my %geneinfo = split /[;=]/, $genefields[8];
  	my $genename = "";
				
	if(exists $geneinfo{"Parent"}) {
		$genename = $geneinfo{"Parent"};
	} elsif(exists $geneinfo{"Name"}) {
		$genename = $geneinfo{"Name"};
	} elsif(exists $geneinfo{"ID"}) {
		$genename = $geneinfo{"ID"};
	}
	
	### If its J.sinuosa then something different
	if($genename =~ /Jalsin_scf/ ) { 
		$genename =~ s/\.[0-9]$//g;
	}
	if($genename =~ /Solyd/) {
		$genename =~ s/\.[0-9]$//g;
	}
	if($genename =~ /Medtr/) {
		$genename =~ s/\..*//;
	}
	# get rid of weird embalishments
	
	$genename =~ s/mRNA\.|CDS//g;
	$genename =~ s/mRNA:|CDS://g;
	$genename =~ s/\.v[0-9].*//;
	$genename =~ s/\.[0-9]\.p$//;
	$genename =~ s/\.[0-9]$//;
	$genename =~ s/\.0[0-9]$//;
				
	if($genefields[2] eq $cdsname) {
		if($genename eq $curgene || $curgene eq "") {
		
			if($curgene eq "") {
				$curgene = $genename;
				@oldgenefields = @genefields;
			}
			if($start == -1 || $genefields[3] < $start) {
				$start = $genefields[3];
			}
			if($end == -1 || $genefields[4] > $end) {
				$end = $genefields[4];
			}
		} else {
			print $oldgenefields[0] . "\t" . $oldgenefields[1] . "\t" . "CDSg" . "\t" . $start . "\t" . $end . "\t" . $oldgenefields[5] . "\t" . $oldgenefields[6] . "\t" . $oldgenefields[7] . "\tName=" . $curgene . "\n";
			@oldgenefields = @genefields;
			$curgene = $genename;
			$start = $genefields[3];
			$end = $genefields[4];
		}
    }
  }
}
