#!/usr/bin/perl

package CREOrtholog;

use POSIX;
use strict;
use warnings;
use List::Util qw(min max);
use lib './scripts';
use ConservatoryUtils;

sub new {
    my $class = shift;
    my ($genomeDB, $genome, $locus, $refLocus,$refUpFastaFileName, $refDownFastaFileName, $laxAlignment) = @_;

    if(!defined $laxAlignment) { $laxAlignment =0; }
    my $tmpDir = $genomeDB->getTemporaryDir();

    mkdir $tmpDir;
	my $tmpUpstreamFastaName ="$tmpDir/$genome.$locus.up.fasta";
	my $tmpDownstreamFastaName = "$tmpDir/$genome.$locus.down.fasta";

	extractFastaToFile($genomeDB->getUpstreamFastaFileName($genome), $locus, $tmpUpstreamFastaName);
    if($refDownFastaFileName ne "") {
	    extractFastaToFile($genomeDB->getDownstreamFastaFileName($genome), $locus, $tmpDownstreamFastaName);
    }

	my $tmpMAFUpFileName = "$tmpDir/$locus.up.maf";
	my $tmpMAFDownFileName = "$tmpDir/$locus.down.maf";

	system(generateLastzCommandLine($genomeDB,$laxAlignment) . " $refUpFastaFileName" . "[multiple] $tmpUpstreamFastaName > $tmpMAFUpFileName");
    unlink($tmpUpstreamFastaName);

    if($refDownFastaFileName ne "") {
	    system(generateLastzCommandLine($genomeDB,$laxAlignment) . " $refDownFastaFileName" . "[multiple] $tmpDownstreamFastaName > $tmpMAFDownFileName");    
        unlink($tmpDownstreamFastaName);        
    }


	my $upQuality = getMAFQuality($tmpMAFUpFileName);
    my $downQuality=0;
    if($refDownFastaFileName ne "") {
	    $downQuality = getMAFQuality($tmpMAFDownFileName);
    }
  
    my $self = {
        _GenomeDB => $genomeDB,
	    _Locus => $locus,
        _Genome => $genome,
        _tmpUpstreamFastaName => $tmpUpstreamFastaName,
        _tmpDownstreamFastaName => $tmpDownstreamFastaName,
        _tmpMAFUpFileName => $tmpMAFUpFileName,
        _tmpMAFDownFileName => $tmpMAFDownFileName,
        _downQuality => $downQuality,
        _upQuality => $upQuality
    };

    bless $self, $class;
    return $self;
}

sub removeTemporaryFiles {
    my ($self) = @_;
#    unlink($self->{_tmpMAFUpFileName});
    unlink($self->{_tmpMAFDownFileName});
}

sub getQuality {
    my ($self) = @_;
    return $self->{_upQuality} + $self->{_downQuality};
}
sub getLocus() {
    my ($self) = @_;
    return $self->{_Locus};  
}

sub getGenome() {
    my ($self) = @_;
    return $self->{_Genome};  
}

sub getMAFFileName {
    my ($self, $upOrDown) = @_;
    if($upOrDown eq "U") {
        return $self->{_tmpMAFUpFileName};
    } elsif($upOrDown eq "D") {
        return $self->{_tmpMAFDownFileName};
    } else {
        die "ERROR: getMAFFileName needs a U or D.\n";
    }
}
sub DESTROY {
    my $self = shift;
    $self->removeTemporaryFiles();
}

sub print {
    my ($self, $outFile) = @_;
    if(!defined $outFile) { $outFile = \*STDOUT; } 
    print $outFile join(",",$self->{_Genome}, $self->{_Locus}, $self->getQuality()) . "\n";
}

1;