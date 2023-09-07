#!/usr/bin/perl

package CREOrthologList;

use POSIX;
use strict;
use warnings;
use List::Util qw(min max);
use lib './scripts';
use ConservatoryUtils;
use CREOrtholog;

sub new {
    my $class = shift;
    my ($genomeDB, $refGenome, $refLocus, $targetGenome, $refUpFastaFileName, $refDownFastaFileName, $laxAlignment, $verbose) = @_;
    my @orthologList;
    my %species;
    if(!defined $laxAlignment) { $laxAlignment=0; }

    if(!defined $verbose) { $verbose=0; }
    if(!defined $genomeDB) { die "ERROR: CREOrthologList constructor must get a GenomeDatabase object.\n"; }

    if(defined $genomeDB && defined $refGenome && defined $refLocus) {
        my $orthologsFileName = $genomeDB->getConservatoryDir() . "/genomes/" . $genomeDB->genomeToFamily($targetGenome) . "/$refGenome.$targetGenome.orthologs.csv";

        if(! -e $orthologsFileName) { return undef; }
        open(my $orthologsFile, $orthologsFileName) || return undef;
        while(my $orthologLine = <$orthologsFile>) {
            chomp $orthologLine;
            my ($curPutativeOrtholog, $curRefLocus, $score) = split /,/, $orthologLine;
            if($curRefLocus eq $refLocus && $curRefLocus ne $curPutativeOrtholog) {
                my $putativeOrtholog = new CREOrtholog($genomeDB, geneToGenome($curPutativeOrtholog), $curPutativeOrtholog, $refLocus,$refUpFastaFileName, $refDownFastaFileName, $laxAlignment);
                if($putativeOrtholog->getQuality > $genomeDB->getMinOrthologQuality() ) {
                    if($verbose) {
                        print "PROGRESS: Found possible CRE-Ortholog " . $putativeOrtholog->getLocus() . " (Quality: " . $putativeOrtholog->getQuality() . ").\n";
                    }
                    push(@orthologList, $putativeOrtholog);
                    $species{ $genomeDB->genomeToSpecies($targetGenome) } = 1;
                } else {
                    $putativeOrtholog->removeTemporaryFiles();
                }
            }
        }
        close($orthologsFile);
    }

    my $self = {
        _GenomeDB=> $genomeDB,
	    _Orthologs => \@orthologList,
        _Species => \%species
    };

    bless $self, $class;
    return $self;
}

sub add {
    my ($self, $toAdd) = @_;

    if(ref($toAdd) eq 'CREOrtholog') {
        push (@{ $self->{_Orthologs} }, $toAdd );
        $self->{_Species}->{ $self->{_GenomeDB}->genomeToSpecies($toAdd->getGenome()) } = 1;
    } elsif(ref($toAdd) eq 'CREOrthologList') {
        my @orthListToAdd = $toAdd->getOrthologs();
        foreach my $curOrtholog (@orthListToAdd) {
            push (@{ $self->{_Orthologs} }, $curOrtholog );
            $self->{_Species}->{ $self->{_GenomeDB}->genomeToSpecies($curOrtholog->getGenome()) } = 1;
        }
    } else {
        die "ERROR: Only CREOrtholog or CREOrthologList can be added to a CREOrthologList object.\n";
    }
}

sub getOrthologs {
    my ($self) = @_;
    return @{ $self->{_Orthologs}};
}

sub getOrthologsNumber {
    my ($self) = @_;
    return scalar @{ $self->{_Orthologs}};
}
sub getSpeciesNumber {
    my ($self) = @_;
    return scalar keys %{ $self->{_Species} };
}

sub getSpecies {
    my ($self) = @_;
    return sort keys %{ $self->{_Species} };
}
#### Returns the orthologs for a particular species 
sub getOrthologsForSpecies {
    my ($self, $species) = @_;
    my @orthologsForSpecies;
    foreach my $curOrtholog (@{ $self->{_Orthologs}}) {
        if($self->{_GenomeDB}->genomeToSpecies( $curOrtholog->getGenome()) eq $species) {
            push @orthologsForSpecies, $curOrtholog;
        }
    }
    return @orthologsForSpecies;
}

#### Returns the orthologs ordered by the alignment quality (highest is first)
sub getOrthologsByQuality {
    my ($self, $species) = @_;
    my @sortedOrthologs = reverse sort { $a->getQuality <=> $b->getQuality } @{ $self->{_Orthologs} };
    return @sortedOrthologs;
}

sub removeTemporaryFiles {
    my ($self) = @_;
    foreach my $curOrtholog (@{ $self->{_Orthologs}}) {
        $curOrtholog->removeTemporaryFiles();
    }
}

#####################################################################
##### Pick the best orthologs

sub pickBestPutativeOrthologs {
    my ($self, $justOneOrtholog, $verbose) = @_;

    my @bestOrthologs;

	foreach my $curSpecies ($self->getSpecies()) {
		# first, check if we have a consistent duplication. if not, pick the best alignment.
		# For consistent duplication it means that the number of genes is a multiple of the number of species.
		# if so, pick the best duplication group

        my @orthologsForSpecies = $self->getOrthologsForSpecies($curSpecies);
        my %genomesOrthologQuality;
        my %orthologsPerGenome;

        #### Calculate the mean overall alignment quality for genome
        foreach my $curOrtholog (@orthologsForSpecies) {
            if(!defined $genomesOrthologQuality{ $curOrtholog->getGenome() } ) { $genomesOrthologQuality{ $curOrtholog->getGenome() }=0; }
            if(!defined $orthologsPerGenome{$curOrtholog->getGenome() }) { $orthologsPerGenome{$curOrtholog->getGenome() } =0 ; }

            $orthologsPerGenome {$curOrtholog->getGenome() }++;
            $genomesOrthologQuality{ $curOrtholog->getGenome() } = $genomesOrthologQuality{ $curOrtholog->getGenome() } + $curOrtholog->getQuality();      
        }
        foreach my $curGenome (keys %orthologsPerGenome) {
            $genomesOrthologQuality{ $curGenome} = $genomesOrthologQuality{ $curGenome} / $orthologsPerGenome{$curGenome};
        }

        if(scalar keys %genomesOrthologQuality >1) {  ### If we have multiple genomes for the same species.
            ## Do we have orthogroups or just one ortholog in different genomes?
            my @genomesSortedByQuality = reverse sort { $genomesOrthologQuality{$a} <=> $genomesOrthologQuality{$b} } keys %genomesOrthologQuality;
            my $alreadyPickedOne=0;

            foreach my $curOrtholog (@orthologsForSpecies) {
                if($curOrtholog->getGenome() eq $genomesSortedByQuality[0] && !($alreadyPickedOne++ && $justOneOrtholog)) {
                    push @bestOrthologs, $curOrtholog;
			        if($verbose) {
                        print "PROGRESS: Picked CRE-ortholog for $curSpecies. Genome: " . $curOrtholog->getGenome(). ". Locus: ". $curOrtholog->getLocus() . ".\n";
                    } 
                } else {
                    $curOrtholog->removeTemporaryFiles();
                }
            }
        } else {
            push @bestOrthologs, @orthologsForSpecies;
        }
    }

    $self->{_Orthologs} = \@bestOrthologs;
    my %species;
    foreach my $curOrtholog (@bestOrthologs) {
        $species{ $self->{_GenomeDB}->genomeToSpecies( $curOrtholog->getGenome()) }=1;
    }
    $self->{_Species} = \%species;
}

sub writeList {
    my ($self, $outputFileName) = @_;
    my $outputFile;

    if(defined $outputFileName) {
        open($outputFile, ">", $outputFileName);
    } else {
        $outputFile = \*STDOUT;
    }
    foreach my $curOrtholog (@{ $self->{_Orthologs}}) {
        $curOrtholog->print($outputFile);
    }
    if(defined $outputFileName) {    
        close($outputFile);
    }
}

sub DESTROY {
    my $self = shift;
}

1;