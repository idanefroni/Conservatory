#!/usr/bin/perl

package GenomeDatabase;

use POSIX;
use strict;
use warnings;
use Cwd 'abs_path';
use lib './scripts';
use ConservatoryUtils;

sub new {
    my $class = shift;
    my ($conservatoryDirectory, $genomeDatabaseFileName) = @_;

    if(!defined $conservatoryDirectory || $conservatoryDirectory eq "") {  ## if not provided, the default file name is genome_database.csv
         $conservatoryDirectory = abs_path(".");
    }

    if(!defined $genomeDatabaseFileName || $genomeDatabaseFileName eq "") {  ## if not provided, the default file name is genome_database.csv
         $genomeDatabaseFileName = "$conservatoryDirectory/genome_database.csv";
    }
     my $self = {
            _conservatoryDirectory => $conservatoryDirectory,
            _genomeDatabaseFileName => $genomeDatabaseFileName
    };

    ## open the genomeDatabase file

    my %genomeToSpeciesHash;
    my %speciesToFamilyHash;    
    my %speciesDatabaseHash;
    my %regUpstreamLengthHash;
    my %regDownstreamLengthHash; 
    my %classificationsHash;

    if(! -e $self->{_genomeDatabaseFileName}) { die "ERROR: Cannot find genome database file $self->{_genomeDatabaseFileName}.\n"; }
    open (my $genomeDatabaseFile, "<", $self->{_genomeDatabaseFileName}) || die "ERROR: Cannot find genome database file $self->{_genomeDatabaseFileName}.\n";
    while(my $curgenomeline = <$genomeDatabaseFile>) {
	    chomp($curgenomeline);
	    if(substr($curgenomeline,0,1) ne "#")  {
		    my ($curGenomeName, $curSpecies,  $curFamily, $curReference, $upstreamLength, $downstreamLength, $geneNameField, $geneProcessingRegEx, $gene2SpeciesIdentifier, $proteinProcessingRegEx, $classification) = split /,/, $curgenomeline;

            $genomeToSpeciesHash{$curGenomeName}= $curSpecies;
            $speciesToFamilyHash{$curSpecies} = $curFamily;
            $speciesDatabaseHash{$curSpecies}{$curGenomeName}=1;

		    $regUpstreamLengthHash{$curGenomeName} = $upstreamLength * 1000;
		    $regDownstreamLengthHash{$curGenomeName} = $downstreamLength * 1000;
		    $classificationsHash{$curGenomeName} = $classification;	
	    }
    }
    close($genomeDatabaseFile);

    $self->{ _genomeToSpeciesHash } = \%genomeToSpeciesHash;
    $self->{ _speciesToFamilyHash } = \%speciesToFamilyHash;
    $self->{ _speciesDatabaseHash } = \%speciesDatabaseHash;
    $self->{ _regUpstreamLengthHash } = \%regUpstreamLengthHash;
    $self->{ _regDownstreamLengthHash } = \%regDownstreamLengthHash;
    $self->{ _classificationsHash } = \%classificationsHash;
    $self->{_TmpDir} = $self->{_conservatoryDirectory} . "/tmp/";
    $self->{_minFamilyThreshold} = 2000;
    $self->{_minFamilyIdentity} = 70;    

    $self->{_minDeepThreshold} = 1600;
    $self->{_minDeepIdentity} = 70;    

    $self->{_minOrthologQuality} = 1000;

    bless $self, $class;
    return $self;
}
sub genomeExists {
    my ($self, $genome) = @_;
    if(defined $self->{_genomeToSpeciesHash}->{$genome}) {
        return 1;
    } else {
        return 0;
    }
}

sub getSpeciesForFamily {
    my ($self, $family) = @_;
    my @speciesForFamily;
    foreach my $curSpecies (keys %{ $self->{_speciesToFamilyHash} }) {
        if($self->{_speciesToFamilyHash}->{$curSpecies} eq $family) {
            push @speciesForFamily, $curSpecies;
        }
    }
    return @speciesForFamily;
}
sub getGenomeNames {
    my ($self) = @_;
    return sort keys %{ $self->{ _genomeToSpeciesHash }};
}

sub getGenomeNamesInFamily {
    my ($self, $family) = @_;
    my @genomesInFamily;
    foreach my $curGenome (keys %{ $self->{ _genomeToSpeciesHash } }) {
        if($self->genomeToFamily($curGenome) eq $family) {
            push @genomesInFamily, $curGenome;
        }
    }
    return sort @genomesInFamily;
}

sub speciesToFamily {
    my ($self, $species) = @_;
    return $self->{_speciesToFamilyHash}->{$species};
}
sub genomeToSpecies {
    my ($self, $genome) = @_;
    return $self->{_genomeToSpeciesHash}->{$genome};
}

sub genomeToFamily {
    my ($self, $genome) = @_;
    return $self->{_speciesToFamilyHash}->{ $self->{_genomeToSpeciesHash}->{$genome} };
}

sub readGenome {
    my ($self, $genomeName) = @_;
    my %genomeSeq;
    if(!defined $self->{_genomeToSpeciesHash}->{$genomeName}) { die "ERROR: Cannot find genome $genomeName in genome database.\n"; }

    my $genomeFastaFileName = $self->getConservatoryDir() . "/genomes/" . $self->genomeToFamily($genomeName) . "/$genomeName.fasta.gz";
    if(! (-e $genomeFastaFileName)) { die "ERROR: Cannot find $genomeName fasta file: $genomeFastaFileName.\n"; }
            
    open (my $genomeFastaFile, "gunzip -c $genomeFastaFileName |");
        my $curChr="";
        while(my $curFastaLine = <$genomeFastaFile>) {
	        chomp($curFastaLine);	
	        if(substr($curFastaLine,0,1) eq ">") {
		        $curChr = cleanChrName($curFastaLine);
	        } else {
		        $genomeSeq{$curChr} .= uc($curFastaLine);
            }
	    }
    close($genomeFastaFile);
    return \%genomeSeq;
}

sub getConservatoryDir {
    my ($self) = @_;
    return $self->{_conservatoryDirectory};
}

sub getTemporaryDir {
    my ($self) = @_;
    return $self->{_TmpDir};
}

sub setTemporaryDir {
    my ($self, $tmpDir) = @_;
    $self->{_TmpDir} = $tmpDir;
}

sub getUpstreamLength {
    my ($self, $genome) = @_;
    return $self->{ _regUpstreamLengthHash }->{$genome};
}
sub getDownstreamLength {
    my ($self, $genome) = @_;
    return $self->{ _regDownstreamLengthHash }->{$genome};
}
sub getClassification {
    my ($self, $genome) = @_;
    return $self->{ _classificationsHash }->{$genome};
}
sub getUpstreamFastaFileName {
    my ($self, $genome) = @_;
    return $self->getConservatoryDir() . "/genomes/" . $self->genomeToFamily($genome) . "/$genome.upstream.fasta.gz";
}
sub getDownstreamFastaFileName {
    my ($self, $genome) = @_;
    return $self->getConservatoryDir() . "/genomes/" . $self->genomeToFamily($genome) . "/$genome.downstream.fasta.gz";
}

sub getLocusFullName {
    my ($self, $genome, $locus) = @_;
    return $self->genomeToSpecies($genome) . "-$genome-$locus";
}
sub getMinFamilyIdentity {
    my ($self) = @_;
    return $self->{_minFamilyIdentity};
}
sub getMinFamilyThreshold {
    my ($self) = @_;
    return $self->{_minFamilyThreshold};
}

sub getMinDeepIdentity {
    my ($self) = @_;
    return $self->{_minDeepIdentity};
}
sub getMinDeepThreshold {
    my ($self) = @_;
    return $self->{_minDeepThreshold};
}

sub setMinFamilyThreshold {
    my ($self, $threshold) = @_;
    $self->{_minFamilyThreshold} = $threshold;
}
sub setMinFamilyIdentity {
    my ($self, $identity) = @_;
    $self->{_minFamilyIdentity} = $identity;
}

sub setMinDeepThreshold {
    my ($self, $threshold) = @_;
    $self->{_minDeepThreshold} = $threshold;
}
sub setMinDeepIdentity {
    my ($self, $identity) = @_;
    $self->{_minDeepIdentity} = $identity;
}

sub getMinOrthologQuality {
    my ($self) = @_;
    return $self->{_minOrthologQuality};
}
1;
