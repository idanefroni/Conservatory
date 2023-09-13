#!/usr/bin/perl

package CNS;

use POSIX;
use strict;
use warnings;
use lib './scripts';
use ConservatoryUtils;
use GenomeDatabase;
use Mapping;

my $minSyntenyBiasToAssignCNS= 3; ### The minimum number of species differences between syntenic coding and non-coding region to delete the less
									# syntenic one
                                    
sub new {
    my $class = shift;
    # we can initialize CNS from either a file, a string, or a parameters
    my ($refGenome, $CNSID, $locus, $pos, $length, $conservationLevel, $supportingSpecies, $ancesteralSeq, $regulatorySeqPart);

    if(scalar @_ == 1) {
        my $firstparam= shift;
        ## if this is a file, just read the next line from it
        if($firstparam->isa('IO::Handle')) {
            $firstparam = <$firstparam>;
            chomp($firstparam);
        }
        ## now parse the string
        ($refGenome, $CNSID, $locus, $pos, $length, $conservationLevel, $supportingSpecies, $ancesteralSeq) = split /,/, $firstparam;
    } else {
        ($refGenome, $CNSID, $locus, $pos, $length, $conservationLevel, $supportingSpecies, $ancesteralSeq, $regulatorySeqPart) = @_;
        if(!defined $ancesteralSeq) { $ancesteralSeq = 'N' x $length; }
    }

    if(!defined $pos) { $pos=""; }
    if(!defined $conservationLevel) { $conservationLevel=""; } 
    if(!defined $supportingSpecies) { $supportingSpecies=0; }     

    if($pos ne "" && !defined $regulatorySeqPart) {
        if($pos<0) { $regulatorySeqPart = "U"; } else { $regulatorySeqPart = "D"; }
    }

    my $self = {
        _RefGenome => $refGenome,
	    _ID => $CNSID,
		_Locus => $locus,
		_Pos => $pos,
		_Level => $conservationLevel,
        _SupSp => $supportingSpecies,
		_AncesteralSeq => $ancesteralSeq,
		_RegulatorySeqPart => $regulatorySeqPart,
        _Alive => 1
       };

    bless $self, $class;

    if($refGenome ne "" && $CNSID ne "") {
        return $self;
    } else {
        return undef;
    }
}

sub getRefGenome {
    my ($self) = @_;
    return $self->{_RefGenome};
}

sub getRegulatorySeqPart {
    my ($self) = @_;
    return $self->{_RegulatorySeqPart};
}
sub setRegulatorySeqPart {
    my ($self, $regulatorySeqPart) = @_;
    $self->{_RegulatorySeqPart}= $regulatorySeqPart;
}

sub getRefChr {
    my ($self) = @_;
    if(!defined $self->{_ReferenceMapping}) { die "ERROR: The CNS $self->{_ID} has no absolute reference coordinates set.\n"; }
    return $self->{_ReferenceMapping}->getAbsChr();
}
sub getRefPos {
    my ($self) = @_;
    if(!defined $self->{_ReferenceMapping}) { die "ERROR: The CNS $self->{_ID} has no absolute reference coordinates set.\n"; }
    return $self->{_ReferenceMapping}->getAbsPos();
}

sub getRefEnd {
    my ($self) = @_;
    if(!defined $self->{_ReferenceMapping}) { die "ERROR: The CNS $self->{_ID} has no absolute reference coordinates set.\n"; }
    return $self->{_ReferenceMapping}->getAbsPos() + $self->{_ReferenceMapping}->getAbsLen();
}

sub getRefLocus {
    my ($self) = @_;
    if(!defined $self->{_ReferenceMapping}) { die "ERROR: The CNS $self->{_ID} has no absolute reference coordinates set.\n"; }
    return $self->{_ReferenceMapping}->getLocus();
}

sub getRefStrand {
    my ($self) = @_;
    if(!defined $self->{_ReferenceMapping}) { die "ERROR: The CNS $self->{_ID} has no absolute reference coordinates set.\n"; }
    my $absStrand = $self->{_ReferenceMapping}->getAbsStrand();
    if( $self->{_ReferenceMapping}->getStrand() eq "-") {
        $absStrand = flipStrand($absStrand);
    }
    return $self->{_ReferenceMapping}->getAbsStrand();
}

sub setComment {
    my ($self, $comment) = @_;
    $self->{_Comment} = $comment;
}
sub getComment {
    my ($self) = @_;
    return $self->{_Comment};
}

sub getID {
    my ($self) = @_;
    return $self->{_ID};
}

sub getLocus {
    my ($self) = @_;
    return $self->{_Locus};
}

sub getSeq {
    my ($self) = @_;
    return $self->{_AncesteralSeq};
}


sub setSeq {
    my ($self, $seq) = @_;
    $self->{_AncesteralSeq} = $seq;
}

sub getConservationLevel {
    my ($self) = @_;
    return $self->{_Level};
}

sub getSupportingSpeciesNumber {
    my ($self) = @_;
    return $self->{_SupSp};
}

sub getLen {
    my ($self) = @_;
    return length($self->{_AncesteralSeq});
}

sub setConservationLevel {
    my ($self,$level) = @_;
    $self->{_Level} = $level;
}

sub getSupportingSpecies {
    my ($self, $mappingDB) = @_;
    if(!defined $mappingDB) { die "ERROR: getSupportingSpecies: no mapping database provided.\n"; }
    my %species;
    if(defined $mappingDB->getMappingsForCNS( $self->{_ID}) ) {
        foreach my $curMapping ( @{ $mappingDB->getMappingsForCNS( $self->{_ID})}) {
            $species{ $curMapping->getSpecies() } =1;
        }
        $self->{_SupSp} = scalar keys %species;
        my @speciesArray = keys %species;
        return \@speciesArray;
    } else {
        $self->{_SupSp} =0;
        return undef;
    }
}

sub findReferenceMapping {
    my ($self, $genomeDB, $MappingDB) = @_;

    foreach my $curMapping ( @{ $MappingDB->getMappingsForCNS( $self->{_ID}) }) {
		if($curMapping->getSpecies() eq $genomeDB->genomeToSpecies( $self->{_RefGenome} ) && geneToLocus( $curMapping->getLocus()) eq $self->{_Locus} )  {
			$self->{_ReferenceMapping} = $curMapping;
		}
	}
    return $self->hasReferenceMapping();
}

sub getReferenceMapping {
    my ($self) = @_;
    if(!defined $self->{_ReferenceMapping}) { die "ERROR: CNS $self->{_ID} has no reference mapping. Call findReferenceMapping to set it up.\n"; }

    return $self->{_ReferenceMapping};
}

sub getPos {
    my ($self) = @_;
    return $self->{_Pos};
}

sub hasReferenceMapping {
    my ($self) = @_;
    if(defined $self->{_ReferenceMapping}) {
        return 1;
    } else { 
        return 0;
    }
}

sub hasConservationLevel {
    my ($self) = @_;
    if(defined $self->{_Level} && $self->{_Level} ne "") {
        return 1;
    } else {
        return 0;
    }
}
sub isMerged {
    my ($self) = @_;
    if(index($self->{_ID}, $superCNSPrefix) != -1) {
        return 1;
    } else {
        return 0;
    }
}

sub isAlive {
    my ($self) = @_;
    return $self->{_Alive};
}

sub kill {
    my ($self) = @_;
    $self->{_Alive} =0;
}


#########################################################################################################
#### Compare synteny of a CNS to two genes
#### format: compareSynteny(CNSOne, CNSTwo)
####
####  Returns a number: 1 - CNS one has deeper synteny
####                    2 - CNS two has deeper synteny
####	                0 - both mappings have similar synteny
####

sub compareSynteny {
	my ($self, $CNSTwo, $mappingDB, $conservatoryTree) = @_;
	### First, compare the level of conservation

	## If the number of species between the two CNS-coding pairs is not large enough to determine synteny bias,
	## return here
	if(abs($self->getSupportingSpeciesNumber() - $CNSTwo->getSupportingSpeciesNumber())< $minSyntenyBiasToAssignCNS) {
		return 0;
	}
    ## if we do not have an assigned conservation level, find it
    if(!$self->hasConservationLevel()) {
        my @species = @{ $self->getSupportingSpecies($mappingDB) };
        $self->setConservationLevel($conservatoryTree->findDeepestCommonNode(\@species));
    }

    if(!$CNSTwo->hasConservationLevel()) {
        my @species = @{ $CNSTwo->getSupportingSpecies($mappingDB) };
        $CNSTwo->setConservationLevel($conservatoryTree->findDeepestCommonNode(\@species));
    }

	### Pick the one with deepest conservation. If conservation level is the same, pick the better supported one by number of species
	### And if they are all the same, return 0.
    if($conservatoryTree->getAgeForNode($self->getConservationLevel()) == -1 || $conservatoryTree->getAgeForNode( $CNSTwo->getConservationLevel() ) == -1 ) {
        print "ERROR: Cannot determine node age for:\n";
        $self->print;
        $CNSTwo->print;
        return 0;
    }

	if($conservatoryTree->getAgeForNode($self->getConservationLevel()) > $conservatoryTree->getAgeForNode( $CNSTwo->getConservationLevel() )) {
		return 1;
	} elsif( $conservatoryTree->getAgeForNode($self->getConservationLevel()) < $conservatoryTree->getAgeForNode( $CNSTwo->getConservationLevel() ) ) {
		return 2;
	} elsif( $self->getSupportingSpeciesNumber() > $CNSTwo->getSupportingSpeciesNumber() ) {
		return 1;
	} elsif($CNSTwo->getSupportingSpeciesNumber() > $self->getSupportingSpeciesNumber() ) {
		return 2;
	} else {
		return 0;   
	}
}

sub print {
    my ($self, $outFile) = @_;
    if(!defined $outFile) { $outFile = \*STDOUT; } 

    if(defined $self->{_Comment}) {
    	print $outFile join(",",
					$self->{_RefGenome},
					$self->{_ID},
					$self->{_Locus},
					$self->{_Pos},
					$self->getLen(),
					$self->{_Level},
					$self->{_SupSp},
					$self->{_AncesteralSeq},
                    $self->{_Comment}) ."\n";
    } else {
    	print $outFile join(",",
					$self->{_RefGenome},
					$self->{_ID},
					$self->{_Locus},
					$self->{_Pos},
					$self->getLen(),
					$self->{_Level},
					$self->{_SupSp},
					$self->{_AncesteralSeq}) ."\n";
    }
}

1;
