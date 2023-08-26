#!/usr/bin/perl

package MappingDatabase;

use POSIX;
use strict;
use warnings;
use lib './scripts';
use ConservatoryUtils;
use Statistics::Basic qw(:all nofill);

use Mapping;

sub new {
    my ($class, $MappingDBFileName, $verbose) = @_;
    if(!defined $verbose) { $verbose =0; }
    my @mappingDB =();
    my %mappingCNSIndex;

    # we can initialize a mapping database from file

    if(defined $MappingDBFileName && $MappingDBFileName ne "") {

        ### now read the CNS DB

        my $curPos=1;
        open(my $MappingInputFile, $MappingDBFileName) || die "ERROR: Cannot open mapping database file $MappingDBFileName.\n";
        while(<$MappingInputFile>) {
            chomp;
	        if($verbose) { print "PROGRESS: Reading Mappings..." . ($curPos++) . ".\r" };

            my $newMapping = new Mapping($_);
            push(@mappingDB, $newMapping);
            push( @{ $mappingCNSIndex{ $newMapping->getCNSID() } } , $newMapping);
        }
        if($verbose) { print "\n"; };
    }

    my $self = {
        _MappingDB => \@mappingDB,
        _mappingCNSIndex => \%mappingCNSIndex
    };

    bless $self, $class;
    return $self;
}

sub getMappingsForCNS {
    my ($self, $CNSID) = @_;
    if(ref($CNSID) eq 'CNS') { $CNSID = $CNSID->getID(); }
    if(defined $self->{_mappingCNSIndex}->{$CNSID}) {
        return $self->{_mappingCNSIndex}->{$CNSID};
    } else {
        my @dummy=();
        return \@dummy;
    }
}
sub getNumberOfMappings {
    my ($self) = @_;
    return scalar @{ $self->{_MappingDB} };
}
sub addMapping {
    my ($self, $newMapping) = @_;
    if(!defined $newMapping->getCNSID() ) { die "ERROR: No CNSID for mapping. Cannot assign to database.\n"; }
    push @{ $self->{_MappingDB} }, $newMapping;
    push @{ $self->{_mappingCNSIndex}->{$newMapping->getCNSID()} }, $newMapping;
}

sub linkMappingToCNS {
    my ($self, $mapping, $CNSID) = @_;
    if(ref($CNSID) eq 'CNS') { $CNSID = $CNSID->getID(); }
    $mapping->setCNSID($CNSID);
    push @{ $self->{_mappingCNSIndex}->{$CNSID} }, $mapping; 
}

sub deleteMapping {
    my ($self, $mappingToDelete) = @_;
    my @newMappingsForCNS = grep { $_ != $mappingToDelete} @{ $self->getMappingsForCNS( $mappingToDelete->getCNSID() ) };
    $self->{_mappingCNSIndex}->{ $mappingToDelete->getCNSID() } = \@newMappingsForCNS; 
    $mappingToDelete->unlink(); ## unlink mapping    
}

sub deleteCNS {
    my ($self, $CNSID) = @_;
    ### accepts a CNS object or a CNSID string
    if(ref($CNSID) eq 'CNS') { $CNSID = $CNSID->getID(); }

    ##unlinkmappings
    foreach my $curMapping (@{ $self->getMappingsForCNS($CNSID) }) {
        if($curMapping->getCNSID() eq $CNSID) {
            $curMapping->unlink();
        }
    }
    delete $self->{_mappingCNSIndex}->{$CNSID};
}

sub orderMappingByPosition {
    my ($self) = @_;
    my @mappingDBDeref = @{ $self->{_MappingDB} };

    @mappingDBDeref = sort { geneToGenome($a->getLocus()) cmp geneToGenome($b->getLocus()) ||
				  		   $a->getAbsChr() cmp $b->getAbsChr() ||
				  		   $a->getAbsPos() <=>  $b->getAbsPos() } @mappingDBDeref;
    $self->{_MappingDB} = \@mappingDBDeref;

}

sub orderMappingByGenome {
    my ($self) = @_;
    my @mappingDBDeref = @{ $self->{_MappingDB} };

    @mappingDBDeref = sort { geneToGenome($a->getLocus()) cmp geneToGenome($b->getLocus())} @mappingDBDeref;
    $self->{_MappingDB} = \@mappingDBDeref;

}

sub getMappingsByOrder {
    my ($self) = @_;
    return $self->{_MappingDB};
}
sub renameMappings {
    my ($self, $verbose) = @_;

    if($verbose) { print "PROGRESS: Sorting mappings..."; }
    ## deref
    $self->orderMappingByPosition();
    if($verbose) { print "Done.\n"; }

    my @mappingDBDeref = @{ $self->{_MappingDB} };
    my $lastMapping = $mappingDBDeref[0];
    my $nameCounter=1;
    my $curPos=2;
    $lastMapping->setName("M_" . $nameCounter);

    foreach my $curMapping ( @{ $self->{_MappingDB} } [1..((scalar @{ $self->{_MappingDB} }) -1 )]) {
	    if($verbose) { print "PROGRESS: Naming..." . ($curPos++) . "/" . (scalar @{ $self->{_MappingDB} }) . ".\r"; };

	    if( geneToGenome($lastMapping->getLocus) ne geneToGenome($curMapping->getLocus()) ) {
		    $nameCounter=1;
	    } elsif (! $lastMapping->overlap($curMapping) ) {
		    $nameCounter++;
	    } 
	    $curMapping->setName("M_" . $nameCounter);
	    $lastMapping = $curMapping;
    }
    if ($verbose) { print "\n"; }
}

sub removeDuplicatesForCNS {
    my ($self, $CNSID)= @_;
    if(ref($CNSID) eq 'CNS') { $CNSID = $CNSID->getID(); }

    my @mappingsForCNS = @{ $self->getMappingsForCNS($CNSID)  };

    @mappingsForCNS = sort { geneToGenome($a->getLocus()) cmp geneToGenome($b->getLocus()) ||
				  		   $a->getAbsChr() cmp $b->getAbsChr() ||
				  		   $a->getAbsPos() <=>  $b->getAbsPos() } @mappingsForCNS;
        
    if(scalar @mappingsForCNS >1) {
        my $lastMapping = $mappingsForCNS[0];
        foreach my $curMapping (@mappingsForCNS[1.. ((scalar @mappingsForCNS) -1)]) {
            if($lastMapping->identical($curMapping)) {
                $self->deleteMapping($curMapping);
            } else {
                $lastMapping = $curMapping;
            }
        }
    }
}

sub mergeOverlappingMappingsForCNS {
	my ($self,$CNSID) = @_;
    if(ref($CNSID) eq 'CNS') { $CNSID = $CNSID->getID(); }

    my @mappingsForCNS = @{ $self->getMappingsForCNS($CNSID)  };

    @mappingsForCNS = sort { geneToGenome($a->getLocus()) cmp geneToGenome($b->getLocus()) ||
				  		   $a->getAbsChr() cmp $b->getAbsChr() ||
				  		   $a->getAbsPos() <=>  $b->getAbsPos() } @mappingsForCNS;

    my $lastMapping = $mappingsForCNS[0];
    if(scalar @mappingsForCNS >1) {
        foreach my $curMapping (@mappingsForCNS[1..(scalar @mappingsForCNS-1)]) {
            if($curMapping->isAlive() && $lastMapping->overlap($curMapping)) {
 #               if($curMapping->getSpecies() eq "Atauschii") {
 #                   print "Merging:\n";
  #                  $lastMapping->print;
   #                 $curMapping->print;
    #            }
                $lastMapping->merge($curMapping);
     #           if($curMapping->getSpecies() eq "Atauschii") {
      #              print "Merged:\n";
       #             $lastMapping->print;
        #        }

                $self->deleteMapping($curMapping);
            } else {
                $lastMapping = $curMapping;
            }
        }
    }
}

#################################################################
#### Accepts two CNSs. Copy the mappings from CNS #2 to CNS #1
sub mergeMappings {
    my ($self, $CNSTarget, $CNSSource, $bias) = @_;

    if(!defined $bias) { $bias=0;}
    foreach my $mapping (@{ $self->getMappingsForCNS( $CNSSource) }) {
        $self->linkMappingToCNS($mapping, $CNSTarget);
	    $mapping->setRRP( $mapping->getRRP() + $bias);
	}
}

############################################################################################################3
#### CNS processing functions
### Find possible breakpoints for a CNS
###
###  sub CNS are created where there is a drop of >3 standard deviations in the number of alignments (and atleast 5) for the nucleotide, as long as the CNS is larger than minCNSLength.
###

sub findBreakPointsForCNS {
	my ($self, $CNS) = @_;

	my @mappings = @{ $self->getMappingsForCNS($CNS) };
	my @CNSCoverage = (0) x $CNS->getLen();
	my @CNSCoverageSmooth = (0) x $CNS->getLen();

	foreach my $curMapping (@mappings) {
		my $start = $curMapping->getRRP();
		my $end = $start + $curMapping->getLen();

		foreach my $pos ($start..$end) {
			my $curNucleotide = substr($curMapping->getSeq(), ($pos-$start),1);
			if($curNucleotide ne "-" && $curNucleotide ne "N") {
				$CNSCoverage[$pos]++;
			}
		}
	}
	### Smooth the coverage vector to avoid spliting on very small gaps
	foreach my $pos (0..($CNS->getLen()-4)) {
		$CNSCoverageSmooth[$pos] = mean(@CNSCoverage[($pos)..($pos+3)])
	}
	## The last positions cannot be smoothed
	foreach my $pos (($CNS->getLen()-4)..($CNS->getLen()-1)) {
		$CNSCoverageSmooth[$pos] = $CNSCoverage[$pos];
	}

	### Change point detection. Our data is generally well behaved with few outliers (if any)
	### we do a simple change detection algorithm. 
	#####  1. Calculate delta between points. 2. Find the stdev of the delta. 3. Identify places where delta is high 
	#####    Don't break if creating too small of a fragment (<$minCNSLength)

	#### Calculate the dX
	my @CNSdelta = (0) x $CNS->getLen();
	foreach my $pos (0..((scalar @CNSCoverageSmooth)-2)) {
		$CNSdelta[$pos+1] = $CNSCoverageSmooth[$pos+1] - $CNSCoverageSmooth[$pos];
	}
	
	my $deltaSD = stddev(@CNSdelta);

	## Now identify breakpoints.
	###
	my @breakpoints;
	my $lastBreakpointPos=0;

	foreach my $pos (0..(scalar @CNSdelta -1)) {
		### split CNS if larger than minCNSlength and high variability. But don't split CNS if its too short, has too many species
		## or if the split point is not well supported
		if((abs($CNSdelta[$pos]) > $deltaSD * $standardDeviationsToSplit) && (abs($CNSdelta[$pos])> $minSpeciesToSplitCNS ) && ($pos + $minCNSLength < $CNS->getLen()) &&
		 ($pos - $lastBreakpointPos > $minCNSLength && ($CNS->getLen() - $pos) > $minCNSLength ) && ($CNSCoverage[$pos+1] >= $minSpeciesForCNS)) {
			push @breakpoints, $pos;
			$lastBreakpointPos = $pos;
		}
	}
	# Add start and end points
	@breakpoints = (0, @breakpoints , $CNS->getLen());

	return @breakpoints;
}


########################################################################################################################################
########
########  Given a set of hits and breakpoints, split the hits so they will be aligned to the breakpoints as closely as possible
########   filter out leftover hits that have low identity
########
########

sub assignMappingsToBreakpoints {
	my ($self, $CNS, $breakpointRef) = @_;
	my @breakpoints = @$breakpointRef;
    my @mappings =  @{ $self->getMappingsForCNS($CNS) };
    my @assignedMappings;

	foreach my $curMapping ( @mappings ) {

		for my $curBreakpoint (0..(scalar @breakpoints - 2) ) {
			### if the hit overlaps the sub CNS

			if( $curMapping->overlap($breakpoints[$curBreakpoint], $breakpoints[$curBreakpoint+1])) {
                my $subsetSeq = $curMapping->getSubsetSeq($breakpoints[$curBreakpoint],$breakpoints[$curBreakpoint+1] );
                ### count internal gaps. First remove leading and trailing
                $subsetSeq =~ s/^(-+)//;
                $subsetSeq =~ s/(-+)$//;                
                my $numOfInternalGaps = () = ($subsetSeq =~ /-/g);
                ### Only include the hit of the coverage is sufficient and if its not too gappy
				if((length($subsetSeq)-$numOfInternalGaps) / ($breakpoints[$curBreakpoint+1] - $breakpoints[$curBreakpoint] ) > $minCNSCoverageAfterSplit && length($subsetSeq) >= $minCNSLength && ($numOfInternalGaps / length($subsetSeq)) < $minSequenceContentInAlignment ) {

                    my $breakpointMapping = $curMapping->createSubSetMapping( $breakpoints[$curBreakpoint],$breakpoints[$curBreakpoint+1] , $CNS->getLen());
                    $breakpointMapping->setBreakpoint($curBreakpoint);
                    $self->linkMappingToCNS($breakpointMapping, $CNS);
                    push(@assignedMappings, $breakpointMapping);
                }
            }
        }
    }

    ### Now delete the original mappings
    foreach my $curMapping (@mappings) {
        $self->deleteMapping($curMapping);
    }
    return @assignedMappings;              
}


sub writeDatabase {
    my ($self, $outputFileName) = @_;
    open(my $outputFile, ">", $outputFileName);
    foreach my $curCNS ( keys %{ $self->{_mappingCNSIndex} }) {
        $self->removeDuplicatesForCNS($curCNS);
        foreach my $curMapping (@{ $self->{_mappingCNSIndex}->{$curCNS} }) {
            $curMapping->print($outputFile);
        }
    }
    close($outputFile);
}


1;
