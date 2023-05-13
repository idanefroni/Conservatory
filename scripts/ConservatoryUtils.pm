package ConservatoryUtils;

use POSIX;
use strict;
use warnings;
use List::Util qw(min max uniq);
use Statistics::Basic qw(:all nofill);
use Exporter;

our @ISA= qw( Exporter );
our @EXPORT = qw (overlap geneToSpecies geneToGenome fullNameToShortName lengthWithoutGaps dropAsterixFromProtein findAll getRandomORFLength getLongestORF getCNSbreakpoints polishCNSAlignments);


##############################################################################
##### Parameters

my $standardDeviationsToSplit = 3; ## The number of standard deviation in number of aligned species to trigger a CNS split
my $minSpeciesToSplitCNS=5;  ### Minimal number of species to consider a split of the CNS
my $maxSpeciesToInitiateCNSSplit = 200; ## Do not split CNS if it is supported by atleast this number of species
									 ## To avoid spliting of very highly conserved CNSs

#############################################################################3
### Compute overlap between two fragments. 
### Parameters:
### overlap(startOne, endOne, StartTwo, endTwo, mode)
### returns relative overlap (0-1) (of the first fragment).
###

#sub overlap {
#	my ($startOne, $endOne, $startTwo, $endTwo, $mode) =@_;

#	if ( ($startOne < $startTwo && $endOne < $startTwo) ||
#		 ($startOne > $endTwo)) {
#		return 0;
#	} else {
		#### return the percent overlap
#		if($endOne < $endTwo) {
#			if($startOne >= $startTwo) {
#				return 1;
#			} else {
#				return ($endOne -$startTwo) /($endOne - $startOne);
#			} 
#		} else {
#			if($startOne <= $startTwo) {
#				return 1;
#			} else {
#				return ($endTwo -$startOne) /($endOne - $startOne);
#			} 			
#		}
#	}
#}

#############################################################################3
### Compute overlap between two fragments. 
### Parameters:
### overlap(startOne, endOne, StartTwo, endTwo)
### returns relative overlap (0-1).
###

sub overlap {
	my ($startOne, $endOne, $startTwo, $endTwo) =@_;
	if ( ($startOne < $startTwo && $endOne < $startTwo) ||
		 ($startOne > $endTwo)) {
		return 0;
	} else {
		return (min($endOne, $endTwo) - max($startOne, $startTwo))/ max( $endOne-$startOne, $endTwo-$startTwo);
	}
}

##################################################################################
# Returns species or genome name from the name of a gene

sub geneToSpecies {
	(my $geneName) = @_;
	my @geneNameComponents = split '-', $geneName;
	return $geneNameComponents[0];
}

sub geneToGenome {
	(my $geneName) = @_;
	my @geneNameComponents = split '-', $geneName;
	return $geneNameComponents[1];
}

###################################################################################
# Returns just the locus name from the name of a gene
sub fullNameToShortName {
	(my $geneName) = @_;
	$geneName =~ s/^.[^-]*-[^-]*-//;
	return $geneName;
}

#############################################################
### Returns the length of a sequence without the gaps

sub lengthWithoutGaps {
	my $seq = $_[0];
	return length($seq =~ s/-//rg);
}

########################################################################################
#### Removes the asterix that some programs use to mark stop and other programs hate

sub dropAsterixFromProtein {
	my ($seqObj) = @_;
	my $seq = $seqObj->seq();
	$seq =~ s/\*//g;
	$seqObj->seq($seq);
	return $seqObj;
}

#################################################################
### Find all occurances of a string in another string

sub findAll {
	my ($seq, $substring) = @_;
	my @occurences;
	my $offset=0;
	
	my $result = index($seq, $substring, $offset);

	while( $result != -1) {
		push (@occurences, $result);
		$offset = $result + 1;
		$result = index($seq, $substring, $offset);
	}
	return @occurences;
}

## Returns the mean value of an array of integers
#sub mean {
#    return sum(@_)/@_;
#}


###################### ORF Processing functions
#### getLongestORF(DNSSequence, Strand) 
#####
#### Strand is optional (+ or -) to limit search to one strand. If not provided, ORFs are searched on both strands.
####

sub getLongestORF {
   my ($seq,$strandToLookAt) = @_;
   
   if(!defined $strandToLookAt) { $strandToLookAt = ""; }
   $seq =~ s/-//g;
   
   my $best=0;
   my ($bests,$beste,$beststrand)=(-1,-1,0);
   my $dna=Bio::Seq->new(-seq => $seq);
   my %strand;
   
   if ($strandToLookAt eq "+") {
   	   %strand=('+'=>$dna->seq);
   } elsif($strandToLookAt eq "-") {
  	   %strand=('-'=>$dna->revcom->seq);
   } else {
   	   %strand=('+'=>$dna->seq,
               '-'=>$dna->revcom->seq);
   }
   
   foreach my $direction (keys %strand) {
      my @starts=();
      my @ends=();
      for (my $frame=0;$frame< length($seq)-3;$frame++) {
         unless ($strand{$direction}=~m/^.{$frame}(taa|tga|tag)/i) {
            push @starts,$frame+1;
         }
      }
 
      while ($strand{$direction}=~m/(taa|tga|tag)/gi) {
         push @ends,pos($strand{$direction})-2;
      }
      push @ends,($dna->length+2,$dna->length+1,$dna->length);


      for my $s (@starts) {
         for my $e (@ends) {
            if ($e%3==$s%3 and $e>$s) {
               my $ORFlength = $e-$s;
               if($e > $dna->length) {
                  $ORFlength -= ($e - $dna->length -1);
                  ### If we have partial codon - remove it
                  $ORFlength -= ($ORFlength%3);
               }
               if ($ORFlength >$best) {
                  $best= $ORFlength;
                  ($bests,$beste,$beststrand)=($s,$e,$direction);
               }
               last
            } else {
               next
            }
         }
      }
   }
   my $bestORF;
   if($beststrand eq "+") {
   	   $bestORF = substr($dna->seq,$bests-1, $best)
   } else {
   	   $bestORF = substr($dna->revcom->seq,$bests-1, $best)
   }
   
   return (length($bestORF), $bestORF, $beststrand ); ## Returns the length of the longest ORF and its sequence	
}


##### determine the expected random ORF based on GC content of cns
### Based on Oliver and Marin, Journal of Molecular Evolution, 1996 43:216-223.
###
sub getRandomORFLength {
	my $LENGTH_FACTOR = 2;
	my $randomORFLength;
	
	my $sequence = $_[0];
	my $seqstr = $sequence;
	my $gcCNS = length($seqstr =~ s/[A|T|N]//rg);
	$seqstr = $sequence;
	my $atCNS = length($seqstr =~ s/[G|C|N]//rg);

	my $q= (1-$gcCNS/($gcCNS+$atCNS))/2;
	my $t= $q**2 - $q**3;
	if($t == 0) {
		$randomORFLength=1000;
	} else {
		$randomORFLength = ceil(1/$t)*3*$LENGTH_FACTOR+1;
	}
	return $randomORFLength;
}

############################################################################################################3
#### CNS processing functions

### Accepts an array of hits hashes and returns the breakpoints

sub getCNSbreakpoints {
	my ($CNSLength, $hitsRef, $minCNSLength) = @_;
	my @hits = @$hitsRef;
	my @CNSCoverage = (0) x $CNSLength;
	my @CNSCoverageSmooth = (0) x $CNSLength;

	foreach my $curHit (@hits) {
		my $start = $curHit->{'ReferenceRelativePosition'};
		my $end = $start + $curHit->{'Length'};

		foreach my $pos ($start..$end) {
			my $curNucleotide = substr($curHit->{'TargetSequence'}, ($pos-$start),1);
			if($curNucleotide ne "-" && $curNucleotide ne "N") {
				$CNSCoverage[$pos]++;
			}
		}
	}
	### Smooth the coverage vector to avoid spliting on very small gaps
	foreach my $pos (0..($CNSLength-4)) {
		$CNSCoverageSmooth[$pos] = mean(@CNSCoverage[($pos)..($pos+3)])
	}
	## The last positions cannot be smoothed
	foreach my $pos (($CNSLength-4)..($CNSLength-1)) {
		$CNSCoverageSmooth[$pos] = $CNSCoverage[$pos];
	}

	### Change point detection. Our data is generally well behaved with few outliers (if any)
	### we do a simple change detection algorithm. 
	#####  1. Calculate delta between points. 2. Find the stdev of the delta. 3. Identify places where delta is high (3xstdev),
	#####    Don't break if creating too small of a fragment (<$minCNSLength)

	#### Calculate the dX
	my @CNSdelta = (0) x $CNSLength;
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
		if((abs($CNSdelta[$pos]) > $deltaSD * $standardDeviationsToSplit) && (abs($CNSdelta[$pos])> $minSpeciesToSplitCNS ) && ($pos + $minCNSLength < $CNSLength) &&
		 ($pos - $lastBreakpointPos > $minCNSLength && ($CNSLength - $pos) > $minCNSLength ) ) {
			push @breakpoints, $pos;
			$lastBreakpointPos = $pos;
		}
	}
	# Add start and end points
	@breakpoints = (0, @breakpoints , $CNSLength);
	return \@breakpoints;
}

sub polishCNSAlignments {
	my ($breakpointRef, $alignmentMapRef, $minCNSConservationAfterSplit) = @_;
	my @breakpoints = @$breakpointRef; 
	my @alignmentMap = @$alignmentMapRef;
	my @splitAndPolishedAlignments;

	foreach my $curHit (@alignmentMap) {
		my $start =  $curHit->{'ReferenceRelativePosition'};
		my $end = $start + $curHit->{'Length'};

		for my $curBreakpoint (0..(scalar @breakpoints - 2) ) {
			### if the hit overlaps the sub CNS
#			print "Testing polish: " . $curHit->{'TargetLocus'} . ": $breakpoints[$curBreakpoint] - " . $breakpoints[$curBreakpoint+1] . " and " . $curHit->{'ReferenceRelativePosition'} . "-" . $curHit->{'Length'}  . "\n";

			if( overlap($breakpoints[$curBreakpoint], $breakpoints[$curBreakpoint+1], $curHit->{'ReferenceRelativePosition'}, $curHit->{'ReferenceRelativePosition'} + $curHit->{'Length'} )) {

				### distance from breakpoint start
				my $positionInCNS = max($breakpoints[$curBreakpoint], $curHit->{'ReferenceRelativePosition'});

			#	print "editing seq: " . $curHit->{'TargetSequence'} . ". start is $start. Cutting from " . $breakpoints[$curBreakpoint] . "\n";
				my $subCNSSeqStart = max(0,$breakpoints[$curBreakpoint] - $start);
				my $subCNSTargetSeq = substr($curHit->{'TargetSequence'}, $subCNSSeqStart , min($breakpoints[$curBreakpoint+1]-$breakpoints[$curBreakpoint], $breakpoints[$curBreakpoint+1] - $curHit->{'ReferenceRelativePosition'}, $curHit->{'Length'} - $breakpoints[$curBreakpoint]+$start ));
				### remove trailing and leading gaps
			#	print "-- Possible hit Seq: $subCNSTargetSeq. Trimmed: ";	
				$subCNSTargetSeq =~ s/(-+)$//;
				my $numOfLeadingGaps = $subCNSTargetSeq =~ s/^(-+)//;
				if($numOfLeadingGaps eq "") { $numOfLeadingGaps = 0; }

				my $numOfInternalGaps = $subCNSTargetSeq =~ tr/-//;
				if($numOfInternalGaps eq "") { $numOfInternalGaps = 0; }

			#	print " $subCNSTargetSeq ($numOfInternalGaps gaps).\n";

				### Only include the hit of the coverage is sufficient
				if((length($subCNSTargetSeq)-$numOfInternalGaps) / ($breakpoints[$curBreakpoint+1] - $breakpoints[$curBreakpoint] ) > $minCNSConservationAfterSplit) {
			#		print "-- hit at $positionInCNS. leading gaps $numOfLeadingGaps. new length " . length($subCNSTargetSeq) . "\n";	
					my $newTargetPosition;  ## update the target coordinates - byt that depends if the hit is on the plus or minus strand
					#### DEBBUG THIS
					if($curHit->{'TargetStrand'} eq "+" ) {
						$newTargetPosition = $curHit->{'TargetPosition'} + $breakpoints[$curBreakpoint] + $numOfLeadingGaps
					} else {
						$newTargetPosition = $curHit->{'TargetPosition'} + $curHit->{'Length'} - length($subCNSTargetSeq) - $subCNSSeqStart;
					}

					push (@splitAndPolishedAlignments, {
							"TargetSpecies" => $curHit->{'TargetSpecies'},
							"TargetLocus" => $curHit->{'TargetLocus'},
							"TargetStrand" => $curHit->{'TargetStrand'},
							"Length" => length($subCNSTargetSeq),
							"ReferenceUpDown" => $curHit->{'ReferenceUpDown'},
							"ReferenceLocus" => $curHit->{'ReferenceLocus'},
							"TargetSequence" => $subCNSTargetSeq,
							"ReferenceRelativePosition" => $positionInCNS+$numOfLeadingGaps,
							"TargetPosition" => $newTargetPosition 
					});
				}
			}
		}
	}
	return \@splitAndPolishedAlignments;
}
