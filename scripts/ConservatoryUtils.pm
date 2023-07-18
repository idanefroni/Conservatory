package ConservatoryUtils;

use POSIX;
use strict;
use warnings;
use List::Util qw(min max uniq);
use Statistics::Basic qw(:all nofill);
use Bio::SimpleAlign;

use Exporter;

our @ISA= qw( Exporter );
our @EXPORT = qw (overlap geneToSpecies geneToGenome isGeneName fullNameToShortName lengthWithoutGaps dropAsterixFromProtein findAll getRandomORFLength getLongestORF getCNSbreakpoints polishCNSAlignments getGappiness shiftTargetCoordinate translateRealtiveToAbsoluteCoordinates getGeneCoordinates);


##############################################################################
##### CNS Splitting Parameters

my $standardDeviationsToSplit = 3; ## The number of standard deviation in number of aligned species to trigger a CNS split
my $minSpeciesToSplitCNS=5;  ### Minimal number of species to consider a split of the CNS
my $maxSpeciesToInitiateCNSSplit = 200; ## Do not split CNS if it is supported by atleast this number of species
									 ## To avoid spliting of very highly conserved CNSs
my $minIdentityToKeepBreakpoint	=50;
my $minSpeciesToKeepBreakpoint	= 5;
my $minSequenceContentInAlignment = 0.5; 

my %footprintDatabase;   ### Database of foot prints for the getGeneCoordiantes function


#############################################################################
### Compute overlap between two fragments. 
### Parameters:
### overlap(startOne, endOne, StartTwo, endTwo)
### returns relative overlap (0-1).
###

sub overlap {
	my ($startOne, $endOne, $startTwo, $endTwo) =@_;
	if ( ($startOne < $startTwo && $endOne < $startTwo) ||
		 ($startOne > $endTwo) || ($startTwo > $endOne)) {
		return 0;
	} else {
		return (min($endOne, $endTwo) - max($startOne, $startTwo))/ max( $endOne-$startOne, $endTwo-$startTwo);
	}
}

##################################################################################
# Returns species or genome name from the name of a gene

sub isGeneName { 
	my ($geneNameToCheck) = @_;
	if($geneNameToCheck =~ /^.+-.+-.+$/) {
		return 1;
	} else {
#		print " --------- $geneNameToCheck is not a gene name.\n";
		return 0;
	}
}
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

##################################################################################
#####
##### Gappiness - return the percentrage of gaps for each position in a simplealign object
sub getGappiness {
	my ($align) = @_;
	my @gapPercentages;

    for (my $pos = 1; $pos <= $align->length; $pos++) {
		my $gap_count = 0;
        
        # Iterate over all sequences at once
        foreach my $seq ($align->each_seq) {
            my $symbol = $seq->subseq($pos, $pos);
            $gap_count++ if $symbol eq '-';
        }
        
        my $gapPercentageForPos = ($gap_count / $align->num_sequences) * 100;
        push @gapPercentages, $gapPercentageForPos;
    }
    
	return @gapPercentages;
}

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

##############################################################################################
##### determine the expected random ORF based on GC content of CNS
###
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
	my $t= $q**2 - $q**3;   ## Olive and Marin, JME
	if($t == 0) {  ### maximal value for ORF length
		$randomORFLength=1000;
	} else {
		$randomORFLength = ceil(1/$t)*3*$LENGTH_FACTOR+1;
	}
	return $randomORFLength;
}

############################################################################################################3
#### CNS processing functions

### Accepts CNS length, an array of hits hashes and returns the breakpoints to sub CNSs.
###
###  sub CNS are created where there is a drop of >3 standard deviations in the number of alignments (and atleast 5) for the nucleotide, as long as the CNS is larger than minCNSLength.
###


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
	#####  1. Calculate delta between points. 2. Find the stdev of the delta. 3. Identify places where delta is high 
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
	return @breakpoints;
}


########################################################################################################################################
########
########  Given a set of hits and breakpoints, split the hits so they will be aligned to the breakpoints as closely as possible
########   filter out leftover hits that have low identity
########
########

sub polishCNSAlignments {
	my ($breakpointRef, $alignmentMapRef, $minCNSConservationAfterSplit, $minCNSLength) = @_;
	my @breakpoints = @$breakpointRef;
	my @alignmentMap = @$alignmentMapRef;
	my @splitAndPolishedAlignments;

	### set up the alignment summary
	my %alignmentsForBreakpoints;
	for my $curBreakpoint (0..(scalar @breakpoints - 2) ) {
		$alignmentsForBreakpoints{$curBreakpoint} = Bio::SimpleAlign->new();
	}

	foreach my $curHit (@alignmentMap) {
		my $start =  $curHit->{'ReferenceRelativePosition'};
		my $end = $start + $curHit->{'Length'};

		for my $curBreakpoint (0..(scalar @breakpoints - 2) ) {
			### if the hit overlaps the sub CNS

			if( overlap($breakpoints[$curBreakpoint], $breakpoints[$curBreakpoint+1], $curHit->{'ReferenceRelativePosition'}, $curHit->{'ReferenceRelativePosition'} + $curHit->{'Length'} )) {


				### distance from breakpoint start
				my $positionInCNS = max($breakpoints[$curBreakpoint], $curHit->{'ReferenceRelativePosition'});
				my $subCNSSeqStart = max(0,$breakpoints[$curBreakpoint] - $start);
				my $subCNSTargetSeq = substr($curHit->{'TargetSequence'}, $subCNSSeqStart , min($breakpoints[$curBreakpoint+1]-$breakpoints[$curBreakpoint], $breakpoints[$curBreakpoint+1] - $curHit->{'ReferenceRelativePosition'}, $curHit->{'Length'} - $breakpoints[$curBreakpoint]+$start ));
				my $newReferenceSeq = substr($curHit->{'ReferenceSequence'}, $subCNSSeqStart , min($breakpoints[$curBreakpoint+1]-$breakpoints[$curBreakpoint], $breakpoints[$curBreakpoint+1] - $curHit->{'ReferenceRelativePosition'}, $curHit->{'Length'} - $breakpoints[$curBreakpoint]+$start ));

				### remove trailing and leading gaps
				$subCNSTargetSeq =~ s/(-+)$//;
 			    my $numOfTrailingGaps = length($1);
       
				$subCNSTargetSeq =~ s/^(-+)//;
				my $numOfLeadingGaps = length($1);

				if(!defined $numOfLeadingGaps) { $numOfLeadingGaps = 0; }
				if(!defined $numOfTrailingGaps) { $numOfTrailingGaps = 0; }

				my $numOfInternalGaps = ($subCNSTargetSeq =~ tr/-//) + 0;

				### Only include the hit of the coverage is sufficient and if its not too gappy
				if((length($subCNSTargetSeq)-$numOfInternalGaps) / ($breakpoints[$curBreakpoint+1] - $breakpoints[$curBreakpoint] ) > $minCNSConservationAfterSplit && length($subCNSTargetSeq) >= $minCNSLength && ($numOfInternalGaps / length($subCNSTargetSeq)) < $minSequenceContentInAlignment ) {
					
					## update the target coordinates
					my $newTargetPosition= shiftTargetCoordinate($curHit->{'TargetPosition'}, $curHit->{'TargetStrand'}, $curHit->{'TargetSequence'}, $curHit->{'ReferenceSequence'}, $subCNSSeqStart, length($subCNSTargetSeq));
					if($curHit->{'TargetSpecies'} eq "Slycopersicum") {
						print "DEBUG: " . $curHit->{'TargetLocus'} . ":" . $curHit->{'TargetPosition'} . ":" . $curHit->{'Length'} . ":" . $curHit->{'TargetSequence'} . "\n";
						print "DEBUG: " . $breakpoints[$curBreakpoint] . ":" . $curHit->{'ReferenceRelativePosition'} . ":" . $curHit->{'ReferenceSequence'} . "\n";
						print "DEBUG: positionInCNS: $positionInCNS subCNSSeqStart: $subCNSSeqStart newsubCNSTargetSet $subCNSTargetSeq\n";
						print "DEBUG: newTarget: $newTargetPosition.\n";
					}


					my $alignedSeq = ('-' x ($positionInCNS-$breakpoints[$curBreakpoint]+$numOfLeadingGaps)) . $subCNSTargetSeq;
          			### pad end of alignments with gaps to make it equal length
		          	$alignedSeq = $alignedSeq . ('-' x ($breakpoints[$curBreakpoint+1] - $breakpoints[$curBreakpoint] - length($alignedSeq)));
          
					### Log the alignment to the breakpointlist
					$alignedSeq =~ tr/-/N/;
					my $seq = Bio::LocatableSeq->new(-seq => $alignedSeq,
													 -start => 1,
													 -end => length($alignedSeq),
												     -id => $curHit->{'TargetLocus'}. $curHit->{'TargetPosition'});
					$alignmentsForBreakpoints{$curBreakpoint}->add_seq($seq);

					push (@splitAndPolishedAlignments, {
							"TargetSpecies" => $curHit->{'TargetSpecies'},
							"TargetLocus" => $curHit->{'TargetLocus'},
							"TargetStrand" => $curHit->{'TargetStrand'},
							"Length" => length($subCNSTargetSeq),
							"ReferenceUpDown" => $curHit->{'ReferenceUpDown'},
							"ReferenceLocus" => $curHit->{'ReferenceLocus'},
							"TargetSequence" => $subCNSTargetSeq,
							"ReferenceSequence" => $newReferenceSeq,
							"ReferenceRelativePosition" => $positionInCNS+$numOfLeadingGaps,
							"TargetPosition" => $newTargetPosition,
							"Breakpoint" => $curBreakpoint
					});
				}
			}
		}
	}

	my %breakpointsToDelete;
	for my $curBreakpoint (0..(scalar @breakpoints - 2) ) {
		if( $alignmentsForBreakpoints{$curBreakpoint}->percentage_identity() < $minIdentityToKeepBreakpoint ||  $alignmentsForBreakpoints{$curBreakpoint}->num_sequences < $minSpeciesToKeepBreakpoint ) {
			$breakpointsToDelete{$curBreakpoint}=1;
		}

	}

	my @splitAndPolishedAlignmentsFiltered;
	foreach my $curAlignment (@splitAndPolishedAlignments) {
		if(!defined $breakpointsToDelete{ $curAlignment->{'Breakpoint'} }) {
			push(@splitAndPolishedAlignmentsFiltered, $curAlignment);
		}
	}

	return @splitAndPolishedAlignmentsFiltered;
}

###################################################################################
####
#### returns the absolute coordinates for a gene in conservatory
####  accepts the conservatorydirectory (to find the footprint files), family name, genome name and full gene name
####
####   returns a hash of values {Chromosome, Strand, Start, End}
####

### The function extract the values from the footprint file which it maintains in memory. The first time a gene is requested
###  from a genome, the whole foot print file is read for that genome. This is to speed up access.

sub getGeneCoordinates {
	my ($conservatoryDir, $family, $genome, $gene) = @_;
	
	if(!defined $footprintDatabase{$genome}) { ### If we didn't see this genome before, load it to memory
		my $footprintFileName = "$conservatoryDir/genomes/$family/$genome.footprint.gff3";
		open (my $footprintFile, "<", $footprintFileName) || die ("ERROR: Cannot open footprint file $footprintFileName.\n");
		while (my $line = <$footprintFile>){
			chomp($line);
			my @array = split /\t/, $line;
			my %fields = split /[;=]/, $array[8];
			my $footprintGeneName = $fields{'Name'};
			$footprintDatabase{$genome}{$footprintGeneName}{'Strand'} = $array[6];
			$footprintDatabase{$genome}{$footprintGeneName}{'Chromosome'} = $array[0];
			$footprintDatabase{$genome}{$footprintGeneName}{'Start'} = $array[3];
			$footprintDatabase{$genome}{$footprintGeneName}{'End'} = $array[4];			
		}
		close($footprintFile);
	}

	return %{ $footprintDatabase{$genome}{$gene} };
}


#######################################
##### Translate a coordiante list of segments to from relative to absolute coordinates.
##### Accepts a hash with "Start" value having the coordiante and "Length" having the length of the segment or an array of such hashes.

sub translateRealtiveToAbsoluteCoordinates {
	my ($CNSCoordinatesRef, $absGeneStart, $absGeneEnd, $strand) = @_;
	my @absCoordinates;
	my @relCoordiantes;
	if(ref($CNSCoordinatesRef) eq "ARRAY") {
		@relCoordiantes = @$CNSCoordinatesRef;
	} elsif(ref($CNSCoordinatesRef) eq "HASH") {
		@relCoordiantes = (%$CNSCoordinatesRef);
	} else {
		die "ERROR: translateRelativeToAbsoluteCoordiante accepts either a HASH or an array of HASHES.\n";
	}

	foreach my $coord (@relCoordiantes) {
		my $absStart;
		if($strand eq "+") {
			if($coord->{'Start'} <0) { ## if it is upstream
				$absStart = $absGeneStart + $coord->{'Start'};
			} else {
				$absStart = $absGeneEnd + $coord->{'Start'} -1;
			}
		} else {
			if($coord->{'Start'} <0) { ## if it is upstream
				$absStart = $absGeneEnd - $coord->{'Start'} - $coord->{'Length'} + 1 ;
			} else {
				$absStart = $absGeneStart - $coord->{'Start'} - $coord->{'Length'} + 1;
			}
		}
		### Now copy the hash
		my %absCoordiantesHash;
		foreach my $curKey (keys %$coord) {
			if($curKey eq "Start") { $absCoordiantesHash{'Start'} = $absStart; }
			elsif ($curKey eq "Length") { $absCoordiantesHash{'Length'} = $coord->{'Length'} -1; }
			else { $absCoordiantesHash{$curKey} = $coord->{$curKey}; }

		}
		push(@absCoordinates, \%absCoordiantesHash);
	}
	if(ref($CNSCoordinatesRef) eq "ARRAY") {
		return @absCoordinates;
	} elsif(ref($CNSCoordinatesRef) eq "HASH") {
		return %{ $relCoordiantes[0] };
	}
}
########################################################################
##### Shift coordiantes

sub shiftTargetCoordinate {
	my ($currentCoordiante, $strand, $targetSequence, $referenceSequence, $shiftBy, $length) = @_;

	my $untrimmedUpSequence = substr($referenceSequence, 0, $shiftBy);
	my $untrimmedDownSequence = substr($referenceSequence, $shiftBy + $length);
	my $hitUpstream = substr($targetSequence, 0, $shiftBy);
	my $hitDownstream = substr($targetSequence, $shiftBy + $length);

	my $upstreamInsertions = ($untrimmedUpSequence =~ tr/Z//) + 0;
	my $downstreamInsertions = ($untrimmedDownSequence =~ tr/Z//) + 0;
	my $upstreamGaps = ($hitUpstream =~ tr/-//) + 0;
	my $downstreamGaps = ($hitDownstream =~ tr/-//) + 0;

	my $newTargetPosition;

	if($strand eq "+") {
		$newTargetPosition = $currentCoordiante + $shiftBy + $upstreamInsertions - $upstreamGaps;
	} else {
		$newTargetPosition = $currentCoordiante - $shiftBy + $downstreamInsertions - $downstreamGaps;
	}
	return $newTargetPosition;	
}

1;
