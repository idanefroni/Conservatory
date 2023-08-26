package ConservatoryUtils;

use POSIX;
use strict;
use warnings;
use List::Util qw(min max);
use File::Temp qw /tempfile/;

use Statistics::Basic qw(:all nofill);
use Bio::SimpleAlign;
use Algorithm::NeedlemanWunsch;

use Exporter;



##############################################################################
##### CNS Splitting Parameters

our $standardDeviationsToSplit = 3; ## The number of standard deviation in number of aligned species to trigger a CNS split
our $minSpeciesToSplitCNS=5;  ### Minimal number of species to consider a split of the CNS
our $maxSpeciesToInitiateCNSSplit = 200; ## Do not split CNS if it is supported by atleast this number of species
									 ## To avoid spliting of very highly conserved CNSs
our $minIdentityToKeepBreakpoint	=50;
our $minSpeciesToKeepBreakpoint	= 5;
our $minSequenceContentInAlignment = 0.5;
our $minCNSLength=8;
our $minSpeciesForCNS=4;
our $minCNSCoverageAfterSplit=0.3; 
our $superCNSPrefix="Super";
my %footprintDatabase;   ### Database of footprints for the getGeneCoordiantes function

our @ISA= qw( Exporter );
our @EXPORT = qw (overlap reverseComplement flipStrand alignPairwise isGeneName geneToSpecies geneToGenome geneToLocus fullNameToShortName lengthWithoutGaps dropAsterixFromProtein findAll getRandomORFLength getLongestORF getCNSbreakpoints polishCNSAlignments 
				  getGappiness shiftTargetCoordinate translateRealtiveToAbsoluteCoordinates getGeneCoordinates clearGeneCooordinateDatabase findDeepestCommonNode multipleAlignment cleanChrName trimTreeToLeaves
				  $minCNSLength $minCNSCoverageAfterSplit  $minSequenceContentInAlignment $minSpeciesForCNS $superCNSPrefix $standardDeviationsToSplit $minSpeciesToSplitCNS $maxSpeciesToInitiateCNSSplit);

#############################################################################
### Returns the reverse complement of a DNA sequence


sub reverseComplement {
	my ($seq) = @_;
	$seq =~ tr/ATCGatcg/TAGCtagc/;
	$seq = reverse $seq;
	return $seq;
}


#############################################################################
### Perform pairwise alignment between two DNA sequences
### Usage: alignPariwise(Sequence1, Sequence2)
###
### Returns a SimpleAlign object with the alignment
###

sub alignPairwise {
	my ($seq1, $seq2) = @_;

	sub score_sub {
		  if (!@_) {  return -0.5;  }  ## gap penalty
  			## mismatch scores -1, match +1
  			return ($_[0] eq $_[1]) ? 1 : -1;
	}
	
	my $nw = Algorithm::NeedlemanWunsch->new(\&score_sub);
	$nw->gap_open_penalty(-5);
	$nw->gap_extend_penalty(-0.25);

	my $seq1Arr = [split//, $seq1];
	my $seq2Arr = [split//, $seq2];
	my (@align1, @align2);

	$nw->align($seq1Arr,$seq2Arr,
	   {
			align   => sub {unshift @align1, $seq1Arr->[shift]; unshift @align2, $seq2Arr->[shift];},
			shift_a => sub {unshift @align1, $seq1Arr->[shift]; unshift @align2,            '-'},
			shift_b => sub {unshift @align1,               '-'; unshift @align2, $seq2Arr->[shift]},
		});

	my $pairwise = Bio::SimpleAlign->new();
	$pairwise->add_seq(Bio::LocatableSeq->new( -seq => join("",@align1), -id=> "seq1", -start=>0, -end=> scalar @align1));
	$pairwise->add_seq(Bio::LocatableSeq->new( -seq => join("",@align2), -id=> "seq2", -start=>0, -end=> scalar @align2));

	return $pairwise;
}

#############################################################################
### Compute overlap between two fragments. 
### Parameters:
### overlap(startOne, endOne, StartTwo, endTwo)
### returns relative overlap (0-1).
###

sub overlap {
	my ($startOne, $endOne, $startTwo, $endTwo) =@_;
	if ( ($startOne <= $startTwo && $endOne <= $startTwo) ||
		 ($startOne >= $endTwo) || ($startTwo >= $endOne)) {
		return 0;
	} else {
		return (min($endOne, $endTwo) - max($startOne, $startTwo))/ max( $endOne-$startOne, $endTwo-$startTwo);
	}
}


##################################################################################
# Flip Strand (turn from - to plus and vice versa)

sub flipStrand {
	my ($strand) = @_;
	if($strand eq "+") {
		return "-";
	} else {
		return "+";
	}
}

##################################################################################
# Returns species or genome name from the name of a gene

sub isGeneName { 
	my ($geneNameToCheck) = @_;
	if($geneNameToCheck =~ /^.+-.+-.+$/) {
		return 1;
	} else {
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

sub geneToLocus { ### accepts full gene name. returns just the locus
	(my $geneName) = @_;
	my @geneNameComponents = split '-', $geneName;
	return $geneNameComponents[2];
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
### Returns a chromosome name that we can work with (without embelshiment and making sure its not a number)

sub cleanChrName {
	my $chr = shift(@_);
	chomp($chr);
	if(substr($chr,0,1) eq ">") {
		$chr = substr($chr,1);
	}
	if(index($chr, " ")!= -1 || index($chr, "\t") != -1) {
		$chr = (split(' ', $chr))[0];
	}
	return "C" . $chr;
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
		my $start = $curHit->{'RRP'};
		my $end = $start + $curHit->{'Len'};

		foreach my $pos ($start..$end) {
			my $curNucleotide = substr($curHit->{'Seq'}, ($pos-$start),1);
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
			print "DEBUG: Break at $pos because " . $CNSdelta[$pos] . ". from " . $CNSCoverageSmooth[$pos-1] . " to " . $CNSCoverageSmooth[$pos] . "to" . $CNSCoverageSmooth[$pos+1]  . "\n";
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
		my $start =  $curHit->{'RRP'};
		my $end = $start + length( $curHit->{'Seq'});

		for my $curBreakpoint (0..(scalar @breakpoints - 2) ) {
			### if the hit overlaps the sub CNS

#      print "DEBUG: polishing. $curBreakpoint  (" . $breakpoints[$curBreakpoint] . "-" . $breakpoints[$curBreakpoint+1] .") ($start - $end).\n";
			if( overlap($breakpoints[$curBreakpoint], $breakpoints[$curBreakpoint+1], $start, $end)) {

				### distance from breakpoint start
				my $positionInCNS = max($breakpoints[$curBreakpoint], $curHit->{'RRP'});
				my $positionInSubCNS = max(0, $curHit->{'RRP'}- $breakpoints[$curBreakpoint]);
				my $subCNSSeqStart = max(0,$breakpoints[$curBreakpoint] - $start);
				my $subCNSSeqLength = min(length($curHit->{'Seq'}), ($breakpoints[$curBreakpoint+1]-$breakpoints[$curBreakpoint])-$positionInSubCNS );
				my $subCNSTargetSeq = substr($curHit->{'Seq'}, $subCNSSeqStart , $subCNSSeqLength);

				my $referenceSequence = $curHit->{'Seq'}; 
				my $newReferenceSeq = $subCNSTargetSeq;  ## Default value for reference sequence is the current sequence
				if(defined $curHit->{'RefSeq'}) {
					$referenceSequence = $curHit->{'RefSeq'};
					$newReferenceSeq = substr($curHit->{'RefSeq'}, $subCNSSeqStart , $subCNSSeqLength);
				}

				### remove trailing and leading gaps
				my $numOfTrailingGaps = $subCNSTargetSeq =~ s/-+$//;
       
				my $numOfLeadingGaps = $subCNSTargetSeq =~ s/^-+//;
				$positionInSubCNS += $numOfLeadingGaps;

				my $numOfInternalGaps = () = ($subCNSTargetSeq =~ /-/g);

				### Only include the hit of the coverage is sufficient and if its not too gappy
				if((length($subCNSTargetSeq)-$numOfInternalGaps) / ($breakpoints[$curBreakpoint+1] - $breakpoints[$curBreakpoint] ) > $minCNSConservationAfterSplit && length($subCNSTargetSeq) >= $minCNSLength && ($numOfInternalGaps / length($subCNSTargetSeq)) < $minSequenceContentInAlignment ) {

					## update the target coordinates
					my $newTargetPosition= shiftTargetCoordinate($curHit->{'Pos'}, $curHit->{'Strand'}, $curHit->{'Seq'}, $referenceSequence, $subCNSSeqStart+$numOfLeadingGaps, length($subCNSTargetSeq));

					my $newAbsolutePosition=0;
					if(defined $curHit->{'AbsPos'}) {
						$newAbsolutePosition = shiftTargetCoordinate($curHit->{'AbsPos'}, $curHit->{'GeneStrand'}, $curHit->{'Seq'}, $referenceSequence, $subCNSSeqStart+$numOfLeadingGaps, length($subCNSTargetSeq));
					}
					# remove internal gaps
					$subCNSTargetSeq =~ s/-//g;
    
					### Log the alignment to the breakpointlist
					push (@splitAndPolishedAlignments, {
							'Species' => $curHit->{'Species'},
							'Locus' => $curHit->{'Locus'},
							'Strand' => $curHit->{'Strand'},
							'Len' => length($subCNSTargetSeq),
							'RefUpDown' => $curHit->{'RefUpDown'},
							'Seq' => $subCNSTargetSeq,
							'RefSeq' => $newReferenceSeq,
							'RRP' => $positionInSubCNS+$numOfLeadingGaps,
							'Pos' => $newTargetPosition,
							'AbsChr' => $curHit->{'AbsChr'},
							'AbsPos' => $newAbsolutePosition,
							'GeneStrand' => $curHit->{'GeneStrand'},
							'Name' => $curHit->{'Name'},
							'Breakpoint' => $curBreakpoint
					});

					my $seq = Bio::LocatableSeq->new(-seq => $subCNSTargetSeq,
													 -start => 1,
													 -end => length($subCNSTargetSeq),
												     -id => $curHit->{'Locus'} . $newTargetPosition);
					$alignmentsForBreakpoints{$curBreakpoint}->add_seq($seq);
				}
			}
		}
	}

	my %breakpointsToDelete;
	for my $curBreakpoint (0..(scalar @breakpoints - 2) ) {

		my $alignedBP = multipleAlignment($alignmentsForBreakpoints{$curBreakpoint});
   
     if(!defined $alignedBP) {
			$breakpointsToDelete{$curBreakpoint}=1;
    } elsif( $alignedBP->percentage_identity() < $minIdentityToKeepBreakpoint || $alignedBP->num_sequences < $minSpeciesToKeepBreakpoint ) {
			$breakpointsToDelete{$curBreakpoint}=1;
		} else {
			## update sequences for BP
			foreach my $curHit (@splitAndPolishedAlignments) {
				if($curHit->{'Breakpoint'} == $curBreakpoint) {
					my $alignedSeq = $alignedBP->get_seq_by_id($curHit->{'Locus'} . $curHit->{'Pos'});
					if(!defined $alignedSeq) { print "ERR: Can't find ID " . $curHit->{'Locus'} . $curHit->{'Pos'} . "\n"; }
					my $newAlignedSeq = uc($alignedSeq->seq());
					### remove leading and trailing gaps
					$newAlignedSeq =~ s/-+$//; # remove trailing
					my $numOfLeadingGaps = $newAlignedSeq =~ s/^-+//;  ## and leading
					$curHit->{'Seq'} = $newAlignedSeq;
					$curHit->{'Len'} = length($newAlignedSeq =~ s/-//gr);
					## and update the positions
					$curHit->{'RRP'} += $numOfLeadingGaps;
					$curHit->{'Pos'} = shiftTargetCoordinate($curHit->{'Pos'}, $curHit->{'Strand'}, $curHit->{'Seq'}, $curHit->{'RefSeq'}, $numOfLeadingGaps, length($curHit->{'Seq'}));
					if(defined $curHit->{'AbsPos'}) {
						$curHit->{'AbsPos'} = shiftTargetCoordinate($curHit->{'AbsPos'}, $curHit->{'Strand'}, $curHit->{'Seq'}, $curHit->{'RefSeq'}, $numOfLeadingGaps, length($curHit->{'Seq'}));
					}					
				}
			}
		}
	}

  print "deleting breakpoints:" . join(",", keys %breakpointsToDelete) . "\n";
  
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
			$footprintDatabase{$genome}{$footprintGeneName}{'Chr'} = $array[0];
			$footprintDatabase{$genome}{$footprintGeneName}{'Start'} = $array[3];
			$footprintDatabase{$genome}{$footprintGeneName}{'End'} = $array[4];			
		}
		close($footprintFile);
	}
	if(! defined $footprintDatabase{$genome}{$gene} ) { die "ERROR: Cannot find gene $gene in genome $genome.\n"; }
	return %{ $footprintDatabase{$genome}{$gene} };
}

###################################################################################
#### clearGeneCooordinateDatabase
#####
#### The gene coordinate (footprint) cache can get quite big. This function removes it from memory.
####
sub clearGeneCooordinateDatabase {
	undef %footprintDatabase;
}

#######################################
##### Translate a coordiante list of segments to from relative to absolute coordinates.
##### Accepts a hash with "Start" value having the coordiante and "Length" having the length of the segment or an array of such hashes.

sub translateRealtiveToAbsoluteCoordinates {
	my ($CNSCoordinatesRef, $absGeneStart, $absGeneEnd, $geneStrand, $segementStrand) = @_;
	my @absCoordinates;
	my @relCoordiantes;
	if(ref($CNSCoordinatesRef) eq "ARRAY") {
		@relCoordiantes = @$CNSCoordinatesRef;
	} elsif(ref($CNSCoordinatesRef) eq "HASH") {
		@relCoordiantes = ($CNSCoordinatesRef);
	} else {
		die "ERROR: translateRelativeToAbsoluteCoordiante accepts either a HASH or an array of HASHES.\n";
	}

	foreach my $coord (@relCoordiantes) {
		my $absStart;
		my $start;
		if(defined $coord->{'Pos'}) { $start = $coord->{'Pos'} }
		elsif(defined $coord->{'Start'}) { $start = $coord->{'Start'}} 
		else {
			die "ERROR: No start or pos key in coordinate hash.\n";
		}

		if($geneStrand eq "+") {
			if($start <0) { ## if it is upstream
				$absStart = $absGeneStart + $start;
			} else {
				$absStart = $absGeneEnd + $start;
			}
		} else {
			if($start <0) { ## if it is upstream
				$absStart = $absGeneEnd - $start - $coord->{'Len'} + 1 ;
			} else {
				$absStart = $absGeneStart - $start - $coord->{'Len'};
			}
		}
		
		if($segementStrand eq "-") {
			if($geneStrand eq "+") {
				$absStart -= $coord->{'Len'};
			} else {
				$absStart += $coord->{'Len'};
			}
		}
		### Now copy the hash
		my %absCoordiantesHash;
		foreach my $curKey (keys %$coord) {
			if($curKey eq "Start") { $absCoordiantesHash{'Start'} = $absStart; }
			elsif ($curKey eq "Len") { $absCoordiantesHash{'Len'} = $coord->{'Len'} -1; }
			else { $absCoordiantesHash{$curKey} = $coord->{$curKey}; }

		}
		push(@absCoordinates, \%absCoordiantesHash);
	}
	if(ref($CNSCoordinatesRef) eq "ARRAY") {
		return @absCoordinates;
	} elsif(ref($CNSCoordinatesRef) eq "HASH") {
		return %{ $absCoordinates[0] };
	}
}
########################################################################
##### Shift coordiantes

sub shiftTargetCoordinate {
	my ($currentCoordinate, $strand, $targetSequence, $referenceSequence, $shiftBy, $length) = @_;

	my $untrimmedUpSequence = substr($referenceSequence, 0, $shiftBy);
	my $hitUpstream = substr($targetSequence, 0, $shiftBy);

	my $upstreamInsertions = ($untrimmedUpSequence =~ tr/Z//) + 0;
	my $upstreamGaps = ($hitUpstream =~ tr/-//) + 0;

	my $newTargetPosition;

	if($strand eq "+") {
		$newTargetPosition = $currentCoordinate + $shiftBy + $upstreamInsertions - $upstreamGaps;
	} else {
		$newTargetPosition = $currentCoordinate - $shiftBy - $upstreamInsertions + $upstreamGaps;
	}

	return $newTargetPosition;	
}

##########################################################################################################
####
####  Given a tree and a set of nodes on the tree, findDeepestCommonNode returns the name of the deepest node
####   shared by all leaves. This can be used to identify the evolutionary origin of a CNS
###

sub findDeepestCommonNode {
	my ($tree, $leavesListRef) = @_;
	my @leavesList = @$leavesListRef;

	my $anchorLeafNode = $tree->find_node( ($leavesList[0] =~ s/\./_/gr) );
	if(!defined $anchorLeafNode) { die "ERROR in findDeepestCommonNode: Cannot find deepest node $leavesList[0] in tree.\n"; }
	my %anchorLeafPath;
	my $node = $anchorLeafNode;
	my $deepestNode = $node->id;
	$anchorLeafPath{$deepestNode}=0;

	my $deepness=1;
	while($node->ancestor) {
		$node = $node->ancestor;
		$anchorLeafPath{$node->id} = $deepness++;
	}

	foreach my $leaf (@leavesList) {
		my $node = $tree->find_node( -id => $leaf=~ s/\./_/gr );
		if(! defined $node) { next;}
		while($node->ancestor) {
			$node = $node->ancestor;
			if(defined $anchorLeafPath{$node->id}) {  ### This is where the path intersect. Check if this is the deepest we found
				if( $anchorLeafPath{$node->id} > $anchorLeafPath{$deepestNode}) {
					$deepestNode = $node->id;
				}
				last;
			}
		}
	}
	return $deepestNode;
}

##########################################################################################################
####
####  multipleAlignment(alignObject)
####
####  Accepts a SimpleAlign object, performs multiple alignment and returns an aligned SimpleAlign
####

sub multipleAlignment {
	my ($inAlign) = @_;
	## make temporary fasta file
	my $pid = $$;
	chomp($pid);
	my $tempFileName = "/dev/shm/$pid.CNS.fasta";

	my $out = Bio::SeqIO->new(-file=> ">$tempFileName", -format => "fasta");
	foreach my $seq ($inAlign->each_seq) { $out->write_seq($seq); }

	## perform the alignment
	system("mafft --quiet --thread -1 --auto --ep 0.3 --op 7 $tempFileName > $tempFileName.align.fasta");

	## Now read the fasta file
	my $alignedCNSFile = Bio::AlignIO->new(-file => "$tempFileName.align.fasta",
										-format => "fasta");
	my $outAlign = $alignedCNSFile->next_aln;
	unlink($tempFileName);
	unlink("$tempFileName.align.fasta");
	return $outAlign;
}


################################################################################################################
#####
#####  trim tree - given a tree and a set of leaves, trim the tree to include just these leaves. 
#####    trim also internal nodes, or nodes with only one child
####
####
sub trimTreeToLeaves {
	my ($tree, $leavesListRef, $startNode) = @_;
	my @leavesList = @$leavesListRef;
	my %leavesMap = map { $_ => 1 } @leavesList;

	if(!defined $startNode) { $startNode = $tree->get_root_node(); }

	if(($startNode == $tree->get_root_node()) && $startNode->is_Leaf()) {
		die "ERROR: From trimTreeToLeaves: Tree is empty.\n";
	}
	## If this is a leaf and not found in our list of leaves to keep, delete
	if($startNode->is_Leaf()) {
		if(!defined $leavesMap{$startNode->id()} ) {
			$tree->remove_Node($startNode);
		} 
	} else {
		### trim the children

		foreach my $curChildNode ($startNode->each_Descendent()) {
			trimTreeToLeaves($tree, $leavesListRef, $curChildNode);
		}

		### If we have no children leaf, delete ourselves. If we have just one, delete it and assume it's id
		my @childrenLeft = $startNode->each_Descendent();

		if(scalar @childrenLeft == 0) {
			if(! ($startNode == $tree->get_root_node() )  ) {
				$tree->remove_Node($startNode);
			} else {
				print "WARNING: From trimmTreeToLeaves: we have deleted all nodes on the tree.\n";
			}

		} elsif(scalar @childrenLeft == 1 && $childrenLeft[0]->is_Leaf()) {

			$startNode->id( $childrenLeft[0]->id());
			$startNode->branch_length( $childrenLeft[0]->branch_length() );
			$tree->remove_Node( $childrenLeft[0]);
		} elsif(scalar @childrenLeft ==1) {
			if($startNode == $tree->get_root_node()) {
				$tree->set_root_node($childrenLeft[0]);
			} else {
				$startNode->ancestor()->add_Descendent($childrenLeft[0]);
				$tree->remove_Node( $startNode);
			}
		}
	}
}

1;


