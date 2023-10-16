package ConservatoryUtils;

use POSIX;
use strict;
use warnings;
use List::Util qw(min max sum);
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
our $minIdentityToKeepBreakpoint=50;
our $minSpeciesToKeepBreakpoint	= 5;
our $minSequenceContentInAlignment = 0.5;
our $minCNSLength=8;
our $minSpeciesForCNS=4;
our $minCNSCoverageAfterSplit=0; 
our $superCNSPrefix="Super";
our $minSequenceContentToConsiderReconstruction = 0.33; ## how many of the sequences in the alignment has to be informative (non-gap)
													  # for the nucleotide to be considered for ancesteral seq reconstruction. 
													  # this is to avoid reconstructing sequences based on small number of samples.

#### Sequence Filtering and parameters
our $homoPolymerFilterLength = 8;         ### Filter homopolymer of this length
our $homoPolymerFuzzyFilterLength = 12;   ### Filter homopolymer with one variant of this length
our $atRichRegionFilterLength = 10;       ### Filter atrich regions of this length
our $atRichRegionFilterFlankLength = 2;   ###   when filtering atrich, leave flanks of this length
our $dimerPolymerFilterLength = 5;        ### Filter repeating dimers of these many repeats
my $filterCharacter = 'X';


my %footprintDatabase;   ### Database of footprints for the getGeneCoordiantes function

our @ISA= qw( Exporter );
our @EXPORT = qw( overlapFragment reverseComplement flipStrand alignPairwise isGeneName geneToSpecies geneToGenome geneToLocus fullNameToShortName lengthWithoutGaps dropAsterixFromProtein findAll getRandomORFLength getLongestORF getCNSbreakpoints polishCNSAlignments 
				  getGappiness shiftTargetCoordinate translateRealtiveToAbsoluteCoordinates getGeneCoordinates clearGeneCooordinateDatabase multipleAlignment cleanChrName trimTreeToLeaves
				  $minCNSLength $minCNSCoverageAfterSplit  $minSequenceContentInAlignment $minSpeciesForCNS $superCNSPrefix $standardDeviationsToSplit $minSpeciesToSplitCNS $maxSpeciesToInitiateCNSSplit $minSequenceContentToConsiderReconstruction
				  $homoPolymerFilterLength $homoPolymerFuzzyFilterLength $atRichRegionFilterLength $atRichRegionFilterFlankLength $dimerPolymerFilterLength
				  generateLastzCommandLine extractFastaToFile getMAFQuality fuzzyFilterHomopolymers fuzzyFilterDimers fuzzyFilterATRich mergeFastaSequences);


######################################################################

sub generateLastzCommandLine {
	my ($genomeDB, $laxAlignment) = @_;
	my $minIdentity, my $minHSPthreshold;

	if($laxAlignment) {
		$minIdentity = $genomeDB->getMinDeepIdentity();
		$minHSPthreshold = $genomeDB->getMinDeepThreshold();
	} else {
		$minIdentity = $genomeDB->getMinFamilyIdentity();
		$minHSPthreshold = $genomeDB->getMinFamilyThreshold();
	}
	 return "lastz --format=maf --gap=200,100 --nochain --noytrim --seed=match4 --gapped --strand=both --step=1 --ambiguous=iupac --identity=$minIdentity --ydrop=1000 --hspthreshold=$minHSPthreshold --gappedthresh=$minHSPthreshold";
}

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
### overlapFragment(startOne, endOne, StartTwo, endTwo)
### returns relative overlap (0-1).
###

sub overlapFragment {
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


###################################################################################
####
#### Returns the overall quality score for a MAF alignment
####  accepts a MAF filename
####
####   returns total alignmetn quality
####

sub getMAFQuality {
	my ($MAFFileName) = @_;

	open (my $m, $MAFFileName);
	my $maf = Bio::AlignIO->new(-fh => $m, -format => 'maf');
	my $quality =0;
	while(my $aln = $maf->next_aln()) {
		$quality = $quality + $aln->score;
	}
	close($m);
	return $quality;
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
				$absStart = $absGeneEnd - $start - $coord->{'Len'} +1;
			} else {
				$absStart = $absGeneStart - $start - $coord->{'Len'} +1;
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

####################################################################################################################################
#### Extract sequence for locus from a compressed fasta file, filter for low-information sequence and output to a seperate fasta file
####    and return the filtered sequence

sub extractFastaToFile {
	(my $fastaFile, my $locus, my $outputFastaFile) = @_;

	open (my $inFasta, "samtools faidx $fastaFile $locus |");
	my $fastaHeader = <$inFasta>;
	chomp($fastaHeader);
	my $seq;
	my @seqLines;
	while(my $curLine = <$inFasta>) {
		chomp($curLine);
		push @seqLines, uc($curLine);
	}
	close($inFasta);
	$seq = join("", @seqLines);

	my $origSeq = $seq;
	### Now mask homopolymers and AT rich regions
	$seq = fuzzyFilterHomopolymers($seq,$homoPolymerFilterLength,0);  
	$seq = fuzzyFilterATRich($seq,$atRichRegionFilterLength,0, $atRichRegionFilterFlankLength); 
	$seq = fuzzyFilterDimers($seq, $dimerPolymerFilterLength,0);

	### and repeating common trimers
	$seq =~ s/([C|G|A]TT){4}/XXXXXXXXXXXX/g;
	$seq =~ s/([C|G|A]AA){4}/XXXXXXXXXXXX/g;
	$seq =~ s/([C|G|A]AT){4}/XXXXXXXXXXXX/g;

	#### and repeating tetramers
	$seq =~ s/(TCTA){4}/XXXXXXXXXXXXXXXX/g;
	### and ambigious sequences
	$seq =~ s/N/$filterCharacter/g;

	## and write to file
	open(my $outFasta, ">$outputFastaFile");
	print $outFasta "$fastaHeader\n$seq\n";
	close($outFasta);

	return ($seq);
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


################################################################################
###### utility filtering functions #######################################################

##### Fuzzy filter homopolymers

### fuzziness is the number of mismatch positions allowed (either 0 or 1)
sub fuzzyFilterHomopolymers {
	my ($inSeq, $polymerLength, $fuzziness) = @_;
	my @inSeqArr = split //, $inSeq;
	my @polymers = ('A' x $polymerLength,
					'T' x $polymerLength,
					'C' x $polymerLength,
					'G' x $polymerLength);

	my $fuzzySearchString="(?=($polymers[1]";
	my @replacementArr = ($filterCharacter) x $polymerLength;

	for my $curPolymer (@polymers) {
		if($fuzziness == 1) {
			for my $curNuc (1..$polymerLength) {
				my $curSearch = $curPolymer;
				substr($curSearch,$curNuc-1,1) = ".";
				$fuzzySearchString .= "|$curSearch"
			}
		} else {
			$fuzzySearchString .= "|$curPolymer";
		}
	}
	$fuzzySearchString .= "))";
	
	my @polymerOccurences;
	$_ = $inSeq;

	while(/$fuzzySearchString/g) {
		push @polymerOccurences, pos();
	}
	foreach my $curPos (@polymerOccurences) {
		@inSeqArr[$curPos..($curPos+($polymerLength-1))] = @replacementArr;
	}

	return join("", @inSeqArr);
}

##### Fuzzy filter for repeating dimers sequences
## Fuzziness is 0 or 1 (number of mismatches accepted)
sub fuzzyFilterDimers {
	my ($inSeq, $polymerLength, $fuzziness) = @_;
	my @inSeqArr = split //, $inSeq;
	my @replacementArr = ($filterCharacter) x ($polymerLength * 2);
	my @dimerPolymers;

	for my $nuca ( ('A','T','C','G')) {
		for my $nucb ( ('A','T','C','G')) {
			if($nuca ne $nucb) {
				push @dimerPolymers, ("$nuca$nucb" x $polymerLength);
			}
		}
	}

	my $fuzzySearchString="(?=($dimerPolymers[1]";

	for my $curPolymer (@dimerPolymers) {
		if($fuzziness==1) {
			for my $curNuc (1..($polymerLength*2-1)) {
				my $curSearch = $curPolymer;
				substr($curSearch,$curNuc-1,1) = ".";
				$fuzzySearchString .= "|$curSearch"
			}
		} else {
				$fuzzySearchString .= "|$curPolymer"
		}
	}
	$fuzzySearchString .= "))";

	## Now collected all the places where there is a match (accepting overlapping regions)
	my @dimerOccurences;
	$_ = $inSeq;

	while(/$fuzzySearchString/g) {
		push @dimerOccurences, pos();
	}

	## Now mask Seq
	foreach my $curPos (@dimerOccurences) {
		@inSeqArr[$curPos..($curPos+($polymerLength*2-1))] = @replacementArr;
	}

	return join("", @inSeqArr);
}

########################################################################################
### Fuzzy filter AT rich regions
###

sub fuzzyFilterATRich {
	my ($inSeq, $ATrichLength, $fuzziness, $flanklength) = @_;
	my @replacementArr = ($filterCharacter) x ($ATrichLength - $flanklength*2);

	my $translatedInSeq = $inSeq;
	$translatedInSeq =~ tr/ATWCGRYSKMBDHVNX/0001111111111111/;

	my @inSeqArr = split //, $inSeq;
	my @inSeqTrArr = split //, $translatedInSeq;
	my $seqLength = (scalar @inSeqArr)-$ATrichLength-1;
	my @ATrichPos;

	## Look for AT rich positions
	foreach my $curPos (0..$seqLength) {
		if(sum(@inSeqTrArr[$curPos..($curPos+$ATrichLength-1)])<=$fuzziness) {
			push @ATrichPos, $curPos;
		}
	}
	## Now mask Seq
	foreach my $curPos (@ATrichPos) {
		@inSeqArr[($curPos+$flanklength)..($curPos+$ATrichLength-1-$flanklength)] = @replacementArr;
	}

	return join("", @inSeqArr);

	return($inSeq);
}

sub mergeFastaSequences {
	my ($originalSequence, $sequenceToAdd, $startCoordinate) = @_;

	$sequenceToAdd =~ s/-/N/g;

	my $mergedSequence = $originalSequence;

	# If we are not overriding anything, just place the sequence in the right place

	if(substr($originalSequence, $startCoordinate, length($sequenceToAdd)) eq ('N' x length($sequenceToAdd) ) ) {
		substr($mergedSequence, $startCoordinate, length($sequenceToAdd)) = $sequenceToAdd;
	} elsif( substr($originalSequence, $startCoordinate, length($sequenceToAdd)) ne $sequenceToAdd ) { ### if there is already sequence there (that is different than the one we want to add), merge the sequences
		for my $curNucleotidePos ($startCoordinate .. ($startCoordinate + length($sequenceToAdd)-1)) {
			if(substr($mergedSequence, $curNucleotidePos, 1) eq "N") { ## If we have no data, put the merged one
				substr($mergedSequence, $curNucleotidePos,1) = substr($sequenceToAdd, $curNucleotidePos - $startCoordinate,1);
			}
		}
	}
	return $mergedSequence;
}

1;


