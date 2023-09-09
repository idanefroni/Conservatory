#!/usr/bin/perl

package Mapping;

use POSIX;
use strict;
use warnings;
use List::Util qw(min max);
use lib './scripts';
use ConservatoryUtils;

sub new {
    my $class = shift;
    # we can initialize mapping from either a file, a string, a fasta alignment or parameters
    my ($CNSID,$species,$locus,$pos, $strand,$referenceRelativePosition,$seq, $refSeq, $absChr, $absPos, $absStrand, $name, );

    if(scalar @_ == 1) {
        my $firstparam= shift;
        ## if this is a file, just read the next line from it
        if($firstparam->isa('IO::Handle')) {
            $firstparam = <$firstparam>;
            chomp($firstparam);
        } elsif($firstparam->isa('Mapping')) { ### If its a mapping, copy
            $CNSID = $firstparam->getCNSID();
            $species = $firstparam->getSpecies();
            $locus = $firstparam->getLocus();
            $pos = $firstparam->getPos();
            $strand = $firstparam->getStrand();
            $referenceRelativePosition = $firstparam->getRRP();
            $seq = $firstparam->getSeq();
            $refSeq = $firstparam->getRefSeq();
            $absChr = $firstparam->getAbsChr();
            $absPos = $firstparam->getAbsPos();
            $absStrand = $firstparam->getAbsStrand();
            $name = $firstparam->getName();
        } elsif($firstparam->isa('Bio::LocatableSeq')) {          
			my ($targetLocus,$targetPosition,$targetStrand,$absChromosome, $geneStrand,$absPosition, $targetName) = split /:/, $firstparam->id();
            $species = geneToSpecies($targetLocus);
            $locus = $targetLocus;
            $pos = $targetPosition;
            $strand = $targetStrand;
            $absChr = $absChromosome;
            $absStrand = $geneStrand;
            $absPos = $absPosition;
            $name = $targetName;
            $CNSID="";
           
            $seq = uc($firstparam->seq);
            my ($leadingGaps) = $seq =~ /^(-+)/;
            if(!defined $leadingGaps) { $leadingGaps= ""; }
            $seq =~ s/^(-+)//;
            $seq =~ s/(-+)$//;
            $referenceRelativePosition = length($leadingGaps);
        } elsif(ref($firstparam) eq '' )  {
            ($CNSID,$species,$locus,$pos, $strand,$referenceRelativePosition,$seq, $refSeq,$absChr, $absPos, $absStrand, $name) = split /,/, $firstparam;
        } else {
            die "ERROR: Cannot instantiate mapping. Unknown input type " . ref($firstparam) . ".\n";
        }
    } else {
        ($CNSID,$species,$locus,$pos, $strand,$referenceRelativePosition,$seq, $refSeq, $absChr, $absPos, $absStrand, $name) = @_;
    }

    if(defined $seq) { 
        $seq = uc($seq);
    }
    if(!defined $absChr) { $absChr="";}
    if(!defined $absPos || $absPos eq "") { $absPos=0;}
    if(!defined $absStrand) { $absStrand="";}        
    if(!defined $name) { $name="";} 
    if(!defined $refSeq) { $refSeq=""; }

    my $self = {
	    _CNS => $CNSID,
        _Species => $species,
        _Locus => $locus,
        _Pos => $pos,
        _Strand => $strand,
        _RRP => $referenceRelativePosition,
        _Seq => $seq,
        _RefSeq => $refSeq,
        _AbsChr => $absChr,
        _AbsPos => $absPos,
        _AbsStrand => $absStrand,
        _Name => $name
    };

    bless $self, $class;
    return $self;
}

sub getLen {
    my ($self) = @_;
    return length($self->{_Seq});

}
sub getAbsLen {
    my ($self) = @_;
    return length($self->{_Seq}) - ($self->{_Seq} =~ tr/-// );
}
sub getAbsChr {
    my ($self) = @_;
    return $self->{_AbsChr};
}

sub setAbsChr {
    my ($self, $absChr) = @_;
    $self->{_AbsPos} = $absChr;
}

sub getAbsPos {
    my ($self) = @_;
    return $self->{_AbsPos};
}

sub setAbsPos {
    my ($self, $absPos) = @_;
    $self->{_AbsPos} = $absPos;
}

sub getPos {
    my ($self) = @_;
    return $self->{_Pos};
}

sub setPos {
    my ($self, $pos) = @_;
    $self->{_Pos} = $pos;
}

sub setRelPosFromAbs {
    my ($self, $absPos, $geneCoordinates) = @_;
    my $newRelativePosition;

    if($geneCoordinates->{'Strand'} eq "+") {
        if($self->getPos()>=0) {
            $newRelativePosition = $absPos - $geneCoordinates->{'End'}; 
        } else {
            $newRelativePosition = $absPos - $geneCoordinates->{'Start'}; 
        }
    } else {
        if($self->getPos()>=0) {
            $newRelativePosition = $geneCoordinates->{'Start'} - $absPos - $self->getAbsLen() + 1;
        } else {
            $newRelativePosition = $geneCoordinates->{'End'} - $absPos - $self->getAbsLen() + 1; 
        }
    }
    $self->{_Pos} = $newRelativePosition;
}

sub getAbsEnd {
    my ($self) = @_;
    return $self->getAbsPos() + $self->getAbsLen();
}

sub getAbsStrand {
    my ($self) = @_;
    return $self->{_AbsStrand};
}

sub setAbsStrand {
    my ($self, $strand) = @_;
    $self->{_AbsStrand} = $strand;
}


sub getStrandInGenome {
    my ($self) = @_;
    if($self->{_AbsStrand} eq "-" && $self->{_Strand} eq "+") {
        return "-";
    } elsif($self->{_AbsStrand} eq "-" && $self->{_Strand} eq "-") {
        return "+";
    } elsif($self->{_AbsStrand} eq "+" && $self->{_Strand} eq "-") {
        return "-";
    } else {
        return "+";
    }
}

sub identical {
    my ($self, $other) = @_;
    if($self->getAbsPos() == $other->getAbsPos() &&
       $self->getAbsChr() eq $other->getAbsChr() &&
       $self->getAbsLen() eq $other->getAbsLen() &&
       $self->getLocus() eq $other->getLocus()) {
        return 1;
    } else {
        return 0;
    }
}

sub getBreakPoint {
    my ($self) = @_;
    return $self->{_Breakpoint};
}

sub setBreakpoint {
    my ($self, $breakpoint) = @_;
    $self->{_Breakpoint} = $breakpoint;
}

sub getSeq {
    my ($self) = @_;
    return $self->{_Seq};
}
#######################################
sub getTargetSeqInGenome {
    my ($self) = @_;

    if($self->getStrandInGenome() eq "-") {
        return reverseComplement($self->getTargetSeq());
    } else {
        return $self->getTargetSeq();
    }
}
#######################################
sub getTargetSeq {
    my ($self) = @_;
    my $targetSeq = $self->{_Seq};
    $targetSeq =~ s/-//g;
    return $targetSeq;
}

### returns the sequence in CNS space.
sub getCNSSeq {
    my ($self) = @_;
    my $genomeSeq = $self->{_Seq};
    my $refSeq = $self->{_RefSeq};
    my $cnsSeq = $genomeSeq;

    if($refSeq ne "") {
        my @insertionInReference = findAll($refSeq, "-");
	   	foreach ( @insertionInReference ) { substr($cnsSeq, $_,1) = "Z"; }
		$cnsSeq =~ s/Z//g;   ## Remove the insertions from the target and reference sequences to keep the frame of reference as the reference sequence only
    }
    return $cnsSeq;
}

sub getRefCNSSeq {
    my ($self) = @_;
    my $refSeq = $self->{_RefSeq};

    $refSeq =~ s/-//g;
    return($refSeq);
}

sub getSeqBP {
    my ($self) = @_;
    my $seqNoGaps = $self->getSeq();
    $seqNoGaps =~ s/-//g;
    return length($seqNoGaps);
}

sub getSubsetSeq {
    my ($self, $start, $end) = @_;

    my $subMapping = $self->createSubSetMapping($start, $end);
    return $subMapping->getCNSSeq();
}

sub setSeq {
    my ($self, $newSeq) = @_;
    $self->{_Seq} = $newSeq;
}

sub getRefSeq {
    my ($self) = @_;
    return $self->{_RefSeq};
}

sub setRefSeq {
    my ($self, $newSeq) = @_;
    $self->{_RefSeq} = $newSeq;
}

sub getLocus {
    my ($self) = @_;
    return $self->{_Locus};
}
sub getGenome {
    my ($self) = @_;
    return geneToGenome($self->{_Locus});
}

sub getStrand {
    my ($self) = @_;
    return $self->{_Strand};
}

sub getCNSID {
    my ($self) = @_;
    return $self->{_CNS};
}
sub setCNSID {
    my ($self, $CNSID) = @_;
    if(ref($CNSID) eq 'CNS') { $CNSID = $CNSID->getID(); }
    $self->{_CNS} = $CNSID;
}
sub unlink {
    my ($self) = @_;
    $self->{_CNS}="CNSDELETED";
}

sub hasAbsCoordiantes {
    my ($self) =@_;
    return ($self->getAbsStrand() ne "");
}

sub getRRP {
    my ($self) = @_;
    return $self->{_RRP};
}
sub setRRP {
    my ($self, $RRP) = @_;
    $self->{_RRP} = $RRP;
}

sub getSpecies {
    my ($self) = @_;
    return $self->{_Species};
}
sub getName {
    my ($self) = @_;
    return $self->{_Name};
}

sub setName {
    my ($self, $name) = @_;
    $self->{_Name} = $name;
}

sub flip {
    my ($self) = @_;
    $self->{_Seq} = reverseComplement($self->{_Seq});
    if($self->getStrand() eq "+") {
        $self->setPos( $self->getPos() + $self->getAbsLen());
    } else {
        $self->setPos( $self->getPos() - $self->getAbsLen());
    }
    $self->{_Strand} = flipStrand($self->{_Strand});

}

sub flipKeepStrand {
    my ($self) = @_;
    $self->{_Seq} = reverseComplement($self->{_Seq});
}

##################################################################
####
#### Is there an overlap between this mapping 
####
#### Can accept a second mapping (and margin) and then the function will check absolute overlap in the genome
#### Alternative, can accept start and end coordiantes and then it will check relative coordiante
####

sub overlap {
    my ($self, $other, $minOverlap) = @_;
    my $startOne, my $endOne, my $startTwo, my $endTwo;

    if(!defined $minOverlap) { $minOverlap=0; }    

    if($other->isa("Mapping")) {
   	    $startOne = $self->getAbsPos();
        $endOne = $startOne + $self->getAbsLen();
        $startTwo = $other->getAbsPos();
        $endTwo = $startTwo + $other->getAbsLen();

	    if ( $self->getAbsChr() ne $other->getAbsChr() ||
            geneToGenome($self->getLocus()) ne geneToGenome($other->getLocus())) {
                return 0;
        }
    } elsif($other->isa("CNS")) {
   	    $startOne = $self->getRRP();
        $endOne = $startOne + $self->getLen();
        $startTwo = $other->getPos();
        $endTwo = $startTwo + $other->getLen();
        if($startTwo>$endTwo) {
            my $tmp = $endTwo;
            $endTwo = $startTwo;
            $startTwo=$tmp;
        }
        if($startOne>$endOne) {
            my $tmp = $endOne;
            $endOne = $startOne;
            $startOne=$tmp;
        }
    } elsif((scalar @_) ==3) {
   	    $startOne = $self->getRRP();
        $endOne = $self->getRRP() + $self->getLen();
        $startTwo = $other;
        $endTwo = $minOverlap;
        $minOverlap = 0;
    } else {
        die "ERROR: Bad parameters to Mapping::overlap.\n";
    }

    if (($startOne <= $startTwo && $endOne <= $startTwo) || ($startOne >= $endTwo) || ($startTwo >= $endOne)) {
		return 0;
	} else {
        my $overlapStart = ($startOne > $startTwo) ? $startOne : $startTwo;
        my $overlapEnd = ($endOne < $endTwo) ? $endOne : $endTwo;
        my $overlap = $overlapEnd - $overlapStart;
        if($overlap < $minOverlap) {
            return 0;
        } else {
            return min(1,$overlap/ ($endTwo-$startTwo));
        }
	}

}

#############################################################################
### Returns the percent identity between the target and reference
### Usage: identity(Sequence1, Sequence2)
###
### Returns a identity (0-100)
###

sub identity {
	my ($self) = @_;
    my $one = $self->getCNSSeq();
    my $two = $self->getRefCNSSeq();

	my ($match, $total) = (0,0);
    if(length($one) != length($two)) {
        $self->print();
	    die "ERROR: In Identity: two sequences must be the same length ($one:$two).\n" unless length($one) == length($two);
    }
	for my $i (0 .. length($one) -1 ){
		my ($chr1,$chr2) = (substr($one,$i,1), substr($two,$i,1));
		next if $chr1 eq '-' || $chr2 eq '-';
		$match += $chr1 eq $chr2;
		$total++;
	}
	return ($total>0) ? ($match/$total) *100:0;

}
#### Remove leading and trailing spaces
sub trim {
    my ($self) = @_;
    my $seq = $self->getSeq();
    $seq =~ s/-*$//;
    my ($leadingGaps) = $seq =~ /^(-+)/;
    if(!defined $leadingGaps) { $leadingGaps= ""; }
    $seq =~ s/^(-+)//;
    $self->transpose(length($leadingGaps));
    $self->{_Seq} = $seq;
}
### strip gaps
sub stripGaps {
    my ($self) = @_;
    my $seq = $self->getSeq();
    $seq =~ s/-//g;
    $self->setSeq($seq);
}

### Move a mapping a given number of bases

sub transpose {
    my ($self, $transposition) = @_;

    if($self->getStrandInGenome() eq "+") {
        $self->setAbsPos( $self->getAbsPos()+$transposition);
    } else {
        $self->setAbsPos( $self->getAbsPos() - $transposition);
    }

    $self->setPos( $self->getPos()+$transposition);
}

###########################################################################################
##### creates a new mapping object with a sequence subset
##### update the coordinates accordingly
sub createSubSetMapping {
    my ($self, $start, $end) = @_;

    my $seqPadding = "-------------------";

    my $startInCNSSpace = max(0,$start - $self->getRRP() );
    my $endInCNSSpace = min( $self->getLen() - 1 ,$end - $self->getRRP() );


    (my $startInTargetSpace, my $endInTargetSpace) = $self->_translateCNSSpaceToTargetSpace($startInCNSSpace, $endInCNSSpace);
    (my $startInTargetGenomeSpace, my $endInTargetGenomeSpace) = $self->_translateCNSSpaceToGenomeSpace($startInCNSSpace, $endInCNSSpace);
    (my $startInReferenceSpace, my $endInReferenceSpace) = $self->_translateCNSSpaceToReferenceSpace($startInCNSSpace, $endInCNSSpace);

    my $subMapping = new Mapping($self);
    my $subSeq = substr(($self->{_Seq} . $seqPadding ), $startInTargetSpace, $endInTargetSpace - $startInTargetSpace);
    $subMapping->setSeq($subSeq);

    if(defined $self->getRefSeq()) {
        my $subRefSeq =substr($self->{_RefSeq} . $seqPadding, $startInReferenceSpace, $endInReferenceSpace - $startInReferenceSpace);
        $subMapping->setRefSeq($subRefSeq);
    }

    $subMapping->setRRP(max($start, $self->getRRP()));
    $subMapping->setPos($subMapping->getPos() + $startInTargetGenomeSpace);

    if($subMapping->hasAbsCoordiantes()) {
        if($subMapping->getAbsStrand eq "+") {
            $subMapping->setAbsPos($subMapping->getAbsPos() + $startInTargetGenomeSpace);
        } else {
            $subMapping->setAbsPos($subMapping->getAbsPos() + $self->getAbsLen() - $endInTargetGenomeSpace );
        }
    }
    return $subMapping;
}

### Takes a coordinate in CNS space and translate to target genome space
sub _translateCNSSpaceToTargetSpace {
    my ($self, $start, $end) = @_;
    if(!defined $start || !defined $end) { die "ERROR: _translateCNSSpaceToTargetSpace must accept start and end coordiantes.\n"; }
    my $startGaps = _countGaps($self->getRefSeq(), $start);
    my $endGaps = _countGaps($self->getRefSeq(), $end); 
    my $targetSpaceStart = $start + $startGaps;
    my $targetSpaceEnd = $end + $endGaps; 
    return ($start, $end);
}

sub _translateCNSSpaceToGenomeSpace { 
    my ($self, $start, $end) = @_;

    ($start,$end) = $self->_translateCNSSpaceToTargetSpace($start,$end);
    my $seqToStart;
    my $seqToEnd;
    if($self->getStrand() eq "+") {
        $seqToStart = substr($self->getSeq(), 0,$start);    
        $seqToEnd = substr($self->getSeq(), 0,$end);
    } else {
        $seqToStart = substr($self->getSeq(), $end);
        $seqToEnd = substr($self->getSeq(), $start);
    }
    my $startGaps = $seqToStart =~ tr/-//;
    my $endGaps = $seqToEnd =~ tr/-//;
    $start -= $startGaps;
    $end -= $endGaps; 

    if($self->getStrand() eq "+") {
        return ($start, $end );
    } else {
        return ($self->getLen() - $end, $self->getLen() - $start);
    }
}

sub _countGaps {
    (my $seq, my $bases) = @_;
    my $gaps=0, my $length=0;
    for my $curbp (split //, $seq) {
        last if $length==$bases;
        if($curbp eq "-") {
            $gaps++;
        } else {
            $length++;
        }
    }
    return $gaps;
}
sub _translateCNSSpaceToReferenceSpace {
    my ($self, $start, $end) = @_;
 #   my $startGaps = _countGaps($self->getRefSeq(), $start);
  #  my $endGaps = _countGaps($self->getRefSeq(), $end); 
   # my $refSpaceStart = $start + $startGaps;
   # my $refSpaceEnd = $end + $endGaps;    
    return $self->_translateCNSSpaceToTargetSpace($start,$end);
}

### Merge two mappings. must be from the same genome.
sub merge {
    my ($self, $other) = @_;
    my $maxBias=3;

    if(geneToGenome($self->getLocus()) ne geneToGenome($other->getLocus()) ) {
        die "ERROR: Trying to merge two mappings from different genomes: " . $self->getLocus() . " and " . $other->getLocus() . "\n";
    }

   if($self->getStrand() ne $other->getStrand()) {
        $other->flip();
   }

    my $mySeq = $self->getTargetSeq();
    my $otherSeq = $other->getTargetSeq();

    my $combinedSequence = '-' x (max( $self->getAbsEnd(), $other->getAbsEnd()) - min($self->getAbsPos(), $other->getAbsPos())+ $maxBias);
    my $combinedSeqStart = min($self->getAbsPos(), $other->getAbsPos());
    if(length($combinedSequence)> 10000) { ## Sanity check
        die "ERROR: merging sequences longer than 10kb:" . $self->getLocus() . " and " . $other->getLocus() . "\n";
    }

    my $bias=0;
    my $rejectMerge=0;
    # Calculate overlap length
    my $minEnd = ($self->getAbsPos + length($mySeq) < $other->getAbsPos + length($otherSeq)) ? $self->getAbsPos + length($mySeq) :  $other->getAbsPos + length($otherSeq);
    my $overlapLen = ($minEnd >= $self->getAbsPos() && $minEnd >= $other->getAbsPos())
                ? $minEnd - max($self->getAbsPos(), $other->getAbsPos())
                : 0;
    if($self->getStrand() eq "+") {
       	substr($combinedSequence,$self->getAbsPos() - $combinedSeqStart, length($mySeq)) = $mySeq;

        if(substr($combinedSequence,$other->getAbsPos() - $combinedSeqStart, $overlapLen) ne substr($otherSeq,0,$overlapLen)) {
            $bias = _findMergeBias($combinedSequence, $other->getAbsPos() - $combinedSeqStart, substr($otherSeq,0,$overlapLen), $maxBias);
        }           
        if(!defined $bias) {
            $rejectMerge=1;
        } else {
    	    substr($combinedSequence,$other->getAbsPos() - $combinedSeqStart + $bias, length($otherSeq)) = $otherSeq;
        }
    } else { 
        $mySeq = reverseComplement($mySeq);
        $otherSeq = reverseComplement($otherSeq);
       	substr($combinedSequence,$self->getAbsPos() - $combinedSeqStart, length($mySeq)) = $mySeq;

        if(substr($combinedSequence,$other->getAbsPos() - $combinedSeqStart, $overlapLen) ne substr($otherSeq,0,$overlapLen)) {
            $bias = _findMergeBias($combinedSequence, $other->getAbsPos() - $combinedSeqStart, substr($otherSeq,0,$overlapLen), $maxBias);
        }           
        if(!defined $bias) {
            $rejectMerge=1;
        } else {
    	    substr($combinedSequence,$other->getAbsPos() - $combinedSeqStart+$bias, length($otherSeq)) = $otherSeq;
            $combinedSequence = reverseComplement($combinedSequence);
        }
    }

    if(!$rejectMerge) {

        $self->setAbsPos($combinedSeqStart);
        $self->setSeq($combinedSequence);
        $self->setPos(min($self->getPos(), $other->getPos));
        return 1;
    } else {
        return 0;
    }
}

sub _findMergeBias {
    my ($combined, $start, $searchSeq, $maxBias) = @_;
    my $bias= 1;
    while( substr($combined,$start + $bias, length($searchSeq)) ne $searchSeq && $bias <= $maxBias ) {
        $bias++;
    }
    if(substr($combined,$start + $bias, length($searchSeq)) eq $searchSeq) {
        return $bias;
    } else {
        return undef;
    }
}

### Get absolute coordiantes for locus, if we don't have them
sub fillAbsoluteCoordiantes {
    my ($self, $genomeDB) = @_;

	my %geneCoordinates = getGeneCoordinates($genomeDB->getConservatoryDir() , $genomeDB->genomeToFamily( $self->getGenome() ),  $self->getGenome(), $self->getLocus());
	my %relativeCoordinate = ('Start' => $self->{_Pos}, 'Len' => $self->getAbsLen());

	my %absGeneCoord = translateRealtiveToAbsoluteCoordinates( \%relativeCoordinate, $geneCoordinates{'Start'}, $geneCoordinates{'End'}, $geneCoordinates{'Strand'}, $self->getStrand() );
	$self->{_AbsChr}  = $geneCoordinates{'Chr'};
	$self->{_AbsStrand} = $geneCoordinates{'Strand'};
	$self->{_AbsPos} = $absGeneCoord{'Start'};
}

sub isAlive {
    my ($self) = @_;
    if($self->{_CNS} eq "CNSDELETED") {
        return 0;
    } else {
        return 1;
    }
}
### Is this mapping to a plastid gene
sub isPlastid {
    my ($self) = @_;

    if($self->{_Species} eq "Athaliana" || $self->{_Species} eq "MtrunA17r5" || $self->{_Species} eq "pgrisea") {
            my $geneName = $self->{_Locus};
            $geneName =~ s/.*-//;
            if (substr($geneName,0,3) eq "ATC" || substr($geneName,0,3) eq "ATM" ||
                substr($geneName,0,9) eq "MtrunA17CP" || substr($geneName,0,9) eq "MtrunA17M" ||
                substr($geneName,0,7) eq "PhygriC" || substr($geneName,0,7) eq "PhygriM") {
                return 1;
            }
    }
    return 0;
}

sub print {
    my ($self, $outFile) = @_;

    if(!defined $outFile) { $outFile = \*STDOUT; } 
	print $outFile join(",",
					$self->{_CNS},
					$self->{_Species},
					$self->{_Locus},
					$self->{_Pos},
					$self->{_Strand},
					$self->{_RRP},                                        
					$self->{_Seq},
					$self->{_RefSeq},                    
					$self->{_AbsChr},
					$self->{_AbsPos},
					$self->{_AbsStrand},
					$self->{_Name} ) . "\n";
}

sub printFasta {
    my ($self, $outFile) = @_;
    if(!defined $outFile) { $outFile = \*STDOUT; } 

	print $outFile ">" . join(":", ($self->{_Locus},
    							    $self->{_Pos},
								    $self->{_Strand},
								    $self->{_AbsChr},
								    $self->{_AbsStrand},
								    $self->{_AbsPos},
								    $self->{_Name}) ) . "\n" . 
								    $self->getCNSSeq() . "\n";
}

1;

