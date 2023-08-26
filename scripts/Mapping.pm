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
    my ($CNSID,$species,$locus,$pos, $strand,$referenceRelativePosition,$len, $seq, $absChr, $absPos, $absStrand, $name, $refSeq);

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
            ($CNSID,$species,$locus,$pos, $strand,$referenceRelativePosition,$len, $seq, $absChr, $absPos, $absStrand, $name) = split /,/, $firstparam;
        } else {
            die "ERROR: Cannot instantiate mapping. Unknown input type " . ref($firstparam) . ".\n";
        }
    } else {
        ($CNSID,$species,$locus,$pos, $strand,$referenceRelativePosition,$len, $seq, $absChr, $absPos, $absStrand, $name) = @_;
    }

    if(defined $seq) { 
        $seq = uc($seq);
    }
    if(!defined $name) { $name="";}

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
        if($self->{_Pos}>=0) {
            $newRelativePosition = $absPos - $geneCoordinates->{'End'} -1; 
        } else {
            $newRelativePosition = $absPos - $geneCoordinates->{'Start'}; 
        }
    } else {
        if($self->{_Pos}>=0) {
            $newRelativePosition = $geneCoordinates->{'Start'} - $absPos - $self->getAbsLen() + 1;
        } else {
            $newRelativePosition = $geneCoordinates->{'End'} - $absPos - $self->getAbsLen() + 1; 
        }
    }
    $self->{_Pos} = $newRelativePosition;
}

sub getAbsEnd {
    my ($self) = @_;
    return $self->{_AbsPos} + $self->getAbsLen();
}

sub getAbsStrand {
    my ($self) = @_;
    return $self->{_AbsStrand};
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
       $self->getAbsLen() eq $other->getAbsLen()) {
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

sub getSubsetSeq {
    my ($self, $start, $end) = @_;
    my $startInSeq = max(0,$start - $self->getRRP() );
    my $endInSeq = min( $self->getLen(),$end - $self->getRRP() ); 
    return substr($self->getSeq(), $startInSeq, $endInSeq - $startInSeq);
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
    my ($self, $other, $margin) = @_;

    if($other->isa("Mapping")) {
        if(!defined $margin) { $margin=0; }

   	    my $startOne = $self->getAbsPos()- $margin;
        my $endOne = $startOne + $self->getLen() + $margin;
        my $startTwo = $other->getAbsPos() - $margin;
        my $endTwo = $startTwo + $other->getLen() + $margin;

	    if ( $self->getAbsChr() ne $other->getAbsChr() ||
            geneToGenome($self->getLocus()) ne geneToGenome($other->getLocus()) ||
            ($startOne <= $startTwo && $endOne <= $startTwo) ||
		    ($startOne >= $endTwo) || ($startTwo >= $endOne)) {
		    return 0;
	    } else {
		    return (min($endOne, $endTwo) - max($startOne, $startTwo))/ max( $endOne-$startOne, $endTwo-$startTwo);
	    }
    } else {
   	    my $startOne = $self->getRRP();
        my $endOne = $self->getRRP() + $self->getLen();
        my $startTwo = $other;
        my $endTwo = $margin;
        if ( ($startOne <= $startTwo && $endOne <= $startTwo) ||
		    ($startOne >= $endTwo) || ($startTwo >= $endOne)) {
		    return 0;
	    } else {
		    return (min($endOne, $endTwo) - max($startOne, $startTwo))/ max( $endOne-$startOne, $endTwo-$startTwo);
	    }
    }
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
    $self->{_Seq} = $seq
}
### strip gaps
sub stripInternalGaps {
    my ($self) = @_;
    my $seq = $self->getSeq();
#    $seq =~ s/(?<=^)-|-(?=$)//g;
    $seq =~ s/-//g;
    $self->setSeq($seq);
}
sub stripGaps {
    my ($self) = @_;
    $self->trim();
    $self->stripInternalGaps();
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

    my $startInSeq; 
    my $endInSeq;
    $startInSeq = max(0,$start - $self->getRRP() );
    $endInSeq = min( $self->getLen() ,$end - $self->getRRP() );

    my $newRRP = max($start, $self->getRRP());
    my $transpositionInCNSSpace = $newRRP - $self->getRRP(); ### movement in CNS space

    my $subSeq;
    my $subRefSeq;
    if($self->getStrand() eq "+") {
        $subSeq = substr($self->getSeq(), $startInSeq, $endInSeq - $startInSeq);
    } else {
        $subSeq = reverseComplement(substr(reverseComplement($self->getSeq()), $startInSeq, $endInSeq - $startInSeq) );
    }

    if(defined $self->getRefSeq()) {
        $subRefSeq = substr($self->getRefSeq(), $startInSeq, $endInSeq - $startInSeq);
    }

    my $upstreamSeq = substr($self->getSeq(), 0, $startInSeq);
    my $transpositionInGenomeSpace = $transpositionInCNSSpace - (($upstreamSeq =~ tr/-//) + 0);

    ### now create the sub mapping
    my $subMapping = new Mapping($self);
    $subMapping->setSeq($subSeq);
    if(defined $self->getRefSeq()) {
        $subMapping->setRefSeq($subRefSeq);
    }
    $subMapping->setRRP($newRRP);
    $subMapping->transpose($transpositionInGenomeSpace);

    return $subMapping;
}


### Merge two mappings. must be from the same genome.
sub merge {
    my ($self, $other) = @_;
    if(geneToGenome($self->getLocus()) ne geneToGenome($other->getLocus()) ) {
        die "ERROR: Trying to merge two mappings from different genomes: " . $self->getLocus() . " and " . $other->getLocus() . "\n";
    }

   if($self->getStrandInGenome() ne $other->getStrandInGenome()) {
        $other->flip();
   }

    my $mySeq = $self->getSeq();
    my $otherSeq = $other->getSeq();

    # remove gaps
    $mySeq =~ s/-//g;
    $otherSeq =~ s/-//g;

    my $combinedSequence = '-' x (max( $self->getAbsEnd(), $other->getAbsEnd()) - min($self->getAbsPos(), $other->getAbsPos()));
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
    if($self->getStrandInGenome() eq "+") {
       	substr($combinedSequence,$self->getAbsPos() - $combinedSeqStart, length($mySeq)) = $mySeq;
        if(substr($combinedSequence,$other->getAbsPos() - $combinedSeqStart, $overlapLen) ne substr($otherSeq,0,$overlapLen)) { $bias=1;}

        if($bias ==1 && substr($combinedSequence,$other->getAbsPos() - $combinedSeqStart + $bias, $overlapLen) ne substr($otherSeq,0,$overlapLen)) { 
            $rejectMerge = 1;
        }

	    substr($combinedSequence,$other->getAbsPos() - $combinedSeqStart + $bias, length($otherSeq)) = $otherSeq;
    } else { 
        $mySeq = reverseComplement($mySeq);
        $otherSeq = reverseComplement($otherSeq);
       	substr($combinedSequence,$self->getAbsPos() - $combinedSeqStart, length($mySeq)) = $mySeq;

        if(substr($combinedSequence,$other->getAbsPos() - $combinedSeqStart, $overlapLen) ne substr($otherSeq,0,$overlapLen)) { $bias= 1;}
        if($bias ==1 && substr($combinedSequence,$other->getAbsPos() - $combinedSeqStart + $bias, $overlapLen) ne substr($otherSeq,0,$overlapLen)) { 
            $rejectMerge = 1;
        }
       
	    substr($combinedSequence,$other->getAbsPos() - $combinedSeqStart+$bias, length($otherSeq)) = $otherSeq;
        $combinedSequence = reverseComplement($combinedSequence);
    }

    if(!$rejectMerge) {
        my $transposed = $self->getAbsPos() - $combinedSeqStart;
        $self->setAbsPos($combinedSeqStart);
        $self->setSeq($combinedSequence);
        $self->setPos(min($self->getPos(), $other->getPos));
        return 1;
    } else {
        return 0;
    }

}

### Get absolute coordiantes for locus, if we don't have them
sub findAbsoluteCoordiantes {
    my ($self, $genomeDB) = @_;
    if($self->{_AbsChr} eq "") {
		my %geneCoordinates = getGeneCoordinates($genomeDB->getConservatoryDir() , $genomeDB->speciesToFamily( $self->{_Species} ),  geneToGenome( $self->{_Locus} ), $self->{_Locus});

		my %relativeCoordinate = ('Start' => $self->{_Pos}, $self->getLen());

		my %absGeneCoord = translateRealtiveToAbsoluteCoordinates( \%relativeCoordinate, $geneCoordinates{'Start'}, $geneCoordinates{'End'}, $geneCoordinates{'Strand'}, $self->{_Strand} );
		$self->{_AbsChr}  = $geneCoordinates{'Chr'};
		$self->{_AbsStrand} = $geneCoordinates{'Strand'};
		$self->{_AbsPos} = $absGeneCoord{'Start'};
	}
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
					$self->getLen(),
					$self->{_Seq},
					$self->{_AbsChr},
					$self->{_AbsPos},
					$self->{_AbsStrand},
					$self->{_Name} ) . "\n";
}

sub printFasta {
    my ($self, $outFile) = @_;
    if(!defined $outFile) { $outFile = \*STDOUT; } 
    my $seq= $self->{_Seq};
    if($self->{_Strand} eq "-") { $seq = reverseComplement($seq); }

	print $outFile ">" . join(":", ($self->{_Locus},
    							    $self->{_Pos},
								    $self->{_Strand},
								    $self->{_AbsChr},
								    $self->{_AbsStrand},
								    $self->{_AbsPos},
								    $self->{_Name}) ) . "\n" . 
								    $seq . "\n";
}

1;

