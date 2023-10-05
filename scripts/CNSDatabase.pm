#!/usr/bin/perl

package CNSDatabase;

use POSIX;
use strict;
use warnings;
use lib './scripts';
use CNS;

sub new {
    my ($class, $CNSDBFileName, $verbose) = @_;
    my %CNSTable =();
    if(!defined $verbose) { $verbose =0; }
    # we can initialize a CNS database from file.

    if(defined $CNSDBFileName && $CNSDBFileName ne "") {
        if(! -e $CNSDBFileName) { die "ERROR: Cannot find CNS database file $CNSDBFileName.\n"; }
        ### now read the CNS DB
        my $CNSLineCount= (qx(wc -l $CNSDBFileName))[0];
        chomp($CNSLineCount);
        $CNSLineCount =~ s/ .*//;
        my $curPos=1;
        open(my $CNSInputFile, "<", $CNSDBFileName) || die "ERROR: Cannot open CNS database file $CNSDBFileName.\n";
        while(<$CNSInputFile>) {
            chomp;
	        if($verbose) { print "PROGRESS: Reading CNS..." . ($curPos++) . "/$CNSLineCount.\r" };
            if(defined $_ && $_ ne "") {
                my $newCNS = new CNS($_);
            
                if(defined $newCNS) {
                    $CNSTable{ $newCNS->getID() } = $newCNS;
                }
            }
        }
        if($verbose) { print "\n"; };
    }

    my @CNSOrder = values %CNSTable;
    my $self = {
        _CNSTable => \%CNSTable,
        _CNSOrder => \@CNSOrder,
        _Ordered => 0
    };

    bless $self, $class;
    return $self;
}

sub createSubsetCNSDatabase {
    my ($self, $subset) = @_;
    my @CNSsubset = @$subset;

    my $subsetDatabase = new CNSDatabase();
    my %newCNSTable;
    foreach my $CNSToCopy (@CNSsubset) {
        $newCNSTable{ $CNSToCopy->getID() } = $CNSToCopy;
    }
    $subsetDatabase->{_CNSTable} = \%newCNSTable;
    $subsetDatabase->{_CNSOrder} = values %{ $subsetDatabase->{_CNSTable} };

    return $subsetDatabase;
}

sub getCNSByOrder {
    my ($self) = @_;
    if(! $self->{_Ordered}) {
        my @CNSOrder = values %{ $self->{_CNSTable} };
        $self->{_CNSOrder} = \@CNSOrder;
    }
    return $self->{_CNSOrder};
}
sub getCNSByID {
    my ($self,$CNSID) = @_;
    if(!defined $self->{_CNSTable}->{$CNSID} ) { die "ERROR: $CNSID not found in CNS database.\n"; }
    return $self->{_CNSTable}->{$CNSID};
}
sub getCNSIDs {
    my ($self) = @_;
    return keys %{ $self->{_CNSTable} };
}

sub findReferenceMappings {
    my ($self, $genomeDB, $mappingDB) = @_;
    foreach my $curCNS ( values %{ $self->{_CNSTable}}) {
        if(!$curCNS->findReferenceMapping($genomeDB, $mappingDB)) {
            $self->deleteCNS($curCNS);
            $mappingDB->deleteCNS($curCNS);
        }
    }
}

sub verify {
    my ($self, $mappingDB, $conservatoryTree) = @_;
    foreach my $curCNS ( values %{ $self->{_CNSTable}}) {
        if(!$conservatoryTree->exists($curCNS->getConservationLevel()) ) {
            $self->deleteCNS($curCNS);
            $mappingDB->deleteCNS($curCNS);
        }
    }
}

sub updateNumberOfSupportingSpecies {
    my ($self, $mappingDB) = @_;
    foreach my $curCNS (values %{ $self->{_CNSTable} }) {
        $curCNS->getSupportingSpecies($mappingDB);
    }
}

sub orderCNSByReference {
    my ($self) = @_;
    my @CNSOrderDeref = @{ $self->getCNSByOrder() };

    @CNSOrderDeref = sort { $a->getRefGenome() cmp $b->getRefGenome() ||
				  $a->getReferenceMapping()->getAbsChr() cmp $b->getReferenceMapping()->getAbsChr() ||
				  $a->getReferenceMapping()->getAbsPos()  <=> $b->getReferenceMapping()->getAbsPos() } @CNSOrderDeref;
    
    $self->{_CNSOrder} = \@CNSOrderDeref;
    $self->{_Ordered} = 1;
}

sub orderCNSByName {
    my ($self) = @_;
    my @CNSOrderDeref = @{ $self->getCNSByOrder() };

    @CNSOrderDeref = sort { $a->getID() cmp $b->getID()} @CNSOrderDeref;
    
    $self->{_CNSOrder} = \@CNSOrderDeref;
    $self->{_Ordered} = 1;
}

sub exists {
    my ($self, $CNSID) = @_;
    if(defined $self->{_CNSTable}->{$CNSID}) { 
        return 1;
    } else {
        return 0;
    }
}
sub deleteCNS {
    my ($self, $CNSID) = @_;
    ### accepts a CNS object or a CNSID string
    if(ref($CNSID) eq 'CNS') { $CNSID = $CNSID->getID(); }
    $self->{_CNSTable}->{$CNSID}->kill();
    delete $self->{_CNSTable}->{$CNSID};
}

sub renameCNS {
    my ($self, $CNS, $newCNSName) = @_;
    $self->{_CNSTable}->{ $newCNSName } = $CNS;
    delete $self->{_CNSTable}->{ $CNS->getID() };

    $CNS->setID($newCNSName);
    $self->{_Ordered} = 0;
}
sub add {
    my ($self, $CNS) = @_;
    $self->{_CNSTable}->{ $CNS->getID() } = $CNS;
    $self->{_Ordered} = 0;
}

sub getNumberOfCNSs {
    my ($self) = @_;
    return scalar keys %{ $self->{_CNSTable} };
}

sub writeDatabase {
    my ($self, $outputFileName) = @_;
    my $outputFile;

    if(defined $outputFileName) {
        open($outputFile, ">", $outputFileName);
    } else {
        $outputFile = \*STDOUT;
    }
    foreach my $curCNS (values %{ $self->{_CNSTable} }) {
        if($curCNS->isAlive()) {
            $curCNS->print($outputFile);
        }
    }
    if(defined $outputFileName) {
        close($outputFile);
    }
}

1;
