#!/usr/bin/perl

package ConservatoryTree;

use POSIX;
use strict;
use warnings;
use lib './scripts';
use ConservatoryUtils;
use Bio::TreeIO;
use Bio::SeqIO;
use Bio::AlignIO;

sub new {
    my ($class, $genomeDB, $treeFileName, $annotatedTreeFileName, $treeNodeAgesFileName) = @_;

    if(!defined $genomeDB) { die "ERROR: no genome database provided.\n"; }
    # we can initialize a CNS database from file.
    ## set up defaults, if not specified
    if(!defined $treeFileName) { $treeFileName = $genomeDB->getConservatoryDir() . "/genomes/Conservatory.tree"; }
    if(!defined $annotatedTreeFileName) { $annotatedTreeFileName = $genomeDB->getConservatoryDir() . "/genomes/Conservatory.Annotated.tree"; }
    if(!defined $treeNodeAgesFileName) { $treeNodeAgesFileName = $genomeDB->getConservatoryDir() . "/genomes/Conservatory.Annotated.NodeAges.csv"; }

    if(! -e $treeFileName) { print "WARNING: Cannot find unannotated conservatory tree ($treeFileName).\n"; }
    if(! -e $annotatedTreeFileName) { print "WARNING: Cannot find annotated conservatory tree ($annotatedTreeFileName).\n"; }
    if(! -e $treeNodeAgesFileName) { print "WARNING: Cannot find node ages files for annotated conservatory tree ($treeNodeAgesFileName).\n"; }

    my $conservatoryAnnotatedTreeio = Bio::TreeIO->new(-file => $annotatedTreeFileName, -format => 'newick');
    my $conservatoryAnnotatedTree = $conservatoryAnnotatedTreeio->next_tree;
    if(!defined $conservatoryAnnotatedTree) { die "ERROR: Cannot read annotated tree file $annotatedTreeFileName. Bad format.\n"; }
    my %speciesInTree = map { ($_)->id => 1} $conservatoryAnnotatedTree->get_leaf_nodes();

    my $conservatoryUnannotatedTreeio = Bio::TreeIO->new(-file => $treeFileName, -format => 'newick');
    my $conservatoryUnannotatedTree = $conservatoryUnannotatedTreeio->next_tree;
    if(!defined $conservatoryUnannotatedTreeio) { die "ERROR: Cannot read unannotated tree file $treeFileName. Bad format.\n"; }

    my %conservatoryTreeNodeAges;
    open(my $nodeAgesFile, $treeNodeAgesFileName) || die ("ERROR: Error reading $treeNodeAgesFileName.\n");
    while(<$nodeAgesFile>) {
	    chomp;
	    my ($node,$age)= split /,/;
	    $conservatoryTreeNodeAges{$node}=$age;
    }
    close($nodeAgesFile);

    my $self = {
        _fastML => $genomeDB->getConservatoryDir() . "/scripts/fastml",
        _genomeDB => $genomeDB,
        _UnannotatedTree => $conservatoryUnannotatedTree,
        _AnnotatedTree => $conservatoryAnnotatedTree,
        _NodeAges => \%conservatoryTreeNodeAges,
        _Species => \%speciesInTree
    };

    bless $self, $class;
    return $self;
}

##########################################################################################################
####
####  Given a tree and a set of nodes on the tree, findDeepestCommonNode returns the name of the deepest node
####   shared by all leaves. This can be used to identify the evolutionary origin of a CNS
###

sub findDeepestCommonNode {
	my ($self, $leavesListRef) = @_;
	my @leavesList = @$leavesListRef;
	if(scalar @leavesList ==0) {
		print "WARNING: No species provided for findDeepestCommonNode.\n"; 
		return "NA";
	}
	
	my $anchorLeafNode = $self->{_AnnotatedTree}->find_node( ($leavesList[0] =~ s/\./_/gr) );
	if(!defined $anchorLeafNode) { 
		#die "ERROR in findDeepestCommonNode: Cannot find deepest node $leavesList[0] in tree.\n";
		return "";
	}
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
		my $node = $self->{_AnnotatedTree}->find_node( -id => $leaf=~ s/\./_/gr );
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

sub getChildren {
	my ($self, $nodeName) = @_;
	my @children=();
	my $treeNode = $self->{_AnnotatedTree}->find_node( $nodeName );
	if(defined $treeNode) { 
		foreach my $curChildNode  ($treeNode->each_Descendent()) {
			push @children, $curChildNode->id;
		}
	}	
	return @children;
}

##########################################################################################################
###############
##### Find the birth of the CNS (one node),
sub getCNSBirth {
	my ($self, $leavesListRef) = @_;
	return $self->findDeepestCommonNode($leavesListRef);
}

##########################################################################################################
###############
##### Find the all deaths of a CNS, given a certain birth point
#####  i.e. all the nodes where a CNS is lost (but not in a single species)
sub getCNSDeaths {
	my ($self, $birthNodeID, $leavesListRef, $deathList) = @_;
	my $thisNode = $self->getNodeByID($birthNodeID);

	if(!defined $thisNode) {
		return 0;
	}
	if($thisNode->is_Leaf) { 
		if(grep /$birthNodeID/, @$leavesListRef ) {
			return 1;
		} else {
			return 0;
		}
	}
	my %childrenHaveCNS;
	my $totalCNSs=0;

	foreach my $curChild ( $thisNode->each_Descendent()) {
		$childrenHaveCNS{$curChild->id} = $self->getCNSDeaths($curChild->id, $leavesListRef, $deathList);
		$totalCNSs += $childrenHaveCNS{$curChild->id};
	}
	foreach my $curChild ( $thisNode->each_Descendent()) {
		my $numOfDescendents = scalar $curChild->get_all_Descendents;

		if($totalCNSs>0 && $numOfDescendents>=0 && $childrenHaveCNS{$curChild->id} == 0) { ## If this child has no CNS, BUT is not a leaf, add a death point
			push @{ $deathList },  $curChild->id;
		}
		### also add the living CNSs
	}

	return $totalCNSs;
}
##########################################################################################################
############
############ Given a CNS, reconstuct the ancesteral sequence for the CNS
###########    Accepts CNS Name (for the temporary files) and CNS sequences with the species being the key

sub getReconstructedSequence {
    my ($self, $CNS, $sequencesHashRef, $minInfoContent) = @_;
    my %sequencesHash = %$sequencesHashRef;
    if(ref($CNS) eq 'CNS') { $CNS = $CNS->getID(); }

	## Set up temporary file names (for MAFFT and fastml)
    my $tmpDir = $self->getTmpDir();

	my $tmpOutputFastaFileName = $tmpDir . "/tmp_sp_$CNS.fasta";
	my $tmpOutputFastaAlignedFileName = $tmpDir . "/tmp_sp_$CNS.align.fasta";
	my $tmpOutputTrimmedTreeFileName = $tmpDir . "/tmp_tree_$CNS.tree";
	my $tmpOutputAncesteralSequenceFileName = $tmpDir . "/tmp_as_$CNS.fasta";

	open(my $tmpOutputFastaFile , ">" , $tmpOutputFastaFileName);

	### Generate fasta file for sequence reconstruction using tree and trim the tree
	foreach my $curSpecies (keys %sequencesHash) {
		print $tmpOutputFastaFile ">$curSpecies\n" . $sequencesHash{$curSpecies} . "\n";
	}
	close($tmpOutputFastaFile);

	## Align using MAFFT
	system("mafft --quiet --auto --thread -1 --ep 0.3 --op 7 $tmpOutputFastaFileName > $tmpOutputFastaAlignedFileName");
#	system("mafft --quiet --auto --thread -1 --auto $tmpOutputFastaFileName > $tmpOutputFastaAlignedFileName");

	## and perform reconstruction

	my $trimmedTree = $self->{_UnannotatedTree}->clone();
	my @trimToTheseLeaves = keys %sequencesHash;
	trimTreeToLeaves($trimmedTree, \@trimToTheseLeaves);

	my $trimmedTreeOutput = Bio::TreeIO->new(-file => ">$tmpOutputTrimmedTreeFileName", -format => 'newick');
	$trimmedTreeOutput->write_tree($trimmedTree);
			
	system( $self->{_fastML} . " -b -qf -mh -s $tmpOutputFastaAlignedFileName -t $tmpOutputTrimmedTreeFileName -y /dev/null/a -d /dev/null/a -x /dev/null/a -k /dev/null/a -d /dev/null/a -e /dev/null/a -j $tmpOutputAncesteralSequenceFileName > /dev/null");
	unlink("log.txt"); ### remove the fastml leftover file
	### Read concensus sequence
	my $inReconstructionFasta = Bio::SeqIO->new(-file => $tmpOutputAncesteralSequenceFileName, -format => "fasta", -alphabet => "dna");
	my $reconstructedSequence;
	while(my $seq = $inReconstructionFasta->next_seq) {
		if($seq->id eq "N1") { $reconstructedSequence = $seq->seq; }
	}

	### Filter the sequence based on information content. When little sequence coverage is available, fastml usually just guesses the nucleotide.
	my $inAlignedFasta = Bio::AlignIO->new(-file => $tmpOutputFastaAlignedFileName,
							  		 	   -format => "fasta");
	my $inAligned = $inAlignedFasta->next_aln;
	my @gappiness = getGappiness($inAligned);

    for my $pos (0..(scalar @gappiness -1 )) {
	    if ($gappiness[$pos] > (100-$minInfoContent)) {
       		substr($reconstructedSequence, $pos, 1) = "-";
   		}
	}

	unlink($tmpOutputAncesteralSequenceFileName);
	unlink($tmpOutputFastaFileName);
	unlink($tmpOutputFastaAlignedFileName);
	unlink($tmpOutputTrimmedTreeFileName);
	return($reconstructedSequence);
}

sub getAgeForNode {
    my ($self, $node) = @_;
    if(defined $self->{_NodeAges}->{$node}) {
        return $self->{_NodeAges}->{$node};
    } else {
        print "WARNING: Cannot find node $node in the conservatory tree.\n";		
		return -1;
    }
}

sub getNodeByID {
	my ($self, $nodeID) = @_;
	return $self->{_AnnotatedTree}->find_node( -id => $nodeID=~ s/\./_/gr );
}
#################################################################
### check if a certain node exists in the tree
sub exists {
    my ($self, $node) = @_;
    if(defined $self->{_AnnotatedTree}->find_node($node) ) {
        return 1;
    } else {
        return 0;
    }
}

sub getTmpDir {
    my ($self) = @_;
    return $self->{_genomeDB}->getTemporaryDir();
}

sub writeTree {
    my ($self, $outputFileName, $annotate) = @_;
    if(defined $annotate && $annotate!=0) { 
        ### Annotated Internal nodes
        my $curNodeNumber=1;
        my @nodes = $self->{_AnnotatedTree}->get_nodes();
        foreach my $curNode (@nodes) {
            if($curNode->id eq "") { 
            	$curNode->id("N$curNodeNumber");
            	$curNodeNumber++;
            }
        }
    }
    my $conservatoryTreeOutput = Bio::TreeIO->new(-file => ">$outputFileName", -format => 'newick');
    $conservatoryTreeOutput->write_tree( $self->{_AnnotatedTree});
}

1;
