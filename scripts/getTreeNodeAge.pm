#!/usr/bin/perl

use strict;

my $treeFileName = $ARGV[0];
if(! -e $treeFileName) { die "ERROR: Cannot file tree file $treeFileName or not tree file name provided.\n"; }

use Bio::TreeIO;

sub get_node_distances {
    my ($tree_file) = @_;

    # Read the phylogenetic tree file
    my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $tree_file);
    my $tree = $treeio->next_tree;

    # Get all the leaf nodes
    my @leaves = $tree->get_leaf_nodes;

    # Calculate the distance from each node to the closest leaf
    my %node_distances;
    foreach my $node ($tree->get_nodes) {
        next if $node->is_Leaf;
        my $max_distance = 0;
        foreach my $child ($node->get_all_Descendents) { 
            my $distance = $tree->distance(-nodes => [$node, $child]);
            $max_distance = $distance if $distance > $max_distance;
        }
        $node_distances{$node->id} = $max_distance;
    }

    return %node_distances;
}

# Example usage
my %distances = get_node_distances($treeFileName);

# Print the distances
foreach my $node_id (keys %distances) {
    print "$node_id,$distances{$node_id}\n";
}

