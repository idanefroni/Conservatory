#!/usr/bin/perl -w
use strict;
use IO::Uncompress::Gunzip;
use Bio::AlignIO;

my $file = "/mnt/dd3/jdutheil/Data/Genomes/UCSC_46ways/maf/chr22.maf.gz";
my $ptr = new IO::Uncompress::Gunzip($file) or die $!;
my $alignio = Bio::AlignIO->new(-fh => $ptr, -format => 'maf');

open (CSV, '>chr22.statistics.perl.csv');
while(my $aln = $alignio->next_aln()){
  my $block_length = $aln->length;
  my $block_size = $aln->num_sequences;
  if ($block_length >= 1000) {
    my @seqs = ($aln->each_seq_with_id("hg19.chr22"));
    if ($#seqs == 0) {
      print "."; $|++;
      print CSV $seqs[0]->start() . "\t" . $seqs[0]->end() . "\t" . "$block_length\t$block_size\n";
    }
  }
}
close (CSV);

#1671.12user 0.19system 27:55.00elapsed 99%CPU (0avgtext+0avgdata 81856maxresident)k
#0inputs+64outputs (0major+5470minor)pagefaults 0swaps

