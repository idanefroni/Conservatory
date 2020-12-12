#!/usr/bin/perl
use strict;

my @lastline;
my @curline;

print "##gff-version 3\n";

while(<>)
{
  if(substr($_, 1,1) ne "#") { 
  	chomp($_);
  	if(scalar @lastline == 0) {
  		@lastline = split(/\t/, $_);
  		$lastline[2]="enhancer";
  	} else {
  		@curline = split(/\t/, $_);
  		$curline[2]="enhancer";
  		
  		if($lastline[4] >= $curline[3] && $lastline[8] eq $curline[8]) {
  			# if it is fully contained, ignore it. If not, combine and 
  			# calculate weighted average for the score
  			if($lastline[4] < $curline[4]) {
  				$lastline[5] = ($lastline[5] * abs($lastline[4] - $lastline[3]) + $curline[5] * abs($curline[4] - $curline[3]))/ (abs($lastline[4] - $lastline[3]) + abs($curline[4] - $curline[3]));
	  			$lastline[4] = $curline[4];
  			}
  		} else {  		
  			print $lastline[0] . "\t" . $lastline[1] . "\t" . $lastline[2] . "\t" . $lastline[3] . "\t" . $lastline[4] . "\t" . int($lastline[5]) . "\t" . $lastline[6] . "\t" . $lastline[7] . "\t" . $lastline[8] ."\n";
  			@lastline = @curline;
  		}
  	}
  }
}

#dump last line
@lastline=@curline;
$lastline[2]="enhancer";
print $lastline[0] . "\t" . $lastline[1] . "\t" . $lastline[2] . "\t" . $lastline[3] . "\t" . $lastline[4] . "\t" . int($lastline[5]) . "\t" . $lastline[6] . "\t" . $lastline[7] . "\t" . $lastline[8] ."\n";

