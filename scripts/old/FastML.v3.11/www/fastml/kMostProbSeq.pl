use strict;

my $FullLogLikeFile=shift;
my $Node=shift;
my $k=shift;
my $OutDir=shift;
my $seqType=shift;
my $isServer=shift;

$OutDir=$OutDir."/" if ($OutDir!~/\/$/);
my $K_MOST_PROB_SEQ="python /bioseq/pupkoSVN/trunk/www/FastML/kMostProbSeq.py";
my $ProbMatrix=$OutDir."$Node.LogLikelihoodMarginalProb.csv";

my %Profile=(); # Lines=AlephBet size; Col:#Pos
open (FULL_PROB,$FullLogLikeFile) || die "Can't open The Full Log Like File '$FullLogLikeFile' $!";
open (OUT,">$ProbMatrix") || die "Can't open Prob matrix file: '$ProbMatrix' $!";
print OUT "site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y\n" if ($seqType eq "aa");
print OUT "site,A,C,G,T\n" if ($seqType eq "nuc");
print OUT "site,AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA,AGC,AGG,AGT,ATA,ATC,ATG,ATT,CAA,CAC,CAG,CAT,CCA,CCC,CCG,CCT,CGA,CGC,CGG,CGT,CTA,CTC,CTG,CTT,GAA,GAC,GAG,GAT,GCA,GCC,GCG,GCT,GGA,GGC,GGG,GGT,GTA,GTC,GTG,GTT,TAC,TAT,TCA,TCC,TCG,TCT,TGC,TGG,TGT,TTA,TTC,TTG,TTT\n" if ($seqType eq "codon");
while (my $line=<FULL_PROB>)
{
	my @line=split(",",$line); # NODE,SITE,PROBS BY AB
	my $CurrNode=shift(@line);
	if ($CurrNode eq "$Node")
	{
		print OUT join(",",@line);
	}
}
close (FULL_PROB);
close (OUT);
my $OutSeq=$OutDir.$Node.".".$k."MostProbSeq.fasta";
my $cmd="$K_MOST_PROB_SEQ -i $ProbMatrix -o $OutSeq -k $k";
system ($cmd);

if ($isServer eq "YES")
{
	# Update the output page
    #######################################

	my $OutPage=$OutDir."output.html";
	if (-e $OutDir."output.php")
	{
		$OutPage=$OutDir."output.php";
	}

	open (OUTPUT,"$OutPage") || die "Can't open '$OutPage' $!";
	my @out=<OUTPUT>;
	close (OUTPUT);
	open (OUTPUT,">$OutPage");
	my $kMostProb_Section=0;
	foreach my $line (@out)
	{
		if ($line=~/most likely ancestral sequences for ancestral node/)
		{
			$kMostProb_Section=1;
			print OUTPUT $line;
		}
		elsif (($line=~/form/) and ($kMostProb_Section==1))
		{
			print OUTPUT $line;
			my $FileNoPath=$Node.".".$k."MostProbSeq.fasta";
			print_message_to_output("<A HREF='$FileNoPath' TARGET=_blank>$k-most likely ancestral sequences for ancestral node $Node ('marginal' reconstruction)</A></p>");
			$kMostProb_Section=0;
		}
		else
		{
			print OUTPUT $line;
		}
	}
	close (OUTPUT);
}

#---------------------------------------------
sub print_message_to_output{
#---------------------------------------------
    my $msg = shift;
    print OUTPUT "\n<ul><li>$msg</li></ul>\n";
}

