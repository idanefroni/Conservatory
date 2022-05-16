use strict;
my $FullProbFile=shift;
my $Node=shift;
my $NumOfSeqToSample=shift;
my $SeqType=shift; # aa | nuc
my $OutFile=shift;
my $isServer=shift;

my @AB=();
my $AB_SIZE;
if ($SeqType eq "nuc")
{
	@AB=qw(A C G T);
	$AB_SIZE=4;
}
if ($SeqType eq "aa")
{
	@AB=qw(A C D E F G H I K L M N P Q R S T V W Y);
	$AB_SIZE=20;
}
if ($SeqType eq "codon")
{
	@AB=qw(AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAC TAT TCA TCC TCG TCT TGC TGG TGT TTA TTC TTG TTT);
	$AB_SIZE=61;
}
my %ProbPerSite=(); # hash of array with prob for each pos
open (PROB_FILE,$FullProbFile) || die "Can't open The Full Prob File '$FullProbFile' $!";
my $SeqLength=0;
my $line=<PROB_FILE>; # header
while ($line=<PROB_FILE>)
{
	chomp ($line);
	my @line=split(",",$line); # NODE,SITE,PROBS BY AB
	my $CurrNode=shift(@line);
	my $CurrPos=shift(@line);
	if ($CurrNode eq $Node)
	{
		$ProbPerSite{$CurrPos}=[@line];
		$SeqLength=$CurrPos if ($CurrPos>$SeqLength);
	}
}
close (PROB_FILE);

open (OUT,">$OutFile") || die "Can't open Out: '$OutFile' $!";
for (my $SeqNum=0;$SeqNum<$NumOfSeqToSample;$SeqNum++)
{
	my $RandomSeq="";
	#if (($SeqType eq "aa") or ($SeqType eq "nuc"))
	#{
		for (my $pos=1;$pos<=$SeqLength;$pos++)
		{
			my $Rand=rand();
			my $i=0;
			my $Size=@{$ProbPerSite{$pos}};
			print "SIZE OF PROB VECTOR at POS $pos:$Size\n" if ($Size<$AB_SIZE);
			while(($Rand+0.0001 >= $ProbPerSite{$pos}[$i]) and ($i<$AB_SIZE-1))
			{
				$Rand=$Rand-$ProbPerSite{$pos}[$i];
				$i++;
			}
			print "UNDIFINED:$i for RAND $Rand and vector ",join (",",@{$ProbPerSite{$pos}}) if (!defined $AB[$i]);
			$RandomSeq=$RandomSeq.$AB[$i];
		}
	#}
	#elsif ($SeqType eq "codon")
	#{
	#	for (my $pos=1;$pos<=($SeqLength/3);$pos++)
	#	{
	#		my $Rand=rand();
	#		my $i=0;
	#		my $Size=@{$ProbPerSite{$pos}};
	#		print "SIZE OF PROB VECTOR at POS $pos:$Size\n" if ($Size<$AB_SIZE);
	#		while(($Rand+0.0001 >= $ProbPerSite{$pos}[$i]) and ($i<$AB_SIZE-1))
	#		{
	#			$Rand=$Rand-$ProbPerSite{$pos}[$i];
	#			$i++;
	#		}
	#		print "UNDIFINED:$i for RAND $Rand and vector ",join (",",@{$ProbPerSite{$pos}}) if (!defined $AB[$i]);
	#		$RandomSeq=$RandomSeq.$AB[$i];
	#	}
	#}
#	print "LENGTH:",length($RandomSeq),"\n";
	print OUT ">",$SeqNum+1,"\n$RandomSeq\n";
}
if ($isServer eq "YES")
{
	# Update the output page
    #######################################

	my $OutDir=getDir($OutFile);
	my $OutPage=$OutDir."output.html";
	if (-e $OutDir."output.php")
	{
		$OutPage=$OutDir."output.php";
	}

	open (OUTPUT,"$OutPage") || die "Can't open '$OutPage' $!";
	my @out=<OUTPUT>;
	close (OUTPUT);
	open (OUTPUT,">$OutPage");
	my $SampledSeq_Section=0;
	foreach my $line (@out)
	{
		if ($line=~/sequences from the posterior distribution for ancestral node/)
		{
			$SampledSeq_Section=1;
			print OUTPUT $line;
		}
		elsif (($line=~/form/) and ($SampledSeq_Section==1))
		{
			print OUTPUT $line;
			my $FileNoPath=getFilename($OutFile);
			print_message_to_output("<A HREF='$FileNoPath' TARGET=_blank>$NumOfSeqToSample sequences sampled from the posterior distribution for ancestral node $Node</A></p>");
			$SampledSeq_Section=0;
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

# Returns the filename without directory
sub getFilename{
	my $fullFile = pop @_;
	if ($fullFile =~ m/.*[\\\/](.*)$/) {
		return $1;
	} else {return $fullFile}

}

sub getDir{
	my $fullFile = pop @_;
	if ($fullFile =~ m/(.*[\\\/]).*$/) {
		return $1;
	} else {return ''}
	
}
