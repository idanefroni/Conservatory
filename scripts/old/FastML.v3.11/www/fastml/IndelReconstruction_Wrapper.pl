use strict;
use Getopt::Long;
use FindBin qw($Bin); # www/FastML_2012/
use File::Copy;
die "USAGE: --MSA_File <MSA_File> --Tree_File <Tree_File> --outDir <outDir> --seqType <aa|nuc|codon>
Optional parameters:
      --indelCutOff <Cutoff for indel vs Char> deafult =0.5

	  --CharsMarginalProb <ProbFile> (deafult=prob.marginal.txt) - prob of ancestral sequences - FASTML output

	  --ML_GapOut <Indels_ML_Prob>         # (deafult: IndelsMarginalProb.txt - IndelsMarginalProb (IndelReconstructOutput)
	  --ML_Ancestral_MSA <Ancestal_ML_MSA> # (deafult: seq.marginal_IndelAndChars.txt) - output for Chars and Gap Ancestral Reconstruction - MSA;
	  --ML_Chars_ML_Gap <AncestralProb>    # (deafult: AncestralMaxMarginalProb_Char_Indel.txt) - File with the max prob of each position on each node

	  --MP_GapOut <Indels_MP_State>        # (deafult: Indels.parsimony) - Indel Satate for each MSA pos by parsimony
	  --ML_Char_MP_Gap <ML_Char_MP_indels> # (deafult: AncestralMaxProbMarginal_Char_Parsimony_Indel.txt) - File with the max prob char of each position on each node and indel parsimony
	  --Ancestral_MSA_MP_GAP <MSA_MP_Gap>  # (deafult: seq.marginal_Chars_ParsimonyIndels.txt) - MSA Output for Chars and Parsimonuis Gap Ancestral Reconstruction;

	  --Debug                              # (deafult: off) printouts debug info
" unless (@ARGV >= 1);

# Assign default
my ($MSA_File,$OutDir,$Tree_File,$IndelsCutoff,$SeqType,$MarginalProb_of_Chars,$GapProb_OutFile,$Ancestral_MSA,$Ancestral_Prob,$GapParsimony_OutFile,$Ancestral_Prob_ParsimonyIndel,$Ancestral_MSA_Parsimony,$DEBUG_F,);
$MSA_File="";
$OutDir="";
$Tree_File="";

$IndelsCutoff=0.5;

my $getoptResult = GetOptions ("MSA_File=s"=>\$MSA_File,  # = means that this parameter is required, s means string
							   "outDir=s"=>\$OutDir,
							   "Tree_File=s"=>\$Tree_File,
							   "seqType=s"=>\$SeqType,                                       # aa|nuc|codon
							   "indelCutOff:f"=>\$IndelsCutoff,

							   "CharsMarginalProb:s"=>\$MarginalProb_of_Chars,               # (prob.marginal.txt) - prob of ancestral sequences - FASTML output

							   "ML_GapOut:s"=>\$GapProb_OutFile,                             # (IndelsMarginalProb.txt) - IndelsMarginalProb (IndelReconstructOutput)
							   "ML_Ancestral_MSA:s"=>\$Ancestral_MSA,                        # (seq.marginal_IndelAndChars.txt) - output for Chars and Gap Ancestral Reconstruction - MSA;
							   "ML_Chars_ML_Gap:s"=>\$Ancestral_Prob,                        # (Ancestral_MaxMarginalProb_Char_Indel.txt) - File with the max prob of each position on each node

							   "MP_GapOut:s"=>\$GapParsimony_OutFile,                        # (Indels.parsimony.txt) - Indel Satate for each MSA pos by parsimony
							   "ML_Char_MP_Gap:s"=>\$Ancestral_Prob_ParsimonyIndel,          # (Ancestral_MaxProb_Marginal_Char_Parsimony_Indel.txt) - File with the max prob char of each position on each node and indel parsimony
							   "Ancestral_MSA_MP_GAP:s"=>\$Ancestral_MSA_Parsimony,          # (seq.marginal_Chars_ParsimonyIndels.txt) - MSA Output for Chars and Parsimonuis Gap Ancestral Reconstruction;

							   "Debug" =>\$DEBUG_F,
							   );

# default file names
if ($OutDir!~/\/$/) {$OutDir.="/";}
$GapProb_OutFile=$OutDir."IndelsMarginalProb.txt"                  if ((!defined $GapProb_OutFile) or ($GapProb_OutFile eq ""));
$MarginalProb_of_Chars=$OutDir."prob.marginal.txt"                 if ((!defined $MarginalProb_of_Chars) or ($MarginalProb_of_Chars eq ""));
$Ancestral_MSA=$OutDir."seq.marginal_IndelAndChars.txt"            if ((!defined $Ancestral_MSA) or ($Ancestral_MSA eq ""));
$Ancestral_Prob=$OutDir."Ancestral_MaxMarginalProb_Char_Indel.txt"  if ((!defined $Ancestral_Prob) or ($Ancestral_Prob eq ""));

# default file names for PARSIMONY BASED OUTPUT
$GapParsimony_OutFile=$OutDir."Indels.parsimony.txt"                                           if ((!defined $GapParsimony_OutFile) or ($GapParsimony_OutFile eq ""));                   # Indel Satate for each MSA pos by parsimony
$Ancestral_Prob_ParsimonyIndel=$OutDir."Ancestral_MaxProb_Marginal_Char_Parsimony_Indel.txt" if ((!defined $Ancestral_Prob_ParsimonyIndel) or ($Ancestral_Prob_ParsimonyIndel eq "")); # File with the max prob char of each position on each node and indel parsimony
$Ancestral_MSA_Parsimony=$OutDir."seq.marginal_Chars_ParsimonyIndels.txt"                  if ((!defined $Ancestral_MSA_Parsimony) or ($Ancestral_MSA_Parsimony eq ""));             # Output for parsimony Chars and Gap Ancestral Reconstruction;

my $DEBUG="NO";
$DEBUG="YES" if ($DEBUG_F);

print "
--MSA_File=$MSA_File
--outDir=$OutDir
--Tree_File=$Tree_File
--seqType=$SeqType
--indelCutOff=$IndelsCutoff

--CharsMarginalProb=$MarginalProb_of_Chars

--ML_GapOut=$GapProb_OutFile
--ML_Ancestral_MSA=$Ancestral_MSA
--ML_Chars_ML_Gap=$Ancestral_Prob

--MP_GapOut=$GapParsimony_OutFile
--ML_Char_MP_Gap=$Ancestral_Prob_ParsimonyIndel
--Ancestral_MSA_MP_GAP=$Ancestral_MSA_Parsimony

--Debug=$DEBUG\n";

#print "WAIT...\n";<STDIN>;


# Constants
my $ParsimonyCostMatrix=2;
my $MSA_Prefix_Name="";
if ($MSA_File=~/([^\/]+?)(.aln|.faa|.mfa|.txt)?$/)
{
	$MSA_Prefix_Name=$1;
}
else
{
	$MSA_Prefix_Name=$MSA_File;
}
$DEBUG=uc($DEBUG);
if (!defined $DEBUG)
{
	$DEBUG="NO";
}

# Programs Path
#my $IndelCoder="/bioseq/FastML/IndelReconstruction/indelCoder";
#my $IndelCoder="/bioseq/FastML/IndelReconstruction/indelCoder.V1.6";
#my $IndelCoder="/bioseq/FastML/IndelReconstruction/indelCoder.V1.71";
my $IndelCoder="$Bin/../../programs/indelCoder/indelCoder";
#my $IndelReconstruction="/bioseq/FastML/IndelReconstruction/gainLoss.V9.9822"; # by gainLoss
#my $IndelReconstruction="/bioseq/FastML/IndelReconstruction/gainLoss.V9.9863"; # by gainLoss
my $IndelReconstruction="$Bin/../../programs/gainLoss/gainLoss"; # by gainLoss

# Globals File Names
$OutDir=$OutDir."/" if ($OutDir!~/\/$/);
my $Indels_Reconstruction_results_Dir=$OutDir."IndelsReconstruction/";
# IndelCoder
my $IndelCoderParamFile="IndelCoderParamFile";
my $indelOutputFastaFile="$Indels_Reconstruction_results_Dir/$MSA_Prefix_Name".".indelOutputFastaFile";
my $indelOutputInfoFile="$Indels_Reconstruction_results_Dir/$MSA_Prefix_Name".".indelOutputInfoFile";
my $nexusFileName="$Indels_Reconstruction_results_Dir/$MSA_Prefix_Name".".indel_nexusFile";
my $indelLogFile="$Indels_Reconstruction_results_Dir/$MSA_Prefix_Name"."IndelCoder.log";

# Indel Reconstruction
my $IndelReconstructionParamFile="IndelReconstructionParamFile";
#my $indelOutputFasta_NO_MISSING_DATA_File="$Indels_Reconstruction_results_Dir/$MSA_Prefix_Name"."_MISING_DATA_TO0.indelOutputFastaFile"; # For now gainLoss don't handle missing data so we replace '?' with 0
my $AncestralReconstructIndelPosterior="$Indels_Reconstruction_results_Dir/RESULTS/AncestralReconstructPosterior.txt"; # The file with ancestral prob of indel
my $AncestralReconstructParsimony="$Indels_Reconstruction_results_Dir/RESULTS/gainLossMP.".$ParsimonyCostMatrix.".AncestralReconstructSankoff.txt";
# Joint character based Ancestral MSA with Indel Reconstruction

mkdir ($Indels_Reconstruction_results_Dir);

my %Species_On_MSA=(); # All species in the MSA - MAYBE TO REMOVE
open (MSA,$MSA_File);
while (my $line=<MSA>)
{
	chomp ($line);
	if ($line=~/^>(.*)/)
	{
		$Species_On_MSA{$1}=1;
	}
}
# Read MSA to Hash
my $MSA_Hash_ref=readMSA($MSA_File);
my %MSA_Hash=%{$MSA_Hash_ref};

# Prepare indel Coder ParamFile
open (INDEL_CODER_PARAMS,">$Indels_Reconstruction_results_Dir$IndelCoderParamFile") || die "IndelReconstruction_Wrapper: Can't open IndelCoderParamFile '$Indels_Reconstruction_results_Dir$IndelCoderParamFile' $!";
print INDEL_CODER_PARAMS "_seqFile $MSA_File\n";
print INDEL_CODER_PARAMS "_indelOutputInfoFile $indelOutputInfoFile\n";
print INDEL_CODER_PARAMS "_indelOutputFastaFile $indelOutputFastaFile\n";
print INDEL_CODER_PARAMS "_nexusFileName $nexusFileName\n";
print INDEL_CODER_PARAMS "_logFile $indelLogFile\n";
print INDEL_CODER_PARAMS "_logValue 9\n";
print INDEL_CODER_PARAMS "_codingType SIC\n";
print INDEL_CODER_PARAMS "_isOmitLeadingAndEndingGaps 0\n";

close (INDEL_CODER_PARAMS);

system ("cd $Indels_Reconstruction_results_Dir; $IndelCoder $IndelCoderParamFile");

if (!-e $indelOutputFastaFile)
{
	die "IndelReconstruction_Wrapper: $indelOutputFastaFile was not created or empty, please have a look on the indel coder log file at: $indelLogFile";
}

# Run indelReconstruction by gainLoss
my $removed_BP_InternalNodeName=remove_InternalNodeName_or_BPvalues($Tree_File,$Tree_File.".Orig");
copy ($Tree_File,"$Tree_File.ForIndelReconstruction");
move ("$Tree_File.Orig",$Tree_File) if (-e "$Tree_File.Orig");
open (INDEL_RECONSTRUCTION_PARAMS,">$Indels_Reconstruction_results_Dir$IndelReconstructionParamFile") || die "Can't open IndelReconstructionParamFile '$Indels_Reconstruction_results_Dir$IndelReconstructionParamFile' $!";
print INDEL_RECONSTRUCTION_PARAMS "_seqFile $indelOutputFastaFile\n";
print INDEL_RECONSTRUCTION_PARAMS "_treeFile $Tree_File.ForIndelReconstruction\n";
print INDEL_RECONSTRUCTION_PARAMS "_isRootFreqEQstationary 1\n";
print INDEL_RECONSTRUCTION_PARAMS "_calculateAncestralReconstruct 1\n";
print INDEL_RECONSTRUCTION_PARAMS "_costMatrixGainLossRatio  2\n";
print INDEL_RECONSTRUCTION_PARAMS "_minNumOfOnes 1\n";
close (INDEL_RECONSTRUCTION_PARAMS);
system ("cd $Indels_Reconstruction_results_Dir; $IndelReconstruction $IndelReconstructionParamFile");
my %MSA_Pos_Species_to_Indel=();
my %MSAtoIndel=();
my ($MSA_Pos_Species_to_Indel,$MSAtoIndel)=Read_MSA_to_Indels_Info($indelOutputInfoFile,\%MSA_Pos_Species_to_Indel,\%MSAtoIndel); # hash1 - key1:MSA_Pos,key2:species; value:IndelMSAPos; hash2 - key: MSA_Pos;value: IndelsMSA_Pos (array)
my %AncestralReconstructIndelPosterior_Hash=();
my $AncestralReconstructIndelPosterior_Reff=Read_Ancestral_Prob_For_Indel($AncestralReconstructIndelPosterior,\%AncestralReconstructIndelPosterior_Hash); # hash = key1:IndelMSA_Pos,key2:species; value Prob for indel

####### HADLE WITH PROB RECONSTRUCTION
%AncestralReconstructIndelPosterior_Hash=%$AncestralReconstructIndelPosterior_Reff;
my %MSA_Pos_Species_AncestorIndelProb=(); # Will hold for each MSA_Pos and Species the vector of IndelPos_ProbOfIndel
print "HADLE WITH PROB RECONSTRUCTION LOOP\n====================================================\n" if ($DEBUG eq "YES");
## MAKE UNIQ
print "+++++++++++++++++DEBUG - PRINT INDEL POS TO INDEL NOT UNIQ ++++++++++++++++++++++\n" if ($DEBUG eq "YES");
foreach my $MSA_Pos (sort {$a<=>$b} keys %$MSAtoIndel)
{
	print "MSA:$MSA_Pos\t",join(",",@{$MSAtoIndel->{$MSA_Pos}}),"\n" if ($DEBUG eq "YES");
	my $tmp_array=uniq_array($MSAtoIndel->{$MSA_Pos});
	$MSAtoIndel->{$MSA_Pos}=[@{$tmp_array}];
}
print "+++++++++++++++++DEBUG - PRINT INDEL POS TO INDEL UNIQ +++++++++++++++++++++++++\n" if ($DEBUG eq "YES");
foreach my $MSA_Pos (sort {$a<=>$b} keys %$MSAtoIndel)
{
	print "MSA:$MSA_Pos\t",join(",",@{$MSAtoIndel->{$MSA_Pos}}),"\n" if ($DEBUG eq "YES");
}
print "+++++++++++++++++ END DEBUG ++++++++++++++++++++++++\n" if ($DEBUG eq "YES");

foreach my $MSA_Pos (sort {$a<=>$b} keys %$MSAtoIndel)
{
	print "MSA:$MSA_Pos," if ($DEBUG eq "YES"); # DEBUG
	foreach my $IndelPos (@{$MSAtoIndel->{$MSA_Pos}})
	{
		print "Indel:$IndelPos - $AncestralReconstructIndelPosterior_Hash{$IndelPos}"  if ($DEBUG eq "YES"); # empty
		foreach my $species (keys %{$AncestralReconstructIndelPosterior_Hash{$IndelPos}})
		{
			if (!exists $Species_On_MSA{$species}) # Ancestral Node # CONSIDER REMOVE
			{
				my $IndelPos_ProbOfIndel=$IndelPos."_".$AncestralReconstructIndelPosterior_Hash{$IndelPos}{$species};
				if (!exists $MSA_Pos_Species_AncestorIndelProb{$MSA_Pos}{$species}){$MSA_Pos_Species_AncestorIndelProb{$MSA_Pos}{$species}=[$IndelPos_ProbOfIndel];}
				else {push @{$MSA_Pos_Species_AncestorIndelProb{$MSA_Pos}{$species}},$IndelPos_ProbOfIndel;}
				print "$MSA_Pos\t$IndelPos\t$species\t$AncestralReconstructIndelPosterior_Hash{$IndelPos}{$species}\n"  if ($DEBUG eq "YES"); # DEBUG
			}
		}
	}
}
open (GAP_PROB,">$GapProb_OutFile") || die "Can't open '$GapProb_OutFile' $!";
print GAP_PROB "Pos\tNode\tProb_Of_Indel\n";
my %MSA_Pos_Node_MaxProbOf_Gap=();
foreach my $MSA_Pos (sort {$a<=>$b} keys %MSA_Pos_Species_AncestorIndelProb)
{
	foreach my $species (sort keys %{$MSA_Pos_Species_AncestorIndelProb{$MSA_Pos}})
	{
		if (!exists $Species_On_MSA{$species}) # Ancestral Node # CONSIDER REMOVE
		{
			print "$MSA_Pos\t$species"  if ($DEBUG eq "YES");
			print GAP_PROB "$MSA_Pos\t$species";
			my $Uniq_Indels_Reff=uniq_array($MSA_Pos_Species_AncestorIndelProb{$MSA_Pos}{$species});
			my @Uniq_Indels=@$Uniq_Indels_Reff;
			my $NumOfIndelCoverMSA_Pos=@Uniq_Indels;
			my @ProbsOfIndel;
			for (my $i=0;$i<$NumOfIndelCoverMSA_Pos;$i++)
			{
				my $Indel_IndelProb=$Uniq_Indels[$i];
				my ($Indel_Pos,$IndelProb)=split("_",$Indel_IndelProb);
				print "\t$Indel_Pos:$IndelProb"  if ($DEBUG eq "YES");
				push (@ProbsOfIndel,$IndelProb);
			}
			my $maxProbOfIndel = (sort { $b <=> $a } @ProbsOfIndel)[0];
			print "\tMAX:$maxProbOfIndel\n"  if ($DEBUG eq "YES");
			print GAP_PROB "\t$maxProbOfIndel\n";
			$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$species}=$maxProbOfIndel;
		}
	}
}
close (GAP_PROB);
my %MSA_Pos_Node_Char_or_Gap=();
# Read the Chars Marginal Prob
my ($MSA_Pos_Node_Char_Marginal_Prob_Reff,$Nodes_Name_Reff,$MSA_Length)=Read_Char_Marginal_Prob($MarginalProb_of_Chars);
print "MSA_Length:$MSA_Length\n"  if ($DEBUG eq "YES");
my @Nodes=@$Nodes_Name_Reff;
open (ANCESTRAL_PROB,">$Ancestral_Prob")|| die "Can't open Ancestral Prob File: '$Ancestral_Prob' $!\n";
print ANCESTRAL_PROB "Pos_on_MSA\tNode\tChar\tCharProb\n";
foreach my $MSA_Pos (sort {$a<=>$b} keys %{$MSA_Pos_Node_Char_Marginal_Prob_Reff})
{
	print "MSA:$MSA_Pos\n"  if ($DEBUG eq "YES");
	foreach my $Node (sort keys %{$MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}})
	{
		my $maxProbChar="NA";
		my $maxProb=0;
		my $Num_Of_1=0;
		foreach my $Char (sort keys %{$MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}->{$Node}})
		{
			
			if (($MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}->{$Node}->{$Char}>$maxProb)&&(defined $MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}->{$Node}->{$Char}))
			{
				$maxProbChar=$Char;
				$maxProb=$MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}->{$Node}->{$Char};
			}
			$Num_Of_1++ if ($MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}->{$Node}->{$Char}==1);
		}

		# Decide what is the most probable char on pos
		if ($Num_Of_1>1) # GAP 
		{
			if ($SeqType eq "codon")
			{
				$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}="---".":1";
			}
			else
			{
				$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}="-".":1";
			}
			
			$maxProbChar="NA";
			$maxProb=0;
		}
		else
		{
			if (!exists $MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node}){$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node}="NA";}
			print "NODE:$Node - $maxProbChar:$maxProb ? -:$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node}\n"  if ($DEBUG eq "YES");#<STDIN>; # DEBUG
			if (($SeqType eq "aa") or ($SeqType eq "nuc"))
			{
				if ($MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node} eq "NA")
				{
					$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}=$maxProbChar.":".$maxProb;
				}
  #			    elsif ($maxProb>=$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node}) # MOST PROBALBE IS THE CHAR
				#elsif ($MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node}<(1-$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node})) # MOST PROBALBE IS THE CHAR
				elsif ($MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node}<$IndelsCutoff) # MOST PROBALBE IS THE CHAR
				{
					$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}=$maxProbChar.":".$maxProb;
				}
				else
				{
					$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}="-".":".$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node};
					#$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}="---".":".$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node} if ($SeqType eq "codon");
				}
			}
			elsif ($SeqType eq "codon")
			{
				# MSA Pos is according to the codon number (i.e ((MSA_Pos-1)/3)+1)
				my $MSA_Pos_GAP=(($MSA_Pos-1)*3)+1; # The real char on the MSA
				if (!exists $MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos_GAP}{$Node}){$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos_GAP}{$Node}="NA";}
				if ($MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos_GAP}{$Node} eq "NA")
				{
					$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}=$maxProbChar.":".$maxProb;
				}
  #			    elsif ($maxProb>=$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node}) # MOST PROBALBE IS THE CHAR
				#elsif ($MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos_GAP}{$Node}<(1-$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos_GAP}{$Node})) # MOST PROBALBE IS THE CHAR
				elsif ($MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos_GAP}{$Node}<$IndelsCutoff) # MOST PROBALBE IS THE CHAR
				{
					$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}=$maxProbChar.":".$maxProb;
				}
				else
				{
					$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}="---".":".$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos_GAP}{$Node};
				}
			}
		}
		my ($CharForPrint,$ProbForPrint)=split(/:/,$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node});
		if ($SeqType eq "codon")
		{
			my $MSA_Pos_GAP=(($MSA_Pos-1)*3)+1; # The real char
			print ANCESTRAL_PROB "$MSA_Pos_GAP\t$Node\t$CharForPrint\t$ProbForPrint\n";#$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}\n";
		}
		else
		{
			print ANCESTRAL_PROB "$MSA_Pos\t$Node\t$CharForPrint\t$ProbForPrint\n";#$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}\n";
		}
	}
} 


### PRINT THE GAP and CHAR Ancestral MSA
open (MSA_OUT,">$Ancestral_MSA") || die "Can't open Output MSA: '$Ancestral_MSA' $!\n";
foreach my $Node (@Nodes)
{
	if (exists $MSA_Hash{$Node}) # Original sequence
	{
		print MSA_OUT ">$Node\n";
		print MSA_OUT "$MSA_Hash{$Node}\n";
	}
	else # Ancestral seq
	{
		print MSA_OUT ">$Node\n";
		for (my $i=1;$i<=$MSA_Length;$i++)
		{
			my ($Char,$Prob)=split(":",$MSA_Pos_Node_Char_or_Gap{$i}{$Node});
			print MSA_OUT $Char;
		}
		print MSA_OUT "\n";
	}
}
### TO HERE

# For Parsimony (COPY OF THE CODE ABOVE...) TO DO: CHANGE IT SOME DAY...
my %AncestralReconstructIndelParsimony_Hash=();
my $AncestralReconstructIndelParsimony_Reff=Read_Ancestral_Parsimony_State($AncestralReconstructParsimony,\%AncestralReconstructIndelParsimony_Hash);     # hash = key1:IndelMSA_Pos,key2:species; value 1 for indel 0 for char
%AncestralReconstructIndelParsimony_Hash=%$AncestralReconstructIndelParsimony_Reff;

my %MSA_Pos_Species_AncestorIndelParsimony=(); # Will hold for each MSA_Pos and Species the vector of IndelPos_ProbOfIndel
foreach my $MSA_Pos (sort {$a<=>$b} keys %$MSAtoIndel)
{
#	print "MSA:$MSA_Pos,";
	foreach my $IndelPos (@{$MSAtoIndel->{$MSA_Pos}})
	{
#		print "Indel:$IndelPos - $AncestralReconstructIndelPosterior_Hash{$IndelPos}"; # empty
		foreach my $species (keys %{$AncestralReconstructIndelParsimony_Hash{$IndelPos}})
		{
			my $IndelPos_ProbOfIndel=$IndelPos."_".$AncestralReconstructIndelParsimony_Hash{$IndelPos}{$species};
			if (!exists $MSA_Pos_Species_AncestorIndelParsimony{$MSA_Pos}{$species}){$MSA_Pos_Species_AncestorIndelParsimony{$MSA_Pos}{$species}=[$IndelPos_ProbOfIndel];}
			else {push @{$MSA_Pos_Species_AncestorIndelParsimony{$MSA_Pos}{$species}},$IndelPos_ProbOfIndel;}
#			print "$MSA_Pos\t$IndelPos\t$species\t$AncestralReconstructIndelPosterior_Hash{$IndelPos}{$species}\n";
		}
	}
}
open (GAP_PARSIMONY,">$GapParsimony_OutFile") || die "Can't open '$GapProb_OutFile' $!";
print GAP_PARSIMONY "Pos\tNode\tGap\n";
my %MSA_Pos_Node_ParsimonyOf_Gap=();
foreach my $MSA_Pos (sort {$a<=>$b} keys %MSA_Pos_Species_AncestorIndelParsimony)
{
	foreach my $species (sort keys %{$MSA_Pos_Species_AncestorIndelParsimony{$MSA_Pos}})
	{
		print "$MSA_Pos\t$species"  if ($DEBUG eq "YES");
		print GAP_PARSIMONY "$MSA_Pos\t$species" if ($species=~/^N\d+$/); # print only ancestral nodes
		my $Uniq_Indels_Reff=uniq_array($MSA_Pos_Species_AncestorIndelParsimony{$MSA_Pos}{$species});
		my @Uniq_Indels=@$Uniq_Indels_Reff;
		my $NumOfIndelCoverMSA_Pos=@Uniq_Indels;
		my @ParsimonyOfIndel;
		for (my $i=0;$i<$NumOfIndelCoverMSA_Pos;$i++)
		{
			my $Indel_IndelParsimony=$Uniq_Indels[$i];
			my ($Indel_Pos,$IndelParsimony)=split("_",$Indel_IndelParsimony);
			print "\t$Indel_Pos:$IndelParsimony"  if ($DEBUG eq "YES");
			push (@ParsimonyOfIndel,$IndelParsimony);
		}
		# my $minProbOfIndel = (sort { $a <=> $b } @ParsimonyOfIndel)[0]; # WE GAVE PRIORITY TO CHAR (used when we had old (<=1.71) indelCoder)
		# print "\tMAX:$minProbOfIndel\n"  if ($DEBUG eq "YES");
		# print GAP_PARSIMONY "\t$minProbOfIndel\n";
		# $MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos}{$species}=$minProbOfIndel;
		my $maxProbOfIndel = (sort { $b <=> $a } @ParsimonyOfIndel)[0]; 
		print "\tMAX:$maxProbOfIndel\n"  if ($DEBUG eq "YES");
		print GAP_PARSIMONY "\t$maxProbOfIndel\n" if ($species=~/^N\d+$/); # print only ancestral nodes;
		$MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos}{$species}=$maxProbOfIndel;
	}
}
close (GAP_PARSIMONY);
my %MSA_Pos_Node_Char_or_Gap_Parsimony=();
open (ANCESTRAL_PROB_PARSIMONY_INDEL,">$Ancestral_Prob_ParsimonyIndel")|| die "IndelReconstruction_Wrapper::Can't open Ancestral Prob Parsimony Indel File: '$Ancestral_Prob_ParsimonyIndel' $!\n";
print ANCESTRAL_PROB_PARSIMONY_INDEL "Pos_on_MSA\tNode\tChar\tCharProb\n";
foreach my $MSA_Pos (sort {$a<=>$b} keys %{$MSA_Pos_Node_Char_Marginal_Prob_Reff})
{
	print "MSA:$MSA_Pos\n"  if ($DEBUG eq "YES");
	foreach my $Node (sort keys %{$MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}})
	{
		my $maxProbChar="NA";
		my $maxProb=0;
		my $Num_Of_1=0;
		foreach my $Char (sort keys %{$MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}->{$Node}})
		{
			
			if (($MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}->{$Node}->{$Char}>$maxProb)&&(defined $MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}->{$Node}->{$Char}))
			{
				$maxProbChar=$Char;
				$maxProb=$MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}->{$Node}->{$Char};
			}
			$Num_Of_1++ if ($MSA_Pos_Node_Char_Marginal_Prob_Reff->{$MSA_Pos}->{$Node}->{$Char}==1);
		}
		# Decide what is the most probable char on pos
		if ($Num_Of_1>1) # GAP ON ORIGINAL SEQ (NOT ANCESTRAL)
		{
			if ($SeqType eq "codon")
			{
				$MSA_Pos_Node_Char_or_Gap_Parsimony{$MSA_Pos}{$Node}="---".":1";
			}
			else
			{
				$MSA_Pos_Node_Char_or_Gap_Parsimony{$MSA_Pos}{$Node}="-".":1";
			}
			$maxProbChar="NA";
			$maxProb=0;
		}
		else
		{
			
			if (($SeqType eq "aa") or ($SeqType eq "nuc"))
			{
#			print "NODE:$Node - $maxProbChar:$maxProb ? -:$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node}\n";#<STDIN>;
				if (!exists $MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos}{$Node})
				{
					$MSA_Pos_Node_Char_or_Gap_Parsimony{$MSA_Pos}{$Node}=$maxProbChar.":".$maxProb;
				}
#			elsif ($maxProb>=$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node}) # MOST PROBALBE IS THE CHAR
				elsif ($MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos}{$Node}==0) # NO GAP BY PARSIMONY - MOST PROBALBE IS THE CHAR
				{
					$MSA_Pos_Node_Char_or_Gap_Parsimony{$MSA_Pos}{$Node}=$maxProbChar.":".$maxProb;
				}
				elsif ($MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos}{$Node}==1)
				{
					$MSA_Pos_Node_Char_or_Gap_Parsimony{$MSA_Pos}{$Node}="-".":"."1";
					$MSA_Pos_Node_Char_or_Gap_Parsimony{$MSA_Pos}{$Node}="---".":"."1" if ($SeqType eq "codon");
				}
			}
			elsif ($SeqType eq "codon")
			{
				# MSA Pos is according to the codon number (i.e ((MSA_Pos-1)/3)+1)
				my $MSA_Pos_GAP=(($MSA_Pos-1)*3)+1; # The real char on the MSA
				if (!exists $MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos_GAP}{$Node}){$MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos_GAP}{$Node}="NA";}
				if ($MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos_GAP}{$Node} eq "NA")
				{
					$MSA_Pos_Node_Char_or_Gap_Parsimony{$MSA_Pos}{$Node}=$maxProbChar.":".$maxProb;
				}
				#			    elsif ($maxProb>=$MSA_Pos_Node_MaxProbOf_Gap{$MSA_Pos}{$Node}) # MOST PROBALBE IS THE CHAR
				#elsif ($MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos_GAP}{$Node}<(1-$MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos_GAP}{$Node})) # MOST PROBALBE IS THE CHAR
				elsif ($MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos_GAP}{$Node}<$IndelsCutoff) # MOST PROBALBE IS THE CHAR
				{
					$MSA_Pos_Node_Char_or_Gap_Parsimony{$MSA_Pos}{$Node}=$maxProbChar.":".$maxProb;
				}
				else
				{
					$MSA_Pos_Node_Char_or_Gap_Parsimony{$MSA_Pos}{$Node}="---".":".$MSA_Pos_Node_ParsimonyOf_Gap{$MSA_Pos_GAP}{$Node};
				}
			}
		}
		my ($CharForPrint,$ProbForPrint)=split(/:/,$MSA_Pos_Node_Char_or_Gap_Parsimony{$MSA_Pos}{$Node});
		if ($SeqType eq "codon")
		{
			my $MSA_Pos_GAP=(($MSA_Pos-1)*3)+1; # The real char on the MSA
			print ANCESTRAL_PROB_PARSIMONY_INDEL "$MSA_Pos_GAP\t$Node\t$CharForPrint\t$ProbForPrint\n";#$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}\n";
		}
		else
		{
			print ANCESTRAL_PROB_PARSIMONY_INDEL "$MSA_Pos\t$Node\t$CharForPrint\t$ProbForPrint\n";#$MSA_Pos_Node_Char_or_Gap{$MSA_Pos}{$Node}\n";
		}
	}
} 


### PRINT THE GAP and CHAR Ancestral MSA
open (MSA_OUT_PARSIMONY,">$Ancestral_MSA_Parsimony") || die "Can't open Output MSA PARSIMONY : '$Ancestral_MSA_Parsimony' $!\n";
foreach my $Node (@Nodes)
{
	if (exists $MSA_Hash{$Node}) # Original sequence
	{
		print MSA_OUT_PARSIMONY ">$Node\n";
		print MSA_OUT_PARSIMONY "$MSA_Hash{$Node}\n";
	}
	else
	{
		print MSA_OUT_PARSIMONY ">$Node\n";
		for (my $i=1;$i<=$MSA_Length;$i++)
		{
			my ($Char,$Prob)=split(":",$MSA_Pos_Node_Char_or_Gap_Parsimony{$i}{$Node});
			print MSA_OUT_PARSIMONY $Char;
		}
		print MSA_OUT_PARSIMONY "\n";
	}
}
close (MSA_OUT_PARSIMONY);


sub Read_MSA_to_Indels_Info
# Will create an hash that map each position on the MSA to the translated indel (or indels)
{
#character number: 0
#Start position relative to MSA: 0
#End position relative to MSA: 1
#Length: 1
#Found in species: DQ373066.PTT Start position relative to genome: 0 Length: 1
#ENDCHARACTER

	print "MAPPING MSA POS TO INDEL\n==============================================================\n"  if ($DEBUG eq "YES");
	my $IndelInfo=shift;
	my $MSA_Pos_Species_to_Indel_Reff=shift;
	my $MSAtoIndel_Reff=shift;
	my %MSA_Pos_Species_to_Indel=%$MSA_Pos_Species_to_Indel_Reff;
	my %MSAtoIndel=%$MSAtoIndel_Reff;
	open (INDELS,$IndelInfo) || die "Can't open IndelInfo File: '$IndelInfo' $!";
	my $IndelPos="";
	my $MSA_Pos="";
	my $Length="";
	while (my $line=<INDELS>)
	{
		chomp ($line);
		if ($line=~/character number: ([0-9]+)/)
		{
			$IndelPos=$1+1; # Indel Pos start from 0
		}
		elsif ($line =~/Start position relative to MSA: ([0-9]+)/)
		{
			$MSA_Pos=$1+1;  # MSA Pos start from 0
		}
		elsif ($line=~/Found in species: (.*?) Start position relative to genome: ([0-9]+) Length: ([0-9]+)/)
		{
			my $Species=$1;
			my $length=$3;
			for (my $i=0;$i<$length;$i++)
			{
				my $tmpPosOnMSA=$MSA_Pos+$i;
				if (exists $MSA_Pos_Species_to_Indel{$tmpPosOnMSA}{$Species}){push (@{$MSA_Pos_Species_to_Indel{$tmpPosOnMSA}{$Species}},$IndelPos);}
				else {$MSA_Pos_Species_to_Indel{$tmpPosOnMSA}{$Species}=[$IndelPos];}

				if (exists $MSAtoIndel{$tmpPosOnMSA}){push (@{$MSAtoIndel{$tmpPosOnMSA}},$IndelPos);}
				else {$MSAtoIndel{$tmpPosOnMSA}=[$IndelPos];}
				print "$tmpPosOnMSA\t",$Species,"\t",join(",",@{$MSAtoIndel{$tmpPosOnMSA}}),"\n"  if ($DEBUG eq "YES"); # QA
			}
		}
		print "===========================\n"  if ($DEBUG eq "YES");
	}
	close (INDELS);
	return (\%MSA_Pos_Species_to_Indel,\%MSAtoIndel);
}
sub Read_Ancestral_Parsimony_State
{
	my $AncestralReconstructParsimony=shift;
	my $AncestralReconstructIndelParsimony_Reff=shift;
	my %AncestralReconstructIndelState=%$AncestralReconstructIndelParsimony_Reff;

	open (ANCESTRAL_INDEL_STATE,$AncestralReconstructParsimony) || die "Can't open AncestralReconstructParsimony: '$AncestralReconstructParsimony' $!";
	my $line=<ANCESTRAL_INDEL_STATE>;
	$line=<ANCESTRAL_INDEL_STATE>;
	$line=<ANCESTRAL_INDEL_STATE>;
	$line=<ANCESTRAL_INDEL_STATE>;
	$line=<ANCESTRAL_INDEL_STATE>;
	$line=<ANCESTRAL_INDEL_STATE>;
# print with MP based on the cost matrix:
#  0->0 =0
#  0->1 =2
#  1->0 =1
#  1->1 =0
#POS     Node    State

	while ($line=<ANCESTRAL_INDEL_STATE>)
	{
		chomp ($line);
		my ($POS,$Node,$State)=split(/\t/,$line);
		if ($State==0) #Char
		{
			$AncestralReconstructIndelState{$POS}{$Node}=0;
		}
		else  # Indel
		{
			$AncestralReconstructIndelState{$POS}{$Node}=1;
		}
	}
	close (ANCESTRAL_INDEL_STATE);
	return \%AncestralReconstructIndelState;
}

sub Read_Ancestral_Prob_For_Indel
{
	my $AncestralReconstructPosterior=shift;
	my $AncestralReconstructIndelPosterior_Reff=shift;
	my %AncestralReconstructIndelPosterior=%$AncestralReconstructIndelPosterior_Reff;
	print "Read_Ancestral_Prob_For_Indel: $AncestralReconstructPosterior $AncestralReconstructIndelPosterior_Reff\n=========================================================================================\n" if ($DEBUG eq "YES"); # DEBUG
	open (ANCESTRAL_INDEL_PROB,$AncestralReconstructPosterior) || die "IndelReconstruction_Wrapper.pl:Can't open AncestralReconstructPosterior: '$AncestralReconstructPosterior' $!";
	my $line=<ANCESTRAL_INDEL_PROB>;
	while ($line=<ANCESTRAL_INDEL_PROB>)
	{
		chomp ($line);
		my ($POS,$Node,$State,$Prob)=split(/\t/,$line);
		$AncestralReconstructIndelPosterior{$POS}{$Node}=$Prob;
		print "AncestralReconstructIndelPosterior{$POS}{$Node}=$Prob\n"  if ($DEBUG eq "YES"); # DEBUG
	}
	close (ANCESTRAL_INDEL_PROB);
	return \%AncestralReconstructIndelPosterior;
}

sub remove_InternalNodeName_or_BPvalues {
   my $IN_treeFile=shift;
   my $OLD_treeFile=shift;
   my $treeFileOneLine;
   open(TREEFILE,"$IN_treeFile") || die "IndelReconstruction_Wrapper.pl:remove_InternalNodeName_or_BPvalues: Can't open TREEFILE for reading '$IN_treeFile' $!";;
   while (<TREEFILE>) {
      my $line = $_;
      chomp($line);
      $treeFileOneLine .= $line;
   }
   close TREEFILE;
   my $changed = "no";
   if ($treeFileOneLine =~ m/\)N[0-9]+:/) {
	   $treeFileOneLine =~ s/\)N[0-9]+:/\):/g;    # remove internal nodes names in the BP palce
	   $changed = "yes";
   }
   if ($treeFileOneLine =~ m/\)N[0-9];/) {
	   $treeFileOneLine =~ s/\)N[0-9];/\);/g;    # remove last internal node names in the BP palce
	   $changed = "yes";
   }
   if ($treeFileOneLine =~ m/\)\d*\.?\d+\:/) {
      $treeFileOneLine =~ s/\)\d*\.?\d+\:/\)\:/g; #replace bootstrap values  which look like this: ((A:0.02,B:0.03)40:0.3);
      $changed = "yes";
   }
   if ($treeFileOneLine =~ m/\d*\.?\d+\[\d*\.?\d+\]/) {
      $treeFileOneLine =~ s/(\d*\.?\d+)\[\d*\.?\d+\]/$1/g;#replace bootstrap values  which look like this:(A:0.4,(B:0.1,C:0.1):0.3[40]);
      $changed = "yes";
   }
   if ($changed eq "yes") {
      rename $IN_treeFile, $OLD_treeFile;
      open (TREE_REMOVED,">$IN_treeFile");
      print TREE_REMOVED $treeFileOneLine."\n";
      close TREE_REMOVED;
   }
   return $changed;
}
sub uniq_array
{
	my $ReffToArray=shift;
	my %hash = ();
	foreach my $item (@$ReffToArray) {
		$hash{$item} = 1;
	}
	my @unique = sort keys(%hash);
	return \@unique;
}
sub Read_Char_Marginal_Prob
{
	my $Chars_MarginalProb_File=shift;
	my %Chars_MarginalProb=(); #Key1: MSA_Pos, Key2:Species, Key3:Char, Value:MarginalProb
	my @Nodes_Name=();
	my $MSA_Length=0;
	open (MARGINAL_PROB,$Chars_MarginalProb_File) || return "Could Not Open the MarginalProb_File: '$Chars_MarginalProb_File' $!";
	my $MSA_Pos="";
	while (my $line=<MARGINAL_PROB>)
	{
		if ($line=~/marginal probabilities at position: ([0-9]+)/)
		{
			$MSA_Pos=$1;
			$MSA_Length++;
#			print "POS:$MSA_Pos\t";
		}
		elsif ($line=~/of node: (.*?): /)
		{
			my $node=$1;
			push (@Nodes_Name,$node) if ($MSA_Pos==1);
#			print "$node\t";
			my @Chars_Prob=$line=~/p\([A-Z]+\)=[0-9\.\-]+/g;
			foreach my $Char_Prob (@Chars_Prob)
			{
				if ($Char_Prob=~/p\(([A-Z]+)\)=([0-9\.\-]+)/)
				{
					my $char=$1;
					my $prob=$2;
					$Chars_MarginalProb{$MSA_Pos}{$node}{$char}=$prob;
#					print "Chars_MarginalProb{$MSA_Pos}{$node}{$char}=$prob\n";
				}
			}
		}
	}
	close (MARGINAL_PROB);
	return (\%Chars_MarginalProb,\@Nodes_Name,$MSA_Length);
}
sub readMSA
{
	# read MSA in FASTA format return hash where key is seq name and value is sequence
	my $MSA=shift;
	my %MSA_Hash=();
	open (my $in, "<",$MSA) || die "IndelReconstruction_Wrapper:readMSA: Can't read the MSA '$MSA' $!";
	## 1.1. Read FASTA header and save it
	my $fastaLine = <$in>;
	while (defined $fastaLine) {
		chomp $fastaLine;
		my $header = substr($fastaLine,1);
		## 1.2. Read seq until next header
		$fastaLine = <$in>;
		my $seq = "";
		while ((defined $fastaLine) and
			   (substr($fastaLine,0,1) ne ">" )) {
			chomp $fastaLine;
			$seq .= $fastaLine;
			$fastaLine = <$in>;
		}
		$MSA_Hash{$header}=$seq;
	}
    # close file
	close ($in); 
	return \%MSA_Hash;
}
