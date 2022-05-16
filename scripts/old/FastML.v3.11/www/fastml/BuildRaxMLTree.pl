use strict;
use FileHandle;
use Bio::SeqIO;
use Bio::AlignIO;

my $MSA=shift;           
my $OutTree=shift;       
my $WorkingDir=shift;

my $Model=shift; #Available AA substitution models: DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, LG, MTART, MTZOA, PMB, HIVB, HIVW, JTTDCMUT, FLU, GTR
                 #NUC:                              GTRCAT

my $MSA_Name=$MSA; # IF WITHOUT PATH
if ($MSA=~/([^\/]+)$/){$MSA_Name=$1;} # NAME WITHOUT PATH

my $OutTree_Suffix=$OutTree; # IF WITHOUT PATH
if ($OutTree=~/([^\/]+)$/){$OutTree_Suffix=$1;} # NAME WITHOUT PATH

$WorkingDir=$WorkingDir."/" if ($WorkingDir!~/\//);
my $Codes2NameIndex=$WorkingDir."$MSA_Name"."Codes2NamesIndex.txt";
my $CodedMSA=$WorkingDir."/$MSA_Name".".coded.aln";
my $CodedMSAPhylip=$WorkingDir."$MSA_Name".".coded.Phylip";
# Convert Names to numbers
my $ans=name2codeFastaFrom1("$MSA",$Codes2NameIndex,$CodedMSA);
#if ($ans ne "ok") {exit_on_error}
# Convert To Phylip
convertMsaFormat($CodedMSA,$CodedMSAPhylip,"fasta","phylip");
#my $convert_cmd="readseq -a -f12 $CodedMSA > $CodedMSAPhylip";
#system ($convert_cmd);
# Run RaxML
$Model="PROTCAT".$Model if ($Model ne "GTRCAT");
my $RaxML_cmd="cd $WorkingDir;raxmlHPC -s $CodedMSAPhylip -n $OutTree_Suffix"." -m $Model";
print "$RaxML_cmd\n";
system ($RaxML_cmd);
# Bring Back names to tree
my $RaxMLTree="RAxML_bestTree.$OutTree_Suffix";
code2nameTree($Codes2NameIndex,$WorkingDir.$RaxMLTree,$WorkingDir."$OutTree_Suffix");


sub name2codeFastaFrom1 {
####################################################################################################################
# Convert the names in a fasta file to numbers, and creates a code file with the names and the codes (running number)
###################################################################################################################
        my $in_fileName = shift;
        my $code_fileName = shift;
        my $out_fileName = shift;
        my $counter_offset=shift; # optional

        my $in_file = Bio::SeqIO->new(-file => $in_fileName , '-format' => 'Fasta');
        my $code_file = new FileHandle(">$code_fileName") or return ("Can't write to $code_fileName $!");
        my $out_file = new FileHandle(">$out_fileName") or return ("Can't write to $out_fileName");
        $counter_offset=1 if (!defined $counter_offset);
        $counter_offset=1 if ($counter_offset==0);
        my $counter = $counter_offset;
        my $i;

        while ( my $seqObj = $in_file->next_seq() ) {
                my $name = $seqObj->display_id();
                $name.= " ".$seqObj->desc()   if ($seqObj->desc());
                print $code_file "$name\t$counter\n";
                my $seq = $seqObj->seq();
                print $out_file ">$counter\n";
                for($i=0;$i<length($seq);$i+=60){
                        print $out_file substr($seq,$i,60) . "\n";
                }
                if($i<length($seq)){
                        print $out_file substr($seq,$i,length($seq)-$i);
                }
                print $out_file "\n";
                $counter++;
        }
        $out_file->close();
        $in_file->close();
        $code_file->close();
        return "ok";
}

sub code2nameTree
{
###############################################################################################################
# Works together (or rather after) the script names2codeFasta.pl. Takes a tree created based on 
# a fasta file with codes, and reverts the codes to the names. Required input is a code file which is created by
# names2codeFasta.pl
# ** very useful for working with all phyml and such, since these programs chop the name to 10 chars
###############################################################################################################


# die "Usage: code2name.pl CODE_FILE TREE_FILE NEW_FILE NAME_LENGTH" if (scalar(@ARGV) < 3);
	my $nameLength = "NA";
	my $code2nameFile = shift;
	my $treeFile = shift;
	my $newFile = shift;
	
	$nameLength = shift;
	if (!defined $nameLength) {
		$nameLength = 30;
	}
	
	
	
	my %names2code;
	my @fields;
	

	open FH, "<$code2nameFile";
	while (my $line=<FH>){
		$line =~ /(.+)\t(\d+)/;
		my $code = $2;
		my $name = $1;
		$name =~ s/[\[\]\,\:\;\(\)]/_/g; #remove characters that are newick format associated
		if ($name =~ m/(.*\|.{$nameLength})/) {
			$name = $1;
		}
		$names2code{$code}=$name;
		print "$code $name\n";
	}
	
	close FH;
	
	open TREE, "<$treeFile";
	open NEWTREE, ">$newFile";
	
	my $full_tree = "";
	my $line2;
	while ($line2 = <TREE>){ # this assumes there are bootstrap values on the input tree
		chomp $line2;
		$full_tree.=$line2;
		
	}
	
	@fields = split(/:/, $full_tree);
	
	foreach my $field (@fields) {
		if ($field =~ /[\,\(](\d+)$/) { # a leaf comes either after a "(" or a ","
			$field =~ s/(\d+)$/$names2code{$1}/;
		}
	
		if ($field !~/;$/) {print NEWTREE "$field:";}
		else {print NEWTREE "$field";}  # Last One
	}
	
	print NEWTREE "\n";	
}

sub convertMsaFormat
{
	my $inFile=shift;
	my $outFile=shift;
	my $inFormat=shift;
	my $outFormat=shift;

	#die "usage: convertMsaFormat.pl <inFile> <outFile> <inFormat> <outFormat>\n"
	
	print "inFile = '$inFile'\n";
	print "outFile = '$outFile'\n";
	print "inFormat = '$inFormat'\n";
	print "outFormat = '$outFormat'\n";
	my $in  = Bio::AlignIO->new( '-format' => $inFormat , -file => $inFile);
	my $out  = Bio::AlignIO->new( '-format' => $outFormat , -file => ">$outFile");
	
	my ($alignObj, $seqStr, $trans);
	while ($alignObj = $in->next_aln()) {
		$alignObj->verbose(1);
		# Otherwise, bioperl adds sequence start/stop values, causing problems 
		# with clustal/bali_score
		$alignObj->set_displayname_flat();
		$out->write_aln($alignObj);
	}
}
