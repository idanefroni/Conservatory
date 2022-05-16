use strict;

use Getopt::Long;


use FindBin qw($Bin); # www/FastML_2012/
use lib "$Bin/../bioSequence_scripts_and_constants/";
#use lib "/bioseq/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;
use POSIX;
use FindBin qw($Bin);
use File::Copy;
use File::Basename;

die "USAGE:FastML_Wrapper.pl --MSA_File <MSA_File> --seqType <aa|nuc|codon> --outDir <FULL_PATH_outDir>
Optional parameters:
  --Tree <phylogenetic tree>
  --TreeAlg <NJ | RAxML> - How to builed tree when tree not provided by user; default=NJ
  --SubMatrix <JTT | LG | mtREV | cpREV | WAG | DAYHOFF > amino acid options, the default is JTT.
              <JC_Nuc | T92 | HKY | GTR> nucleotide options, the default is JC_Nuc.
              <yang | empiriCodon> codon options, the default is yang.
  --OptimizeBL <yes | no> default: yes
  --UseGamma   <yes | no> default: yes
#  --OptAlpha   <yes | no> default: no (relevant only when UseGamma==yes)
  --Alpha      <User provide alpha>   (relevant only when UseGamma==yes) 
                                       user alpha parameter of the gamma distribution [if alpha is not given, alpha and branches will be evaluated from the data]
  --jointReconstruction <yes | no> default: yes
  --indelReconstruction <PARSIMONY|ML|BOTH> - which method is used for indel reconstruction
  --indelCutOff <Cutoff for indel vs Char> deafult =0.5
" unless (@ARGV >= 1);
my @ARGV_forPrint=@ARGV;
my %VARS=(); # FOR PROGRAM VARS
my %FORM=(); # FOR USER INPUTS

# Assign default
$FORM{MSA_File}="";
$FORM{outDir}="";
$FORM{TreeAlg}="NA";
$FORM{Tree}="NA";
$FORM{OptimizeBL}="YES";
$FORM{UseGamma}="YES";
#$FORM{OptAlpha}="NO";
$FORM{Alpha}="";
$VARS{RunNumber}="NA";
$VARS{isServer}="NO";
$FORM{JointReconstruction}="YES";

$FORM{IndelReconstructionMethod}="BOTH";
$FORM{IndelsCutoff}=0.5;
$FORM{DEBUG}="NO";
my $getoptResult = GetOptions ("MSA_File=s"=>\$FORM{MSA_File},  # = means that this parameter is required, s means string
							   "outDir=s"=>\$FORM{outDir},
							   "seqType=s"=>\$FORM{seqType},
							   "Tree:s"=>\$FORM{Tree},
							   "TreeAlg:s"=>\$FORM{TreeAlg},     # NJ | RaxML
							   "SubMatrix:s"=>\$FORM{SubMatrix},
							   "OptimizeBL:s"=>\$FORM{OptimizeBL},
							   "UseGamma:s"=>\$FORM{UseGamma},
#							   "OptAlpha:s"=>\$FORM{OptAlpha},
							   "Alpha:i"=>\$FORM{Alpha},
							   "jointReconstruction:s"=>\$FORM{JointReconstruction},
							   "indelReconstruction:s"=>\$FORM{IndelReconstructionMethod}, #Parsimony|ML
							   "RunNum:i"=>\$VARS{RunNumber},  # RELEVANT FOR SERVER ONLY
							   "isServer:s"=>\$VARS{isServer}, # RELEVANT FOR SERVER ONLY
							   "indelCutOff:f"=>\$FORM{IndelsCutoff},
							   "DEBUG:s"=>\$FORM{DEBUG} # YES | NO  
							   );

$FORM{JointReconstruction}=uc($FORM{JointReconstruction});
$FORM{UseGamma}=uc($FORM{UseGamma});
$FORM{OptimizeBL}=uc($FORM{OptimizeBL});
$FORM{TreeAlg}=uc($FORM{TreeAlg});
$FORM{DEBUG}=uc($FORM{DEBUG});
$FORM{IndelReconstructionMethod}=uc($FORM{IndelReconstructionMethod});

die "ERROR: No path for output\n" if ($FORM{outDir} eq "");
die "ERROR: MSA_File is requiered\n" if ($FORM{MSA_File} eq "");

$FORM{seqType}=lc ($FORM{seqType});
die "ERROR: seqType must be aa or nuc or codon - NOT $FORM{seqType}\n" if (($FORM{seqType} ne "aa") and ($FORM{seqType} ne "codon") and ($FORM{seqType} ne "nuc"));
unless ($FORM{outDir} =~ m/\/$/) {
	$FORM{outDir} .= "/";
}
print "outDir: $FORM{outDir}\n";
unless (-e $FORM{outDir}) {
	mkdir ($FORM{outDir});
}

if (!defined $FORM{SubMatrix}) # assign default
{
	if ($FORM{seqType} eq "aa") {$FORM{SubMatrix}="JTT"; print "SubMatrix=JTT (default)\n";}
	elsif ($FORM{seqType} eq "nuc") {$FORM{SubMatrix}="JC_Nuc"; print "SubMatrix=JC_Nuc (default)\n";}
	elsif ($FORM{seqType} eq "codon") {$FORM{SubMatrix}="yang"; print "SubMatrix=yang (default)\n";}
}

if (($FORM{Tree} ne "NA") and ($FORM{TreeAlg} ne "NA"))
{
	die "ERROR: Notice, only --Tree or --TreeAlg should be provided, not both...\n";
}
if (($FORM{Tree} ne "NA") and (!-e $FORM{Tree}))
{
	die "ERROR: The tree file '$FORM{Tree}' does not exists...\n";
}

if (($FORM{IndelsCutoff}<0) or ($FORM{IndelsCutoff}>1))
{
	die "ERROR: The --indelCutOff must be between 0 and 1...\n";
}
if (($FORM{IndelReconstructionMethod} ne "BOTH") and ($FORM{IndelReconstructionMethod} ne "PARSIMONY") and ($FORM{IndelReconstructionMethod} ne  "ML"))
{
	die "ERROR: The --indelReconstruction must be ML or PARSIMONY or BOTH Only...\n";
}
# Assign other defaults
$VARS{Aln_format}="FASTA";
$FORM{TreeAlg}="NJ" if ($FORM{TreeAlg} eq "NA");
###### here are the name of the result files.

###### tree file output in Newick format: 
$VARS{tree_newick} = "tree.newick.txt";

###### ree file output in ANCESTOR format: 
$VARS{tree_ancestor} = "tree.ancestor.txt";

###### joint sequences output file: 
$VARS{seq_joint} = "seq.joint.txt";
###### marginal sequences output file: 
$VARS{seq_marginal} = "seq.marginal.txt";
###### joint probabilities output file: 
$VARS{prob_joint} = "prob.joint.txt";
 
###### marginal probabilities output file: 
$VARS{prob_marginal} = "prob.marginal.txt";
$VARS{prob_marginal_csv} = "prob.marginal.csv";
$VARS{log_likelihood_prob_marginal_csv}="LogLikelihood_prob.margianl.csv";

# Indel Reconstructions
# Likelihood
$VARS{marginal_seq_chars_and_indel}="seq.marginal_IndelAndChars.txt";
$VARS{marginal_prob_chars_and_indel}="Ancestral_MaxMarginalProb_Char_Indel.txt";
$VARS{marginal_indel_prob}="IndelsMarginalProb.txt";
# Parsimony
$VARS{marginal_prob_chars_and_parsimony_indels}="Ancestral_MaxProb_Marginal_Char_Parsimony_Indel.txt";
$VARS{marginal_seq_chars_and_parsimony_indels}="seq.marginal_Chars_ParsimonyIndels.txt";
$VARS{parsimony_indels}="Indels.parsimony.txt";
###### JalView Ouputs
$VARS{JalViewMarginalFeaturesFile}="JalView_Features_Marginal_Prob";
$VARS{seq_marginal_JalView}="seq.marginal_NO_IndelReconstruction_JalView.$VARS{Aln_format}".".aln";
$VARS{Tree_JalView}="tree.JalView.newick";
$VARS{JalView_Marginal_Reconstruction}="JalViewMarginal_Seq_Reconstruction_NO_IndelReconstruction.html" if ($VARS{isServer} eq "YES");
$VARS{JalView_Marginal_Reconstruction}="JalViewMarginal_Seq_Reconstruction_NO_IndelReconstruction.jnlp" if ($VARS{isServer} eq "NO");
##Chars and Indels
# ML BASED
$VARS{JalViewMarginal_Chars_and_Indels_FeaturesFile}="JalView_Features_CharsAndIndels_Marginal_Prob";
$VARS{seq_marginal_Chars_and_Indels_JalView}="seq.marginal_CharsAndIndels_JalView.$VARS{Aln_format}".".aln";
if ($VARS{isServer} eq "YES")
{
	$VARS{JalView_Marginal_Chars_and_Indel_Reconstruction}="JalViewMarginal_CharsAndIndels_Reconstruction.html";
}
else
{
	$VARS{JalView_Marginal_Chars_and_Indel_Reconstruction}="JalViewMarginal_CharsAndIndels_Reconstruction.jnlp";
}
# ML CHARS PARSIMONY INDELS
$VARS{seq_marginal_chars_and_parsimony_indels_JalView}="seq.marginal_Chars_ParsimonyIndels_JalView.$VARS{Aln_format}".".aln";
$VARS{JalViewMarginal_Chars_and_Parsimony_Indels_FeaturesFile}="JalView_Features_Marginal_Prob_Chars_And_Parsimony_Indels";
if ($VARS{isServer} eq "YES")
{
	$VARS{JalView_Marginal_Chars_and_Parsimony_Indel_Reconstruction}="JalViewMarginal_Chars_And_Parsimony_Indels_Reconstruction.html";
}
else
{
	$VARS{JalView_Marginal_Chars_and_Parsimony_Indel_Reconstruction}="JalViewMarginal_Chars_And_Parsimony_Indels_Reconstruction.jnlp";
}
# Joint reconstruction
$VARS{JalViewJointAnnotationGraphFile}="JalView_Annotation_Graph_Joint_Prob";
$VARS{seq_joint_JalView}="seq.joint_JalView.$VARS{Aln_format}".".aln";
if ($VARS{isServer} eq "YES")
{
	$VARS{JalView_Joint_Reconstruction}="JalViewJoint_Reconstruction.html";
}
else
{
	$VARS{JalView_Joint_Reconstruction}="JalViewJoint_Reconstruction.jnlp";
}

###### here we set the html output file (where links to all files will be)
if ($VARS{isServer} eq "NO")
{
	$VARS{OutHtmlFile} = "output.html";
}
else
{
	$VARS{OutHtmlFile} = "output.php"; 
}

#TO DO

# Convert sequence names to num to avoid problems with RAxML and LIB

if ($VARS{isServer} eq "NO")
# Copy input files to the running dir and work on them from now on
{
	copy ($FORM{MSA_File},$FORM{outDir});
	my ($MSA_FileName,$MSA_dir)=fileparse($FORM{MSA_File});
	$FORM{MSA_File}=$FORM{outDir}.$MSA_FileName;
	print "Copy and analyse MSA: $FORM{MSA_File}\n";
	if (-e $FORM{Tree})
	{
		copy ($FORM{Tree},$FORM{outDir});
		my ($Tree_FileName,$Tree_dir)=fileparse($FORM{Tree});
		$FORM{Tree}=$FORM{outDir}.$Tree_FileName;
		print "Copy and analyse tree: $FORM{Tree}\n";
	}
}
my %SeqNamesToCode=();
my %CodeToSeqName=();
my ($SeqNamesToCode,$CodeToSeqName)=MSASeqNamesToCode($FORM{MSA_File},$FORM{outDir});
TreeNamesToCodes ($FORM{Tree},$SeqNamesToCode) if (-e $FORM{Tree});
%CodeToSeqName=%$CodeToSeqName;
%SeqNamesToCode=%$SeqNamesToCode;

################

if ($FORM{Tree} ne "NA")
{
	$VARS{UserProvideTree}="YES";
}
else
{
	$VARS{UserProvideTree}="NO";
	if ($FORM{TreeAlg} eq "RAXML")
	{
		$VARS{RAxML_Tree}="RAxML_tree.newick";
	}
}
if ($VARS{isServer} eq "YES")
{
	$VARS{All_Outputs_Zip}="FASTML_run_".$VARS{RunNumber}.".zip"; # All Outputs ZIP
	$VARS{logs_dir} = GENERAL_CONSTANTS::SERVERS_LOGS_DIR."fastml/" if ($VARS{isServer} eq "YES");
	$VARS{OutLogFile} = $VARS{logs_dir}.$VARS{RunNumber}.".log";
    ###### WWWdir is where the web=page is.
	$VARS{WWWdir} = GENERAL_CONSTANTS::FASTML_URL."results/" .$VARS{RunNumber}. "/"; #XMXMXMXMX
	$VARS{run_url} = $VARS{WWWdir}.$VARS{OutHtmlFile};
    ###### here we set the reload interval (in seconds).
	$VARS{reload_interval} = 30;
    ###### here we set the email of the server - for problems...
	$VARS{DEVELOPER_MAIL} = GENERAL_CONSTANTS::ADMIN_EMAIL;
	$VARS{UserMailFile}=$FORM{outDir}."user_email.txt";
	$VARS{DevMail} = "\"mailto:$VARS{DEVELOPER_MAIL}?subject=Fastml%20Run%20No.:%20$VARS{RunNumber}\"";
	$VARS{ContactDef} = "\n<H3><center>For assistance please <a href=$VARS{DevMail}>contact us</a> and mention this number: $VARS{RunNumber}</H3>\n";
	###### this are the name of the program.
#	$VARS{fastml} = "/bioseq/pupkoSVN/tags/fastml.v2.05/programs/fastml/fastml";                # TO DO
#	$VARS{fastml} = "/groups/pupko/haim/pupkoSVN/trunk/programs/fastml/fastml";                 # TO DO
	$VARS{fastml} = "/bioseq/FastML/fastml";
	$VARS{Indel_Reconstruction} = "/bioseq/FastML/IndelReconstruction/IndelReconstruct.pl";     # TO DO
	$VARS{RAxML} = "/bioseq/FastML/BuildRaxMLTree.pl";                                          # TO DO
	###### Send mail Global VARS
	$VARS{send_email_dir} = GENERAL_CONSTANTS::SEND_EMAIL_DIR_IBIS;
	$VARS{smtp_server} = GENERAL_CONSTANTS::SMTP_SERVER;
	$VARS{userName} = GENERAL_CONSTANTS::ADMIN_USER_NAME;
	$VARS{userPass} = GENERAL_CONSTANTS::ADMIN_PASSWORD;
	my $estimated_run_time=estimate_run_time($FORM{MSA_File},$FORM{seqType},$VARS{UserProvideTree},$FORM{UseGamma});
	
	# UPDATE STATE
	open OUTPUT, "$FORM{outDir}$VARS{OutHtmlFile}" || exit_on_error("sys_error","Can't open output page: '$FORM{outDir}$VARS{OutHtmlFile}' $!");
	my @OUTPUT=<OUTPUT>;
	close (OUTPUT);
	my $currentTime=time;
	print "CURRENT TIME:$currentTime\n";#<STDIN>;
	open (SUBMITING_TIME,">$FORM{outDir}SUBMISSION_TIME");
	print SUBMITING_TIME $currentTime;
	close (SUBMITING_TIME);
	open (STATUS,">$FORM{outDir}QUEUE_STATUS");
	print STATUS "Running";
	close (STATUS);
	open (OUTPUT, ">$FORM{outDir}$VARS{OutHtmlFile}") || exit_on_error("sys_error","Can't open output page: '$FORM{outDir}$VARS{OutHtmlFile}' $!");
	foreach my $line (@OUTPUT)
	{
		if ($line=~/QUEUED/)
		{
			$line=~s/QUEUED/RUNNING/;
			print OUTPUT $line;
		}
		elsif ($line=~/The time that passed since submitting the query is:/)
		{
			$line=~s/The time that passed since submitting the query is:/Running time is:/;
			print OUTPUT "$line";
		}
		elsif ($line=~/\<!-- HERE WILL COME ESTIMATED TIME --\>/)
		{
			print OUTPUT "<font size=\"4\">Estimated running time:<b> $estimated_run_time</b></font><br>\n";
		}
		else
		{
			print OUTPUT $line;
		}
	}
	close (OUTPUT);
}
else
{
	$VARS{logs_dir} = $FORM{outDir};
	$VARS{OutLogFile} = $FORM{outDir}."FastML_log.log";
	###### this are the name of the program
#	$VARS{fastml} = "/bioseq/pupkoSVN/tags/fastml.v2.05/programs/fastml/fastml";
	$VARS{fastml} = "$Bin/../../programs/fastml/fastml";
#	$VARS{fastml} = "/groups/pupko/haim/pupkoSVN/trunk/programs/fastml/fastml";
	$VARS{Indel_Reconstruction} = "$Bin/IndelReconstruction_Wrapper.pl";
#	$VARS{Indel_Reconstruction} = "/bioseq/FastML/IndelReconstruction/IndelReconstruct.pl";
	$VARS{RAxML} = "$Bin/BuildRaxMLTree.pl";
#	$VARS{RAxML} = "/bioseq/FastML/BuildRaxMLTree.pl";
	$VARS{DEVELOPER_MAIL} = GENERAL_CONSTANTS::ADMIN_EMAIL;
	$VARS{DevMail} = "\"mailto:$VARS{DEVELOPER_MAIL}?subject=FastML\"";

	# VALIDATION FOR COMMAND LINE
	removeEndLineExtraChars($FORM{MSA_File});
	## TO DO - ADD MORE FROM THE CGI
}

###### here we set the error definitions.

$VARS{ErrorDef} = "<font size=+3 color='red'>ERROR!  FASTML session has been terminated: </font>";
$VARS{SysErrorDef} = "<p><font size=+3 color='red'>SYSTEM ERROR - FASTML session has been terminated!</font><br><b>Please wait for a while and try to run FASTML again</b></p>\n";
print "LOG: $VARS{OutLogFile}\n";
open (LOG,">>$VARS{OutLogFile}") || exit_on_error('sys_error', "Can't open Log File: $VARS{OutLogFile} $!");
print LOG "\n\n========================================= NEW FASTML RUN STARTED ===========================================\n";
print LOG "COMMAND: perl $0 "."@ARGV_forPrint"."\n";
print LOG "FULL RUNNING PARAMETERS (INCLUDING DEFAULTS):\n--MSA_File $FORM{MSA_File} --outDir $FORM{outDir} --seqType $FORM{seqType} --Tree $FORM{Tree} --TreeAlg $FORM{TreeAlg} --SubMatrix $FORM{SubMatrix} --OptimizeBL $FORM{OptimizeBL} --UseGamma $FORM{UseGamma} --Alpha $FORM{Alpha} --jointReconstruction $FORM{JointReconstruction} --indelReconstruction $FORM{IndelReconstructionMethod} --indelCutOff $FORM{IndelsCutoff}\n";

open OUTPUT, ">>$FORM{outDir}$VARS{OutHtmlFile}" || exit_on_error("sys_error","Can't open output page: '$FORM{outDir}$VARS{OutHtmlFile}' $!");
print OUTPUT "<h4><font face=Verdana><u>Running Messages:</u></h4></font>\n";
close OUTPUT;
if (($VARS{UserProvideTree} eq "NO") and ($FORM{TreeAlg} eq "RAXML"))
{
	if ($VARS{isServer} eq "YES")
	{
		open OUTPUT, ">>$FORM{outDir}$VARS{OutHtmlFile}" || exit_on_error("sys_error","Can't open output page: '$FORM{outDir}$VARS{OutHtmlFile}' $!");
		print_message_to_output("Generating the phylogenetic tree using RAxML");
		close (OUTPUT);
	}
	RAxML();
	$FORM{Tree}="$FORM{outDir}$VARS{RAxML_Tree}";
}
RunFastML();

$FORM{Tree}="$FORM{outDir}$VARS{tree_newick}" if ($FORM{Tree} eq "NA");
# Check if there are indels to reconstruct
$VARS{AreThereIndels}=AreThereIndels($FORM{MSA_File});
if ($VARS{AreThereIndels} eq "YES")
{
	open OUTPUT, ">>$FORM{outDir}$VARS{OutHtmlFile}" || exit_on_error("sys_error","Can't open output page: '$FORM{outDir}$VARS{OutHtmlFile}' $!");
	print_message_to_output("Ancestral reconstruction of indels");
	close (OUTPUT);
	IndelReconstruction($FORM{MSA_File},$FORM{Tree},$FORM{outDir},$FORM{IndelsCutoff});
}


#### BRING BACK NAMES
#TreeCodesToNames($FORM{Tree},$CodeToSeqName);
#TreeCodesToNames("$FORM{outDir}/$VARS{tree_newick}",$CodeToSeqName);
TreeCodesToNamesShort($FORM{Tree},$CodeToSeqName);
TreeCodesToNamesShort("$FORM{outDir}/$VARS{tree_newick}",$CodeToSeqName);
MSACodesToNames($FORM{MSA_File},$CodeToSeqName);

MSACodesToNames("$FORM{outDir}/$VARS{seq_marginal}",$CodeToSeqName);
MSACodesToNames("$FORM{outDir}/$VARS{marginal_seq_chars_and_indel}",$CodeToSeqName) if ($VARS{AreThereIndels} eq "YES");
if((($FORM{IndelReconstructionMethod} eq "PARSIMONY") or ($FORM{IndelReconstructionMethod} eq "BOTH")) and ($VARS{AreThereIndels} eq "YES"))
{
#	print "MSACodesToNames($FORM{outDir}/$VARS{marginal_seq_chars_and_parsimony_indels},$CodeToSeqName);\n";<STDIN>;
	MSACodesToNames("$FORM{outDir}/$VARS{marginal_seq_chars_and_parsimony_indels}",$CodeToSeqName);
	TabDelFileCodesToNames("$FORM{outDir}/$VARS{marginal_prob_chars_and_parsimony_indels}",1,$CodeToSeqName);
	TabDelFileCodesToNames("$FORM{outDir}/$VARS{parsimony_indels}",1,$CodeToSeqName);
}

#TabDelFileCodesToNames();
TabDelFileCodesToNames("$FORM{outDir}/$VARS{marginal_prob_chars_and_indel}",1,$CodeToSeqName) if ((($FORM{IndelReconstructionMethod} eq "ML") or ($FORM{IndelReconstructionMethod} eq "BOTH")) and ($VARS{AreThereIndels} eq "YES"));

#CommaDelFileCodesToNames("$FORM{outDir}/$VARS{prob_marginal_csv}",0,$CodeToSeqName);
TabDelFileCodesToNames("$FORM{outDir}/$VARS{marginal_indel_prob}",1,$CodeToSeqName) if ((($FORM{IndelReconstructionMethod} eq "ML") || ($FORM{IndelReconstructionMethod} eq "BOTH")) and ($VARS{AreThereIndels} eq "YES"));

# print_message_to_output ("<A HREF= $VARS{marginal_prob_chars_and_indel} TARGET=_blank> The probabilities of the marginal reconstruction (including ancestral reconstruction of indels)</A>") if (($FORM{IndelReconstructionMethod} eq "ML") and ($VARS{AreThereIndels} eq "YES"));
# print_message_to_output ("<A HREF= $VARS{prob_marginal} TARGET=_blank> The probabilities of the marginal reconstruction (without ancestral reconstruction of indels)</A>");
# print_message_to_output ("<A HREF= $VARS{marginal_indel_prob} TARGET=_blank> The probabilities of the marginal reconstruction for indels</A>") if (($FORM{IndelReconstructionMethod} eq "ML") and ($VARS{AreThereIndels} eq "YES"));

if ($FORM{JointReconstruction} eq "YES"){
	MSACodesToNames("$FORM{outDir}/$VARS{seq_joint}",$CodeToSeqName);
}

AncestorFileCodesToNames("$FORM{outDir}/$VARS{tree_ancestor}",$CodeToSeqName);
# print_message_to_output ("<A HREF= $VARS{tree_ancestor} TARGET=_blank> Tree in Ancestor format</A>");
	



#### OUTPUTS
MakeJalViewOutputs();
ExtractAncestralProbPerNodePerSite("$FORM{outDir}$VARS{prob_marginal}","$FORM{outDir}$VARS{prob_marginal_csv}",$FORM{seqType});
ExtractAncestralLogLikelihoodPerNodePerSite("$FORM{outDir}$VARS{prob_marginal}","$FORM{outDir}$VARS{log_likelihood_prob_marginal_csv}",$FORM{seqType});
ZipAllOutputs() if ($VARS{isServer} eq "YES");
OrganizeOutputs();
#print to output
open OUTPUT, ">>$FORM{outDir}$VARS{OutHtmlFile}" || exit_on_error("sys_error","Can't open output page: '$FORM{outDir}$VARS{OutHtmlFile}' $!");
print OUTPUT "\n<H1><center><a name=finish> FastML has finished.</a></center></H1><br><br>\n";
print OUTPUT "<p><A HREF= $VARS{JalView_Marginal_Chars_and_Indel_Reconstruction} TARGET=_blank> The sequences of the marginal reconstruction colored by probabilities with tree (with reconstruction of indels)</A></br>\n" if ((($FORM{IndelReconstructionMethod} eq "ML") || ($FORM{IndelReconstructionMethod} eq "BOTH"))and ($VARS{AreThereIndels} eq "YES"));
print OUTPUT "<p><A HREF= $VARS{JalView_Marginal_Chars_and_Parsimony_Indel_Reconstruction} TARGET=_blank> The sequences when using marginal probabilities to reconstruct characters and maximum parsimony to reconstruct indels colored by probabilities with tree (with reconstruction of indels)</A></br>\n" if ((($FORM{IndelReconstructionMethod} eq "PARSIMONY") ||($FORM{IndelReconstructionMethod} eq "BOTH")) and ($VARS{AreThereIndels} eq "YES"));
print OUTPUT "<p><A HREF= $VARS{JalView_Marginal_Reconstruction} TARGET=_blank> The sequences of the marginal reconstruction colored by probabilities with Tree (without reconstruction of indels)</A></br>\n";
print OUTPUT "<p><A HREF= $VARS{JalView_Joint_Reconstruction} TARGET=_blank> The sequences of the joint reconstruction with probabilities and Tree (without reconstruction of indels)</A></br>\n" if ($FORM{JointReconstruction} eq "YES");

print OUTPUT "<br><br>\n";
print OUTPUT "<br><br><h4><u>Output Files:</u></h4>\n";
print_message_to_output ("<B>Download all FastML outputs in a <A HREF='$VARS{All_Outputs_Zip}'>click!</A><br><br></B>") if ($VARS{isServer} eq "YES");
print OUTPUT "<span class=\"PrintOutH\"><div id=\"PrintOutH\">Marginal reconstruction</span></div>\n";
print OUTPUT "<br><span class=\"PrintOutH_3\"><div id=\"PrintOutH_3\">Sequences<br></span></div>\n";
print_message_to_output ("<A HREF= $VARS{marginal_seq_chars_and_indel} TARGET=_blank> The sequences of the marginal reconstruction (including ancestral reconstruction of indels)</A>") if ((($FORM{IndelReconstructionMethod} eq "ML") or ($FORM{IndelReconstructionMethod} eq "BOTH")) and ($VARS{AreThereIndels} eq "YES"));
print_message_to_output ("<A HREF= $VARS{marginal_seq_chars_and_parsimony_indels} TARGET=_blank> The sequences when using marginal probability to reconstruct characters and maximum parsimony to reconstruct indels</A>") if((($FORM{IndelReconstructionMethod} eq "PARSIMONY") or ($FORM{IndelReconstructionMethod} eq "BOTH")) and ($VARS{AreThereIndels} eq "YES"));
print_message_to_output ("<A HREF= $VARS{seq_marginal} TARGET=_blank> The sequences of the marginal reconstruction (without ancestral reconstruction of indels)</A>");
if ($VARS{isServer} eq "YES")
{
	print_message_to_output ("<form enctype=\"multipart/form-data\" action=\"http://guidance.tau.ac.il/make_logo.php\" method=\"post\">\n Create a logo of the posterior probability at ancestral node ".print_make_logo_selection_box("$FORM{outDir}$VARS{prob_marginal_csv}",$FORM{seqType}));
	print_message_to_output ("<form enctype=\"multipart/form-data\" action=\"http://guidance.tau.ac.il/generate_kMostProbSeq.php\" method=\"post\">\nGenerate <select name=\"k\" id=\"k\"><option value=\"10\">10</option><option value=\"25\">25</option><option value=\"50\">50</option><option value=\"100\">100</option><option value=\"250\">250</option><option value=\"500\">500</option><option value=\"1000\">1000</option></select> most likely ancestral sequences for ancestral node ".print_make_kMostProbSeq_selection_box("$FORM{outDir}$VARS{log_likelihood_prob_marginal_csv}",$FORM{seqType}));
	print_message_to_output ("<form enctype=\"multipart/form-data\" action=\"http://guidance.tau.ac.il/SampleSeqFromProb.php\" method=\"post\">\n Sample <select name=\"NumOfSeq\" id=\"NumOfSeq\"><option value=\"10\">10</option><option value=\"25\">25</option><option value=\"50\">50</option><option value=\"100\">100</option><option value=\"250\">250</option><option value=\"500\">500</option><option value=\"1000\">1000</option><option value=\"2000\">2000</option><option value=\"3000\">3000</option><option value=\"4000\">4000</option><option value=\"5000\">5000</option></select> sequences from the posterior distribution for ancestral node ".print_SampleSeq_selection_box("$FORM{outDir}$VARS{prob_marginal_csv}",$FORM{seqType}));
}

print OUTPUT "<br><span class=\"PrintOutH_3\"><div id=\"PrintOutH_3\">Probabilities<br></span></div>\n";
print_message_to_output ("<A HREF= $VARS{marginal_prob_chars_and_indel} TARGET=_blank> The probabilities of the marginal reconstruction (including ancestral reconstruction of indels)</A>") if ((($FORM{IndelReconstructionMethod} eq "ML") or ($FORM{IndelReconstructionMethod} eq "BOTH")) and ($VARS{AreThereIndels} eq "YES"));
print_message_to_output ("<A HREF= $VARS{prob_marginal_csv} TARGET=_blank> The probabilities of the marginal reconstruction (without ancestral reconstruction of indels) - csv format</A>");
print_message_to_output ("<A HREF= $VARS{marginal_indel_prob} TARGET=_blank> The probabilities of the marginal reconstruction for indels</A>") if ((($FORM{IndelReconstructionMethod} eq "ML") or ($FORM{IndelReconstructionMethod} eq "BOTH")) and ($VARS{AreThereIndels} eq "YES"));

if ($FORM{JointReconstruction} eq "YES"){
	print OUTPUT "<span class=\"PrintOutH\"><div id=\"PrintOutH\">Joint reconstruction</span></div>\n";
	print_message_to_output ("<A HREF= $VARS{seq_joint} TARGET=_blank> The sequences of the joint reconstruction</A>");
	print_message_to_output ("<A HREF= $VARS{prob_joint} TARGET=_blank> The probabilities of the joint reconstruction</A>");
}

print OUTPUT "<span class=\"PrintOutH\"><div id=\"PrintOutH\">Phylogenetic tree</span></div>\n";
print_message_to_output ("<A HREF= $VARS{tree_newick} TARGET=_blank> Tree in Newick format</A>");
print_message_to_output ("<A HREF= $VARS{tree_ancestor} TARGET=_blank> Tree in Ancestor format</A>");
	
print OUTPUT "<br><br>\n";

print OUTPUT "\n<br><br><p><center>Please <a href= $VARS{DevMail}>report any problem</a> in case of need.</center></p>\n";
open (END_OK,">$FORM{outDir}"."FASTML_$VARS{RunNumber}".".END_OK");
close (END_OK);
close LOG;
close OUTPUT;

if ($VARS{isServer} eq "YES")
{
	sleep(25);
	
# UPDATE HEADER that the job has finished
	open OUTPUT, "<$FORM{outDir}$VARS{OutHtmlFile}";
	my @output = <OUTPUT>;
	close OUTPUT;
	open OUTPUT, ">$FORM{outDir}$VARS{OutHtmlFile}";
	foreach my $line (@output){
		if ($line =~ /FastML job status page/i){ #finds the phrase "FATML" job status page, case-insensitive
		    print OUTPUT "<H1 align=center>FastML Job Status Page - <font color='red'>FINISHED</font></h1>\n";
		}
		elsif ($line =~/Queued/)
		{
			$line=~s/Queued/Finished/;
			print OUTPUT $line;
		}
		elsif ($line =~ /REFRESH/ or $line =~ /NO-CACHE/){next;}
		else {
			print OUTPUT $line; 
		}
	}
	close OUTPUT;
	unlink ("$FORM{outDir}QUEUE_STATUS");
	&send_finish_email_to_user();
}

################################################ SUBROUTINES #############################################

#---------------------------------------------
sub OrganizeOutputs
{
	my ($Tree_FileName,$Tree_dir)=fileparse($FORM{Tree});
	my ($MSA_FileName,$MSA_dir)=fileparse($FORM{MSA_File});
	my @DebugFiles=($Tree_FileName.".ORIG_NAMES",$Tree_FileName.".ForIndelReconstruction",$Tree_FileName.".CODES",$VARS{tree_newick}.".CODES",$MSA_FileName.".ORIG_NAMES",$MSA_FileName.".CODES","SeqCodes",$VARS{seq_marginal}.".CODES",$VARS{marginal_seq_chars_and_indel}.".CODES",$VARS{marginal_prob_chars_and_indel}.".CODES",$VARS{marginal_indel_prob}.".CODES",$VARS{seq_joint}.".CODES",$VARS{tree_ancestor}.".CODES",$VARS{seq_marginal_JalView}.".CODES",$VARS{seq_marginal_Chars_and_Indels_JalView}.".CODES",$VARS{Tree_JalView}.".CODES",$VARS{seq_joint_JalView}.".CODES",$VARS{marginal_prob_chars_and_parsimony_indels}.".CODES",$VARS{marginal_seq_chars_and_parsimony_indels}.".CODES",$VARS{parsimony_indels}.".CODES",$VARS{seq_marginal_chars_and_parsimony_indels_JalView}.".CODES");
	my @JalViewFiles=($VARS{JalViewMarginalFeaturesFile},$VARS{seq_marginal_JalView},$VARS{Tree_JalView},$VARS{JalViewMarginal_Chars_and_Indels_FeaturesFile},$VARS{seq_marginal_Chars_and_Indels_JalView},$VARS{JalViewJointAnnotationGraphFile},$VARS{seq_joint_JalView},$VARS{JalView_Marginal_Chars_and_Indel_Reconstruction},$VARS{JalView_Joint_Reconstruction},$VARS{JalView_Marginal_Reconstruction},$VARS{seq_marginal_chars_and_parsimony_indels_JalView},$VARS{JalViewMarginal_Chars_and_Parsimony_Indels_FeaturesFile},$VARS{JalView_Marginal_Chars_and_Parsimony_Indel_Reconstruction});
	
	my @FilesWithParsimonyIndels=($VARS{marginal_prob_chars_and_parsimony_indels},$VARS{marginal_seq_chars_and_parsimony_indels},$VARS{parsimony_indels});
	my @FilesWithMLIndels=($VARS{marginal_seq_chars_and_indel},$VARS{marginal_prob_chars_and_indel},$VARS{marginal_indel_prob},$VARS{JalViewMarginal_Chars_and_Indels_FeaturesFile},$VARS{seq_marginal_Chars_and_Indels_JalView},$VARS{JalView_Marginal_Chars_and_Indel_Reconstruction});
	my @Logs_files=("IndelReconstruction.log","log.txt");

	if ($FORM{IndelReconstructionMethod} eq "PARSIMONY")
	{
		foreach my $file (@FilesWithMLIndels)
		{
			unlink ("$FORM{outDir}$file");
		}
	}
	if ($FORM{IndelReconstructionMethod} eq "ML")
	{
		foreach my $file (@FilesWithParsimonyIndels)
		{
			unlink ("$FORM{outDir}$file");
		}
	}
	if ($FORM{DEBUG} eq "YES")
	{
		my $ForDebugDir="$FORM{outDir}/FilesForDebug/";
		mkdir ("$ForDebugDir") if (!-e "$ForDebugDir");
		foreach my $file (@DebugFiles)
		{
			move ("$FORM{outDir}$file",$ForDebugDir) if (-e "$FORM{outDir}$file");
		}
	}
	else # delete files
	{
		foreach my $file (@DebugFiles)
		{
			unlink ("$FORM{outDir}$file");
		}
	}
	if ($VARS{isServer} eq "NO")
	{
		# JALVEIW RELATED FILES
		my $ForJalViewDir="$FORM{outDir}/FilesForJalView/";
		mkdir ($ForJalViewDir) if (!-e $ForJalViewDir);
		foreach my $file (@JalViewFiles)
		{
			move ("$FORM{outDir}$file",$ForJalViewDir) if (-e "$FORM{outDir}$file");
		}
	}
	my $test_gzip=`which gzip`;
	if (($test_gzip !~/not found/) and ($test_gzip!~/^\s+$/))
	{
		foreach my $file (@Logs_files)
		{
			system ("gzip $FORM{outDir}$file");
		}
	}
}
#---------------------------------------------
sub ZipAllOutputs
{
	
	#Marginal reconstruction
	system ("cp $VARS{marginal_seq_chars_and_indel} $FORM{outDir}sequences_of_the_marginal_reconstruction_including_indels.fas");
	system ("zip -r $FORM{outDir}$VARS{All_Outputs_Zip} sequences_of_the_marginal_reconstruction_including_indels.fas; rm -f $FORM{outDir}sequences_of_the_marginal_reconstruction_including_indels.fas");

	if(($FORM{IndelReconstructionMethod} eq "PARSIMONY") and ($VARS{AreThereIndels} eq "YES"))
	{
		system ("cp $VARS{marginal_prob_chars_and_parsimony_indels} $FORM{outDir}sequences_of_marginal_reconstruction_characters_and_parsimony_indels.fas");
		system ("zip -r $FORM{outDir}$VARS{All_Outputs_Zip} sequences_of_marginal_reconstruction_characters_and_parsimony_indels.fas; rm -f $FORM{outDir}sequences_of_marginal_reconstruction_characters_and_parsimony_indels.fas");
	}

	system ("cp $VARS{seq_marginal} $FORM{outDir}sequences_of_the_marginal_reconstruction_without_reconstruction_of_indels.fas");
	system ("zip -r $FORM{outDir}$VARS{All_Outputs_Zip} sequences_of_the_marginal_reconstruction_without_reconstruction_of_indels.fas; rm -f $FORM{outDir}sequences_of_the_marginal_reconstruction_without_reconstruction_of_indels.fas");

	
    #Probabilities
	if (($FORM{IndelReconstructionMethod} eq "ML") and ($VARS{AreThereIndels} eq "YES"))
	{
		system ("cp $VARS{marginal_prob_chars_and_indel} $FORM{outDir}Max_probabilities_of_marginal_reconstruction_with_indels.txt");
		system ("zip -r $FORM{outDir}$VARS{All_Outputs_Zip} Max_probabilities_of_marginal_reconstruction_with_indels.txt; rm -f $FORM{outDir}Max_probabilities_of_marginal_reconstruction_with_indels.txt");
	}
	
	system ("cp $VARS{prob_marginal} $FORM{outDir}probabilities_of_the_marginal_reconstruction_without_indels.txt");
	system ("zip -r $FORM{outDir}$VARS{All_Outputs_Zip} probabilities_of_the_marginal_reconstruction_without_indels.txt; rm -f $FORM{outDir}probabilities_of_the_marginal_reconstruction_without_indels.txt");

	if (($FORM{IndelReconstructionMethod} eq "ML") and ($VARS{AreThereIndels} eq "YES"))
	{
		system ("cp $VARS{marginal_indel_prob} $FORM{outDir}probabilities_of_the_marginal_reconstruction_for_indels.txt");
		system ("zip -r $FORM{outDir}$VARS{All_Outputs_Zip} probabilities_of_the_marginal_reconstruction_for_indels.txt; rm -f $FORM{outDir}probabilities_of_the_marginal_reconstruction_for_indels.txt");
	}

	# Joint reconstruction
	if ($FORM{JointReconstruction} eq "YES")
	{
		system ("cp $VARS{seq_joint} $FORM{outDir}sequences_of_the_joint_reconstruction.fas");
		system ("zip -r $FORM{outDir}$VARS{All_Outputs_Zip} sequences_of_the_joint_reconstruction.fas; rm -f $FORM{outDir}sequences_of_the_joint_reconstruction.fas");
		
		system ("cp  $VARS{prob_joint} $FORM{outDir}probabilities_of_the_joint_reconstruction.txt");
		system ("zip -r $FORM{outDir}$VARS{All_Outputs_Zip} probabilities_of_the_joint_reconstruction.txt; rm -f $FORM{outDir}probabilities_of_the_joint_reconstruction.fas");
	}

    #Phylogenetic tree
	system ("cp $VARS{tree_newick} $FORM{outDir}Tree.txt");
	system ("zip -r $FORM{outDir}$VARS{All_Outputs_Zip} Tree.txt; rm -f $FORM{outDir}Tree.txt");

	system ("cp $VARS{tree_ancestor} $FORM{outDir}Tree_in_Ancestor_format.txt");
	system ("zip -r $FORM{outDir}$VARS{All_Outputs_Zip} Tree_in_Ancestor_format.txt; rm -f $FORM{outDir}Tree_in_Ancestor_format.txt");

}

#---------------------------------------------
sub print_message_to_output{
#---------------------------------------------
    my $msg = shift;
    print OUTPUT "\n<ul><li>$msg</li></ul>\n";
}

sub print_make_logo_selection_box
{
	my $Marginal_ProbFile_CSV=shift;
	my $SeqType=shift;
	open (MARGINAL_PROB,$Marginal_ProbFile_CSV) || exit_on_error("sys_error","print_make_logo_selection_box:Can't open MARGINAL_PROB: '$Marginal_ProbFile_CSV' $!");
	my $return_str="";
	$return_str.="<input type=\"hidden\" name=\"run_num\" value=\"$VARS{RunNumber}\">\n";
    $return_str.="<input type=\"hidden\" name=\"MarginalProbFile\" value=\"$Marginal_ProbFile_CSV\"/>\n";
	$return_str.="<input type=\"hidden\" name=\"SeqType\" value=\"$SeqType\"/>\n";
    $return_str.="<select name=\"Node\" id=\"Node\">\n";
	
	my %Nodes=();
	while (my $line=<MARGINAL_PROB>)
	{
		if ($line=~/N([0-9]+)/)
		{
			$Nodes{$1}=1;
		}
	}
	for my $Node ( sort {$a<=>$b} keys %Nodes) {
		$return_str.="<option value=\"N$Node\">N$Node</option>";
	}
	$return_str.="<input type=\"submit\" value=\"Make Logo\"/>\n";
    $return_str.="</form>\n";
    close (MARGINAL_PROB);
    return $return_str;
}

sub print_make_kMostProbSeq_selection_box
{
	my $LLMarginal_ProbFile_CSV=shift;
	my $SeqType=shift;
	open (MARGINAL_PROB,$LLMarginal_ProbFile_CSV) || exit_on_error("sys_error","print_make_kMostProbSeq_selection_box:Can't open MARGINAL_PROB: '$LLMarginal_ProbFile_CSV' $!");
	my $return_str="";
	$return_str.="<input type=\"hidden\" name=\"run_num\" value=\"$VARS{RunNumber}\">\n";
    $return_str.="<input type=\"hidden\" name=\"LLMarginalProbFile\" value=\"$LLMarginal_ProbFile_CSV\"/>\n";
	$return_str.="<input type=\"hidden\" name=\"SeqType\" value=\"$SeqType\"/>\n";
    $return_str.="<select name=\"Node\" id=\"Node\">\n";
	
	my %Nodes=();
	while (my $line=<MARGINAL_PROB>)
	{
		if ($line=~/N([0-9]+)/)
		{
			$Nodes{$1}=1;
		}
	}
	for my $Node ( sort {$a<=>$b} keys %Nodes) {
		$return_str.="<option value=\"N$Node\">N$Node</option>";
	}
	$return_str.="<input type=\"submit\" value=\"Generate Sequences\"/>\n";
    $return_str.="</form>\n";
    close (MARGINAL_PROB);
    return $return_str;
}
sub send_finish_email_to_user
{
# send mail
	if (-s $VARS{UserMailFile}){
		open USR_MAIL,$VARS{UserMailFile};
		my $recipient=<USR_MAIL>;
		chomp ($recipient);
		close USR_MAIL;
        my $email_subject;
        my $HttpPath=$VARS{run_url};
		$email_subject = "'Your FastML results for run number $VARS{RunNumber} are ready'";
        my $email_message = "'Hello,\\n\\nThe results for your FastML run are ready at:\\n".$HttpPath."\\n\\nRunning Parameters:\\n";
		if (-e "$FORM{outDir}JOB_TITLE")
		{
			open (JOB_TITLE,"$FORM{outDir}JOB_TITLE");
			my $JOB_TITLE_STR=<JOB_TITLE>;
			chomp ($JOB_TITLE_STR);
			$email_message.="Job Title:$JOB_TITLE_STR\\n";
			close (JOB_TITLE);
		}
        $email_message.="\\nPlease note: the results will be kept on the server for three months.\\n\\nThanks\\nFastML Team'";

        my $msg = './sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t \''.$recipient.'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message;
#       my $msg = "ssh bioseq\@lecs \"cd $VARS{send_email_dir}; ".'perl sendEmail.pl -f \''.GENERAL_CONSTANTS::ADMIN_EMAIL.'\' -t \''.$FORM{user_mail}.'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message."\""; # TO ACTIVATE IF THE NODES IN CLUSTER FAILS TO COMMUNICATE WITH NET
        #if ($attach ne ''){$msg.=" -a $attach"; print LOG "sending $msg\n";}
		open LOG, ">>$VARS{OutLogFile}";
        print LOG "MESSAGE:$email_message\nCOMMAND:$msg\n";
        chdir $VARS{send_email_dir};
        my $email_system_return = `$msg`;
        unless ($email_system_return =~ /successfully/)    {
            print LOG "send_mail: The message was not sent successfully. system returned: $email_system_return\n";
		}
		close (LOG);
	}
}

sub ValidateInput
{
}
sub RAxML
{
	my $SubModelRaxML=uc($FORM{SubMatrix});
	$SubModelRaxML="GTRCAT" if (($FORM{SubMatrix} eq "yang") or ($FORM{SubMatrix} eq "goldman_yang") or ($FORM{SubMatrix} eq "empiriCodon") or ($FORM{SubMatrix} eq "JC_Nuc") or($FORM{SubMatrix} eq "T92") or($FORM{SubMatrix} eq "HKY") or($FORM{SubMatrix} eq "GTR"));
	my $RAxML_cmd="cd $FORM{outDir};perl $VARS{RAxML} $FORM{MSA_File} $VARS{RAxML_Tree} $FORM{outDir} $SubModelRaxML"; 
	print LOG "RAxML: running: $RAxML_cmd\n";
	print "$RAxML_cmd\n";# <STDIN>;# DEBUG
	system($RAxML_cmd);

	if (!-e "$FORM{outDir}$VARS{RAxML_Tree}")
	{
		exit_on_error("sys_error","RAxML: failed to create the output file: '$FORM{outDir}$VARS{RAxML_Tree}'");
	}
		
}
sub RunFastML
{
	my %MatrixHash = (JTT => '-mj',
					  LG => '-ml',
					  mtREV => '-mr',
					  cpREV => '-mc',
					  WAG => '-mw',
					  Dayhoff => '-md',
					  yang => '-my',
					  goldman_yang => '-mg',
					  empiriCodon => '-me',
					  JC_Nuc => '-mn',
					  JC_AA => '-ma',
					  T92 => '-mt',
					  HKY => '-mh',
					  GTR => '-mg'
					  );

	open OUTPUT, ">>$FORM{outDir}$VARS{OutHtmlFile}" || exit_on_error("sys_error","Can't open output page: '$FORM{outDir}$VARS{OutHtmlFile}' $!");
	print_message_to_output("Ancestral reconstruction of characters");
	close (OUTPUT);
	my $fastml_comm ="cd $FORM{outDir}; ".$VARS{fastml}." -s $FORM{MSA_File} $MatrixHash{$FORM{SubMatrix}} -qf";
	if ($VARS{UserProvideTree} eq "YES") { # was $treeFilePath ne "") {
		$fastml_comm = $fastml_comm." -t $FORM{Tree}";
	}
	
	elsif (($VARS{UserProvideTree} eq "NO") and ($FORM{TreeAlg} eq "RAXML")) # Tree build by RAxML
	{
		$fastml_comm = $fastml_comm." -t $FORM{outDir}$VARS{RAxML_Tree}";
	}
	if ($FORM{OptimizeBL} eq "NO") {
			$fastml_comm = $fastml_comm." -b";
	}
	if (($VARS{UserProvideTree} eq "NO") and ($FORM{TreeAlg} eq "RAXML")) # Building Tree with RAxML - NO BL OPTIMIZATION
	{
		$fastml_comm = $fastml_comm." -b";
	}
	if ($FORM{UseGamma} eq "YES") {
		$fastml_comm = $fastml_comm." -g";
		if ($FORM{Alpha} ne "") {
			$fastml_comm = $fastml_comm." -p $FORM{Alpha}";
		}
	}
	if ($FORM{JointReconstruction} eq "NO") {
		$fastml_comm = $fastml_comm." -f";
	}	
	
	$fastml_comm = $fastml_comm." > $FORM{outDir}"."fastml.std";
    print LOG "\nRunFastML: running $fastml_comm\n";
	print "$fastml_comm\n";#<STDIN>;#DEBUG
	system ($fastml_comm);
	
	# CHECK FOR ERRORS - TO DO VALIDATE FILE NAME
	my $found_error=0;
	if ((-e $FORM{outDir}."fastml.std") and ((-s $FORM{outDir}."fastml.std") > 0)){
		open STD, $FORM{outDir}."fastml.std";
        while (<STD>){
            if (/error/i){
                $found_error = 1;                
                print LOG "\nAn error was found in file \"$FORM{outDir}fastml.std\":$_ \n";
				exit_on_error("sys_error","An error was found in file \"$FORM{outDir}fastml.std\"");
                last;
            }
        }
		close STD;
	}
	if(!-e $FORM{outDir}.$VARS{seq_marginal} or -z $FORM{outDir}.$VARS{seq_marginal} or !-e $FORM{outDir}.$VARS{prob_marginal}  or -z $FORM{outDir}.$VARS{prob_marginal}){
		print LOG "\nThe file $VARS{seq_marginal} or $VARS{prob_marginal} was not created.\n";
		exit_on_error("sys_error","FASTML failed to run, one of the files $VARS{seq_marginal} or $VARS{prob_marginal} were not created...\n");
	}
}
	
sub IndelReconstruction
{
	my $MSA=shift;
	my $Tree=shift;
	my $OutDir=shift;
	my $IndelsCutoff=shift;

	## TO DO: TO THINK IF WE WANT TO RECONSTRUCT INDELS ALWAYS WITH USER/RaxML provided tree or better to do it with FASTML output and than need to change it from  $VARS{treePathOnServer}
	my $IndelReconstruction_cmd="cd $OutDir; perl $VARS{Indel_Reconstruction} --MSA_File $MSA --Tree_File $Tree --outDir $OutDir --seqType $FORM{seqType} --indelCutOff $IndelsCutoff ";
	$IndelReconstruction_cmd=$IndelReconstruction_cmd." --Debug" if ($FORM{DEBUG} eq "YES");
	$IndelReconstruction_cmd=$IndelReconstruction_cmd." > $OutDir/IndelReconstruction.log";

#	my $IndelReconstruction_cmd="cd $OutDir;perl $VARS{Indel_Reconstruction} $OutDir $MSA $Tree $OutDir/IndelsMarginalProb.txt $OutDir/prob.marginal.txt $OutDir/seq.marginal_IndelAndChars.txt $OutDir/AncestralMaxMarginalProb_Char_Indel.txt $OutDir/Indels.parsimony $OutDir/AncestralMaxProbMarginal_Char_Parsimony_Indel.txt $OutDir/seq.marginal_Chars_ParsimonyIndels.txt $FORM{seqType} $IndelsCutoff > $OutDir/IndelReconstruction.log";
	print LOG "IndelReconstruction: $IndelReconstruction_cmd\n";
	print "$IndelReconstruction_cmd\n";# <STDIN>; #DEBUG
	system($IndelReconstruction_cmd);
	# CHECK FOR ERRORS
	
}

sub MakeJalViewOutputs # TO DO
{
	# Prepare JalView Outputs
#	system ("cp $FORM{outDir}/$VARS{tree_newick}.CODES $FORM{outDir}/$VARS{Tree_JalView}");
#	TreeCodesToNamesShort("$FORM{outDir}/$VARS{Tree_JalView}",$CodeToSeqName);
#	FixTreeNamesForJalView("$FORM{outDir}/$VARS{Tree_JalView}");

	my $ans=RemoveN_FormAncestralNames("$FORM{outDir}/$VARS{seq_marginal}.CODES","$FORM{outDir}/$VARS{seq_marginal_JalView}",$VARS{Aln_format},"$FORM{outDir}/$VARS{tree_newick}.CODES","$FORM{outDir}/$VARS{Tree_JalView}");
	if ($ans ne "OK"){exit_on_error("sys_error","RemoveN_FormAncestralNames($FORM{outDir}/$VARS{seq_marginal}.CODES,$FORM{outDir}/$VARS{seq_marginal_JalView},$VARS{Aln_format},$FORM{outDir}/$VARS{tree_newick}.CODES,$FORM{outDir}/$VARS{Tree_JalView}): FAILED: $ans");}
	MSACodesToNamesShort("$FORM{outDir}/$VARS{seq_marginal_JalView}",$CodeToSeqName);
	TreeCodesToNamesShort("$FORM{outDir}/$VARS{Tree_JalView}",$CodeToSeqName);
	# Fix JalView Species names - If species name only numbers add to it ID_
    # FixAlnNamesForJalView("$FORM{outDir}/$VARS{seq_marginal_JalView}");
	
	$ans=make_Jalview_Features_MarginalProb("$FORM{outDir}/$VARS{prob_marginal}","$FORM{outDir}/$VARS{JalViewMarginalFeaturesFile}","$FORM{outDir}/$VARS{seq_marginal_JalView}");
	if ($ans ne "OK"){exit_on_error("sys_error","make_Jalview_Features_MarginalProb($FORM{outDir}/$VARS{prob_marginal},$FORM{outDir}/$VARS{JalViewMarginalFeaturesFile},$FORM{outDir}/$VARS{seq_marginal_JalView}): FAILED: $ans");}
	
	if ($VARS{isServer} eq "YES")
	{
		$ans=PrepareJalView($VARS{Tree_JalView},$VARS{seq_marginal_JalView},$VARS{WWWdir},"$FORM{outDir}/$VARS{JalView_Marginal_Reconstruction}","NA",$VARS{JalViewMarginalFeaturesFile});
		if ($ans ne "OK"){exit_on_error("sys_error","PrepareJalView($VARS{Tree_JalView},$VARS{seq_marginal_JalView},$VARS{WWWdir},$FORM{outDir}/$VARS{JalView_Marginal_Reconstruction},NA,$VARS{JalViewMarginalFeaturesFile}): FAILED: $ans");}
	}
	elsif ($VARS{isServer} eq "NO")
	{
		$ans=PrepareJalViewJNLP($VARS{Tree_JalView},$VARS{seq_marginal_JalView},"$FORM{outDir}/$VARS{JalView_Marginal_Reconstruction}","NA",$VARS{JalViewMarginalFeaturesFile});
		if ($ans ne "OK"){exit_on_error("sys_error","PrepareJalViewJNLP($VARS{Tree_JalView},$VARS{seq_marginal_JalView},$FORM{outDir}/$VARS{JalView_Marginal_Reconstruction},NA,$VARS{JalViewMarginalFeaturesFile}): FAILED: $ans");}
	}
	if ($VARS{AreThereIndels} eq "YES")
	{
		if (($FORM{IndelReconstructionMethod} eq "ML") or ($FORM{IndelReconstructionMethod} eq "BOTH"))
		{
			# Prepare JalView outputs for united indels and chars reconstruction
			$ans=RemoveN_FormAncestralNames("$FORM{outDir}/$VARS{marginal_seq_chars_and_indel}.CODES","$FORM{outDir}/$VARS{seq_marginal_Chars_and_Indels_JalView}",$VARS{Aln_format},"$FORM{outDir}/$VARS{tree_newick}.CODES","$FORM{outDir}/$VARS{Tree_JalView}");
			if ($ans ne "OK"){exit_on_error("sys_error","RemoveN_FormAncestralNames($FORM{outDir}/$VARS{marginal_seq_chars_and_indel}.CODES,$FORM{outDir}/$VARS{seq_marginal_Chars_and_Indels_JalView},$VARS{Aln_format},$FORM{outDir}/$VARS{tree_newick}.CODES,$FORM{outDir}/$VARS{Tree_JalView}): FAILED: $ans");}
			MSACodesToNamesShort("$FORM{outDir}/$VARS{seq_marginal_Chars_and_Indels_JalView}",$CodeToSeqName);
			TreeCodesToNamesShort("$FORM{outDir}/$VARS{Tree_JalView}",$CodeToSeqName);
#		FixAlnNamesForJalView("$FORM{outDir}/$VARS{seq_marginal_Chars_and_Indels_JalView}");
			$ans=make_Jalview_Features_MarginalProb_IndelsChar("$FORM{outDir}/$VARS{marginal_prob_chars_and_indel}","$FORM{outDir}/$VARS{JalViewMarginal_Chars_and_Indels_FeaturesFile}","$FORM{outDir}/$VARS{seq_marginal_Chars_and_Indels_JalView}");
			if ($ans ne "OK"){exit_on_error("sys_error","make_Jalview_Features_MarginalProb_IndelsChar($FORM{outDir}/$VARS{marginal_prob_chars_and_indel},$FORM{outDir}/$VARS{JalViewMarginal_Chars_and_Indels_FeaturesFile},$FORM{outDir}/$VARS{seq_marginal_Chars_and_Indels_JalView}): FAILED: $ans");}
			my $ans="";
			if ($VARS{isServer} eq "YES")
			{
				$ans=PrepareJalView($VARS{Tree_JalView},$VARS{seq_marginal_Chars_and_Indels_JalView},$VARS{WWWdir},"$FORM{outDir}/$VARS{JalView_Marginal_Chars_and_Indel_Reconstruction}","NA",$VARS{JalViewMarginal_Chars_and_Indels_FeaturesFile});
				if ($ans ne "OK"){exit_on_error("sys_error","PrepareJalView($VARS{Tree_JalView},$VARS{seq_marginal_Chars_and_Indels_JalView},$VARS{WWWdir},$FORM{outDir}/$VARS{JalView_Marginal_Chars_and_Indel_Reconstruction},NA,$VARS{JalViewMarginal_Chars_and_Indels_FeaturesFile}): FAILED: $ans");}
			}
			else
			{
				$ans=PrepareJalViewJNLP($VARS{Tree_JalView},$VARS{seq_marginal_Chars_and_Indels_JalView},"$FORM{outDir}/$VARS{JalView_Marginal_Chars_and_Indel_Reconstruction}","NA",$VARS{JalViewMarginal_Chars_and_Indels_FeaturesFile});
				if ($ans ne "OK"){exit_on_error("sys_error","PrepareJalViewJNLP($VARS{Tree_JalView},$VARS{seq_marginal_Chars_and_Indels_JalView},$FORM{outDir}/$VARS{JalView_Marginal_Chars_and_Indel_Reconstruction},NA,$VARS{JalViewMarginal_Chars_and_Indels_FeaturesFile}): FAILED: $ans");}
			}
		}
		if(($FORM{IndelReconstructionMethod} eq "PARSIMONY") or ($FORM{IndelReconstructionMethod} eq "BOTH"))
		{
			# Prepare JalView outputs for ML chars and Parsimony indels reconstruction
			$ans=RemoveN_FormAncestralNames("$FORM{outDir}/$VARS{marginal_seq_chars_and_parsimony_indels}.CODES","$FORM{outDir}/$VARS{seq_marginal_chars_and_parsimony_indels_JalView}",$VARS{Aln_format},"$FORM{outDir}/$VARS{tree_newick}.CODES","$FORM{outDir}/$VARS{Tree_JalView}");
			if ($ans ne "OK"){exit_on_error("sys_error","RemoveN_FormAncestralNames($FORM{outDir}/$VARS{marginal_seq_chars_and_parsimony_indels}.CODES,$FORM{outDir}/$VARS{seq_marginal_chars_and_parsimony_indels_JalView},$VARS{Aln_format},$FORM{outDir}/$VARS{tree_newick}.CODES,$FORM{outDir}/$VARS{Tree_JalView}): FAILED: $ans");}
			MSACodesToNamesShort("$FORM{outDir}/$VARS{seq_marginal_chars_and_parsimony_indels_JalView}",$CodeToSeqName);
			TreeCodesToNamesShort("$FORM{outDir}/$VARS{Tree_JalView}",$CodeToSeqName) if (!-e "$FORM{outDir}/$VARS{Tree_JalView}.CODES"); # convert if not already converted
#		FixAlnNamesForJalView("$FORM{outDir}/$VARS{seq_marginal_Chars_and_Indels_JalView}");
			$ans=make_Jalview_Features_MarginalProb_IndelsChar("$FORM{outDir}/$VARS{marginal_prob_chars_and_parsimony_indels}","$FORM{outDir}/$VARS{JalViewMarginal_Chars_and_Parsimony_Indels_FeaturesFile}","$FORM{outDir}/$VARS{seq_marginal_chars_and_parsimony_indels_JalView}");
			if ($ans ne "OK"){exit_on_error("sys_error","make_Jalview_Features_MarginalProb_IndelsChar($FORM{outDir}/$VARS{marginal_prob_chars_and_parsimony_indels},$FORM{outDir}/$VARS{JalViewMarginal_Chars_and_Parsimony_Indels_FeaturesFile},$FORM{outDir}/$VARS{seq_marginal_chars_and_parsimony_indels_JalView}): FAILED: $ans");}
			my $ans="";
			if ($VARS{isServer} eq "YES")
			{
				$ans=PrepareJalView($VARS{Tree_JalView},$VARS{seq_marginal_chars_and_parsimony_indels_JalView},$VARS{WWWdir},"$FORM{outDir}/$VARS{JalView_Marginal_Chars_and_Parsimony_Indel_Reconstruction}","NA",$VARS{JalViewMarginal_Chars_and_Parsimony_Indels_FeaturesFile});
				if ($ans ne "OK"){exit_on_error("sys_error","PrepareJalView($VARS{Tree_JalView},$VARS{seq_marginal_chars_and_parsimony_indels_JalView},$VARS{WWWdir},$FORM{outDir}/$VARS{JalView_Marginal_Chars_and_Parsimony_Indel_Reconstruction},NA,$VARS{JalViewMarginal_Chars_and_Parsimony_Indels_FeaturesFile}): FAILED: $ans");}
			}
			else
			{
				$ans=PrepareJalViewJNLP($VARS{Tree_JalView},$VARS{seq_marginal_chars_and_parsimony_indels_JalView},"$FORM{outDir}/$VARS{JalView_Marginal_Chars_and_Parsimony_Indel_Reconstruction}","NA",$VARS{JalViewMarginal_Chars_and_Parsimony_Indels_FeaturesFile});
				if ($ans ne "OK"){exit_on_error("sys_error","PrepareJalViewJNLP($VARS{Tree_JalView},$VARS{seq_marginal_chars_and_parsimony_indels_JalView},$FORM{outDir}/$VARS{JalView_Marginal_Chars_and_Parsimony_Indel_Reconstruction},NA,$VARS{JalViewMarginal_Chars_and_Parsimony_Indels_FeaturesFile}): FAILED: $ans");}
			}
		}

	}
	# FOR JOINT
	if ($FORM{JointReconstruction} eq "YES")
	{
		$ans=Make_JalView_AnnotGraph_JointProb("$FORM{outDir}/$VARS{prob_joint}","$FORM{outDir}/$VARS{JalViewJointAnnotationGraphFile}","Joint log likelihood");
		if ($ans ne "OK"){exit_on_error("sys_error","Make_JalView_AnnotGraph_JointProb($FORM{outDir}/$VARS{prob_joint},$FORM{outDir}/$VARS{JalViewJointAnnotationGraphFile},Joint log likelihood): FAILED: $ans");}
		#system ("cp $FORM{outDir}/$VARS{tree_newick} $FORM{outDir}/$VARS{Tree_JalView}");
		#FixTreeNamesForJalView("$FORM{outDir}/$VARS{Tree_JalView}");
		$ans=RemoveN_FormAncestralNames("$FORM{outDir}/$VARS{seq_joint}.CODES","$FORM{outDir}/$VARS{seq_joint_JalView}",$VARS{Aln_format},"$FORM{outDir}/$VARS{tree_newick}.CODES","$FORM{outDir}/$VARS{Tree_JalView}");
		if ($ans ne "OK"){exit_on_error("sys_error","RemoveN_FormAncestralNames($FORM{outDir}/$VARS{seq_joint}.CODES,$FORM{outDir}/$VARS{seq_joint_JalView},$VARS{Aln_format},$FORM{outDir}/$VARS{tree_newick}.CODES,$FORM{outDir}/$VARS{Tree_JalView}): Failed: $ans");}
		MSACodesToNamesShort("$FORM{outDir}/$VARS{seq_joint_JalView}",$CodeToSeqName);
		TreeCodesToNamesShort("$FORM{outDir}/$VARS{Tree_JalView}",$CodeToSeqName);

#		FixTreeNamesForJalView("$FORM{outDir}/$VARS{Tree_JalView}");
#		FixAlnNamesForJalView("$FORM{outDir}/$VARS{seq_joint_JalView}");
		if ($VARS{isServer} eq "YES")
		{
			$ans=PrepareJalView($VARS{Tree_JalView},$VARS{seq_joint_JalView},$VARS{WWWdir},"$FORM{outDir}/$VARS{JalView_Joint_Reconstruction}",$VARS{JalViewJointAnnotationGraphFile},"NA");
			if ($ans ne "OK"){exit_on_error("sys_error","PrepareJalView($VARS{Tree_JalView},$VARS{seq_joint_JalView},$VARS{WWWdir},$FORM{outDir}/$VARS{JalView_Joint_Reconstruction},$VARS{JalViewJointAnnotationGraphFile},NA): FAILED: $ans");}
		}
		elsif ($VARS{isServer} eq "NO")
		{
			$ans=PrepareJalViewJNLP($VARS{Tree_JalView},$VARS{seq_joint_JalView},"$FORM{outDir}/$VARS{JalView_Joint_Reconstruction}",$VARS{JalViewJointAnnotationGraphFile},"NA");
			if ($ans ne "OK"){exit_on_error("sys_error","PrepareJalViewJNLP($VARS{Tree_JalView},$VARS{seq_joint_JalView},$FORM{outDir}/$VARS{JalView_Joint_Reconstruction},$VARS{JalViewJointAnnotationGraphFile},NA): FAILED: $ans");}
		}
	}
}

sub FixTreeNamesForJalView
{
	my $TreeFile=shift;
	open (TREE,$TreeFile);
	my $Tree=<TREE>;
	close (TREE);
	open (TREE,">$TreeFile");
	if ($Tree !~/:/)              # With No distances
	{
		$Tree=~s/\(([0-9]+),/\(ID_$1,/g;  # between ( and ,
		$Tree=~s/,([0-9]+)\),/,ID_$1,/g;  # between , and ,
		$Tree=~s/,([0-9]+)\)/,ID_$1\)/g;  # between , and )
#		$Tree=~s/([0-9]+),/ID_$1,/g;
#		$Tree=~s/,([0-9]+)\)/,ID_$1\)/g;
		print TREE $Tree,"\n";
	}
	else # (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);
	{
#		$Tree=~s/([0-9]+):(.*?),/ID_$1:$2,/g;
#		$Tree=~s/,([0-9]+):(.*?)\)/,ID_$1:$2\)/g;
		$Tree=~s/\(([0-9]+):(.*?),/\(ID_$1:$2,/g;
		$Tree=~s/,([0-9]+):(.*?),/,ID_$1:$2,/g;
		$Tree=~s/,([0-9]+):(.*?)\)/,ID_$1:$2\)/g;
		print TREE $Tree,"\n";
	}
	close (TREE);
}
sub FixAlnNamesForJalView
{
	my $MSA=shift;
	open (MSA,$MSA) || die "FixAlnNamesForJalView: Can't open MSA for reading: '$MSA' $!";
	my @MSA=<MSA>;
	close (MSA);
	open (MSA,">$MSA") || die "FixAlnNamesForJalView: Can't open MSA for writing: '$MSA' $!";
	foreach my $line (@MSA)
	{
		if ($line=~/^>(.*)/)
		{
			chomp ($line);
			my $ID=$1;
			if ($ID=~/^([0-9]+)$/)
			{
				$ID="ID_".$1;
			}
			print MSA ">$ID\n";
		}
		else
		{
			print MSA $line;
		}
	}
	close (MSA);
}
## Handle ERRORS
# HANDLE EXIT
sub exit_on_error{
    my $which_error = shift;
    my $error_msg = shift;
    
    my $error_definition = "<font size=+2 color='red'>ERROR! FastML session has been terminated:</font><br />\n";
    my $syserror = "<font size=+1 color='red'>A SYSTEM ERROR OCCURRED!</font><br />Plesae try to run FastML again in a few minutes.<br />We apologize for the inconvenience.<br />\n";
    
    if ($which_error eq 'user_error'){
        open LOG, ">>$VARS{OutLogFile}";
        print LOG "\n\t EXIT on error:\n$error_msg\n";
		close LOG;
		if ($VARS{isServer} eq "YES")
		{
			if (-e "$FORM{outDir}$VARS{OutHtmlFile}") # OUTPUT IS OPEN
			{
				open (OUTPUT,">>$FORM{outDir}$VARS{OutHtmlFile}");
				print OUTPUT  $error_definition."$error_msg";
				close (OUTPUT);
				print "$error_definition.$error_msg\n";
			}
			else # OUTPUT WAS NOT CREATED
			{
				print "Content-type: text/html\n\n";
				print "<html>\n";
				print "<head>\n";
				print "<title>ERROR has occurred</title>\n";
				print "</head>\n";
				print "<body>\n";
				print $error_definition."$error_msg";
			}
			# print $error_msg to the screen
		}
    }
    elsif ($which_error eq 'sys_error')
    {
		open LOG, ">>$VARS{OutLogFile}";
		if ($VARS{isServer} eq "YES")
		{
			send_administrator_mail_on_error ($error_msg) if ($VARS{isServer} eq "YES");
			if (-e "$FORM{outDir}$VARS{OutHtmlFile}") # OUTPUT IS OPEN
			{
				#open LOG, ">>$FORM{outDir}$VARS{OutLogFile}";
				print LOG "\n$error_msg\n";
				open (OUTPUT,">>$FORM{outDir}$VARS{OutHtmlFile}");
				print OUTPUT $syserror;
				close OUTPUT;
				print "\n$error_msg\n";
			}
			else  # Output not open
			{
				print "Content-type: text/html\n\n";
				print "<html>\n";
				print "<head>\n";
				print "<title>ERROR has occurred</title>\n";
				print "</head>\n";
				print "<body>\n";
				print $syserror;
			}
			#print $error_msg to the log file
		}
		else
		{
			print STDERR "\n\tFASTML EXIT ON ERROR:\n$error_msg\n";
			print "\n\tFASTML EXIT ON ERROR:\n$error_msg\n";
			print LOG "\n\tFASTML EXIT ON ERROR:\n$error_msg\n";
		}
		close LOG;
	}
    close OUTPUT;
    
    if (($VARS{isServer} eq "YES") and (-e $VARS{UserMailFile}))
    {
		open (EMAIL,$VARS{UserMailFile});
		$VARS{user_email}=<EMAIL>;
		chomp($VARS{user_email});
		send_mail_on_error();
		update_output_that_run_failed();
	}
    open LOG, ">>$VARS{OutLogFile}";
    print LOG "\nExit Time: ".(BIOSEQUENCE_FUNCTIONS::printTime)."\n";
    close LOG;
    exit;
}

########################################################################################
sub send_mail_on_error
{
	my $email_subject;
	my $HttpPath=$VARS{run_url};
	$email_subject = "'Your FastML run $VARS{RunNumber} FAILED'";
	my $JOB_TITLE_STR="";
	if (-e "$FORM{outDir}JOB_TITLE")
	{
		open (JOB_TITLE,"$FORM{outDir}JOB_TITLE");
		$JOB_TITLE_STR=<JOB_TITLE>;
		chomp ($JOB_TITLE_STR);
		close (JOB_TITLE);
	}
	my $email_message = "'Hello,\\n\\nUnfortunately your FastML run (number ".$VARS{RunNumber}.") has failed.\\n";
	if ($JOB_TITLE_STR ne "")
	{
		$email_message.="Job Title:$JOB_TITLE_STR\\n";
	}
	$email_message.="Please have a look at ".$HttpPath." for further details\\n\\nSorry for the inconvenience\\nFastML Team'";
	
	my $msg = "ssh bioseq\@lecs \"cd $VARS{send_email_dir}; ".'./sendEmail.pl -f \'TAU BioSequence <bioSequence@tauex.tau.ac.il>\' -t \''.$VARS{user_email}.'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message.'"';
	#if ($attach ne ''){$msg.=" -a $attach"; print LOG "sending $msg\n";}
	open LOG, ">>$VARS{OutLogFile}";
	print LOG "MESSAGE:$email_message\nCOMMAND:$msg\n";
	chdir $VARS{send_email_dir};
	my $email_system_return = `$msg`;
	unless ($email_system_return =~ /successfully/)    {
		print LOG "send_mail: The message was not sent successfully. system returned: $email_system_return\n";
	}
	close LOG;
}

####################################################################################
sub send_administrator_mail_on_error
{
	my $message=shift;
	my $email_subject;
	chomp ($message);
	$email_subject = "'System ERROR has occurred on FastML: $VARS{run_url}'";
	my $email_message = "'Hello,\\n\\nUnfortunately a system System ERROR has occurred on FastML: $VARS{run_url}.\\nERROR: $message.'";
	my $Admin_Email=GENERAL_CONSTANTS::ADMIN_EMAIL;
	my $msg = "ssh bioseq\@lecs \"cd $VARS{send_email_dir}; ".'./sendEmail.pl -f \'bioSequence@tauex.tau.ac.il\' -t \''."bioSequence\@tauex.tau.ac.il".'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message.'"';
	#if ($attach ne ''){$msg.=" -a $attach"; print LOG "sending $msg\n";}
	print LOG "MESSAGE:$email_message\nCOMMAND:$msg\n";
	chdir $VARS{send_email_dir};
	my $email_system_return = `$msg`;
	print LOG "RESULT: $email_system_return\n";
}
sub update_output_that_run_failed
{
	open (STATUS,">$FORM{outDir}QUEUE_STATUS");
	print STATUS "Failed";
	close (STATUS);
	close OUTPUT;
    # finish the output page
    open OUTPUT, "$FORM{outDir}$VARS{OutHtmlFile}";
    my @output = <OUTPUT>;
    close OUTPUT;
    # remove the refresh commands from the output page
    open OUTPUT, ">$FORM{outDir}$VARS{OutHtmlFile}";
    foreach my $line (@output){
		if (($line=~/TTP-EQUIV="REFRESH"/) or ($line=~/CONTENT="NO-CACHE"/))
		{
			next;
		}
		elsif ($line=~/(.*)RUNNING(.*)/)
		{
			print OUTPUT $1."FAILED".$2;
		}
		elsif ($line=~/Queued/)
		{
			$line=~s/Queued/Failed/;
			print OUTPUT $line;
		}
		else {
			print OUTPUT $line;
		}
	}
    print OUTPUT "<h4 class=footer align=\"center\">Questions and comments are welcome! Please <span class=\"admin_link\"><a href=\"mailto:bioSequence\@tauex.tau.ac.il\?subject=FastML\%20Run\%20Number\%20$VARS{RunNumber}\">contact us</a></span></h4>";
    print OUTPUT "</body>\n";
    print OUTPUT "</html>\n";
    close OUTPUT;
	unlink ("$FORM{outDir}QUEUE_STATUS");
}
sub stop_reload
{    
    open LOG, ">>$FORM{outDir}$VARS{OutLogFile}";
    print LOG "\nEnd time: ".BIOSEQUENCE_FUNCTIONS::printTime();     
    close LOG;
	
    sleep ($VARS{reload_interval});
	
    open OUTPUT, "<$FORM{outDir}$VARS{OutHtmlFile}";
    my @output = <OUTPUT>;
    close OUTPUT;
    
    open OUTPUT, ">$FORM{outDir}$VARS{OutHtmlFile}";
    foreach my $line (@output){
		unless ($line =~ /REFRESH/ or $line =~ /NO-CACHE/){
			print OUTPUT $line;
		}
    }
    close OUTPUT;
}



### MAKE JalView


sub Make_JalView_AnnotGraph_JointProb
{
	my $JointProb_File=shift;
	my $OutJalviewAnnotFile=shift;
	my $Y_label=shift;           # The Y label
	my $last_x = 0;
	open (OUT,">$OutJalviewAnnotFile") || return ("Can't open outAnnotationsFile '$OutJalviewAnnotFile': $!");
	print OUT "JALVIEW_ANNOTATION\n";
	print OUT "BAR_GRAPH\t$Y_label\t";
	open (DATA_FILE,$JointProb_File) || die ("make_Jalview_AnnotationGraph: Can't open data file '$JointProb_File' $!");
	my $line=<DATA_FILE>; # header
	while ($line=<DATA_FILE>)
	{
#Joint log likelihood of position 147: -24.7398
		chomp ($line);
		if ($line=~/Joint log likelihood of position ([0-9]+): ([0-9\-e\.]+)/)
		{
			my $Pos=$1;
			my $Score=$2;
			print OUT "$Score,$Score,$Score|";
		}
	}
	print OUT "\n";
	close (OUT);
	return ("OK");
}


sub make_Jalview_Features_MarginalProb_IndelsChar
{
	my $MarginalProb_Chars_and_Indels_File=shift;
	my $outJalviewFeaturesFile=shift;
	my $aln=shift; # The naming of the seq in the features file is according to their place in the alignment

	my %NodeId_To_Seq_Num=();
	open (ALN,$aln) or return ("Can't open ALN for reading: '$aln': $!");
	my $seq_number=-1;
	while (my $line=<ALN>)
	{

		chomp ($line);
		if ($line=~/^>(.*)/)
		{
			$seq_number++;
			$NodeId_To_Seq_Num{$1}=$seq_number;
		}
	}
	close (ALN);

	open JALVIEW_FEATURES, ">$outJalviewFeaturesFile" or return ("Can't open outFeaturesFile '$outJalviewFeaturesFile': $!");
	print JALVIEW_FEATURES "Score0\t8E0152\n"; # When rounding down
	print JALVIEW_FEATURES "Score1\t8E0152\n"; 
	print JALVIEW_FEATURES "Score2\tC51B7D\n"; 
	print JALVIEW_FEATURES "Score3\tDE77AE\n"; 
	print JALVIEW_FEATURES "Score4\tF1B6DA\n"; 
	print JALVIEW_FEATURES "Score5\tFDE0EF\n";
	print JALVIEW_FEATURES "Score6\tE6F5D0\n";
	print JALVIEW_FEATURES "Score7\tB8E186\n";
	print JALVIEW_FEATURES "Score8\t7FBC41\n";
	print JALVIEW_FEATURES "Score9\t4D9221\n";
	print JALVIEW_FEATURES "Score10\t276419\n"; 
	
	print JALVIEW_FEATURES "STARTGROUP\tMarginalProb\n";
	my $pos=0;
	open (MARGINAL_PROB_CHAR_AND_INDELS,$MarginalProb_Chars_and_Indels_File) || return "make_Jalview_Color_MSA: Could Not Open the MarginalProb_Chars_and_Indels_File: '$MarginalProb_Chars_and_Indels_File' $!";
	my $line=<MARGINAL_PROB_CHAR_AND_INDELS>; # header

## SEQ ID IS NOT IDENTIFIED CORRECTLY CHANGE THE ANNOTATION TO BE ACCORDING SEQ NUMBER... (LIKE IN GUIDANCE)
## http://www.jalview.org/help/html/features/featuresFormat.html
## description	sequenceId	sequenceIndex	start	end	featureType	score (optional)
## This format allows two alternate ways of referring to a sequence, either by its text ID, or its index in an associated alignment. Normally, sequence features are associated with sequences rather than alignments, and the sequenceIndex field is given as "-1". In order to specify a sequence by its index in a particular alignment, the sequenceId should be given as "ID_NOT_SPECIFIED", otherwise the sequenceId field will be used in preference to the sequenceIndex field.

	while (my $line=<MARGINAL_PROB_CHAR_AND_INDELS>)
	{
		chomp ($line);
		my ($Pos_on_MSA,$Node,$Char,$CharProb)=split(/\t/,$line);
		if ($Node=~/^([0-9]+)$/){$Node="ID_".$1;}		
		if ($Node=~/^N([0-9]+)/)
		{
			$Node=~s/^N([0-9]+)$/$1/;
			my $Color_Class="Score".sprintf("%.0f",$CharProb*10);
#		print JALVIEW_FEATURES "$CharProb\t$Node\t-1\t",$Pos_on_MSA,"\t",$Pos_on_MSA,"\t$Color_Class\t$CharProb\n";
			print JALVIEW_FEATURES "$CharProb\tID_NOT_SPECIFIED\t$NodeId_To_Seq_Num{$Node}\t",$Pos_on_MSA,"\t",$Pos_on_MSA,"\t$Color_Class\t$CharProb\n";#FEATURE ACCORDING TO SEQ NUMBER
		}
	}
	print JALVIEW_FEATURES "ENDGROUP\tMarginalProb\n";
	close (JALVIEW_FEATURES);
	close (MARGINAL_PROB_CHAR_AND_INDELS);
	return ("OK");
}

sub make_Jalview_Features_MarginalProb
{
	my $MarginalProb_File=shift;
	my $outJalviewFeaturesFile=shift;
	my $aln=shift; # The naming of the seq in the features file is according to their place in the alignment

	my %NodeId_To_Seq_Num=();
	# print "ALN: $aln\n"; #QA
	open (ALN,$aln) or return ("Can't open ALN for reading '$aln': $!");
	my $seq_number=-1;
	while (my $line=<ALN>)
	{

		chomp ($line);
		if ($line=~/^>(.*)/)
		{
			$seq_number++;
			$NodeId_To_Seq_Num{$1}=$seq_number;
		}
	}
	close (ALN);

	
	# Global VARS
	
	open JALVIEW_FEATURES, ">$outJalviewFeaturesFile" or return ("Can't open outFeaturesFile '$outJalviewFeaturesFile': $!");
	# print "JALVIEW OUTPUT WRITE TO: $outJalviewFeaturesFile\n"; #QA
	print JALVIEW_FEATURES "Score0\t8E0152\n"; # When rounding down
	print JALVIEW_FEATURES "Score1\t8E0152\n"; 
	print JALVIEW_FEATURES "Score2\tC51B7D\n"; 
	print JALVIEW_FEATURES "Score3\tDE77AE\n"; 
	print JALVIEW_FEATURES "Score4\tF1B6DA\n"; 
	print JALVIEW_FEATURES "Score5\tFDE0EF\n";
	print JALVIEW_FEATURES "Score6\tE6F5D0\n";
	print JALVIEW_FEATURES "Score7\tB8E186\n";
	print JALVIEW_FEATURES "Score8\t7FBC41\n";
	print JALVIEW_FEATURES "Score9\t4D9221\n";
	print JALVIEW_FEATURES "Score10\t276419\n"; 
	
	print JALVIEW_FEATURES "STARTGROUP\tMarginalProb\n";
	my $pos=0;
	open (MARGINAL_PROB,$MarginalProb_File) || return "make_Jalview_Color_MSA: Could Not Open the MarginalProb_File: '$MarginalProb_File' $!";
	while (my $line=<MARGINAL_PROB>)
	{
		if ($line=~/marginal probabilities at position: ([0-9]+)/)
		{
			$pos=$1;
#			print "POS:$pos\n";
		}
		elsif ($line=~/of node: (.*?): p\([A-Z]\)=([0-9.]+)/)
		{
			my $prob=$2;
			my $Seq_ID=$1;
			if ($Seq_ID=~/^([0-9]+)$/){$Seq_ID="ID_".$1;}
			if ($Seq_ID=~/^N/) # prob only for the ancestral nodes
			{
				$Seq_ID=~s/^N([0-9]+)$/$1/;
				my $Color_Class="Score".sprintf("%.0f",$prob*10);
#			print JALVIEW_FEATURES "$prob\t$Seq_ID\t-1\t",$pos,"\t",$pos,"\t$Color_Class\t$prob\n"; # FEATURE ACCORDING TO SEQ ID (DOES PROBLEMS)
				print "[WARNNING] MISSING: NodeId_To_Seq_Num{$Seq_ID}\n" if (!defined $NodeId_To_Seq_Num{$Seq_ID});
				print JALVIEW_FEATURES "$prob\tID_NOT_SPECIFIED\t$NodeId_To_Seq_Num{$Seq_ID}\t",$pos,"\t",$pos,"\t$Color_Class\t$prob\n"; #FEATURE ACCORDING TO SEQ NUMBER
			}
		}
	}
	print JALVIEW_FEATURES "ENDGROUP\tMarginalProb\n";
	close (JALVIEW_FEATURES);
	close (MARGINAL_PROB);
	return ("OK");
}


sub RemoveN_FormAncestralNames
{
	my $Aln=shift;
	my $Aln_No_N=shift;
	my $Aln_Format=shift;
	
	my $Tree=shift;
	my $Tree_No_N=shift;
	
	open (ALN,$Aln) || return ("RemoveN_FormAncestralNames: Could not open '$Aln' $!");
	open (ALN_NO_N,">$Aln_No_N") || return ("RemoveN_FormAncestralNames: Could not open '$Aln_No_N' $!");
	
	if ($Aln_Format eq "FASTA")
	{
		while (my $line=<ALN>)
		{
			if ($line=~/^>N([0-9]+)$/)
			{
				print ALN_NO_N ">$1\n";
			}
			elsif ($line=~/^>([0-9]+)$/)
			{
				print ALN_NO_N ">ID_$1\n";
			}
			else
			{
				print ALN_NO_N $line;
			}
		}
	}
	elsif ($Aln_Format eq "CLUSTALW")
	{
		while (my $line=<ALN>)
		{
			if ($line=~/^N([0-9]+)(.*)/)
			{
				print ALN_NO_N "$1 $2\n";
			}
			else
			{
				print ALN_NO_N $line;
			}
		}
	}
	close (ALN);
	close (ALN_NO_N);
	
	open (TREE,$Tree) || return "RemoveN_FormAncestralNames: Could not open file: '$Tree' $!";
	open (TREE_NO_N,">$Tree_No_N") || return "RemoveN_FormAncestralNames: Could not open '$Tree_No_N' $!";
	
	my $Tree_String=<TREE>;
	$Tree_String=~s/\)N([0-9]+)/\)$1/g;
	print TREE_NO_N $Tree_String;
	close TREE;
	close TREE_NO_N;
	return ("OK");
}
sub PrepareJalView
{
	my $tree=shift;
	my $MSA=shift;
	my $http=shift;
	my $JalViewPage=shift;
	my $Jalview_AnnotFile=shift; # Optional, otherwise NA
	my $JalviewFeaturesFile=shift; # Optional Otherwise NA
	#print "$Jalview_AnnotFile\n$JalviewFeaturesFile\n";
	open (JALVIEW,">$JalViewPage") || return ("PrepareJalView: Can't open $JalViewPage for writing $!");
	
	print JALVIEW "<HTML>\n";
	print JALVIEW "<applet  CODEBASE=\"http://fastml.tau.ac.il/\"\n";
	print JALVIEW "CODE=\"jalview.bin.JalviewLite\" width=100% height=100%\n";
	print JALVIEW "ARCHIVE=\"jalviewApplet.jar\">\n";
	print JALVIEW "<param name=\"file\" value=\"$http"."$MSA\">\n";
	print JALVIEW "<param name=\"tree\" value=\"$http"."$tree\">\n";
	print JALVIEW "<param name=\"features\" value=\"$http".$JalviewFeaturesFile."\">\n" if ($JalviewFeaturesFile ne "NA");
	print JALVIEW "<param name=\"annotations\" value=\"$http".$Jalview_AnnotFile."\">\n" if ($Jalview_AnnotFile ne "NA");
	print JALVIEW "<param name=\"application_url\" value=\"http://www.jalview.org/services/launchApp\">\n";
	print JALVIEW "<param name=\"showbutton\" VALUE=\"false\">\n";
	print JALVIEW "<param name=\"showConservation\" VALUE=\"false\">\n";
	print JALVIEW "<param name=\"showQuality\" VALUE=\"false\">\n";
	print JALVIEW "<param name=\"showConsensus\" VALUE=\"false\">\n";
	
	print JALVIEW "<param name=\"showTreeBootstraps\"  VALUE=\"true\">\n";
	print JALVIEW "<param name=\"showTreeDistances\" VALUE=\"false\">\n";
	
	print JALVIEW "</APPLET>\n";
	print JALVIEW "</HTML>\n";
	close (JALVIEW);
	return ("OK");
	
}

sub ExtractAncestralLogLikelihoodPerNodePerSite
{
	my $FastMLProbFile=shift;
	my $OutLL=shift;
	my $SeqType=shift;
	my @AB=();
	if ($SeqType eq "nuc")
	{
		@AB=qw(A C G T);
	}
	if ($SeqType eq "aa")
	{
		@AB=qw(A C D E F G H I K L M N P Q R S T V W Y);
	}
	if ($SeqType eq "codon")
	{
		@AB=qw(AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAC TAT TCA TCC TCG TCT TGC TGG TGT TTA TTC TTG TTT);
	}
	
	open (FASTML,$FastMLProbFile) || return ("ExtractAncestralProbPerNodePerSite: Failed to open FastML Prob File: '$FastMLProbFile' $!");
	open (OUT,">$OutLL") || return ("ExtractAncestralProbPerNodePerSite: Failed to open OutLL File: '$OutLL' $!");
	print OUT "Ancestral Node,Pos,",join (",",@AB),"\n";
	my @observervedChars=();
	my $AB_Size=0;
	my %Data=();
	my $Positions=1;
	my $line=<FASTML>;
	while (($line!~/([\+]+) marginal log likelihood ([\+]+)/) and ($line)){$line=<FASTML>;} # GO TILL marginal LL SECTION
	while ($line=<FASTML>)
	{
		if ($line=~/^node/)
		{
			chomp ($line);
			@observervedChars=split(",",$line);
			$AB_Size=@observervedChars;
			$AB_Size=$AB_Size-2; # For Node and site
		}
		else
		{
			if ($line=~/([0-9]+)/)
			{
				chomp ($line);
				my @line=split(",",$line);
				for (my $i=2;$i<$AB_Size+2;$i++) # Array start with node and pos
				{
					# print "Data{$line[0]}{$line[1]}{$observervedChars[$i]}=$line[$i]\n";<STDIN>;
					$Data{$line[0]}{$line[1]}{$observervedChars[$i]}=$line[$i]; #Key1:Ancestral Node,key2:pos,key3:char value: prob
				}
				$Positions=$line[1];
			}
		}
	}
	# print "Pos:$Positions\n";
	foreach my $node (sort keys(%Data))
	{
		for (my $pos=1;$pos<=$Positions;$pos++)
		{
			print OUT "$node,$pos,";
			foreach my $Char (@AB)
			{
				if ($Char ne $AB[-1]) # NOT THE LAST POSITION
				{
					if (defined $Data{$node}{$pos}{$Char})
					{
						print OUT "$Data{$node}{$pos}{$Char},";
					}
					else
					{
						print OUT "0,"
					}
				}
				else
				{
					if (defined $Data{$node}{$pos}{$Char})
					{
						print OUT "$Data{$node}{$pos}{$Char}";
					}
					else
					{
						print OUT "0";
					}
				}
			}
			print OUT "\n";
		}
	}
}

sub ExtractAncestralProbPerNodePerSite
{
	my $FastMLProbFile=shift;
	my $OutProb=shift;
	my $SeqType=shift;
	my @AB=();
	if ($SeqType eq "nuc")
	{
		@AB=qw(A C G T);
	}
	if ($SeqType eq "aa")
	{
		@AB=qw(A C D E F G H I K L M N P Q R S T V W Y);
	}
	if ($SeqType eq "codon")
	{
		@AB=qw(AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAC TAT TCA TCC TCG TCT TGC TGG TGT TTA TTC TTG TTT);
	}

	open (FASTML,$FastMLProbFile) || return ("ExtractAncestralProbPerNodePerSite: Failed to open FastML Prob File: '$FastMLProbFile' $!");
	open (OUT,">$OutProb") || return ("ExtractAncestralProbPerNodePerSite: Failed to open OutProb File: '$OutProb' $!");
	print OUT "Ancestral Node,Pos,".join(",",@AB),"\n";
	my @observervedChars=();
	my $AB_Size=0;
	my %Data=();
	my $Positions=1;
	my $line=<FASTML>;
	while (($line!~/([\+]+) marginal probs ([\+]+)/) and ($line)){$line=<FASTML>;} # GO TILL marginal probs SECTION
	while ($line=<FASTML>)
	{
		if ($line=~/^node/)
		{
			chomp ($line);
			@observervedChars=split(",",$line);
			$AB_Size=@observervedChars;
			$AB_Size=$AB_Size-2; # For Node and site
		}
		elsif ($line=~/([\+]+) marginal log likelihood ([\+]+)/){last;} # finished
		else
		{
			if ($line=~/([0-9]+)/)
			{
				chomp ($line);
				my @line=split(",",$line);
				for (my $i=2;$i<$AB_Size+2;$i++) # Array start with node and pos
				{
					# print "Data{$line[0]}{$line[1]}{$observervedChars[$i]}=$line[$i]\n";#<STDIN>;
					$Data{$line[0]}{$line[1]}{$observervedChars[$i]}=$line[$i]; #Key1:Ancestral Node,key2:pos,key3:char value: prob
				}
				$Positions=$line[1];
			}
		}
	}
	# print "Pos:$Positions\n";
	foreach my $node (sort keys(%Data))
	{
		for (my $pos=1;$pos<=$Positions;$pos++)
		{
			print OUT "$node,$pos,";
			foreach my $Char (@AB)
			{
				if ($Char ne $AB[-1]) # NOT THE LAST POSITION
				{
					if ((defined $Data{$node}{$pos}{$Char}) and ($Data{$node}{$pos}{$Char}>0.0001))
					{
						print OUT "$Data{$node}{$pos}{$Char},";
					}
					else
					{
						print OUT "0,"
					}
				}
				else
				{
					if ((defined $Data{$node}{$pos}{$Char}) and ($Data{$node}{$pos}{$Char}>0.0001))
					{
						print OUT "$Data{$node}{$pos}{$Char}";
					}
					else
					{
						print OUT "0";
					}
				}
			}
			print OUT "\n";
		}
	}
}

sub print_SampleSeq_selection_box
{
	my $Marginal_ProbFile_CSV=shift;
	my $SeqType=shift;
	open (MARGINAL_PROB,$Marginal_ProbFile_CSV) or exit_on_error ("sys_error',print_SampleSeq_selection_box: Can't open MARGINAL_PROB: '$Marginal_ProbFile_CSV' $!");
	my $return_str="";
	$return_str.="<input type=\"hidden\" name=\"run_num\" value=\"$VARS{RunNumber}\">\n";
    $return_str.="<input type=\"hidden\" name=\"MarginalProbFile\" value=\"$Marginal_ProbFile_CSV\"/>\n";
	$return_str.="<input type=\"hidden\" name=\"SeqType\" value=\"$SeqType\"/>\n";
    $return_str.="<select name=\"Node\" id=\"Node\">\n";
	
	my %Nodes=();
	while (my $line=<MARGINAL_PROB>)
	{
		if ($line=~/N([0-9]+)/)
		{
			$Nodes{$1}=1;
		}
	}
	for my $Node ( sort {$a<=>$b} keys %Nodes) {
		$return_str.="<option value=\"N$Node\">N$Node</option>";
	}
	$return_str.="<input type=\"submit\" value=\"Sample sequences\"/>\n";
    $return_str.="</form>\n";
    close (MARGINAL_PROB);
    return $return_str;
}

sub AreThereIndels
{
	my $MSA=shift; # Fasta Format
	my $AreThereIndels="NO";
	open (MSA,$MSA) || exit_on_error ("sys_error", "AreThereIndels: Can't open the MSA file: '$MSA' $!");
	while (my $line=<MSA>)
	{
		if ($line!~/^>/)
		{
			if ($line=~/-/)
			{
				$AreThereIndels="YES";
				close (MSA);
				return ($AreThereIndels);
			}
		}
	}
	close (MSA);
	return ($AreThereIndels);
									  
}


sub removeEndLineExtraChars{
    # remove extra chars on end of lines (^M,spaces);
    my $inputFile = shift;
    my @lines;
    if (open FILE, $inputFile){
		@lines=<FILE>;
        close (FILE);
    }
    if (open (NEWFILE,">$inputFile")){
        my $line;
        foreach $line (@lines){
			#     $line=~s/(\r)$/\n/;
			$line=~s/(\s+)$//;
			print NEWFILE "$line\n";
        }
        close NEWFILE;
    }
}

sub estimate_run_time
{
	my $MSA=shift;
	my $SeqType=shift; # aa|nuc|codon
	my $IsTree=shift;
	my $IsGamma=shift;
	my %CODONS_350=(
					50 => '120',
					100 => '330',
					150 => '870',
					200 => '1230',
					250 => '1920',
					300 => '2700'
					);
	my %NUC_350=(
				 50 => '2',
				 100 => '5',
				 150 => '15',
				 200 => '35',
				 250 => '80',
				 300 => '130'
				 );
	my %AA_350=(
				50 => '2',
				100 => '6',
				150 => '25',
				200 => '60',
				250 => '180',
				300 => '300'
				);
	my $MSA_Length=0;
	my $NumOfSeq=0;
	open (MSA,$MSA);
	my $Seq="";
	while (my $line=<MSA>)
	{
		chomp ($line);
		if ($line=~/^>/)
		{
			$NumOfSeq++;
			$MSA_Length=length($Seq);
			$Seq=$Seq.$line;
		}
		else
		{
			$Seq=$line;
		}
	}
	close (MSA);
	my $Time=0;
	if ($NumOfSeq<=50)
	{
		if ($SeqType eq "aa"){$Time=$AA_350{50};}
		elsif ($SeqType eq "nuc"){$Time=$NUC_350{50};}
		elsif ($SeqType eq "codon"){$Time=$CODONS_350{50};}
	}
	if (($NumOfSeq>50)and($NumOfSeq<=100))
	{
		if ($SeqType eq "aa"){$Time=$AA_350{100};}
		elsif ($SeqType eq "nuc"){$Time=$NUC_350{100};}
		elsif ($SeqType eq "codon"){$Time=$CODONS_350{100};}
	}
	if (($NumOfSeq>100)and($NumOfSeq<=150))
	{
		if ($SeqType eq "aa"){$Time=$AA_350{150};}
		elsif ($SeqType eq "nuc"){$Time=$NUC_350{150};}
		elsif ($SeqType eq "codon"){$Time=$CODONS_350{150};}
	}
	if (($NumOfSeq>150)and($NumOfSeq<=200))
	{
		if ($SeqType eq "aa"){$Time=$AA_350{200};}
		elsif ($SeqType eq "nuc"){$Time=$NUC_350{200};}
		elsif ($SeqType eq "codon"){$Time=$CODONS_350{200};}
	}
	if (($NumOfSeq>200)and($NumOfSeq<=250)) 
	{
		if ($SeqType eq "aa"){$Time=$AA_350{2500};}
		elsif ($SeqType eq "nuc"){$Time=$NUC_350{250};}
		elsif ($SeqType eq "codon"){$Time=$CODONS_350{250};}
	}
	if (($NumOfSeq>250)and($NumOfSeq<=300))
	{
		if ($SeqType eq "aa"){$Time=$AA_350{300};}
		elsif ($SeqType eq "nuc"){$Time=$NUC_350{300};}
		elsif ($SeqType eq "codon"){$Time=$CODONS_350{300};}
	}
	if ($NumOfSeq>300)
	{
		if ($SeqType eq "aa"){$Time=$AA_350{300}*2;}
		elsif ($SeqType eq "nuc"){$Time=$NUC_350{300}*2;}
		elsif ($SeqType eq "codon"){$Time=$CODONS_350{300}*3;}
	}
	my $Seq_Factor=$MSA_Length/350;
	print "SeqFactor:$Seq_Factor\n";
	$Time=$Time*$Seq_Factor;
	
	if ($IsTree eq "NO")
	{
		$Time=$Time*3;
	}
	if ($IsGamma eq "YES")
	{
		$Time=$Time*18;
	}
	$Time=int($Time);
	print "TIME:$Time\n";#<STDIN>;
	my $RetTime=BIOSEQUENCE_FUNCTIONS::time_in_days_from_minutes($Time);
	

	my @Time=split(/:/,$RetTime);
	my $ElementInTime=@Time;
	if ($ElementInTime==2){$RetTime=$RetTime." minutes";}
	if ($ElementInTime==3){$RetTime=$RetTime." hours";}
	print "TIME:$Time\tRet:$RetTime\n";#<STDIN>;
	return ($RetTime);
	
}

sub MSASeqNamesToCode
{
	my $MSA=shift;
	my $OutDir=shift;
	
	my %SeqNamesToCode=();
	my %CodeToSeqName=();
	copy($MSA,$MSA.".ORIG_NAMES");
	open (MSA,$MSA) || exit_on_error("sys_error","MSASeqNamesToCode: Can't open MSA: '$MSA' $!");
	my @MSA=<MSA>;
	close (MSA);
	open (NEW_MSA,">$MSA") || exit_on_error("sys_error","MSASeqNamesToCode: Can't open NEW MSA: '$MSA' for writing $!");
	open (CODES,">$OutDir/SeqCodes") || exit_on_error("sys_error","MSASeqNamesToCode: Can't open SeqCodes '$OutDir/SeqCodes' $!");
	my $SeqCounter=0;
	foreach my $line (@MSA)
	{
		chomp ($line);
		if ($line=~/^>(.*)/)
		{
			my $SeqName=$1;
			my $SeqName_Short=$SeqName;
			my $SeqNameOrig=$SeqName;
			if ((length ($SeqName)>50) and ($SeqName=~/\|/))
			{
				my @SeqName=split(/\s/,$SeqName);
				$SeqName_Short=$SeqName[0];
			}
			else
			{
				$SeqName_Short=~s/[_\[\]\)\(\s\?\:\-\/\=]+/_/g;
			}
			$SeqCounter++;
			my $SeqCode="S".$SeqCounter;
			$SeqNamesToCode{$SeqNameOrig}=$SeqCode;            # Take the orig name a key
			$CodeToSeqName{$SeqCode}{'short'}=$SeqName_Short;
			$CodeToSeqName{$SeqCode}{'full'}=$SeqNameOrig;
			print NEW_MSA ">$SeqCode\n";
			print CODES "$SeqCode\t$CodeToSeqName{$SeqCode}{'short'}\t$CodeToSeqName{$SeqCode}{'full'}\n";
		}
		else
		{
			print NEW_MSA "$line\n";
		}
	}
	close (NEW_MSA);
	return (\%SeqNamesToCode,\%CodeToSeqName);
}

sub TreeNamesToCodes
{
	my $Tree=shift;
	my $SeqNamesToCode_HashRef=shift;

	copy($Tree,"$Tree".".ORIG_NAMES");
	open (TREE,$Tree) || exit_on_error("sys_error","TreeNamesToCodes: Can't open Tree: '$Tree' $!");
	my @Tree=<TREE>;
	close (TREE);
	open (TREE,">$Tree") || exit_on_error("sys_error","TreeNamesToCodes: Can't open Tree for writing: '$Tree' $!");
	foreach my $line (@Tree)
	{
		chomp $line;
		my @TreeSplit=split(/,/,$line);
		foreach my $elem (@TreeSplit)
		{
			my @elem=split(/:/,$elem);
			foreach my $elem1 (@elem)
			{
				if ($elem1 =~/([^()]+)/)
				{
					my $TaxID=$1;
					# Comment out Haim 24/2/13 beacuse after we change the name it is no longer match the line one...
					#if ((length ($Tax/ID)>50) and ($TaxID=~/\|/))
					#{
					#	my @TaxID=split(/\s/,$TaxID);
					#	$TaxID=$TaxID[0];
					#}
					#else
					#{
					#	$TaxID=~s/[_\[\]\)\(\s\?\:\-\/\=]+/_/g;
					#}
					#$TaxID=~s/[_\[\]\)\(\s\?\:\-\/\=]+/_/g;
					# print "$elem1\t$TaxID\n"; # QA
					
					if (exists $SeqNamesToCode_HashRef->{$TaxID})
					{
						$line=~s/$TaxID/$SeqNamesToCode_HashRef->{$TaxID}/;
					}
				}
			}
		}
		print TREE "$line";
	}
	print TREE "\n";
	close (TREE);
}


sub TreeCodesToNames
{
	my $Tree=shift;
	my $CodeToSeqNames_HashRef=shift;

	copy($Tree,"$Tree".".CODES");
	open (TREE,$Tree) || exit_on_error("sys_error","TreeCodesToNames:Can't open Tree: '$Tree' $!");
	my @Tree=<TREE>;
	close (TREE);
	open (TREE,">$Tree") || exit_on_error("sys_error","TreeCodesToNames:Can't open Tree for writing: '$Tree' $!"); 
	foreach my $line (@Tree)
	{
		chomp $line;
		my @TreeSplit=split(/,/,$line);
		foreach my $elem (@TreeSplit)
		{
			my @elem=split(/:/,$elem);
			foreach my $elem1 (@elem)
			{
				if ($elem1 =~/([^()]+)/)
				{
					my $TaxIDCode=$1;
					# print "$elem1\t$TaxIDCode\n"; # QA
					
					if (exists $CodeToSeqNames_HashRef->{$TaxIDCode}{'full'})
					{
						$line=~s/$TaxIDCode/$CodeToSeqNames_HashRef->{$TaxIDCode}{'full'}/;
					}
				}
			}
		}
		print TREE "$line";
	}
	print TREE "\n";
	close (TREE);
}
sub TreeCodesToNamesShort
{
	my $Tree=shift;
	my $CodeToSeqNames_HashRef=shift;

	copy($Tree,"$Tree".".CODES");
	open (TREE,$Tree) || exit_on_error("sys_error","TreeCodesToNames:Can't open Tree: '$Tree' $!");
	my @Tree=<TREE>;
	close (TREE);
	open (TREE,">$Tree") || exit_on_error("sys_error","TreeCodesToNames:Can't open Tree for writing: '$Tree' $!"); 
	foreach my $line (@Tree)
	{
		chomp $line;
		my @TreeSplit=split(/,/,$line);
		foreach my $elem (@TreeSplit)
		{
			my @elem=split(/:/,$elem);
			foreach my $elem1 (@elem)
			{
				if ($elem1 =~/([^()]+)/)
				{
					my $TaxIDCode=$1;
					# print "$elem1\t$TaxIDCode\n"; # QA
					
					if (exists $CodeToSeqNames_HashRef->{$TaxIDCode}{'short'})
					{
						$line=~s/$TaxIDCode/$CodeToSeqNames_HashRef->{$TaxIDCode}{'short'}/;
					}
				}
			}
		}
		print TREE "$line";
	}
	print TREE "\n";
	close (TREE);
}
sub MSACodesToNames
{
	my $MSA=shift;
	my $CodeToSeqNames_HashRef=shift;
	copy($MSA,"$MSA".".CODES");
	open (MSA,$MSA) ||  exit_on_error("sys_error","MSACodesToNames:Can't open MSA: '$MSA' $!");
	my @MSA=<MSA>;
	close (MSA);
	open (MSA,">$MSA") || exit_on_error("sys_error","MSACodesToNames:Can't open MSA for writing: '$MSA' $!");
	foreach my $line (@MSA)
	{
		chomp ($line);
		if ($line=~/^>(.*)/)
		{
			if (exists $CodeToSeqNames_HashRef->{$1}{'full'})
			{
				print MSA ">$CodeToSeqNames_HashRef->{$1}{'full'}\n";
			}
			else
			{
				print MSA "$line\n";
			}
		}
		else
		{
			print MSA "$line\n";
		}
	}
	close (MSA);
}
sub MSACodesToNamesShort
{
	my $MSA=shift;
	my $CodeToSeqNames_HashRef=shift;
	copy($MSA,"$MSA".".CODES");
	open (MSA,$MSA) ||  exit_on_error("sys_error","MSACodesToNames:Can't open MSA: '$MSA' $!");
	my @MSA=<MSA>;
	close (MSA);
	open (MSA,">$MSA") || exit_on_error("sys_error","MSACodesToNames:Can't open MSA for writing: '$MSA' $!");
	foreach my $line (@MSA)
	{
		chomp ($line);
		if ($line=~/^>(.*)/)
		{
			if (exists $CodeToSeqNames_HashRef->{$1}{'short'})
			{
				print MSA ">$CodeToSeqNames_HashRef->{$1}{'short'}\n";
			}
			else
			{
				print MSA "$line\n";
			}
		}
		else
		{
			print MSA "$line\n";
		}
	}
	close (MSA);
}
sub TabDelFileCodesToNames
{
	my $File=shift;
	my $Col=shift;
	my $CodeToSeqNames_HashRef=shift;
	
	copy($File,"$File".".CODES");
	open (FILE,$File) || exit_on_error("sys_error","TabDelFileCodesToNames:Can't open File: '$File' $!");
	my @File=<FILE>;
	close (FILE);
	open (FILE,">$File") || exit_on_error ("sys_error", "TabDelFileCodesToNames: Can't open File for writing: '$File' $!");
	foreach my $line (@File)
	{
		chomp ($line);
		my @line=split(/\t/,$line);
		if (exists $CodeToSeqNames_HashRef->{$line[$Col]}{'full'})
		{
			$line[$Col]=$CodeToSeqNames_HashRef->{$line[$Col]}{'full'};
		}
		my $line=join("\t",@line);
		print FILE "$line\n";
	}
	close FILE;
}

sub CommaDelFileCodesToNames
{
	my $File=shift;
	my $Col=shift;
	my $CodeToSeqNames_HashRef=shift;
	
	copy($File,"$File".".CODES");
	open (FILE,$File) || exit_on_error("sys_error","CommaDelFileCodesToNames:Can't open File: '$File' $!");
	my @File=<FILE>;
	close (FILE);
	open (FILE,">$File") || exit_on_error ("sys_error", "CommaDelFileCodesToNames: Can't open File for writing: '$File' $!");
	foreach my $line ($File)
	{
		chomp ($line);
		my @line=split(/,/,$line);
		if (exists $CodeToSeqNames_HashRef->{$line[$Col]}{'full'})
		{
			$line[$Col]=$CodeToSeqNames_HashRef->{$line[$Col]}{'full'};
		}
		print FILE join (",",@line),"\n";
	}
	close FILE;
}

sub AncestorFileCodesToNames
{
	my $File=shift;
	my $CodeToSeqNames_HashRef=shift;
	copy($File,"$File".".CODES");
	open (FILE,$File) || exit_on_error("sys_error","AncestorFileCodesToNames:Can't open File: '$File' $!");
	my @File=<FILE>;
	close (FILE);
	open (FILE,">$File") || exit_on_error ("sys_error", "AncestorFileCodesToNames: Can't open File for writing: '$File' $!");
	foreach my $line (@File)
	{
		chomp ($line);
		my @line=split(/\s+/,$line);
		foreach my $elem (@line)
		{
			if (exists $CodeToSeqNames_HashRef->{$elem}{'full'})
			{
				$line=~s/\b$elem\b/$CodeToSeqNames_HashRef->{$elem}{'full'}/;
			}
		}
		print FILE "$line\n";
	}
	close (FILE);
}
sub PrepareJalViewJNLP # Prepare JalView JNLP for Desktop app (not server)
{
	my $tree=shift;
	my $MSA=shift;
#	my $http=shift;
	my $JalViewJNLP=shift;
	my $Jalview_AnnotFile=shift; # Optional, otherwise NA
	my $JalviewFeaturesFile=shift; # Optional Otherwise NA
	#print "$Jalview_AnnotFile\n$JalviewFeaturesFile\n";
	open (JALVIEW,">$JalViewJNLP") || return ("PrepareJalViewJNLP: Can't open $JalViewJNLP for writing $!");
	print JALVIEW <<JALVIEWDESKTOP;
<!--

If you have downloaded this file after pressing "Launch Full Application" from Jalview on a web page and you don't know what to do with this file, you must install Java from http://www.java.sun.com then try opening this file again.

	JNLP generated by /jalviewServlet/launchApp
	JNLP generated from http://www.jalview.org/webstart/jalview.jnlp
Available servlet parameters (please URLEncode):
	open=<alignment file URL>
	jvm-max-heap=heap size in M or G
	features maps to	'-features'
	treeFile maps to	'-tree'
	tree maps to	'-tree'
	annotations maps to	'-annotations'
	colour maps to	'-colour'


-->
<?xml version="1.0" encoding="UTF-8"?>
<jnlp spec="1.0+" codebase="http://www.jalview.org/webstart">
	<information>
		<title>Jalview</title>
		<vendor>The Barton Group</vendor>
		<homepage href="http://www.jalview.org"/>
		<description>Jalview Multiple Alignment Editor</description>
		<description kind="short">Jalview</description>
		<icon href="JalviewLogo_big.png"/>
		<offline-allowed/>	
		<association mime-type="application-x/ext-file" extensions="fa"/>
        <association mime-type="application-x/ext-file" extensions="fasta"/>
        <association mime-type="application-x/ext-file" extensions="fastq"/>
        <association mime-type="application-x/ext-file" extensions="blc"/>
        <association mime-type="application-x/ext-file" extensions="msf"/>
        <association mime-type="application-x/ext-file" extensions="pfam"/>
        <association mime-type="application-x/ext-file" extensions="aln"/>
        <association mime-type="application-x/ext-file" extensions="pir"/>
        <association mime-type="application-x/ext-file" extensions="amsa"/>
        <association mime-type="application-x/ext-file" extensions="stk"/>
        <association mime-type="application-x/ext-file" extensions="jar"/>
	</information>
	<security>
		<all-permissions/>
	</security>
	<resources>
		<j2se version="1.6+" initial-heap-size="10M"/>
		<jar href="jalview.jar"/>
		<jar href="JGoogleAnalytics_0.3.jar"/>
		<jar href="Jmol-12.2.4.jar"/>
		<jar href="VARNAv3-9-dev.jar"/>
		<jar href="activation.jar"/>
		<jar href="apache-mime4j-0.6.jar"/>
		<jar href="axis.jar"/>
		<jar href="castor-1.1-cycle-xml.jar"/>
		<jar href="commons-codec-1.3.jar"/>
		<jar href="commons-discovery.jar"/>
		<jar href="commons-logging-1.1.1.jar"/>
		<jar href="groovy-all-1.8.2.jar"/>
		<jar href="httpclient-4.0.3.jar"/>
		<jar href="httpcore-4.0.1.jar"/>
		<jar href="httpmime-4.0.3.jar"/>
		<jar href="jaxrpc.jar"/>
		<jar href="jdas-1.0.4.jar"/>
		<jar href="jhall.jar"/>
		<jar href="jswingreader-0.3.jar"/>
		<jar href="log4j-1.2.8.jar"/>
		<jar href="mail.jar"/>
		<jar href="miglayout-4.0-swing.jar"/>
		<jar href="min-jaba-client-2.0.jar"/>
		<jar href="regex.jar"/>
		<jar href="saaj.jar"/>
		<jar href="spring-core-3.0.5.RELEASE.jar"/>
		<jar href="spring-web-3.0.5.RELEASE.jar"/>
		<jar href="vamsas-client.jar"/>
		<jar href="wsdl4j.jar"/>
		<jar href="xercesImpl.jar"/>
		<jar href="xml-apis.jar"/>
		<property name="jalview.version" value="2.8"/>
	</resources>
	<application-desc main-class="jalview.bin.Jalview">

JALVIEWDESKTOP



 	print JALVIEW "<argument>-open</argument>\n";
print JALVIEW "<argument>$MSA</argument>\n";
if ($JalviewFeaturesFile ne "NA")
{
	
	print JALVIEW "<argument>-features</argument>\n";
	print JALVIEW "<argument>$JalviewFeaturesFile</argument>\n";
}
if ($Jalview_AnnotFile ne "NA")
{
	print JALVIEW "<argument>-annotations</argument>\n";
	print JALVIEW "<argument>$Jalview_AnnotFile</argument>\n"
}
print JALVIEW "<argument>-tree</argument>\n";
print JALVIEW "<argument>$tree</argument>\n";
close (JALVIEW);
return ("OK");

}
