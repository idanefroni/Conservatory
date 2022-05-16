#!/usr/bin/perl

package GENERAL_CONSTANTS; #don't forget: a package must end with a return value (1; in the end)!!!!!

# constants to use when sending e-mails using the server admin's email address.
use constant ADMIN_EMAIL =>                 "TAU BioSequence \<bioSequence\@tauex.tau.ac.il\>";
use constant ADMIN_USER_NAME =>             "";
use constant ADMIN_PASSWORD =>              "";
#use constant SMTP_SERVER =>                 "";
use constant SMTP_SERVER =>                 "";

# the name of the list of all running processes
use constant QUEUING_JOBS =>                "/bioseq/bioSequence_scripts_and_constants/queuing_jobs.list";
use constant RUNNING_JOBS =>                "/bioseq/bioSequence_scripts_and_constants/running_jobs.list";
use constant SUBMITTED_JOBS =>              "/bioseq/bioSequence_scripts_and_constants/submitted_jobs.list";
use constant JOBS_ON_BIOSEQ_NODE =>         "/bioseq/bioSequence_scripts_and_constants/jobs_on_bioc.01_node.list";
use constant JOBS_WAITING_BIOSEQ_NODE =>         "/bioseq/bioSequence_scripts_and_constants/jobs_waiting_bioc.01_node.list";
use constant CONSURF_RUNNING_JOBS =>        "/bioseq/bioSequence_scripts_and_constants/consurf_running_jobs.list";
use constant SELECTON_RUNNING_JOBS =>       "/bioseq/bioSequence_scripts_and_constants/selecton_running_jobs.list";
use constant CONSEQ_RUNNING_JOBS =>       "/bioseq/bioSequence_scripts_and_constants/conseq_running_jobs.list";
use constant PEPITOPE_RUNNING_JOBS =>       "/bioseq/bioSequence_scripts_and_constants/pepitope_running_jobs.list";

# Databases urls
use constant PROTEOPEDIA =>       "http://proteopedia.org/wiki/index.php/";
use constant PDB_DB =>            "http://www.rcsb.org/pdb/explore/explore.do?structureId=";
use constant RCSB_WGET=>	  "wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/";
use constant RCSB => 		  "http://www.rcsb.org/";
use constant PISA_WGET => 	  "wget http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/multimer.pdb?";


# CGIs paths
use constant CONSURF_CGI_DIR =>             "/var/www/cgi-bin/ConSurf"; 

#general paths
use constant SERVERS_RESULTS_DIR =>         "/bioseq/data/results/";
use constant SERVERS_LOGS_DIR =>            "/bioseq/data/logs/";
#use constant SEND_EMAIL_DIR =>              "/db1/Local/src/sendEmail"; # path on biocluster
use constant SEND_EMAIL_DIR =>              "/bioseq/bioSequence_scripts_and_constants/sendEmail";
use constant SEND_EMAIL_DIR_IBIS =>         "/bioseq/bioSequence_scripts_and_constants/sendEmail";  # path on ibis
use constant DAEMON_LOG_FILE =>             "/bioseq/bioSequence_scripts_and_constants/daemon.log";
use constant UPDATE_RUN_TIME_LOG_FILE =>    "/bioseq/bioSequence_scripts_and_constants/update_runTime.log";
use constant CONSURF_CGI =>                 "/var/www/cgi-bin/ConSurf";  #on ibis
use constant BIOSEQ_TEMP =>                 "/bioseq/temp/";

# servers urls:
use constant SELECTON_URL =>     "http://selecton.tau.ac.il";
use constant CONSEQ_URL =>       "http://conseq.tau.ac.il/";
use constant CONSURF_URL =>      "http://consurf.tau.ac.il/";
use constant NEW_CONSURF_URL =>  "http://consurf.tau.ac.il/"; #"http://consurftest.tau.ac.il/";
use constant EPITOPIA_URL =>     "http://epitopia.tau.ac.il/";
use constant PEPITOPE_URL =>     "http://pepitope.tau.ac.il/";
use constant QMF_URL =>          "http://quasimotifinder.tau.ac.il/";
use constant PATCHFINDER_URL =>  "http://patchfinder.tau.ac.il/";
#use constant FASTML_URL =>       "http://ibis.tau.ac.il/fastml/";
use constant FASTML_URL =>       "http://fastml.tau.ac.il/";
use constant RECONST_URL =>      "http://fastml.tau.ac.il/reconst/";
use constant GAIN_LOSS_URL =>    "http://gloome.tau.ac.il/";
use constant CONSURF_DB_URL =>   "http://consurfdb.tau.ac.il/";
#use constant GILAD_SERVER_URL => "http://consurftest.tau.ac.il/Gilad/";
use constant GILAD_SERVER_URL => "http://mud.tau.ac.il/";
use constant MCPep_URL => "http://bental.tau.ac.il/MCPep/";
use constant GUIDANCE_URL =>     "http://guidance.tau.ac.il/";
use constant GUIDANCE_INDELS_URL =>     "http://guidance.tau.ac.il/indels/";
use constant SPECBOOST_URL =>    "http://bental.tau.ac.il/specBoost/";
use constant PROMAYA_URL =>    "http://bental.tau.ac.il/ProMaya/";
use constant HOMOLOGY_SEARCH_URL =>    "http://fastml.tau.ac.il/HomologySearch/";
use constant COPAP_URL =>    "http://copap.tau.ac.il/";

#servers logs:
use constant CONSURF_LOG =>      "/bioseq/ConSurf_old/consurf.log"; 
use constant CONSURF_NEW_LOG =>  "/bioseq/ConSurf/consurf.log";
use constant SELECTON_LOG =>     "/bioseq/Selecton/selecton.log"; 
use constant EPITOPIA_LOG =>     "/bioseq/epitopia/epitopia.log"; 
use constant CONSEQ_LOG =>       "/bioseq/ConSeq/conseq.log";  
use constant PEPITOPE_LOG =>     "/bioseq/pepitope/pepitope.log";
use constant RECONST_LOG =>      "/bioseq/ReConst_Server/reconst.log";
use constant MCPep_LOG => 	 "/bioseq/MCPep/mcpep.log";
use constant Guidance_LOG =>     "/bioseq/Guidance/guidance.log";
use constant Guidance_Indels_LOG =>     "/bioseq/GuidanceIndels/guidance_Indels.log";
use constant MuD_LOG =>          "/bioseq/Gilad_Server/MuD.log";
use constant FASTML_LOG =>       "/bioseq/FastML/fastml.log";
use constant SPECBOOST_LOG =>    "/bioseq/specBoost/specBoost.log";
use constant GAIN_LOSS_LOG =>    "/bioseq/GainLoss/GainLoss.log";
use constant PROMAYA_LOG =>      "/bioseq/ProMaya/ProMaya.log";
use constant COPAP_LOG =>      "/bioseq/CoPAP/CoPAP.log";

#servers results urls:
# servers urls:
use constant SELECTON_RESULTS_URL =>        SELECTON_URL."/results/";

#external databases
#use constant PQS=>                  "/bioseq/data/results/PQS/";
use constant PQS=>                  "/biodb/PQS/";
use constant PDB_DIVIDED =>         "/biodb/PDB/data/structures/divided/pdb/";
use constant SWISSPROT_DB =>        "/biodb/BLAST/Proteins/swissprot";
use constant UNIPROT_DB =>          "/biodb/BLAST/Proteins/uniprot";
use constant CLEAN_UNIPROT_DB =>    "/biodb/BLAST/Proteins/clean_uniprot";
use constant UNIREF90_DB =>         "/biodb/BLAST/Proteins/uniref90";#"/groups/bioseq.home/HAIM/UNIREF90/uniref90";
use constant PDBAA_NCBI=>           "/biodb/BLAST/Proteins/pdbaa";
use constant CULLED_PDB =>          "/groups/bioseq.home/HAIM/PDBAA/pdbaaent";  # TO CHANGE TO: /biodb/BLAST/dunbrack.fccc.edu/Guoli/culledpdb/pdbaaent_dun
use constant PDB_DUNBRACK =>        "/groups/bioseq.home/HAIM/PDBAA/pdbaa";     # TO CHANGE TO: /biodb/BLAST/dunbrack.fccc.edu/Guoli/culledpdb/pdbaa_dun
use constant NR_PROT_DB =>          "/biodb/BLAST/Proteins/nr";
use constant NR_NUC_DB =>           "/biodb/BLAST/Nucleotides/nt";
use constant UNIPROT_DAT_INDEX =>   "/bioseq/data/results/GB_CDS/uniprot.dat.bp_index";
use constant PDB_TO_UNIPROT =>      "/bioseq/data/results/PDB_to_UNIPROT/idmapping_PDB_UNIPROTKB.dat";#"/biodb/idmapping_PDB_UNIPROTKB.dat";
use constant PDB_TO_UNIPROT_test =>      "/biodb/idmapping_PDB_UNIPROTKB.dat";
#internal databases
use constant EPITOPIA_DATA => "/bioseq/epitopia/data";

#external programs
use constant BLASTALL =>            "/opt/bio/ncbi/bin/blastall"; #"/opt/Bio/ncbi/bin/blastall"; # on the lecs
use constant BLASTPGP =>            "blastpgp"; # "/opt/Bio/ncbi/bin/blastpgp"; # on the lecs
use constant CS_BLAST =>            "/share/apps/csblast-2.1.0-linux64/csblast_static"; # on the lecs
use constant MUSCLE_LECS =>         "/share/apps/bin/muscle";  # on the lecs
use constant MUSCLE =>              "/usr/local/bin/muscle";  # on the biocluster
use constant MUSCLE_3_6 =>          "/bioseq/Programs/muscle_3.6_from_BIOCLUSTER/muscle3.6/muscle";  # for servers who came from biocluster (Selecton?, old ConSurf, ConSeq)
use constant CLUSTALW_LECS =>       "/share/apps/bin/clustalw"; # on the lecs
use constant CLUSTALW =>            "/usr/local/bin/clustalw"; # on the biocluster
use constant CLUSTALW_1_82 =>       "/bioseq/Programs/ClustalW_1.82/clustalw1.82/clustalw"; # for servers who came from biocluster (Selecton?, old ConSurf, ConSeq)
use constant CLUSTALW_1_81 =>       "/bioseq/Programs/ClustalW_1.81/clustalw1.81/clustalw"; # for servers who came from biocluster (Selecton?, old ConSurf, ConSeq)
use constant CLUSTALW_2_0_10 =>       "/bioseq/Programs/ClustalW_2.0.10/clustalw-2.0.10-linux-i386-libcppstatic/clustalw2"; # for servers who came from biocluster (Selecton?, old ConSurf, ConSeq)

use constant MAFFT_LINSI =>	    "/usr/local/bin/mafft-linsi"; # on the biocluster
use constant MAFFT =>         	    "/usr/local/bin/mafft"; # on the biocluster
#use constant MAFFT_GUIDANCE =>      "/groups/pupko/privmane/bin/mafft"; #v6.711b
#use constant MAFFT_LINSI_GUIDANCE => 	    "/groups/pupko/privmane/bin/mafft --localpair --maxiterate 1000"; #v6.711b
#use constant MAFFT_GUIDANCE =>      "/bioseq/Programs/MAFFT_6.711b/mafft"; #v6.711b
use constant MAFFT_GUIDANCE =>       "/bioseq/Programs/MAFFT_6.833/bin/mafft"; #v6.833
#use constant MAFFT_GUIDANCE =>       "/bioseq/Programs/MAFFT_6.857/bin/mafft"; #v6.857 !!! make sure: 'setenv MAFFT_BINARIES /bioseq/Programs/MAFFT_6.857/mafft-6.857-with-extensions/binaries' BEFORE
#use constant MAFFT_LINSI_GUIDANCE => 	    "/bioseq/Programs/MAFFT_6.711b/mafft --localpair --maxiterate 1000"; #v6.711b
use constant MAFFT_LINSI_GUIDANCE => "/bioseq/Programs/MAFFT_6.833/bin/mafft --localpair --maxiterate 1000"; #v6.833
#use constant MAFFT_LINSI_GUIDANCE => "/bioseq/Programs/MAFFT_6.857/bin/mafft --localpair --maxiterate 1000"; #v6.857 !!! make sure: 'setenv MAFFT_BINARIES /bioseq/Programs/MAFFT_6.857/mafft-6.857-with-extensions/binaries' BEFORE
use constant PRANK_LECS =>          "/share/apps/bin/prank"; # on the lecs
use constant PRANK =>		    "/usr/local/bin/prank"; # on the biocluster
use constant T_COFFEE =>            "/share/apps/T-COFFEE-8.47/bin/binaries/linux/t_coffee"; # requiers setenv PATH /share/apps/T-COFFEE-8.47/bin/binaries/linux:$PATH  
use constant PAGAN_LECS =>          "/share/apps/pagan-msa/bin/pagan"; # requires:  "module load gcc/gcc461" before!!

use constant TREE_VIEWER_DIR =>     "/bioseq/ConSurf_old/treeViewer/";
use constant PACC_path =>           "/bioseq/ConSeq/external_scripts/PACC/";
use constant RATE4SITE_BIOC_VER =>  "/bioseq/rate4site/BioCluster_Nov_06_dev/rate4site.exe";
use constant RATE4SITE_SLOW_BIOC_VER => "/bioseq/rate4site/BioCluster_Nov_06_dev/rate4siteSlow.exe";
use constant RATE4SITE =>           "/db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4site.exe";
use constant RATE4SITE_SLOW =>      "/db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4siteSlow.exe";
use constant RATE4SITE_SLOW_LECS => "/share/apps/bin/rate4site_slow";
use constant RATE4SITE_LOCAL =>     "/bioseq/rate4site/rate4site";
use constant RATE4SITE_SLOW_LOCAL =>"/bioseq/rate4site/rate4site.doubleRep";
use constant RATE4SITE_WITH_LG =>   "/bioseq/rate4site/With_LG/rate4site";
use constant RATE4SITE_WITH_LG_SLOW => "/bioseq/rate4site/With_LG/rate4site.doubleRep";
use constant RUBY => "/share/apps/bin/ruby"; #"/usr/bin/ruby";
#use constant CD_HIT_DIR =>          "/db1/Local/src/cd-hit_redundency/";
use constant CD_HIT_DIR =>          "/bioseq/cd_hit/";
use constant PREDICT_PACC =>	    "/bioseq/ConSeq/external_scripts/PACC/run.sh";
use constant MSA_to_HSSP =>	    "/bioseq/ConSeq/external_scripts/PACC/MSA2hssp.pl";
#use constant SEMPHY =>              "/groups/pupko/privmane/alignment/run/semphy"; #on Biocluster
use constant SEMPHY =>              "/bioseq/Programs/Semphy/semphy.doubleRep";

#internal programs
use constant EPITOPIA_EXECUTABLES => "/bioseq/epitopia/executables";

# constant values
use constant BLAST_MAX_HOMOLOGUES_TO_DISPLAY => 500;
use constant BLAST_PDB_MAX_HOMOLOGUES_TO_DISPLAY => 25;
use constant CONSURF_PIPE_FORM =>               "/bioseq/ConSurf_old/consurf_pipe.form";
use constant SELECTON_MAX_NUCLEOTIDE =>         15000;
use constant MAX_WALLTIME => "96:00:00";

# Queue Details
use constant BIOSEQ_NODE =>             "bioc01.tau.ac.il"; #Node on BioCluster dedicated to Bioseq runs (Not part of the queue)
#use constant MAX_QUEUE_RUNS =>          60;
use constant MAX_QUEUE_RUNS =>          999;

# external links
use constant RCSB_WEB =>                "http://www.rcsb.org/";
use constant PYMOL_WEB =>               "http://pymol.sourceforge.net/";
use constant CHIMERA_WEB =>             'http://www.rbvi.ucsf.edu/chimera/';
use constant CHIMERA_SAVING_FIGURE =>   'http://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/print.html';
use constant CHIMERA_DOWNLOAD =>        CHIMERA_WEB."download.html";
use constant MSA_CONVERT =>             'http://www.ebi.ac.uk/cgi-bin/readseq.cgi';
use constant MSA_FORMATS =>             'http://www.ebi.ac.uk/help/formats.html';

# redirect pages
use constant CONSURF_REDIRECT_PAGE => CONSURF_URL."too_many_runs.html";
use constant SELECTON_REDIRECT_PAGE => SELECTON_URL."/too_many_runs.html";
use constant CONSEQ_REDIRECT_PAGE => CONSEQ_URL."too_many_runs.html";
use constant PEPITOPE_REDIRECT_PAGE => PEPITOPE_URL."too_many_runs.html";

#faq pages
use constant CONSURF_TREE_FAQ => CONSURF_URL.'quick_help.html#note5';

#Files Name Conventions
use constant TEMPLATES_LIST_FILE=>"List_of_Templates";
use constant PISA_ERRORS_FILE=>"PISA_Errors";


#---------------------------------------------
sub print_to_output{
    my $OutHtmlFile = shift;
    my $server_name = shift;
    my $run_name = shift;
    my $recipient = shift;    
    
    open OUTPUT, ">>$OutHtmlFile";
    flock OUTPUT, 2;
    print OUTPUT "\n<p><font size=+3 color='red'>ERROR! $server_name session has been terminated: </font>\n<br><b>A system error occured during the calculation. Please try to run $server_name again in a few minutes.</b>\n</p>\n";
    print OUTPUT "<H3><center>For assistance please <a href=\"mailto:".ADMIN_EMAIL."?subject=".$server_name."%20Run%20No:%20".$run_name."\">contact us</a> and mention this number: $run_name</H3>\n";
    flock OUTPUT, 8;
    close OUTPUT;    
    &send_mail($server_name, $recipient, $run_name, "error","error") if ($recipient ne "NO");
    &stop_reload($OutHtmlFile);
}
#---------------------------------------------

# in case the desired mail report on error: the vars $email_subject and $email_message should be 'error'
sub send_mail { # to user
    my $server_name = shift;
    my $recipient = shift;
    my $run_name = shift;
    my $email_subject= shift;
    my $email_message = shift;
    my $email_attach = shift;
    my $from_server = "";
    $from_server = shift;
    my $OutputURL;
    my $mail;
    
    if ($server_name eq "Selecton") {$OutputURL = SELECTON_URL."/results/$run_name"."/output.html";}
    elsif ($server_name eq "ConSeq") {$OutputURL = CONSEQ_URL."results/$run_name"."/output.html";}
    elsif ($server_name eq "Epitopia") {$OutputURL = EPITOPIA_URL."results/$run_name"."/output.html";}
    elsif ($server_name eq "pepitope") {$OutputURL = PEPITOPE_URL."results/$run_name"."/output.html";}
    elsif ($server_name eq "ConSurf") {$OutputURL = CONSURF_URL."results/$run_name"."/output.html";}
    elsif ($server_name eq "QuasiMotiFinder") {$OutputURL = QMF_URL."results/$run_name"."/output.html";}
    elsif ($server_name eq "fastml") {$OutputURL = FASTML_URL."results/$run_name"."/output.html";}
        
    $email_subject = "Error in $server_name running" if $email_subject eq "error";
    $email_message = "Hello!\n\nUnfortunately there was an error while running the $server_name server.\nPlease click on the following link to see more details\nWe apologize for the inconvenience\n\n$OutputURL\n" if $email_message eq "error"; 
    chdir SEND_EMAIL_DIR;
    chdir SEND_EMAIL_DIR_IBIS if ($from_server eq "ibis");
    $mail ='perl sendEmail.pl -f \''.ADMIN_EMAIL.'\' -t \''.$recipient.'\' -u \''.$email_subject.'\' -s '.SMTP_SERVER.' -m \''.$email_message."\'";
    #$mail ='perl sendEmail.pl -f \''.ADMIN_EMAIL.'\' -t \''.$recipient.'\' -u \''.$email_subject.'\' -xu '.ADMIN_USER_NAME.' -xp '.ADMIN_PASSWORD.' -s '.SMTP_SERVER.' -m \''.$email_message."\'";
    if ($email_attach ne '') {$mail.=" -a $email_attach";}
    `$mail`;
}
#---------------------------------------------
sub stop_reload{
    my $OutHtmlFile = shift;
    
    sleep 10;
    open OUTPUT, "<$OutHtmlFile";
    my @output = <OUTPUT>;
    close OUTPUT;   
    open OUTPUT, ">$OutHtmlFile";
    foreach my $line (@output){  # we remove the refresh lines and the button which codes for Selecton cancelled job
        unless ($line =~ /REFRESH/i or $line =~ /NO-CACHE/i or $line =~ /ACTION=\"\/cgi\/kill_process.cgi/ or
            $line =~ /VALUE=\"Cancel Selecton Job\"/ or $line =~ /TYPE=hidden NAME=\"pid\"/ or 
            $line =~ /TYPE=hidden NAME=\"selecton_http\"/ or $line =~ /TYPE=hidden NAME=\"run_no\"/ or
            $line =~ /<!--job_/){
                print OUTPUT $line;
        }
    }
    close OUTPUT;
}
#---------------------------------------------
sub print_Q_status_in_html{
    my $html_file = shift;
    my $_status = shift;
    my $_time = shift;
    my $_estimated_run_time = shift;
    
    my ($line, $line1, $line2);
    my $out = "/bioseq/ELANA/from_GENERAL_CONST.txt";
        
    $_time = "" if ($_time eq "no");
    unless (open HTML, "+>>".$html_file) {
        return "print_Q_status_in_html : Could not open file $html_file to update the status. Status is: $_status  reason: $!\n";}
    else{
        flock HTML, 2;
        seek HTML, 0, 0; #rewind the pointer to the beginning
        my @html_lines = <HTML>; # read the contents into the array
        truncate HTML, 0; # remove all the information, The 0 represents the size of the file that we want
        foreach (@html_lines){            
            if(/<!--job_stat--><.+>Your job status is:<\/a> (.+)<br>/){
                if ($_status ne ""){
                    s/$1/$_status/;
                }
            }
            elsif(/<!--job_pass-->The time that passed since submitting the query is: (.+)<br>/){
                if($_time ne ""){
                    s/$1/$_time/;
                }
            }
            elsif(/<!--(job_time--)Estimated run time is: (-->)/ and $_estimated_run_time ne "none"){
                $line = $_;
                $line1 = $1;
                $line2 = $2;
                if ($_estimated_run_time =~ m/\d+:\d+:\d+:\d+/) {
                    $_estimated_run_time .= " days";
                }
                elsif ($_estimated_run_time =~ m/\d+:\d+:\d+/) {
                    $_estimated_run_time .= " hours";
                }
                elsif($_estimated_run_time =~ m/\d+:\d+/){
                    $_estimated_run_time .= " minutes";
                }
                $_ = $line; # since we make another RE comparison, the original values of $_ and $1 are changing, therefore we must save them at the beginning and change them back here.
                s/$line2/$_estimated_run_time<br>/; # the reason we first substitue the second part, is that the first part creates an expression --> which might be wrongly replaced with this value
                s/$line1/$line1>/;                
            }
        }
        print HTML $_ foreach (@html_lines);
        flock HTML, 8;
        close HTML;
        return "OK";
    }
}


# in case the desired mail report on error: the vars $email_subject and $email_message should be 'error'
sub send_mail2 { # to user
    my $server_name = shift;
    my $recipient = shift;
    my $run_name = shift;
    my $email_subject= shift;
    my $email_message = shift;
    my $email_attach = shift;
    my $from_server = shift;
    my $OutputURL;
    my $mail;
    
    if ($server_name eq "Selecton") {$OutputURL = SELECTON_URL."/results/$run_name"."/output.html";}
    elsif ($server_name eq "ConSeq") {$OutputURL = CONSEQ_URL."results/$run_name"."/output.html";}
    elsif ($server_name eq "Epitopia") {$OutputURL = EPITOPIA_URL."results/$run_name"."/output.html";}
    elsif ($server_name eq "pepitope") {$OutputURL = PEPITOPE_URL."results/$run_name"."/output.html";}
    elsif ($server_name eq "ConSurf") {$OutputURL = CONSURF_URL."results/$run_name"."/output.html";}
    elsif ($server_name eq "QuasiMotiFinder") {$OutputURL = QMF_URL."results/$run_name"."/output.html";}
    elsif ($server_name eq "fastml") {$OutputURL = FASTML_URL."results/$run_name"."/output.html";}
        
    $email_subject = "Error in $server_name running" if $email_subject eq "error";
    $email_message = "Hello!\n\nUnfortunately there was an error while running the $server_name server.\nPlease click on the following link to see more details\nWe apologize for the inconvenience\n\n$OutputURL\n" if $email_message eq "error"; 
    chdir SEND_EMAIL_DIR;
    chdir SEND_EMAIL_DIR_IBIS if ($from_server eq "ibis");
    $mail ='perl sendEmail.pl -f \''.ADMIN_EMAIL.'\' -t \''.$recipient.'\' -u \''.$email_subject.'\' -s '.SMTP_SERVER.' -m \''.$email_message."\'";
    #$mail ='perl sendEmail.pl -f \''.ADMIN_EMAIL.'\' -t \''.$recipient.'\' -u \''.$email_subject.'\' -xu '.ADMIN_USER_NAME.' -xp '.ADMIN_PASSWORD.' -s '.SMTP_SERVER.' -m \''.$email_message."\'";
    if ($email_attach ne '') {$mail.=" -a $email_attach";}
    $mail = 'sh -c \' $mail 2>/dev/null\'';
    `$mail`;
}
1;
