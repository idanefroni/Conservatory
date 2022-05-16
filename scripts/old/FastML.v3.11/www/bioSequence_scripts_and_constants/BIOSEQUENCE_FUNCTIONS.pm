#!/usr/bin/perl

package BIOSEQUENCE_FUNCTIONS; #don't forget: a package must end with a return value (1; in the end)!!!!!

use strict;
use GENERAL_CONSTANTS;

#------------------------------------------------------------------------------------
sub subtract_time_from_now{
# receieves the begin time in format of: HH:MN:SS DD-MO-YEAR
# returns the the time (in hours) passed from the time of calculation to the begin time.
# if an error was found during calculation: returns "no"
# error will be found in case the time that passed is more than 1 month different.
    
    my $begin_time = shift;
    $begin_time .= " ".shift;
    my %date1;
    my %date2;
    my $date1_ref;
    my $date2_ref;
    my @time_difference;
    my $dir_counter = 0;
    
    $begin_time =~ m/(\d+):(\d+):(\d+) (\d+)-(\d+)-(\d+)/;
    %date1 = (Year => $6, Month => $5, Day => $4, Hour => $1, Minute => $2, Second => $3);
    %date2 = (Year => "", Month => "", Day => "", Hour => "", Minute => "", Second => "");
    &convert_currentTime(\%date2);
    
    @time_difference = &compare_time(\%date1, \%date2);
    #if ($time_difference[0] eq "no") {
    #    return "no";
    #}
    if ($time_difference[0] =~ m/error/) {
        return $time_difference[0];
    }
    else{
        return $time_difference[1];
    }
}
#------------------------------------------------------------------------------------
# the routine converts the "Begin/End" time line from Selecton's log files to a numeric string.
# it insertes the new values to the hash' reference .
sub convertTime
{
    my $inputTimeString = $_[0];
    my $answer = $_[1]; #reference to hash
    my %months =
    (   Jan => "01", Feb => "02", Mar => "03", Apr => "04", May => "05", Jun => "06",
        Jul => "07",Aug => "08", Sep => "09", Oct => "10", Nov => "11", Dec => "12");

    if ($inputTimeString =~ m/(\d+):(\d+):(\d+),\s+\w+\s(\w+)\s(\d+),\s(\d+)/)
    {
        my $HH = &convertNum($1);
        my $MN = &convertNum($2);
        my $SS = &convertNum($3);
        my $MM = $months{$4};
        my $DD = &convertNum($5);
        my $YYYY = $6;
        
        $answer->{Year} = $YYYY;
        $answer->{Month} = $MM;
        $answer->{Day} = $DD;
        $answer->{Hour} = $HH;
        $answer->{Minute} = $MN;
        $answer->{Second} = $SS;
    }
}#convertTime
#__________________________________________________________
# converts a number from one digit to 2 digits
sub convertNum 
{
    my $input_num = shift;
    if ($input_num < 10)
        {return "0".$input_num;}
    else
        {return $input_num;}
}

#__________________________________________________________
# calculates the time differences by comparing seperately months, days, minutes and seconds.
# this functions assumes that the year is the same year.
# input: references to 2 hashs with time's details
# output: string with time difference, messured by hours:minutes:seconds

sub compare_time()
{
    my $time1 = $_[0]; #refernce to the time array
    my $time2 = $_[1]; #refernce to the time array
    my $time_difference;
    my $no_of_Days_passed;
    my $no_of_hours_passed;
    my %days_each_month = ('01' => '31', '02' => '28', '03' => '31', '04' => '30', '05' => '31', '06' => '30',
                        '07' => '31', '08' => '31', '09' => '30', '10' => '31', '11' => '30', '12' => '31');
    
    if ($time1->{Month} eq $time2->{Month}) {#same month
        if ($time1->{Day} eq $time2->{Day}) {#same day
            if ($time2->{Hour} >= $time1->{Hour}) {#compare hour: h2>h1
                $time_difference = &calculate_time_difference($time1->{Hour}, $time2->{Hour}, $time1->{Minute}, $time2->{Minute}, $time1->{Second}, $time2->{Second}, 0);
            }
            else{
                #return("no");
                 return("error: H1 is: $time1->{Hour} H2 is: $time2->{Hour} it is the same day, therefor it is impossible that H1>H2. \n");
            }
        }
        else  {# different day
            if ($time2->{Day} >= $time1->{Day}){
                $no_of_Days_passed = ($time2->{Day}-$time1->{Day});                
                $time_difference = &calculate_time_difference($time1->{Hour}, $time2->{Hour}, $time1->{Minute}, $time2->{Minute}, $time1->{Second}, $time2->{Second}, $no_of_Days_passed);
            }
            else{
                #return("no");
                 return("error: D1 is: $time1->{Day} D2 is: $time2->{Day}, it is impossible in the same month that D1>D2.\n");
            }
        }
    }
    else {#different month
        #if ($time2->{Month} >= $time1->{Month}){
            if (($time2->{Month} - $time1->{Month})>1 or ($time2->{Month} - $time1->{Month})<0){
                #return("no");
                 return("error: M1 is: $time1->{Month}, M2 is: $time2->{Month}. The program doesn't allow a difference bigger than 1 month.\n");
            }
            else {# 1 month difference
                $no_of_Days_passed = ($time2->{Day} + $days_each_month{$time1->{Month}} - $time1->{Day});               $time_difference = &calculate_time_difference($time1->{Hour}, $time2->{Hour}, $time1->{Minute}, $time2->{Minute}, $time1->{Second}, $time2->{Second}, $no_of_Days_passed);
            }            
        #}
        #else{
            #return("no");#, "error: M1 is: $time1->{Month}, M2 is: $time2->{Month}. It is impossible for M1 to be bigger within the same year\n");
        #}
    }
    return ("yes", $time_difference);    
} # finish: compare_time()

#__________________________________________________________
# does the part of calculating minutes and seconds difference.
# input: hours difference (just for formating the string output) M1, M2, D1, D2
# output: string output, sent to the compare_time() function for display
sub calculate_time_difference()
{
    my $hour1 = $_[0];
    my $hour2= $_[1];    
    my $minute1 = $_[2];
    my $minute2 = $_[3];
    my $second1 = $_[4];
    my $second2 = $_[5];
    my $days_passed = $_[6];
    my $minutes_passed;
    my $seconds_passed;
    my $hours_passed;
    my $reduce_minute = "no";
    my $reduce_hour = "no";
    my $reduce_day = "no";
    
    # seconds
    if ($second2>=$second1)
        {$seconds_passed = $second2-$second1;}
    else
        {$seconds_passed = 60+$second2-$second1;
         $reduce_minute = "yes";}
    #minutes    
    if ($minute2>=$minute1)
        {$minutes_passed = $minute2-$minute1;}
    else
        {$minutes_passed = 60+$minute2-$minute1;
         $reduce_hour = "yes";}
    if ($reduce_minute eq "yes")
    {
        if ($minutes_passed == 0)
           {$minutes_passed = 59;}
        else
           {$minutes_passed -=1;}
    }
    #hours
    if ($hour2>=$hour1)
        {$hours_passed = $hour2-$hour1;}
    else
        {$hours_passed = 24+$hour2-$hour1;
         $reduce_day = "yes";}
    if ($reduce_hour eq "yes")
    {
        if($hours_passed == 0)
            {$hours_passed = 23;}
        else
            {$hours_passed -=1;}
    }
    #days
    if ($days_passed > 0)
    {
        if($reduce_day eq "yes")
            {$days_passed-=1;}
        $hours_passed += 24*$days_passed;        
    }
    $hours_passed = &convertNum($hours_passed);
    $minutes_passed = &convertNum($minutes_passed);
    $seconds_passed = &convertNum($seconds_passed);
    return  "$hours_passed:$minutes_passed:$seconds_passed";
}
#------------------------------------------------------------------------------------
sub convert_currentTime {
    my $answer = shift; #reference to hash
   my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
   my $year = 1900 + $yearOffset;
    $second = &convertNum($second);
    $minute = &convertNum($minute);
    $hour = &convertNum($hour);
    $month = &convertNum($month+1);
    $dayOfMonth = &convertNum($dayOfMonth);

    $answer->{Year} = $year;
    $answer->{Month} = $month;
    $answer->{Day} = $dayOfMonth;
    $answer->{Hour} = $hour;
    $answer->{Minute} = $minute;
    $answer->{Second} = $second;    

#print "Current time is: ".$answer->{Hour}.":".$answer->{Minute}.":".$answer->{Second}."  ".$answer->{Day}."-".$answer->{Month}."-".$answer->{Year}."\n";    
        
}
#---------------------------------------------
sub check_if_user_is_allowed{
    
    my $server_name = shift;
    my $user_ip = shift;
    my $user_email = shift;
    
    my $file_to_open;
    
    my %ip_total = ();
    my ($ip, $_mail, $redirect_html);
    
    if ($server_name eq "consurf"){
        $redirect_html = GENERAL_CONSTANTS::CONSURF_REDIRECT_PAGE;
        $file_to_open = GENERAL_CONSTANTS::CONSURF_RUNNING_JOBS;        
    }
    elsif ($server_name eq "selecton"){
        $redirect_html = GENERAL_CONSTANTS::SELECTON_REDIRECT_PAGE;
        $file_to_open = GENERAL_CONSTANTS::SELECTON_RUNNING_JOBS;        
    }
    elsif ($server_name eq "conseq"){
        $redirect_html = GENERAL_CONSTANTS::CONSEQ_REDIRECT_PAGE;
        $file_to_open = GENERAL_CONSTANTS::CONSEQ_RUNNING_JOBS;        
    }
    elsif ($server_name eq "pepitope"){
        $redirect_html = GENERAL_CONSTANTS::PEPITOPE_REDIRECT_PAGE;
        $file_to_open = GENERAL_CONSTANTS::PEPITOPE_RUNNING_JOBS;        
    }
    if (-e $file_to_open and !(-z $file_to_open)){
        open RUN_LIST, $file_to_open;
        flock RUN_LIST, 2;
        while (<RUN_LIST>){
            chomp;
            if(/^(null_)?\d+ (.+) (.+)$/){
                $ip = $2;
                $_mail = $3;
                if (exists $ip_total{$ip}){
                    $ip_total{$ip}++;}
                else{
                    $ip_total{$ip} = 1;}
                if (exists $ip_total{$_mail}){
                    $ip_total{$_mail}++;}
                else{
                    $ip_total{$_mail} = 1;}
            }
            #redirects unwanted visitors to the site
            if ($ip =~ /66\.232\.100\.62/ or $ip =~ /83\.97\.\177\.107/ or $ip =~ /91\.74\.160\.18/){
                #print "Location: http://www.tau.ac.il/lifesci/\n\n";
                exit;
            }
        }
        close RUN_LIST;
        if ((exists $ip_total{$user_ip} && $ip_total{$user_ip} >=7) or (exists $ip_total{$user_email} && $ip_total{$user_email} >= 7)){
        # output a message to the user that he cannot continue the run
            print "Location: $redirect_html\n\n";
            exit;        
        }
    }    
}
#---------------------------------------------
# the values for this statistics were determined in a statistical test we did on November 2007,
# on Selecton seccsful runs for 3 months on the bioinfo machine
#sub selecton_estimated_run_time1{
#    my $seq_times_length = shift;
#    my $model = shift;
#
#    my ($time_in_minutes, $time_in_hours, $time_in_days);
#    # set the time according to each model's parameters
#    $time_in_minutes = $seq_times_length*0.0251 +  20.345 if ($model eq "M8");
#    $time_in_minutes = $seq_times_length*0.0256 + 17.391 if ($model eq "MEC");
#    # to be on the safe side - we add 20% for the time
#    $time_in_minutes = int($time_in_minutes*1.2);
#    # calculate time in DD:HH:MM:SS format
#    $time_in_minutes = int($time_in_minutes); # remove numbers after the "."
#    
#    return(&time_in_days_from_minutes($time_in_minutes));
#}
#---------------------------------------------
# the values for this statistics were determined in a statistical test we did on October 2009, on Selecton seccsful runs for a few month on biocluster.
# the file can be found at: /bioseq/Selecton/total_models_statistics.csv
sub selecton_estimated_run_time{
    my $seq_length = shift;
    my $num_of_seq = shift;
    my $model = shift;
    
    my ($time_in_minutes, $time_in_hours, $time_in_days);
    # set the time according to each model's parameters
    if ($model eq "MEC"){
        $time_in_minutes = $seq_length*$num_of_seq*0.0035 + 12.677 ;
    }
    elsif ($model eq "M8"){
        if($num_of_seq<11){
            $time_in_minutes = $seq_length*$num_of_seq*0.022 + 3.5198;
        }
        elsif($num_of_seq>10 and $num_of_seq<21){
            $time_in_minutes = $seq_length*$num_of_seq*0.0025 + 14.82;
        }        
        elsif($num_of_seq>20 and $num_of_seq<31){
            $time_in_minutes = $seq_length*$num_of_seq*0.0021 + 35.153;
        }
        elsif($num_of_seq>30 and $num_of_seq<41){
            $time_in_minutes = $seq_length*$num_of_seq*0.0026 + 48.412;
        }
        elsif($num_of_seq>40 and $num_of_seq<51){
            $time_in_minutes = $seq_length*$num_of_seq*0.0024 + 65.947;
        }
        else{
            $time_in_minutes = $seq_length*$num_of_seq*0.003 + 91.341;
        }
    }
    
    # to be on the safe side - we triple the time
    $time_in_minutes = int($time_in_minutes*3);
    # calculate time in DD:HH:MM:SS format
    $time_in_minutes = int($time_in_minutes); # remove numbers after the "."
    
    return(&time_in_days_from_minutes($time_in_minutes));
}
#---------------------------------------------
# input: int represents sum of minutes
# output: time in format: HH:MM:SS (maybe change in the future to time in format: DD:HH:MM:SS)
sub time_in_days_from_minutes{
    my $minutes = shift;
    my $hours = 0;
    my $days = 0;
    my $ret = "";
    
    if($minutes <=59){
        $ret = $minutes.":00"; 
    }
    elsif ($minutes >59){
        $hours = int($minutes/60);
        $minutes = $minutes%60;
        $minutes = new_num($minutes);
        # ---- if the format needed inculdes only hours
        $hours = new_num($hours);
        $ret = $hours.":".$minutes.":00"; 
        ## --- if the format needed inculdes days in seperate
        #if($hours <= 23){
        #    $hours = new_num($hours);
        #    $ret = $hours.":".$minutes.":00"; 
        #}
        #else{
        #    $days = int($hours/24);
        #    $hours = $hours%24;
        #    $hours = new_num($hours);
        #    $days = new_num($days);
        #    $ret = $days.":".$hours.":".$minutes.":00";
        #}
    }    
    return $ret;
}
#---------------------------------------------
# gives the number in minimum 2 digits
sub new_num{
    my $num = shift;
    ($num < 10) ? return "0".$num : return $num;
}
#---------------------------------------------
# returns the time in format hh:mm:ss dd:mn:yyy
sub printTime {
   my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
   my $year = 1900 + $yearOffset;
   
   $second = &new_num($second);
   $minute = &new_num($minute);
   $hour = &new_num($hour);
   $month = &new_num($month+1);
   $dayOfMonth = &new_num($dayOfMonth);
   
   return "$hour:$minute:$second $dayOfMonth-".$month."-$year";
}
#---------------------------------------------
sub printYear {
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	return $year;
}
#---------------------------------------------
sub printMonth {
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    # localtime returns array. in its 5th cell (4 when coutin from 0) is the number denotin the current month minus 1
    # example: in December, $time[4] = 11. So in the above @months array, $months[11] = Dec.
    my @time = localtime();
    return $months[$time[4]];
}
#---------------------------------------------
# input: the server name and run_name
# the routine will remove this run_name from the list of running jobs
# please note: the var $server should be spelled: "Selecton", "ConSurf"
sub remove_job_from_running_log{
    
    my $server = shift;
    my $run_name = shift;
    my $log;
    
    if($server eq "Selecton") {
        $log = GENERAL_CONSTANTS::SELECTON_RUNNING_JOBS;}
    elsif($server eq "ConSurf"){
        $log = GENERAL_CONSTANTS::CONSURF_RUNNING_JOBS;}
    elsif($server eq "ConSeq"){
        $log = GENERAL_CONSTANTS::CONSEQ_RUNNING_JOBS;}
    elsif($server eq "pepitope"){
        $log = GENERAL_CONSTANTS::PEPITOPE_RUNNING_JOBS;}
    
    # remove the job from the running jobs list
    open LIST, "+>>".$log;
    flock LIST, 2;
    seek LIST, 0, 0; #rewind the pointer to the beginning
    my @all_lines_in_list = <LIST>; # read the contents into the array
    truncate LIST, 0; # remove all the information, The 0 represents the size of the file that we want
    foreach (@all_lines_in_list){
        chomp;
        unless(/$run_name/){
            print LIST $_."\n";
        }
    }
    flock LIST, 8;
    close LIST;
}
#---------------------------------------------
# prints the job in the queuing jobs list
sub enqueue_job{
    my $job_num = shift;
    my $server = shift;
    my $run_name = shift;
    my $ret = "ok";
    
    unless (open LIST, ">>".GENERAL_CONSTANTS::QUEUING_JOBS){
       $ret = "Could not open file ".GENERAL_CONSTANTS::QUEUING_JOBS.". Reason: $!\nThe job was not listed in the queuing_jobs list.\n".printTime();
    }
    else{
       flock LIST, 2; # locks the list, so no other process will write to it. On the same time - if the list is currently locked by another process - it waits until the list file is realeased. The "2" and "8" are the operation symbols for "lock" and "unlock".
       print LIST "$job_num $server $run_name ".printTime()."\n";
       flock LIST, 8;
       close LIST;
    }
    return $ret;
}

#------------------------------------------------------
# prints the job in the bioseq node running jobs list
sub enqueue_job_to_bioseq_node{
    my $job_num = shift;
    my $server = shift;
    my $run_name = shift;
    my $ret = "ok";
    
    unless (open LIST, ">>".GENERAL_CONSTANTS::JOBS_ON_BIOSEQ_NODE){
       $ret = "Could not open file ".GENERAL_CONSTANTS::JOBS_ON_BIOSEQ_NODE.". Reason: $!\nThe job was not listed in the bioseq node running job list.\n".printTime();
    }
    else{
       flock LIST, 2; # locks the list, so no other process will write to it. On the same time - if the list is currently locked by another process - it waits until the list file is realeased. The "2" and "8" are the operation symbols for "lock" and "unlock".
       print LIST "$job_num $server $run_name ".printTime()."\n";
       flock LIST, 8;
       close LIST;
    }
    return $ret;
}
#------------------------------------------------------
# prints the job in the bioseq node waiting jobs list
sub waiting_jobs_for_bioseq_node{
     my $server = shift;
    my $run_name = shift;
    my $ret = "ok";
    
    unless (open LIST, ">>".GENERAL_CONSTANTS::JOBS_WAITING_BIOSEQ_NODE){
       $ret = "Could not open file ".GENERAL_CONSTANTS::JOBS_WAITING_BIOSEQ_NODE.". Reason: $!\nThe job was not listed in the bioseq node waiting job list.\n".printTime();
    }
    else{
       flock LIST, 2; # locks the list, so no other process will write to it. On the same time - if the list is currently locked by another process - it waits until the list file is realeased. The "2" and "8" are the operation symbols for "lock" and "unlock".
       print LIST "$server $run_name ".printTime()."\n";
       flock LIST, 8;
       close LIST;
    }
    return $ret;
}
#------------------------------------------------------
# remove the job from the bioseq node waiting jobs list
sub remove_job_from_bioseq_node_waiting_list{
    my $server = shift;
    my $run_name = shift;
    my $ret = "ok";
    
    unless (open LIST, "+>>".GENERAL_CONSTANTS::JOBS_WAITING_BIOSEQ_NODE){
       $ret = "Could not open file ".GENERAL_CONSTANTS::JOBS_WAITING_BIOSEQ_NODE.". Reason: $!\nThe job was not listed in the bioseq node waiting job list.\n".printTime();
    }
    else{
       flock LIST, 2;
       seek LIST, 0, 0; #rewind the pointer to the beginning
       my @all_lines_in_list = <LIST>; # read the contents into the array
       truncate LIST, 0; # remove all the information, The 0 represents the size of the file that we want
       foreach my $line (@all_lines_in_list){
       		chomp;
		if (($line=~/$run_name/) and ($line=~/$server/))
		{
			$line = ""; # removing this line from the lines array
		}
		elsif ($line =~/([A-Za-z0-9])+/)
		{
			print LIST "$line\n";
		}
	}
    flock LIST, 8;
    close LIST;
    }
   return $ret; 
}
#---------------------------------------------
# input:  path to pdb file
# output: 3 options:
# 1. --PDB_NOT_OPEN if couldn't open the pdb file
# 2. --NO_CHAINS if no chain was founded in column 22
# 3. string with all the chains founded in this pdb.

sub which_chain_in_pdb_and_seqres{
    my $input_pdb = shift;
    my $chain_founded;
    my %all_chains;
    my @ret;
    my $seqres_found = "--SEQRES_no";
    
    unless (open PDB, $input_pdb){
        @ret = ("--PDB_NOT_OPEN $input_pdb $!");
        return \@ret;}
    while (<PDB>){
        if (/^ATOM/){
            $chain_founded = substr $_, 21, 1;
            if (!(exists $all_chains{$chain_founded})){
                $all_chains{$chain_founded} = 1;
            }
        }
        if ($seqres_found eq "--SEQRES_no" && /^SEQRES/){
            $seqres_found = "--SEQRES_yes";
        }
    }
    close PDB;    
    $chain_founded = "";
    foreach my $key (keys %all_chains){
        $chain_founded.=$key;
    }
    if($chain_founded !~ /\S/){
        @ret = ("--NO_CHAINS", $seqres_found);}
    else{
        @ret = ($chain_founded, $seqres_found);}
    return \@ret;
}
#---------------------------------------------
# input : 1. path to a pdb file, where there is no chain identifier in the 22 column of ATOM and 12 column of SEQRES
#         2. one letter denotes a chain identifier to add
# output : the same file, in the same path, where the letter given as input is added to the previously empty 22 column.
sub add_chain_to_pdb{
    my $input_pdb = shift;
    my $chain_id_to_add = shift;
    
    my ($beg_line, $end_line, $line);
    
    open PDB_IN, "+>>".$input_pdb;
    seek PDB_IN, 0, 0;
    my @all_lines_in_pdb = <PDB_IN>;
    truncate PDB_IN, 0;
    foreach(@all_lines_in_pdb){
        if (/^ATOM/){
            $line = $_;
            $beg_line = substr $line, 0, 21;
            $end_line = substr $line, 22, length($line);
            $_ = $beg_line.$chain_id_to_add.$end_line;
        }
        elsif (/^SEQRES/){
            $line = $_;
            $beg_line = substr $line, 0, 11;
            $end_line = substr $line, 12, length($line);
            $_ = $beg_line.$chain_id_to_add.$end_line;
        }
        print PDB_IN $_;
    }    
    close PDB_IN;    
}
#---------------------------------------------
sub convertNewline{
    # runs dos2unix, the program that converts plain text files in DOS/MAC format to UNIX format.
    my $inputFilePath = shift;
    my $WorkingDir = shift;
    my $dos2unix="cd $WorkingDir;dos2unix -q $inputFilePath";        
    system "$dos2unix";
    # if the input file was in mac format, the simple dos2unix will not work.
    # read the file - if it is only one line, it might mean that the new line characters
    # are not read well (for example: ^M). Trying to run dos2unix again, saying the format is mac
    $WorkingDir.='/' unless $WorkingDir =~ /\/$/;
    if (open FILE, $WorkingDir.$inputFilePath){
        my $num_of_lines = 0;
        while (<FILE>){
            $num_of_lines++;
        }
        close FILE;
        if ($num_of_lines==1){
            $dos2unix="cd $WorkingDir;dos2unix -c mac $inputFilePath -q ";
            system "$dos2unix";
        }
    }
    
}
#---------------------------------------------
sub removeEndLineExtraChars{
    # remove extra chars on end of lines (^M,spaces);
    my $inputFilePath = shift;
    my $WorkingDir = shift;
    $WorkingDir.='/' unless $WorkingDir =~ /\/$/;
    my @lines;	
    if (open FILE, $WorkingDir.$inputFilePath){
       @lines=<FILE>; 
	close (FILE);
    }
    if (open (NEWFILE,">$WorkingDir$inputFilePath")){
	my $line;
	foreach $line (@lines){
	   #     $line=~s/(\r)$/\n/;
		$line=~s/(\s+)$//;
          	print NEWFILE "$line\n";
        }
        close NEWFILE;
    }
}
#---------------------------------------------
sub check_file_type{
    
    my $FileName=shift;
    my $Type="PLAIN_TEXT";
    if (-e "$FileName")
    {   
        #$Type="Executable" if (-x $FileName); #Executable
        $Type="Binary" if (-c $FileName); #Contains Special Chars;
        $Type="Binary" if (-B $FileName); #Binary
        
        if (-T $FileName and $Type ne "BINARY") # Potentially Text File but maybe not: The first block or so of the file is examined for odd characters such as strange control codes or characters with the high bit set. If too many strange characters (>30%) are found, it's a -B  file; otherwise it's a -T  file...
        {
           unless (open FILE,$FileName){
                return ("ERR", "check_file_type : cannot open the file $FileName for reading $!");
           }
           my $line=<FILE>;
           close (FILE);
           if ($line=~/%PDF-/){
                $Type="PDF";
            }
           elsif ($line=~/\\rtf/){
                $Type="RTF";
            }
           
        }
    }
    else
    {
        return ("ERR", "check_file_type : the file $FileName was not found");
    }
    return ("OK", $Type);    
}
#---------------------------------------------
1;
