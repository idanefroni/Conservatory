#! /bin/sh
#######################################################################################################################
Help()
{
   # Display Help
   echo
   echo "This script performs a blastx query for CNS in an input file"
   echo "against a sample DB and filters CNS that are likely to be ORFs/TEs."
   echo "Takes sequences that are more than 30 bp and less than 4000 from START codon"
   echo "filters automatically sequences over 200 bp."
   echo "Divides sequences into 4 files during processing."
   echo	"Outputs:"
   echo "blast results"
   echo "csv of filtered CNS files"
   echo "csv with all entries that were dropped"
   echo
   echo "Syntax: scriptTemplate [-i -t -o]"
   echo "Arguments:"
   echo "i     Input file name."
   echo "t     Number of threads to use in blast."
   echo	"o     Output path (optional)."
   echo
}
######################################################################################################################

while getopts i:t:o::h flag; do
        case "${flag}" in
                i) input_file=${OPTARG};;
                t) threads=${OPTARG};;
                o) output_path=${OPTARG};;
                h) Help
                        exit;;
        esac
done

output_path=${output_path:-""}
file_name=${input_file##*/}
blast_file_name=$output_path"${file_name%.*}".blast.tsv

echo 'Intializing'
# Create fasta file of suspected ORFs
SEQS=$(mktemp /tmp/seqs-XXXX.tsv)
FASTA1=$(mktemp /tmp/orfs-XXXX.fasta)
FASTA2=$(mktemp /tmp/orfs-XXXX.fasta)
FASTA3=$(mktemp /tmp/orfs-XXXX.fasta)
FASTA4=$(mktemp /tmp/orfs-XXXX.fasta)
BLAST1=$(mktemp /tmp/blast-XXXX.out)
BLAST2=$(mktemp /tmp/blast-XXXX.out)
BLAST3=$(mktemp /tmp/blast-XXXX.out)
BLAST4=$(mktemp /tmp/blast-XXXX.out)


trap "rm -f $SEQS $FASTA1 $FASTA2 $FASTA3 $FASTA4 $BLAST1 $BLAST2 $BLAST3 $BLAST4" EXIT

echo 'Pulling sequences' 
awk -F "," 'function abs(v) {return v < 0 ? -v : v} length($15) > 30 && abs($8) < 4000 {print $2,"\t",$15}' $input_file > $SEQS
awk -v f1="$FASTA1" -v f2="$FASTA2" -v f3="$FASTA3" -v f4="$FASTA4" '
        {if (NR%4 == 0)
                {printf ">%s\n%s\n",$1,$2  > f1}
        else if (NR%4 == 1)
                {printf ">%s\n%s\n",$1,$2 > f2}
        else if (NR%4 == 2)
                {printf ">%s\n%s\n",$1,$2  > f3}
	else if (NR%4 == 3)
                {printf ">%s\n%s\n",$1,$2  > f4}}' 	$SEQS 

echo 'Running Blastx'
# Run Blastx
blastx -db /data/c/workspace2/adar/blast_filter/new_db/whole_db -query $FASTA1 -num_threads $threads -outfmt 6 -out $BLAST1 &
blastx -db /data/c/workspace2/adar/blast_filter/new_db/whole_db -query $FASTA2 -num_threads $threads -outfmt 6 -out $BLAST2 &
blastx -db /data/c/workspace2/adar/blast_filter/new_db/whole_db -query $FASTA3 -num_threads $threads -outfmt 6 -out $BLAST3 &
blastx -db /data/c/workspace2/adar/blast_filter/new_db/whole_db -query $FASTA4 -num_threads $threads -outfmt 6 -out $BLAST4 &


wait
# consolidate blast results
cat $BLAST1 >> $blast_file_name
cat $BLAST2 >> $blast_file_name
cat $BLAST3 >> $blast_file_name
cat $BLAST4 >> $blast_file_name

echo 'Processing results'
# Process results using R script
Rscript /data/c/workspace2/adar/blast_filter/process_blast.R $input_file $output_path
rm $FASTA1 $FASTA2 $FASTA3 $FASTA4 $SEQS $BLAST1 $BLAST2 $BLAST3 $BLAST4

echo 'Done'
