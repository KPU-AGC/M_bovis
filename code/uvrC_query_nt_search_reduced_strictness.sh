#!/bin/bash
#This script will run BLAST searches for user inputed fasta files
#It was originally built for searching alignments of a specific primer set and probe for the *M. bovis* project
#It can readily be used for any BLAST searches using desired query sequences of interest
#To define queries just put your fasta files of interest in the "probe_primers_ref" directory
#BLAST parameters can easily be changed too by changing the flags used
#Refer to Table C1 in the Appendicies of the BLAST manual for all options and flags
#(https://www.ncbi.nlm.nih.gov/books/NBK279690/)
#-----------------------------------------------

#Load the conda blast environment
#--------------------------------

eval "$(conda shell.bash hook)"
conda activate ncbi-blast+

#Directories and files that need to be created or referred to
#------------------------------------------------------------
ref=~/Documents/paddy/m_bovis/ref
#nt_ref=~/Documents/ref/nt_database
uvrC_query=$ref/AF003959.1_uvrC_1kbp_up_500bp_dn_rev_comp.fasta
results=~/Documents/paddy/m_bovis/results
output_alignment=~/Documents/paddy/m_bovis/results/uvrC_output_alignment
fasta_format_alignment_long=~/Documents/paddy/m_bovis/results/uvrC_fasta_format_alignment_long
final_fasta=~/Documents/paddy/m_bovis/results/uvrC_final_fasta

uvrC_list=$ref/uvrC_list.txt
raw_alignment_list=$ref/raw_alignment_list.txt
fasta_list=$ref/fasta_list.txt
alignment_list=$ref/alignment_list.txt

#Make folders that may not already exist
#---------------------------------------
[ ! -d $output_alignment ] && mkdir -p $output_alignment
[ ! -d $fasta_format_alignment_long ] && mkdir -p $fasta_format_alignment_long
[ ! -d $final_fasta ] && mkdir -p $final_fasta

#BLAST search
#------------
blastn -db "nt" -query $uvrC_query -num_threads 12 -max_target_seqs 5000 -parse_deflines -evalue 5e-2 -perc_identity 75 -outfmt '6 sscinames scomnames staxids sallgi sacc sseqid evalue qcovs pident sstart send sseq' -out $output_alignment/uvrC_q_nt_s_reduced_strictness.txt

#Transform the tab delimited BLAST output into fasta formatted files
#-------------------------------------------------------------------
#Generate a list of output files to format
#-----------------------------------------
ls $output_alignment | grep .txt | sed "s/.txt//g" | uniq > $raw_alignment_list

#Format into fasta files
#-----------------------
while IFS= read -r prefix;
do
awk -F "\t" 'BEGIN { OFS = "\n"} {print ">" $1 "_" $3 ":" $6, $12}' $output_alignment/"$prefix".txt > $fasta_format_alignment_long/"$prefix"_long.fasta
fold -w 70 -s $fasta_format_alignment_long/"$prefix"_long.fasta > $final_fasta/"$prefix".fasta
done < $raw_alignment_list

#Remove duplicate results in fasta files (Clustal won't work with files that have duplicates)
#--------------------------------------------------------------------------------------------
#Generate a list of fasta files to work on
#-----------------------------------------
ls $final_fasta | grep .fasta | sed "s/.fasta//g" | uniq > $fasta_list

#Remove duplicated results in fasta files
#----------------------------------------
while IFS= read -r prefix;
do
seqkit rmdup -n $final_fasta/"$prefix".fasta -o $final_fasta/"$prefix"_clean.fasta
done < $fasta_list

#Replace spaces in identifiers with "_"
#aligners just will not handle them
#---------------------------------- 
sed -i 's/\s/_/g' $final_fasta/*_clean.fasta

#Clustal alignments
#---------------------------
#Generate a list of files for the alignment
#------------------------------------------
ls $final_fasta | grep _clean.fasta | sed "s/.fasta//g" | uniq > $alignment_list

#Finally, the last thing, the actual Clustal alignment
#-----------------------------------------------------
while IFS= read -r prefix;
do
clustalo -i $final_fasta/"$prefix".fasta -o $final_fasta/"$prefix".aln --outfmt=clu --seqtype=DNA
done < $alignment_list