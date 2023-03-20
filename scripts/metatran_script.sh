#!/bin/bash

################################################
    ######### VARIABLES LIST##########

GuppyBinary="https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_5.0.11_linux64.tar.gz"
FAST5s=/home/ekariuki/ONTSeq/analysis/basecalling/raw_f5s/combf5s.zip
ONTdata=/home/ekariuki/ONTSeq/analysis/basecalling/raw_f5s
FASTQs=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/fq_files
mergfqs=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/fq_files/merg_combfqs
guppy_bc=./ont-guppy/bin/guppy_basecaller
lis=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/lis.txt
list=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/list.txt
cfg_file=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/dna_r9.4.1_450bps_hac_mod.cfg
SEQ_SUMMARY_2=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/sequencing_summary_run2.txt
SEQ_SUMMARY_1=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/sequencing_summary_run1.txt
porechop=/home/ekariuki/eanbitRT21/Porechop/results
RES=/home/ekariuki/ONTSeq/analysis/analysis2/results
NANOPLOT=/home/ekariuki/ONTSeq/analysis2/results/nanoplot
isonclust=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/fq_files/merg_combfqs/isonclust
isondir=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/fq_files/ison_correct
edfq=/home/ekariuki/ONTSeq/analysis/edlib_output
mergcor=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/fq_files/ison_correct/merg_corrected
rrnadb=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/rRNA_databases
silvref=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/data
sortmerna=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/sortmerna
sortmerna_uncor=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/sortmerna_uncor
genome=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/BSF_genome/GCF_905115235.1_iHerIll2.2.curated.20191125_genomic.fa
ind_genome=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/BSF_genome/GCF_905115235.1_iHerIll2.2.curated.20191125_genomic.fa.gmidx
gtf=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/BSF_genome/GCF_905115235.1_iHerIll2.2.curated.20191125_genomic.gtf
ref_trans=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/BSF_genome/GCF_905115235.1_iHerIll2.2.curated.20191125_rna.fna
mmap2=/home/ekariuki/ONTSeq/analysis/analysis2/results/minimap2
unm_fqs=/home/ekariuki/ONTSeq/analysis/analysis2/results/minimap2/corrected_folder/unmapped_fqs
mmap2_cor=/home/ekariuki/ONTSeq/analysis/analysis2/results/minimap2/corrected_folder
mmap2_uncor=/home/ekariuki/ONTSeq/analysis/analysis2/results/minimap2/uncorrected_folder
genome_mmap2=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/BSF_genome/GCF_905115235.1_iHerIll2.2.curated.20191125_genomic.fa
database=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/setup_and_test/RefSeq_bac.fa
diamond_database=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/setup_and_test/RefSeq_bac.fa
subsys_db=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/setup_and_test/subsys_db.fa
subsys_database=/home/ekariuki/samsa2/setup_and_test/subsys_db.fa
daa=/home/ekariuki/ONTSeq/analysis/daafiles
fast_list=/home/ekariuki/ONTSeq/analysis/analysis2/results/minimap2/corrected_folder/unmapped_fqs/list.txt
diam=/home/ekariuki/ONTSeq/analysis/analysis2/results/Diamond
py_scripts=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/python_scripts
R_scripts=/home/ekariuki/samsa2/R_scripts
PAT=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/path_list.txt
dbcan_fas=/home/ekariuki/ONTSeq/analysis/eanbitRT_22/blastn_out
dbase=/home/ekariuki/ONTSeq/analysis/run_dbcan/db
dbcan2=/home/ekariuki/ONTSeq/analysis/run_dbcan/results

	######## First things first; Please activate conda environment where all packages are installed .... #######
echo 'Activate conda package manager... '

conda activate transcriptomics
	######## Ready to go :) ########
#########################################################
#echo 'Downloading Guppy'

wget $GuppyBinary
#########################################################
#echo 'Unpack the guppy binary files'

tar -xvzf ont-guppy_5.0.11_linux64.tar.gz
#########################################################
#echo 'unzipping fast5 files ...'

unzip -q $FAST5s -d $ONTdata
#########################################################
# To perform the basecalling step you need to know the flowcell and ONT kit used to generate your
# fast5 files and select the appropriate config file
#########################################################
#echo 'performing basecalling ...'

#CPU (Slower)
guppy_basecaller --compress_fastq -i $ONTdata/ \
-s $FASTQs --cpu_threads_per_caller 16 --num_callers 4 \
-c $cfg_file

#GPU
#$guppy_bc -i $ONTdata -s $FASTQs  \
--recursive \
--config $cfg_file \
--gpu_runners_per_device 16 \
--cpu_threads_per_caller 2 \
--device cuda:0
--barcode_kits "SQK-PCB109"

########################################################
#echo 'Demultiplexing with Porechop  ...'

porechop -i input_reads.fastq.gz -b output_dir
########################################################
#echo 'Running pycoQC ...'

pycoQC -f $SEQ_SUMMARY_1 -o $RES/PyqoQc_run1.html

#echo 'Running pycoQC ...'

pycoQC -f $SEQ_SUMMARY_2 -o $RES/PyqoQc_run2.html
########################################################
#echo 'Convert basecalled fastqs.gz into one .gz file ...'

gunzip -c $FASTQs/*.gz | gzip > $RES/reads.fastq.gz
########################################################
#echo 'Running nanofilt ...'

gunzip -c $RES/reads.fastq.gz | NanoFilt -q 8 | gzip > $RES/filtered-reads.fastq.gz
########################################################
#echo 'Merging the fastq files into one fastq.gz file per barcode ...'

for ff in `cat $lis`
        do echo ==================
        find -L $FASTQs/comb_fqs/${ff} -type f -name "*.fastq.gz" -exec cat {} \; > $FASTQs/merg_combfqs/${ff}_merg.fastq.gz
done
########################################################
#echo "unzip the merged fastq files ..."
gunzip $FASTQs/merg_combfqs/*fastq.gz 
########################################################
#echo "Running Pychopper (edlib_backend)..."
for fs in `cat $lis`
do 
	cdna_classifier.py -t 32 -m edlib -r $RES/${fs}_report.pdf -u $mergfqs/${fs}_unclassified.fq -w $mergfqs/${fs}_rescued.fq $mergfqs/${fs}_merg.fastq $mergfqs/${fs}_full_length_output.fq
done
###############################################
#echo merge all fastqs

find -L $mergfqs/temp_fqs/ -type f -name "*.fastq" -exec cat {} \; > $mergfqs/temp_fqs/all_merg.fastq


echo "Combined pychopper stats..."
cdna_classifier.py -t 32 -m edlib -r $RES/allmerged_report.pdf -u $mergfqs/temp_fqs/allmerged_unclassified.fq -w $mergfqs/temp_fqs/allmerged_rescued.fq $mergfqs/temp_fqs/all_merg.fastq $mergfqs/temp_fqs/all_full_length_output.fq

###############################################
#echo 'Running Pychopper ...'

#echo 'phmm backend ...'

#cdna_classifier.py -t 24 -r report.pdf -u unclassified.fq -w rescued.fq input.fq full_length_output.fq

#echo 'edlib/parasail backend ...'
cdna_classifier.py -t 32 -m edlib -k PCB109 -r $FASTQs/report.pdf -u $FASTQs/unclassified.fq -w $FASTQs/rescued.fq $FASTQs/pass/unclassified/*.fastq $FASTQs/full_length_output.fq
###############################################
#echo 'Running Porechop...'

echo 'Adapter trimming with porechop... '
for fs in `cat $lis`
do
	porechop -i $mergfqs/${fs}_full_length_output.fq -o $mergfqs/${fs}_adfree.fastq.gz --threads 32
done
###############################################
#echo 'Unzip adapter trimmed data...'

for fs in `cat $lis`
do
	gunzip $mergfqs/${fs}_adfree.fastq.gz
done
###############################################
#echo "Nanoplot on uncorrected reads... "

mkdir -p $RES/nanoplot_uncorrected
for f in `cat $list`
do
        NanoPlot -t 32 --fastq $mergfqs/*adfree.fastq --N50 --maxlength 40000 -o $RES/nanoplot_uncorrected --plots dot --legacy hex
done

#echo "Finished Nanoplot"
###############################################
      ##### Clustering of Isoforms ######
#echo 'Running isonclust ...'

for fs in `cat $lis`
do
        isONclust --ont --t 32 --fastq $mergfqs/${fs}_adfree.fastq --outfolder $isonclust/${fs}_clust.fastq
done

#echo 'obtain separate cluster fastq files ...'

for fs in `cat $list` 
do
	isONclust write_fastq --clusters $isonclust/${fs}_clust.fastq/final_clusters.tsv --fastq $mergfqs/${fs}_adfree.fastq --outfolder $isonclust/${fs}_clust.fastq/fastq_clusters --N 1
done
###############################################
 ####### Correcting of clustered reads ######

echo "Running isONcorrect ..."
for fs in `cat $list`
do
	clusfqs=$isonclust/${fs}_clust.fastq/fastq_clusters 
	run_isoncorrect --t 32 --fastq_folder $clusfqs --outfolder $isondir/${fs}_correction/
done

echo "Finished isONcorrect ... "

###############################################
#echo "Merging corrected fastqs ... "

mkdir -p $isondir/merg_corrected

for fs in `cat $list`
do
	cor=$isondir/"${fs}"_correction
	touch $isondir/$mergcor/"${fs}"_corrected_reads.fastq

	for f in $cor/*/corrected_reads.fastq
		do
		cat "${f}" >> $mergcor/"${fs}"_corrected_reads.fastq
	done
done
echo 'Reads merged ...'
###############################################
#echo "Nanoplot Quality stats on corrected reads... "
mkdir -p $RES/nanoplot.merged
for f in `cat $list`
do
 	NanoPlot -t 32 --fastq $mergcor/*.fastq --N50 --maxlength 40000 -o $RES/nanoplot.merged --plots dot --legacy hex
done

echo "Finished Nanoplot"
###############################################
echo "Running sortmerna to remove rRNA... "

######### SORTMERNA DEFAULT FAST PARAMETERS ###############

#sortmerna -ref ./data/silva-bac-16s-database-id85.fasta \
#  -ref ./data/silva-arc-16s-database-id95.fasta \
#  -reads /home/ekariuki/ONTSeq/isONcorrect/set2/"${r}"_corrected_reads.fastq \
#  -sam -fastx -blast 1 -num_alignments 1 -v

########## CUSTOMIZED PARAMETERS ################
mkdir -p $sortmerna
for r in `cat $list`
do
	sortmerna -ref $silvref/silva-bac-16s-database-id85.fasta \
	-ref $silvref/silva-arc-16s-database-id95.fasta \
	-ref $rrnadb/silva-euk-28s-id98.fasta \
	-ref $rrnadb/silva-euk-18s-id95.fasta \
  	-ref $rrnadb/silva-arc-23s-id98.fasta \
	-ref $rrnadb/rfam-5.8s-database-id98.fasta \
 	-reads $mergcor/"${r}"_corrected_reads.fastq \
	--aligned $sortmerna/"${r}"_rna --other $sortmerna/"${r}"_clean --threads 32 \
	-sam -fastx -blast 1 -num_alignments 1 -v

	#Empty the kdvb folder after each run
	rm -rfv /home/ekariuki/sortmerna/run/kvdb 

done
########################
#This step is necessary to evaluate the efficacy of the isONcorrect error correction step
#by comparing the output with that of clustered corrected reads.

echo 'running sortmerna on unclustered/uncorrected reads...'
mkdir -p $sortmerna_uncor
for r in `cat $list`
do
       sortmerna -ref $silvref/silva-bac-16s-database-id85.fasta \
       -ref $silvref/silva-arc-16s-database-id95.fasta \
       -ref $rrnadb/silva-euk-28s-id98.fasta \
       -ref $rrnadb/silva-euk-18s-id95.fasta \
       -ref $rrnadb/silva-arc-23s-id98.fasta \
       -ref $rrnadb/rfam-5.8s-database-id98.fasta \
       -reads $mergfqs/${r}_adfree.fastq \
       --aligned $sortmerna_uncor/"${r}"_rna --other $sortmerna_uncor/"${r}"_clean --threads 32 \
       -sam -fastx -blast 1 -num_alignments 1 -v

        #Empty the kdvb folder after each run
       rm -rfv /home/ekariuki/sortmerna/run/kvdb

done

##################################################
echo 'change filename extension ...'

for fs in `cat $list`
do
	mv $sortmerna/${fs}_clean.fq $sortmerna/${fs}_clean.fastq
done 

############# RUNNING MINIMAP2 ALIGNER ###########
echo 'Aligning corrected reads with Minimap2 ...'

for fs in `cat $list`
do
        minimap2 -t 16 -G 500k -k 13 -w 5 -ax splice $genome_mmap2 $sortmerna/${fs}_clean.fastq > $mmap2/${fs}_cor.sam
done
####################
echo 'Aligning unclustered/uncorrected reads with Minimap2 ...'

for fs in `cat $list`
do
        minimap2 -t 16 -G 500k -k 13 -w 5 -ax splice $genome_mmap2 $sortmerna_uncor/"${fs}"_clean.fq > $mmap2/${fs}_uncor.sam
done

echo 'alignment complete!'

##################################################

echo 'Samtools suite on Minimap2 Results... '

for b in `cat $list`
do
	echo ' Generating .bam files ... '
	
	samtools view -b $mmap2/${b}_cor.sam > $mmap2_cor/"${b}"_cor.bam  #Corrected
	samtools view -b $mmap2/${b}_uncor.sam > $mmap2_uncor/"${b}"_uncor.bam  #Uncorrected

	echo 'Sorting the .bam file...'

	samtools sort $mmap2_cor/"${b}"_cor.bam -o $mmap2_cor/"${b}"_sorted_cor.bam #Corrected
	samtools sort $mmap2_uncor/"${b}"_uncor.bam -o $mmap2_uncor/"${b}"_sorted_uncor.bam #Uncorrected

	echo 'Indexing the sorted bam file...'

	samtools index $mmap2_cor/"${b}"_sorted_cor.bam #Corrected
	samtools index $mmap2_uncor/"${b}"_sorted_uncor.bam #Uncorrected

	echo 'Run flagstat ... '

	samtools flagstat $mmap2_cor/"${b}"_sorted_cor.bam  >> $mmap2_cor/cor.tsv
	samtools flagstat $mmap2_uncor/"${b}"_sorted_uncor.bam  >> $mmap2_cor/uncor.tsv

	echo 'Obtain alignment statistics ...'

	samtools stats $mmap2_cor/"${b}"_sorted_cor.bam #Corrected
	samtools stats $mmap2_uncor/"${b}"_sorted_uncor.bam #Uncorrected

	echo 'Obtain coverage statistics ...'

	samtools coverage $mmap2_cor/"${b}"_sorted_cor.bam -o $mmap2_cor/${b}_cor_coverage.out # Corrected
	samtools coverage $mmap2_uncor/"${b}"_sorted_uncor.bam -o $mmap2_uncor/${b}_uncor_coverage.out # Uncorrected

	echo 'Getting unmapped reads in .bam file...'

	samtools view -b -f 4 $mmap2_cor/"${b}"_sorted_cor.bam > $mmap2_cor/"${b}"_unmapped_cor.bam

	echo 'Indexing the unmapped reads'

	samtools index $mmap2_cor/"${b}"_unmapped_cor.bam

	echo 'Getting only the mapped reads in .bam format... '

	samtools view -b -F 4 $mmap2_cor/"${b}"_sorted_cor.bam > $mmap2_cor/"${b}"_mapped_cor.bam

	echo 'Indexing the mapped reads... '

	samtools index $mmap2_cor/"${b}"_mapped_cor.bam

	echo 'Sorting the bam files... '
	samtools sort $mmap2_cor/"${b}"_unmapped_cor.bam -o $mmap2_cor/"${b}"_unm_sorted_cor.bam

	samtools sort $mmap2_cor/"${b}"_mapped_cor.bam -o $mmap2_cor/"${b}"_map_sorted_cor.bam

        ########## CONVERTING THE READS TO FASTQ ############
	echo 'Using samtools to convert bam files to fastq format... '

	samtools bam2fq $mmap2_cor/"${b}"_unm_sorted_cor.bam > $mmap2_cor/"${b}"_Unmapped_cor.fastq

	samtools bam2fq $mmap2_cor/"${b}"_map_sorted_cor.bam > $mmap2_cor/"${b}"_Mapped_cor.fastq
done

####################################################################
# DIAMOND_example_script.bash by Sam Westreich, github.com/transcript
(Westreich, 2018)
####################################################################
#
# This bash script shows how to structure commands for the following
# purposes:
#
#       1. Converting a reference database to DIAMOND-searchable format;
#       2. Performing an annotation search, running an input file against
#          a DIAMOND-searchable database;
#       3. Converting the results of an annotation search to BLAST m8
#      format.
####################################################################
echo 'QUERY DATABASE AND CREATE .daa FILES FOR INPUT FASTQS ...'
###################### REFSEQ DATABASE #############################
echo 'creating .daa files for RefSeq database ...'
for d in `cat $list`
do
	
	diamond blastx --query $mmap2_cor/unmapped_fqs/"${d}"_Unmapped_cor.fastq \
	--db $database --daa $daa/"${d}"_reads_unm.daa
done
####################### SUBSYSTEMS DATABASE#########################
echo ' creating .daa files for subsystems database ...'
for d in `cat $list`
do

	diamond blastx --query $mmap2_cor/"${d}"_Unmapped_cor.fastq \
	--db $subsys_db --daa $daa/"${d}"_subsys_reads.daa
done
####################################################################
# DATABASE CREATION

echo "Creating database ... " 

	######################## REFSEQ #######################
diamond makedb --in $database --db $database

	###################### SUBSYSTEMS #####################
diamond makedb --in $subsys_db --db $subsys_db
	#######################################################
# explanation of settings:

#       --in    starting file that will be used to create DIAMOND-searchable db
#       --db    name of created DIAMOND database (automatically appended with .daa)

####################################################################

# ANNOTATION SEARCH
	####################################
echo "Performing annotation search on  $filename against RefSeq database ... "

for d in `cat $list`
do 
	filename=$mmap2_cor/unmapped_fqs/"${d}"_Unmapped_cor.fastq
	diamond_output=/home/ekariuki/ONTSeq/analysis/daafiles/"${d}"_unmreads_refseq.daa
	diamond blastx --db $diamond_database -q $filename -a $diamond_output --sensitive -t ./ -k 1
done
	####################################
echo "Performing annotation search on  $filename  against  Subsystems database ... "

for d in `cat $list`
do
       filename=/home/ekariuki/ONTSeq/analysis/analysis2/results/minimap2/corrected_folder/"${d}"_Unmapped_cor.fastq
       diamond_output2=/home/ekariuki/ONTSeq/analysis/daafiles/"${d}"_subsys_reads.daa
       diamond blastx --db $subsys_database -q $filename -a $diamond_output2 --sensitive -t ./ -k 1
done
	###################################
# explanation of settings:
#
#       --db    sets database (must be in DIAMOND-readable form)
#       -q              query file name
#       -a              name of results file (in DIAMOND format)
#       -t              sets temporary directory locations
#       -k              number of hits above cutoff to return (if run without specifying,
#                       default is 25)

###############################################################

# CONVERTING RESULTS

echo "Converting file " $diamond_output " to readable format"
        ######################## REFSEQ #######################

for d in `cat $list`
do 
	diamond_output=$daa/"${d}"_unmreads_refseq.daa
	final_output=$diam/Refseq/"${d}"_BLAST_results_unmapped_RefSeq.m8
	diamond view -a $diamond_output -o $final_output -f tab
done

	 ###################### SUBSYSTEMS #####################
for d in `cat $list`
do
       diamond_output2=/home/ekariuki/ONTSeq/analysis/daafiles/"${d}"_subsys_reads.daa
       final_output2=$diam/Subsystems/"${d}"_BLAST_results_unmapped_subsys.m8
       diamond view -a $diamond_output2 -o $final_output2 -f tab
done
	########################################################
# explanation of settings:
#
#       --daa(-a)       name of output file in .daa format, from the annotation search
#       -o              output file name
#       -f              separator of values in final output

#################################################################

####################### AGGREGATION STEP ########################
echo 'Aggregation step ... '
### ADD TABS TO THE DATABASES IF YOU ENCOUNTER AN INDEX ERROR ###

#sed -e 's/[:blank:] +/\t//g' $REF_DB > $edited.db.fa

#################### SUBSYSTEMS DATABASE ########################

echo 'Running subsystems database aggregation ...'

for d in `cat $list`
do 
	final_output2=$diam/Subsystems/"${d}"_BLAST_results_unmapped_subsys.m8
	python $py_scripts/DIAMOND_subsystems_analysis_counter.py \
	-I $final_output2 -D $subsys_db -O $diam/Subsystems/"${d}"_hierarchy -P $diam/Subsystems/"${d}"_unmap_receipt
done
################### SUMMARIZING SUBSYS. OUTPUT ###################
echo 'Reducing identical annotations... '

for d in `cat $list`
o
	#This quick program reduces down identical hierarchy annotations
	python $py_scripts/subsys_reducer.py -I $diam/Subsystems/"${d}"_hierarchy
done
####################### REFSEQ DATABASE ##########################
echo 'aggregating the RefSeq database counts by organism... '

for d in `cat $list`
do
	final_output=$diam/Refseq/"${d}"_BLAST_results_unmapped_RefSeq.m8
	python $py_scripts/standardized_DIAMOND_analysis_counter.py \
	-I $final_output -D $database -O
done

echo 'aggregating the RefSeq database counts by function... '

for d in `cat $list`
do
	final_output=$diam/Refseq/"${d}"_BLAST_results_unmapped_RefSeq.m8
	python $py_scripts/standardized_DIAMOND_analysis_counter.py \
	-I $final_output -D $database -F
done
###############################################################
echo 'aggregating counts per specific organism... '

for d in `cat $list`
do
        final_output=$diam/Refseq/"${d}"_BLAST_results_unmapped_RefSeq.m8
	python $py_scripts/DIAMOND_specific_organism_retriever.py -I \
	$final_output -SO Dysgonomonas -D $database
done
	### This can be done on any specific organism ###
for d in `cat $list`
do
        final_output=$diam/Refseq/"${d}"_BLAST_results_unmapped_RefSeq.m8
        python $py_scripts/DIAMOND_specific_organism_retriever.py -I \
        $final_output -SO Bacteroides -D $database
################################################################
echo 'obtaining raw read counts ...'

for d in `cat $list`
do
	
	python $py_scripts/raw_read_counter.py -I $mmap2_cor/unmapped_fqs/"${d}"_Unmapped_cor.fastq -O \
	$RES/rawcounts.txt
done
################################################################
echo 'Running DIAMOND Refseq Organism Annotation on pooled fasta files from the same condition.. '

for d in `cat $fast_list`
do
	diamond blastx -d $database -q $mmap2_cor/unmapped_fqs/${d}.fa -outfmt \
	'6 qseqid sseqid sseq evalue' -o $diam/${d}_blast.tsv
	diamond blastx -d $database -q $mmap2_cor/unmapped_fqs/${d}.fa -o $diam/${d}_matches.tsv
done
################################################################
echo 'Running DIAMOND subsystems Annotation on pooled fasta files from the same condition.. '

for d in `cat $fast_list`
do
	diamond blastx -d $subsys_db -q $mmap2_cor/unmapped_fqs/${d}.fa -outfmt \
	'6 qseqid sseqid sseq evalue' -o $diam/${d}_sub_blast.tsv
	diamond blastx -d $subsys_db -q $mmap2_cor/unmapped_fqs/${d}.fa -o $diam/${d}_subsys.tsv
done
################################################################
echo 'Aggregating counts by function per condition i.e. dietary group... '

for d in `cat $fast_list`
do
       final_output=$diam/"${d}"_matches.tsv
       python $py_scripts/standardized_DIAMOND_analysis_counter.py \
       -I $final_output -D $database -F
done

################################################################
echo 'Aggregating counts by organism per condition i.e. dietary group... '

for d in `cat $fast_list`
do
       final_output=$diam/"${d}"_matches.tsv
       python $py_scripts/standardized_DIAMOND_analysis_counter.py \
       -I $final_output -D $database -O
done

#################### SUBSYSTEMS DATABASE ######################

echo 'Running subsystems database aggregation ...'

for d in `cat $fast_list`
do 
	final_output2=$diam/${d}_subsys.tsv
	python $py_scripts/DIAMOND_subsystems_analysis_counter.py \
	-I $final_output2 -D $subsys_db -O $diam/Subsystems/"${d}"_merg_hierarchy -P $diam/Subsystems/"${d}"_merg_receipt
done
################### SUMMARIZING SUBSYS. OUTPUT ###################
echo 'Reducing identical annotations... '

for d in `cat $fast_list`
do
	#This quick program reduces down identical hierarchy annotations
	python $py_scripts/subsys_reducer.py -I $diam/Subsystems/"${d}"_merg_hierarchy
done

################################################################
##### Next: Statistical Analysis in R #####

# Find Scripts in R folder... :)

################################################################
echo 'Running dbCAN to find and annotate CAZymes ...'
#This step can alternatively be run on the dbCAN GUI/web server
#This step assumed you have already pre-installed dbCAN2 and its dependencies
#The step also assumes you have merged all samples from one condition and converted fastqs to .fa files

for f in `cat $PAT`
do
	# Make directories if they don't already exist
	mkdir -p $dbcan2/${f}

	# Run dbCAN using the Hotpep backend  
	run_dbcan.py $dbcan_fas/${f}.fa meta --db_dir $dbase --dbCANFile dbCAN-HMMdb-V8.txt -t hotpep \
	--hotpep_cpu 20 --out_dir $dbcan2/${f}/
done
#################################################################

# The end... For now ;)
