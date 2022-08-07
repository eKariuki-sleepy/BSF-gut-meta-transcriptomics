
# vsearch-Usearch pipeline

 ## Quality filtering

  vsearch -fastq_filter merged-90-85.fq --fastq_stripleft 0 --fastq_stripright 0 -fastq_maxee 600 --fastaout QCd_merged.fa
 
 --ratained all sequences
 
 ## Dereplication
 
  vsearch --derep_fulllength QCd_merged.fa -sizeout -relabel Uniq -output unique_seqs.fa
 
 --73 sequences discarded.
 -- 286667 sequences retained.
 
 ## Clustering otus
 
  usearch -unoise3 unique_seqs.fa -zotus ASVs.fa -minsize 2
 
 ## Generating count table
 
  sed -i.tmp 's/Zotu/ASV_/' ASVs.fa
  rm ASVs.fa.tmp
  vsearch -usearch_global QCd_merged.fa --db ASVs.fa --id 0.85 --otutabout ASV_counts.txt

 ## Assigning taxonomy

  usearch -sintax ASVs.fa -db ../rdp_16s_v18.fa -tabbedout ASV_tax_raw.txt -strand both -sintax_cutoff 0.1

  sed -i.tmp 's/#OTU ID//' ASV_counts.txt

  cp ASV_counts.txt R_working_dir/

 ## This script is used to convert the output by remving the digit values.

cp ASV_tax.txt R_working_dir/

```
nano convert_usearch_tax.sh

#!/bin/bash

cut -f1 ASV_tax_raw.txt > ASV_IDs
cut -f4 ASV_tax_raw.txt > raw_tax

cut -d : -f2 raw_tax | cut -d , -f1 | sed 's/\"//g' > domain
cut -d : -f3 raw_tax | cut -d , -f1 | sed 's/\"//g' > phylum
cut -d : -f4 raw_tax | cut -d , -f1 | sed 's/\"//g' > class
cut -d : -f5 raw_tax | cut -d , -f1 | sed 's/\"//g' > order
cut -d : -f6 raw_tax | cut -d , -f1 | sed 's/\"//g' > family
cut -d : -f7 raw_tax | cut -d , -f1 | sed 's/\"//g' > genus
cut -d : -f8 raw_tax | cut -d , -f1 | sed 's/\"//g' > species

paste ASV_IDs domain phylum class order family genus species > tax_temp

echo -e "\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" > header

cat header tax_temp > ASV_tax.txt

rm ASV_IDs raw_tax domain phylum class order family genus species tax_temp header

```

