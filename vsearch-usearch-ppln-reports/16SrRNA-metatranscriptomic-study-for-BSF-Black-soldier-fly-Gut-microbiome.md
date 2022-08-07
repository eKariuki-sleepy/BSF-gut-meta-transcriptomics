# 16SrRNA metatranscriptomic study for BSF(Black soldier fly) Gut microbiome.

:::info
#### NOTES
<ins>Questions</ins>
1.Is it 16s only?
2.Word count 65046, report 270k.?
3.more data in new data than old data?

  *** 
- Long read qc -https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7144081/
- OTU picking options for sortmerna.Give a try:
http://manpages.ubuntu.com/manpages/bionic/man1/sortmerna.1.html
- Generate two pipelines 
Pipeline A - contains samples with read count greater that 50,000 i.e. Samples 8.1, 9.1, 7.1, 11.1 sorted by a 90% identity on the silver database and samples with read counts less than 50,000 i.e. the rest of the samples sorted by an 85% identity on the silver database. The two sets of samples will be merged and analysed using the Qiime Vsearch pipeline

    Pipeline B - Contains samples all the samples all of them sorted using a silver database that has an identity of 90%. This dataset will be analysed using the dada2 pipeline. The two results will be compared to see how they differ.
-Merge all replicates from one sample if the first two approaches fail
-Use sortme rna for OTU picking as an alternative to vsearch.




:::

:::info
#### WORKFLOW
:::
#### 1. Attain the data of 5 samples and its replicates

#### 2. Validated metadata using keemei

#### 3. Drafted a schematic workflow of the pipeline



#### 4. Creating and activating a working conda environment.
```
conda create metatranscriptomics
conda activate metatranscriptomics
```
#### 5. Filtering and sorting rRNA Sequences with SortMeRNA
Description

"SortMeRNA is a program tool for filtering, mapping and OTU-picking NGS reads in metatranscriptomic and metagenomic data. The core algorithm is based on approximate seeds and allows for fast and sensitive analyses of nucleotide sequences. The main application of SortMeRNA is filtering ribosomal RNA from metatranscriptomic data."

Once we have removed low quality sequences and remove any adapter contamination, we can then proceed to an additional (and optional) step to remove rRNA sequences from the samples. If your samples were not prepared with an rRNA depletion protocol before library preparation, it is reccomended to run this step to computational remove any rRNA sequence contiamation that may otheriwse take up a majority of the aligned sequences."

##### 5.1 Installation
```
conda install sortmerna
```
##### 5.2 Preparation
Before we can run the sortmerna command, we must first download the eukaryotic, archeal and bacterial rRNA databases. The sortmerna_db/ folder will be the location that we will keep the files necessary to run SortMeRNA. These databases only need to be created once, so any future RNAseq experiements can use these files.

##### 5.3 Download the sortmerna package (~2min) into sortmerna_db folder
```
    wget -P sortmerna_db https://github.com/biocore/sortmerna/archive/2.1b.zip

    # Decompress folder 
    unzip sortmerna_db/2.1b.zip -d sortmerna_db

    # Move the database into the correct folder
    mv sortmerna_db/sortmerna-2.1b/rRNA_databases/ sortmerna_db/

    # Remove unnecessary folders
    rm sortmerna_db/2.1b.zip
    rm -r sortmerna_db/sortmerna-2.1b

```
##### 5.4 Sorting 
```
**id85- < 50k reads
#This script takes in our fastq files and sorts the rRNA based on ribosomal subunits.

mkdir bac-16S
mkdir aligned

#The loop generates mapped rRNA reads ,unmapped rRNA reads with respect to each sample and a log file.
#Removes the kvdb file to allows next run.(sortmerna guidlines)

for file in ./*.fq;
do

        name=$(basename ${file} .fq)
       # echo $name
        sortmerna --ref ../sortmerna/sortmerna_db/rRNA_databases/silva-bac-16s-database-id85.fasta   --reads $file \
        --aligned ./aligned/bac-16S-85-mapped$name  --other ./bac-16S-85-umapped$name \
        --fastx ./bac-16S-85-unmapped$name  --threads 9
             rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
             mv ./aligned/bac-16S* bac-16S

done


**id90 > 50k reads.
#This script takes in our fastq files and sorts the rRNA based on ribosomal subunits.

mkdir bac-16S
mkdir aligned

#The loop generates mapped rRNA reads ,unmapped rRNA reads with respect to each sample and a log file.
#Removes the kvdb file to allows next run.(sortmerna guidlines)

for file in ./*.fq;
do

        name=$(basename ${file} .fq)
       # echo $name
        sortmerna --ref ../sortmerna/sortmerna_db/rRNA_databases/silva-bac-16s-database-id85.fasta   --reads $file \
        --aligned ./aligned/bac-16S-85-mapped$name  --other ./bac-16S-85-umapped$name \
        --fastx ./bac-16S-85-unmapped$name  --threads 9
             rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
             mv ./aligned/bac-16S* bac-16S

done
```


##### 5.4 Renaming identifier lines 

 ```
 awk '{if( (NR-1)%4 ) print; else printf("@CF3-%d\n",cnt++)}' bac-16S-85-mappedbarcode04_rna.fq > new.fastq
 ```


##### 5.5 merging data
```
cat *.fq > merged.fq
```



:::info
## vsearch-Usearch pipeline
**Quality filtering

 **  vsearch -fastq_filter merged-90-85.fq --fastq_stripleft 0 --fastq_stripright 0 -fastq_maxee 600 --fastaout QCd_merged.fa
 
 --ratained all sequences
 
 **Dereplication
 
 ** vsearch --derep_fulllength QCd_merged.fa -sizeout -relabel Uniq -output unique_seqs.fa
 
 --73 sequences discarded.
-- 286667 sequences retained.
 
 **clustering otus
 
 ** usearch -unoise3 unique_seqs.fa -zotus ASVs.fa -minsize 2
 
 **generating count table
 
 ** sed -i.tmp 's/Zotu/ASV_/' ASVs.fa
rm ASVs.fa.tmp
vsearch -usearch_global QCd_merged.fa --db ASVs.fa --id 0.85 --otutabout ASV_counts.txt

**Assigning taxonomy

** usearch -sintax ASVs.fa -db ../rdp_16s_v18.fa -tabbedout ASV_tax_raw.txt -strand both -sintax_cutoff 0.1

** sed -i.tmp 's/#OTU ID//' ASV_counts.txt
cp ASV_counts.txt R_working_dir/
#this script is used to convert the output by remving the digit values.
bash convert_usearch_tax.sh
cp ASV_tax.txt R_working_dir/

```
convert_usearch_tax.sh

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




:::

:::info
### R-RESULTS.

#### Taxonomy
##### Initial plots
![](https://i.imgur.com/8xTSmak.png)
![](https://i.imgur.com/WEIROL6.png)
#### New plot subset(Bacteroidales)
![](https://i.imgur.com/aySqzdf.png)
![](https://i.imgur.com/9DvuH0o.png)
#### New bacteriodales
![](https://i.imgur.com/qHXxywV.png)

#### Clostridales(Order)
![](https://i.imgur.com/Z9L4ZOO.png)
#### Clostridia(Class)
![](https://i.imgur.com/dTi8Hg1.png)

#### Updated taxa
![](https://i.imgur.com/k06qrgR.png)


#### Bact
![](https://i.imgur.com/ke17tdN.png)
#### BSG -bact
![](https://i.imgur.com/PRocDMS.png)

#### GH51
![](https://i.imgur.com/SKE9hm6.png)
#### clostdium
![](https://i.imgur.com/5iwxhc8.png)
#### clostdia
![](https://i.imgur.com/5Xl7YkA.png)
###


#### Alpha diversity
![](https://i.imgur.com/sloU2vb.png)
** Shannon * try
![](https://i.imgur.com/YLPpAde.png)

** Betterplots
![](https://i.imgur.com/We8dAdO.png)

** new plots
![](https://i.imgur.com/NuAgd6V.png)

#### PcoA
##### Phylum
![](https://i.imgur.com/n9mm2NC.png)
##### phylum new
![](https://i.imgur.com/eyYfssf.png)

#### Phylogeny
##### bacteroidales
![](https://i.imgur.com/1rceoWY.png)


#### Beta diversity
Bray–Curtis: The sum of lesser counts for species present in both communities divided by the sum of all counts in both communities. This can be thought of as a quantitative version of the Sørensen index.


![](https://i.imgur.com/BmsIgJh.png)
![](https://i.imgur.com/AYLiCaM.png)
![](https://i.imgur.com/ioSRDQf.png)
![](https://i.imgur.com/zTVVYzO.png)
#### new plot
![](https://i.imgur.com/ayuF8nB.png)



##### Notes
Weighted Unifrac: The fraction of the phylogenetic tree branch lengths shared by the two communities, weighted by the counts of organisms, so more abundant organisms have a greater influence.

Non-metric multidimensional scaling (NMDS) is an indirect gradient analysis approach which produces an ordination based on a distance or dissimilarity matrix. Unlike methods which attempt to maximise the variance or correspondence between objects in an ordination, NMDS attempts to represent, as closely as possible, the pairwise dissimilarity between objects in a low-dimensional space. Any dissimilarity coefficient or distance measure may be used to build the distance matrix used as input.
NMDS is a rank-based approach. This means that the original distance data is substituted with ranks. Thus, rather than object A being 2.1 units distant from object B and 4.4 units distant from object C, object C is the "first" most distant from object A while object C is the "second" most distant. While information about the magnitude of distances is lost, rank-based methods are generally more robust to data which do not have an identifiable distribution.

NMDS is a robust technique. It can:
* Tolerate missing pairwise distances
* Be applied to a (dis)similarity matrix built with any (dis)similarity measure and
* use quantitative, semi-quantitative, qualitative, or mixed variables
:::



:::info
## R-Scripts

```
# 16sRRNA analysis.

### 1.Alpha diversity.
### 2.Beta diversity.
### 3.PcoA.
### 4.Taxonomic ranking.

# 1. Set working directory.

setwd("~/analysis_16s/90-85-16s-data/usearch/vsearch/r_work2")

# 2.  loading Packages

library("phyloseq")
library("permute")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("tidyverse")

# 3. Reading data
count_tab <- read.table("ASV_counts.txt", header=T, row.names=1, check.names=F)

tax_tab <- as.matrix(read.table("ASV_tax_new.txt", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

sample_info_tab <- read.table("sample_info1.txt", header=T, row.names=1, check.names=F)
sample_info_tab

##modified sample data
sample_info_tab1 <- read.table("sample_info1_new.txt", header=T, row.names=1, check.names=F)




# 4.Take in

# first we need to get a sum for each ASV across all 4 blanks and all 16 samples
#blank_ASV_counts <- rowSums(count_tab[,1:4])
#sample_ASV_counts <- rowSums(count_tab[,5:15])

# now we normalize them, here by dividing the samples' total by 4 – as there are 4x as many samples (16) as there are blanks (4)
#norm_sample_ASV_counts <- sample_ASV_counts/2

# here we're getting which ASVs are deemed likely contaminants based on the threshold noted above:
#blank_ASVs <- names(blank_ASV_counts[blank_ASV_counts * 10 > norm_sample_ASV_counts])
#length(blank_ASVs) # this approach identified about 50 out of ~1550 that are likely to have orginated from contamination

# looking at the percentage of reads retained for each sample after removing these presumed contaminant ASVs shows that the blanks lost almost all of their sequences, while the samples, other than one of the bottom water samples, lost less than 1% of their sequences, as would be hoped
#colSums(count_tab[!rownames(count_tab) %in% blank_ASVs, ]) / colSums(count_tab) * 100

# now that we've used our extraction blanks to identify ASVs that were likely due to contamination, we're going to trim down our count table by removing those sequences and the blank samples from further analysis
#filt_count_tab <- count_tab[!rownames(count_tab) %in% blank_ASVs, -c(1:1)]
# and here make a filtered sample info table that doesn't contain the blanks
filt_sample_info_tab<-sample_info_tab1[-c(1:2), ]

filt_sample_info_tab1<-sample_info_tab1[-c(1:2),]

# and let's add some colors to the sample info table that are specific to sample types and characteristics that we will use when plotting things
# we'll color the water samples blue: 
filt_sample_info_tab$color[filt_sample_info_tab$diet == "water_hyacinth"] <- "blue"
# the biofilm sample a darkgreen:
filt_sample_info_tab$color[filt_sample_info_tab$diet == "chicken_manure"] <- "darkgreen"
# the basalts with highly altered, thick outer rinds (>1 cm) brown ("chocolate4" is the best brown I can find...):
filt_sample_info_tab$color[filt_sample_info_tab$diet == "chicken_feed"] <- "chocolate4"
# the basalts with smooth, glassy, thin exteriors black:
filt_sample_info_tab$color[filt_sample_info_tab$diet == "feed_mix"] <- "black"
# and the calcified carbonate sample an ugly yellow:
filt_sample_info_tab$color[filt_sample_info_tab$diet == "brewers_spent_grain"] <- "darkkhaki"
# and the calcified carbonate sample an red:
filt_sample_info_tab$color[filt_sample_info_tab$diet == "water_hyacinth"] <- "red"

# and now looking at our filtered sample info table we can see it has an addition column for color
filt_sample_info_tab




# 5.Alpha diversity.# working

## Phyloseq
## first we need to create a phyloseq object using our un-transformed count table
count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
sample_info_tab_phy <- sample_data(filt_sample_info_tab)
treefilename = system.file("clusters.txt",  package="phyloseq")

##modified sample data
sample_info_tab_phy1 <- sample_data(sample_info_tab1)

#AS_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
A_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy1)
physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy1, treefilename)
# and now we can call the plot_richness() function on our phyloseq object
#plot_richness(AS_physeq, color="diet", measures=c("Chao1", "Shannon")) + 
  #scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$diet)])) +
  #theme(legend.title = element_blank())

## alpha boxplot
library(phylosmith)
### Shannon *try
Shannon <-alpha_diversity_graph(AS_physeq, index = 'shannon',
                      treatment = c( 'diet'), subset = NULL, colors = 'default') + labs(x="Samples", title="shannon diversity")
### Simpson*try
alpha_diversity_graph(AS_physeq, index = 'simpson',
                      treatment = c( 'diet'), subset = NULL, colors = 'default')
## BEST plots#working
plot_richness(A_physeq, x="diet", measures=c("Observed", "Shannon")) + geom_boxplot()


# 6. PCOA 

# not working
# generating and visualizing the PCoA with phyloseq
#vst_pcoa <- ordinate(AS_physeq, method="MDS", distance="euclidean")
#eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
## plotting
#plot_ordination(AS_physeq, vst_pcoa, color="diet") + 
  #labs(col="Pylum") + geom_point(size=1) + 
  #geom_text(aes(label=rownames(filt_sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  #coord_fixed(sqrt(eigen_vals[12]/eigen_vals[1])) + ggtitle("PCoA") + 
  #scale_color_manual(values=unique(filt_sample_info_tab$color[order(filt_sample_info_tab$diet)])) + 
  #theme(legend.position="none")

### Working better

physeq.ord <- ordinate(AS_physeq, "RDA", "bray" )

plot_ordination(AS_physeq, physeq.ord, type="split", color="Genus",
                shape="diet", title="PcoA", label = "station") +
  geom_point(size=3)



# 7. Beta diversity.
library(ape)
library(ggforce)

ps_ord <- ordinate(A_physeq, method = "RDA", distance = "jsd")

p=plot_ordination(A_physeq, ps_ord, type = "samples", color = "diet")

##Polygons
p + geom_point(size = 5) + geom_polygon(aes(fill = diet))

### ellipse 
## gforce
p + geom_mark_ellipse() 



# 8. Taxonomy.
#poor plots.
plot_bar(AS_physeq, fill = "Order") +
  geom_bar(aes(color=Order), stat="identity", position="stack")


## NEW APPROACH.

# how you want to break things down will depend on your data and your questions, as usual
# for now let's just generate a table of proportions of each phylum, and breakdown the Proteobacteria to classes

phyla_counts_tab <- otu_table(tax_glom(AS_physeq, taxrank="Phylum")) # using phyloseq to make a count table that has summed all ASVs that were in the same phylum
phyla_tax_vec <- as.vector(tax_table(tax_glom(AS_physeq, taxrank="Phylum"))[,2]) # making a vector of phyla names to set as row names

rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

# we also have to account for sequences that weren't assigned any taxonomy even at the phylum level 
# these came into R as 'NAs' in the taxonomy table, but their counts are still in the count table
# so we can get that value for each sample by substracting the column sums of this new table (that has everything that had a phylum assigned to it) from the column sums of the starting count table (that has all representative sequences)
unclassified_tax_counts <- colSums(filt_count_tab) - colSums(phyla_counts_tab)
# and we'll add this row to our phylum count table:
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

# now we'll remove the Proteobacteria, so we can next add them back in broken down by class
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% " Proteobacteria", ]

class_counts_tab <- otu_table(tax_glom(AS_physeq, taxrank="Phylum")) # making count table broken down by class (contains classes beyond the Proteobacteria too at this point)

class_tax_phy_tab <- tax_table(tax_glom(AS_physeq, taxrank="Phylum")) # making a table that holds the phylum and class level info
phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
#
#rows_tmp <- row.names(class_tax_phy_tab)
#class_tax_tab <- data.frame("Phylum"=phy_tmp_vec, "Phylum"=class_tmp_vec, row.names = rows_tmp)
#proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$Class== "Proteobacteria", "Class"]) # making a vector of just the Proteobacteria classes
#proteo_classes_vec <- proteo_classes_vec[proteo_classes_vec != "NA"] # removing the NA (these are Proteobacteria not assigned anything at the class level), we'll add the counts back in for these in a second. leaving them here is not an option as there are multiple ones in our Class column (because they can come from any Phylum)

#rownames(class_counts_tab) <- as.vector(class_tax_tab$Class) # changing the row names like above so that they correspond to the taxonomy, rather than an ASV identifier
#proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] # making a table of the counts of the Proteobacterial classes

# there are also likely some some sequences that were resolved to the level of Proteobacteria, but not any further, and therefore would be missing from our class table
# we can find the sum of them by subtracting the proteo class count table from just the Proteobacteria row from the original phylum-level count table
#proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

# now combining the tables:
#major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_NA"=proteo_no_class_annotated_counts)
# and to check we didn't miss any other sequences, we can compare the column sums to see if they are the same:
#identical(colSums(major_taxa_counts_tab), colSums(filt_count_tab)) # if "TRUE", we know nothing fell through the cracks

# now we'll generate a proportions table for summarizing:
#major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

# if we check the dimensions of this table at this point
#dim(major_taxa_proportions_tab)
# we see there are currently 25 rows, which might be a little busy for a summary figure
# many of these taxa make up a very small percentage, so we're going to filter some out
# this is a completely arbitrary decision solely to ease visualization and intepretation, entirely up to your data and you
# here, we'll only keep rows (taxa) that make up greater than 5% in any individual sample
#temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 3, ])
# checking how many we have that were above this threshold
#dim(temp_filt_major_taxa_proportions_tab) # now we have 13, much more manageable for an overview figure

# though each of the filtered taxa made up less than 5% alone, together they can be more
# so we're going to add a row called "Other" that keeps track of how much we filtered out (which will also keep our totals at 100%)
#filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
#filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

## part 2
# first let's make a copy of our table that's safe for manipulating
#filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab
# and add a column of the taxa names so that it is within the table, rather than just as row names (this makes working with ggplot easier)
#filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)

# now we'll transform the table into narrow, or long, format (also makes plotting easier)
#filt_major_taxa_proportions_tab_for_plot.g <- gather(filt_major_taxa_proportions_tab_for_plot, Sample, Proportion, -Major_Taxa)

# take a look at the new table and compare it with the old one
#head(filt_major_taxa_proportions_tab_for_plot.g)
#head(filt_major_taxa_proportions_tab_for_plot)
# manipulating tables like this is something you may need to do frequently in R

# now we want a table with "color" and "characteristics" of each sample to merge into our plotting table so we can use that more easily in our plotting function
# here we're making a new table by pulling what we want from the sample information table
#sample_info_for_merge<-data.frame("Sample"=row.names(filt_sample_info_tab), "diet"=filt_sample_info_tab$diet, "color"=filt_sample_info_tab$color, stringsAsFactors=F)
# and here we are merging this table with the plotting table we just made (this is an awesome function!)
#filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)

# and now we're ready to make some summary figures with our wonderfully constructed table
## a good color scheme can be hard to find, i included the viridis package here because sometimes it's been really helpful for me, though this is not demonstrated in all of the following :/ 

# one common way to look at this is with stacked bar charts for each taxon per sample:
#ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  #geom_bar(width=0.6, stat="identity") +
  #scale_fill_viridis(discrete=TRUE) +
  #theme_bw() +
  #theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  #labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")

## Boxplots of the same.
#ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Major_Taxa, Proportion)) +
  #geom_jitter(aes(color=factor(diet)), size=2, width=0.15, height=0) +
  #scale_color_manual(values=unique(filt_major_taxa_proportions_tab_for_plot.g2$color[order(filt_major_taxa_proportions_tab_for_plot.g2$diet)])) +
  #geom_boxplot(fill=NA, outlier.color=NA) +
  #theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
 # labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples")

### Subseting Bacteroidales
gp1.ch = subset_taxa(A_physeq, Order == "Bacteroidales")
plot_bar(gp1.ch, fill="Genus")
plot_bar(gp1.ch, x="diet", fill="Genus")

gp2.ch = subset_taxa(A_physeq, Order =="Clostridiales")
gp3.ch = subset_taxa(A_physeq, Class =="Clostridia")
plot_bar(gp2.ch, fill="Genus")
plot_taxa_bar(gp2.ch, x="diet", fill="Genus")
plot_taxa_bar(gp3.ch, x="diet", fill="Genus")


##new phylogenetics
library(ape)
random_tree = rtree(ntaxa(A_physeq), rooted=TRUE, tip.label=taxa_names(A_physeq))
plot(random_tree)
##merge
physeq1 = merge_phyloseq(A_physeq,random_tree)
physeq1
gp4.ch = subset_taxa(physeq1, Order =="Bacteroidales")
plot_tree(gp4.ch, color="diet", shape="Order", label.tips="Genus",size="abundance", ladderize="left", plot.margin=0.3)

plot_tree(physeq1, color="diet", shape="Family", label.tips="Order", ladderize="left",size="abundance", plot.margin=0.3)

#size="abundance",shape="Genus"color="diet

## Taxonomy take three(best)
library(RColorBrewer)
library("ggsci")
library("colorspace")
library("dichromat")
library(ggsci)
library(scales)
#display.brewer.pal(n = 48, name = 'Paired')

cv3 <- plot_bar(A_physeq, fill = "Order") +
  geom_bar(aes(color=Order), colour=NA, stat="identity", position="stack") +
  scale_fill_igv()
cv3


```

:::


:::info
#### References
1.https://astrobiomike.github.io/amplicon/workflow_ex#analysis-in-r
2.https://drive5.com/usearch/manual/ex_hmp.html
3.https://www.sciencedirect.com/science/article/pii/S2001037019303745
:::

usearch -calc_distmx ASVs.fa -tabbedout mx.txt -maxdist 0.1 -termdist 0.1

#### palette
https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6405277/
