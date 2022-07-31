# make_DESeq_PCA.R
# Adapted from SAMSA (Westreich, 2016) by Eric G Kariuki, & Eric Njiraini (Oct 2021)
# Run with --help flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input"), type="character", default="A:/Analy2is/R_analysis/counts/Organism",
              help="Input directory", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="PCA_plot.tab",
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("USAGE: $ run_DESeq_stats.R -I working_directory/ -O save.filename -L level (1,2,3,4)")

# check for necessary specs
if (is.null(opt$input)) {
  print ("WARNING: No working input directory specified with '-I' flag.")
  stop()
} else {  cat ("Working directory is ", opt$input, "\n")
  wd_location <- opt$input
  setwd(wd_location)  }

cat ("Saving results as ", opt$out, "\n")
save_filename <- opt$out

cat ("Calculating DESeq results for hierarchy level ", opt$level, "\n")

# import other necessary packages
suppressPackageStartupMessages({
  library(DESeq2)
  library("pheatmap")
  library(ggplot2)
  library(factoextra)
  library(ggforce)
  library(ggrepel)
})

# GET FILE NAMES
control_files <- list.files(
  pattern = "control_*", full.names = T, recursive = FALSE)
control_names = ""
for (name in control_files) {
  control_names <- c(control_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
control_names <- control_names[-1]
control_names_trimmed = ""
for (name in control_names) {
  control_names_trimmed <- c(control_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
control_names_trimmed <- control_names_trimmed[-1]

exp_files <- list.files(
  pattern = "experimental_*", full.names = T, recursive = FALSE)
exp_names = ""
for (name in exp_files) {
  exp_names <- c(exp_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
exp_names <- exp_names[-1]
exp_names_trimmed = ""
for (name in exp_names) {
  exp_names_trimmed <- c(exp_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
exp_names_trimmed <- exp_names_trimmed[-1]

# READ IN FILES
# loading the control table
y <- 0
for (x in control_files) {
  y <- y + 1
  if (y == 1) {
    control_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(control_table) = c("DELETE", x, "V3")
    control_table <- control_table[,c(2,3)]      }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    temp_table <- temp_table[,c(2,3)]
    control_table <- merge(control_table, temp_table, by = "V3", all = T)  }
}
control_table[is.na(control_table)] <- 0
rownames(control_table) = control_table$V3
control_table_trimmed <- control_table[,-1]

# loading the experimental table
y <- 0
for (x in exp_files) {
  y <- y + 1
  if (y == 1) {
    exp_table <- read.table(file = x, header=F, quote = "", sep = "\t")
    colnames(exp_table) = c("DELETE", x, "V3")
    exp_table <- exp_table[,c(2,3)]  }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    exp_table <- merge(exp_table, temp_table[,c(2,3)], by = "V3", all = T)  }
}
exp_table[is.na(exp_table)] <- 0
rownames(exp_table) = exp_table$V3
exp_table_trimmed <- exp_table[,-1]

# getting the column names simplified
colnames(control_table_trimmed) = control_names_trimmed
colnames(exp_table_trimmed) = exp_names_trimmed

complete_table <- merge(control_table_trimmed, exp_table_trimmed, by=0, all = TRUE)
complete_table[is.na(complete_table)] <- 1
rownames(complete_table) <- complete_table$Row.names
complete_table <- complete_table[,-1]
completeCondition <- data.frame(condition=factor(c(
  rep(paste("control", 1:length(control_files), sep=".")),
  rep(paste("experimental", 1:length(exp_files), sep=".")))))
completeCondition1 <- t(completeCondition)
colnames(complete_table) <- completeCondition1
completeCondition2 <- data.frame(condition=factor(c(
  rep("control", length(control_files)),
  rep("experimental", length(exp_files)))))

#Add an extra field for SampleID_ EK

completeCondition2 = data.frame(row.names = c("CF1","CF3","CF4","BSG1","BSG2","BSG3", "CM1", "CM2", "CM3","FM1","FM2","FM3","WH1","WH2","WH3"),completeCondition2)
###convert rowname to column
completeCondition2 <- tibble::rownames_to_column(completeCondition2, "samples")
### Add new column with condition names
completeCondition2$conditions <- c("CF","CF","CF","BSG","BSG","BSG","CM","CM","CM","FM","FM","FM","WH","WH","WH")
dds <- DESeqDataSetFromMatrix(complete_table, completeCondition2, ~conditions)

dds <- DESeq(dds)
transformed_data <- rlog(dds, blind=FALSE)

# making the PCA plot

# calculate euclidean distances from the variance-stabilized data
dists <- dist(t(assay(transformed_data)))
PCAplot <- plotPCA(transformed_data, intgroup = "conditions", returnData = TRUE)
percentVar <- round(100 * attr(PCAplot, "percentVar"))


pcplot <- ggplot(PCAplot, aes(PC1, PC2, color=conditions)) +
    geom_point(size=2.5) + ##### Adding ellipses using ggforce ######
    ggforce::geom_mark_ellipse(aes(fill = conditions,
                               color = conditions)) +
    theme(legend.position = 'bottom') +
    geom_label_repel(aes(label=completeCondition2$samples), fill = "white", size=2.8, 
                     max.overlaps = Inf, fontface = 'bold') +
#    geom_text(aes(label=name), hjust=1, vjust=-1) +
    ggtitle("PCA Plot of control vs. experimental organism data") +
    theme(legend.position = "bottom") + xlim(-12,5)+ ylim(-7.5,6) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
  
 # coord_equal()
  
pcplot

###saving and finishing up####
cat ("Saving PCA plot as ", save_filename, " now.\n")
pdf(file = paste(save_filename,".pdf",sep = ""), width=10, height=7)
dev.off()
