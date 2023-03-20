
###### PAIRED T-TEST FOR ERROR CORRECTION STATISTICS
###### by Eric G. Kariuki 17/11/21

library(dplyr)
library(ggpubr)

PATH <- "A:/Thesis/error_correction_table.csv"
df <- read.csv(PATH) 
# Calculating the mean and SD
library("dplyr")
group_by(df, Correction_status) %>%
  summarise(
    count = n(),
    mean = mean(Mapping_percentage, na.rm = TRUE),
    sd = sd(Mapping_percentage, na.rm = TRUE)
  )
#Plotting with Boxplots
library("ggpubr")
ggboxplot(df, x = "Correction_status", y = "Mapping_percentage", 
          color = "Correction_status", palette = c("#00AFBB", "#E7B800"),
          order = c("uncorrected", "corrected"),
          ylab = "Mapping percentage", xlab = "Correction status")

library(PairedData)
# Subset weight data before treatment
Uncorrected <- subset(df,  Correction_status == "uncorrected", Mapping_percentage,
                      drop = TRUE)
# subset weight data after treatment
Corrected <- subset(df,  Correction_status == "corrected", Mapping_percentage,
                    drop = TRUE)
# Plot paired data
library(PairedData)
pd <- paired(Uncorrected, Corrected)
plot(pd, type = "profile") + theme_bw()

#Alternative for paired data plot
ggpaired(df, x = "Correction_status", y = "Mapping_percentage",
         color = "Correction_status", line.color = "gray", line.size = 0.4,
         ylab = "Mapping percentage", xlab = "Correction status",
         palette = "jco") 
stat_compare_means(paired = TRUE)

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", paired = TRUE)


#Check for Normality using the Shapiro Test since n<30 (no. of samples)

# compute the difference
d <- with(df, 
          Mapping_percentage[Correction_status == "uncorrected"] - Mapping_percentage[Correction_status == "corrected"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.6141

# Computing t-test
res <- t.test(Mapping_percentage ~ Correction_status, data = df, paired = TRUE)
res

################################ MAPPING QUALITY PAIRED T-TEST #######################################

library(dplyr)
library(ggpubr)

PATH <- "A:/Thesis/error_correction_table.csv"
df <- read.csv(PATH) 
# Calculating the mean and SD
library("dplyr")
group_by(df, Correction_status) %>%
  summarise(
    count = n(),
    mean = mean(Mapping_qual., na.rm = TRUE),
    sd = sd(Mapping_qual., na.rm = TRUE)
  )
#Plotting with Boxplots
library("ggpubr")
ggboxplot(df, x = "Correction_status", y = "Mapping_qual.", 
          color = "Correction_status", palette = c("#00AFBB", "#E7B800"),
          order = c("uncorrected", "corrected"),
          ylab = "Mapping Quality", xlab = "Correction status")

library(PairedData)
# Subset weight data before treatment
uncorrected_mapq <- subset(df,  Correction_status == "uncorrected", Mapping_qual.,
                      drop = TRUE)
# subset weight data after treatment
corrected_mapq <- subset(df,  Correction_status == "corrected", Mapping_qual.,
                    drop = TRUE)
# Plot paired data
library(PairedData)
pd <- paired(uncorrected_mapq, corrected_mapq)
plot(pd, type = "profile") + theme_bw()

#Alternative for paired data plot
ggpaired(df, x = "Correction_status", y = "Mapping_qual.",
         font.label = list(size = 18, face = "bold", color = "black"),
         color = "Correction_status", line.color = "gray", line.size = 0.4,
         ylab = "Mapping Quality", xlab = "Correction status",
         palette = "jco") 
stat_compare_means(paired = TRUE)

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", paired = TRUE)


#Check for Normality using the Shapiro Test since n<30 (no. of samples)

# compute the difference
d <- with(df, 
          Mapping_qual.[Correction_status == "uncorrected"] - Mapping_qual.[Correction_status == "corrected"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.6141

# Computing t-test
res2 <- t.test(Mapping_qual. ~ Correction_status, data = df, paired = TRUE)
res2

# saving and finishing up
