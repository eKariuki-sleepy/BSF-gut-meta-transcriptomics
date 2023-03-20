library(ggplot2) # loading the library

df <- data.frame("feature.csv") # read in the data as data frame
df2$Sample.Id <- as.character(df2$Sample.Id) # Turn your 'Sample.Id' column into a character vector

df2$Sample.Id <- factor(df2$Sample.Id, levels=unique(df2$Sample.Id)) # Then turn it back into a factor with the levels in the correct order

ggplot(data = df2,mapping = aes(x = Sample.Id ,y = Feature.counts, fill = Sample.Id )) + geom_bar(stat = "identity") + ggtitle("Feature Counts graph") + theme(plot.title = element_text(hjust = 0.5))  # Generates the Barplot


