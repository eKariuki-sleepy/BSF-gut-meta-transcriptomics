library(microbiome)
library(ggpubr)
#calculatin an even sampling depth for all the samples
summary(sample_sums(physeq))#
#rarefied <- rarefy_even_depth(carbom, sample.size = 1648)
ps.meta <- meta(physeq)
tab <- microbiome::alpha(physeq, index = "diversity_shannon")
ps.meta$diversity_Shannon <- tab$diversity_shannon
tab1 <- microbiome::alpha(physeq, index = "Chao1")
ps.meta$Chao1 <- tab1$chao1#plotting the boxplots for the shannon index data

p4 <- ggboxplot(ps.meta, "DIET","Chao1",
                color = "DIET", palette = c("#F5EB02
", "#D36318
", "#120004
", "#FF00FF
", "#16A3A3
", "#0921FF
", "#FF0000
", "#19BA24
"), add = "jitter", linetype = "solid", Family = "arial", add.params = list(),
                error.plot = "pointrange", legand = NULL, size = NULL, width = 0.7, notch = FALSE, outlier.shape = 20, facet.by = NULL,
                panel.labs = NULL, short.panel.labs = TRUE,bxp.errorbar = FALSE, bxp.errorbar.width = 0.4, ggtheme = theme_pubr())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic", colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))
p4 + stat_compare_means()+aes(x = fct_inorder(DIET)) + theme(legend.position = "none") + xlab("DIET") + ylab("Chao1 diversity") + labs(tag = "E", plot.tag.position = c(0.2, -0.1))
