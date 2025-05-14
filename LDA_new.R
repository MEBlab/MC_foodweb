options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:/tap/R/MC/CO1/97%/new/new/total")
set.seed(123)


library("phyloseq")
library("tidyverse")
library("ggplot2")
library("xlsx")
library("vegan")
library("RColorBrewer")
library("microbiomeMarker")
library("extrafont")
library("grid")


otu_mat <- read.xlsx(file = "Total decontam.xlsx", 2)
tax_mat <- read.xlsx(file = "Total decontam.xlsx", 4)
samples_df <- read.xlsx(file = "Total decontam.xlsx", 3)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
MC_total <- phyloseq(OTU, TAX, samples)


MC_Lefse <- run_lefse(MC_total, wilcoxon_cutoff = 0.05,  norm = "CPM", group = "Cove", taxa_rank = "Genus", kw_cutoff = 0.05, multigrp_strat = TRUE, lda_cutoff = 3 )
write.csv(file = "MC_Lefse.csv", MC_Lefse@marker_table)


MC_Lefse2 <- run_lefse(MC_total, wilcoxon_cutoff = 0.05,  norm = "CPM", group = "Cove", taxa_rank = "Genus", kw_cutoff = 0.05, multigrp_strat = TRUE, lda_cutoff = 3 )
write.csv(file = "MC_Lefse2.csv", MC_Lefse2@marker_table)

MC_LDA <- read.xlsx(file = "MC_Lefse2.xlsx", 3)


LDA_plot <- ggplot(data = MC_LDA) + 
  geom_bar( aes(x= reorder(Genus, -A), y= Score, fill = Cove), 
            color= "black", stat= "identity", alpha = 0.8) +
  theme_minimal() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.5), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18), 
        axis.text.y = element_text(face = "italic"), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = 3, linewidth = 0.5, color = "black"), 
        axis.ticks.x = element_line(), 
        legend.title = element_text(family = "Tahoma", face = "bold", size = 20), 
        legend.text = element_text(family = "Tahoma", size = 18)) + 
  coord_flip() + 
  xlab("Genus") + 
  scale_y_continuous(limits = c(0,5)) + 
  ylab(expression(bold("LDA score" ~  (log[10])))) + 
  scale_fill_manual(values = c("#D62728","#1F77B4"))

ggsave("LDA_plot -minor.png", dpi=300, dev="png", height=8, width=12, units="in")
