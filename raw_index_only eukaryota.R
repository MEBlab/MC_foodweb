options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:/tap/R/MC/CO1/97%/new/new/total/beta/raw/newfig")
set.seed(123)


library("phyloseq")
library("xlsx")
library("tidyverse")
library("extrafont")
library("ggplot2")
library("ggsignif") 
library("patchwork")
library("ggforce")
library("ggnewscale")



otu_mat <- read.xlsx(file = "Fig1_raw_index.xlsx", 1)
tax_mat <- read.xlsx(file = "Fig1_raw_index.xlsx", 3)
samples_df <- read.xlsx(file = "Fig1_raw_index.xlsx", 2)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
MC_raw_index <- phyloseq(OTU, TAX, samples)


MC_ord <- ordinate(MC_raw_index, "PCoA", "bray")
plot_ord <- plot_ordination(MC_raw_index, MC_ord,  type = "samples") +
  geom_point(size = 9, aes(shape = Depth2, color = Station)) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 22),
        legend.title = element_text(family = "Tahoma", face = "bold", size = 20),
        legend.text = element_text(family = "Tahoma", size = 18))  +
  geom_mark_ellipse(aes(fill = Cove), alpha = 0.2,  con.color = "none") +
  scale_fill_manual(values=c("#D62728", "#1F77B4")) +
  scale_color_manual(labels =c("MC19", "MC03", "MC07", "MC10", "MC13"),
                     values = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")) +
  scale_shape_manual(labels = c("Surface", "Middle", "Bottom"),
                     values = c("Surface" = 19, "Middle" = 17, "Bottom" = 15 )) 
p1 <- plot_ord + guides(color = guide_legend(title= "Station"), shape = guide_legend(title = "Depth"))

p1

plot_ord <- plot_ordination(MC_raw_index, MC_ord,  type = "samples") +
  geom_point(size = 9, aes(shape = Depth2, fill = Station), colour = "black", stroke = 1) +
  scale_fill_manual(values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"),
                    labels =c("MC19", "MC03", "MC07", "MC10", "MC13")) + 
  scale_shape_manual(labels = c("Surface", "Middle", "Bottom"),
                     values = c("Surface" = 21, "Middle" = 24, "Bottom" = 22 )) + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 22),
        legend.title = element_text(family = "Tahoma", face = "bold", size = 20),
        legend.text = element_text(family = "Tahoma", size = 18))  +
  ggnewscale::new_scale_fill()  +
  geom_mark_ellipse(aes(fill = Cove), alpha = 0.2,  con.color = "none") +
  scale_fill_manual(values=c("#D62728", "#1F77B4")) 
 
p1 <- plot_ord + guides(color = guide_legend(title= "Station"), shape = guide_legend(title = "Depth"))

p1

# 18S
otu_mat <- read.xlsx(file = "Fig1_raw_index.xlsx", 4)
tax_mat <- read.xlsx(file = "Fig1_raw_index.xlsx", 5)
samples_df <- read.xlsx(file = "Fig1_raw_index.xlsx", 2)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
MC_18S_decon <- phyloseq(OTU, TAX, samples)
MC_18S_decon.melt = psmelt(MC_18S_decon)

# COI

otu_mat <- read.xlsx(file = "Fig1_raw_index.xlsx", 6)
tax_mat <- read.xlsx(file = "Fig1_raw_index.xlsx", 7)
samples_df <- read.xlsx(file = "Fig1_raw_index.xlsx", 2)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
MC_COI_decon <- phyloseq(OTU, TAX, samples)
MC_COI_decon.melt = psmelt(MC_COI_decon)

# 12S (MiF) 

otu_mat <- read.xlsx(file = "Fig1_raw_index.xlsx", 8)
tax_mat <- read.xlsx(file = "Fig1_raw_index.xlsx", 9)
samples_df <- read.xlsx(file = "Fig1_raw_index.xlsx", 2)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
MC_12S_decon <- phyloseq(OTU, TAX, samples)
MC_12S_decon.melt = psmelt(MC_12S_decon)


MC_18S_decon.melt = psmelt(MC_18S_decon)
MC_COI_decon.melt = psmelt(MC_COI_decon)
MC_12S_decon.melt = psmelt(MC_12S_decon)

estimate_richness(MC_18S_decon, measures = c("Observed", "shannon"))
estimate_richness(MC_COI_decon, measures = c("Observed", "shannon"))
estimate_richness(MC_12S_decon, measures = c("Observed", "shannon"))

MC_18S_decon.melt$Station <- factor(MC_18S_decon.melt$Station, levels = c("St. 19", "St. 03", "St. 07", "St. 10", "St. 13"))
MC_COI_decon.melt$Station <- factor(MC_COI_decon.melt$Station, levels = c("St. 19", "St. 03", "St. 07", "St. 10", "St. 13"))
MC_12S_decon.melt$Station <- factor(MC_12S_decon.melt$Station, levels = c("St. 19", "St. 03", "St. 07", "St. 10", "St. 13"))

MC_18S_decon.melt$Depth2 <- factor(MC_18S_decon.melt$Depth2, levels = c("Surface", "Middle", "Bottom"))
MC_COI_decon.melt$Depth2 <- factor(MC_COI_decon.melt$Depth2, levels = c("Surface", "Middle", "Bottom"))
MC_12S_decon.melt$Depth2 <- factor(MC_12S_decon.melt$Depth2, levels = c("Surface", "Middle", "Bottom"))




######### Plot  ######################
test <- read.xlsx(file = "Fig1_raw_index.xlsx", 14)
test$Station <- factor(test$Station, levels = c("MC19", "MC03", "MC07", "MC10", "MC13"))
test$Depth2 <- factor(test$Depth2, levels = c("Surface", "Middle", "Bottom"))


plot_alpha_COI <- ggplot(data = test, aes(x = Station, y = R_COI, fill=Station)) +
  geom_point(aes(fill = Station, shape = factor(Depth2)),  size = 5, stroke = 1) + 
  theme_bw() +
  ggtitle("COI") + 
  theme(text = element_text(family = "Tahoma"), 
        axis.text = element_text(size= 22), 
        axis.text.x = element_text(angle = 90, vjust = -0.01),
        axis.title = element_text(size = 25, face = "bold", family = "Tahoma"), 
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -0.5),
        axis.line = element_line(),
        legend.title = element_text(size = 20, face = "bold", family = "Tahoma"), 
        legend.text = element_text(size = 18, family = "Tahoma"), 
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold.italic", family = "Tahoma", hjust = 0.5)) +
  scale_fill_manual("Station", values = c("#F8766D", "#ABA300", "#0CB702", "#00A9FF", "#ED35ED")) +
  scale_shape_manual("Depth", values = c(22, 24, 21)) +
 # geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(0,300)) +
  ylab("Observed") + 
  xlab("")



plot_alpha_18S <- ggplot(data = test, aes(x = Station, y = R_18S, fill=Station)) +
  geom_point(aes(fill = Station, shape = factor(Depth2)),  size = 5, stroke = 1) + 
  theme_bw() +
  ggtitle("18S rRNA") + 
  theme(text = element_text(family = "Tahoma"), 
        axis.text = element_text(size= 22), 
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 25, face = "bold", family = "Tahoma"), 
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -0.5),
        axis.line = element_line(),
        legend.title = element_text(size = 20, face = "bold", family = "Tahoma"), 
        legend.text = element_text(size = 18, family = "Tahoma"), 
        legend.position = "none", 
        plot.title = element_text(size = 20, face = "bold.italic", family = "Tahoma", hjust = 0.5)) +
  scale_fill_manual("Station", values = c("#F8766D", "#ABA300", "#0CB702", "#00A9FF", "#ED35ED")) +
  scale_shape_manual("Depth", values = c(22, 24, 21)) +
#  geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(0,310)) +
  ylab("") 


plot_alpha_12S <- ggplot(data = test, aes(x = Station, y = R_12S, fill=Station)) +
  geom_point(aes(fill = Station, shape = factor(Depth2)),  size = 5, stroke = 1) + 
  theme_bw() +
  ggtitle("12S rRNA") + 
  theme(text = element_text(family = "Tahoma"), 
        axis.text = element_text(size= 22), 
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 25, face = "bold", family = "Tahoma"), 
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -0.5),
        axis.line = element_line(),
        legend.title = element_text(size = 20, face = "bold", family = "Tahoma"), 
        legend.text = element_text(size = 18, family = "Tahoma"),
        plot.title = element_text(size = 20, face = "bold.italic", family = "Tahoma", hjust = 0.5)) +
  scale_fill_manual("Station", values = c("#F8766D", "#ABA300", "#0CB702", "#00A9FF", "#ED35ED")) +
  scale_shape_manual("Depth", values = c(22, 24, 21)) +
#  geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(0,80), breaks = c(0,20,40,60,80)) +
  ylab("") + 
  xlab("")

# For legend
plot_alpha_12S <- ggplot(data = test, aes(x = Station, y = R_12S, fill=Station)) +
  geom_point(aes( shape = factor(Depth2)),  size = 2, stroke = 1) + 
  theme_bw() +
  ggtitle("12S rRNA") + 
  theme(text = element_text(family = "Tahoma"), 
        axis.text = element_text(size= 22), 
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 25, face = "bold", family = "Tahoma"), 
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -0.5),
        axis.line = element_line(),
        legend.title = element_text(size = 20, face = "bold", family = "Tahoma"), 
        legend.text = element_text(size = 18, family = "Tahoma"),
        plot.title = element_text(size = 20, face = "bold.italic", family = "Tahoma", hjust = 0.5)) +
  scale_shape_manual("Depth", values = c(21, 24, 22)) +
 # geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(0,80)) +
  ylab("") + 
  xlab("")

p2 <- plot_alpha_COI+plot_alpha_18S+plot_alpha_12S 

p2

p2/p1

ggsave("Abu_alpha diversity222.png", dpi=300, dev="png", height=12, width=12, units="in")

######## Shannon #####

plot_alpha_COI_S <- ggplot(data = test, aes(x = Station, y = S_COI, fill=Station)) +
  geom_point(aes(fill = Station, shape = factor(Depth2)),  size = 5, stroke = 1) + 
  theme_bw() +
  ggtitle("COI") + 
  theme(text = element_text(family = "Tahoma"), 
        axis.text = element_text(size= 22), 
        axis.text.x = element_text(angle = 90, vjust = -0.01),
        axis.title = element_text(size = 25, face = "bold", family = "Tahoma"), 
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -0.5),
        axis.line = element_line(),
        legend.title = element_text(size = 20, face = "bold", family = "Tahoma"), 
        legend.text = element_text(size = 18, family = "Tahoma"), 
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold.italic", family = "Tahoma", hjust = 0.5)) +
  scale_fill_manual("Station", values = c("#F8766D", "#ABA300", "#0CB702", "#00A9FF", "#ED35ED")) +
  scale_shape_manual("Depth", values = c(22, 24, 21)) +
  # geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(2,4)) +
  ylab("Shannon index") + 
  xlab("")



plot_alpha_18S_S <- ggplot(data = test, aes(x = Station, y = S_18S, fill=Station)) +
  geom_point(aes(fill = Station, shape = factor(Depth2)),  size = 5, stroke = 1) + 
  theme_bw() +
  ggtitle("18S rRNA") + 
  theme(text = element_text(family = "Tahoma"), 
        axis.text = element_text(size= 22), 
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 25, face = "bold", family = "Tahoma"), 
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -0.5),
        axis.line = element_line(),
        legend.title = element_text(size = 20, face = "bold", family = "Tahoma"), 
        legend.text = element_text(size = 18, family = "Tahoma"), 
        legend.position = "none", 
        plot.title = element_text(size = 20, face = "bold.italic", family = "Tahoma", hjust = 0.5)) +
  scale_fill_manual("Station", values = c("#F8766D", "#ABA300", "#0CB702", "#00A9FF", "#ED35ED")) +
  scale_shape_manual("Depth", values = c(22, 24, 21)) +
  #  geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(3,4)) +
  ylab("") 


plot_alpha_12S_S <- ggplot(data = test, aes(x = Station, y = S_12S, fill=Station)) +
  geom_point(aes(fill = Station, shape = factor(Depth2)),  size = 5, stroke = 1) + 
  theme_bw() +
  ggtitle("12S rRNA") + 
  theme(text = element_text(family = "Tahoma"), 
        axis.text = element_text(size= 22), 
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 25, face = "bold", family = "Tahoma"), 
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -0.5),
        axis.line = element_line(),
        legend.title = element_text(size = 20, face = "bold", family = "Tahoma"), 
        legend.text = element_text(size = 18, family = "Tahoma"),
        plot.title = element_text(size = 20, face = "bold.italic", family = "Tahoma", hjust = 0.5)) +
  scale_fill_manual("Station", values = c("#F8766D", "#ABA300", "#0CB702", "#00A9FF", "#ED35ED")) +
  scale_shape_manual("Depth", values = c(22, 24, 21)) +
  #  geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(0,3)) +
  ylab("") + 
  xlab("")

# FOr legend
plot_alpha_12S_S <- ggplot(data = test, aes(x = Station, y = S_12S, fill=Station)) +
  geom_point(aes( shape = factor(Depth2)),  size = 2, stroke = 1) + 
  theme_bw() +
  ggtitle("12S rRNA") + 
  theme(text = element_text(family = "Tahoma"), 
        axis.text = element_text(size= 22), 
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 25, face = "bold", family = "Tahoma"), 
        axis.title.y = element_text(vjust = 2), 
        axis.title.x = element_text(vjust = -0.5),
        axis.line = element_line(),
        legend.title = element_text(size = 20, face = "bold", family = "Tahoma"), 
        legend.text = element_text(size = 18, family = "Tahoma"),
        plot.title = element_text(size = 20, face = "bold.italic", family = "Tahoma", hjust = 0.5)) +
  scale_shape_manual("Depth", values = c(21, 24, 22)) +
  # geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(0,80)) +
  ylab("") + 
  xlab("")

p2_S <- plot_alpha_COI_S+plot_alpha_18S_S+plot_alpha_12S_S

p2_S

p2/p1

ggsave("Abu_alpha diversity222_S.png", dpi=300, dev="png", height=8, width=12, units="in")
