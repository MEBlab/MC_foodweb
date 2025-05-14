options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:/tap/R/MC/CO1/97%/new/new/total/beta")
set.seed(123)


library("phyloseq")
library("xlsx")
library("tidyverse")
library("extrafont")
library("ggplot2")
library("ggsignif") 
library("patchwork")
library("ggforce")

#Data input (Beta diversity combined loci)

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

MC_ord <- ordinate(MC_total, "PCoA", "bray")
plot_ord <- plot_ordination(MC_total, MC_ord,  type = "samples") +
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
  scale_color_manual(labels =c("St. 19", "St. 03", "St. 07", "St. 10", "St. 13"), 
                     values = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")) +
  scale_shape_manual(labels = c("Surface", "Middle", "Bottom"), 
                     values = c("Surface" = 19, "Middle" = 17, "Bottom" = 15 ))
p1 <- plot_ord + guides(color = guide_legend(title= "Station"), shape = guide_legend(title = "Depth"))

p1

ggsave("Index_beta diversity.png", dpi=300, dev="png", height=8, width=12, units="in")


###Each primer data input for alpha diversity (FIG 2A)

# 18S
otu_mat <- read.xlsx(file = "Total input.xlsx", 1)
tax_mat <- read.xlsx(file = "Total input.xlsx", 3)
samples_df <- read.xlsx(file = "Total input.xlsx", 2)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
MC_18S_decon_p <- phyloseq(OTU, TAX, samples)
MC_18S_decon_p.melt = psmelt(MC_18S_decon_p)

# COI

otu_mat <- read.xlsx(file = "Total input.xlsx", 4)
tax_mat <- read.xlsx(file = "Total input.xlsx", 5)
samples_df <- read.xlsx(file = "Total input.xlsx", 2)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
MC_COI_decon_p <- phyloseq(OTU, TAX, samples)
MC_COI_decon_p.melt = psmelt(MC_COI_decon_p)

# 12S (MiF) (Diatom included)

otu_mat <- read.xlsx(file = "Total input.xlsx", 6)
tax_mat <- read.xlsx(file = "Total input.xlsx", 7)
samples_df <- read.xlsx(file = "Total input.xlsx", 2)
otu_mat <- otu_mat %>% tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
MC_12S_decon_p <- phyloseq(OTU, TAX, samples)
MC_12S_decon_p.melt = psmelt(MC_12S_decon_p)


MC_18S_decon_p.melt = psmelt(MC_18S_decon_p)
MC_COI_decon_p.melt = psmelt(MC_COI_decon_p)
MC_12S_decon_p.melt = psmelt(MC_12S_decon_p)

estimate_richness(MC_18S_decon_p, measures = c("Observed"))
estimate_richness(MC_COI_decon_p, measures = c("Observed"))
estimate_richness(MC_12S_decon_p, measures = c("Observed"))

MC_18S_decon_p.melt$Station <- factor(MC_18S_decon_p.melt$Station, levels = c("St. 19", "St. 03", "St. 07", "St. 10", "St. 13"))
MC_COI_decon_p.melt$Station <- factor(MC_COI_decon_p.melt$Station, levels = c("St. 19", "St. 03", "St. 07", "St. 10", "St. 13"))
MC_12S_decon_p.melt$Station <- factor(MC_12S_decon_p.melt$Station, levels = c("St. 19", "St. 03", "St. 07", "St. 10", "St. 13"))

MC_18S_decon_p.melt$Depth2 <- factor(MC_18S_decon_p.melt$Depth2, levels = c("Surface", "Middle", "Bottom"))
MC_COI_decon_p.melt$Depth2 <- factor(MC_COI_decon_p.melt$Depth2, levels = c("Surface", "Middle", "Bottom"))
MC_12S_decon_p.melt$Depth2 <- factor(MC_12S_decon_p.melt$Depth2, levels = c("Surface", "Middle", "Bottom"))



######### 여기 아래부터가 그림@@@@######################
test <- read.xlsx(file = "Total input.xlsx", 8)
test$Station <- factor(test$Station, levels = c("St. 19", "St. 03", "St. 07", "St. 10", "St. 13"))
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
  geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(0,150)) +
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
  geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(0,260)) +
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
  geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10)) +
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
  geom_boxplot(aes(fill = Station), alpha = 0.2, color = "black") +
  scale_y_continuous(limits = c(0,75)) +
  ylab("") + 
  xlab("")

p2 <- plot_alpha_COI+plot_alpha_18S+plot_alpha_12S 

p2

p2/p1

ggsave("Abu_alpha diversity22.png", dpi=300, dev="png", height=12, width=12, units="in")


#ggsave("Abu_alpha diversity.png", dpi=300, dev="png", height=8, width=12, units="in")
#p1

Station_col <- c("St. 19" = "#F8766D", "St. 03"= "#ABA300", "St. 07" = "#0CB702", "St. 10" = "#00A9FF", "St. 13" = "#ED35ED") 
Depth_shp <- c("Surface" = 21, )
scale_fill_manual("", values = c("#F8766D", "#ABA300", "#0CB702", "#00A9FF", "#ED35ED"))


