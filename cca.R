options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:/tap/R/MC/CO1/97%/new/new/total/PCA")
set.seed(123)

library("phyloseq")
library("dplyr")
library("xlsx")
library("extrafont")
library("ggplot2")
library("ggfortify")
library("vegan")
library("ggrepel")
library("ggord")
library("vegan")

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

MC_total_melt <- psmelt(MC_total)
write.csv(file = "MC_total_decontam.csv", MC_total_melt)

#####
abu <- read.xlsx(file = "PCA_LDA.xlsx",2)
sm <- read.xlsx(file = "PCA_LDA.xlsx",3)

Sp <- abu[,2:ncol(abu)]
Env <- sm[,2:10]

cca_result <- cca(Sp, Env)
plot(cca_result) 
veg_1 = as.data.frame(cca_result$CCA$biplot)
veg_1$env <- c("Temperature", "Salinity", "Fucoxanthin", "Chl-a", "Phosphate", "Nitrite+Nitrate", "Silicate", "Ammonium", "SPM")
veg_1["env"] = row.names(veg_1)
veg_2 = as.data.frame(cca_result$CCA$v)
veg_2["genus"] = row.names(veg_2)

write.xlsx(file = "CCA_ggplot_input.xlsx", veg_2)
write.xlsx(file = "CCA_ggplot_input2.xlsx", veg_1)
cca_input <- read.xlsx(file = "CCA_ggplot_input.xlsx", 2)
cca_input2 <- read.xlsx(file = "CCA_ggplot_input.xlsx", 3)
plot = ggplot() + 
  geom_point(data = cca_input, aes(x = CCA1, y = CCA2, color = Cove), size = 3) + 
  geom_point(data = cca_input2, aes(x = CCA1, y = CCA2), size = 0) +
  scale_color_manual(values = c("#D62728","#1F77B4")) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) 

#  stat_ellipse(data = cca_input, aes(color = Cove))


plot + geom_text_repel(data = veg_2,aes(x = CCA1, y = CCA2, label = veg_2$genus),nudge_y = -0.05, size = 5) +
  theme_bw() +
  geom_segment(data = cca_input2,aes(x = 0,y = 0,xend = CCA1,yend = CCA2),color ="blue" , arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel(data = cca_input2,aes(x = CCA1, y = CCA2, label = Sample),nudge_y = -0.1,color = "blue",size = 5) +
  theme(axis.text = element_text( family = "Tahoma", size = 18),
        axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = -0.5),
        legend.title = element_text(size = 20, family = "Tahoma",  face = "bold"), 
        legend.text = element_text(size = 18, family = "Tahoma"),
        panel.grid = element_blank())

ggsave("CCA_LDA genus .png", dpi=300, dev="png", height=8, width=12, units="in")


###이거 공통된 샘플의 개수가 부족해서 불가능####
st <- read.xlsx(file = "PCA_LDA.xlsx", 6)
sm2 <- read.xlsx(file = "PCA_LDA.xlsx",7)
St <- st[,2:ncol(st)]
Env2 <- sm2[,2:10] 

cca_result_st <- cca(St, Env2)
plot(cca_result_st) #결과 확인
veg_3 = as.data.frame(cca_result$CCA$biplot)
veg_4["env"] = row.names(veg_3)
veg_3 = as.data.frame(cca_result$CCA$v)
veg_4["Station"] = row.names(veg_4)
#####



cc_result <- cancor(Sp,Env)
cc_result$xcoef
cc_result$ycoef

cc_result$cor

cc1_X <- as.matrix(X) %*% cc_result$xcoef[,1]
cc1_Y <- as.matrix(Y) %*% cc_result$ycoef[,1]

cc2_X <- as.matrix(X) %*% cc_result$xcoef[,2]
cc2_Y <- as.matrix(Y) %*% cc_result$ycoef[,2]

cor(cc1_X, cc1_Y)

assertthat::are_equal(cc_result$cor[1],
                      cor(cc1_X,cc1_Y)[1])

cca_df <- ggplot(aes(x=cc1_X, y=cc1_Y)) +
  geom_point()





#RDA
ord <- rda(Sp, Env)
ggord(ord)
ggord(ord, ptslab =TRUE, size = -1, addsize = 3, addcol = "black",veccol = "blue", labcol = "blue", parse = TRUE, ext=1.05, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))

#NMDS
m_Sp <- as.matrix(Sp)
nmds = metaMDS(m_Sp, distance = "bray")
en = envfit(nmds, Env, permutations = 999, na.rm =TRUE)

write.csv(file = "genus_NMDS.csv", nmds$species)
genus_nmds <- read.xlsx(file = "CCA_ggplot_input.xlsx", 4)
#data.scores_nmds = as.data.frame(scores(nmds)$sites)
#data.scores_nmds$Cove2 = sm$Cove2

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)
write.csv(file = "Sample_NMDS.csv", en_coord_cont)
sample_nmds <- read.xlsx(file = "CCA_ggplot_input.xlsx", 5)

gg = ggplot() + 
  geom_point(data = genus_nmds, aes(x = NMDS1, y = NMDS2, color = Cove), size = 3) + 
  scale_color_manual(values = c("#D62728","#1F77B4")) + 
  geom_text_repel(data = genus_nmds, aes(x = NMDS1, y = NMDS2, label = genus_nmds$Genus), size = 5) +
  theme_bw() +
  theme(axis.text = element_text( family = "Tahoma", size = 18),
        axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = -0.5),
        legend.title = element_text(size = 20, family = "Tahoma",  face = "bold"), 
        legend.text = element_text(size = 18, family = "Tahoma"),
        panel.grid = element_blank()) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) 

gg
gg  + geom_text_repel(data = sample_nmds,aes(x = NMDS1, y = NMDS2, label = sample_nmds$Sample), nudge_y = 0.01, color = "blue", size = 5) +
  geom_segment(data = sample_nmds, aes(x = 0,y = 0,xend = NMDS1, yend = NMDS2),color ="blue" , arrow = arrow(length = unit(0.25, "cm"))) 


ggsave("NMDS_LDA genus .png", dpi=300, dev="png", height=8, width=12, units="in")


