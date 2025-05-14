options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:/tap/R/MC/CO1/97%/new/new/total/metaweb/new")
set.seed(123)

library("phyloseq")
library("tidyverse")
library("ggplot2")
library("xlsx")
library("graphlayouts")
library("igraph")
library("ggraph")
library("colorspace")
library("fluxweb")
library("dplyr")
library("brainGraph")
library("bipartite")
library("rlang")
library("patchwork")
library("webr")


nodes <- read.xlsx("metawb_new2.xlsx", 2)
links <- read.xlsx("metawb_new2.xlsx", 1)
mc_out <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
colFG <- c("salmon","saddlebrown","turquoise4","steelblue1","burlywood","maroon","violet","tan1","green","gray","green3","pink2" )
rownames(nodes) <- nodes[,1]
nodes$fg <- factor(nodes$fg)
levels(nodes$fg)
nodes$colfg <- colFG[as.numeric(nodes$fg)]
V(mc_out)$shape <- ifelse(V(mc_out)$FFG2 == "Solitary", "square", "circle")
plotfw(mc_out, col = nodes$colfg, size=nodes$Range*2, vertex.frame.color="black", vertex.frame.width = 2 , edge.width=0.3, edge.arrow.size=0.3, vertex.label=NA, vertex.shapes = mc_out$shape)
legend("topright", legend=levels(nodes$fg),  pch=19,col=colFG)


tlnodes <- trophiclevels(mc_out)

write.csv(file = "tlnodes_out.csv", tlnodes)


png("MC_out_foodweb_final.png", width=2800, height=2600, res=300)
plotfw(mc_out, col = nodes$colfg, size=nodes$Range*2, vertex.frame.color="black", vertex.frame.width = 2 , edge.width=0.3, edge.arrow.size=0.3, vertex.label=NA, vertex.shapes = mc_out$shape)
legend("topright", legend=levels(nodes$fg),  pch=19,col=colFG)
dev.off()

png("MC_out_foodweb_final_legendx.png", width=2800, height=2600, res=300)
plotfw(mc_out, col = nodes$colfg, size=nodes$Range*2, vertex.frame.color="black", vertex.frame.width = 2 , edge.width=0.3, edge.arrow.size=0.3, vertex.label=NA, vertex.shapes = mc_out$shape)

dev.off()


nodes2 <- read.xlsx("metawb_new2.xlsx", 4) 
links2 <- read.xlsx("metawb_new2.xlsx", 3)
mc_in <- graph_from_data_frame(d=links2, vertices=nodes2, directed=T)
colFG2 <- c("salmon","steelblue1","burlywood","maroon","violet","tan1","green","gray","green3","pink2")
rownames(nodes2) <- nodes2[,1]
nodes2$fg <- factor(nodes2$fg)
levels(nodes2$fg)
nodes2$colfg <- colFG2[as.numeric(nodes2$fg)]
V(mc_in)$shape <- ifelse(V(mc_in)$FFG2 == "Solitary", "square", "circle")
plotfw(mc_in, col = nodes2$colfg, size=nodes2$Range*2, vertex.frame.color="black", vertex.frame.width = 2 , edge.width=0.3, edge.arrow.size=0.3, vertex.label=NA, vertex.shapes = mc_in$shape)
legend("topright", legend=levels(nodes2$fg),  pch=19,col=colFG2)

tlnodes2 <- trophiclevels(mc_in)

write.csv(file = "tlnodes_in.csv", tlnodes2)

png("MC_in_foodweb_final.png", width=2800, height=2600, res=300)
plotfw(mc_in, col = nodes2$colfg, size=nodes2$Range*2, vertex.frame.color="black", vertex.frame.width = 2 , edge.width=0.3, edge.arrow.size=0.3, vertex.label=NA, vertex.shapes = mc_in$shape)
legend("topright", legend=levels(nodes2$fg),  pch=19,col=colFG2)
dev.off()


png("MC_in_foodweb_final_legendx.png", width=2800, height=2600, res=300)
plotfw(mc_in, col = nodes2$colfg, size=nodes2$Range*2, vertex.frame.color="black", vertex.frame.width = 2 , edge.width=0.3, edge.arrow.size=0.3, vertex.label=NA, vertex.shapes = mc_in$shape)

dev.off()

#### Simulated annealing algorithm ####
spinglass_communities_out <- cluster_spinglass(mc_out)
spinglass_communities_in <- cluster_spinglass(mc_in)

membership_in <- membership(spinglass_communities_in)
membership_out <- membership(spinglass_communities_out)
membership_in_list <- as.list(membership_in)
membership_out_list <- as.list(membership_out)
write.csv(file = "membership_in.csv", membership_in_list)
write.csv(file = "membership_out.csv", membership_out_list)


#  total link (L) 
L_out <- ecount(mc_out)
L_in <- ecount(mc_in)
# Number of module (N_M)
N_M_out <- length(unique(membership_out))
N_M_in <- length(unique(membership_in))

# igraph modularity 
modularity(mc_out,membership = membership_out) #0.3747972
modularity(mc_in,membership = membership_in) #0.3067098


# Plot network 
#plot(mc_out, vertex.color = membership(spinglass_communities_out), vertex.label=NA, main = "Community Detection with Spin-glass Model")
#plot(mc_in, vertex.color = membership(spinglass_communities_in), vertex.label=NA, main = "Community Detection with Spin-glass Model")


#Sp Topological roles
Topo_role_in <- read.xlsx(file = "PC2.xlsx", 1)
Topo_role_out <- read.xlsx(file = "PC2.xlsx", 2)


Out_fill <- c("Apex predator" = "salmon", "Omnivorous fish"= "violet", "Omnivorous scavenger" = "tan1", "Invertebrate predator" = "maroon", "Filter feeder" = "steelblue1", "Grazer" = "burlywood", "Sessile filterer" = "pink2", "Decomposer" = "saddlebrown", "Deposit feeder" = "turquoise4", "Phytoplankton" = "green", "Resource" = "gray", "Seaweed" = "green3") #2

In_fill <- c("Apex predator" = "salmon",  "Omnivorous fish"= "violet", "Omnivorous scavenger" = "tan1", "Invertebrate predator" = "maroon", "Filter feeder" = "steelblue1", "Grazer" = "burlywood", "Sessile filterer" = "pink2", "Phytoplankton" = "green", "Resource" = "gray", "Seaweed" = "green3") #2

plot_topo_in <- ggplot()  +
  geom_point(data = Topo_role_in , aes(x = PC, y = Zscore, fill = FFG), 
             pch = 21, color = "black", size = 6, stroke = 1) +
  labs(x = "Among module connectivity (PC)", y = expression(paste(bold("within module degree" ~(italic("z")))))) + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15)) +
  scale_y_continuous(limits = c(-2,4), breaks = c(-2,-1,0,1,2,3,4)) +
  scale_x_continuous(limits = c(0,0.81), breaks = c(0,0.2,0.4,0.6,0.8)) +
  geom_vline(xintercept = c(0.625), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(2.5), color = "grey70", linetype = 2) + 
  scale_fill_manual(name = "Functional feeding group" ,values = In_fill)
plot_topo_in
ggsave("Z-score, PC _ in.png", dpi=300, dev="png", height=8, width=12, units="in")

plot_topo_out <- ggplot()  +
  geom_point(data = Topo_role_out , aes(x = PC, y = Zscore, fill = FFG), 
             pch = 21, color = "black", size = 6, stroke = 1) +
  labs(x = "Among module connectivity (PC)", y = expression(paste(bold("within module degree" ~(italic("z")))))) + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15)) +
  scale_y_continuous(limits = c(-2,3), breaks = c(-2,-1,0,1,2,3)) +
  scale_x_continuous(limits = c(0,0.81), breaks = c(0,0.2,0.4,0.6,0.8)) +
  geom_vline(xintercept = c(0.625), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(2.5), color = "grey70", linetype = 2) + 
  scale_fill_manual(name = "Functional feeding group", values = Out_fill)

plot_topo_out
ggsave("Z-score, PC _ out.png", dpi=300, dev="png", height=8, width=12, units="in")

legend_size <- c(1,3,5,7,9) 
plot_topo_in <- ggplot()  +
  geom_point(data = Topo_role_in , aes(x = PC, y = Zscore, fill = FFG, size = Index), 
             pch = 21, color = "black", stroke = 1) +
  labs(x = "Among module connectivity (PC)", y = expression(paste(bold("within module degree" ~(italic("z")))))) + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15)) +
  scale_y_continuous(limits = c(-2,4), breaks = c(-2,-1,0,1,2,3,4)) +
  scale_x_continuous(limits = c(0,0.81), breaks = c(0,0.2,0.4,0.6,0.8)) +
  geom_vline(xintercept = c(0.625), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(2.5), color = "grey70", linetype = 2) + 
  scale_fill_manual(name = "Functional feeding group" ,values = In_fill) +
  scale_size(breaks = c(1:5), range = c(2,10)) + 
  guides(fill = guide_legend(), 
         size = guide_legend(override.aes = list(size = legend_size)))


plot_topo_in 

ggsave("Z-score, PC _ in(size).png", dpi=300, dev="png", height=8, width=12, units="in")


plot_topo_out <- ggplot()  +
  geom_point(data = Topo_role_out , aes(x = PC, y = Zscore, fill = FFG, size = Index), 
             pch = 21, color = "black", stroke = 1) +
  labs(x = "Among module connectivity (PC)", y = expression(paste(bold("within module degree" ~(italic("z")))))) + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15)) +
  scale_y_continuous(limits = c(-2,3), breaks = c(-2,-1,0,1,2,3)) +
  scale_x_continuous(limits = c(0,0.81), breaks = c(0,0.2,0.4,0.6,0.8)) +
  geom_vline(xintercept = c(0.625), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(2.5), color = "grey70", linetype = 2) + 
  scale_fill_manual(name = "Functional feeding group", values = Out_fill) +
  scale_size(breaks = c(1:5), range = c(2,10)) + 
  guides(fill = guide_legend(), 
         size = guide_legend(override.aes = list(size = legend_size)))

plot_topo_out

ggsave("Z-score, PC _ out(size).png", dpi=300, dev="png", height=8, width=12, units="in")

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 1000 randomized network Zscore and PC @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  mc_out3 <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
num_nodes_out <- vcount(mc_out3)
num_randomizations <- 1000


modularity_values_out <- numeric(num_randomizations)
z_scores_matrix_out <- matrix(0, nrow = num_nodes_out, ncol = num_randomizations)
pc_matrix_out <- matrix(0, nrow = num_nodes_out, ncol = num_randomizations)

# 1000 times Z,PC
for (i in 1:num_randomizations) {
  community_out <- cluster_spinglass(mc_out3)
  membership_out <- membership(community_out)
  modularity_values_out[i] <- modularity(community_out)
  z_scores_matrix_out[, i] <- within_module_deg_z_score(mc_out3, membership_out)
  pc_matrix_out[, i] <- part_coeff(mc_out3, membership_out)
  
  
}

modularity_mean_out <- mean(modularity_values_out) #0.3240228
modularity_sd_out <- sd(modularity_values_out) #0.0008002786



modularity_empirical_out <- 0.375
# Calculate Z-score
z_score_out <- (modularity_empirical_out - modularity_mean_out) / modularity_sd_out

# Calculate p-value based on z-score
p_value_out<- 2 * pnorm(-abs(z_score_out))

p_value_out

write.csv(z_scores_matrix_out, "Z_scores_matrix_out.csv", row.names = FALSE)
write.csv(pc_matrix_out, "PC_matrix_out.csv", row.names = FALSE)

Topo_role_out1000 <- read.xlsx(file = "PC1000.xlsx", 1)


plot_topo_out1000 <- ggplot()  +
  geom_point(data = Topo_role_out1000 , aes(x = PC, y = Zscore, fill = FFG, size = Index), 
             pch = 21, color = "black", stroke = 1) +
  geom_errorbar(data = Topo_role_out1000 , aes(x = PC, y = Zscore, xmin = PC-PC_sd, xmax = PC+PC_sd), width = 0.05, color = "black") +
  geom_errorbar(data = Topo_role_out1000 , aes(x = PC, y = Zscore, ymin = Zscore-Z_sd, ymax = Zscore+Z_sd), width = 0.01, color = "black") +
  labs(x = "Among module connectivity (PC)", y = expression(paste(bold("Within module degree" ~(italic("z")))))) + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15)) +
  scale_y_continuous(limits = c(-2.3,3), breaks = c(-2,-1,0,1,2,3)) +
  scale_x_continuous(limits = c(-0.1,0.81), breaks = c(0,0.2,0.4,0.6,0.8)) +
  geom_vline(xintercept = c(0.625), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(2.5), color = "grey70", linetype = 2) + 
  scale_fill_manual(name = "Functional feeding group", values = Out_fill) +
  scale_size(breaks = c(1:5), range = c(2,10)) + 
  guides(fill = guide_legend(), 
         size = guide_legend(override.aes = list(size = legend_size)))

plot_topo_out1000

ggsave("Z-score, PC _ out(size) 1000 randomized_v1.png", dpi=300, dev="png", height=8, width=12, units="in")


### Inner @@@

mc_in3 <- graph_from_data_frame(d=links2, vertices=nodes2, directed=F)

num_nodes_in <- vcount(mc_in3)
num_randomizations <- 1000

modularity_values_in <- numeric(num_randomizations)
z_scores_matrix_in <- matrix(0, nrow = num_nodes_in, ncol = num_randomizations)
pc_matrix_in <- matrix(0, nrow = num_nodes_in, ncol = num_randomizations)

for (i in 1:num_randomizations) {
  community_in <- cluster_spinglass(mc_in3)
  membership_in <- membership(community_in)
  modularity_values_in[i] <- modularity(community_in)
  z_scores_matrix_in[, i] <- within_module_deg_z_score(mc_in3, membership_in)
  pc_matrix_in[, i] <- part_coeff(mc_in3, membership_in)
  
}
mean(modularity_values_in) # 0.2839162
sd(modularity_values_in) #0.0007865634
write.csv(z_scores_matrix_in, "Z_scores_matrix_in.csv", row.names = FALSE)
write.csv(pc_matrix_in, "PC_matrix_in.csv", row.names = FALSE)


modularity_empirical_in <- 0.307
modularity_mean_in <- mean(modularity_values_in) # 0.2839162
modularity_sd_in <- sd(modularity_values_in) #0.0007865634


# Calculate Z-score
z_score_in <- (modularity_empirical_in - modularity_mean_in) / modularity_sd_in

# Calculate p-value based on z-score
p_value_in <- 2 * pnorm(-abs(z_score_in))

p_value_in


Topo_role_in1000 <- read.xlsx(file = "PC1000.xlsx", 2)


plot_topo_in1000 <- ggplot()  +
  geom_point(data = Topo_role_in1000 , aes(x = PC, y = Zscore, fill = FFG, size = Index), 
             pch = 21, color = "black", stroke = 1) +
  geom_errorbar(data = Topo_role_in1000 , aes(x = PC, y = Zscore, xmin = PC-PC_sd, xmax = PC+PC_sd), width = 0.05, color = "black") +
  geom_errorbar(data = Topo_role_in1000 , aes(x = PC, y = Zscore, ymin = Zscore-Z_sd, ymax = Zscore+Z_sd), width = 0.01, color = "black") +
  labs(x = "Among module connectivity (PC)", y = expression(paste(bold("Within module degree" ~(italic("z")))))) + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15)) +
  scale_y_continuous(limits = c(-2.3,3.5), breaks = c(-2,-1,0,1,2,3)) +
  scale_x_continuous(limits = c(-0.1,0.81), breaks = c(0,0.2,0.4,0.6,0.8)) +
  geom_vline(xintercept = c(0.625), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(2.5), color = "grey70", linetype = 2) + 
  scale_fill_manual(name = "Functional feeding group", values = Out_fill) +
  scale_size(breaks = c(1:5), range = c(2,10)) + 
  guides(fill = guide_legend(), 
         size = guide_legend(override.aes = list(size = legend_size)))

plot_topo_in1000

ggsave("Z-score, PC _ in (size) 1000 randomized_v1.png", dpi=300, dev="png", height=8, width=12, units="in")


#Indegree ###
degree_data <- read.xlsx("Degree_new.xlsx", 1)

# In-degree & Cove data
degree_values <- as.numeric(degree_data$Indegree)
cove_values <- degree_data$Cove

# Divide Inner, Outer 
inner_data <- degree_values[cove_values == "Inner"]
outer_data <- degree_values[cove_values == "Outer"]

# Cumulative distribution
calculate_cumulative <- function(degree) {
  degree <- sort(degree, decreasing = TRUE)
  cumulative <- cumsum(degree)
  cumulative_distribution <- cumulative / sum(degree)
  return(cumulative_distribution)
}

cumulative_inner <- calculate_cumulative(inner_data)
cumulative_outer <- calculate_cumulative(outer_data)

# data frame 
df_inner <- data.frame(Degree = sort(inner_data, decreasing = TRUE), Cumulative = cumulative_inner, Cove = "Inner")
df_outer <- data.frame(Degree = sort(outer_data, decreasing = TRUE), Cumulative = cumulative_outer, Cove = "Outer")

# Combine
df_combined_indegree <- rbind(df_inner, df_outer)


# Plot
p <- ggplot()  +
  geom_point(data = df_combined_indegree, aes(x = Degree, y = Cumulative, fill = Cove), 
             pch = 21, colour = "black", size = 6, stroke = 1) +
  labs(x = "In-degree", y = "Log cumulative distribution") + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15), 
        legend.position = c(0.9, 0.9), legend.margin = unit(-0.6,"mm"))  +
  scale_fill_manual(values = c("Inner" = "#D62728", "Outer" = "#1F77B4")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))

p
ggsave("In-degree_Log Cumltatuve distribution_new.png", dpi=300, dev="png", height=8, width=12, units="in")

#Outdegree###
degree_data <- read.xlsx("Degree_new.xlsx", 2)

degree_values <- as.numeric(degree_data$Outdegree)
cove_values <- degree_data$Cove

inner_data <- degree_values[cove_values == "Inner"]
outer_data <- degree_values[cove_values == "Outer"]

calculate_cumulative <- function(degree) {
  degree <- sort(degree, decreasing = TRUE)
  cumulative <- cumsum(degree)
  cumulative_distribution <- cumulative / sum(degree)
  return(cumulative_distribution)
}

cumulative_inner <- calculate_cumulative(inner_data)
cumulative_outer <- calculate_cumulative(outer_data)


df_inner <- data.frame(Degree = sort(inner_data, decreasing = TRUE), Cumulative = cumulative_inner, Cove = "Inner")
df_outer <- data.frame(Degree = sort(outer_data, decreasing = TRUE), Cumulative = cumulative_outer, Cove = "Outer")

df_combined_outdegree <- rbind(df_inner, df_outer)


p1 <- ggplot()  +
  geom_point(data = df_combined_outdegree , aes(x = Degree, y = Cumulative, fill = Cove), 
             pch = 21, color = "black", size = 6, stroke = 1) +
  labs(x = "Out-degree", y = "Log cumulative distribution") + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15), 
        legend.position = c(0.9, 0.9), legend.margin = unit(-0.6,"mm"))  +
  scale_fill_manual(values = c("Inner" = "#D62728", "Outer" = "#1F77B4")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))
p1

ggsave("Out-degree_Log Cumltatuve distribution_new.png", dpi=300, dev="png", height=8, width=12, units="in")


###################
####Total degree###
###################
degree_total <- read.xlsx("Degree_new.xlsx", 3)

total_degree_values <- as.numeric(degree_total$Degree)
cove_values <- degree_total$Cove

inner_data <- total_degree_values[cove_values == "Inner"]
outer_data <- total_degree_values[cove_values == "Outer"]


calculate_cumulative <- function(degree) {
  degree <- sort(degree, decreasing = TRUE)
  cumulative <- cumsum(degree)
  cumulative_distribution <- cumulative / sum(degree)
  return(cumulative_distribution)
}

cumulative_inner <- calculate_cumulative(inner_data)
cumulative_outer <- calculate_cumulative(outer_data)

df_inner <- data.frame(Degree = sort(inner_data, decreasing = TRUE), Cumulative = cumulative_inner, Cove = "Inner")
df_outer <- data.frame(Degree = sort(outer_data, decreasing = TRUE), Cumulative = cumulative_outer, Cove = "Outer")

df_combined_outdegree <- rbind(df_inner, df_outer)


p2 <- ggplot()  +
  geom_point(data = df_combined_outdegree , aes(x = Degree, y = Cumulative, fill = Cove), 
             pch = 21, color = "black", size = 6, stroke = 1) +
  labs(x = "Degree", y = "Log cumulative distribution") + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15), 
        legend.position = c(0.9, 0.9), legend.margin = unit(-0.6,"mm"))  +
  scale_fill_manual(values = c("Inner" = "#D62728", "Outer" = "#1F77B4")) +
  scale_y_continuous(limits = c(0,1)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))
p2

ggsave("Total_degree_Log Cumltatuve distribution_new.png", dpi=300, dev="png", height=8, width=12, units="in")


#Topological role of species


FFG_Topo_role_in1000 <- read.xlsx(file = "PC1000.xlsx", 4)


plot_FFG_topo_in1000 <- ggplot()  +
  geom_point(data = FFG_Topo_role_in1000 , aes(x = PC, y = Zscore, fill = FFG), 
             pch = 21, color = "black", stroke = 1, size = 7 ) +
  geom_errorbar(data = FFG_Topo_role_in1000 , aes(x = PC, y = Zscore, xmin = PC-PC_sd, xmax = PC+PC_sd), width = 0.05, color = "black") +
  geom_errorbar(data = FFG_Topo_role_in1000 , aes(x = PC, y = Zscore, ymin = Zscore-Z_sd, ymax = Zscore+Z_sd), width = 0.01, color = "black") +
  labs(x = "Among module connectivity (PC)", y = expression(paste(bold("within module degree" ~(italic("z")))))) + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15)) +
  scale_y_continuous(limits = c(-2.3,3.5), breaks = c(-2,-1,0,1,2,3)) +
  scale_x_continuous(limits = c(-0.1,0.81), breaks = c(0,0.2,0.4,0.6,0.8)) +
  geom_vline(xintercept = c(0.625), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(2.5), color = "grey70", linetype = 2) + 
  scale_fill_manual(name = "Functional feeding group", values = Out_fill) +
  scale_size(breaks = c(1:5), range = c(2,10)) + 
  guides(fill = guide_legend(), 
         size = guide_legend(override.aes = list(size = legend_size)))

plot_FFG_topo_in1000

ggsave("FFG-Z-score, PC _ in (size) 1000 randomized_v1.png", dpi=300, dev="png", height=8, width=12, units="in")


FFG_Topo_role_out1000 <- read.xlsx(file = "PC1000.xlsx", 3)


plot_FFG_topo_out1000 <- ggplot()  +
  geom_point(data = FFG_Topo_role_out1000 , aes(x = PC, y = Zscore, fill = FFG), 
             pch = 21, color = "black", stroke = 1 , size = 7) +
  geom_errorbar(data = FFG_Topo_role_out1000 , aes(x = PC, y = Zscore, xmin = PC-PC_sd, xmax = PC+PC_sd), width = 0.05, color = "black") +
  geom_errorbar(data = FFG_Topo_role_out1000 , aes(x = PC, y = Zscore, ymin = Zscore-Z_sd, ymax = Zscore+Z_sd), width = 0.01, color = "black") +
  labs(x = "Among module connectivity (PC)", y = expression(paste(bold("within module degree" ~(italic("z")))))) + 
  theme_bw() + 
  theme(axis.title = element_text(face = "bold", family = "Tahoma", size = 25),
        axis.title.x = element_text(vjust = -0.15), 
        axis.title.y = element_text(vjust = 2),
        axis.text = element_text(family = "Tahoma", size= 18),
        panel.grid = element_blank() ,
        legend.title = element_text(family = "Tahoma", face = "bold", size = 18), 
        legend.text = element_text(family = "Tahoma", size = 15)) +
  scale_y_continuous(limits = c(-2.3,3), breaks = c(-2,-1,0,1,2,3)) +
  scale_x_continuous(limits = c(-0.1,0.81), breaks = c(0,0.2,0.4,0.6,0.8)) +
  geom_vline(xintercept = c(0.625), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(2.5), color = "grey70", linetype = 2) + 
  scale_fill_manual(name = "Functional feeding group", values = Out_fill) +
  scale_size(breaks = c(1:5), range = c(2,10)) + 
  guides(fill = guide_legend(), 
         size = guide_legend(override.aes = list(size = legend_size)))

plot_FFG_topo_out1000

ggsave("FFG-Z-score, PC _ out (size) 1000 randomized_v1.png", dpi=300, dev="png", height=8, width=12, units="in")








