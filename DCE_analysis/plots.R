######BEFORE YOU START PLOTTING!!! 
###Create a folder named "svg_plots" in your working directory so that all the plots can be saved in .svg format in the same folder
dir.create("svg_plots")

library(ggplot2) #data visualization
library(ggpubr) #data visualization
library(gginnards) #data visualization
library(ggsci) #data visualization
library(plyr) #data and plot wrangling
library(pheatmap) #create heatmaps
library(tidyverse) #data wrangling
library(karyoploteR) #create chromosome plots
library(TxDb.Hsapiens.UCSC.hg38.knownGene) #annotation data
library(readxl) #import excel files (.xlsx)

#create a vector named "my_comparisons" to compare the data of statistical significance
my_comparisons <- list(c("NP", "EBV"), c("NP", "CIN"), c("NP","MSI"), c("NP","GS"))

#do the same with the levels for factorizing

my_levels <- c("NP", "EBV", "CIN", "MSI", "GS")

size <- 20 #fontsize value

########Average Binsignal########

average_binsignal <- 
  apply(X = sapply(bin_signals, function(x) sapply(x, function(y) mean(y[,1]))), 
        MARGIN = 2, 
        FUN = mean)

average_binsignal %<>% 
  as.data.frame %>% 
  rownames_to_column('group')
colnames(average_binsignal)[2] <- "Average_binsignal"
average_binsignal$group <- factor(average_binsignal$group, levels = my_levels)


svg("svg_plots/average_binsignal.svg", width=10, height=8)
ggplot(average_binsignal, aes(x = group, y = Average_binsignal, fill = group)) +
  ggtitle("Average Binsignal for each GC subtype")+
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = pal_uchicago()(8)[c(3,2,1,5,7,8)]) +
  ylab(label = "Average binsignal") +
  theme(text = element_text(size = size), 
        axis.title = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 30, hjust = 0.5),
        legend.position = "none")+
  ylab("Average Binsignal")
dev.off()

######Average binsignal with significance#####
bsignals <- data.frame(
average_binsignal,
  sd = c(sd(unlist(bin_signals[["NP"]])),
         sd(unlist(bin_signals[["EBV"]])),
         sd(unlist(bin_signals[["CIN"]])),
         sd(unlist(bin_signals[["MSI"]])),
         sd(unlist(bin_signals[["GS"]])))
)

bsignals$group <- factor(bsignals$group, levels = my_levels)

svg("svg_plots/average_binsignal_signif.svg", width=10, height=8)
ggplot(bsignals, aes(x=group, y=Average_binsignal, fill=group)) + 
  ggtitle("Average Binsignal for each GC subtype")+
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values = pal_uchicago()(8)[c(3,2,1,5,7,8)]) +
  geom_errorbar(aes(ymin=Average_binsignal-sd, ymax=Average_binsignal+sd), width=.2,
                position=position_dodge(.9)) +
  stat_compare_means(method = "kruskal.test", label.y = 0.43, size = 5, label.x = 3.5) + 
  stat_compare_means(label = "p.signif",comparisons = my_comparisons) +
  theme(text = element_text(size = size), 
        axis.title = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 30, hjust = 0.5),
        legend.position = "none")+
  ylab("Average Binsignal")
dev.off()


rm(bsignals)


#####Number of bins with non-zero expression (all of them, not just in DCE)


no_nz_bins_NP <- t(rbind.data.frame(no_nz_bins[["NP"]])) 
colnames(no_nz_bins_NP) <- "NP"

no_nz_bins_CIN <- t(rbind.data.frame(no_nz_bins[["CIN"]])) 
colnames(no_nz_bins_CIN) <- "CIN"

no_nz_bins_EBV <- t(rbind.data.frame(no_nz_bins[["EBV"]])) 
colnames(no_nz_bins_EBV) <- "EBV"

no_nz_bins_MSI <- t(rbind.data.frame(no_nz_bins[["MSI"]])) 
colnames(no_nz_bins_MSI) <- "MSI"

no_nz_bins_GS <- t(rbind.data.frame(no_nz_bins[["GS"]])) 
colnames(no_nz_bins_GS) <- "GS"

#all together
no_nz_bins_merged <- cbind.data.frame(
  sum(no_nz_bins_NP),
  sum(no_nz_bins_CIN),
  sum(no_nz_bins_EBV),
  sum(no_nz_bins_MSI),
  sum(no_nz_bins_GS)
)

colnames(no_nz_bins_merged) <- c("NP", "CIN", "EBV", "MSI", "GS")
row.names(no_nz_bins_merged) <- "sum of nz bins"

no_nz_bins_merged <- rbind.data.frame(
  no_nz_bins_merged,
  c("NP", "CIN", "EBV", "MSI", "GS")
)

row.names(no_nz_bins_merged) <- c("sum", "subtype")
no_nz_bins_merged <- t(no_nz_bins_merged)
no_nz_bins_merged <- as.data.frame(no_nz_bins_merged)
no_nz_bins_merged[,1] <- as.numeric(no_nz_bins_merged[,1])
no_nz_bins_merged$subtype <- factor(no_nz_bins_merged$subtype, levels = my_levels)

#barplot
svg("svg_plots/number_of_nz_bins.svg", width=10, height=8)
ggplot(no_nz_bins_merged, aes(x = subtype, y = sum, fill = subtype, color = "black")) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = pal_uchicago()(8)[c(3,2,1,5,7,8)]) +
  theme(legend.position = "NONE",
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  labs(title="Number of non zero bins in each GC subtype", 
       y = "Number of bins")
dev.off()

########Number of bins found inside a DCE########
nob <- function(cods){
  n_bins <- list()
  for (group in names(cods)){
    for (chr in names(cods[[names(cods)[1]]])){
      x <- cods[[group]][[chr]]
      n <- x[,"END"]-x[,"START"]+1
      n_bins[[group]][[chr]] <- unname(n)
    }
  }
  return(n_bins)
}
n_bins <- nob(CODs)
sapply(n_bins, function(x) summary(unlist(x , use.names = F)))


n_bins <- lapply(n_bins, unlist, use.names = F)
n_bins <- data.frame(no_bins =c(n_bins$NP, n_bins$CIN, n_bins$EBV, n_bins$MSI, n_bins$GS), 
                     group = c( rep("NP", length(n_bins$NP)), 
                                rep("CIN", length(n_bins$CIN)), 
                                rep("EBV", length(n_bins$EBV)), 
                                rep("MSI", length(n_bins$MSI)), 
                                rep("GS", length(n_bins$GS))))

n_bins$group <- factor(n_bins$group, levels = my_levels)
perc1 <- data.frame(percentage = as.numeric(perc), group = rep(c(rownames(perc)), ncol(perc))) #put them in the order of the "perc" matrix rownames!
perc1$group <- factor(perc1$group, levels = my_levels)

svg("svg_plots/perc_of_bins_in_dces.svg", width=10, height=8)
ggplot(perc1, aes(x=group, y=percentage, fill=group)) +
  geom_violin(trim=FALSE)+ 
  scale_fill_manual(values = pal_uchicago()(8)[c(3,2,1,5,7,8)]) +
  geom_hline(yintercept = mean(perc1$percentage), linetype = 2) +
  geom_boxplot(width=0.1) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = 60, label.x = 4, size = 5) +
  theme(legend.position = "NONE",
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  labs(title="Percentage of bins found inside a DCE", 
       y = "Percentage (%)")
dev.off()

#######DCE sizes##########
s_CODs1$group <- factor(s_CODs1$group, levels = my_levels)


# ggplot(s_CODs1, aes(x=group, y=size, fill=group)) +
#   geom_violin(trim=FALSE)+
#   geom_boxplot(width=0.1, outlier.shape = NA) +
#   geom_hline(yintercept = mean(s_CODs1$size), linetype = 2) +
#   stat_compare_means(ref.group = "NP", label = "p.signif", label.y = 1.5) +
#   stat_compare_means(method = "kruskal.test", label.y = 1.25, label.x = 1.25, size = 5)+
#   #scale_y_continuous(limits = c(-0.2, 1.5)) +
#   theme(legend.position = "NONE",
#         plot.title = element_text(size = 25, hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 20),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15)) +
#   scale_fill_manual(values = pal_uchicago()(8)[c(3,2,1,5,7,8)]) +
#   labs(title="DCE sizes of GC subtypes", 
#        y = "DCE size (log10(bps))")
  

p <- ggplot(s_CODs1, aes(x=group, y=size, fill=group)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = mean(s_CODs1$size), linetype = 2) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  stat_compare_means(method = "kruskal.test", size = 5, label.y = 7.5, label.x = 4)+
  #scale_y_continuous(limits = c(-0.2, 1.5)) +
  theme(legend.position = "NONE",
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_fill_manual(values = pal_uchicago()(8)[c(3,2,1,5,7,8)]) +
  labs(title="DCE sizes of GC subtypes", 
       y = "DCE size (log10(bps))")


#p[["layers"]][[4]][["geom"]][["default_aes"]][["size"]] <- 6

svg("svg_plots/dce_sizes.svg", width=10, height=8)
p
dev.off()


######Number of DCEs in each subtype###########
sum_cods <- cbind.data.frame(
  length(unlist(s_CODs[["NP"]])),
  length(unlist(s_CODs[["CIN"]])),
  length(unlist(s_CODs[["EBV"]])),
  length(unlist(s_CODs[["MSI"]])),
  length(unlist(s_CODs[["GS"]]))
)

colnames(sum_cods) <- c("NP", "CIN", "EBV", "MSI", "GS")
row.names(sum_cods) <- "sum of dce sizes"

sum_cods <- rbind.data.frame(
  sum_cods,
  c("NP", "CIN", "EBV", "MSI", "GS")
)

row.names(sum_cods) <- c("sum", "subtype")
sum_cods <- t(sum_cods)
sum_cods <- as.data.frame(sum_cods)
sum_cods[,1] <- as.numeric(sum_cods[,1])
sum_cods$subtype <- factor(sum_cods$subtype, levels = my_levels)


#barplot
svg("svg_plots/number_of_dces.svg", width=10, height=8)
ggplot(sum_cods, aes(x = subtype, y = sum, fill = subtype, color = "black")) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = pal_uchicago()(8)[c(3,2,1,5,7,8)]) +
  theme(legend.position = "NONE",
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  labs(title="Number of DCEs in each gastric cancer subtype",
       x ="Subtypes", 
       y = "Number of DCEs")
dev.off()

###########Domainograms##########
#load the chromosomes dataframe (provided in the github repo as an xlsx file)
chromosomes_hg38p14 <- read_excel("chromosomes.hg38p14.xlsx", 
                                  col_names = FALSE, 
                                  col_types = c("numeric","text", "numeric", "numeric"))
View(chromosomes_hg38p14)


DCE_NP <- read.table("DCEs_NP", sep = '\t', header = T)
DCE_CIN <- read.table("DCEs_CIN", sep = '\t', header = T)
DCE_EBV <- read.table("DCEs_EBV", sep = '\t', header = T)
DCE_MSI <- read.table("DCEs_MSI", sep = '\t', header = T)
DCE_GS <- read.table("DCEs_GS", sep = '\t', header = T)

source("plot_domainograms.R")
#par(mfrow = c(1,1))


svg("svg_plots/domainogram_NP.svg", height = 8, width = 10)
p2 <- human.genome.plot(chromosomes= chromosomes_hg38p14, data= DCE_NP, dataname= "NP")
dev.off()

svg("svg_plots/domainogram_CIN.svg", height = 8, width = 10)
p3 <- human.genome.plot(chromosomes= chromosomes_hg38p14, data= DCE_CIN, dataname= "CIN")
dev.off()

svg("svg_plots/domainogram_EBV.svg", height = 8, width = 10)
p4 <- human.genome.plot(chromosomes= chromosomes_hg38p14, data= DCE_EBV, dataname= "EBV")
dev.off()

svg("svg_plots/domainogram_MSI.svg", height = 8, width = 10)
p5 <- human.genome.plot(chromosomes= chromosomes_hg38p14, data= DCE_MSI, dataname= "MSI")
dev.off()

svg("svg_plots/domainogram_GS.svg", height = 8, width = 10)
p6 <- human.genome.plot(chromosomes= chromosomes_hg38p14, data= DCE_GS, dataname= "GS")
dev.off()


####MAKE HEATMAPS WITH THE NUMBER OF DCEs AMONG CHROMOSOMES FOR EACH SUBTYPE################

df <- data.frame("Number_of_cods" = numeric(), "study_group" = character(), "type_of_rearrangement" = character(), "chromosome" = character(), 
                 stringsAsFactors = F)
for (chr in names(rearrangements[[1]])){
  for (group in names(rearrangements)){
    Number_of_cods <- c(length(rearrangements[[group]][[chr]][["depleted"]]), length(rearrangements[[group]][[chr]][["emerged"]]),
                        length(unique(rearrangements[[group]][[chr]][["split"]][, "control"])),
                        length(unique(rearrangements[[group]][[chr]][["merged"]][, "control"])),
                        length(rearrangements[[group]][[chr]][["intact"]]),
                        nrow(rearrangements[[group]][[chr]][["expanded"]]),
                        nrow(rearrangements[[group]][[chr]][["contracted"]]),
                        nrow(rearrangements[[group]][[chr]][["shift"]]))
    
    study_group <- rep(group, length(Number_of_cods))
    type_of_rearrangement <- c("depleted", "emerged", "split", "merged", "intact", "expanded", "contracted", "shift")
    chromosome <- rep(chr, length(Number_of_cods))
    df <- rbind(df, cbind(Number_of_cods, study_group, type_of_rearrangement, chromosome))
  }
}
#print(df)
df$Number_of_cods <- as.integer(as.character(df$Number_of_cods))
df$type_of_rearrangement <- factor(as.character(df$type_of_rearrangement), 
                                   levels =  c("depleted", "emerged", "split", "merged", "intact", "expanded", "contracted", "shift"))

length(df[df$study_group == "EBV",])
#length(df[df$study_group == "NP",])
length(df[df$study_group == "CIN",])
length(df[df$study_group == "MSI",])
length(df[df$study_group == "GS",])

num_rear_chrs_EBV <- df[df$study_group == "EBV",]
#num_rear_chrs_NP <- df[df$study_group == "NP",]
num_rear_chrs_CIN <- df[df$study_group == "CIN",]
num_rear_chrs_MSI <- df[df$study_group == "MSI",]
num_rear_chrs_GS <- df[df$study_group == "GS",]


#par(mfrow = c(2,3))

pheatmap_cods_EBV <- num_rear_chrs_EBV %>% 
  pivot_wider(
    values_from = Number_of_cods,
    id_cols = chromosome,
    names_from = type_of_rearrangement
  )

#pheatmap_cods_EBV <- as.data.frame(pheatmap_cods_EBV)
pheatmap_cods_EBV <- pheatmap_cods_EBV %>% remove_rownames %>% column_to_rownames(var = "chromosome")
pheatmap(t(pheatmap_cods_EBV),
         fontsize_row = 14,
         fontsize_col = 14,
         main = "   DCE reorganization patterns among chromosomes (EBV)
         ",
         fontsize = 14)

#NP
pheatmap_cods_NP <- num_rear_chrs_NP %>% 
  pivot_wider(
    values_from = Number_of_cods,
    id_cols = chromosome,
    names_from = type_of_rearrangement
  )

#pheatmap_cods_EBV <- as.data.frame(pheatmap_cods_EBV)
pheatmap_cods_NP <- pheatmap_cods_NP %>% remove_rownames %>% column_to_rownames(var = "chromosome")
pheatmap(t(pheatmap_cods_NP),
         fontsize_row = 14,
         fontsize_col = 14,
         main = "DCE reorganization patterns among chromosomes (NP)
         ",
         fontsize = 14)

#CIN
pheatmap_cods_CIN <- num_rear_chrs_CIN %>% 
  pivot_wider(
    values_from = Number_of_cods,
    id_cols = chromosome,
    names_from = type_of_rearrangement
  )

#pheatmap_cods_EBV <- as.data.frame(pheatmap_cods_EBV)
pheatmap_cods_CIN <- pheatmap_cods_CIN %>% remove_rownames %>% column_to_rownames(var = "chromosome")
pheatmap(t(pheatmap_cods_CIN),
         fontsize_row = 14,
         fontsize_col = 14,
         main = "   DCE reorganization patterns among chromosomes (CIN)
         ",
         fontsize = 14)

#MSI

pheatmap_cods_MSI <- num_rear_chrs_MSI %>% 
  pivot_wider(
    values_from = Number_of_cods,
    id_cols = chromosome,
    names_from = type_of_rearrangement
  )

#pheatmap_cods_EBV <- as.data.frame(pheatmap_cods_EBV)
pheatmap_cods_MSI <- pheatmap_cods_MSI %>% remove_rownames %>% column_to_rownames(var = "chromosome")
pheatmap(t(pheatmap_cods_MSI),
         fontsize_row = 14,
         fontsize_col = 14,
         main = "   DCE reorganization patterns among chromosomes (MSI)
         ",
         fontsize = 14)

#GS
pheatmap_cods_GS <- num_rear_chrs_GS %>% 
  pivot_wider(
    values_from = Number_of_cods,
    id_cols = chromosome,
    names_from = type_of_rearrangement
  )

#pheatmap_cods_EBV <- as.data.frame(pheatmap_cods_EBV)
pheatmap_cods_GS <- pheatmap_cods_GS %>% remove_rownames %>% column_to_rownames(var = "chromosome")
pheatmap(t(pheatmap_cods_GS),
         fontsize_row = 14,
         fontsize_col = 14,
         main = "   DCE reorganization patterns among chromosomes (GS)
         ",
         fontsize = 14)




##########plot the rearrangement occurence for each subtype compared to NP
df <- tapply(X = as.numeric(as.character(df$Number_of_cods)), INDEX = list(df$study_group, df$type_of_rearrangement), FUN = sum)
class(df)
df <- t(df)

rear <- as.data.frame(df)

rear <- cbind.data.frame( "Type of rearrangement" = c("depleted",
                                                      "emerged",
                                                      "split",
                                                      "merged",
                                                      "intact",
                                                      "expanded",
                                                      "contracted",
                                                      "shift"), df)

rear <- rear %>% 
  gather(subtype, value, 
         -`Type of rearrangement`)



rear$`Type of rearrangement` <- factor(rear$`Type of rearrangement`, levels = c("depleted",
                                                                                "emerged",
                                                                                "split",
                                                                                "merged",
                                                                                "intact",
                                                                                "expanded",
                                                                                "contracted",
                                                                                "shift"))



rearrangements_g <- rear
rearrangements_g$subtype <- factor(rearrangements_g$subtype, levels = c("CIN",
                                                                        "EBV",
                                                                        "MSI",
                                                                        "GS"))


rearrangements_g$`Type of rearrangement` <- factor(rearrangements_g$`Type of rearrangement`, levels = c("depleted",
                                                                                                        "emerged",
                                                                                                        "split",
                                                                                                        "merged",
                                                                                                        "intact",
                                                                                                        "expanded",
                                                                                                        "contracted",
                                                                                                        "shift"))


svg("svg_plots/rearrangements_absolute.svg", width=10, height=8)
ggplot(rearrangements_g, aes(x = subtype, y = value, fill = rearrangements_g$`Type of rearrangement`)) + 
  geom_bar(stat='identity', position = "stack") +
  ggtitle("Number of rearrangement occurence
per subtype compared to NP")+
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15))+
  ylab("Rearrangement occurence")+
  guides(fill=guide_legend(title="Type of rearrangement"))
dev.off()


#rearrangements heatmaps
svg("svg_plots/rearrangements_heatmaps.svg", width=10, height=8)

as.ggplot(
  function(){
    pheatmap::pheatmap(fontsize = size,
                       apply(df, 2, function(x) x/sum(x)), scale = "column", cluster_cols = F,
                       clustering_method = "ward.D2", cutree_rows = 3,
                       display_numbers = apply(apply(df, 2, function(x) x/sum(x)), 2, function(x) as.character(round(x,3))))
  }
) +
  theme(text = element_text( size = size),
        plot.tag = element_text(size = tag_size, face = "bold"), legend.position = "none")

dev.off()



#plot a stacked barplot with the percentage of rearrangements occuring in each subtype
rear <- as.data.frame(df)
rearprop <- cbind.data.frame(
  rear$CIN/sum(rear$CIN)*100,
  rear$EBV/sum(rear$EBV)*100,
  rear$MSI/sum(rear$MSI)*100,
  rear$GS/sum(rear$GS)*100
)


colnames(rearprop) <- c("CIN",
                        "EBV",
                        "MSI",
                        "GS")

rownames(rearprop) <- c("depleted",
                        "emerged",
                        "split",
                        "merged",
                        "intact",
                        "expanded",
                        "contracted",
                        "shift")



rearprop <- cbind.data.frame( "Type of rearrangement" = c("depleted",
                                                          "emerged",
                                                          "split",
                                                          "merged",
                                                          "intact",
                                                          "expanded",
                                                          "contracted",
                                                          "shift"), rearprop)




rearprop <- rearprop %>% 
  gather(subtype, value, 
         -`Type of rearrangement`)

rearprop$`Type of rearrangement` <- factor(rearprop$`Type of rearrangement`, levels = c("depleted",
                                                                                        "emerged",
                                                                                        "split",
                                                                                        "merged",
                                                                                        "intact",
                                                                                        "expanded",
                                                                                        "contracted",
                                                                                        "shift"))


rearprop$subtype <- factor(rearprop$subtype, levels = c("NP",
                                                        "EBV",
                                                        "CIN",
                                                        "MSI",
                                                        "GS"))



svg("svg_plots/rearrangements_percentages_GC.svg", width=10, height=8)
ggplot(rearprop, aes(x = subtype, y = value, fill = rearprop$`Type of rearrangement`)) + 
  geom_bar(stat='identity', position = "stack") +
  ggtitle("   DCE pattern reorganization 
  occurence per GC subtype")+
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15))+
  ylab("Pattern reorganization occurence (%)")+
  guides(fill=guide_legend(title="Type of reorganization"))
dev.off()

rearprop <- rearprop %>% remove_rownames() %>% column_to_rownames(var = 'Type of rearrangement')

pheatmap(rearprop,
         display_numbers = TRUE,
         main = "Rearrangements proportions
         ",
         number_color = "black", 
         fontsize_number = 8
)

######################PLOT GENE DENSITY#########################################

coordinates <- read.csv("tcga_annotation_data_GC_ensids.txt", sep = ",", header = T)
colnames(coordinates) <- c("chromosome", "start", "end", "strand","gene_name", "gene_type","ensid", "ensid_version", "gene_description")
coordinates$chromosome <- paste("chr", coordinates$chromosome, sep = '')

#plot the genes to infer gene density (use the genes from the normal GTEx data)
# coordinates <- read.delim("GTEx_normal_gc_gene_information.txt", sep = "\t", header = T)
# coordinates$Chromosome.scaffold.name <- paste("chr", coordinates$Chromosome.scaffold.name, sep = '')
#kp <- plotKaryotype(genome="hg38", chromosomes = c("chr13","chr14","chr15"))

genes <- makeGRangesFromDataFrame(data.frame(
  chr = coordinates$chromosome,
  start = coordinates$start,
  end = coordinates$end
))

svg("svg_plots/gene_density_genome_tcga.svg", width = 10, height = 10)
kp <- plotKaryotype(genome="hg38") #chromosomes = c("chr13","chr14","chr15")
kpPlotDensity(kp,genes)
kpAddMainTitle(kp, main = "Gene density", cex=2)
dev.off()

svg("svg_plots/gene_density_TCGA_chr13_14_15_22.svg", width = 8, height = 8)
kp <- plotKaryotype(genome="hg38", chromosomes = c("chr13","chr14","chr15","chr22"))
kpPlotDensity(kp,genes)
kpAddMainTitle(kp, main = "Gene density: Chromosomes 13,14,15,22", cex=2)
dev.off()


########Find intra DCE coexpression########
####WARNING!!!!!! Before you run the scripts above, make sure you have run the script for the intra DCE coexpression in the DCE_statistics.R script!!!
intra_cod1 <- data.frame(intra_cod_correlation = numeric(), group = character(), chromosome = character())

for (chr in names(intra_cod$NP)){
  
  intra_cod2 <- list()
  for (group in names(intra_cod)){
    intra_cod2[[group]] <- intra_cod[[group]][[chr]]
  }
  intra_cod1 <- rbind(intra_cod1, 
                      data.frame(intra_cod_correlation =c(intra_cod2$NP, intra_cod2$EBV, intra_cod2$CIN, intra_cod2$MSI, intra_cod2$GS), 
                                 group = c(rep("NP", length(intra_cod2$NP)), rep("EBV", length(intra_cod2$EBV)), 
                                           rep("CIN", length(intra_cod2$CIN)), rep("MSI", length(intra_cod2$MSI)), rep("GS", length(intra_cod2$GS))),
                                 chromosome = rep(chr, length(unlist(intra_cod2)))))
  
}

#This is a script for each chromosome separately
# ggplot(intra_cod1, aes(x=group, y=intra_cod_correlation, fill=group)) +
#   geom_violin(trim=FALSE)+
#   geom_boxplot(width=0.1) + 
#   stat_compare_means(ref.group = "NP", label = "p.signif", size = 20, label.y.npc = "top", col = "firebrick") + 
#   facet_wrap(~ chromosome, ncol = 6, scales = "free") +
#   theme(text = element_text(size = 50), axis.text.x = element_text(size = 20), panel.spacing = unit(3, "lines"), legend.text = element_text(face = "bold"),
#         legend.title = element_text(face = "bold"), legend.key.size = unit(x = 15, units = "mm"), legend.box.background = element_rect()) 
# dev.off()

intra_cod <- lapply(intra_cod, unlist, use.names = F)
intra_cod1 <- data.frame(intra_cod_correlation =c(intra_cod$NP, intra_cod$EBV, intra_cod$CIN, intra_cod$MSI, intra_cod$GS),
                         group = c(rep("NP", length(intra_cod$NP)), rep("EBV", length(intra_cod$EBV)),
                                   rep("CIN", length(intra_cod$CIN)), rep("MSI", length(intra_cod$MSI)), rep("GS", length(intra_cod$GS)) ))
intra_cod1$group <- factor(intra_cod1$group, levels = my_levels)

library(ggsci)
svg("svg_plots/Intra_DCEs_g.svg", height = 8, width = 8)
ggplot(intra_cod1, aes(x=group, y=intra_cod_correlation, fill=group)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1) +
  scale_fill_manual(values = pal_uchicago()(8)[c(3,2,1,5,7,8)]) +
  geom_hline(yintercept = mean(intra_cod1$intra_cod_correlation), linetype = 2) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  stat_compare_means(method = "kruskal.test", size = 5, label.y = 1.1, label.x = 4)+
theme(legend.position = "none",
      plot.title = element_text(size = 25, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15)) +
  labs(title="Intra-DCE correlation in each GC subtype", 
       y = "Correlation")
dev.off()


########Find inter DCE coexpression########
####WARNING!!!!!! Before you run the scripts above, make sure you have run the script for the intra DCE coexpression in the DCE_statistics.R script!!!
inter_cod1 <- data.frame(inter_cod_correlation = numeric(), group = character(), chromosome = character())

for (chr in names(inter_cod$NP)){
  
  inter_cod2 <- list()
  for (group in names(inter_cod)){
    inter_cod2[[group]] <- inter_cod[[group]][[chr]]
  }
  inter_cod1 <- rbind(inter_cod1, 
                      data.frame(inter_cod_correlation =c(inter_cod2$NP, inter_cod2$EBV, inter_cod2$CIN, inter_cod2$MSI, inter_cod2$GS), 
                                 group = c(rep("NP", length(inter_cod2$NP)), rep("EBV", length(inter_cod2$EBV)), 
                                           rep("CIN", length(inter_cod2$CIN)), rep("MSI", length(inter_cod2$MSI)), rep("GS", length(inter_cod2$GS))),
                                 chromosome = rep(chr, length(unlist(inter_cod2)))))
  
}
inter_cod1$group <- factor(inter_cod1$group, levels = my_levels)

#This is a script for each chromosome separately

# p <- ggplot(inter_cod1, aes(x=group, y=inter_cod_correlation, fill=group)) +
#   geom_violin(trim=FALSE)
# p + geom_boxplot(width=0.1) + 
#   stat_compare_means(ref.group = "NP", label = "p.signif", size = 15, label.y.npc = "top", col = "firebrick") + 
#   facet_wrap(~ chromosome, ncol = 6, scales = "free") +
#   theme(text = element_text(size = 50), axis.text.x = element_text(size = 20), panel.spacing = unit(3, "lines"), legend.text = element_text(face = "bold"),
#         legend.title = element_text(face = "bold"), legend.key.size = unit(x = 15, units = "mm"), legend.box.background = element_rect()) 
# 
# dev.off()

inter_cod <- lapply(inter_cod, unlist, use.names = F)
inter_cod1 <- data.frame(inter_cod_correlation =c(inter_cod$NP, inter_cod$EBV, inter_cod$CIN, inter_cod$MSI, inter_cod$GS),
                         group = c(rep("NP", length(inter_cod$NP)), rep("EBV", length(inter_cod$EBV)),
                                   rep("CIN", length(inter_cod$CIN)), rep("MSI", length(inter_cod$MSI)), rep("GS", length(inter_cod$GS))))
inter_cod1$group <- factor(inter_cod1$group, levels = my_levels)


svg("svg_plots/Inter_DCEs_g.svg", height = 8, width = 8)
ggplot(inter_cod1, aes(x=group, y=inter_cod_correlation, fill=group)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1) +
  scale_fill_manual(values = pal_uchicago()(8)[c(3,2,1,5,7,8)]) +
  geom_hline(yintercept = mean(inter_cod1$inter_cod_correlation), linetype = 2) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  stat_compare_means(method = "kruskal.test", size = 5, label.y = 1, label.x = 4)+
  theme(legend.position = "none",
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  labs(title="Inter-DCE correlation in each GC subtype", 
       y = "Correlation")
dev.off()


##########GENOME COVERAGE BY DCEs (HAPLOID)############# 
#you will need the sum of chromosome length data (the current version is hg38p14) (provided as an excel file)
s_CODs <- soCODs(CODs)
sapply(s_CODs, function(x)summary(unlist(x, use.names = F)))
s_CODs <- lapply(s_CODs, unlist, use.names = F)
gen_cover_perc <- data.frame()
#gen_cover_perc <- t(gen_cover_perc)
for (i in names(s_CODs)){
  x <- data.frame(sum(s_CODs[[i]])/sum(chromosomes_hg38p14$...3))*100

gen_cover_perc <- rbind(gen_cover_perc, x)
}
gen_cover_perc$subtypes <- names(s_CODs) 
colnames(gen_cover_perc) <- c("genome_cover","subtypes")
gen_cover_perc$subtypes <- factor(gen_cover_perc$subtypes, levels = my_levels)

#barplot
svg("svg_plots/DCE_haploid_genome_coverage.svg", width=8, height=8)
ggplot(gen_cover_perc, aes(x = subtypes, y = genome_cover, fill = subtypes, color = "black")) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = pal_uchicago()(8)[c(3,2,1,5,7,8)]) +
  theme(legend.position = "NONE",
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  labs(title="DCE genome coverage (haploid)
in each gastric cancer subtype",
       x ="Subtypes", 
       y = "Coverage Percentage (%)")
dev.off()


######DCE Chromosome coverage ##########
s_CODs <- soCODs(CODs)
sapply(s_CODs, function(x)summary(unlist(x, use.names = F)))

chr1 <- list( 
  s_CODs[["NP"]][["chr1"]],
  s_CODs[["CIN"]][["chr1"]],
  s_CODs[["EBV"]][["chr1"]], 
  s_CODs[["MSI"]][["chr1"]],
  s_CODs[["GS"]][["chr1"]])
names(chr1) <- c("NP", "CIN", "EBV", "MSI","GS")

chr2 <- list(
  s_CODs[["NP"]][["chr2"]],
  s_CODs[["CIN"]][["chr2"]],
  s_CODs[["EBV"]][["chr2"]], 
  s_CODs[["MSI"]][["chr2"]],
  s_CODs[["GS"]][["chr2"]])
names(chr2) <- c( "NP", "CIN", "EBV", "MSI", "GS")

chr3 <- list(
  s_CODs[["NP"]][["chr3"]],
  s_CODs[["CIN"]][["chr3"]],
  s_CODs[["EBV"]][["chr3"]], 
  s_CODs[["MSI"]][["chr3"]],
  s_CODs[["GS"]][["chr3"]])
names(chr3) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr4 <- list(
  s_CODs[["NP"]][["chr4"]],
  s_CODs[["CIN"]][["chr4"]],
  s_CODs[["EBV"]][["chr4"]], 
  s_CODs[["MSI"]][["chr4"]],
  s_CODs[["GS"]][["chr4"]])
names(chr4) <- c("NP", "CIN", "EBV", "MSI","GS")


chr5 <- list( 
  s_CODs[["NP"]][["chr5"]],
  s_CODs[["CIN"]][["chr5"]],
  s_CODs[["EBV"]][["chr5"]], 
  s_CODs[["MSI"]][["chr5"]],
  s_CODs[["GS"]][["chr5"]])
names(chr5) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr6 <- list(s_CODs[["NP"]][["chr6"]], 
             s_CODs[["CIN"]][["chr6"]],
             s_CODs[["EBV"]][["chr6"]], 
             s_CODs[["MSI"]][["chr6"]],
             s_CODs[["GS"]][["chr6"]])
names(chr6) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr7 <- list(s_CODs[["NP"]][["chr7"]],
             s_CODs[["CIN"]][["chr7"]],
             s_CODs[["EBV"]][["chr7"]], 
             s_CODs[["MSI"]][["chr7"]],
             s_CODs[["GS"]][["chr7"]])
names(chr7) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr8 <- list(s_CODs[["NP"]][["chr8"]],
             s_CODs[["CIN"]][["chr8"]],
             s_CODs[["EBV"]][["chr8"]], 
             s_CODs[["MSI"]][["chr8"]],
             s_CODs[["GS"]][["chr8"]])
names(chr8) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr9 <- list(s_CODs[["NP"]][["chr9"]],
             s_CODs[["CIN"]][["chr9"]],
             s_CODs[["EBV"]][["chr9"]], 
             s_CODs[["MSI"]][["chr9"]],
             s_CODs[["GS"]][["chr9"]])
names(chr9) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr10 <- list(s_CODs[["NP"]][["chr10"]],
              s_CODs[["CIN"]][["chr10"]],
              s_CODs[["EBV"]][["chr10"]], 
              s_CODs[["MSI"]][["chr10"]],
              s_CODs[["GS"]][["chr10"]])
names(chr10) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr11 <- list(s_CODs[["NP"]][["chr11"]],
              s_CODs[["CIN"]][["chr11"]],
              s_CODs[["EBV"]][["chr11"]], 
              s_CODs[["MSI"]][["chr11"]],
              s_CODs[["GS"]][["chr11"]])
names(chr11) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr12 <- list(s_CODs[["NP"]][["chr12"]],
              s_CODs[["CIN"]][["chr12"]],
              s_CODs[["EBV"]][["chr12"]], 
              s_CODs[["MSI"]][["chr12"]],
              s_CODs[["GS"]][["chr12"]])
names(chr12) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr13 <- list(s_CODs[["NP"]][["chr13"]],
              s_CODs[["CIN"]][["chr13"]],
              s_CODs[["EBV"]][["chr13"]], 
              s_CODs[["MSI"]][["chr13"]],
              s_CODs[["GS"]][["chr13"]])
names(chr13) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr14 <- list(s_CODs[["NP"]][["chr14"]],
              s_CODs[["CIN"]][["chr14"]],
              s_CODs[["EBV"]][["chr14"]], 
              s_CODs[["MSI"]][["chr14"]],
              s_CODs[["GS"]][["chr14"]])
names(chr14) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr15 <- list(s_CODs[["NP"]][["chr15"]],
              s_CODs[["CIN"]][["chr15"]],
              s_CODs[["EBV"]][["chr15"]], 
              s_CODs[["MSI"]][["chr15"]],
              s_CODs[["GS"]][["chr15"]])
names(chr15) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr16 <- list(s_CODs[["NP"]][["chr16"]],
              s_CODs[["CIN"]][["chr16"]],
              s_CODs[["EBV"]][["chr16"]], 
              s_CODs[["MSI"]][["chr16"]],
              s_CODs[["GS"]][["chr16"]])
names(chr16) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr17 <- list(s_CODs[["NP"]][["chr17"]],
              s_CODs[["CIN"]][["chr17"]],
              s_CODs[["EBV"]][["chr17"]], 
              s_CODs[["MSI"]][["chr17"]],
              s_CODs[["GS"]][["chr17"]])
names(chr17) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr18 <- list(s_CODs[["NP"]][["chr18"]],
              s_CODs[["CIN"]][["chr18"]],
              s_CODs[["EBV"]][["chr18"]], 
              s_CODs[["MSI"]][["chr18"]],
              s_CODs[["GS"]][["chr18"]])
names(chr18) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr19 <- list(s_CODs[["NP"]][["chr19"]],
              s_CODs[["CIN"]][["chr19"]],
              s_CODs[["EBV"]][["chr19"]], 
              s_CODs[["MSI"]][["chr19"]],
              s_CODs[["GS"]][["chr19"]])
names(chr19) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr20 <- list(s_CODs[["NP"]][["chr20"]],
              s_CODs[["CIN"]][["chr20"]],
              s_CODs[["EBV"]][["chr20"]], 
              s_CODs[["MSI"]][["chr20"]],
              s_CODs[["GS"]][["chr20"]])
names(chr20) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr21 <- list(s_CODs[["NP"]][["chr21"]],
              s_CODs[["CIN"]][["chr21"]],
              s_CODs[["EBV"]][["chr21"]], 
              s_CODs[["MSI"]][["chr21"]],
              s_CODs[["GS"]][["chr21"]])
names(chr21) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr22 <- list(s_CODs[["NP"]][["chr22"]],
              s_CODs[["CIN"]][["chr22"]],
              s_CODs[["EBV"]][["chr22"]], 
              s_CODs[["MSI"]][["chr22"]],
              s_CODs[["GS"]][["chr22"]])
names(chr22) <- c("NP", "CIN", "EBV", "MSI", "GS")


chrX <- list(s_CODs[["NP"]][["chrX"]],
             s_CODs[["CIN"]][["chrX"]],
             s_CODs[["EBV"]][["chrX"]], 
             s_CODs[["MSI"]][["chrX"]],
             s_CODs[["GS"]][["chrX"]])
names(chrX) <- c("NP", "CIN", "EBV", "MSI", "GS")

chrY <- list(s_CODs[["NP"]][["chrY"]],
             s_CODs[["CIN"]][["chrY"]],
             s_CODs[["EBV"]][["chrY"]], 
             s_CODs[["MSI"]][["chrY"]],
             s_CODs[["GS"]][["chrY"]])
names(chrY) <- c("NP", "CIN", "EBV", "MSI", "GS")


chr_dce_coverage <- list(chr1,chr2,chr3,
                         chr4,chr5,chr6,
                         chr7,chr8,chr9,
                         chr10,chr11,chr12,
                         chr13,chr14,chr15,
                         chr16,chr17,chr18,
                         chr19,chr20,chr21,
                         chr22,chrX,chrY)

rm(chr1,chr2,chr3,
   chr4,chr5,chr6,
   chr7,chr8,chr9,
   chr10,chr11,chr12,
   chr13,chr14,chr15,
   chr16,chr17,chr18,
   chr19,chr20,chr21,
   chr22,chrX,chrY)


chrordered <- c("chr1","chr2","chr3",
                "chr4","chr5","chr6",
                "chr7","chr8","chr9",
                "chr10","chr11","chr12",
                "chr13","chr14","chr15",
                "chr16","chr17","chr18",
                "chr19","chr20","chr21",
                "chr22","chrX","chrY")

names(chr_dce_coverage) <- chrordered
rm(chrordered)



#ALL together

#NP-like
dce_sizes_NP <- cbind.data.frame(
  sum(chr_dce_coverage[["chr1"]][["NP"]]),
  sum(chr_dce_coverage[["chr2"]][["NP"]]),
  sum(chr_dce_coverage[["chr3"]][["NP"]]),
  sum(chr_dce_coverage[["chr4"]][["NP"]]),
  sum(chr_dce_coverage[["chr5"]][["NP"]]),
  sum(chr_dce_coverage[["chr6"]][["NP"]]),
  sum(chr_dce_coverage[["chr7"]][["NP"]]),
  sum(chr_dce_coverage[["chr8"]][["NP"]]),
  sum(chr_dce_coverage[["chr9"]][["NP"]]),
  sum(chr_dce_coverage[["chr10"]][["NP"]]),
  sum(chr_dce_coverage[["chr11"]][["NP"]]),
  sum(chr_dce_coverage[["chr12"]][["NP"]]),
  sum(chr_dce_coverage[["chr13"]][["NP"]]),
  sum(chr_dce_coverage[["chr14"]][["NP"]]),
  sum(chr_dce_coverage[["chr15"]][["NP"]]),
  sum(chr_dce_coverage[["chr16"]][["NP"]]),
  sum(chr_dce_coverage[["chr17"]][["NP"]]),
  sum(chr_dce_coverage[["chr18"]][["NP"]]),
  sum(chr_dce_coverage[["chr19"]][["NP"]]),
  sum(chr_dce_coverage[["chr20"]][["NP"]]),
  sum(chr_dce_coverage[["chr21"]][["NP"]]),
  sum(chr_dce_coverage[["chr22"]][["NP"]]),
  sum(chr_dce_coverage[["chrX"]][["NP"]]),
  sum(chr_dce_coverage[["chrY"]][["NP"]]))



colnames(dce_sizes_NP) <- c("chr1","chr2","chr3",
                            "chr4","chr5","chr6",
                            "chr7","chr8","chr9",
                            "chr10","chr11","chr12",
                            "chr13","chr14","chr15",
                            "chr16","chr17","chr18",
                            "chr19","chr20","chr21",
                            "chr22","chrX","chrY")
row.names(dce_sizes_NP) <- "sum of dce sizes per chr"


#NP
dce_sizes_NP <- cbind.data.frame(
  sum(chr_dce_coverage[["chr1"]][["NP"]]),sum(chr_dce_coverage[["chr2"]][["NP"]]),sum(chr_dce_coverage[["chr3"]][["NP"]]),
  sum(chr_dce_coverage[["chr4"]][["NP"]]),sum(chr_dce_coverage[["chr5"]][["NP"]]),sum(chr_dce_coverage[["chr6"]][["NP"]]),
  sum(chr_dce_coverage[["chr7"]][["NP"]]),sum(chr_dce_coverage[["chr8"]][["NP"]]),sum(chr_dce_coverage[["chr9"]][["NP"]]),
  sum(chr_dce_coverage[["chr10"]][["NP"]]),sum(chr_dce_coverage[["chr11"]][["NP"]]),sum(chr_dce_coverage[["chr12"]][["NP"]]),
  sum(chr_dce_coverage[["chr13"]][["NP"]]),sum(chr_dce_coverage[["chr14"]][["NP"]]),sum(chr_dce_coverage[["chr15"]][["NP"]]),
  sum(chr_dce_coverage[["chr16"]][["NP"]]),sum(chr_dce_coverage[["chr17"]][["NP"]]),sum(chr_dce_coverage[["chr18"]][["NP"]]),
  sum(chr_dce_coverage[["chr19"]][["NP"]]),sum(chr_dce_coverage[["chr20"]][["NP"]]),sum(chr_dce_coverage[["chr21"]][["NP"]]),
  sum(chr_dce_coverage[["chr22"]][["NP"]]),sum(chr_dce_coverage[["chrX"]][["NP"]]),sum(chr_dce_coverage[["chrY"]][["NP"]]))



colnames(dce_sizes_NP) <- c("chr1","chr2","chr3",
                            "chr4","chr5","chr6",
                            "chr7","chr8","chr9",
                            "chr10","chr11","chr12",
                            "chr13","chr14","chr15",
                            "chr16","chr17","chr18",
                            "chr19","chr20","chr21",
                            "chr22","chrX","chrY")
row.names(dce_sizes_NP) <- "sum of dce sizes per chr"

#CIN
dce_sizes_CIN <- cbind.data.frame(
  sum(chr_dce_coverage[["chr1"]][["CIN"]]),sum(chr_dce_coverage[["chr2"]][["CIN"]]),sum(chr_dce_coverage[["chr3"]][["CIN"]]),
  sum(chr_dce_coverage[["chr4"]][["CIN"]]),sum(chr_dce_coverage[["chr5"]][["CIN"]]),sum(chr_dce_coverage[["chr6"]][["CIN"]]),
  sum(chr_dce_coverage[["chr7"]][["CIN"]]),sum(chr_dce_coverage[["chr8"]][["CIN"]]),sum(chr_dce_coverage[["chr9"]][["CIN"]]),
  sum(chr_dce_coverage[["chr10"]][["CIN"]]),sum(chr_dce_coverage[["chr11"]][["CIN"]]),sum(chr_dce_coverage[["chr12"]][["CIN"]]),
  sum(chr_dce_coverage[["chr13"]][["CIN"]]),sum(chr_dce_coverage[["chr14"]][["CIN"]]),sum(chr_dce_coverage[["chr15"]][["CIN"]]),
  sum(chr_dce_coverage[["chr16"]][["CIN"]]),sum(chr_dce_coverage[["chr17"]][["CIN"]]),sum(chr_dce_coverage[["chr18"]][["CIN"]]),
  sum(chr_dce_coverage[["chr19"]][["CIN"]]),sum(chr_dce_coverage[["chr20"]][["CIN"]]),sum(chr_dce_coverage[["chr21"]][["CIN"]]),
  sum(chr_dce_coverage[["chr22"]][["CIN"]]),sum(chr_dce_coverage[["chrX"]][["CIN"]]),sum(chr_dce_coverage[["chrY"]][["CIN"]]))



colnames(dce_sizes_CIN) <- c("chr1","chr2","chr3",
                             "chr4","chr5","chr6",
                             "chr7","chr8","chr9",
                             "chr10","chr11","chr12",
                             "chr13","chr14","chr15",
                             "chr16","chr17","chr18",
                             "chr19","chr20","chr21",
                             "chr22","chrX","chrY")
row.names(dce_sizes_CIN) <- "sum of dce sizes per chr"

#EBV
dce_sizes_EBV <- cbind.data.frame(
  sum(chr_dce_coverage[["chr1"]][["EBV"]]),sum(chr_dce_coverage[["chr2"]][["EBV"]]),sum(chr_dce_coverage[["chr3"]][["EBV"]]),
  sum(chr_dce_coverage[["chr4"]][["EBV"]]),sum(chr_dce_coverage[["chr5"]][["EBV"]]),sum(chr_dce_coverage[["chr6"]][["EBV"]]),
  sum(chr_dce_coverage[["chr7"]][["EBV"]]),sum(chr_dce_coverage[["chr8"]][["EBV"]]),sum(chr_dce_coverage[["chr9"]][["EBV"]]),
  sum(chr_dce_coverage[["chr10"]][["EBV"]]),sum(chr_dce_coverage[["chr11"]][["EBV"]]),sum(chr_dce_coverage[["chr12"]][["EBV"]]),
  sum(chr_dce_coverage[["chr13"]][["EBV"]]),sum(chr_dce_coverage[["chr14"]][["EBV"]]),sum(chr_dce_coverage[["chr15"]][["EBV"]]),
  sum(chr_dce_coverage[["chr16"]][["EBV"]]),sum(chr_dce_coverage[["chr17"]][["EBV"]]),sum(chr_dce_coverage[["chr18"]][["EBV"]]),
  sum(chr_dce_coverage[["chr19"]][["EBV"]]),sum(chr_dce_coverage[["chr20"]][["EBV"]]),sum(chr_dce_coverage[["chr21"]][["EBV"]]),
  sum(chr_dce_coverage[["chr22"]][["EBV"]]),sum(chr_dce_coverage[["chrX"]][["EBV"]]),sum(chr_dce_coverage[["chrY"]][["EBV"]]))



colnames(dce_sizes_EBV) <- c("chr1","chr2","chr3",
                             "chr4","chr5","chr6",
                             "chr7","chr8","chr9",
                             "chr10","chr11","chr12",
                             "chr13","chr14","chr15",
                             "chr16","chr17","chr18",
                             "chr19","chr20","chr21",
                             "chr22","chrX","chrY")
row.names(dce_sizes_EBV) <- "sum of dce sizes per chr"

#MSI

dce_sizes_MSI <- cbind.data.frame(
  sum(chr_dce_coverage[["chr1"]][["MSI"]]),sum(chr_dce_coverage[["chr2"]][["MSI"]]),sum(chr_dce_coverage[["chr3"]][["MSI"]]),
  sum(chr_dce_coverage[["chr4"]][["MSI"]]),sum(chr_dce_coverage[["chr5"]][["MSI"]]),sum(chr_dce_coverage[["chr6"]][["MSI"]]),
  sum(chr_dce_coverage[["chr7"]][["MSI"]]),sum(chr_dce_coverage[["chr8"]][["MSI"]]),sum(chr_dce_coverage[["chr9"]][["MSI"]]),
  sum(chr_dce_coverage[["chr10"]][["MSI"]]),sum(chr_dce_coverage[["chr11"]][["MSI"]]),sum(chr_dce_coverage[["chr12"]][["MSI"]]),
  sum(chr_dce_coverage[["chr13"]][["MSI"]]),sum(chr_dce_coverage[["chr14"]][["MSI"]]),sum(chr_dce_coverage[["chr15"]][["MSI"]]),
  sum(chr_dce_coverage[["chr16"]][["MSI"]]),sum(chr_dce_coverage[["chr17"]][["MSI"]]),sum(chr_dce_coverage[["chr18"]][["MSI"]]),
  sum(chr_dce_coverage[["chr19"]][["MSI"]]),sum(chr_dce_coverage[["chr20"]][["MSI"]]),sum(chr_dce_coverage[["chr21"]][["MSI"]]),
  sum(chr_dce_coverage[["chr22"]][["MSI"]]),sum(chr_dce_coverage[["chrX"]][["MSI"]]),sum(chr_dce_coverage[["chrY"]][["MSI"]]))



colnames(dce_sizes_MSI) <- c("chr1","chr2","chr3",
                             "chr4","chr5","chr6",
                             "chr7","chr8","chr9",
                             "chr10","chr11","chr12",
                             "chr13","chr14","chr15",
                             "chr16","chr17","chr18",
                             "chr19","chr20","chr21",
                             "chr22","chrX","chrY")
row.names(dce_sizes_MSI) <- "sum of dce sizes per chr"

#GS
dce_sizes_GS <- cbind.data.frame(
  sum(chr_dce_coverage[["chr1"]][["GS"]]),sum(chr_dce_coverage[["chr2"]][["GS"]]),sum(chr_dce_coverage[["chr3"]][["GS"]]),
  sum(chr_dce_coverage[["chr4"]][["GS"]]),sum(chr_dce_coverage[["chr5"]][["GS"]]),sum(chr_dce_coverage[["chr6"]][["GS"]]),
  sum(chr_dce_coverage[["chr7"]][["GS"]]),sum(chr_dce_coverage[["chr8"]][["GS"]]),sum(chr_dce_coverage[["chr9"]][["GS"]]),
  sum(chr_dce_coverage[["chr10"]][["GS"]]),sum(chr_dce_coverage[["chr11"]][["GS"]]),sum(chr_dce_coverage[["chr12"]][["GS"]]),
  sum(chr_dce_coverage[["chr13"]][["GS"]]),sum(chr_dce_coverage[["chr14"]][["GS"]]),sum(chr_dce_coverage[["chr15"]][["GS"]]),
  sum(chr_dce_coverage[["chr16"]][["GS"]]),sum(chr_dce_coverage[["chr17"]][["GS"]]),sum(chr_dce_coverage[["chr18"]][["GS"]]),
  sum(chr_dce_coverage[["chr19"]][["GS"]]),sum(chr_dce_coverage[["chr20"]][["GS"]]),sum(chr_dce_coverage[["chr21"]][["GS"]]),
  sum(chr_dce_coverage[["chr22"]][["GS"]]),sum(chr_dce_coverage[["chrX"]][["GS"]]),sum(chr_dce_coverage[["chrY"]][["GS"]]))



colnames(dce_sizes_GS) <- c("chr1","chr2","chr3",
                            "chr4","chr5","chr6",
                            "chr7","chr8","chr9",
                            "chr10","chr11","chr12",
                            "chr13","chr14","chr15",
                            "chr16","chr17","chr18",
                            "chr19","chr20","chr21",
                            "chr22","chrX","chrY")
row.names(dce_sizes_GS) <- "sum of dce sizes per chr"



dce_sizes_per_chr <- rbind.data.frame(
  dce_sizes_NP,
  dce_sizes_CIN,
  dce_sizes_EBV,
  dce_sizes_MSI,
  dce_sizes_GS,
  chromosomes_hg38p14$...3
)

row.names(dce_sizes_per_chr) <- c("NP", "CIN", "EBV", "MSI", "GS", "chromosome_length")

dce_colsums_lenghts <- rbind.data.frame(
  (colSums(dce_sizes_per_chr)/chromosomes_hg38p14$...3)*100
) 

colnames(dce_colsums_lenghts) <- c("chr1","chr2","chr3",
                                   "chr4","chr5","chr6",
                                   "chr7","chr8","chr9",
                                   "chr10","chr11","chr12",
                                   "chr13","chr14","chr15",
                                   "chr16","chr17","chr18",
                                   "chr19","chr20","chr21",
                                   "chr22","chrX","chrY")

dce_colsums_lenghts <- rbind.data.frame(dce_colsums_lenghts,
                                        c("chr1","chr2","chr3",
                                          "chr4","chr5","chr6",
                                          "chr7","chr8","chr9",
                                          "chr10","chr11","chr12",
                                          "chr13","chr14","chr15",
                                          "chr16","chr17","chr18",
                                          "chr19","chr20","chr21",
                                          "chr22","chrX","chrY"))

row.names(dce_colsums_lenghts) <- c("Percentage","chromosome")

dce_colsums_lenghts <- t()

rm(dce_colsums_lenghts)


#calculate the percentage of DCEs covering each chromosome
perc_dce_sizes_per_chr <- data.frame(
  rbind((dce_sizes_NP/chromosomes_hg38p14$...3)*100,
        (dce_sizes_CIN/chromosomes_hg38p14$...3)*100,
        (dce_sizes_EBV/chromosomes_hg38p14$...3)*100,
        (dce_sizes_MSI/chromosomes_hg38p14$...3)*100,
        (dce_sizes_GS/chromosomes_hg38p14$...3)*100
  )
)

row.names(perc_dce_sizes_per_chr) <- c("NP", "CIN", "EBV", "MSI", "GS")

perc_dce_sizes_per_chr <- as.data.frame(t(perc_dce_sizes_per_chr))
#perc_dce_sizes_per_chr <-as.matrix(perc_dce_sizes_per_chr)
#pheatmap(t(perc_dce_sizes_per_chr))

pheatmap(perc_dce_sizes_per_chr,
         main = "Chromosome coverage (%) by DCEs in each gastric cancer subtype")



