

library(tidyr)
library(stringr)
library(seqinr)
library(Biostrings)
library(ggplot2)
library(Rsubread)
library(DESeq2)
library(utils)
library(pheatmap)
library(ggplot2)
library(calibrate)
library(Biobase)
library(BiocGenerics)
library(dplyr)
library(tidyverse)
library(ggforce)
library(ggfortify)
library(AnnotationDbi)
library(plyranges)
library(readxl)

#Read in counts data files
setwd("/Users/qkrueger/Documents/terminal_env/vcvd_nvec/All Genes + GO terms/qak")

pr <- readDNAStringSet("~/Documents/terminal_env/reference_genomes/nvec_uk/GCF_932526225.1/rna.fna")
pr_trans <- data.frame(full=names(pr))
remove(pr)
rexp <- "^(\\w+)\\s?(.*)$"

pr_trans$gene <- str_extract(string=pr_trans$full, pattern = "(?<=\\()([^()]*?)(?=\\)[^()]*$)")                  
pr_trans$name <- word(pr_trans$full, 1) 
#write.table(pr_trans,"genetable.txt", sep = "\t")


###############split for population types##################
countdata <- read.csv("counts_table.csv",row.names = 1)
condition <- data.frame(condition=(substr(colnames(countdata),1,nchar(colnames(countdata))-1)))
colData <- data.frame(row.names = colnames(countdata), condition)

FL <- countdata[,grep("FL", colnames(countdata), value = T)]
condition <- substr(colnames(FL),1,nchar(colnames(FL))-1)
colData <- data.frame(row.names = colnames(FL), condition)
countdata <- as.matrix(FL)

NS <- countdata[,grep("NS", colnames(countdata), value = T)]
condition <- substr(colnames(NS),1,nchar(colnames(NS))-1)
colData <- data.frame(row.names = colnames(NS), condition)
countdata <- data.matrix(NS)

NC <- countdata[,grep("NC", colnames(countdata), value = T)]
condition <- substr(colnames(NC),1,nchar(colnames(NC))-1)
colData <- data.frame(row.names = colnames(NC), condition)
countdata <- as.matrix(NC)


##################continue#####################
countdata <- countdata[rowSums(countdata) >= 20,]
condition
colData

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = colData, design = ~condition)
rld <- vst(dds)
rld <- DESeqTransform(rld)
dds <- DESeq(dds)

dists <- dist(t(assay(rld)))


write.csv(as.matrix(dists),"NV_VCVD_dist.csv")

plot(hclust(dists))
pdf("clust_all.pdf", height = 4, width = 8)
plot(hclust(dists))
dev.off()

plotPCA(rld, intgroup = "condition",ntop=nrow(countdata))

pdf("PCA_all.pdf", height = 4, width = 8)
plotPCA(rld, intgroup = "condition",ntop=nrow(countdata))
dev.off()

svg("PCA_all_NC.svg", height = 4, width = 8)
plotPCA(rld, intgroup = "condition",ntop=nrow(countdata))
dev.off()

pca <- plotPCA(rld, intgroup = "condition", 
               ntop=nrow(countdata), returnData = T)
pca
write.csv(pca,"pca_NC.csv")


#immune heatmap
library(circlize)
library(ComplexHeatmap)
immune <- read_excel("list_of_56_immune_genes_nve_with_uk.xlsx")
immune <- rename(immune,"Gene"="UK_gene_model")

#############INDIVIDUAL###################
matrix <- left_join(immune,rld)
matrix <- dist(t()) 

matrix <- assay(rld)
matrix_scaled <- data.frame(t(apply(matrix, 1, scale)))
colnames(matrix_scaled) <- colnames(matrix)
matrix_scaled$Gene <- rownames(matrix_scaled)


matrix_scaled <- na.omit(left_join(immune,matrix_scaled))


#######AVERAGE##########
matrix <- assay(rld)
matrix_scaled <- data.frame(t(apply(matrix, 1, scale)))
colnames(matrix_scaled) <- colnames(matrix)
matrix_scaled$Gene <- rownames(matrix_scaled)
matrix_scaled <- na.omit(matrix_scaled)
tmat <- matrix_scaled %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Value") %>%
  mutate(Group = substr(Sample, 1, 4))
matrix_scaled_means <- tmat %>%
  group_by(Gene, Group) %>%
  summarise(Mean = mean(Value, na.rm = TRUE)) %>%
  pivot_wider(names_from = Group, values_from = Mean)
matrix_scaled <- na.omit(left_join(immune,matrix_scaled_means))
matrix_scaled <- matrix_scaled[order(matrix_scaled$Pathway), ]





immdists <- dist(t(matrix_scaled[5:ncol(matrix_scaled)]))
plot(hclust(immdists))
pdf("clust_immune.pdf", height = 4, width = 8)
plot(hclust(immdists))
dev.off()

pathway_colors <- data.frame(Pathway=c(
  "Apoptosis",
  "cGAS-STING",
  "Interferon-like",
  "RNA Interference",
  "Antibacterial",
  "Other"           
), Color = c(
  "#E41A1C",
  "#377EB8",
  "#4DAF4A",
  "#984EA3",
  "#FF7F00",
  "#999999"
))

matrix_scaled <- left_join(pathway_colors,matrix_scaled)

########################

range(matrix_scaled[ , -1:-5])

c_zscore <- colorRamp2(c(min(matrix_scaled[,1:ncol(matrix_scaled)-1]),
                         0,
                         max(matrix_scaled[,1:ncol(matrix_scaled)-1])), 
                       c("#33CCFF", 
                         "#363636", 
                         "#FFFF00"))

#Immune Heatmap
heat1 <- Heatmap(matrix_scaled[,-1:-5], 
                 col = c_zscore,
                 cluster_rows = T,
                 row_labels = matrix_scaled$Genes_Name,
                 row_names_gp = gpar(col = matrix_scaled$Color, fontsize = 10),
                 column_labels = colnames(matrix_scaled[,-1:-5]), name = "Z-score",
                 column_names_centered = TRUE,
                 column_order = c(1,4,7,3,6,9,2,5,8),
                 row_names_max_width = unit(6, "cm"),
                 show_heatmap_legend = TRUE,
                 column_names_rot = 90,
                 )
heat1
row_ha <- rowAnnotation(
  Group = matrix_scaled$Pathway,
  col = list(Group = c(matrix_scaled$Pathway, matrix_scaled$Color))
)


pdf("heatmap_immune.pdf", width = 6, height = 12)
heat1
dev.off()

############SEPARATE Immune HEATMAPS###################

row_height_unit = unit(8, "mm")
column_width = unit(((1.75*row_height_unit)*length(colnames(matrix_scaled[,-1:-5]))),"mm")

for (f in 1:length(unique(matrix_scaled$Pathway))) {

  dynamic_height = row_height_unit * nrow(matrix_scaled[grep(pathway_colors[f,1], matrix_scaled$Pathway), -1:-5])
  
  heat1 <- Heatmap(matrix_scaled[grep(pathway_colors[f,1], matrix_scaled$Pathway), -1:-5],
                   heatmap_height = dynamic_height,
                   heatmap_width = column_width,
                   col = c_zscore,
                   cluster_rows = T,
                   row_labels = matrix_scaled[grep(pathway_colors[f,1], matrix_scaled$Pathway), 3],
                   row_names_gp = gpar(col = matrix_scaled[grep(pathway_colors[f,1], matrix_scaled$Pathway), 2],
                                       fontsize = 10),
                   column_labels = colnames(matrix_scaled[,-1:-5]), 
                   name = "Z-score",
                   column_names_centered = TRUE,
                   column_order = c(1,4,7,3,6,9,2,5,8),
                   row_names_max_width = unit(6, "cm"),
                   show_heatmap_legend = T,
                   column_names_rot = 90,)
  heat1

  pdf(paste("heatmap_immune_meanZ_sort_",
            unique(matrix_scaled[grep(pathway_colors[f,1], matrix_scaled$Pathway), 1]),
            ".pdf",sep = ""), width = 6, height = 11)
  print(heat1)
  dev.off()
  
}

pca_res <- prcomp(t(matrix_scaled[,-1:-5]), scale. = F) # convert to pca object
str(pca_res)
summary(pca_res)

tcountdata <- t(countdata)

#otherwise continue
cond <- rownames(colData)
df <- NA
df$PC1 <- as.numeric(pca_res$x[,1]) / (pca_res$sdev[1] * sqrt(nrow(tcountdata)))
df$PC2 <- as.numeric(pca_res$x[,2]) / (pca_res$sdev[2] * sqrt(nrow(tcountdata)))

df <- as.data.frame(df)
df <- df[-1]

df <- cbind(df,cond)
df$Condition <- substr(df$cond,1,nchar(df$cond)-1)
df$Genotype <- substr(df$cond,1,2)
df$Bacteria <- substr(df$cond,3,nchar(df$cond)-1)
propvar <-summary(pca_res)
pervar <- propvar$importance

write.csv(df,"pca_immune.csv")

pgg1 <- ggplot(data = df, aes(x=PC1, y=PC2)) + 
  geom_encircle(aes(group = Condition, fill = Genotype, color = Genotype), 
                s_shape=1, expand=0, alpha = 0.65, size = 3) +
  #geom_encircle(aes(group = Genotype), 
  #              s_shape=1, expand=0, alpha = 0.15, size = 3) +
  geom_point(aes(shape = Bacteria), size = 3, alpha = 0.95) +
  #scale_fill_manual(values=c("blue","yellow","red"))+
  #scale_color_manual(values=c("blue","yellow","red"), aesthetics = "color")+
  theme_bw() + #labs(fill="Time", color = "Time", shape = "Temp") + 
  xlab(paste("PC1 - ",round((pervar[2,1]*100), 2),"%", sep = "")) +
  ylab(paste("PC2 - ",round((pervar[2,2]*100), 2),"%", sep = "")) +
  ggtitle("4H Thermal Exposure") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=22,face="bold"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22,face = "bold"),
        plot.title = element_text(color="black", size=24),
        #legend.position = "none"
        )

pgg1


pdf("PCA_immune_all.pdf", width = 12, height = 8)  
pgg1
dev.off()


#pairwise
comploop <- read.csv("comps.csv")
comploop

gene <- read.delim("genetable.txt")
gr <- data.frame(read_gff("/Users/qkrueger/Documents/terminal_env/reference_genomes/nvec_uk/GCF_932526225.1/genomic.gff")) %>% 
  dplyr::select(Name, product)

gr <- na.omit(gr)
gene <- left_join(gene,gr,by="Name")

rm(gr)

f <- 1
dl <- as.list(NA)

degs <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(degs)<-c("cond1","cond2",
                  "1.5all","1.5up","1.5down")

for (f in 1:nrow(comploop)){
  
  countdata <- read.csv("counts_table.csv",row.names = 1)
  countdata <- data.frame(countdata[,grep(comploop[f,3],
                                       colnames(countdata))])
  condition <- data.frame(condition=(substr(colnames(countdata),1,nchar(colnames(countdata))-1)))
  colData <- data.frame(row.names = colnames(countdata), condition)
  countdata <- countdata[rowSums(countdata) >= 10,]
  dds <- DESeqDataSetFromMatrix(countData = countdata, colData = colData, design = ~condition)
  rld <- vst(dds)
  rld <- DESeqTransform(rld)
  dds <- DESeq(dds)
  
  comp <- as.data.frame(results(dds, 
                            contrast = c("condition",comploop[f,2],comploop[f,1])))
  comp$padj[is.na(comp$padj)] <- 1
  
  comp <- na.omit(comp)
  comp$gene <- rownames(comp)
  countdata$gene <- rownames(countdata)
  comp <- left_join(comp,gene)
  comp <- left_join(comp,countdata)
  comp <- comp %>%
    distinct(gene,.keep_all = T)
  
  sigcomp <- comp[comp$padj <0.05,]
  sigcomp <- sigcomp[(abs(sigcomp$log2FoldChange) > 1.5),]
  #write.csv(sigcomp, file = paste(comploop[f,1], "_", comploop[f,2], ".csv", sep=""))
  #write.csv(comp, file = paste(comploop[f,1],"_", comploop[f,2], "_all.csv", sep = ""))
  
  comp$diffexpressed <- "Not Significant"
  comp$diffexpressed[comp$log2FoldChange < -1.5 & comp$padj < 0.05] <- "Down-regulated"
  comp$diffexpressed[comp$log2FoldChange > 1.5 & comp$padj < 0.05] <- "Up-regulated"
  
  table(comp$diffexpressed)
  comp$diffexpressed <- factor(comp$diffexpressed, 
                               levels = c("Up-regulated", "Down-regulated", "Not Significant"))
  
  volc <- ggplot(data = comp, 
    aes(x=log2FoldChange, y=-log10(padj), col = diffexpressed)) + 
    geom_vline(xintercept = c(-1.5, 1.5), col = "red") + 
    geom_hline(yintercept = -log10(0.05), col = "red") + 
    geom_point() + theme_minimal() + 
    scale_color_manual(values = c("red", "blue", "black")) + 
    ggtitle(paste(comploop[f,1]," ", comploop[f,2], sep = "")) + 
    guides(color = guide_legend(title = "Differential Expression"))
  
  
  #pdf(paste(comploop[f,1],"_", comploop[f,2], ".pdf", sep = ""), 
  #    height = 4, width = 8)
  #print(volc)
  #dev.off()
  
  up15 <- sigcomp[sigcomp$log2FoldChange > 1.5,]
  down15 <- sigcomp[sigcomp$log2FoldChange < -1.5,]

  degs[f,] <- cbind(comploop[f,1], comploop[f,2],
                    nrow(sigcomp),
                    nrow(up15),
                    nrow(down15))

  up15 <- data.frame(sigcomp[sigcomp$log2FoldChange > 1.5,7])
  names(up15) <- paste(comploop[f,1],comploop[f,2],"up1.5",sep = "_")
  down15 <- data.frame(sigcomp[sigcomp$log2FoldChange < -1.5,7])
  names(down15) <- paste(comploop[f,1],comploop[f,2],"down1.5",sep = "_")
  comp15 <- data.frame(sigcomp[,7])
  names(comp15) <- paste(comploop[f,1],comploop[f,2],"1.5",sep = "_")
  
  dl[[f]] <- c(as.list(comp15),
               as.list(up15),
               as.list(down15))
  
}

write.csv(degs,"DEGs_VCVD_4h_QAK.csv")


capture.output(print(dl), file = "upset_DEGs.csv")
#write.csv(dl,"upset_DEGs.csv")

alldegs <- do.call(rbind,dl)
alldegs <- data.frame(Gene = unique(unlist(alldegs, use.names = F)))

matrix_scaled <- left_join(alldegs,matrix_scaled)

#FULL Heatmap
heat1 <- Heatmap(as.matrix(matrix_scaled[,2:ncol(matrix_scaled)]), 
                 col = c_zscore,
                 cluster_rows = T, show_row_names = F,
                 column_labels = colnames(matrix_scaled[2:ncol(matrix_scaled)]), 
                 name = "Z-score",
                 column_names_centered = TRUE,
                 column_order = rev(colnames(matrix_scaled[,2:ncol(matrix_scaled)])),
                 cluster_columns = F, cluster_column_slices = F,
                 show_heatmap_legend = TRUE,
                 column_names_rot = 90, use_raster = F
                
)


pdf("heatmap_alldegs.pdf", width = 12, height = 24)
heat1
dev.off()

library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(purrr)
library(reshape)
library(ggVennDiagram)

all <- map(dl, 1)
up <- map(dl, 2)
down <- map(dl, 3)

# Retrieve the name of the first element from each sub-list
alln<- sapply(dl, function(sublist) names(sublist)[1])
names(all) <- alln
all<- all[grep(".*NT", names(all))]
cm <- make_comb_mat(all,mode="distinct")
Reduce(intersect, all)

upn<- sapply(dl, function(sublist) names(sublist)[2])
names(up) <- upn
up<- up[grep(".*NT", names(up))]
cm <- make_comb_mat(up,mode="distinct")


downn<- sapply(dl, function(sublist) names(sublist)[3])
names(down) <- downn
down<- down[grep(".*NT", names(down))]
cm <- make_comb_mat(down,mode="distinct")
Reduce(intersect, down)




cm<-cm[comb_size(cm) >= 10]


ss = set_size(cm)
cs = comb_size(cm)

pdf("upsetupall.pdf", height = 4, width = 8)
ht = UpSet(cm, 
           #set_order = order(ss),
           set_order = order(rownames(cm), decreasing = T),
           comb_order = order(comb_degree(cm), -cs),
           top_annotation = HeatmapAnnotation(
             "Gene Intersections" = anno_barplot(cs, 
                                                 ylim = c(0, 600),
                                                 border = FALSE, 
                                                 gp = gpar(fill = "red"), 
                                                 height = unit(5, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           right_annotation = rowAnnotation(
             "Total Genes" = anno_barplot(ss, 
                                          baseline = 0,
                                          add_numbers = TRUE,
                                          ylim = c(0,900),
                                          #axis_param = list(
                                          #  at = c(0, 100, 200, 300),
                                          #  labels = c(0, 100,200,300),
                                          #  labels_rot = 0),
                                          border = FALSE, 
                                          gp = gpar(fill = "red"), 
                                          width = unit(1.5, "cm")
             )
             
           ), 
           left_annotation = rowAnnotation(
             set_name = anno_text(set_name(cm), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(cm)) + unit(4, "mm"))
           ),
           show_row_names = FALSE)

ht = draw(ht)
od = column_order(ht)
yeet <- decorate_annotation("Gene Intersections", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(15, "pt"), 
            default.units = "native", just = "center", 
            gp = gpar(fontsize = 12, col = "red"), rot = 0)
})
dev.off()


#Check filtering counts

setwd("~/Documents/terminal_env/vcvd_nvec/All Genes + GO terms/")

files <- list.files(".", pattern = "*1.5.csv", full.names = F)

counts <- read.csv("counts_table.csv")
counts <- rename(counts,Transcript_ID=X)
rownames(counts) <- counts$Transcript_ID
counts <- counts [,-1]

keep1 <- counts[rowSums(counts) >= 1*(ncol(counts)),]
keep3 <- counts[rowSums(counts) >= 3*(ncol(counts)),]
nrow(counts)
nrow(keep1)
nrow(keep3)


filesdf <- data.frame(comp=files)

filesdf <- data.frame(comp=t(rbind(
  filesdf[grep(paste(c("FL.*FL","NS.*NS","NC.*NC"),collapse="|"),filesdf$comp),]
)))

fname <- filesdf %>% mutate(comp = str_replace(comp, "\\_", "|")) %>% 
  separate(comp, into = c("F1", "F2"), sep = "\\|")
fname$F2 <- str_extract(string=fname$F2, pattern = ".*(?=\\_1)")

f <- 1
ftg <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(ftg)<-c("cond1","cond2",
                 "1.5_pre-filter","1-sample","3-sample")

for (f in 1:nrow(filesdf)){
  
  comp <- read.csv(filesdf[f,])
  fname[f,1]
  fname[f,2]
  
  comp1 <- comp[rowSums(comp[,13:ncol(comp)]>=(1*ncol(comp)-12)),]
  comp3 <- comp[rowSums(comp[,13:ncol(comp)]>=(3*ncol(comp)-12)),]
  
  
  ftg[f,] <- cbind(fname[f,1], fname[f,2],
                   nrow(comp),
                   nrow(comp1),
                   nrow(comp3))
  
  
}

write.csv(ftg,"postadmera_1.5_filter_VCVD_4h.csv")

############ASV COUNTS################

setwd("~/Documents/terminal_env/vcvd_nvec/")
laura <- read_xlsx("laura_blast_table.xlsx")

laura %>%
  group_by(Location) %>%
  summarize(sum.average = sum(average))


###########Immune plot#############
library(data.table)
immune <- read_excel("list_of_56_immune_genes_nve_with_uk.xlsx")
immune <- rename(immune,"gene"="UK_gene_model")

files <- list.files(".", pattern = "*all.csv", full.names = F)
files <- files[grep("NT",files)]
files <- files[-grep("NT.*NT",files)]
files
imm <- as.list(NA)
for (f in 1:length(files)){
  
  imm[[f]] <- read.csv(files[f])
  imm[[f]] <- left_join(immune,imm[[f]])
  imm[[f]] <- imm[[f]][1:11]

}
names(imm) <- files

Limm <- rbindlist(imm, idcol = "Condition")
Limm$Condition <- substr(Limm$Condition,1,nchar(Limm$Condition)-8)

dodge = 0.5

ggimmune <- ggplot(Limm, aes(x=Genes_Name, y=`log2FoldChange`, colour=`Condition`)) +
  #facet_grid(cols = vars(immunetemp$`Condition 2`)) +
  xlab(Limm$Short) +
  #scale_y_continuous(breaks=seq(-6,6,by=0.5)) +
  #geom_abline(data = hlines, aes(intercept = s, slope = 0)) +
  geom_hline(yintercept = rep(1.5,nrow(Limm)), col = "#919191") +
  geom_hline(yintercept = rep(-1.5,nrow(Limm)), col = "#919191") +
  geom_hline(yintercept = rep(0,nrow(Limm)), col = "#D3D3D3") +
  #scale_color_manual(values = c("#000CFF","yellow","red"))+
  geom_linerange(aes(ymin=`log2FoldChange`-lfcSE, ymax=`log2FoldChange`+lfcSE), position = position_dodge(width = dodge)) +
  geom_point(position=position_dodge(width = dodge), stat="identity", size = 2) +
  ggtitle("Time")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

ggimmune



pdf("immune.pdf", width = 24, height = 6)
ggimmune
dev.off()

write.csv(Limm,"immune_l2fc.csv")

