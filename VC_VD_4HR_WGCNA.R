
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(BiocManager)
BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
BiocManager::install("Biostrings")
library(Biostrings)
BiocManager::install("WGCNA")
library(WGCNA)
library(flashClust)
library(dynamicTreeCut)
library(SummarizedExperiment)
library(magrittr)
library(purrr)

#Read in counts data files
setwd("~/work/projects/Thesis/vcvd_4hr/")
files <- list.files(
  path = "raw_counts/",
  pattern = NULL,
  full.names = TRUE,
  recursive = FALSE
)


#Read files into a list
install.packages("ragg", type = "source")
install.packages(c("tidyverse", "magrittr", "purrr"))
library(tidyverse)
library(magrittr)
library(purrr)
library(tidyr)
file_list <- files %>%
  map(~ {
    # read the file
    df <- read.delim(.x, header = FALSE, sep = "\t")
    # make a safe column name out of just the filename (no path, no .txt)
    sample_name <- tools::file_path_sans_ext(basename(.x))
    # set names: first column always "Transcript_ID", second is the sample name
    colnames(df) <- c("Transcript_ID", sample_name)
    df
  })

# 2. Read them into a named list of tibbles
dfs <- files %>% 
  set_names(~ tools::file_path_sans_ext(basename(.x))) %>% 
  map(~ read_tsv(.x, col_names = c("Transcript_ID", basename(.x))))
# 3. Merge by Transcript_ID
merged_counts <- purrr::reduce(dfs, dplyr::full_join, by = "Transcript_ID")

#Filter out tRNA
merged_counts <- merged_counts %>% 
  filter(!str_detect(Transcript_ID, regex("tRNA", ignore_case = TRUE)))

#Move to a matrix with proper rownames
library(tibble)

counts_matrix <- merged_counts %>% 
  column_to_rownames("Transcript_ID") %>% 
  as.matrix()

#Auto-generate metadata from sample names
metadata <- tibble(sample = colnames(counts_matrix)) %>% 
  mutate(
    Condition = sample %>% 
      str_extract("^[A-Z]{2,4}") %>%      # e.g. “FLNT”, “NCVD”, etc.
      factor(levels = unique(.))
  ) %>% 
  column_to_rownames("sample")

# move the Transcript_ID column into row-names
merged_counts <- merged_counts %>% 
  column_to_rownames(var = "Transcript_ID")

#Final check, should be true
all(rownames(metadata) == colnames(counts_matrix))

#Keep only FL samples
merged_fl <- merged_counts[,grep("FL", colnames(merged_counts))]
metadata_fl <- data.frame(Condition = colnames(merged_fl))
rownames(metadata_fl) <- metadata_fl$Condition
metadata_fl$Condition <- substr(metadata_fl$Condition, 1, nchar(metadata_fl$Condition)-1)

#Run deg and transform just FL data
dds_fl <- DESeqDataSetFromMatrix(countData = merged_fl,
                                 colData = metadata_fl,
                                 design = ~ Condition)
rld_fl <- vst(dds_fl, fitType = "mean")
fl_pca <- plotPCA(rld_fl, intgroup = "Condition", returnData = TRUE)
fl_pca_object <- prcomp(t(assay(rld_fl)))
fl_percentvar <- fl_pca_object$sdev^2 / sum(fl_pca_object$sdev^2) * 100
custom_fl_pca_colors <- c("FLNT" = "firebrick1", "FLVC" = "firebrick4", "FLVD" = "orange")
custom_fl_pca_shapes <- c("FLNT" = 15, "FLVC" = 15, "FLVD" = 15)
fl_pca <- ggplot(fl_pca, aes(x = PC1, y = PC2, color = Condition, shape = Condition)) +
  geom_point(size = 4, alpha = 1) +  # Adjust size & transparency
  scale_color_manual(values = custom_fl_pca_colors ) +  # Apply custom colors
  scale_shape_manual(values = custom_fl_pca_shapes) +  # Apply custom shapes
  labs(title = "Florida Genotypes",
       x = paste0("PC1 (", round(fl_percentvar[1], 2), "%)"),
       y = paste0("PC2 (", round(fl_percentvar[2], 2), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("FloridaPCAPlot.pdf", plot = fl_pca, width = 6, height = 5, dpi = 300)

#Deseq and pca for NC
merged_nc <- merged_counts[,grep("NC", colnames(merged_counts))]
metadata_nc <- data.frame(Condition = colnames(merged_nc))
rownames(metadata_nc) <- metadata_nc$Condition
metadata_nc$Condition <- substr(metadata_nc$Condition, 1, nchar(metadata_nc$Condition)-1)

dds_nc <- DESeqDataSetFromMatrix(countData = merged_nc,
                                 colData = metadata_nc,
                                 design = ~ Condition)
rld_nc <- vst(dds_nc, fitType = "mean")
nc_pca <- plotPCA(rld_nc, intgroup = "Condition", returnData = TRUE)
nc_pca_object <- prcomp(t(assay(rld_nc)))
nc_percentvar <- nc_pca_object$sdev^2 / sum(nc_pca_object$sdev^2) * 100
custom_nc_pca_colors <- c("NCNT" = "olivedrab3", "NCVC" = "darkolivegreen", "NCVD" = "seagreen3")
custom_nc_pca_shapes <- c("NCNT" = 16, "NCVC" = 16, "NCVD" = 16)
nc_pca <- ggplot(nc_pca, aes(x = PC1, y = PC2, color = Condition, shape = Condition)) +
  geom_point(size = 4, alpha = 1) +  # Adjust size & transparency
  scale_color_manual(values = custom_nc_pca_colors ) +  # Apply custom colors
  scale_shape_manual(values = custom_nc_pca_shapes) +  # Apply custom shapes
  labs(title = "North Carolina Genotypes",
       x = paste0("PC1 (", round(nc_percentvar[1], 2), "%)"),
       y = paste0("PC2 (", round(nc_percentvar[2], 2), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("NC_PCAPlot.pdf", plot = nc_pca, width = 6, height = 5, dpi = 300)

#Deseq and pca for NS
merged_ns <- merged_counts[,grep("NS", colnames(merged_counts))]
metadata_ns <- data.frame(Condition = colnames(merged_ns))
rownames(metadata_ns) <- metadata_ns$Condition
metadata_ns$Condition <- substr(metadata_ns$Condition, 1, nchar(metadata_ns$Condition)-1)

dds_ns <- DESeqDataSetFromMatrix(countData = merged_ns,
                                 colData = metadata_ns,
                                 design = ~ Condition)
rld_ns <- vst(dds_ns, fitType = "mean")
ns_pcas <- plotPCA(rld_ns, intgroup = "Condition", returnData = TRUE)
ns_pca_object <- prcomp(t(assay(rld_ns)))
ns_percentvar <- ns_pca_object$sdev^2 / sum(ns_pca_object$sdev^2) * 100
custom_ns_pca_colors <- c("NSNT" = "dodgerblue", "NSVC" = "turquoise2", "NSVD" = "navyblue")
custom_ns_pca_shapes <- c("NSNT" = 17, "NSVC" = 17, "NSVD" = 17)
ns_pca <- ggplot(ns_pcas, aes(x = PC1, y = PC2, color = Condition, shape = Condition)) +
  geom_point(size = 4, alpha = 1) +  # Adjust size & transparency
  scale_color_manual(values = custom_ns_pca_colors ) +  # Apply custom colors
  scale_shape_manual(values = custom_ns_pca_shapes) +  # Apply custom shapes
  labs(title = "Nova Scotia Genotypes",
       x = paste0("PC1 (", round(ns_percentvar[1], 2), "%)"),
       y = paste0("PC2 (", round(ns_percentvar[2], 2), "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("NS_PCAPlot.pdf", plot = ns_pca, width = 6, height = 5, dpi = 300)

#Create deseq object
library(BiocManager)
BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = metadata,
                              design = ~ Condition)
rld <- vst(dds, fitType = "mean")

# the three “ntop” values you want
ntop_values <- c(500, 1000, 10000)

# custom aesthetics
colors <- c(
  "FLNT"="firebrick1","FLVC"="firebrick4","FLVD"="orange",
  "NCNT"="olivedrab3","NCVC"="darkolivegreen","NCVD"="seagreen3",
  "NSNT"="dodgerblue","NSVC"="turquoise2","NSVD"="navyblue"
)
shapes <- c(
  "FLNT"=15,"FLVC"=15,"FLVD"=15,
  "NCNT"=16,"NCVC"=16,"NCVD"=16,
  "NSNT"=17,"NSVC"=17,"NSVD"=17
)

output_dir <- "figures/pcas/"


for (ntop in ntop_values) {
  # extract the PCA data
  pca_data <- plotPCA(rld,
                      intgroup="Condition",
                      ntop=ntop,
                      returnData=TRUE)
  
  # pull out % variances
  var_pct <- attr(pca_data, "percentVar")
  
  # make the ggplot
  p <- ggplot(pca_data, aes(PC1, PC2, color=Condition, shape=Condition)) +
    geom_point(size=2) +
    scale_color_manual(values=colors) +
    scale_shape_manual(values=shapes) +
    labs(
      title   = paste0("PCA (top ", ntop, " most variable genes)"),
      x       = sprintf("PC1 (%.1f%%)", var_pct[1]*100),
      y       = sprintf("PC2 (%.1f%%)", var_pct[2]*100)
    ) +
    theme_minimal(base_size=14) +
    theme(plot.title=element_text(hjust=0.5))
  
  # save
  ggsave(filename = paste0("PCA_ntop", ntop, ".pdf"),
         plot     = p,
         path     = output_dir, 
         width    = 6, height = 5, dpi = 300)
}
# a function to make one PCA plot for a given ntop
make_pca_plot <- function(ntop){
  # get the PCA data
  pca_data <- plotPCA(rld, intgroup="Condition", ntop=ntop, returnData=TRUE)
  var_pct  <- attr(pca_data, "percentVar")
  
  # build the ggplot
  ggplot(pca_data, aes(PC1, PC2, color=Condition, shape=Condition)) +
    geom_point(size=2) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(
      title = paste0("Top ", ntop, " most variable genes"),
      x     = sprintf("PC1 (%.1f%%)", var_pct[1]*100),
      y     = sprintf("PC2 (%.1f%%)", var_pct[2]*100)
    ) +
    theme_minimal(base_size=14) +
    theme(plot.title = element_text(hjust=0.5))
}

# generate all three plots
pca_plots <- map(ntop_values, make_pca_plot)

# stitch into a 1×3 collage
collage <- wrap_plots(pca_plots, ncol=3) +
  plot_annotation(title = "PCA comparison at different ntop thresholds")

# save
ggsave("figures/PCA_collage.pdf", collage,
       width=16, height=5, dpi=300)

#Normalize data and run deseq analysis
dds <- DESeq(dds)

BiocManager::install("Biostrings")
library(Biostrings)
library(stringr)
library(tibble)

# 1) Read fasta and get headers
uk_genome  <- readDNAStringSet("~/work/projects/Thesis/uk_genome/ukgenome/GCF_932526225.1/cds_from_genomic.fna")
uk_headers <- names(uk_genome)

# 2) Pull out “gene=” and “protein=” with str_extract()
transcript_ids       <- str_extract(uk_headers, "(?<=gene=)[^\\]]+")
protein_annotations  <- str_extract(uk_headers, "(?<=protein=)[^\\]]+")

# 3) Assemble into a tibble (no loop, no growing vectors)
uk_genes_protein <- tibble(
  Transcript_ID       = transcript_ids,
  Protein_Annotation  = protein_annotations
)

# 4) (Optional) Replace any NA strings with actual NA
uk_genes_protein <- mutate(uk_genes_protein,
                           Transcript_ID      = if_else(is.na(Transcript_ID), NA_character_, Transcript_ID),
                           Protein_Annotation = if_else(is.na(Protein_Annotation), NA_character_, Protein_Annotation)
)

uk_genes_protein


#Get all pairwise comparisons based on condition
library(BiocParallel)

# 1. Define exactly the comparisons you need
desired <- c(
  "FLVCvsFLNT","NCVCvsNCNT","NSVCvsNSNT",
  "FLVCvsFLVD","NCVCvsNCVD","NSVCvsNSVD",
  "FLVDvsFLNT","NCVDvsNCNT","NSVDvsNSNT"
)
pair_list <- strsplit(desired, "vs")
pairwise_results <- setNames(
  lapply(pair_list, function(pair) {
    results(dds, contrast = c("Condition", pair[1], pair[2]))
  }),
  desired
)

# 2. Register a multicore backend (use one fewer than total cores to leave 
#    one free for the system)
register(MulticoreParam(workers = parallel::detectCores() - 1))

# 3. Split each “A vs B” string into c("A","B")
pair_list <- strsplit(desired, "vs")

# 4. Run results() in parallel
pairwise_results <- bplapply(pair_list, function(pair, dds_obj) {
  results(dds_obj, contrast = c("Condition", pair[1], pair[2]))
}, dds, BPPARAM = bpparam())

# 5. Name the list elements for easy lookup
names(pairwise_results) <- desired
condition_pairs <- combn(levels(metadata$Condition), 2, simplify = FALSE)
pairwise_results <- list()
for(pair in condition_pairs) {
  res <- results(dds, contrast = c("Condition", pair[1], pair[2]))
  pairwise_results[[paste(pair[1], "vs", pair[2])]] <- res
}
# 1) Define exactly the comparisons you need
desired <- c(
  "FLVCvsFLNT","NCVCvsNCNT","NSVCvsNSNT",
  "FLVCvsFLVD","NCVCvsNCVD","NSVCvsNSVD",
  "FLVDvsFLNT","NCVDvsNCNT","NSVDvsNSNT"
)

# 2) Iterate over them, splitting on “vs” to feed into DESeq2
pairwise_results <- lapply(desired, function(comp) {
  parts <- strsplit(comp, "vs")[[1]]
  # contrast = c("Condition", numerator, denominator)
  results(dds, contrast = c("Condition", parts[1], parts[2]))
})

library(BiocParallel)
register(MulticoreParam(4))  

pairwise_results <- bplapply(desired, function(comp) {
  parts <- strsplit(comp, "vs")[[1]]
  results(dds, contrast = c("Condition", parts[1], parts[2]))
})

# 3) Name your list so you can pull them back out by the same labels
names(pairwise_results) <- desired
#For each comparison make three files, one being all the genes which will have all
#standard deseq information plus product and transcript; the 2nd will be only the 
#differentially expressed features but same information; 3rd will be all genes but#
#just two columns, gene id and logfoldchange

#Read in UK genome annotations
go_terms <- read_delim("~/work/projects/Thesis/uk_genome/ukgenome/UK_GO (1).tsv")

go_terms <- go_terms %>%
  mutate(GENEID = str_trim(GENEID),
         GENEID = paste0("LOC", GENEID))
go_terms <- go_terms %>%
  rename("Transcript_ID" = "GENEID")

# 0) Build a GO-map that collapses multiple terms per gene
go_map <- go_terms %>%
  group_by(GENEID) %>% 
  summarise(
    GOTERMS = paste(unique(GOTERMS), collapse = ";"),
    .groups = "drop"
  )

library(dplyr)
#Loop through all pairwise comparisons and make a volcano plot for each one
pairwise_sel <- pairwise_results[desired]

# thresholds you want
lfc_thresholds <- c(3, 1.5)

# make an output directory
outdir <- "figures/volcanoes/"
if (!dir.exists(outdir)) dir.create(outdir)

for(contrast_name in names(pairwise_sel)) {
  message("Plotting: ", contrast_name)
  res_df <- as.data.frame(pairwise_sel[[contrast_name]]) %>%
    rownames_to_column("Transcript_ID") %>%
    left_join(uk_genes_protein, by = "Transcript_ID") %>%
    mutate(log10_pvalue = -log10(padj))
  
  plot_df <- res_df %>%
    mutate(
      Category = case_when(
        log2FoldChange >=  3 & padj < 0.05 ~ "Strong Up",
        log2FoldChange <= -3 & padj < 0.05 ~ "Strong Down",
        log2FoldChange >=  1.5 & padj < 0.05 ~ "Moderate Up",
        log2FoldChange <= -1.5 & padj < 0.05 ~ "Moderate Down",
        TRUE                                 ~ "Nonsignificant"
      ),
      Category = factor(
        Category,
        levels = c("Strong Up","Moderate Up","Moderate Down","Strong Down","Nonsignificant")
      )
    )
  
  p <- ggplot(plot_df, aes(x = log2FoldChange, y = log10_pvalue, color = Category)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c(
      "Strong Up"       = "red3",
      "Moderate Up"     = "orange",
      "Moderate Down"   = "skyblue",
      "Strong Down"     = "blue3",
      "Nonsignificant"  = "gray"
    )) +
    # draw both cutoffs
    geom_vline(xintercept = c(-1.5,  1.5), linetype = "dotted") +
    geom_vline(xintercept = c(-3.0,  3.0), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title    = contrast_name,
      x        = "Log₂ Fold Change",
      y        = "-Log₁₀ Adjusted P-value",
      color    = "Category"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title    = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  out_file <- file.path(outdir, paste0(contrast_name, "_volcano_bothLFC.pdf"))
  ggsave(out_file, p, width = 6, height = 5, dpi = 300)
}
library(readr)
# Your named list of DESeq2 results
pairwise_sel <- pairwise_results[desired]

# Output directory (will be created if needed)
outdir <- "spreadsheets_new/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for(contrast_name in desired) {
  message("Writing CSVs for: ", contrast_name)
  
  # 1) Base DESeq2 table with your protein annotations
  df_base <- as.data.frame(pairwise_sel[[contrast_name]]) %>%
    rownames_to_column("Transcript_ID") %>%
    left_join(uk_genes_protein, by = "Transcript_ID") %>%
    mutate(log10_pvalue = -log10(padj))
  
  # 2) Attach collapsed GO terms
  df_base <- df_base %>%
    left_join(go_map, by = c("Transcript_ID" = "GENEID"))
  
  # 3) Annotate strength categories
  df_all <- df_base %>%
    mutate(
      Category = case_when(
        log2FoldChange >=  3 & padj < 0.05   ~ "Strong Up",
        log2FoldChange <= -3 & padj < 0.05   ~ "Strong Down",
        log2FoldChange >=  1.5 & padj < 0.05 ~ "Moderate Up",
        log2FoldChange <= -1.5 & padj < 0.05 ~ "Moderate Down",
        TRUE                                  ~ "Nonsignificant"
      )
    )
  
  # 4) Significant = any non-“Nonsignificant”
  df_sig <- df_all %>% 
    filter(Category != "Nonsignificant")
  
  # 5) LFC only
  df_lfc <- df_all %>%
    select(Transcript_ID, log2FoldChange)
  
  # 6) Write CSVs
  write_csv(
    df_all,
    file.path(outdir, paste0(contrast_name, "_all_genes.csv"))
  )
  write_csv(
    df_sig,
    file.path(outdir, paste0(contrast_name, "_significant_genes.csv"))
  )
  write_csv(
    df_lfc,
    file.path(outdir, paste0(contrast_name, "_log2FC_only.csv"))
  )
}
  # Check if the necessary columns exist
  #if(!all(c("log2FoldChange", "padj") %in% colnames(volcano_data))) {
    stop("Missing necessary columns: log2FoldChange or padj in the results data")

  # Add transcript IDs (assuming row names contain them)
  volcano_data$Transcript_ID <- rownames(volcano_data)
  
  #Make csv with all genes + GO terms
  #res_go <- merge(volcano_data, go_terms, by = "Transcript_ID", all.x = TRUE)
  #res_go_simple <- res_go[, c("Transcript_ID", "GOTERMS")]
  #write.csv(res_go_simple, file = paste0(pair, "_all_GO.csv"), row.names = FALSE)
  # Merge with annotation data
  volcano_data <- merge(volcano_data, uk_genes_protein, by = "Transcript_ID", all.x = TRUE)
  
  #Add -log10 pvalue column
  volcano_data$log10_pvalue <- -log10(volcano_data$padj)
  
  # Remove rows with NA in log2FoldChange or log10_pvalue
  volcano_data <- volcano_data[!is.na(volcano_data$log2FoldChange) & !is.na(volcano_data$log10_pvalue), ]
  
  #Create columns for not significant, upragulated, and downregulated
  volcano_data$Category <- "Nonsignificant"  # Default category
  volcano_data$Category[volcano_data$log2FoldChange < -3 & volcano_data$padj < 0.05] <- "Downregulated"
  volcano_data$Category[volcano_data$log2FoldChange > 3 & volcano_data$padj < 0.05] <- "Upregulated"
  
  #Create volcano plots
  p <- ggplot(volcano_data, aes(x = log2FoldChange, y = log10_pvalue, colour = Category)) +
    geom_point(alpha = 0.6) + 
    scale_color_manual(values = c("Nonsignificant" = "gray", 
                                  "Upregulated" = "red", 
                                  "Downregulated" = "blue")) +  # Custom color mapping
    labs(title = paste("Pairwise Comparison:", pair),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = c(-3,3), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
  
  #Create spreadsheets of the comparisons and save
  #significant_genes <- volcano_data[volcano_data$Category != "Nonsignificant", ]
  #write.csv(significant_genes, file = paste0(pair, ".csv"), row.names = FALSE)

  #write .csv volcano data
  #write.csv(volcano_data, file = paste0(pair, "_all.csv"), row.names = FALSE)
  #genes_plus_go <- volcano_data %>%
    #left_join(go_terms, by = "Transcript_ID")
  #genes_plus_go <- genes_plus_go %>%
    #distinct(Transcript_ID, .keep_all = TRUE)
  
  #write.csv(genes_plus_go, file = paste0(pair, "_all_GO.csv"), row.names = FALSE)
  
  #3rd table:gene names, and logfold change from volcano data
  #gene_plus_logfold <- volcano_data %>%
    select(Transcript_ID, log2FoldChange)
  #write_delim(gene_plus_logfold,delim = "\t", file = paste0(pair, "_mwu.tsv"))
  #Save plot to computer
  dir.create("Plots1", showWarnings = FALSE)
  ggsave(paste0("Plots1/", pair, "_volcano_plot.pdf"), plot = p, width = 8, height = 6)

#Make volcano collage
if (!requireNamespace("patchwork", quietly=TRUE)) {
    install.packages("patchwork")
}
library(patchwork)
# 2) Re-generate your nine volcano plots as a named list of ggplot objects
volcano_plots <- lapply(names(pairwise_sel), function(contrast_name) {
  # build the data
  df <- as.data.frame(pairwise_sel[[contrast_name]]) %>%
    rownames_to_column("Transcript_ID") %>%
    left_join(uk_genes_protein, by="Transcript_ID") %>%
    mutate(
      log10_pvalue = -log10(padj),
      Category = case_when(
        log2FoldChange >=  3 & padj < 0.05   ~ "Strong Up",
        log2FoldChange <= -3 & padj < 0.05   ~ "Strong Down",
        log2FoldChange >=  1.5 & padj < 0.05 ~ "Moderate Up",
        log2FoldChange <= -1.5 & padj < 0.05 ~ "Moderate Down",
        TRUE                                  ~ "Nonsignificant"
      )
    )
  
  # create the ggplot
  p <- ggplot(df, aes(x=log2FoldChange, y=log10_pvalue, color=Category)) +
    geom_point(alpha=0.6) +
    scale_color_manual(values = c(
      "Strong Up"      = "red3",
      "Moderate Up"    = "orange",
      "Moderate Down"  = "skyblue",
      "Strong Down"    = "blue3",
      "Nonsignificant" = "gray"
    )) +
    geom_vline(xintercept = c(-1.5,  1.5), linetype="dotted") +
    geom_vline(xintercept = c(-3.0,  3.0), linetype="dashed") +
    geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    labs(
      title    = contrast_name,
      subtitle = "|LFC| dotted=1.5, dashed=3",
      x        = "Log₂ Fold Change",
      y        = "-Log₁₀ Adj. P-value"
    ) +
    theme_minimal(base_size=10) +
    theme(
      legend.position = "none",
      plot.title      = element_text(size=10, hjust=0.5),
      plot.subtitle   = element_text(size=8, hjust=0.5)
    )
  
  return(p)
})

# give them nice names
names(volcano_plots) <- names(pairwise_sel)

# 3) Assemble into a 3×3 grid
collage <- wrap_plots(volcano_plots, ncol = 3)

# 4) Save the collage
ggsave(
  filename = "figures/volcanoes/volcano_collage.pdf",
  plot     = collage,
  width    = 10,   # adjust to taste
  height   = 10,
  dpi      = 300
)  
# 1. Define exactly the row‐by‐row order you want:
ordered_contrasts <- c(
  # FL row
  "FLVCvsFLNT", "FLVCvsFLVD", "FLVDvsFLNT",
  # NC row
  "NCVCvsNCNT", "NCVCvsNCVD", "NCVDvsNCNT",
  # NS row
  "NSVCvsNSNT", "NSVCvsNSVD", "NSVDvsNSNT"
)

# 2. Reorder your existing list of ggplots:
volcano_plots_ordered <- volcano_plots[ordered_contrasts]

# 3. Lay them out in a 3×3 grid:
collage <- wrap_plots(volcano_plots_ordered, ncol = 3)

# 4. Save the new collage
ggsave(
  filename = "volcanoes/volcano_collage_by_genotype.pdf",
  plot     = collage,
  width    = 10,    # tweak as needed
  height   = 10,
  dpi      = 300
)

#Upset Plot
expr_mat <- assay(rld)
# for each genotype (or exposure) get the genes that are >0 in ANY replicate
genes_by_genotype <- metadata %>%
  rownames_to_column("sample") %>%
  pull(sample, name = Condition) %>%
  split(names(.)) %>%               # a list: names = each Condition
  map(~ rownames(expr_mat)[rowMeans(expr_mat[, .x] > 0) > 0]) 

if (!requireNamespace("UpSetR", quietly=TRUE)) {
  install.packages("UpSetR")
}
library(UpSetR)
if (!requireNamespace("ComplexHeatmap", quietly=TRUE)) {
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)
library(circlize)

# 1) Build a binary membership matrix (genes × sets)
membership <- fromList(genes_by_genotype)

# 2) Make the combination matrix
cm <- make_comb_mat(membership)

# 3) Draw a basic UpSet plot
ComplexHeatmap::UpSet(cm)
  
# only keep the contrasts you care about
sel <- pairwise_results[desired]

library(purrr)
ordered_sets <- c(
  # FL up
  "FLVCvsFLNT_up", "FLVCvsFLVD_up", "FLVDvsFLNT_up",
  # NC up
  "NCVCvsNCNT_up", "NCVCvsNCVD_up", "NCVDvsNCNT_up",
  # NS up
  "NSVCvsNSNT_up", "NSVCvsNSVD_up", "NSVDvsNSNT_up",
  # FL down
  "FLVCvsFLNT_down", "FLVCvsFLVD_down", "FLVDvsFLNT_down",
  # NC down
  "NCVCvsNCNT_down", "NCVCvsNCVD_down", "NCVDvsNCNT_down",
  # NS down
  "NSVCvsNSNT_down", "NSVCvsNSVD_down", "NSVDvsNSNT_down"
)

sets_list <- imap(sel, function(res, nm) {
  up   <- rownames(res)[ res$padj < 0.05 & res$log2FoldChange >=  1.5 ]
  down <- rownames(res)[ res$padj < 0.05 & res$log2FoldChange <= -1.5 ]
  # return a small list with two named elements
  setNames(list(up, down),
           c(paste0(nm, "_up"), paste0(nm, "_down")))
}) %>% 
  flatten()

membership <- fromList(sets_list)
cm         <- make_comb_mat(membership)

# 3) Your custom set‐order
ordered_sets <- c(
  # FL up…
  "FLVCvsFLNT_up", "FLVCvsFLVD_up", "FLVDvsFLNT_up",
  # NC up…
  "NCVCvsNCNT_up", "NCVCvsNCVD_up", "NCVDvsNCNT_up",
  # NS up…
  "NSVCvsNSNT_up", "NSVCvsNSVD_up", "NSVDvsNSNT_up",
  # FL down…
  "FLVCvsFLNT_down","FLVCvsFLVD_down","FLVDvsFLNT_down",
  # NC down…
  "NCVCvsNCNT_down","NCVCvsNCVD_down","NCVDvsNCNT_down",
  # NS down…
  "NSVCvsNSNT_down","NSVCvsNSVD_down","NSVDvsNSNT_down"
)

# 4) Build a RowAnnotation for the right margin that:
#    • Shows set‐size bars
#    • Then writes the set names again
right_anno <- rowAnnotation(
  `Set size` = anno_barplot(
    set_size(cm),
    border = FALSE,
    gp     = gpar(fill = "darkgray")
  ),
  ` ` = anno_text(
    rownames(membership),
    location = 0.5,    # center the text in each row
    just     = "left",
    gp       = gpar(fontsize = 10)
  ),
  annotation_name_side = "left",
  width = unit.c(unit(4, "cm"), unit(4, "cm"))
)

# 5) Draw the UpSet with:
#    – intersection sizes on top
#    – set sizes + labels on the right
#    – row names on the left as usual
up <- UpSet(
  cm,
  set_order       = ordered_sets,
  comb_order      = order(comb_size(cm), decreasing = TRUE),
  row_names_side  = "left",
  top_annotation  = HeatmapAnnotation(
    `Intersection\nsize` = anno_barplot(
      comb_size(cm),
      border = FALSE,
      gp     = gpar(fill = "steelblue")
    ),
    annotation_name_side = "left"
  ),
  right_annotation = right_anno,
  pt_size         = unit(3, "mm")
)
#============================================================================================

#=======================WGCNA analysis=======================================================

#Clean and tidy data
col_sel <- names(merged_counts[,-1]) #Get all but first column names

mdata <- merged_counts %>%
  tidyr::pivot_longer(
    .,
    col = all_of(col_sel)
  ) %>%  
  mutate(
    group = sub("\\d+$", "", name)
  )

#Plot groups to identify outliers
(
  p <- mdata %>%
    ggplot(., aes(x = name, y = value)) +
    geom_violin() +
    geom_point(alpha = 0.2) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)          # Rotate treatment text
    ) +
    labs(x = "Treatment Groups", y = "RNA Seq Counts") +
    facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")      # Facet by hour
)

#Normalize counts with DEseq

library(DESeq2)

# 1) Pull out the matrix of raw counts (genes × samples)
count_mat <- as.matrix( merged_counts[ , -1 ] )
rownames(count_mat) <- merged_counts$Transcript_ID

# 2) Build a colData data.frame with rownames = sample names
sample_names <- colnames(count_mat)
cond <- gsub("\\d+$", "", sample_names)   # strip off trailing digits
colData <- data.frame(
  Condition = factor(cond),
  row.names  = sample_names
)

# 3) Check that everything matches
stopifnot( all(rownames(colData) == colnames(count_mat)) )

# 4) Finally call DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),  # integer matrix
  colData   = colData,           # data.frame with rownames = sample names
  design    = ~ Condition
)

dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds)
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

q75_wpn <- quantile( rowVars(wpn_vsd), .75)
q95_wpn <- quantile( rowVars(wpn_vsd), .95)
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )

input_mat <- t(expr_normalized)

library(WGCNA)
allowWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 11, to = 20, by = 1))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 10
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")




# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

#Examine expression analysis

# pick out a few modules of interest here
modules_of_interest = c("yellow", "red", "pink")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)
row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized[1:5,1:10]

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")

#Generate and export networks

genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )
# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")
#Make Module membership csv

modules <- unique(module_df$colors)
transcript_ids <- module_df$gene_id
transposed_MEs <- t(MEs0)
transposed_MEs <- transposed_MEs[-nrow(transposed_MEs), ]
mm_matrix <- matrix(nrow = nrow(expr_normalized), ncol = nrow(transposed_MEs))

for (gene in 1:nrow(expr_normalized)) {
  for (module in 1:nrow(transposed_MEs)) {
    # Calculate correlation between the gene's expression profile and the module eigengene
    mm_matrix[gene, module] <- cor(expr_normalized[gene, ], transposed_MEs[module, ])
  }
}
mm_df <- as.data.frame(mm_matrix)
rownames(mm_df) <- rownames(expr_normalized)
colnames(mm_df) <- rownames(transposed_MEs)
write.csv(mm_df, file = "modulemembership.csv")
