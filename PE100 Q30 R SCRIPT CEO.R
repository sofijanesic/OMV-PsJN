##############################################
# Final Adjusted DESeq2 Workflow (PE100 OMV vs CELLS, Q30 filtered)
##############################################

setwd("C:/Users/sofija.nesic/Documents/Sof R 30")

# Load required libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)
library(matrixStats)

# -------------------------------
# 1. Read featureCounts output (Q30)
# -------------------------------
counts_PE100 <- read.delim("counts_PE100_unique30.txt", comment.char="#")
rownames(counts_PE100) <- counts_PE100$Geneid
countData_PE100 <- counts_PE100[ , -(1:6) ]

# -------------------------------
# 2. Define sample metadata
# -------------------------------
colData_PE100 <- data.frame(
  row.names = colnames(countData_PE100),
  Group = c("CELLS","CELLS","OMV","OMV")
)
colData_PE100$Group <- factor(colData_PE100$Group, levels = c("CELLS","OMV"))

# -------------------------------
# 3. Run DESeq2
# -------------------------------
dds_PE100 <- DESeqDataSetFromMatrix(countData = countData_PE100,
                                    colData = colData_PE100,
                                    design = ~ Group)
dds_PE100 <- DESeq(dds_PE100)

res_raw_PE100 <- results(dds_PE100, contrast = c("Group","OMV","CELLS"))
res_PE100_shrink <- lfcShrink(dds_PE100,
                              coef = "Group_OMV_vs_CELLS",
                              type = "apeglm")

res_df_PE100 <- as.data.frame(res_raw_PE100)
res_df_PE100$log2FoldChange <- res_PE100_shrink$log2FoldChange
res_df_PE100$Locus_ID <- rownames(res_df_PE100)
res_df_PE100$minus_log10_pvalue <- -log10(res_df_PE100$pvalue)
res_df_PE100$comparison <- "OMV_vs_CELLS"

# -------------------------------
# 4. MA plot
# -------------------------------
plotMA(res_PE100_shrink, main="PE100")

# -------------------------------
# 5. PCA plot
# -------------------------------
vsd_PE100 <- vst(dds_PE100, blind=FALSE)
pcaData_PE100 <- plotPCA(vsd_PE100, intgroup="Group", returnData=TRUE)
percentVar_PE100 <- round(100 * attr(pcaData_PE100, "percentVar"))

# Clean sample names
pcaData_PE100$name <- gsub("Q30.|_PE100_bwa_dedup_unique30.bam", "", pcaData_PE100$name)

ggplot(pcaData_PE100, aes(PC1, PC2, color=Group, label=name)) +
  geom_point(size=4) +
  geom_text(nudge_x = 15, nudge_y = 2) +
  labs(
    title = "PCA Plot (PE100)",
    x = paste0("PC1: ", percentVar_PE100[1], "% variance"),
    y = paste0("PC2: ", percentVar_PE100[2], "% variance")
  ) +
  theme_minimal()

library(ggrepel)

ggplot(pcaData_PE100, aes(PC1, PC2, color=Group, label=name)) +
  geom_point(size=4) +
  geom_text_repel(max.overlaps = 20, nudge_x = 10) +
  labs(
    title = "PCA Plot (PE100)",
    x = paste0("PC1: ", percentVar_PE100[1], "% variance"),
    y = paste0("PC2: ", percentVar_PE100[2], "% variance")
  ) +
  theme_minimal()
# -------------------------------
# 6. Import annotation and biotypes
# -------------------------------
gtf <- import("genomic.gtf")
gtf_genes <- gtf[gtf$type == "gene"]

anno_PE100 <- data.frame(
  Locus_ID = gtf_genes$locus_tag,
  biotype  = gtf_genes$gene_biotype
)

res_anno_PE100 <- merge(res_df_PE100, anno_PE100, by="Locus_ID", all.x=TRUE)

# -------------------------------
# 7. Summarize by biotype
# -------------------------------
summary_by_biotype_PE100 <- res_anno_PE100 %>%
  group_by(biotype) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange, na.rm=TRUE),
    sig_genes   = sum(padj < 0.05, na.rm=TRUE),
    total_genes = n()
  )
print(summary_by_biotype_PE100)

# -------------------------------
# 8. Export results per biotype
# -------------------------------
dir.create("biotype_results_PE100_Q30", showWarnings = FALSE)
biotypes_PE100 <- unique(na.omit(res_anno_PE100$biotype))

for (bt in biotypes_PE100) {
  safe_bt <- gsub("[^A-Za-z0-9_]", "_", bt)
  subset_df <- res_anno_PE100[res_anno_PE100$biotype == bt,
                              c("Locus_ID","log2FoldChange","pvalue","padj","biotype")]
  write.csv(subset_df,
            file = paste0("biotype_results_PE100_Q30/", safe_bt, "_results_PE100_Q30.csv"),
            row.names = FALSE)
}

# -------------------------------
# 9. VolcaNoseR-ready export (PE100 Q30)
# -------------------------------
res_df_PE100$contrast   <- "Group: OMV vs CELLS"
res_df_PE100$reference  <- "CELLS"
res_df_PE100$condition  <- "OMV"

write.csv(res_df_PE100[order(res_df_PE100$padj), ],
          file = "DESeq2_results_OMV_vs_CELLS_PE100_unique30.csv",
          row.names = TRUE)

res_out_PE100 <- res_df_PE100[, c("Locus_ID","baseMean","log2FoldChange","lfcSE","stat",
                                  "pvalue","padj","minus_log10_pvalue","comparison")]
Sys.setlocale("LC_NUMERIC","C")
write.csv(res_out_PE100, file="Volcano_input_PE100_Q30.csv",
          row.names=FALSE, quote=FALSE)

res_minimal_PE100 <- res_out_PE100[, c("Locus_ID","log2FoldChange","pvalue","padj","minus_log10_pvalue")]
write.csv(res_minimal_PE100, "Volcano_input_minimal_PE100_Q30.csv", row.names=FALSE)

# -------------------------------
# 10. Volcano plots per biotype (PE100 Q30)
# -------------------------------
plot_top10_volcano <- function(res_df, subtitle_text) {
  keyvals <- ifelse(is.na(res_df$padj) | is.na(res_df$log2FoldChange), '#949494',
                    ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 0.6, '#FB9528',
                           ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -0.6, '#38C0FF',
                                  '#949494')))
  names(keyvals) <- rownames(res_df)
  names(keyvals)[keyvals == '#FB9528'] <- 'Upregulated'
  names(keyvals)[keyvals == '#38C0FF'] <- 'Downregulated'
  names(keyvals)[keyvals == '#949494'] <- 'Not significant'
  
  res_clean <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
  res_ordered <- res_clean[order(res_clean$padj), ]
  
  up_candidates <- rownames(res_ordered[res_ordered$padj < 0.05 & res_ordered$log2FoldChange > 0.6, ])
  down_candidates <- rownames(res_ordered[res_ordered$padj < 0.05 & res_ordered$log2FoldChange < -0.6, ])
  
  top_up <- head(up_candidates, min(5, length(up_candidates)))
  top_down <- head(down_candidates, min(5, length(down_candidates)))
  
  select_genes <- c(top_up, top_down)
  
  EnhancedVolcano(res_df,
                  lab = rownames(res_df),
                  selectLab = select_genes,
                  x = 'log2FoldChange',
                  y = 'padj',
                  pCutoff = 0.05,
                  FCcutoff = 0.6,
                  colCustom = keyvals,
                  title = 'OMV vs CELLS (PE100)',
                  subtitle = subtitle_text,
                  caption = 'DESeq2 results',
                  drawConnectors = TRUE,
                  arrowheads = FALSE,
                  lengthConnectors = 0.05,
                  max.overlaps = Inf)
}

# Load per-biotype CSVs (PE100 Q30)
res_protein_PE100 <- read.csv("biotype_results_PE100_Q30/protein_coding_results_PE100_Q30.csv", row.names = 1)
res_rrna_PE100    <- read.csv("biotype_results_PE100_Q30/rRNA_results_PE100_Q30.csv", row.names = 1)
res_trna_PE100    <- read.csv("biotype_results_PE100_Q30/tRNA_results_PE100_Q30.csv", row.names = 1)
res_pseudo_PE100  <- read.csv("biotype_results_PE100_Q30/pseudogene_results_PE100_Q30.csv", row.names = 1)
res_miscRNA_PE100 <- read.csv("biotype_results_PE100_Q30/misc_RNA_results_PE100_Q30.csv", row.names = 1)

# Plot per biotype volcano plots
plot_top10_volcano(res_protein_PE100, "Protein coding derived sRNAs (PE100)")
plot_top10_volcano(res_rrna_PE100,    "rRNA derived sRNAs (PE100)")
plot_top10_volcano(res_trna_PE100,    "tRNA derived sRNAs (PE100)")
plot_top10_volcano(res_pseudo_PE100,  "Pseudogene derived sRNAs (PE100)")
plot_top10_volcano(res_miscRNA_PE100, "Misc RNA derived sRNAs (PE100)")

# -------------------------------
# 11. General Volcano Plot (all genes, PE100 Q30)
# -------------------------------
keyvals <- ifelse(is.na(res_df_PE100$padj) | is.na(res_df_PE100$log2FoldChange), '#949494',
                  ifelse(res_df_PE100$padj < 0.05 & res_df_PE100$log2FoldChange > 0.6, '#FB9528',
                         ifelse(res_df_PE100$padj < 0.05 & res_df_PE100$log2FoldChange < -0.6, '#38C0FF',
                                '#949494')))
names(keyvals) <- rownames(res_df_PE100)
names(keyvals)[keyvals == '#FB9528'] <- 'Upregulated'
names(keyvals)[keyvals == '#38C0FF'] <- 'Downregulated'
names(keyvals)[keyvals == '#949494'] <- 'Not significant'

res_clean_PE100 <- res_df_PE100[!is.na(res_df_PE100$padj) & !is.na(res_df_PE100$log2FoldChange), ]
res_ordered_PE100 <- res_clean_PE100[order(res_clean_PE100$padj), ]

up_candidates_PE100 <- rownames(res_ordered_PE100[res_ordered_PE100$padj < 0.05 & res_ordered_PE100$log2FoldChange > 0.6, ])
down_candidates_PE100 <- rownames(res_ordered_PE100[res_ordered_PE100$padj < 0.05 & res_ordered_PE100$log2FoldChange < -0.6, ])

top_up_PE100 <- head(up_candidates_PE100, min(10, length(up_candidates_PE100)))
top_down_PE100 <- head(down_candidates_PE100, min(10, length(down_candidates_PE100)))

select_genes_PE100 <- c(top_up_PE100, top_down_PE100)

EnhancedVolcano(res_df_PE100,
                lab = rownames(res_df_PE100),
                selectLab = select_genes_PE100,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.6,
                colCustom = keyvals,
                title = 'OMV vs CELLS (PE100)',
                subtitle = "DESeq2 results",
                caption = 'DESeq2 results',
                drawConnectors = TRUE,
                arrowheads = FALSE,
                lengthConnectors = 10,
                max.overlaps = Inf)

# -------------------------------
# 12. Heatmap of top variable genes (PE100 Q30)
# -------------------------------
vsd_PE100 <- vst(dds_PE100, blind=FALSE)
mat_PE100 <- assay(vsd_PE100)

# Select top 20 most variable genes
topVarGenes_PE100 <- head(order(rowVars(mat_PE100), decreasing=TRUE), 20)

# Scale rows for visualization
mat_scaled_PE100 <- mat_PE100[topVarGenes_PE100, ]
mat_scaled_PE100 <- t(scale(t(mat_scaled_PE100)))

# Rename columns to simplified sample names
colnames(mat_scaled_PE100) <- c("CELLS1","CELLS3","OMV1","OMV3")
rownames(colData_PE100) <- c("CELLS1","CELLS3","OMV1","OMV3")  # ensure annotation matches

# Plot heatmap
pheatmap(mat_scaled_PE100,
         annotation_col = colData_PE100,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         main = "Top 20 variable genes (PE100 - DESeq2)")

# 13. Heatmap of top variable genes + biotzpe(Q30)
# -------------------------------

vsd_PE100 <- vst(dds_PE100, blind=FALSE)   # still using dds_SE50 object, but naming output PE100
mat_PE100 <- assay(vsd_PE100)

# Select top 20 most variable genes
topVarGenes_PE100 <- head(order(rowVars(mat_PE100), decreasing=TRUE), 20)

# Scale rows for visualization
mat_scaled_PE100 <- mat_PE100[topVarGenes_PE100, ]
mat_scaled_PE100 <- t(scale(t(mat_scaled_PE100)))

# Rename columns to simplified sample names
colnames(mat_scaled_PE100) <- c("CELLS1","CELLS3","OMV1","OMV3")

# Fix annotation rownames to match simplified sample names
rownames(colData_PE100) <- c("CELLS1","CELLS3","OMV1","OMV3")

# -------------------------------
# Add biotype annotation for rows
# -------------------------------
row_anno <- data.frame(Biotype = res_anno_PE100$biotype[match(rownames(mat_scaled_PE100),
                                                             res_anno_PE100$Locus_ID)])
rownames(row_anno) <- rownames(mat_scaled_PE100)

biotype_colors <- c("protein_coding"="blue",
                    "lncRNA"="green",
                    "pseudogene"="red")

# Now plot heatmap with both column and row annotations
pheatmap(mat_scaled_PE100,
         annotation_col = colData_PE100,
         annotation_row = row_anno,
         annotation_colors = list(Biotype = biotype_colors),
         show_rownames = TRUE,
         cluster_cols = TRUE,
         main = "Top 20 variable genes (PE100 - DESeq2)")