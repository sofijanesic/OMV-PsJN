##############################################
# Final Adjusted DESeq2 Workflow (SE50 OMV vs CELLS, Q30 filtered)
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
counts_SE50 <- read.delim("counts_SE50_unique30.txt", comment.char="#")
rownames(counts_SE50) <- counts_SE50$Geneid
countData_SE50 <- counts_SE50[ , -(1:6) ]

# -------------------------------
# 2. Define sample metadata
# -------------------------------
colData_SE50 <- data.frame(
  row.names = colnames(countData_SE50),
  Group = c("CELLS","CELLS","OMV","OMV")
)
colData_SE50$Group <- factor(colData_SE50$Group, levels = c("CELLS","OMV"))

# -------------------------------
# 3. Run DESeq2
# -------------------------------
dds_SE50 <- DESeqDataSetFromMatrix(countData = countData_SE50,
                                   colData = colData_SE50,
                                   design = ~ Group)
dds_SE50 <- DESeq(dds_SE50)

res_raw <- results(dds_SE50, contrast = c("Group","OMV","CELLS"))
res_SE50_shrink <- lfcShrink(dds_SE50,
                             coef = "Group_OMV_vs_CELLS",
                             type = "apeglm")

res_df_SE50 <- as.data.frame(res_raw)
res_df_SE50$log2FoldChange <- res_SE50_shrink$log2FoldChange
res_df_SE50$Locus_ID <- rownames(res_df_SE50)
res_df_SE50$minus_log10_pvalue <- -log10(res_df_SE50$pvalue)
res_df_SE50$comparison <- "OMV_vs_CELLS"

# -------------------------------
# 4. MA plot
# -------------------------------
plotMA(res_SE50_shrink, main="SE50")

# -------------------------------
# 5. PCA plot
# -------------------------------
vsd_SE50 <- vst(dds_SE50, blind=FALSE)
pcaData_SE50 <- plotPCA(vsd_SE50, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData_SE50, "percentVar"))

# Clean sample names
pcaData_SE50$name <- gsub("_SE50_bwa_dedup_unique30.bam", "", pcaData_SE50$name)
# Remove both the suffix and the Q30 prefix
pcaData_SE50$name <- gsub("Q30.|_SE50_bwa_dedup_unique30.bam$", "", pcaData_SE50$name)
ggplot(pcaData_SE50, aes(PC1, PC2, color=Group, label=name)) +
  geom_point(size=4) +
  geom_text(nudge_x = 10, nudge_y = 2) +
  labs(
    title = "PCA Plot (SE50)",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_minimal()

# -------------------------------
# 6. Import annotation and biotypes
# -------------------------------
gtf <- import("genomic.gtf")
gtf_genes <- gtf[gtf$type == "gene"]

anno_SE50 <- data.frame(
  Locus_ID = gtf_genes$locus_tag,
  biotype  = gtf_genes$gene_biotype
)

res_anno_SE50 <- merge(res_df_SE50, anno_SE50, by="Locus_ID", all.x=TRUE)

# -------------------------------
# 7. Summarize by biotype
# -------------------------------
summary_by_biotype_SE50 <- res_anno_SE50 %>%
  group_by(biotype) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange, na.rm=TRUE),
    sig_genes   = sum(padj < 0.05, na.rm=TRUE),
    total_genes = n()
  )
print(summary_by_biotype_SE50)

# -------------------------------
# 8. Export results per biotype
# -------------------------------
dir.create("biotype_results_SE50_Q30", showWarnings = FALSE)
biotypes_SE50 <- unique(na.omit(res_anno_SE50$biotype))

for (bt in biotypes_SE50) {
  safe_bt <- gsub("[^A-Za-z0-9_]", "_", bt)
  subset_df <- res_anno_SE50[res_anno_SE50$biotype == bt,
                             c("Locus_ID","log2FoldChange","pvalue","padj","biotype")]
  write.csv(subset_df,
            file = paste0("biotype_results_SE50_Q30/", safe_bt, "_results_SE50_Q30.csv"),
            row.names = FALSE)
}

# -------------------------------
# 9. VolcaNoseR-ready export
# -------------------------------
res_df_SE50$contrast   <- "Group: OMV vs CELLS"
res_df_SE50$reference  <- "CELLS"
res_df_SE50$condition  <- "OMV"

write.csv(res_df_SE50[order(res_df_SE50$padj), ],
          file = "DESeq2_results_OMV_vs_CELLS_SE50_unique30.csv",
          row.names = TRUE)

res_out <- res_df_SE50[, c("Locus_ID","baseMean","log2FoldChange","lfcSE","stat",
                           "pvalue","padj","minus_log10_pvalue","comparison")]
Sys.setlocale("LC_NUMERIC","C")
write.csv(res_out, file="Volcano_input_SE50_Q30.csv",
          row.names=FALSE, quote=FALSE)

res_minimal <- res_out[, c("Locus_ID","log2FoldChange","pvalue","padj","minus_log10_pvalue")]
write.csv(res_minimal, "Volcano_input_minimal_SE50_Q30.csv", row.names=FALSE)

# -------------------------------
# 10. Volcano plots per biotype (Q30)
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
                  title = 'OMV vs CELLS (SE50)',
                  subtitle = subtitle_text,
                  caption = 'DESeq2 results',
                  drawConnectors = TRUE,
                  arrowheads = FALSE,
                  lengthConnectors = 0.5,
                  max.overlaps = 20,
                  #ylim = c(0, 6), #samo za tRNA
                  #ylim = c(0, 15), #samo za pseudo
                  )
}

# Load per-biotype CSVs (Q30)
res_protein_SE50 <- read.csv("biotype_results_SE50_Q30/protein_coding_results_SE50_Q30.csv", row.names = 1)
res_rrna_SE50    <- read.csv("biotype_results_SE50_Q30/rRNA_results_SE50_Q30.csv", row.names = 1)
res_trna_SE50    <- read.csv("biotype_results_SE50_Q30/tRNA_results_SE50_Q30.csv", row.names = 1)
res_pseudo_SE50  <- read.csv("biotype_results_SE50_Q30/pseudogene_results_SE50_Q30.csv", row.names = 1)
res_miscRNA_SE50 <- read.csv("biotype_results_SE50_Q30/misc_RNA_results_SE50_Q30.csv", row.names = 1)

# Plot per biotype
plot_top10_volcano(res_protein_SE50, "Protein coding derived sRNAs")
plot_top10_volcano(res_rrna_SE50,    "rRNA derived sRNAs")
plot_top10_volcano(res_trna_SE50,    "tRNA derived sRNAs")
plot_top10_volcano(res_pseudo_SE50,  "Pseudogene derived sRNAs")
plot_top10_volcano(res_miscRNA_SE50, "Misc RNA derived sRNAs")

# -------------------------------
# 11. General Volcano Plot (all genes, Q30)
# -------------------------------
keyvals <- ifelse(is.na(res_df_SE50$padj) | is.na(res_df_SE50$log2FoldChange), '#949494',
                  ifelse(res_df_SE50$padj < 0.05 & res_df_SE50$log2FoldChange > 0.6, '#FB9528',
                         ifelse(res_df_SE50$padj < 0.05 & res_df_SE50$log2FoldChange < -0.6, '#38C0FF',
                                '#949494')))
names(keyvals) <- rownames(res_df_SE50)
names(keyvals)[keyvals == '#FB9528'] <- 'Upregulated'
names(keyvals)[keyvals == '#38C0FF'] <- 'Downregulated'
names(keyvals)[keyvals == '#949494'] <- 'Not significant'

res_clean_SE50 <- res_df_SE50[!is.na(res_df_SE50$padj) & !is.na(res_df_SE50$log2FoldChange), ]
res_ordered_SE50 <- res_clean_SE50[order(res_clean_SE50$padj), ]

up_candidates_SE50 <- rownames(res_ordered_SE50[res_ordered_SE50$padj < 0.05 & res_ordered_SE50$log2FoldChange > 0.6, ])
down_candidates_SE50 <- rownames(res_ordered_SE50[res_ordered_SE50$padj < 0.05 & res_ordered_SE50$log2FoldChange < -0.6, ])

top_up_SE50 <- head(up_candidates_SE50, min(10, length(up_candidates_SE50)))
top_down_SE50 <- head(down_candidates_SE50, min(10, length(down_candidates_SE50)))

select_genes_SE50 <- c(top_up_SE50, top_down_SE50)

EnhancedVolcano(res_df_SE50,
                lab = rownames(res_df_SE50),
                selectLab = select_genes_SE50,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.6,
                colCustom = keyvals,
                title = 'OMV vs CELLS (SE50)',
                subtitle = "DESeq2 results",
                caption = 'DESeq2 results',
                drawConnectors = TRUE,
                arrowheads = FALSE,
                lengthConnectors = 10,
                max.overlaps = Inf)

# -------------------------------
# 12. Heatmap of top variable genes (Q30)
# -------------------------------
vsd_SE50 <- vst(dds_SE50, blind=FALSE)
mat_SE50 <- assay(vsd_SE50)

# Select top 20 most variable genes
topVarGenes_SE50 <- head(order(rowVars(mat_SE50), decreasing=TRUE), 20)

# Scale rows for visualization
mat_scaled_SE50 <- mat_SE50[topVarGenes_SE50, ]
mat_scaled_SE50 <- t(scale(t(mat_scaled_SE50)))

# Rename columns to simplified sample names
colnames(mat_scaled_SE50) <- c("CELLS1","CELLS3","OMV1","OMV3")

# Fix annotation rownames to match simplified sample names
rownames(colData_SE50) <- c("CELLS1","CELLS3","OMV1","OMV3")

# Now plot heatmap
pheatmap(mat_scaled_SE50,
         annotation_col = colData_SE50,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         main = "Top 20 variable genes (SE50 - DESeq2)")

# 13. Heatmap of top variable genes +biotzpe (Q30)
# -------------------------------

vsd_SE50 <- vst(dds_SE50, blind=FALSE)
mat_SE50 <- assay(vsd_SE50)

# Select top 20 most variable genes
topVarGenes_SE50 <- head(order(rowVars(mat_SE50), decreasing=TRUE), 20)

# Scale rows for visualization
mat_scaled_SE50 <- mat_SE50[topVarGenes_SE50, ]
mat_scaled_SE50 <- t(scale(t(mat_scaled_SE50)))

# Rename columns to simplified sample names
colnames(mat_scaled_SE50) <- c("CELLS1","CELLS3","OMV1","OMV3")

# Fix annotation rownames to match simplified sample names
rownames(colData_SE50) <- c("CELLS1","CELLS3","OMV1","OMV3")

# -------------------------------
# Add biotype annotation for rows
# -------------------------------
# Match biotype info to selected genes
row_anno <- data.frame(Biotype = res_anno_SE50$biotype[match(rownames(mat_scaled_SE50),
                                                             res_anno_SE50$Locus_ID)])
rownames(row_anno) <- rownames(mat_scaled_SE50)

# Define colors for biotypes (optional)
biotype_colors <- c("protein_coding"="blue",
                    "lncRNA"="green",
                    "pseudogene"="red")

# Now plot heatmap with both column and row annotations
pheatmap(mat_scaled_SE50,
         annotation_col = colData_SE50,
         annotation_row = row_anno,
         annotation_colors = list(Biotype = biotype_colors),
         show_rownames = TRUE,
         cluster_cols = TRUE,
         main = "Top 20 variable genes (SE50 - DESeq2)")
