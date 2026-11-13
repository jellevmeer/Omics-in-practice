library(DESeq2)
library(data.table)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

## Data Transformation ----
# Importing counts matrix, metadata and conversion dataframes
counts_df <- fread("data/gene_counts.tsv", check.names = F, data.table = FALSE)
metadata <- fread("metadata/sample_metadata.tsv", data.table = F)

gene_conversion <- fread("../results/starsolo/Solo.out/Gene/raw/features.tsv", header = FALSE, data.table = FALSE)

barcode_conversion <- fread("metadata/barcode_samples.tsv",header = FALSE, data.table = FALSE)
barcode_conversion <- barcode_conversion[-1,]

# Converting to gene ids and sample numbers instead of barcodes 
gene_ids <- counts_df[[1]]
counts <- as.matrix(counts_df[, -1, drop = FALSE])

new_cols <- barcode_conversion[[2]][match(colnames(counts), barcode_conversion[[1]])]
colnames(counts) <- ifelse(is.na(new_cols), colnames(df), new_cols)


new_rows <- gene_conversion[[2]][match(gene_ids, gene_conversion[[1]])]
gene_ids <- ifelse(is.na(new_rows), gene_ids, new_rows)


rownames(counts) <- gene_ids
mode(counts) <- "numeric"

# Identifying and collapsing duplicate rownames
duplicaterows <- counts[duplicated(rownames(counts)) | duplicated(rownames(counts), fromLast = TRUE), ]
counts <- rowsum(counts, group = rownames(counts))

# Check if the transformations have gone correctly
stopifnot(!anyDuplicated(rownames(counts)))
stopifnot(!anyDuplicated(metadata$sample))
stopifnot(all(counts >= 0))
stopifnot(all(counts == round(counts)))
stopifnot(setequal(colnames(counts), metadata$sample))

metadata <- metadata[match(colnames(counts), metadata$sample), ]
mode(metadata$sample) <- "character"
stopifnot(identical(colnames(counts), metadata$sample))
rownames(metadata) <- metadata$sample

# Initial summary
dim(counts)
summary(colSums(counts))
table(metadata$condition)
head(metadata)

# Saving counts matrix containing sample names and gene ids
fwrite(counts,file = "data/gene_sample_counts_.tsv",
sep = "\t",
quote = FALSE,
row.names = TRUE,
col.names = TRUE
)

# Library size before normalization
lib_df <- data.frame(sample = colnames(counts), library_size = colSums(counts), condition = metadata$condition)

p_lib <- ggplot(lib_df, aes(sample, library_size, fill = condition)) +
geom_col() +
scale_y_continuous(labels = scales::comma) +
labs(x = NULL, y = "Total raw counts", title = "Library sizes") +
theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/library_sizes.pdf", p_lib, width = 8, height = 5)


## Creating DESeq2 dataset----
metadata$condition <- factor(metadata$condition)
metadata$condition <- relevel(metadata$condition, ref = "Control")

dds <- DESeqDataSetFromMatrix(countData = counts,colData = metadata,
design = ~ condition
)

# Filtering low expressed genes
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
cat("Genes retained:", nrow(dds), "\n")

saveRDS(dds, "results/dds_prefiltered.rds")


## PCA 
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

pca_data <- plotPCA(vsd, intgroup = c("condition", "treatment"), returnData = TRUE)

percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
geom_point(size = 3) +
geom_text_repel(show.legend = FALSE, max.overlaps = Inf) +
xlab(paste0("PC1: ", percent_var[1], "% variance")) +
ylab(paste0("PC2: ", percent_var[2], "% variance")) +
coord_equal() + theme_bw() +
labs(title = "PCA of variance-stabilized counts")

ggsave("figures/PCA_samples.pdf", p_pca, width = 7, height = 5.5)
fwrite(pca_data, "results/PCA_coordinates.tsv", sep = "\t")

## Sample-to-tsample distance heatmap
sample_dists <- dist(t(vsd_mat))
sample_dist_mat <- as.matrix(sample_dists)
rownames(sample_dist_mat) <- colnames(vsd_mat)
colnames(sample_dist_mat) <- colnames(vsd_mat)

ann_col <- metadata[, c("condition"), drop = FALSE]

pdf("figures/sample_distance_heatmap.pdf", width = 7, height = 6)
pheatmap(sample_dist_mat, annotation_col = ann_col, annotation_row = ann_col,main = "Sample-to-sample distances")
dev.off()


## Heatmap (z-score) of genes of interest
selected_genes <- scan("metadata/selected_genes.txt", what = character(), quiet = TRUE)

present <- intersect(selected_genes, rownames(vsd_mat))

missing <- setdiff(selected_genes, rownames(vsd_mat))

writeLines(missing, "results/selected_genes_not_found.txt")

if (length(present) < 2) stop("Fewer than two selected genes were found")

heat_mat <- vsd_mat[present, , drop = FALSE]
heat_z <- t(scale(t(heat_mat)))
heat_z <- heat_z[complete.cases(heat_z), , drop = FALSE]

ann_col <- metadata[, c("condition"), drop = FALSE]

pdf("figures/selected_gene_heatmap_zscore.pdf", width = 8, height = 7)
pheatmap(heat_z, annotation_col = ann_col,
cluster_rows = TRUE, cluster_cols = TRUE,
show_colnames = TRUE,
main = "Selected genes: row z-scores")

dev.off()

# Heatmap with ordered samples
sample_order <- metadata$sample[order(metadata$condition)]
heat_z_ordered <- heat_z[, sample_order, drop = FALSE]

pdf("figures/selected_gene_heatmap_ordered.pdf", width = 8, height = 7)
pheatmap(heat_z_ordered, annotation_col = ann_col[sample_order, , drop = FALSE],
cluster_rows = TRUE, cluster_cols = FALSE,
main = "Selected genes by experimental group")

dev.off()

# Heatmap without row scaling
pdf("figures/selected_gene_heatmap_VST.pdf", width = 8, height = 7)
pheatmap(heat_mat, annotation_col = ann_col, scale = "none",
main = "Selected genes: VST expression")

dev.off()



## Differential Expression Analysis
dds <- DESeq(dds)
resultsNames(dds)

treated_condition <- "Treatment_6h"

res <- results(dds, contrast = c("condition", treated_condition, "Control"), alpha = 0.05) # Contrast requires the specific input name of the sample!

summary(res)

# Shrinking Log2FC for ranking and visualization
res_shrunk <- lfcShrink(
dds, contrast = c("condition", treated_condition, "Control"),
res = res,
type = "apeglm"
)
# If above function does not run, coef = integer, where int is the position of the comparison in resultsNames(dds) + remove contrast = x statment.

res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
filename <- paste("results/DE_", treated_condition, "_vs_Control.tsv")
fwrite(res_df, filename, sep = "\t", na = "NA")


## DESeq2 results interpretation
padj_cutoff <- 0.05
lfc_cutoff <- 1

res_df$status <- "Not significant"
res_df$status[!is.na(res_df$padj) & res_df$padj < padj_cutoff &
res_df$log2FoldChange >= lfc_cutoff] <- "Higher in Treatment"
res_df$status[!is.na(res_df$padj) & res_df$padj < padj_cutoff &
res_df$log2FoldChange <= -lfc_cutoff] <- "Lower in Treatment"

table(res_df$status)

# MA plot (Log2FC vs avg normalized expression)
pdf("figures/MA_Treatment_vs_Control.pdf", width = 7, height = 5)
plotMA(res_shrunk, alpha = 0.05, ylim = c(-5, 5),
main = "Treatment versus Control")
dev.off()


## Volcano plot
plot_df <- res_df
plot_df$minus_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
plot_df$minus_log10_padj[is.na(plot_df$padj)] <- NA_real_


label_pool <- subset(plot_df, status != "Not significant" & !is.na(padj))
label_genes <- head(label_pool[order(label_pool$padj), ], 12)

p_volcano <- ggplot(plot_df, aes(log2FoldChange, minus_log10_padj, color = status)) +
geom_point(alpha = 0.65, size = 1.4, na.rm = TRUE) +
geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2) +
geom_hline(yintercept = -log10(padj_cutoff), linetype = 2) +
scale_color_manual(values = c("Lower in Treatment" = "#377EB8",
"Not significant" = "grey70", "Higher in Treatment" = "#E41A1C")) +
geom_text_repel(data = label_genes, aes(label = gene_id),
size = 3, max.overlaps = Inf, show.legend = FALSE) +
labs(title = "Treatment 6h versus Control",
x = "Shrunken log2 fold change",
y = expression(-log[10](adjusted~p)), color = NULL) +
theme_bw()

ggsave("figures/volcano_Treatment_6h_vs_Control.pdf", p_volcano,
width = 7.5, height = 6)

ggsave("figures/volcano_Treatment_6h_vs_Control.png", p_volcano,
width = 7.5, height = 6, dpi = 300)


## DE expression validation
gene_to_plot <- "ZBTB33"
count_df <- plotCounts(dds, gene = gene_to_plot,
intgroup = c("condition"), returnData = TRUE)
p_count <- ggplot(count_df, aes(condition, count, color =condition)) +
geom_point(position = position_jitter(width = 0.08), size = 3) +
scale_y_log10() +
theme_bw() +
labs(title = gene_to_plot, x = "Condition", y = "Normalized count (log scale)") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 


ggsave(paste0("figures/counts_", gene_to_plot, ".pdf"),
p_count, width = 5, height = 4.5)


writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
