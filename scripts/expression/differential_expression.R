#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop(paste(
        "Usage: differential_expression.R COUNTS SAMPLE_TABLE OUT_PREFIX",
        "DESIGN CONTRAST [SUBSET]",
        "\n  Example DESIGN='~ sex + genotype' CONTRAST='genotype,APOE4,APOE3'",
        "SUBSET='cell_type=astrocytes;age=40 weeks'"
    ))
}

counts_file <- args[[1]]
sample_file <- args[[2]]
out_prefix <- args[[3]]
design_text <- args[[4]]
contrast_parts <- strsplit(args[[5]], ",", fixed = TRUE)[[1]]
subset_text <- if (length(args) >= 6) args[[6]] else ""
if (length(contrast_parts) != 3) stop("CONTRAST must be factor,alternative,reference")

metadata <- read.delim(sample_file, check.names = FALSE, stringsAsFactors = FALSE)
if (nzchar(subset_text)) {
    clauses <- strsplit(subset_text, ";", fixed = TRUE)[[1]]
    for (clause in clauses) {
        pieces <- strsplit(clause, "=", fixed = TRUE)[[1]]
        if (length(pieces) != 2 || !pieces[[1]] %in% colnames(metadata)) {
            stop(paste("Invalid subset clause:", clause))
        }
        metadata <- metadata[metadata[[pieces[[1]]]] == pieces[[2]], , drop = FALSE]
    }
}

fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
count_matrix <- as.matrix(fc[, 7:ncol(fc), drop = FALSE])
rownames(count_matrix) <- fc$Geneid
sample_names <- sub("_Aligned\\.out\\.bam$", "", basename(colnames(count_matrix)))
colnames(count_matrix) <- sample_names
metadata <- metadata[match(sample_names, metadata$sample), , drop = FALSE]
if (anyNA(metadata$sample)) stop("BAM columns and sample table do not match")
rownames(metadata) <- metadata$sample

design_formula <- as.formula(design_text)
design_vars <- all.vars(design_formula)
missing_vars <- setdiff(design_vars, colnames(metadata))
if (length(missing_vars)) stop(paste("Missing design columns:", paste(missing_vars, collapse = ", ")))
for (variable in design_vars) {
    if (any(!nzchar(metadata[[variable]]))) stop(paste("Empty values in design column", variable))
    metadata[[variable]] <- factor(metadata[[variable]])
}

factor_name <- contrast_parts[[1]]
alternative <- contrast_parts[[2]]
reference <- contrast_parts[[3]]
if (!factor_name %in% design_vars) stop("Contrast factor must occur in DESIGN")
observed <- levels(metadata[[factor_name]])
if (!all(c(alternative, reference) %in% observed)) {
    stop(paste("Contrast levels not present. Observed:", paste(observed, collapse = ", ")))
}

dds <- DESeqDataSetFromMatrix(round(count_matrix), metadata, design_formula)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds)
result <- results(dds, contrast = c(factor_name, alternative, reference), alpha = 0.05)
result_df <- data.frame(gene_id = rownames(result), as.data.frame(result), check.names = FALSE)
result_df <- result_df[order(result_df$padj, na.last = TRUE), ]

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
write.table(result_df, paste0(out_prefix, ".all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
significant <- result_df[!is.na(result_df$padj) & result_df$padj < 0.05 & abs(result_df$log2FoldChange) >= 1, ]
write.table(significant, paste0(out_prefix, ".significant.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene_id = rownames(dds), counts(dds, normalized = TRUE), check.names = FALSE),
            paste0(out_prefix, ".normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(paste0(out_prefix, ".MA.pdf"), width = 6, height = 6)
plotMA(result, alpha = 0.05, ylim = c(-5, 5))
dev.off()

vsd <- vst(dds, blind = FALSE)
pca <- prcomp(t(assay(vsd)))
percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
colors <- as.integer(metadata[[factor_name]])
pdf(paste0(out_prefix, ".PCA.pdf"), width = 7, height = 6)
plot(pca$x[, 1], pca$x[, 2], col = colors, pch = 19,
     xlab = paste0("PC1: ", percent[[1]], "%"), ylab = paste0("PC2: ", percent[[2]], "%"))
text(pca$x[, 1], pca$x[, 2], labels = rownames(metadata), pos = 3, cex = 0.6)
legend("topright", legend = levels(metadata[[factor_name]]), col = seq_along(levels(metadata[[factor_name]])), pch = 19)
dev.off()

writeLines(capture.output(sessionInfo()), paste0(out_prefix, ".sessionInfo.txt"))
