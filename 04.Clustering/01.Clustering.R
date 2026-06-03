BASE_DIR <- "/path/to/project"   # SET THIS

quant_dir <- file.path(BASE_DIR, "02.Quantification")
tpm_dir   <- file.path(BASE_DIR, "03.Normalization")
out_dir   <- file.path(BASE_DIR, "04.Clustering")
dir.create(out_dir, showWarnings = FALSE)

tpm     <- as.matrix(read.table(file.path(tpm_dir, "merged_matrix_TPM.txt"),
                                  sep = "\t", header = TRUE, row.names = 1,
                                  check.names = FALSE))
coldata <- read.table(file.path(quant_dir, "coldata.txt"),
                      sep = "\t", header = TRUE, row.names = 1,
                      stringsAsFactors = FALSE)
common  <- intersect(rownames(coldata), colnames(tpm))
tpm     <- tpm[, common]; coldata <- coldata[common, ]

# log2(TPM+1), top 5,000 variable genes
log_tpm <- log2(tpm + 1)
log_top <- log_tpm[order(apply(log_tpm, 1, var), decreasing = TRUE)[1:5000], ]

# PCA
pca    <- prcomp(t(log_top), center = TRUE, scale. = FALSE)
pct    <- round(summary(pca)$importance["Proportion of Variance", ] * 100, 1)
pca_df <- cbind(as.data.frame(pca$x[, 1:10]),
                coldata[rownames(pca$x), c("Sample_ID","Group","Day")])
write.table(pca_df, file.path(out_dir, "PCA_coords.txt"),
            sep = "\t", quote = FALSE)
cat(sprintf("PCA: PC1=%.1f%% PC2=%.1f%%\n", pct[1], pct[2]))

# Hierarchical clustering (ward.D2, k=2)
norm_mat <- t(scale(t(log_top)))
norm_mat[!is.finite(norm_mat)] <- 0
hc       <- hclust(dist(t(norm_mat)), method = "ward.D2")
clust_df <- data.frame(SRR     = hc$labels,
                       cluster = cutree(hc, k = 2),
                       stringsAsFactors = FALSE)[hc$order, ]
clust_df <- cbind(clust_df, coldata[clust_df$SRR, c("Sample_ID","Group","Day")])
write.table(clust_df, file.path(out_dir, "Cluster_assignments.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Clustering done\n")
