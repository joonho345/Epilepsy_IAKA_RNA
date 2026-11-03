# Revision 1 (Reviewer 1, Major comment 7): IAKA vs IHKA correlations computed separately
# for each IHKA dataset, because hippocampal subregion is confounded with dataset identity.
#   CA1 - GSE99577  (raw counts, single dataset), ~ Treatment, phases HA/AC/IM/CR
#   DG  - GSE213393 (raw counts, single dataset), ~ Treatment, phases AC/IM/CR
# IAKA log2FC vectors are read from DEG_Correlation.R output (DEG_IAKA_*.csv).
# Output: Pearson_correlation_by_dataset.csv (Fig. 4c, d)
library(DESeq2); library(dplyr)
set.seed(42)

BASE_DIR <- "/path/to/project"    # SET THIS
GO_PATH  <- "/path/to/GO.txt"     # SET THIS
N_BOOT <- 1000; N_PERM <- 1000

quant_dir <- file.path(BASE_DIR, "07.CrossDataset", "01.Quantification")
deg_dir   <- file.path(BASE_DIR, "07.CrossDataset", "03.DEGCorrelation")
out_dir   <- file.path(BASE_DIR, "07.CrossDataset", "04.IHKAByDataset")
dir.create(out_dir, showWarnings = FALSE)

lines     <- readLines(GO_PATH)
epi_genes <- unique(trimws(unlist(lapply(
  lines[grepl("^[A-Z][a-z].*\tn=", lines)],
  function(l) strsplit(sub(".*\tn=[0-9]+\t", "", l), ", ")[[1]]))))

col <- read.table(file.path(quant_dir, "coldata_ihka.txt"), sep = "\t", header = TRUE,
                  row.names = 1, stringsAsFactors = FALSE)
cnt <- round(data.matrix(read.table(file.path(quant_dir, "merged_matrix_ihka.txt"), sep = "\t",
                                      header = TRUE, row.names = 1, check.names = FALSE)))
k <- intersect(rownames(col), colnames(cnt)); col <- col[k, ]; cnt <- cnt[, k]

run_deseq2 <- function(batch, phase, label) {
  sel   <- rownames(col)[col$Batch == batch & col$Phase == phase]
  col_s <- col[sel, ]; col_s$Treatment <- factor(col_s$Treatment, c("Saline", "KA"))
  dds   <- DESeqDataSetFromMatrix(cnt[, sel], col_s, ~ Treatment)
  n_min <- min(table(col_s$Treatment))
  dds   <- dds[rowSums(counts(dds) >= 10) >= max(2L, n_min %/% 2), ]
  dds   <- DESeq(dds, quiet = TRUE)
  res   <- as.data.frame(results(dds, contrast = c("Treatment", "KA", "Saline")))
  res$gene <- rownames(res)
  write.csv(res[order(res$padj, na.last = TRUE), ],
            file.path(out_dir, paste0("DEG_", label, ".csv")), row.names = FALSE, quote = FALSE)
  setNames(res$log2FoldChange, rownames(res))
}

by_ds <- list(
  "CA1 (GSE99577)" = list(HA = run_deseq2("GSE99577", "HA", "IHKA_CA1_HA"),
                          AC = run_deseq2("GSE99577", "AC", "IHKA_CA1_AC"),
                          IM = run_deseq2("GSE99577", "IM", "IHKA_CA1_IM"),
                          CR = run_deseq2("GSE99577", "CR", "IHKA_CA1_CR")),
  "DG (GSE213393)" = list(AC = run_deseq2("GSE213393", "AC", "IHKA_DG_AC"),
                          IM = run_deseq2("GSE213393", "IM", "IHKA_DG_IM"),
                          CR = run_deseq2("GSE213393", "CR", "IHKA_DG_CR")))

load_lfc <- function(label) {
  df <- read.csv(file.path(deg_dir, paste0("DEG_", label, ".csv")))
  setNames(df$log2FoldChange, df$gene)
}
iaka_lfc <- list(AC = load_lfc("IAKA_AC"), IM = load_lfc("IAKA_IM"), CR = load_lfc("IAKA_CR"))

pearson_stats <- function(v1, v2, filter = NULL) {
  nm <- intersect(names(v1), names(v2))
  if (!is.null(filter)) nm <- intersect(nm, filter)
  v1 <- v1[nm]; v2 <- v2[nm]
  ok <- is.finite(v1) & is.finite(v2); v1 <- v1[ok]; v2 <- v2[ok]; n <- length(v1)
  r  <- cor(v1, v2)
  ci <- quantile(replicate(N_BOOT, { i <- sample(n, n, TRUE); cor(v1[i], v2[i]) }),
                 c(0.025, 0.975), names = FALSE)
  pv <- mean(abs(replicate(N_PERM, cor(v1, sample(v2)))) >= abs(r))   # 0 means p < 1/N_PERM
  c(r = r, ci_lo = ci[1], ci_hi = ci[2], pval = pv, n_genes = n)
}

rows <- list()
for (gs in c("Whole transcriptome", "Epileptogenesis gene set")) {
  gf <- if (gs == "Whole transcriptome") NULL else epi_genes
  for (ds in names(by_ds)) for (pa in names(iaka_lfc)) for (pb in names(by_ds[[ds]])) {
    s <- pearson_stats(iaka_lfc[[pa]], by_ds[[ds]][[pb]], gf)
    rows[[length(rows) + 1]] <- data.frame(gene_set = gs, dataset = ds, iaka = pa,
                                           comparison = pb, t(s))
  }
}
write.csv(bind_rows(rows), file.path(out_dir, "Pearson_correlation_by_dataset.csv"),
          row.names = FALSE, quote = FALSE)
cat("Saved: Pearson_correlation_by_dataset.csv\n")
