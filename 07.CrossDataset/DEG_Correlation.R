library(DESeq2); library(dplyr)
set.seed(42)

BASE_DIR <- "/path/to/project"   # SET THIS
N_BOOT <- 1000; N_PERM <- 1000

quant_dir <- file.path(BASE_DIR, "07.CrossDataset", "01.Quantification")
batch_dir <- file.path(BASE_DIR, "07.CrossDataset", "02.BatchCorrection")
out_dir   <- file.path(BASE_DIR, "07.CrossDataset", "03.DEGCorrelation")
dir.create(out_dir, showWarnings = FALSE)

load_data <- function(m, corrected = FALSE) {
  col <- read.table(file.path(quant_dir, paste0("coldata_", m, ".txt")),
                    sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  cnt <- round(data.matrix(read.table(
    file.path(if (corrected) batch_dir else quant_dir,
              paste0(if (corrected) "counts_corrected_" else "merged_matrix_", m, ".txt")),
    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)))
  k <- intersect(rownames(col), colnames(cnt))
  list(col = col[k, ], cnt = cnt[, k])
}

run_deseq2 <- function(cnt, col, groups, case_t, ctrl_t, label) {
  sel   <- rownames(col)[col$Group %in% groups]
  col_s <- col[sel, ]; col_s$Treatment <- factor(col_s$Treatment, c(ctrl_t, case_t))
  dds   <- DESeqDataSetFromMatrix(cnt[, sel], col_s, ~ Treatment)
  n_min <- min(table(col_s$Treatment))
  dds   <- dds[rowSums(counts(dds) >= 10) >= max(2L, n_min %/% 2), ]
  dds   <- DESeq(dds, quiet = TRUE)
  res   <- as.data.frame(results(dds, contrast = c("Treatment", case_t, ctrl_t)))
  res$gene <- rownames(res)
  write.csv(res[order(res$padj, na.last = TRUE), ],
            file.path(out_dir, paste0("DEG_", label, ".csv")),
            row.names = FALSE, quote = FALSE)
  setNames(res$log2FoldChange, rownames(res))
}

pearson_stats <- function(v1, v2, filter = NULL) {
  nm <- intersect(names(v1), names(v2))
  if (!is.null(filter)) nm <- intersect(nm, filter)
  v1 <- v1[nm]; v2 <- v2[nm]
  ok <- is.finite(v1) & is.finite(v2); v1 <- v1[ok]; v2 <- v2[ok]; n <- length(v1)
  r  <- cor(v1, v2)
  ci <- quantile(replicate(N_BOOT, { i <- sample(n, n, TRUE); cor(v1[i], v2[i]) }),
                 c(0.025, 0.975), names = FALSE)
  pv <- mean(abs(replicate(N_PERM, cor(v1, sample(v2)))) >= abs(r))
  c(r = r, ci_lo = ci[1], ci_hi = ci[2], pval = pv, n = n)
}

cor_table <- function(lfc_a, lfc_b, model, filter = NULL)
  bind_rows(lapply(names(lfc_a), function(pa)
    bind_rows(lapply(names(lfc_b), function(pb) {
      s <- pearson_stats(lfc_a[[pa]], lfc_b[[pb]], filter)
      as.data.frame(as.list(c(iaka = pa, comparison = pb, model = model, s)))
    }))))

# ── DESeq2 ────────────────────────────────────────────────────────────────────
iaka <- load_data("iaka")
iaka_lfc <- list(
  AC = run_deseq2(iaka$cnt, iaka$col, c("IAKA_AC","IASL_AC"), "KA","Saline","IAKA_AC"),
  IM = run_deseq2(iaka$cnt, iaka$col, c("IAKA_IM","IASL_IM"), "KA","Saline","IAKA_IM"),
  CR = run_deseq2(iaka$cnt, iaka$col, c("IAKA_CR","IASL_CR"), "KA","Saline","IAKA_CR"))

ihka <- load_data("ihka", corrected = TRUE)
ihka_lfc <- list(
  HA = run_deseq2(ihka$cnt, ihka$col, c("IHKA_HA","IHSL_HA"), "KA","Saline","IHKA_HA"),
  AC = run_deseq2(ihka$cnt, ihka$col, c("IHKA_AC","IHSL_AC"), "KA","Saline","IHKA_AC"),
  IM = run_deseq2(ihka$cnt, ihka$col, c("IHKA_IM","IHSL_IM"), "KA","Saline","IHKA_IM"),
  CR = run_deseq2(ihka$cnt, ihka$col, c("IHKA_CR","IHSL_CR"), "KA","Saline","IHKA_CR"))

pilo <- load_data("pilo", corrected = TRUE)
pilo_lfc <- list(
  HA = run_deseq2(pilo$cnt, pilo$col, c("PILO_HA","PSLC_HA","PSLC_GSE72402"), "PILO","Saline","PILO_HA"),
  AC = run_deseq2(pilo$cnt, pilo$col, c("PILO_AC","PSLC_AC"),                 "PILO","Saline","PILO_AC"),
  IM = run_deseq2(pilo$cnt, pilo$col, c("PILO_IM","PSLC_GSE72402"),           "PILO","Saline","PILO_IM"),
  CR = run_deseq2(pilo$cnt, pilo$col, c("PILO_CR","PSLC_GSE72402"),           "PILO","Saline","PILO_CR"))

# ── Pearson r (whole transcriptome) ───────────────────────────────────────────
cor_all <- bind_rows(cor_table(iaka_lfc, ihka_lfc, "IHKA"),
                     cor_table(iaka_lfc, pilo_lfc, "PILO"))
write.csv(cor_all, file.path(out_dir, "Pearson_correlation_summary.csv"),
          row.names = FALSE, quote = FALSE)
cat("Saved: Pearson_correlation_summary.csv\n")
