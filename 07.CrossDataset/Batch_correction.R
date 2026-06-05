library(sva)

BASE_DIR <- "/path/to/project"   # SET THIS

quant_dir <- file.path(BASE_DIR, "07.CrossDataset", "01.Quantification")
out_dir   <- file.path(BASE_DIR, "07.CrossDataset", "02.BatchCorrection")
dir.create(out_dir, showWarnings = FALSE)

load_model <- function(m) {
  col <- read.table(file.path(quant_dir, paste0("coldata_", m, ".txt")),
                    sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  cnt <- data.matrix(read.table(file.path(quant_dir, paste0("merged_matrix_", m, ".txt")),
                                  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  k <- intersect(rownames(col), colnames(cnt))
  list(col = col[k, ], cnt = cnt[, k])
}

save_counts <- function(cnt, label)
  write.table(cnt, file.path(out_dir, paste0("counts_corrected_", label, ".txt")),
              sep = "\t", quote = FALSE)

# IAKA — single batch, no correction
iaka <- load_model("iaka")
save_counts(iaka$cnt, "iaka")

# IHKA — 2 batches (GSE99577, GSE213393), Treatment as protected covariate
ihka <- load_model("ihka")
save_counts(ComBat_seq(ihka$cnt, ihka$col$Batch,
                        model.matrix(~ ihka$col$Treatment)[, -1, drop = FALSE]), "ihka")

# PILO — 4 batches, Treatment as protected covariate
pilo <- load_model("pilo")
save_counts(ComBat_seq(pilo$cnt, pilo$col$Batch,
                        model.matrix(~ pilo$col$Treatment)[, -1, drop = FALSE]), "pilo")

cat("Done: counts_corrected_iaka/ihka/pilo.txt\n")
