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

# Revision 1:
#  - IAKA: single sequencing batch, no correction.
#  - IHKA: no count-level correction. Both datasets contain kainate and saline animals
#          in every phase they cover, so dataset is blocked in the DESeq2 design instead
#          (DEG_Correlation.R). Sample labels follow the GEO records (coldata/coldata_ihka.txt).
#  - PILO: cases and controls of some phases come from different datasets, so count-level
#          batch effects are corrected with ComBat-seq, with treatment as a protected covariate.
pilo <- load_model("pilo")
write.table(ComBat_seq(pilo$cnt, pilo$col$Batch,
                       model.matrix(~ pilo$col$Treatment)[, -1, drop = FALSE]),
            file.path(out_dir, "counts_corrected_pilo.txt"), sep = "\t", quote = FALSE)
cat("Done: counts_corrected_pilo.txt\n")
