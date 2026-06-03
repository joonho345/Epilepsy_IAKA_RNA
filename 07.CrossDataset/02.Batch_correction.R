library(dplyr); library(sva); library(stringr)
set.seed(42)

BASE_DIR <- "/path/to/project"                         # SET THIS
GTF_PATH <- "/path/to/Mus_musculus.GRCm39.112.gtf"    # SET THIS

quant_dir <- file.path(BASE_DIR, "07.CrossDataset", "01.Quantification")
out_dir   <- file.path(BASE_DIR, "07.CrossDataset", "02.BatchCorrection")
dir.create(out_dir, showWarnings = FALSE)

# Gene lengths (non-overlapping exons)
gtf <- read.table(GTF_PATH, skip = 5, sep = "\t", quote = "", comment.char = "",
                  col.names = c("chr","source","feature","start","end",
                                "score","strand","frame","attribute"))
gene_lengths <- gtf %>%
  filter(feature == "exon") %>%
  mutate(gn = str_replace_all(str_extract(attribute, 'gene_name "[^"]+"'),
                               'gene_name |"', "")) %>%
  filter(!is.na(gn)) %>% group_by(gn) %>% arrange(start, end) %>%
  mutate(adj = ifelse(row_number() == 1, start, pmax(start, lag(end) + 1)),
         len = ifelse(adj <= end, end - adj + 1, 0)) %>%
  summarise(gene_length = sum(len, na.rm = TRUE), .groups = "drop")

load_model <- function(m) {
  col <- read.table(file.path(quant_dir, paste0("coldata_", m, ".txt")),
                    sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  cnt <- data.matrix(read.table(file.path(quant_dir, paste0("merged_matrix_", m, ".txt")),
                                  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  k <- intersect(rownames(col), colnames(cnt))
  list(col = col[k, ], cnt = cnt[, k])
}

counts_to_tpm <- function(cnt) {
  lk  <- gene_lengths$gene_length[match(rownames(cnt), gene_lengths$gn)] / 1000
  ok  <- !is.na(lk) & lk > 0
  rpk <- cnt[ok, ] / lk[ok]
  t(t(rpk) / colSums(rpk, na.rm = TRUE)) * 1e6
}

save_outputs <- function(cnt, label) {
  write.table(cnt, file.path(out_dir, paste0("counts_corrected_", label, ".txt")),
              sep = "\t", quote = FALSE)
  write.table(counts_to_tpm(cnt), file.path(out_dir, paste0("TPM_corrected_", label, ".txt")),
              sep = "\t", quote = FALSE)
  cat(label, "saved\n")
}

# IAKA — single batch, no correction
iaka <- load_model("iaka")
save_outputs(iaka$cnt, "iaka")

# IHKA — 2 batches, ComBat_seq
ihka <- load_model("ihka")
save_outputs(ComBat_seq(ihka$cnt, ihka$col$Batch,
                         model.matrix(~ ihka$col$Treatment)[, -1, drop = FALSE]), "ihka")

# PILO — 4 batches, ComBat_seq
pilo <- load_model("pilo")
save_outputs(ComBat_seq(pilo$cnt, pilo$col$Batch,
                         model.matrix(~ pilo$col$Treatment)[, -1, drop = FALSE]), "pilo")
