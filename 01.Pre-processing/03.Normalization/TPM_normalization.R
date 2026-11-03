library(dplyr); library(stringr)

BASE_DIR <- "/path/to/project"                         # SET THIS
GTF_PATH <- "/path/to/Mus_musculus.GRCm39.112.gtf"    # SET THIS

quant_dir <- file.path(BASE_DIR, "02.Quantification")
out_dir   <- file.path(BASE_DIR, "03.Normalization")
dir.create(out_dir, showWarnings = FALSE)

# Non-overlapping exon lengths from GTF
gtf <- read.table(GTF_PATH, skip = 5, sep = "\t", quote = "", comment.char = "",
                  col.names = c("chr","source","feature","start","end",
                                "score","strand","frame","attribute"))
gene_lengths <- gtf %>%
  filter(feature == "exon") %>%
  mutate(gn = str_replace_all(str_extract(attribute, 'gene_name "[^"]+"'),
                               'gene_name |"', "")) %>%
  filter(!is.na(gn)) %>%
  group_by(gn) %>% arrange(start, end) %>%
  mutate(adj = ifelse(row_number() == 1, start, pmax(start, lag(end) + 1)),
         len = ifelse(adj <= end, end - adj + 1, 0)) %>%
  summarise(gene_length = sum(len, na.rm = TRUE), .groups = "drop")

counts  <- as.matrix(read.table(file.path(quant_dir, "merged_matrix.txt"),
                                  sep = "\t", header = TRUE, row.names = 1,
                                  check.names = FALSE))
len_kb  <- gene_lengths$gene_length[match(rownames(counts), gene_lengths$gn)] / 1000
valid   <- !is.na(len_kb) & len_kb > 0
rpk     <- counts[valid, ] / len_kb[valid]
tpm     <- t(t(rpk) / colSums(rpk, na.rm = TRUE)) * 1e6

write.table(tpm, file.path(out_dir, "merged_matrix_TPM.txt"),
            sep = "\t", quote = FALSE)
cat("TPM saved:", nrow(tpm), "genes x", ncol(tpm), "samples\n")
