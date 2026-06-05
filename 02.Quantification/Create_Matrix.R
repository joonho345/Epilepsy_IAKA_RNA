library(dplyr)

# SET THESE PATHS
BASE_DIR  <- "/path/to/project"              # output root
HTSEQ_DIR <- "/path/to/HTseq/output"         # directory with *.htseq.count.txt
CSV_FILE  <- "/path/to/sample_sheet.csv"     # sample metadata (columns: Sample_ID, SRR)
GTF_PATH  <- "/path/to/Mus_musculus.GRCm39.112.gtf"

out_dir <- file.path(BASE_DIR, "02.Quantification")
dir.create(out_dir, showWarnings = FALSE)

# Sample mapping
meta    <- read.csv(CSV_FILE, check.names = FALSE, stringsAsFactors = FALSE)
mapping <- meta[grepl("^IAKA|^IASL", meta$Sample_ID), c("Sample_ID","SRR")]
mapping$Group <- ifelse(grepl("^IAKA", mapping$Sample_ID), "IAKA", "IASL")
mapping$Day   <- as.integer(regmatches(mapping$Sample_ID, regexpr("[0-9]+", mapping$Sample_ID)))
mapping <- mapping[order(mapping$Group, mapping$Day), ]

# Merge HTSeq count files
read_htseq <- function(srr)
  read.table(file.path(HTSEQ_DIR, paste0(srr, ".htseq.count.txt")),
             sep = "\t", col.names = c("Gene", srr))
merged <- Reduce(function(x, y) merge(x, y, by = "Gene"),
                 lapply(mapping$SRR, read_htseq))
merged <- merged[!grepl("^__", merged$Gene), ]
rownames(merged) <- merged$Gene; merged <- merged[, -1]

write.table(merged, file.path(out_dir, "merged_matrix_id.txt"),
            sep = "\t", quote = FALSE)

# gene_id → gene_name from GTF
gtf_raw  <- read.table(GTF_PATH, skip = 5, sep = "\t", quote = "",
                        comment.char = "")[, c(3, 9)]
gtf_gene <- gtf_raw[gtf_raw[, 1] == "gene", 2]
exfield  <- function(v9, field)
  sub(paste0('.*', field, ' "(\\S+)".*'), "\\1",
      regmatches(v9, regexpr(paste0(field, ' "\\S+"'), v9)))
gtf_map <- data.frame(gene_id   = exfield(gtf_gene, "gene_id"),
                      gene_name = exfield(gtf_gene, "gene_name"),
                      stringsAsFactors = FALSE)
gtf_map <- gtf_map[!is.na(gtf_map$gene_name) & !duplicated(gtf_map$gene_id), ]

merged$gene_id <- rownames(merged)
merged <- merge(merged, gtf_map, by = "gene_id")
merged <- merged[!is.na(merged$gene_name), setdiff(colnames(merged), "gene_id")]
samp_cols <- setdiff(colnames(merged), "gene_name")
merged[samp_cols] <- lapply(merged[samp_cols], as.numeric)
merged <- aggregate(. ~ gene_name, data = merged, FUN = sum, na.rm = TRUE)
rownames(merged) <- merged$gene_name; merged <- merged[, -1]

write.table(merged, file.path(out_dir, "merged_matrix.txt"),
            sep = "\t", quote = FALSE)

coldata <- mapping[, c("SRR","Sample_ID","Group","Day")]
coldata$Day <- paste0("Day", coldata$Day)
write.table(coldata, file.path(out_dir, "coldata.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Done:", nrow(merged), "genes x", ncol(merged), "samples\n")
