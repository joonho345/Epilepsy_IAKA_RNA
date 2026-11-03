library(DESeq2)

BASE_DIR <- "/path/to/project"   # SET THIS

quant_dir <- file.path(BASE_DIR, "02.Quantification")
out_dir   <- file.path(BASE_DIR, "05.DESeq")
dir.create(out_dir, showWarnings = FALSE)

coldata    <- read.table(file.path(quant_dir, "coldata.txt"),
                         sep = "\t", header = TRUE, row.names = 1,
                         stringsAsFactors = FALSE)
counts_raw <- as.matrix(read.table(file.path(quant_dir, "merged_matrix.txt"),
                                    sep = "\t", header = TRUE, row.names = 1,
                                    check.names = FALSE))
common     <- intersect(rownames(coldata), colnames(counts_raw))
counts_raw <- counts_raw[, common]; coldata <- coldata[common, ]

for (day in c("Day3","Day7","Day14")) {
  sel     <- rownames(coldata)[coldata$Day == day]
  cnt_sub <- round(counts_raw[, sel])
  col_sub <- coldata[sel, ]
  col_sub$Group <- factor(col_sub$Group, levels = c("IASL","IAKA"))

  dds <- DESeqDataSetFromMatrix(cnt_sub, col_sub, ~ Group)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("Group","IAKA","IASL")))
  res$gene      <- rownames(res)
  res$Direction <- "NS"
  res$Direction[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange >  1] <- "Up"
  res$Direction[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1] <- "Down"

  write.csv(res[order(res$padj, na.last = TRUE), ],
            file.path(out_dir, paste0("DESeq2_IAKA_vs_IASL_", day, "_all.csv")),
            row.names = FALSE, quote = FALSE)
  write.csv(res[res$Direction != "NS", ],
            file.path(out_dir, paste0("DEG_IAKA_vs_IASL_", day, ".csv")),
            row.names = FALSE, quote = FALSE)
  cat(sprintf("%s — Up:%d  Down:%d\n", day,
              sum(res$Direction == "Up"), sum(res$Direction == "Down")))
}
