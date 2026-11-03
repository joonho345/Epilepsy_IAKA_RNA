# Gene-set enrichment of DEGs in the nine epileptogenesis-related categories (Fig. 3c,
# Supplementary Data 3).
# Revision 1 (Reviewer 1, Minor comment 5):
#  - Background: all genes with a valid DESeq2 test (non-missing p-value) at each time point,
#    taken from DESeq2_IAKA_vs_IASL_<day>_all.csv (05.DESeq/DESeq.R).
#  - DEGs: padj < 0.05 and |log2FC| > 1 (Direction column from DESeq.R).
#  - Combined analysis: genes in more than one category are counted once.
#  - Two-sided Fisher's exact tests, Benjamini-Hochberg adjustment within each table.
library(dplyr)

BASE_DIR       <- "/path/to/project"      # SET THIS
GENESET_SCRIPT <- "/path/to/Geneset.R"    # SET THIS
source(GENESET_SCRIPT)                     # provides M_gene_sets (list of 9 categories)

deseq_dir <- file.path(BASE_DIR, "05.DESeq")
out_dir   <- file.path(BASE_DIR, "06.EpilepsyGeneSet")
dir.create(out_dir, showWarnings = FALSE)

union_set <- unique(unlist(M_gene_sets))
fisher_row <- function(deg, set, bg) {
  set <- intersect(set, bg); a <- length(intersect(deg, set))
  mat <- matrix(c(a, length(deg) - a, length(set) - a,
                  length(bg) - length(set) - length(deg) + a), 2)
  ft  <- fisher.test(mat)
  data.frame(background = length(bg), set_in_background = length(set),
             pct_background = 100 * length(set) / length(bg), n_DEG = length(deg),
             DEG_in_set = a, pct_DEG = 100 * a / length(deg),
             odds_ratio = unname(ft$estimate), p = ft$p.value)
}

union_rows <- list(); cat_rows <- list()
for (day in c("Day3", "Day7", "Day14")) {
  res <- read.csv(file.path(deseq_dir, paste0("DESeq2_IAKA_vs_IASL_", day, "_all.csv")))
  bg  <- res$gene[!is.na(res$pvalue)]
  for (dir in c("Up", "Down")) {
    deg <- res$gene[res$Direction == dir]
    union_rows[[length(union_rows) + 1]] <- cbind(Day = day, Direction = dir,
                                                  fisher_row(deg, union_set, bg))
    for (cat in names(M_gene_sets))
      cat_rows[[length(cat_rows) + 1]] <- cbind(Day = day, Direction = dir, Category = cat,
                                                fisher_row(deg, M_gene_sets[[cat]], bg))
  }
}
union_df <- bind_rows(union_rows); union_df$padj_BH <- p.adjust(union_df$p, "BH")
cat_df   <- bind_rows(cat_rows);   cat_df$padj_BH   <- p.adjust(cat_df$p, "BH")
write.csv(union_df, file.path(out_dir, "GeneSet_Enrichment_combined.csv"), row.names = FALSE, quote = FALSE)
write.csv(cat_df,   file.path(out_dir, "GeneSet_Enrichment_by_category.csv"), row.names = FALSE, quote = FALSE)
cat("Done: GeneSet_Enrichment_combined.csv, GeneSet_Enrichment_by_category.csv\n")
