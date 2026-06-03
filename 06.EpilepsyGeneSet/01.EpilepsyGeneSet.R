library(dplyr)

BASE_DIR       <- "/path/to/project"      # SET THIS
GENESET_SCRIPT <- "/path/to/00.Geneset.R" # SET THIS

source(GENESET_SCRIPT)

tpm_dir   <- file.path(BASE_DIR, "03.Normalization")
deseq_dir <- file.path(BASE_DIR, "05.DESeq")
out_dir   <- file.path(BASE_DIR, "06.EpilepsyGeneSet")
dir.create(out_dir, showWarnings = FALSE)

tpm <- as.matrix(read.table(file.path(tpm_dir, "merged_matrix_TPM.txt"),
                              sep = "\t", header = TRUE, row.names = 1,
                              check.names = FALSE))
bg  <- rownames(tpm)[rowSums(tpm > 0) >= 1]

degs <- bind_rows(lapply(c("Day3","Day7","Day14"), function(d)
  read.csv(file.path(deseq_dir, paste0("DEG_IAKA_vs_IASL_", d, ".csv"))) %>%
    filter(padj < 0.05) %>%
    mutate(Day = d, Direction = ifelse(log2FoldChange > 0, "Up", "Down"))))

# Fisher's exact test: 9 categories x 3 phases x 2 directions
enrich_df <- bind_rows(lapply(c("Day3","Day7","Day14"), function(day)
  bind_rows(lapply(c("Up","Down"), function(dir) {
    deg_sub <- intersect(degs$gene[degs$Day == day & degs$Direction == dir], bg)
    n_bg <- length(bg); n_deg <- length(deg_sub)
    bind_rows(lapply(names(M_gene_sets), function(cat) {
      cg <- intersect(M_gene_sets[[cat]], bg); a <- length(intersect(deg_sub, cg))
      mat <- matrix(c(a, n_deg - a, length(cg) - a,
                      n_bg - a - (length(cg) - a) - (n_deg - a)), 2)
      data.frame(Day = day, Direction = dir, Category = cat, n_overlap = a,
                 fold_enrich = if (length(cg) > 0 && n_deg > 0)
                   (a / length(cg)) / (n_deg / n_bg) else NA,
                 pval = fisher.test(mat)$p.value)
    }))
  }))))

enrich_df$padj <- p.adjust(enrich_df$pval, method = "BH")
write.csv(enrich_df[order(enrich_df$Day, enrich_df$Direction, enrich_df$padj), ],
          file.path(out_dir, "GeneSet_Enrichment.csv"), row.names = FALSE, quote = FALSE)

coverage <- degs %>%
  mutate(in_set = gene %in% M_gene_set_W) %>%
  group_by(Day, Direction) %>%
  summarise(total = n(), in_set = sum(in_set),
            pct = round(100 * in_set / n(), 1), .groups = "drop")
write.csv(coverage, file.path(out_dir, "GeneSet_Coverage.csv"),
          row.names = FALSE, quote = FALSE)
cat("Done: GeneSet_Enrichment.csv, GeneSet_Coverage.csv\n")
