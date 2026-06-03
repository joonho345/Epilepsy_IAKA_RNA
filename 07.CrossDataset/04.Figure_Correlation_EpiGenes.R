library(dplyr)
set.seed(42)

BASE_DIR <- "/path/to/project"    # SET THIS
GO_PATH  <- "/path/to/GO.txt"     # SET THIS: gene list file for epileptogenesis gene set
N_BOOT <- 1000; N_PERM <- 1000

deg_dir <- file.path(BASE_DIR, "07.CrossDataset", "03.DEGCorrelation")

# Parse epileptogenesis gene set from GO.txt
lines     <- readLines(GO_PATH)
epi_genes <- unique(trimws(unlist(lapply(
  lines[grepl("^[A-Z][a-z].*\tn=", lines)],
  function(l) strsplit(sub(".*\tn=[0-9]+\t", "", l), ", ")[[1]]))))
cat("Epileptogenesis gene set:", length(epi_genes), "genes\n")

# Load log2FC vectors
load_lfc <- function(label) {
  df <- read.csv(file.path(deg_dir, paste0("DEG_", label, ".csv")))
  setNames(df$log2FoldChange, df$gene)
}
iaka_lfc <- list(AC = load_lfc("IAKA_AC"), IM = load_lfc("IAKA_IM"), CR = load_lfc("IAKA_CR"))
ihka_lfc <- list(HA = load_lfc("IHKA_HA"), AC = load_lfc("IHKA_AC"),
                 IM = load_lfc("IHKA_IM"), CR = load_lfc("IHKA_CR"))
pilo_lfc <- list(HA = load_lfc("PILO_HA"), AC = load_lfc("PILO_AC"),
                 IM = load_lfc("PILO_IM"), CR = load_lfc("PILO_CR"))

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

cor_table <- function(lfc_a, lfc_b, model, filter = NULL) {
  bind_rows(lapply(names(lfc_a), function(pa)
    bind_rows(lapply(names(lfc_b), function(pb) {
      s <- pearson_stats(lfc_a[[pa]], lfc_b[[pb]], filter)
      cat(sprintf("  IAKA-%s vs %s-%s  r=%.3f [%.3f,%.3f]  p=%s\n",
                  pa, model, pb, s["r"], s["ci_lo"], s["ci_hi"],
                  if (s["pval"] < 0.001) "<0.001" else sprintf("%.3f", s["pval"])))
      data.frame(iaka = pa, comparison = pb, model = model,
                 r = s["r"], ci_lo = s["ci_lo"], ci_hi = s["ci_hi"],
                 pval = s["pval"], n_genes = s["n"])
    }))))
}

cat("=== IAKA vs IHKA (epileptogenesis gene set) ===\n")
cor_ihka <- cor_table(iaka_lfc, ihka_lfc, "IHKA", epi_genes)
cat("=== IAKA vs PILO (epileptogenesis gene set) ===\n")
cor_pilo <- cor_table(iaka_lfc, pilo_lfc, "PILO", epi_genes)

write.csv(bind_rows(cor_ihka, cor_pilo),
          file.path(deg_dir, "Pearson_correlation_summary_EpiGenes.csv"),
          row.names = FALSE, quote = FALSE)
cat("Saved: Pearson_correlation_summary_EpiGenes.csv\n")
