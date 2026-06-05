library(dplyr)

GENESET_DIR <- "/path/to/BioMart_GO"  # SET THIS: directory with BioMart GO *.txt files

upper_categories <- c("IC","NT","SP","GG","NG","MF","NI","ND","TNF")

M_gene_sets <- setNames(lapply(upper_categories, function(cat) {
  files <- list.files(GENESET_DIR, pattern = cat, full.names = TRUE)
  unique(unlist(lapply(files, function(f) {
    dat <- read.table(f, header = TRUE, sep = "\t", fill = TRUE,
                      stringsAsFactors = FALSE, quote = "")
    dat$Mouse.gene.name[dat$Mouse.homology.type == "ortholog_one2one" &
                        dat$Mouse.gene.name != ""]
  })))
}), upper_categories)

M_gene_set_W <- unique(unlist(M_gene_sets))

for (cat in names(M_gene_sets))
  cat(sprintf("  %-4s: %d genes\n", cat, length(M_gene_sets[[cat]])))
cat(sprintf("  Union: %d genes\n", length(M_gene_set_W)))
