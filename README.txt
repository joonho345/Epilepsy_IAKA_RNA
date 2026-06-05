================================================================================
IAKA RNA-SEQ ANALYSIS PIPELINE
================================================================================

Hippocampal bulk RNA-seq across epileptogenic phases in the
intra-amygdala kainate mouse model
Kim J, Kim B, Cho S, Kim W-J, Kim S — Scientific Data (2026)

  GEO: GSE319769  |  SRA: PRJNA1425071
  36 samples (IAKA n=18, IASL n=18)  |  Illumina NovaSeq X Plus, 101 bp PE

--------------------------------------------------------------------------------
STRUCTURE
--------------------------------------------------------------------------------
  settings.sh                   Path settings for shell scripts
  Geneset.R                     Epileptogenesis GO gene sets (9 categories)
  01.Pre-processing/
    fastp_paired.sh             QC trimming (FASTP v0.21.0)
    STAR_index.sh               Genome index (STAR v2.7.10a, GRCm39)
    STAR_align.sh               Read alignment
  02.Quantification/
    HTseq.sh                    Gene counting (HTSeq v0.12.4)
    Create_Matrix.R             Merge counts + gene name conversion
  03.Normalization/
    TPM_normalization.R         TPM calculation
  04.Clustering/
    Clustering.R                PCA + hierarchical clustering
  05.DESeq/
    DESeq.R                     DESeq2 (IAKA vs IASL, Fig. 3)
  06.EpilepsyGeneSet/
    EpilepsyGeneSet.R           Fisher enrichment in gene sets (Fig. 3)
  07.CrossDataset/
    Batch_correction.R          ComBat-seq per model
    DEG_Correlation.R           DESeq2 + Pearson r, whole transcriptome (Fig. 4a)
    Correlation_EpiGenes.R      Pearson r, epileptogenesis gene set (Fig. 4b)

--------------------------------------------------------------------------------
SETUP
--------------------------------------------------------------------------------
1. Edit settings.sh: set DATA_DIR (pipeline data) and REF_DIR (GRCm39 reference)
2. Edit each R script: set BASE_DIR to your project root
   - Geneset.R: also set GENESET_DIR (BioMart GO *.txt files)
   - Correlation_EpiGenes.R: also set GO_PATH
3. Reference files required (Ensembl release 112):
     Mus_musculus.GRCm39.dna.primary_assembly.fa
     Mus_musculus.GRCm39.112.gtf

