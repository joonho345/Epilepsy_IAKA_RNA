================================================================================
IAKA RNA-SEQ ANALYSIS PIPELINE
================================================================================

Time-course hippocampal transcriptomes after intra-amygdala kainate
injection in mice
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
    EpilepsyGeneSet.R           Fisher enrichment in gene sets (Fig. 3c, Supplementary Data 3)
  07.CrossDataset/
    coldata/                    Sample sheets for IAKA, IHKA and PILO (labels verified against GEO)
    Batch_correction.R          ComBat-seq for the pilocarpine datasets
    DEG_Correlation.R           DESeq2 + Pearson r, whole transcriptome (Fig. 4a)
                                (IHKA: raw counts, ~ Batch + Treatment)
    Correlation_EpiGenes.R      Pearson r, epileptogenesis gene set (Fig. 4b)
    IHKA_by_dataset.R           Pearson r for each IHKA dataset, CA1 and DG (Fig. 4c, d)
    Figure4_dotline.py          Fig. 4 dot-line plots (Python 3, matplotlib)

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


--------------------------------------------------------------------------------
REVISION 1 (2026-11) - CHANGES
--------------------------------------------------------------------------------
- 07.CrossDataset/coldata/coldata_ihka.txt: kainate and saline labels of 12 of the
  24 GSE213393 samples corrected according to the GEO records.
- Batch_correction.R / DEG_Correlation.R: IHKA is no longer ComBat-seq corrected.
  Dataset is included as a blocking factor in the DESeq2 design (~ Batch + Treatment).
  The hyperacute phase (GSE99577 only) uses ~ Treatment. PILO still uses ComBat-seq.
- IHKA_by_dataset.R (new): correlations with each IHKA dataset separately.
- EpilepsyGeneSet.R: background = genes with a valid DESeq2 test at each time point,
  combined test with genes counted once, category-level tests with BH adjustment.
- Figure4_dotline.py (new): Fig. 4 redrawn as dot-line plots with 95% CIs.
- Permutation p-values: 0 means that none of the 1,000 permuted correlations reached
  the observed value (p < 0.001).
