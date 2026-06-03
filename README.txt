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
  settings.sh             Path settings for shell scripts (edit DATA_DIR, REF_DIR)
  00.Geneset.R            Epileptogenesis GO gene sets (9 categories, mouse orthologs)
  01.Pre-processing/      FASTP → STAR index → STAR align
  02.Quantification/      HTSeq count → count matrix + coldata
  03.Normalization/       TPM normalization
  04.Clustering/          PCA + hierarchical clustering
  05.DESeq/               DESeq2 differential expression (Fig. 3)
  06.EpilepsyGeneSet/     Gene set enrichment — Fisher's exact test (Fig. 3)
  07.CrossDataset/        Batch correction + Pearson r vs IHKA/PILO (Fig. 4a, 4b)

--------------------------------------------------------------------------------
SETUP
--------------------------------------------------------------------------------
1. Edit settings.sh: set DATA_DIR (pipeline data) and REF_DIR (GRCm39 reference)
2. Edit each R script: set BASE_DIR to your project root
   - 00.Geneset.R: also set GENESET_DIR (BioMart GO *.txt files)
   - 07.CrossDataset/01.Quantification.py: also set NATCOM_MTX and ANIMAL_SHEET
     (NatCom matrix: Kim et al., Nat Commun, 2025)
   - 07.CrossDataset/04.Figure_Correlation_EpiGenes.R: also set GO_PATH
3. Reference files required (Ensembl release 112):
     Mus_musculus.GRCm39.dna.primary_assembly.fa
     Mus_musculus.GRCm39.112.gtf

--------------------------------------------------------------------------------
EXECUTION ORDER
--------------------------------------------------------------------------------
  0.  00.Geneset.R
  1.  bash 01.Pre-processing/01.fastp_paired.sh  <ID> settings.sh
  2.  bash 01.Pre-processing/02.STAR_index.sh    settings.sh
  3.  bash 01.Pre-processing/03.STAR_align.sh    <ID> settings.sh
  4.  bash 02.Quantification/01.HTseq.sh         <ID> settings.sh
  5.  Rscript 02.Quantification/02.Create_Matrix.R
  6.  Rscript 03.Normalization/01.TPM_normalization.R
  7.  Rscript 04.Clustering/01.Clustering.R
  8.  Rscript 05.DESeq/01.DESeq_Volcano.R
  9.  Rscript 06.EpilepsyGeneSet/01.EpilepsyGeneSet.R
  10. python3 07.CrossDataset/01.Quantification.py
  11. Rscript 07.CrossDataset/02.Batch_correction.R
  12. Rscript 07.CrossDataset/03.DEG_Correlation.R        (Fig. 4a)
  13. Rscript 07.CrossDataset/04.Figure_Correlation_EpiGenes.R  (Fig. 4b)

  Steps 1, 3, 4 are run once per sample. Steps 2, 5–13 are run once.

--------------------------------------------------------------------------------
SOFTWARE
--------------------------------------------------------------------------------
  FASTP v0.21.0  |  STAR v2.7.10a  |  HTSeq v0.12.4  |  samtools v1.17
  R: DESeq2 v1.38.3 | sva | dplyr | stringr
  Python 3.8+ (standard library only)

================================================================================
