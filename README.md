# TCGA Multi-Omics Data Processing Pipeline

This repository contains organized TCGA multi-omics data and processing scripts for immune-focused analysis.

## Directory Structure

```
├── data/
│   ├── raw_data/
│   │   ├── rnaseq/
│   │   │   ├── tcga_acc_rnaseq_rsem_log2.cct
│   │   │   └── tcga_acc_rnaseq_rsem_log2_compressed.cct.gz
│   │   ├── mirna/
│   │   │   └── tcga_acc_mirna_rpm_log2.cct.gz
│   │   ├── methylation/
│   │   │   └── tcga_acc_methylation_450k.cct.gz
│   │   ├── rppa/
│   │   │   ├── tcga_acc_rppa_gene_level.cct
│   │   │   └── tcga_acc_rppa_analyte_level.cct
│   │   └── copy_number/
│   │       ├── tcga_acc_copy_number_gene_level.cct.gz
│   │       └── tcga_acc_copy_number_focal_level.cct
│   └── processed_data/
│       └── (processed output files)
└── scripts/
    ├── process_rnaseq.py
    ├── process_mirna.py (to be created)
    ├── process_methylation.py (to be created)
    ├── process_rppa.py (to be created)
    ├── process_copy_number.py (to be created)
    └── integrate_multiomics.py (to be created)
```

## Data Types

### 1. RNA-seq Data (Gene Expression)
- **File**: `tcga_acc_rnaseq_rsem_log2.cct`
- **Description**: Log2-transformed RSEM gene expression values
- **Samples**: 79 TCGA samples
- **Features**: ~19,000 genes

### 2. miRNA-seq Data (MicroRNA Expression)
- **File**: `tcga_acc_mirna_rpm_log2.cct.gz`
- **Description**: Log2-transformed RPM (Reads Per Million) miRNA expression values
- **Platform**: Illumina HiSeq miRNA sequencing

### 3. Methylation Data (DNA Methylation)
- **File**: `tcga_acc_methylation_450k.cct.gz`
- **Description**: DNA methylation beta values from Illumina 450K arrays
- **Platform**: Illumina Infinium HumanMethylation450 BeadChip

### 4. RPPA Data (Protein Expression)
- **Files**: 
  - `tcga_acc_rppa_gene_level.cct` (gene-level protein data)
  - `tcga_acc_rppa_analyte_level.cct` (analyte-level with phosphorylation states)
- **Description**: Reverse Phase Protein Array quantification
- **Features**: ~200 proteins and phospho-proteins

### 5. Copy Number Data (Somatic Copy Number Alterations)
- **Files**:
  - `tcga_acc_copy_number_gene_level.cct.gz` (gene-level GISTIC2 scores)
  - `tcga_acc_copy_number_focal_level.cct` (focal alteration regions)
- **Description**: GISTIC2-processed copy number alteration data
- **Platform**: Affymetrix SNP 6.0 arrays

## Processing Scripts

### Current Scripts
- **process_rnaseq.py**: Extracts immune marker genes from RNA-seq data and merges with WSI metadata

### Available Scripts
- **process_rnaseq.py**: Extracts immune marker genes from RNA-seq data and merges with WSI metadata
- **process_mirna.py**: Extracts immune-related miRNAs (checkpoint regulation, T cell function, etc.)
- **process_methylation.py**: Focuses on immune gene promoter methylation patterns
- **process_rppa.py**: Extracts immune pathway proteins and calculates pathway scores
- **process_copy_number.py**: Focuses on immune-related genomic regions and calculates instability metrics
- **integrate_multiomics.py**: Master integration script that runs all processing and creates unified dataset

## File Naming Convention

All data files use the `.cct` extension as they contain transcriptomics and related omics data. The naming pattern is:
`tcga_acc_[datatype]_[processing_method].cct[.gz]`

## Usage

### Option 1: Run Master Integration Script (Recommended)
```bash
python scripts/integrate_multiomics.py
```
This will automatically run all individual processing scripts and create an integrated multi-omics dataset.

### Option 2: Run Individual Scripts
```bash
python scripts/process_rnaseq.py
python scripts/process_mirna.py
python scripts/process_methylation.py
python scripts/process_rppa.py
python scripts/process_copy_number.py
```

### Output Files
All processed outputs are saved to `data/processed_data/`:
- Individual data type files: `[datatype]_immune_markers_with_metadata.csv`
- Integrated dataset: `tcga_acc_multiomics_integrated.csv`
- Summary report: `multiomics_integration_report.md`

## Key Features

### Multi-Omics Integration
- **Immune-focused analysis**: All scripts focus on immune-related genes, proteins, and pathways
- **Standardized processing**: Consistent sample ID handling and metadata integration
- **WSI integration**: Merges with whole slide image (WSI) metadata when available
- **Quality control**: Handles missing data and provides summary statistics

### Immune Markers Analyzed
- **RNA-seq**: CD274, CTLA4, PDCD1, LAG3, TIGIT, CD8A, FOXP3, GZMB
- **miRNA**: Let-7 family, miR-155, miR-146a, miR-21, miR-200 family, checkpoint regulators
- **Methylation**: Immune gene promoters, CpG sites in immune pathways
- **RPPA**: Checkpoint proteins, cytokines, transcription factors, kinase pathways
- **Copy Number**: Immune gene regions, MHC locus, checkpoint gene loci

### Analysis Outputs
- **Individual data matrices**: Immune markers × samples for each data type
- **Pathway scores**: Calculated for key immune pathways (T cell, checkpoint, etc.)
- **Multi-omics scores**: Integrated immune scores across all data types
- **Sample classification**: Hot/Cold/Intermediate immune phenotypes
- **Summary statistics**: Comprehensive analysis reports

## Dependencies

- pandas
- numpy
- json (for metadata handling)
- gzip (for compressed file handling)
- subprocess (for script orchestration)
- pathlib (for file path handling)
