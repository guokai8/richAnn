# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

richAnn is a Python package for gene set enrichment analysis with publication-quality visualizations. It provides GO enrichment, KEGG pathway analysis, and GSEA with permutation testing, inspired by the R package richR.

## Common Commands

```bash
# Install for development
pip install -e ".[dev]"

# Run all tests
pytest

# Run specific test file
pytest tests/test_enrichment.py

# Run with verbose output
pytest -v

# Run with coverage
pytest --cov=richAnn

# Code formatting
black richAnn/

# Linting
flake8 richAnn/
```

## Architecture

### Core Module Structure

- **`richAnn/core.py`**: `EnrichResult` class - the central data structure returned by all enrichment functions. Wraps a pandas DataFrame with methods for filtering (`.filter()`), getting top results (`.top()`), extracting gene-term details (`.detail()`), and export (`.to_csv()`, `.to_excel()`).

- **`richAnn/enrichment.py`**: Main enrichment functions:
  - `richGO()` - GO enrichment using hypergeometric test
  - `richKEGG()` - KEGG pathway enrichment
  - `richGSEA()` - Gene Set Enrichment Analysis with permutation testing
  - Internal helper `_convert_padj_method()` converts R-style method names (BH, BY) to statsmodels format (fdr_bh, fdr_by)

- **`richAnn/clustering.py`**: `richCluster()` - clusters enrichment terms using kappa statistics based on gene overlap

- **`richAnn/vis.py`**: Visualization functions (all return matplotlib Figure objects):
  - `ggbar()`, `ggdot()` - bar/dot plots for enrichment results
  - `ggnetplot()` - term-gene bipartite network
  - `ggnetwork()` - term-term similarity network (Jaccard)
  - `ggnetmap()` - combined network from multiple databases
  - `ggheatmap()`, `comparedot()` - cross-sample comparison plots
  - `ggGSEA()` - GSEA enrichment score plot

- **`richAnn/data.py`**: Data loading utilities for GMT, GAF, and KEGG mapping files

- **`richAnn/utils.py`**: Internal helpers for gene normalization, annotation validation, effect size calculation

### Key Data Structures

**EnrichResult DataFrame columns:**
- `Annot`: Annotation ID (GO:0006281, hsa03440)
- `Term`: Human-readable name
- `Pvalue`, `Padj`: Raw and adjusted p-values
- `Count`: Number of query genes in term
- `GeneID`: Semicolon-separated gene list
- `RichFactor`: Enrichment ratio (observed/expected)
- `GeneRatio`, `BgRatio`: Query and background ratios

**GSEA-specific columns:** `ES`, `NES`, `LeadingEdge`

### Test Fixtures

`tests/conftest.py` provides reusable fixtures:
- `sample_genes`, `sample_gene_scores` - test gene lists
- `sample_go_data`, `sample_kegg_data`, `sample_geneset_db` - annotation DataFrames
- `sample_enrich_result` - pre-built EnrichResult object

## pathwaydb Integration

richAnn can directly use annotation data from the `pathwaydb` package via converter functions:

```python
import richAnn as ra
from pathwaydb import GOAnnotationDB, KEGGAnnotationDB

# GO annotations
go_db = GOAnnotationDB('go_human.db')
go_data = ra.from_pathwaydb_go(go_db, ontology="BP")
result = ra.richGO(genes, go_data, ontology="BP")

# KEGG annotations
kegg_db = KEGGAnnotationDB('kegg_human.db')
kegg_data = ra.from_pathwaydb_kegg(kegg_db)
result = ra.richKEGG(genes, kegg_data)
```

**Column mapping:**
- GO: `TERM` → `GOterm`, `Aspect` (P/F/C) → `Ontology` (BP/MF/CC), `term_name` → `GOname`
- KEGG: `PATH` → `Pathway`, `Annot` → `PathwayName`

## Key Implementation Details

- Gene matching is case-insensitive by default (`case_sensitive=False`)
- Multiple testing correction uses statsmodels `multipletests()` with automatic R-to-Python method name conversion
- Visualization functions handle edge cases where all p-values are identical
- GSEA uses weighted Kolmogorov-Smirnov running sum statistic (weight=1.0 default matches GSEApy)
