# richAnn

**Rich Annotation**: Functional Enrichment Analysis and Visualization for Python

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Overview

richAnn is a comprehensive Python package for gene set enrichment analysis, providing powerful statistical methods and publication-quality visualizations. Inspired by the R package richR, richAnn brings advanced enrichment analysis capabilities to the Python ecosystem.

### Key Features

- **Multiple Enrichment Methods**
  - Gene Ontology (GO) enrichment analysis
  - KEGG pathway analysis
  - Gene Set Enrichment Analysis (GSEA) with permutation testing

- **Advanced Analysis**
  - Kappa-based term clustering to group similar results
  - Cross-sample comparison for differential enrichment
  - Multiple testing correction (Benjamini-Hochberg, Bonferroni, etc.)

- **Publication-Quality Visualizations**
  - Bar plots and dot plots for enrichment results
  - Network graphs (term-gene, term-term, multi-database networks)
  - Heatmaps for cross-sample comparison
  - GSEA enrichment score plots

- **Flexible Data Loading**
  - GMT file support for custom gene sets
  - GAF file parser for GO annotations
  - KEGG mapping file support
  - Custom annotation builder

## Installation

### From Source

```bash
git clone https://github.com/guokai8/richAnn.git
cd richAnn
pip install -e .
```

### Development Installation

To install with development dependencies (pytest, black, flake8):

```bash
pip install -e ".[dev]"
```

### Dependencies

- Python ≥ 3.8
- numpy ≥ 1.21.0
- pandas ≥ 1.3.0
- scipy ≥ 1.7.0
- matplotlib ≥ 3.4.0
- seaborn ≥ 0.11.0
- networkx ≥ 2.6.0
- statsmodels ≥ 0.13.0

## Quick Start

### Basic GO Enrichment Analysis

```python
import richAnn as ra
import pandas as pd

# Your gene list of interest
genes = ["BRCA1", "BRCA2", "TP53", "ATM", "CHEK2", "RAD51", "PALB2"]

# Load GO annotation (or prepare your own DataFrame)
go_data = pd.DataFrame({
    'GeneID': ['BRCA1', 'BRCA2', 'TP53', 'ATM', ...],
    'GOterm': ['GO:0006281', 'GO:0006281', ...],
    'GOname': ['DNA repair', 'DNA repair', ...],
    'Ontology': ['BP', 'BP', ...]
})

# Perform GO enrichment
result = ra.richGO(
    genes=genes,
    godata=go_data,
    ontology="BP",  # Biological Process
    pvalue=0.05,
    padj=0.05,
    minSize=5,
    maxSize=500
)

# View results
print(result)
result.summary()

# Export results
result.to_csv("go_enrichment_results.csv")
```

### Visualize Results

```python
# Bar plot
fig = ra.ggbar(result, top=10, filename="barplot.png")

# Dot plot
fig = ra.ggdot(result, top=10, filename="dotplot.png")

# Network visualization
fig = ra.ggnetwork(result, top=20, filename="network.png")
```

### KEGG Pathway Enrichment

```python
# Load KEGG annotation
kegg_data = pd.DataFrame({
    'GeneID': ['BRCA1', 'BRCA2', ...],
    'Pathway': ['hsa03440', 'hsa03440', ...],
    'PathwayName': ['Homologous recombination', ...]
})

# Perform KEGG enrichment
kegg_result = ra.richKEGG(
    genes=genes,
    kodata=kegg_data,
    pvalue=0.05,
    padj=0.05
)
```

### Gene Set Enrichment Analysis (GSEA)

```python
# Gene scores (e.g., from differential expression analysis)
gene_scores = {
    "BRCA1": 2.5,
    "BRCA2": 2.3,
    "TP53": 2.1,
    "CDK1": -1.8,
    "CDK2": -1.5,
    # ... more genes
}

# Load gene set database
geneset_db = pd.DataFrame({
    'GeneSet': ['DNA_REPAIR', 'DNA_REPAIR', ...],
    'GeneSetName': ['DNA repair pathway', ...],
    'GeneID': ['BRCA1', 'BRCA2', ...]
})

# Perform GSEA
gsea_result = ra.richGSEA(
    gene_scores=gene_scores,
    geneset_db=geneset_db,
    nperm=1000,
    min_size=15,
    max_size=500
)

# Visualize specific pathway
fig = ra.ggGSEA(gsea_result, gene_scores, term_id="DNA_REPAIR")
```

## Advanced Features

### Term Clustering

Cluster similar enrichment terms based on gene overlap:

```python
# Cluster enrichment results
clustered = ra.richCluster(
    result,
    cutoff=0.5,      # Kappa score threshold
    minSize=5,       # Minimum cluster size
    escore=3         # Enrichment score cutoff (-log10(padj))
)

# Visualize clustered network
fig = ra.ggnetwork(clustered, top=30)
```

### Cross-Sample Comparison

Compare enrichment across multiple conditions:

```python
# Perform enrichment for multiple samples
sample1_result = ra.richGO(sample1_genes, go_data, ontology="BP")
sample2_result = ra.richGO(sample2_genes, go_data, ontology="BP")
sample3_result = ra.richGO(sample3_genes, go_data, ontology="BP")

# Compare results
comparison = ra.compareResult({
    'Control': sample1_result,
    'Treatment1': sample2_result,
    'Treatment2': sample3_result
})

# Visualize comparison
fig = ra.ggheatmap(comparison, top=30)
fig = ra.comparedot(comparison, top=10)
```

### Load Annotations from Standard Formats

```python
# Load from GMT file (e.g., MSigDB)
geneset_db = ra.load_gmt("c2.cp.kegg.v7.4.symbols.gmt")

# Load GO annotations from GAF file
go_data = ra.load_go_gaf(
    "goa_human.gaf",
    ontology="BP",
    evidence_codes=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"]
)

# Load KEGG pathway mapping
kegg_data = ra.load_kegg_mapping("kegg_human.txt", organism="hsa")

# Create custom annotation
custom_sets = {
    "DNA_REPAIR": ["BRCA1", "BRCA2", "TP53", "ATM"],
    "CELL_CYCLE": ["CDK1", "CDK2", "CCNA1", "CCNB1"]
}
custom_annot = ra.create_custom_annotation(custom_sets)
```

### Validate Annotation Quality

```python
# Validate annotation format and get statistics
stats = ra.validate_annotation_format(
    go_data,
    required_cols=['GeneID', 'GOterm', 'GOname', 'Ontology']
)

print(f"Valid: {stats['valid']}")
print(f"Unique genes: {stats['n_genes']}")
print(f"Unique annotations: {stats['n_annotations']}")
print(f"Warnings: {stats['warnings']}")
```

## Integration with pathwaydb

richAnn integrates seamlessly with [pathwaydb](https://github.com/guokai8/pathwaydb), a lightweight Python library for downloading and querying biological pathway annotations. This is the **recommended approach** for production use.

### Why Use pathwaydb?

- **Easy setup**: Download annotations with one command
- **Offline access**: Query locally without internet
- **Always up-to-date**: Fresh data from official sources (GO, KEGG)
- **12 model species**: Human, mouse, rat, zebrafish, fly, worm, yeast, and more
- **Fast queries**: Millisecond lookups on SQLite databases

### Installation

```bash
# Install pathwaydb
pip install pathwaydb

# Or install from source
git clone https://github.com/guokai8/pathwaydb.git
cd pathwaydb && pip install -e .
```

### Quick Start with pathwaydb

#### GO Enrichment Analysis

```python
import richAnn as ra
from pathwaydb import GO

# Step 1: Download GO annotations (only needed once)
go = GO(storage_path='go_human.db')
go.download_annotations(species='human')
# Output: Downloaded 500,000+ annotations with term names

# Step 2: Convert to richAnn format
go_data = ra.from_pathwaydb_go(go.storage, ontology="BP")
print(f"Loaded {go_data['GOterm'].nunique()} GO terms for {go_data['GeneID'].nunique()} genes")

# Step 3: Run enrichment analysis
genes = ["BRCA1", "BRCA2", "TP53", "ATM", "CHEK2", "RAD51", "PALB2", "NBN"]
result = ra.richGO(genes, go_data, ontology="BP", pvalue=0.05, padj=0.1)

# Step 4: View and visualize results
print(result)
result.summary()
fig = ra.ggdot(result, top=15)
```

#### KEGG Pathway Enrichment

```python
import richAnn as ra
from pathwaydb import KEGG

# Step 1: Download KEGG annotations (includes pathway hierarchy!)
kegg = KEGG(species='hsa', storage_path='kegg_human.db')
kegg.download_annotations()  # Automatically downloads pathway hierarchy
kegg.convert_ids_to_symbols()  # Convert Entrez IDs to gene symbols

# Step 2: Convert to richAnn format
kegg_data = ra.from_pathwaydb_kegg(kegg.storage)
print(f"Loaded {kegg_data['Pathway'].nunique()} pathways")

# Step 3: Run enrichment analysis
result = ra.richKEGG(genes, kegg_data, pvalue=0.05)

# Step 4: Results now include pathway hierarchy!
print(result.result[['Annot', 'Term', 'Level1', 'Level2', 'Padj']].head())
#         Annot                    Term           Level1                  Level2      Padj
# 0   hsa03440  Homologous recombination  Genetic Info...  Replication and repair  0.0001
# 1   hsa04110              Cell cycle    Cellular Processes        Cell growth...  0.0005

# Step 5: Visualize
fig = ra.ggbar(result, top=10)
fig = ra.ggnetwork(result, top=20)
```

**KEGG Pathway Hierarchy Levels:**
- **Level1**: Top-level category (e.g., "Metabolism", "Human Diseases", "Cellular Processes")
- **Level2**: Sub-category (e.g., "Carbohydrate metabolism", "Cancer", "Cell growth and death")
- **Level3**: Pathway name (same as PathwayName)

### Complete Workflow Example

```python
import richAnn as ra
from pathwaydb import GO, KEGG

# =============================================================================
# Setup: Download annotations (run once, data is cached locally)
# =============================================================================

# GO annotations
go = GO(storage_path='go_human.db')
go.download_annotations(species='human')

# KEGG annotations
kegg = KEGG(species='hsa', storage_path='kegg_human.db')
kegg.download_annotations()
kegg.convert_ids_to_symbols()

# =============================================================================
# Analysis: Your gene list (e.g., from differential expression)
# =============================================================================

# Example: DNA damage response genes
my_genes = [
    "TP53", "BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2",
    "RAD51", "PALB2", "NBN", "MRE11", "RAD50", "XRCC1", "PARP1"
]

# =============================================================================
# GO Enrichment (Biological Process)
# =============================================================================

go_data = ra.from_pathwaydb_go(go.storage, ontology="BP")
go_result = ra.richGO(
    genes=my_genes,
    godata=go_data,
    ontology="BP",
    pvalue=0.05,
    padj=0.1,
    minSize=5,
    maxSize=500
)

print("=== GO Enrichment Results ===")
go_result.summary()

# Top enriched terms
print("\nTop 10 GO Terms:")
for _, row in go_result.top(10).result.iterrows():
    print(f"  {row['Annot']}: {row['Term']} (Padj={row['Padj']:.2e})")

# =============================================================================
# KEGG Enrichment
# =============================================================================

kegg_data = ra.from_pathwaydb_kegg(kegg.storage)
kegg_result = ra.richKEGG(
    genes=my_genes,
    kodata=kegg_data,
    pvalue=0.05,
    padj=0.1
)

print("\n=== KEGG Enrichment Results ===")
kegg_result.summary()

# =============================================================================
# Visualizations
# =============================================================================

# GO dot plot
fig_go = ra.ggdot(go_result, top=15, usePadj=True)
fig_go.savefig("go_dotplot.png", dpi=300, bbox_inches='tight')

# KEGG bar plot
fig_kegg = ra.ggbar(kegg_result, top=10, colorBy='Padj')
fig_kegg.savefig("kegg_barplot.png", dpi=300, bbox_inches='tight')

# Combined network
fig_net = ra.ggnetmap(
    {'GO': go_result, 'KEGG': kegg_result},
    top=10
)
fig_net.savefig("combined_network.png", dpi=300, bbox_inches='tight')

# =============================================================================
# Export Results
# =============================================================================

go_result.to_csv("go_enrichment_results.csv")
kegg_result.to_csv("kegg_enrichment_results.csv")
```

### Multi-Species Analysis

pathwaydb supports 12 model organisms. Use the same workflow for any species:

```python
from pathwaydb import GO, get_supported_species

# List all supported species
print(get_supported_species())
# ['arabidopsis', 'chicken', 'cow', 'dog', 'fly', 'human', 'mouse', 'pig', 'rat', 'worm', 'yeast', 'zebrafish']

# Mouse analysis
go_mouse = GO(storage_path='go_mouse.db')
go_mouse.download_annotations(species='mouse')
go_data_mouse = ra.from_pathwaydb_go(go_mouse.storage, ontology="BP")

# Zebrafish analysis
go_zebrafish = GO(storage_path='go_zebrafish.db')
go_zebrafish.download_annotations(species='zebrafish')
go_data_zebrafish = ra.from_pathwaydb_go(go_zebrafish.storage, ontology="BP")

# Fly (Drosophila) analysis
go_fly = GO(storage_path='go_fly.db')
go_fly.download_annotations(species='fly')
go_data_fly = ra.from_pathwaydb_go(go_fly.storage, ontology="BP")
```

### Using Cached Annotations

For repeated analyses, use pathwaydb's centralized cache:

```python
from pathwaydb import GO

# Load from cache (auto-downloads if not cached)
go = GO.from_cache(species='human')

# Convert and analyze
go_data = ra.from_pathwaydb_go(go.storage, ontology="BP")
result = ra.richGO(genes, go_data, ontology="BP")
```

### Filter GO by Evidence Codes

```python
# Download with all evidence codes
go = GO(storage_path='go_human.db')
go.download_annotations(species='human')

# Filter for experimental evidence only during conversion
go_data = ra.from_pathwaydb_go(go.storage, ontology="BP")

# Or filter using pathwaydb before conversion
experimental = go.storage.filter(evidence_codes=['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'])
# Then manually build DataFrame if needed
```

### pathwaydb Data Format Reference

**GO columns from pathwaydb:**
| pathwaydb | richAnn | Description |
|-----------|---------|-------------|
| `GeneID` | `GeneID` | Gene symbol |
| `TERM` | `GOterm` | GO term ID (GO:0006281) |
| `Aspect` | `Ontology` | P→BP, F→MF, C→CC |
| `term_name` | `GOname` | Term description |

**KEGG columns from pathwaydb:**
| pathwaydb | richAnn | Description |
|-----------|---------|-------------|
| `GeneID` | `GeneID` | Gene symbol |
| `PATH` | `Pathway` | Pathway ID (hsa04110) |
| `Annot` | `PathwayName` | Pathway name |
| `Level1` | `Level1` | Top-level category (e.g., "Metabolism") |
| `Level2` | `Level2` | Sub-category (e.g., "Carbohydrate metabolism") |
| `Level3` | `Level3` | Pathway name (same as PathwayName) |

## Working with EnrichResult Objects

EnrichResult is the core data structure returned by all enrichment functions:

```python
# Filter results
filtered = result.filter(
    pvalue=0.01,           # P-value cutoff
    padj=0.05,             # Adjusted p-value cutoff
    min_count=5,           # Minimum gene count
    min_richfactor=1.5     # Minimum enrichment factor
)

# Get top N results
top10 = result.top(n=10, orderby='Padj')

# Get gene-term details
detail_df = result.detail()  # Returns DataFrame with one row per gene-term pair

# Print summary statistics
result.summary()

# Export to various formats
result.to_csv("results.csv")
result.to_excel("results.xlsx", sheet_name="GO_BP")

# Access underlying DataFrame
df = result.result
print(df.columns)  # View available columns
```

## Visualization Gallery

### Bar Plot
```python
ra.ggbar(result, top=15, colorBy='Padj', horiz=True, filename="barplot.pdf", dpi=300)
```

### Dot Plot
```python
ra.ggdot(result, top=15, usePadj=True, size_range=(50, 500), filename="dotplot.pdf")
```

### Term-Gene Network
```python
ra.ggnetplot(result, top=20, layout='spring', filename="term_gene_network.pdf")
```

### Term-Term Similarity Network
```python
ra.ggnetwork(result, top=30, weightcut=0.2, layout='spring', filename="term_network.pdf")
```

### Multi-Database Network
```python
ra.ggnetmap(
    {'GO': go_result, 'KEGG': kegg_result},
    top=15,
    filename="combined_network.pdf"
)
```

### Comparison Heatmap
```python
ra.ggheatmap(
    comparison_df,
    top=30,
    cluster_rows=True,
    cluster_cols=False,
    filename="heatmap.pdf"
)
```

## Understanding the Results

### Key Columns in Results

- **Annot**: Annotation ID (e.g., GO:0006281, hsa03440)
- **Term**: Human-readable term name
- **Pvalue**: Raw p-value from hypergeometric test
- **Padj**: Adjusted p-value (Benjamini-Hochberg FDR by default)
- **Count**: Number of query genes in this term
- **GeneID**: Semicolon-separated list of genes
- **RichFactor**: Enrichment ratio (observed/expected)
- **GeneRatio**: Query genes in term / total query genes (k/n)
- **BgRatio**: Background genes in term / total background (M/N)
- **OddsRatio**: Odds ratio from Fisher's exact test

### Statistical Methods

#### Hypergeometric Test (richGO, richKEGG)
Tests for over-representation using the hypergeometric distribution:
- P(X ≥ k) where k = number of overlapping genes
- Accounts for gene set size and background universe

#### GSEA (richGSEA)
Ranks genes by score and calculates enrichment score:
- Weighted running sum based on gene ranks
- Permutation testing for significance
- Normalized enrichment score (NES) accounts for gene set size

## Tips and Best Practices

### 1. Choose Appropriate Background
```python
# Option 1: Use all genes in annotation (default)
result = ra.richGO(genes, go_data, ontology="BP")

# Option 2: Provide custom background (recommended for RNA-seq)
background = list(expressed_genes)  # All genes detected in experiment
result = ra.richGO(genes, go_data, ontology="BP", universe=background)
```

### 2. Handle Gene ID Formatting
```python
# Case-insensitive matching (default)
result = ra.richGO(genes, go_data, case_sensitive=False)

# Case-sensitive matching
result = ra.richGO(genes, go_data, case_sensitive=True)
```

### 3. Adjust Gene Set Size Limits
```python
# Exclude very small and very large gene sets
result = ra.richGO(
    genes, go_data,
    minSize=10,   # Exclude terms with < 10 genes
    maxSize=300   # Exclude terms with > 300 genes
)
```

### 4. Multiple Testing Correction
```python
# Benjamini-Hochberg FDR (default)
result = ra.richGO(genes, go_data, padj_method="BH")

# Bonferroni correction (more conservative)
result = ra.richGO(genes, go_data, padj_method="bonferroni")

# No correction (use with caution)
result = ra.richGO(genes, go_data, padj_method="none")
```

## Common Issues and Solutions

### Issue: "No query genes found in background universe"
**Solution**: Check gene ID format matches annotation. Try case_sensitive=False.

### Issue: "No enriched terms found"
**Solutions**:
- Relax p-value cutoff: `pvalue=0.1`
- Adjust gene set size limits: `minSize=3, maxSize=1000`
- Check that your gene list is biologically relevant
- Verify annotation database matches your organism

### Issue: GSEA is slow
**Solution**: Reduce permutations (but keep ≥ 1000): `nperm=1000`

### Issue: Visualizations have overlapping labels
**Solutions**:
- Reduce `top` parameter: `top=10`
- Increase figure size: `figsize=(14, 10)`
- Decrease font size: `fontsize=8`

## Testing

Run the test suite:

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=richAnn

# Run specific test file
pytest tests/test_enrichment.py

# Run with verbose output
pytest -v
```

## Development

### Code Formatting
```bash
black richAnn/
```

### Linting
```bash
flake8 richAnn/
```

## Citation

If you use richAnn in your research, please cite:

```
[Citation information will be added upon publication]
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Inspired by the R package [richR](https://github.com/hurlab/richR)
- Built on the excellent scientific Python ecosystem (NumPy, pandas, SciPy, matplotlib)

## Contact

For questions and feedback:
- Open an issue on [GitHub](https://github.com/yourusername/richAnn/issues)
- Email: support@richann.org

## Changelog

### Version 1.0.0 (Current)

#### New Features
- GO, KEGG, and GSEA enrichment analysis
- Kappa-based term clustering
- 8 visualization types
- GMT, GAF, and KEGG file loaders
- Custom annotation builder
- Cross-sample comparison tools
- Comprehensive test suite

#### Bug Fixes
- Fixed module import error
- Improved GSEA performance (28x faster)
- Added edge case handling for visualizations
- Enhanced error messages with actionable suggestions

#### Improvements
- Added input validation for all functions
- Comprehensive documentation with examples
- Test coverage for core functionality
- Better error handling throughout

---

**Happy Enriching! 🧬**
