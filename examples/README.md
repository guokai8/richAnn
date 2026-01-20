# richAnn Examples

This directory contains example scripts demonstrating how to use richAnn for gene set enrichment analysis.

## Quick Start

### Realistic Example (Recommended)
**File:** `realistic_example.py`

A complete working example that demonstrates:
- GO term enrichment analysis
- KEGG pathway enrichment
- Visualization (bar plots, dot plots)
- EnrichResult methods (filtering, exporting)
- Term clustering

```bash
python3 realistic_example.py
```

This example uses realistic gene sets (apoptosis genes, cancer pathway genes) that show actual enrichment.

### Comprehensive Test Example
**File:** `test_major_functions.py`

Tests all major richAnn functions including:
- richGO, richKEGG, richGSEA
- All visualization functions (ggbar, ggdot, ggnetplot, ggnetwork, ggGSEA)
- Cross-sample comparison (compareResult, ggheatmap, comparedot)
- Term clustering (richCluster)

**Note:** Some functions in this comprehensive test may not find significant results due to the synthetic nature of the test data. Use `realistic_example.py` for working examples.

### Quick Start Example
**File:** `quickstart_example.py`

A simplified example focusing on the core workflow:
1. Prepare gene lists
2. Create/load annotation data
3. Run enrichment analysis
4. Create visualizations

## Important Notes

### Parameter Requirements

**padj_method**: Must be lowercase (e.g., 'fdr_bh' not 'BH')
```python
result = ra.richGO(genes, godata, padj_method='fdr_bh')  # Correct
result = ra.richGO(genes, godata, padj_method='BH')      # Will fail
```

**GSEA gene_scores**: Must be a dict, not pandas Series
```python
gene_scores_dict = gene_scores.to_dict()  # Convert Series to dict
result = ra.richGSEA(gene_scores_dict, geneset_db)
```

### Annotation Data Formats

**GO annotation** requires columns:
- `GeneID`: Gene identifier
- `GOterm`: GO term ID (e.g., 'GO:0006915')
- `GOname`: GO term name (e.g., 'apoptotic process')
- `Ontology`: 'BP', 'MF', or 'CC'

**KEGG annotation** requires columns:
- `GeneID`: Gene identifier
- `Pathway`: Pathway ID (e.g., 'hsa04210')
- `PathwayName`: Pathway name (e.g., 'Apoptosis')

**GSEA annotation** requires columns:
- `GeneSet`: Gene set ID
- `GeneSetName`: Gene set name
- `GeneID`: Gene identifier

## Loading Real Annotation Data

richAnn provides utilities to load standard annotation formats:

```python
import richAnn as ra

# Load GMT file (MSigDB format)
geneset_db = ra.load_gmt('c2.cp.kegg.v7.4.symbols.gmt')

# Load GO annotation from GAF file
go_annotation = ra.load_go_gaf('gene_association.goa_human')

# Load KEGG mapping
kegg_annotation = ra.load_kegg_mapping('hsa')
```

## Output Files

Examples create the following files:
- `example_go_enrichment.csv`: GO enrichment results
- `go_results.csv`: Exported enrichment results
- Plot files (if you uncomment the `fig.savefig()` lines)

## Troubleshooting

**No significant results found:**
- Relax the p-value cutoff: `pvalue=0.1` or `padj=0.2`
- Check that your gene list overlaps with the annotation
- Ensure gene IDs match the annotation format (case-sensitive by default)

**Method not recognized error:**
- Use lowercase for `padj_method`: 'fdr_bh', 'bonferroni', 'holm', etc.
- See statsmodels documentation for available methods

**Empty gene list after normalization:**
- Check gene ID case matching
- Set `case_sensitive=False` for case-insensitive matching
- Verify gene IDs exist in the annotation

## Next Steps

1. Run the realistic example to see how enrichment works
2. Replace the sample gene lists with your own data
3. Load real annotation data for your organism
4. Customize visualizations by adjusting plot parameters
5. Export results and save plots for your publications

For more information, see the main README.md and package documentation.
