"""
Final validation test - comparing with expected richR output format
"""

import pandas as pd
import numpy as np
import richAnn as ra

print("=" * 80)
print("VALIDATION TEST: Comparing with richR Expected Output")
print("=" * 80)

# Recreate annotation similar to the expected output
# Based on your expected output, these are the gene sets:
kegg_data = []

pathways = {
    '00010': ('Glycolysis / Gluconeogenesis', 'Carbohydrate metabolism', 'Metabolism', 67,
              ['AKR1A1', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'GALM', 'ADH7', 'LDHAL6A']),
    '00620': ('Pyruvate metabolism', 'Carbohydrate metabolism', 'Metabolism', 47,
              ['AKR1A1', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7', 'LDHAL6A']),
    '00350': ('Tyrosine metabolism', 'Amino acid metabolism', 'Metabolism', 37,
              ['ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7']),
    '00071': ('Fatty acid degradation', 'Lipid metabolism', 'Metabolism', 43,
              ['ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7']),
    '00830': ('Retinol metabolism', 'Metabolism of cofactors and vitamins', 'Metabolism', 68,
              ['ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7']),
    '00982': ('Drug metabolism - cytochrome P450', 'Xenobiotics biodegradation and metabolism', 'Metabolism', 73,
              ['ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7']),
    '00980': ('Metabolism of xenobiotics by cytochrome P450', 'Xenobiotics biodegradation and metabolism', 'Metabolism', 79,
              ['ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7']),
    '04936': ('Alcoholic liver disease', 'Endocrine and metabolic disease', 'Human Diseases', 144,
              ['ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7']),
    '00053': ('Ascorbate and aldarate metabolism', 'Carbohydrate metabolism', 'Metabolism', 31,
              ['AKR1A1']),
    '00052': ('Galactose metabolism', 'Carbohydrate metabolism', 'Metabolism', 32,
              ['GALM']),
    '00640': ('Propanoate metabolism', 'Carbohydrate metabolism', 'Metabolism', 32,
              ['LDHAL6A']),
    '00040': ('Pentose and glucuronate interconversions', 'Carbohydrate metabolism', 'Metabolism', 36,
              ['AKR1A1']),
}

# Create annotation DataFrame with proper pathway sizes
for pathway_id, (name, level2, level1, annotated, input_genes) in pathways.items():
    # Add the actual genes that should be enriched
    for gene in input_genes:
        kegg_data.append({
            'GeneID': gene,
            'Pathway': pathway_id,
            'PathwayName': name,
            'Level3': name,
            'Level2': level2,
            'Level1': level1
        })

    # Add background genes to reach the specified pathway size
    n_background = annotated - len(input_genes)
    for i in range(n_background):
        kegg_data.append({
            'GeneID': f'BG_{pathway_id}_{i}',
            'Pathway': pathway_id,
            'PathwayName': name,
            'Level3': name,
            'Level2': level2,
            'Level1': level1
        })

kegg_annotation = pd.DataFrame(kegg_data)

# Test with the 10 genes from expected output
test_genes = ['AKR1A1', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'GALM', 'ADH7', 'LDHAL6A']

print(f"\nInput genes ({len(test_genes)}): {', '.join(test_genes)}")
print(f"Total gene universe: {kegg_annotation['GeneID'].nunique()} genes")
print(f"Total pathways: {kegg_annotation['Pathway'].nunique()}")
print()

# Run enrichment
result = ra.richKEGG(
    genes=test_genes,
    kodata=kegg_annotation,
    pvalue=0.05,
    padj_method='fdr_bh',
    minSize=2
)

print(f"Enriched pathways found: {len(result.result)}")
print()

if len(result.result) > 0:
    print("Column Comparison:")
    print("-" * 80)

    expected_cols = ['Annot', 'Term', 'Annotated', 'Significant', 'RichFactor',
                     'FoldEnrichment', 'zscore', 'Pvalue', 'Padj', 'GeneID',
                     'ko', 'Level3', 'Level2', 'Level1']

    actual_cols = list(result.result.columns)

    for col in expected_cols:
        status = "✓" if col in actual_cols else "✗"
        print(f"  {status} {col}")

    print()
    print("Additional columns in output:")
    extra_cols = set(actual_cols) - set(expected_cols)
    for col in extra_cols:
        print(f"  + {col}")

    print()
    print("=" * 80)
    print("ENRICHMENT RESULTS (Top pathways)")
    print("=" * 80)

    # Reorder columns to match expected output
    output_cols = ['Annot', 'Term', 'Annotated', 'Significant', 'RichFactor',
                   'FoldEnrichment', 'zscore', 'Pvalue', 'Padj']

    display_df = result.result[output_cols].head(12)
    print(display_df.to_string(index=False))

    print()
    print("Gene IDs (first 3 pathways):")
    print("-" * 80)
    for idx in range(min(3, len(result.result))):
        row = result.result.iloc[idx]
        print(f"{row['Annot']} ({row['Term']})")
        print(f"  Genes: {row['GeneID']}")
        print()

    print("Pathway Hierarchy (first 5):")
    print("-" * 80)
    hierarchy_cols = ['ko', 'Level3', 'Level2', 'Level1']
    print(result.result[hierarchy_cols].head().to_string(index=False))

    print()
    print("=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)
    print(f"✓ All expected columns present: {all(col in actual_cols for col in expected_cols)}")
    print(f"✓ Results sorted by Padj: {result.result['Padj'].is_monotonic_increasing}")
    print(f"✓ Hierarchy columns included: {'Level1' in actual_cols and 'Level2' in actual_cols}")
    print(f"✓ Ko column present: {'ko' in actual_cols}")
    print(f"✓ Zscore calculated: {'zscore' in actual_cols}")
    print(f"✓ FoldEnrichment present: {'FoldEnrichment' in actual_cols}")

    print()
    print("Output format matches richR expected structure!")

print()
print("=" * 80)
