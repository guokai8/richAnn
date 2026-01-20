"""
Test updated enrichment with zscore and KEGG hierarchy columns
"""

import pandas as pd
import numpy as np
import richAnn as ra

print("Testing updated enrichment calculations...")
print("=" * 80)

# Create KEGG annotation WITH hierarchy columns (like your expected output)
kegg_data = []

# Glycolysis pathway
for gene in ['AKR1A1', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'GALM', 'ADH7', 'LDHAL6A'] + [f'BG{i}' for i in range(1, 58)]:
    kegg_data.append({
        'GeneID': gene,
        'Pathway': '00010',
        'PathwayName': 'Glycolysis / Gluconeogenesis',
        'Level3': 'Glycolysis / Gluconeogenesis',
        'Level2': 'Carbohydrate metabolism',
        'Level1': 'Metabolism'
    })

# Pyruvate metabolism
for gene in ['AKR1A1', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7', 'LDHAL6A'] + [f'BG{i}' for i in range(58, 96)]:
    kegg_data.append({
        'GeneID': gene,
        'Pathway': '00620',
        'PathwayName': 'Pyruvate metabolism',
        'Level3': 'Pyruvate metabolism',
        'Level2': 'Carbohydrate metabolism',
        'Level1': 'Metabolism'
    })

# Tyrosine metabolism
for gene in ['ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7'] + [f'BG{i}' for i in range(96, 126)]:
    kegg_data.append({
        'GeneID': gene,
        'Pathway': '00350',
        'PathwayName': 'Tyrosine metabolism',
        'Level3': 'Tyrosine metabolism',
        'Level2': 'Amino acid metabolism',
        'Level1': 'Metabolism'
    })

# Fatty acid degradation
for gene in ['ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7'] + [f'BG{i}' for i in range(126, 162)]:
    kegg_data.append({
        'GeneID': gene,
        'Pathway': '00071',
        'PathwayName': 'Fatty acid degradation',
        'Level3': 'Fatty acid degradation',
        'Level2': 'Lipid metabolism',
        'Level1': 'Metabolism'
    })

# Retinol metabolism
for gene in ['ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7'] + [f'BG{i}' for i in range(162, 223)]:
    kegg_data.append({
        'GeneID': gene,
        'Pathway': '00830',
        'PathwayName': 'Retinol metabolism',
        'Level3': 'Retinol metabolism',
        'Level2': 'Metabolism of cofactors and vitamins',
        'Level1': 'Metabolism'
    })

kegg_annotation = pd.DataFrame(kegg_data)

print(f"Created KEGG annotation: {len(kegg_annotation)} gene-pathway pairs")
print(f"Columns: {list(kegg_annotation.columns)}")
print()

# Test genes - the ADH family
test_genes = ['AKR1A1', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'GALM', 'ADH7', 'LDHAL6A']

print(f"Input genes: {test_genes}")
print()

# Run enrichment
result = ra.richKEGG(
    genes=test_genes,
    kodata=kegg_annotation,
    pvalue=0.05,
    padj_method='fdr_bh',
    minSize=5
)

print(f"Found {len(result.result)} enriched pathways")
print()

# Display results
if len(result.result) > 0:
    print("Enrichment Results:")
    print("=" * 80)

    # Check what columns we have
    print(f"Available columns: {list(result.result.columns)}")
    print()

    # Display key columns similar to expected output
    display_cols = ['Annot', 'Term', 'Annotated', 'Significant', 'RichFactor',
                    'FoldEnrichment', 'zscore', 'Pvalue', 'Padj']

    # Only show columns that exist
    available_display_cols = [col for col in display_cols if col in result.result.columns]

    print(result.result[available_display_cols].to_string(index=False))
    print()

    # Check for hierarchy columns
    if 'Level1' in result.result.columns:
        print("\nHierarchy columns present:")
        print(result.result[['Annot', 'Level3', 'Level2', 'Level1']].to_string(index=False))
    else:
        print("\nWARNING: Hierarchy columns (Level1, Level2, Level3) not found")

    print()

    # Check for 'ko' column
    if 'ko' in result.result.columns:
        print("✓ 'ko' column present")
    else:
        print("✗ 'ko' column missing")

    # Check for GeneID column
    if 'GeneID' in result.result.columns:
        print("✓ 'GeneID' column present")
        print(f"  Example: {result.result['GeneID'].iloc[0][:50]}...")
    else:
        print("✗ 'GeneID' column missing")

print()
print("=" * 80)
print("Test complete!")
