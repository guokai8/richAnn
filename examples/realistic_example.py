"""
Realistic Working Example for richAnn Package

This example uses realistic gene sets that will show actual enrichment.
"""

import pandas as pd
import numpy as np
import richAnn as ra

print("=" * 80)
print("richAnn - Realistic Enrichment Example")
print("=" * 80)

# ============================================================================
# Setup: Create realistic annotation with apoptosis genes
# ============================================================================
print("\n[Setup] Creating gene annotation database...")

# Define gene universe (background)
all_genes = [
    # Apoptosis pathway genes
    'TP53', 'BCL2', 'BAX', 'CASP3', 'CASP9', 'APAF1', 'CYCS', 'BID', 'BAK1',
    'CASP8', 'FADD', 'TRADD', 'BCL2L1', 'MCL1', 'BCL2L11', 'BAD', 'PMAIP1',
    # Cell cycle genes
    'CDK4', 'CDK6', 'CCND1', 'E2F1', 'CDKN1A', 'CDKN1B', 'RB1', 'CDKN2A',
    # MAPK pathway
    'MAPK1', 'MAPK3', 'RAF1', 'MAP2K1', 'BRAF', 'KRAS', 'NRAS', 'EGFR',
    # PI3K-AKT pathway
    'PIK3CA', 'AKT1', 'PTEN', 'MTOR', 'TSC1', 'TSC2', 'FOXO1',
    # DNA repair
    'BRCA1', 'ATM', 'CHEK2', 'RAD51', 'XRCC1',
    # Inflammation
    'IL6', 'TNF', 'IFNG', 'IL1B', 'NFKB1',
    # Angiogenesis
    'VEGFA', 'FGF2', 'PDGFA', 'TGFB1',
    # Random background genes
] + [f'BACKGROUND{i}' for i in range(1, 101)]

print(f"   Total gene universe: {len(all_genes)} genes")

# Create GO annotation - realistic gene-term assignments
go_annotation = pd.DataFrame([
    # Apoptosis genes
    *[{'GeneID': g, 'GOterm': 'GO:0006915', 'GOname': 'apoptotic process', 'Ontology': 'BP'}
      for g in ['TP53', 'BCL2', 'BAX', 'CASP3', 'CASP9', 'APAF1', 'CYCS', 'BID', 'BAK1', 'CASP8', 'FADD', 'BCL2L1']],
    # Cell proliferation
    *[{'GeneID': g, 'GOterm': 'GO:0008283', 'GOname': 'cell proliferation', 'Ontology': 'BP'}
      for g in ['CDK4', 'CDK6', 'CCND1', 'E2F1', 'EGFR', 'MAPK1', 'AKT1', 'MTOR']],
    # DNA repair
    *[{'GeneID': g, 'GOterm': 'GO:0006281', 'GOname': 'DNA repair', 'Ontology': 'BP'}
      for g in ['TP53', 'BRCA1', 'ATM', 'CHEK2', 'RAD51', 'XRCC1']],
    # Signal transduction
    *[{'GeneID': g, 'GOterm': 'GO:0007165', 'GOname': 'signal transduction', 'Ontology': 'BP'}
      for g in ['EGFR', 'KRAS', 'BRAF', 'MAPK1', 'MAPK3', 'PIK3CA', 'AKT1', 'PTEN', 'MTOR']],
    # Inflammatory response
    *[{'GeneID': g, 'GOterm': 'GO:0006954', 'GOname': 'inflammatory response', 'Ontology': 'BP'}
      for g in ['IL6', 'TNF', 'IFNG', 'IL1B', 'NFKB1']],
    # Add background genes to various terms
    *[{'GeneID': f'BACKGROUND{i}', 'GOterm': 'GO:0008283', 'GOname': 'cell proliferation', 'Ontology': 'BP'}
      for i in range(1, 21)],
    *[{'GeneID': f'BACKGROUND{i}', 'GOterm': 'GO:0007165', 'GOname': 'signal transduction', 'Ontology': 'BP'}
      for i in range(21, 46)],
])

# Create KEGG annotation
kegg_annotation = pd.DataFrame([
    # Apoptosis pathway
    *[{'GeneID': g, 'Pathway': 'hsa04210', 'PathwayName': 'Apoptosis'}
      for g in ['TP53', 'BCL2', 'BAX', 'CASP3', 'CASP9', 'BID', 'BAK1', 'CASP8', 'FADD']],
    # PI3K-Akt pathway
    *[{'GeneID': g, 'Pathway': 'hsa04151', 'PathwayName': 'PI3K-Akt signaling pathway'}
      for g in ['EGFR', 'PIK3CA', 'AKT1', 'PTEN', 'MTOR', 'TSC1', 'TSC2', 'FOXO1', 'VEGFA']],
    # Pathways in cancer
    *[{'GeneID': g, 'Pathway': 'hsa05200', 'PathwayName': 'Pathways in cancer'}
      for g in ['TP53', 'KRAS', 'BRAF', 'PIK3CA', 'EGFR', 'RB1', 'CDKN2A', 'BCL2', 'VEGFA']],
    # MAPK pathway
    *[{'GeneID': g, 'Pathway': 'hsa04010', 'PathwayName': 'MAPK signaling pathway'}
      for g in ['KRAS', 'BRAF', 'NRAS', 'EGFR', 'MAPK1', 'MAPK3', 'RAF1', 'MAP2K1']],
    # Add background genes
    *[{'GeneID': f'BACKGROUND{i}', 'Pathway': 'hsa04151', 'PathwayName': 'PI3K-Akt signaling pathway'}
      for i in range(1, 16)],
])

print(f"   GO annotation: {len(go_annotation)} gene-term pairs ({go_annotation['GOterm'].nunique()} unique terms)")
print(f"   KEGG annotation: {len(kegg_annotation)} gene-pathway pairs ({kegg_annotation['Pathway'].nunique()} unique pathways)")

# ============================================================================
# Test 1: GO Enrichment with Apoptosis Genes
# ============================================================================
print("\n[Test 1] GO Enrichment - Apoptosis Gene Set")
print("-" * 80)

# Select mostly apoptosis genes (should show enrichment for GO:0006915)
apoptosis_genes = ['TP53', 'BCL2', 'BAX', 'CASP3', 'CASP9', 'APAF1', 'CYCS', 'BID']

go_result = ra.richGO(
    genes=apoptosis_genes,
    godata=go_annotation,
    pvalue=0.05,
    padj_method='fdr_bh',
    minSize=3
)

print(f"Input: {len(apoptosis_genes)} genes")
print(f"Result: {len(go_result.result)} significant terms\n")

if len(go_result.result) > 0:
    print("Enriched GO terms:")
    display_cols = ['Term', 'Count', 'Pvalue', 'Padj', 'RichFactor', 'GeneRatio']
    print(go_result.result[display_cols].to_string(index=False))

    # Create visualization
    print("\nCreating bar plot...")
    fig = ra.ggbar(go_result, top=5, colorBy='Padj')
    print("✓ Visualization created (uncomment savefig to save)")
    # fig.savefig('example_go_barplot.png', dpi=300, bbox_inches='tight')
else:
    print("No significant enrichment found")

# ============================================================================
# Test 2: KEGG Enrichment with Cancer Pathway Genes
# ============================================================================
print("\n[Test 2] KEGG Enrichment - Cancer Pathway Genes")
print("-" * 80)

# Select genes involved in cancer pathways
cancer_genes = ['TP53', 'KRAS', 'BRAF', 'PIK3CA', 'EGFR', 'RB1', 'BCL2']

kegg_result = ra.richKEGG(
    genes=cancer_genes,
    kodata=kegg_annotation,
    pvalue=0.05,
    padj_method='fdr_bh',
    minSize=3
)

print(f"Input: {len(cancer_genes)} genes")
print(f"Result: {len(kegg_result.result)} significant pathways\n")

if len(kegg_result.result) > 0:
    print("Enriched KEGG pathways:")
    display_cols = ['Term', 'Count', 'Pvalue', 'Padj', 'RichFactor', 'GeneRatio']
    print(kegg_result.result[display_cols].to_string(index=False))

    # Create dot plot
    print("\nCreating dot plot...")
    fig = ra.ggdot(kegg_result, top=5)
    print("✓ Visualization created")
else:
    print("No significant enrichment found")

# ============================================================================
# Test 3: EnrichResult Methods
# ============================================================================
print("\n[Test 3] EnrichResult Object Methods")
print("-" * 80)

if len(go_result.result) > 0:
    # Filtering
    strict_filter = go_result.filter(pvalue=0.01)
    print(f"Original results: {len(go_result.result)} terms")
    print(f"After filter (p<0.01): {len(strict_filter.result)} terms")

    # Top results
    top_terms = go_result.top(n=2)
    print(f"Top 2 results: {len(top_terms.result)} terms")

    # Export to CSV
    go_result.to_csv('example_go_enrichment.csv')
    print("✓ Results exported to 'example_go_enrichment.csv'")

    # Access parameters
    print(f"\nAnalysis parameters:")
    for key, value in go_result.parameters.items():
        print(f"   {key}: {value}")

# ============================================================================
# Test 4: Term Clustering
# ============================================================================
print("\n[Test 4] Term Clustering")
print("-" * 80)

if len(go_result.result) >= 2:
    clustered = ra.richCluster(
        go_result,
        cutoff=0.3,
        minSize=2
    )

    if 'Cluster' in clustered.result.columns:
        n_clusters = clustered.result['Cluster'].nunique()
        print(f"Clustered {len(clustered.result)} terms into {n_clusters} clusters")
        print("\nCluster assignments:")
        print(clustered.result[['Term', 'Cluster']].to_string(index=False))
    else:
        print("No clusters formed (terms too dissimilar)")
else:
    print("Not enough terms for clustering")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 80)
print("EXAMPLE COMPLETED SUCCESSFULLY!")
print("=" * 80)
print(f"✓ GO enrichment: {len(go_result.result)} terms found")
print(f"✓ KEGG enrichment: {len(kegg_result.result)} pathways found")
print(f"✓ Visualizations created")
print(f"✓ Results exported to CSV")
print("\nKey findings:")
print(f"- Apoptosis genes showed enrichment for '{go_result.result.iloc[0]['Term'] if len(go_result.result) > 0 else 'N/A'}'")
print(f"- Cancer genes showed enrichment for '{kegg_result.result.iloc[0]['Term'] if len(kegg_result.result) > 0 else 'N/A'}'")
print("\nFiles created:")
print("- example_go_enrichment.csv")
print("\nNext steps:")
print("- Uncomment fig.savefig() lines to save plots")
print("- Try with your own gene lists")
print("- Load real annotation data using richAnn.load_gmt() or load_go_gaf()")
print("=" * 80)
