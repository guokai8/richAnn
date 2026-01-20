"""
Quick Start Example for richAnn Package

This demonstrates the core functionality with working examples.
"""

import pandas as pd
import numpy as np
import richAnn as ra

print("=" * 80)
print("richAnn Quick Start Example")
print("=" * 80)

# ============================================================================
# 1. PREPARE SAMPLE DATA
# ============================================================================
print("\n1. Setting up sample data...")

# Create a larger gene universe
np.random.seed(42)
all_genes = [
    'TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS', 'PTEN', 'AKT1', 'PIK3CA',
    'BRAF', 'NRAS', 'ERBB2', 'RB1', 'CDKN2A', 'ATM', 'CHEK2', 'MDM2',
    'BCL2', 'BAX', 'CASP3', 'CASP9', 'APAF1', 'CYCS', 'BID', 'BAK1',
    'VEGFA', 'FGF2', 'PDGFA', 'TGFB1', 'IL6', 'TNF', 'IFNG', 'IL1B',
    'STAT3', 'JAK2', 'SMAD4', 'NOTCH1', 'WNT1', 'CTNNB1', 'GSK3B', 'AXIN1',
    'CDK4', 'CDK6', 'CCND1', 'E2F1', 'CDKN1A', 'CDKN1B', 'SKP2', 'CUL1',
    'MAPK1', 'MAPK3', 'RAF1', 'MAP2K1', 'MAP2K2', 'SOS1', 'GRB2', 'SHC1',
    'MTOR', 'RICTOR', 'RAPTOR', 'TSC1', 'TSC2', 'RHEB', 'AKT2', 'AKT3',
] + [f'GENE{i}' for i in range(1, 200)]

# Differentially expressed genes (subset of universe)
de_genes = [
    'TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS', 'PTEN', 'AKT1', 'PIK3CA',
    'BRAF', 'NRAS', 'ERBB2', 'RB1', 'CDKN2A', 'ATM', 'CHEK2', 'MDM2',
    'BCL2', 'BAX', 'CASP3', 'CASP9', 'APAF1', 'CYCS', 'BID', 'BAK1',
    'VEGFA', 'FGF2', 'PDGFA', 'TGFB1', 'IL6', 'TNF', 'IFNG', 'IL1B'
]

print(f"   Universe: {len(all_genes)} genes")
print(f"   DE genes: {len(de_genes)} genes")

# ============================================================================
# 2. CREATE ANNOTATION DATA
# ============================================================================
print("\n2. Creating annotation data...")

# GO annotation
go_terms = {
    'GO:0006915': ('apoptotic process', 'BP'),
    'GO:0008283': ('cell proliferation', 'BP'),
    'GO:0006281': ('DNA repair', 'BP'),
    'GO:0007165': ('signal transduction', 'BP'),
    'GO:0006954': ('inflammatory response', 'BP'),
    'GO:0045765': ('regulation of angiogenesis', 'BP'),
}

go_mapping = {
    'GO:0006915': ['TP53', 'BCL2', 'BAX', 'CASP3', 'CASP9', 'APAF1', 'CYCS', 'BID', 'BAK1'] + [f'GENE{i}' for i in range(1, 12)],
    'GO:0008283': ['MYC', 'EGFR', 'ERBB2', 'KRAS', 'BRAF', 'PIK3CA', 'AKT1', 'CDK4', 'CDK6', 'CCND1'] + [f'GENE{i}' for i in range(12, 20)],
    'GO:0006281': ['TP53', 'BRCA1', 'ATM', 'CHEK2', 'PTEN'] + [f'GENE{i}' for i in range(20, 30)],
    'GO:0007165': ['EGFR', 'KRAS', 'BRAF', 'PIK3CA', 'AKT1', 'PTEN', 'MAPK1', 'MAPK3'] + [f'GENE{i}' for i in range(30, 42)],
    'GO:0006954': ['IL6', 'TNF', 'IFNG', 'IL1B', 'TGFB1'] + [f'GENE{i}' for i in range(42, 53)],
    'GO:0045765': ['VEGFA', 'FGF2', 'PDGFA', 'TGFB1'] + [f'GENE{i}' for i in range(53, 63)],
}

go_data = []
for go_id, genes in go_mapping.items():
    term_name, ontology = go_terms[go_id]
    for gene in genes:
        go_data.append({
            'GeneID': gene,
            'GOterm': go_id,
            'GOname': term_name,
            'Ontology': ontology
        })
go_annotation = pd.DataFrame(go_data)

# KEGG annotation
kegg_pathways = {
    'hsa04210': 'Apoptosis',
    'hsa04151': 'PI3K-Akt signaling pathway',
    'hsa05200': 'Pathways in cancer',
    'hsa04010': 'MAPK signaling pathway',
}

kegg_mapping = {
    'hsa04210': ['TP53', 'BCL2', 'BAX', 'CASP3', 'CASP9', 'BID', 'BAK1'] + [f'GENE{i}' for i in range(70, 82)],
    'hsa04151': ['EGFR', 'PIK3CA', 'AKT1', 'PTEN', 'VEGFA', 'FGF2', 'MTOR', 'TSC1'] + [f'GENE{i}' for i in range(82, 95)],
    'hsa05200': ['TP53', 'KRAS', 'BRAF', 'PIK3CA', 'EGFR', 'ERBB2', 'MYC', 'RB1', 'CDK4'] + [f'GENE{i}' for i in range(95, 108)],
    'hsa04010': ['KRAS', 'BRAF', 'NRAS', 'EGFR', 'MAPK1', 'MAPK3', 'RAF1'] + [f'GENE{i}' for i in range(108, 120)],
}

kegg_data = []
for pathway_id, genes in kegg_mapping.items():
    pathway_name = kegg_pathways[pathway_id]
    for gene in genes:
        kegg_data.append({
            'GeneID': gene,
            'Pathway': pathway_id,
            'PathwayName': pathway_name
        })
kegg_annotation = pd.DataFrame(kegg_data)

print(f"   GO annotation: {len(go_annotation)} gene-term pairs")
print(f"   KEGG annotation: {len(kegg_annotation)} gene-pathway pairs")

# ============================================================================
# 3. GO ENRICHMENT ANALYSIS
# ============================================================================
print("\n3. Running GO enrichment...")

go_result = ra.richGO(
    genes=de_genes,
    godata=go_annotation,
    pvalue=0.1,
    padj=0.2,
    padj_method='fdr_bh',  # Note: lowercase required
    minSize=3
)

print(f"   Found {len(go_result.result)} significant terms")
if len(go_result.result) > 0:
    print("\n   Top results:")
    print(go_result.result[['Term', 'Count', 'Pvalue', 'Padj', 'RichFactor']].head().to_string(index=False))

# ============================================================================
# 4. KEGG PATHWAY ENRICHMENT
# ============================================================================
print("\n4. Running KEGG enrichment...")

kegg_result = ra.richKEGG(
    genes=de_genes,
    kodata=kegg_annotation,
    pvalue=0.1,
    padj=0.2,
    padj_method='fdr_bh',
    minSize=3
)

print(f"   Found {len(kegg_result.result)} significant pathways")
if len(kegg_result.result) > 0:
    print("\n   Top results:")
    print(kegg_result.result[['Term', 'Count', 'Pvalue', 'Padj', 'RichFactor']].head().to_string(index=False))

# ============================================================================
# 5. VISUALIZATIONS
# ============================================================================
print("\n5. Creating visualizations...")

if len(go_result.result) > 0:
    # Bar plot
    fig = ra.ggbar(go_result, top=5, colorBy='Padj')
    print("   ✓ Bar plot created")
    # fig.savefig('go_barplot.png', dpi=300, bbox_inches='tight')

    # Dot plot
    fig = ra.ggdot(go_result, top=5, colorBy='Padj')
    print("   ✓ Dot plot created")

if len(kegg_result.result) > 0:
    # Network plot
    fig = ra.ggnetplot(kegg_result, top=5)
    print("   ✓ Network plot created")

# ============================================================================
# 6. ENRICHRESULT METHODS
# ============================================================================
print("\n6. Testing EnrichResult methods...")

if len(go_result.result) > 0:
    # Filter
    filtered = go_result.filter(pvalue=0.05)
    print(f"   Filtered (p<0.05): {len(filtered.result)} terms")

    # Get top
    top3 = go_result.top(n=3)
    print(f"   Top 3: {len(top3.result)} terms")

    # Export
    go_result.to_csv('go_results.csv')
    print("   ✓ Exported to CSV")

# ============================================================================
# 7. TERM CLUSTERING
# ============================================================================
print("\n7. Testing term clustering...")

if len(go_result.result) >= 2:
    clustered = ra.richCluster(go_result, cutoff=0.3, minSize=2)
    if 'Cluster' in clustered.result.columns:
        n_clusters = clustered.result['Cluster'].nunique()
        print(f"   Formed {n_clusters} clusters")
    else:
        print("   No clusters formed")

print("\n" + "=" * 80)
print("Example completed successfully!")
print("=" * 80)
print("\nNext steps:")
print("- Uncomment fig.savefig() lines to save plots")
print("- Check 'go_results.csv' for exported results")
print("- Try with your own gene lists and annotation data")
print("=" * 80)
