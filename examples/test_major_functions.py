"""
Comprehensive example demonstrating major richAnn functions

This script tests:
1. Gene set enrichment analysis (GO, KEGG, GSEA)
2. Term clustering
3. Multiple visualization types
4. Cross-sample comparison

Requirements:
- richAnn installed (pip install -e .)
- Sample annotation data (GO, KEGG, or GMT files)
"""

import pandas as pd
import numpy as np
import richAnn as ra
import logging

# Set up logging to see what's happening
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

print("=" * 80)
print("richAnn Major Functions Test")
print("=" * 80)

# ============================================================================
# 1. PREPARE SAMPLE DATA
# ============================================================================
print("\n1. Preparing sample data...")

# Create a larger gene universe first
np.random.seed(42)
all_genes_universe = [
    'TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS', 'PTEN', 'AKT1', 'PIK3CA',
    'BRAF', 'NRAS', 'ERBB2', 'RB1', 'CDKN2A', 'ATM', 'CHEK2', 'MDM2',
    'BCL2', 'BAX', 'CASP3', 'CASP9', 'APAF1', 'CYCS', 'BID', 'BAK1',
    'VEGFA', 'FGF2', 'PDGFA', 'TGFB1', 'IL6', 'TNF', 'IFNG', 'IL1B',
    'STAT3', 'JAK2', 'SMAD4', 'NOTCH1', 'WNT1', 'CTNNB1', 'GSK3B', 'AXIN1',
    'CDK4', 'CDK6', 'CCND1', 'E2F1', 'CDKN1A', 'CDKN1B', 'SKP2', 'CUL1',
    'MAPK1', 'MAPK3', 'RAF1', 'MAP2K1', 'MAP2K2', 'SOS1', 'GRB2', 'SHC1',
    'MTOR', 'RICTOR', 'RAPTOR', 'TSC1', 'TSC2', 'RHEB', 'AKT2', 'AKT3',
    'FOXO1', 'FOXO3', 'GSK3A', 'IRS1', 'IRS2', 'INSR', 'IGF1R', 'PDGFRA',
    'PDGFRB', 'KIT', 'FLT3', 'CSF1R', 'MET', 'RET', 'ALK', 'ROS1',
    'HRAS', 'RRAS', 'MRAS', 'RAC1', 'CDC42', 'RHOA', 'PAK1', 'PAK2',
    'JUN', 'FOS', 'ATF2', 'CREB1', 'ELK1', 'SRF', 'MYB', 'ETS1',
    'NFKB1', 'NFKB2', 'RELA', 'RELB', 'IKBKA', 'IKBKB', 'IKBKG', 'CHUK',
    'TRAF2', 'TRAF6', 'TRADD', 'FADD', 'RIPK1', 'CASP8', 'CFLAR', 'BIRC2',
    'XIAP', 'BIRC3', 'DIABLO', 'HTRA2', 'BCL2L1', 'MCL1', 'BCL2L11', 'BMF',
    'PUMA', 'NOXA', 'HRK', 'BIK', 'BAD', 'BIM', 'PMAIP1', 'BBC3',
] + [f'GENE{i}' for i in range(1, 151)]

# Sample gene list for enrichment (subset of universe - differentially expressed)
sample_genes = [
    'TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS', 'PTEN', 'AKT1', 'PIK3CA',
    'BRAF', 'NRAS', 'ERBB2', 'RB1', 'CDKN2A', 'ATM', 'CHEK2', 'MDM2',
    'BCL2', 'BAX', 'CASP3', 'CASP9', 'APAF1', 'CYCS', 'BID', 'BAK1',
    'VEGFA', 'FGF2', 'PDGFA', 'TGFB1', 'IL6', 'TNF', 'IFNG', 'IL1B'
]

# Sample ranked gene list for GSEA (gene: score pairs)
# Higher scores = upregulated, lower scores = downregulated
gene_scores = pd.Series(
    np.random.randn(len(all_genes_universe)),
    index=all_genes_universe
).sort_values(ascending=False)

print(f"  - Sample gene list: {len(sample_genes)} genes")
print(f"  - Gene universe: {len(all_genes_universe)} genes")
print(f"  - Ranked gene list for GSEA: {len(gene_scores)} genes")

# ============================================================================
# 2. CREATE SAMPLE ANNOTATION DATA
# ============================================================================
print("\n2. Creating sample annotation data...")

# Create sample GO annotation
go_terms = {
    'GO:0006915': ('apoptotic process', 'BP'),
    'GO:0008283': ('cell proliferation', 'BP'),
    'GO:0006281': ('DNA repair', 'BP'),
    'GO:0007165': ('signal transduction', 'BP'),
    'GO:0006954': ('inflammatory response', 'BP'),
    'GO:0045765': ('regulation of angiogenesis', 'BP'),
    'GO:0071456': ('cellular response to hypoxia', 'BP'),
    'GO:0006468': ('protein phosphorylation', 'BP'),
}

go_gene_mapping = {
    'GO:0006915': ['TP53', 'BCL2', 'BAX', 'CASP3', 'CASP9', 'APAF1', 'CYCS', 'BID', 'BAK1',
                   'CASP8', 'FADD', 'TRADD', 'BCL2L1', 'MCL1', 'BCL2L11', 'BAD', 'BIM', 'PMAIP1', 'BBC3'],
    'GO:0008283': ['MYC', 'EGFR', 'ERBB2', 'KRAS', 'BRAF', 'PIK3CA', 'AKT1',
                   'CDK4', 'CDK6', 'CCND1', 'E2F1', 'MTOR', 'FOXO1', 'GENE1', 'GENE2', 'GENE3'],
    'GO:0006281': ['TP53', 'BRCA1', 'ATM', 'CHEK2', 'PTEN', 'GENE4', 'GENE5', 'GENE6'],
    'GO:0007165': ['EGFR', 'KRAS', 'BRAF', 'PIK3CA', 'AKT1', 'PTEN', 'MAPK1', 'MAPK3', 'RAF1',
                   'MAP2K1', 'GRB2', 'SOS1', 'GENE7', 'GENE8'],
    'GO:0006954': ['IL6', 'TNF', 'IFNG', 'IL1B', 'TGFB1', 'NFKB1', 'RELA', 'IKBKA', 'TRAF2', 'GENE9', 'GENE10'],
    'GO:0045765': ['VEGFA', 'FGF2', 'PDGFA', 'TGFB1', 'PDGFRA', 'PDGFRB', 'GENE11', 'GENE12'],
    'GO:0071456': ['VEGFA', 'EGFR', 'AKT1', 'MTOR', 'GENE13', 'GENE14', 'GENE15'],
    'GO:0006468': ['EGFR', 'AKT1', 'BRAF', 'ATM', 'CHEK2', 'MAPK1', 'MAPK3', 'CDK4', 'CDK6', 'GENE16', 'GENE17'],
}

# Build GO annotation DataFrame
go_data = []
for go_id, genes in go_gene_mapping.items():
    term_name, ontology = go_terms[go_id]
    for gene in genes:
        go_data.append({
            'GeneID': gene,
            'GOterm': go_id,
            'GOname': term_name,
            'Ontology': ontology
        })

go_annotation = pd.DataFrame(go_data)
print(f"  - GO annotation: {len(go_annotation)} gene-term pairs")

# Create sample KEGG annotation
kegg_pathways = {
    'hsa04210': 'Apoptosis',
    'hsa04151': 'PI3K-Akt signaling pathway',
    'hsa05200': 'Pathways in cancer',
    'hsa04010': 'MAPK signaling pathway',
    'hsa04060': 'Cytokine-cytokine receptor interaction',
}

kegg_gene_mapping = {
    'hsa04210': ['TP53', 'BCL2', 'BAX', 'CASP3', 'CASP9', 'BID', 'BAK1', 'CASP8', 'FADD', 'TRADD',
                 'BCL2L1', 'GENE20', 'GENE21'],
    'hsa04151': ['EGFR', 'PIK3CA', 'AKT1', 'PTEN', 'VEGFA', 'FGF2', 'PDGFA', 'MTOR', 'RICTOR',
                 'TSC1', 'TSC2', 'FOXO1', 'IRS1', 'GENE22', 'GENE23', 'GENE24'],
    'hsa05200': ['TP53', 'KRAS', 'BRAF', 'PIK3CA', 'EGFR', 'ERBB2', 'MYC', 'RB1', 'CDK4', 'CDKN2A',
                 'VEGFA', 'GENE25', 'GENE26', 'GENE27'],
    'hsa04010': ['KRAS', 'BRAF', 'NRAS', 'EGFR', 'FGF2', 'TGFB1', 'MAPK1', 'MAPK3', 'RAF1', 'MAP2K1',
                 'JUN', 'FOS', 'GENE28', 'GENE29'],
    'hsa04060': ['IL6', 'TNF', 'IFNG', 'IL1B', 'VEGFA', 'FGF2', 'PDGFA', 'TGFB1', 'CSF1R', 'KIT',
                 'GENE30', 'GENE31', 'GENE32'],
}

kegg_data = []
for pathway_id, genes in kegg_gene_mapping.items():
    pathway_name = kegg_pathways[pathway_id]
    for gene in genes:
        kegg_data.append({
            'GeneID': gene,
            'Pathway': pathway_id,
            'PathwayName': pathway_name
        })

kegg_annotation = pd.DataFrame(kegg_data)
print(f"  - KEGG annotation: {len(kegg_annotation)} gene-pathway pairs")

# ============================================================================
# 3. GENE ONTOLOGY ENRICHMENT (richGO)
# ============================================================================
print("\n3. Running GO enrichment analysis...")

try:
    go_result = ra.richGO(
        genes=sample_genes,
        godata=go_annotation,
        pvalue=0.1,
        padj=0.2,
        padj_method='fdr_bh'  # Note: must be lowercase
    )

    print(f"  - Significant GO terms: {len(go_result.result)}")
    if len(go_result.result) > 0:
        print("\n  Top 3 enriched GO terms:")
        print(go_result.result[['Term', 'Pvalue', 'Padj', 'Count', 'RichFactor']].head(3).to_string(index=False))
except Exception as e:
    import traceback
    print(f"  ERROR in richGO: {e}")
    traceback.print_exc()
    go_result = None

# ============================================================================
# 4. KEGG PATHWAY ENRICHMENT (richKEGG)
# ============================================================================
print("\n4. Running KEGG pathway enrichment...")

try:
    kegg_result = ra.richKEGG(
        genes=sample_genes,
        kodata=kegg_annotation,
        pvalue=0.1,
        padj=0.2,
        padj_method='fdr_bh'
    )

    print(f"  - Significant KEGG pathways: {len(kegg_result.result)}")
    if len(kegg_result.result) > 0:
        print("\n  Top 3 enriched pathways:")
        print(kegg_result.result[['Term', 'Pvalue', 'Padj', 'Count', 'RichFactor']].head(3).to_string(index=False))
except Exception as e:
    print(f"  ERROR in richKEGG: {e}")
    kegg_result = None

# ============================================================================
# 5. GENE SET ENRICHMENT ANALYSIS (richGSEA)
# ============================================================================
print("\n5. Running GSEA...")

# Create gene set annotation for GSEA (using GO terms)
gsea_data = []
for go_id, genes in go_gene_mapping.items():
    term_name, _ = go_terms[go_id]
    for gene in genes:
        gsea_data.append({
            'GeneSet': go_id,
            'GeneSetName': term_name,
            'GeneID': gene
        })

gsea_annotation = pd.DataFrame(gsea_data)

try:
    gsea_result = ra.richGSEA(
        gene_scores=gene_scores.to_dict(),
        geneset_db=gsea_annotation,
        nperm=1000,
        min_size=3
    )

    print(f"  - Significant gene sets: {len(gsea_result.result)}")
    if len(gsea_result.result) > 0:
        print("\n  Top 3 enriched gene sets:")
        print(gsea_result.result[['Term', 'ES', 'NES', 'Pvalue', 'Padj']].head(3).to_string(index=False))
except Exception as e:
    print(f"  ERROR in richGSEA: {e}")
    gsea_result = None

# ============================================================================
# 6. TERM CLUSTERING (richCluster)
# ============================================================================
print("\n6. Running term clustering...")

if go_result and len(go_result.result) > 0:
    try:
        clustered_result = ra.richCluster(
            go_result,
            cutoff=0.3,
            minSize=2
        )

        if 'Cluster' in clustered_result.result.columns:
            n_clusters = clustered_result.result['Cluster'].nunique()
            print(f"  - Number of clusters: {n_clusters}")
            print(f"  - Clustered terms: {len(clustered_result.result)}")
        else:
            print("  - No clusters formed (terms too dissimilar)")
    except Exception as e:
        print(f"  ERROR in richCluster: {e}")
        clustered_result = go_result
else:
    print("  - Skipping (no GO results)")
    clustered_result = None

# ============================================================================
# 7. VISUALIZATIONS
# ============================================================================
print("\n7. Creating visualizations...")

# 7a. Bar plot
if go_result and len(go_result.result) > 0:
    try:
        fig = ra.ggbar(go_result, top=5, colorBy='Padj')
        print("  - Bar plot created successfully")
        # Uncomment to save: fig.savefig('go_barplot.png', dpi=300, bbox_inches='tight')
    except Exception as e:
        print(f"  ERROR in ggbar: {e}")

# 7b. Dot plot
if kegg_result and len(kegg_result.result) > 0:
    try:
        fig = ra.ggdot(kegg_result, top=5, colorBy='Padj')
        print("  - Dot plot created successfully")
        # Uncomment to save: fig.savefig('kegg_dotplot.png', dpi=300, bbox_inches='tight')
    except Exception as e:
        print(f"  ERROR in ggdot: {e}")

# 7c. Network plot (term-gene bipartite network)
if go_result and len(go_result.result) > 0:
    try:
        fig = ra.ggnetplot(go_result, top=5)
        print("  - Network plot (term-gene) created successfully")
        # Uncomment to save: fig.savefig('go_network.png', dpi=300, bbox_inches='tight')
    except Exception as e:
        print(f"  ERROR in ggnetplot: {e}")

# 7d. Term-term network (based on gene overlap)
if go_result and len(go_result.result) > 0:
    try:
        fig = ra.ggnetwork(go_result, top=5, weightcut=0.2)
        print("  - Term-term network created successfully")
        # Uncomment to save: fig.savefig('go_term_network.png', dpi=300, bbox_inches='tight')
    except Exception as e:
        print(f"  ERROR in ggnetwork: {e}")

# 7e. GSEA enrichment plot
if gsea_result and len(gsea_result.result) > 0:
    try:
        # Plot the top enriched gene set
        top_geneset = gsea_result.result.iloc[0]['Annot']
        fig = ra.ggGSEA(gsea_result, geneSet=top_geneset)
        print("  - GSEA enrichment plot created successfully")
        # Uncomment to save: fig.savefig('gsea_plot.png', dpi=300, bbox_inches='tight')
    except Exception as e:
        print(f"  ERROR in ggGSEA: {e}")

# ============================================================================
# 8. CROSS-SAMPLE COMPARISON
# ============================================================================
print("\n8. Testing cross-sample comparison...")

# Create a second sample with slightly different genes
sample_genes_2 = sample_genes[:20] + ['STAT3', 'JAK2', 'SMAD4', 'NOTCH1']

try:
    go_result_2 = ra.richGO(
        genes=sample_genes_2,
        godata=go_annotation,
        pvalue=0.1,
        padj=0.2,
        padj_method='fdr_bh'
    )

    # Compare results
    if go_result and go_result_2:
        compared = ra.compareResult(
            [go_result, go_result_2],
            sample_names=['Sample1', 'Sample2']
        )
        print(f"  - Comparison table created: {compared.shape}")

        # Heatmap comparison
        try:
            fig = ra.ggheatmap(compared, usePadj=True)
            print("  - Comparison heatmap created successfully")
            # Uncomment to save: fig.savefig('comparison_heatmap.png', dpi=300, bbox_inches='tight')
        except Exception as e:
            print(f"  ERROR in ggheatmap: {e}")

        # Dot plot comparison
        try:
            fig = ra.comparedot(compared, usePadj=True)
            print("  - Comparison dot plot created successfully")
            # Uncomment to save: fig.savefig('comparison_dotplot.png', dpi=300, bbox_inches='tight')
        except Exception as e:
            print(f"  ERROR in comparedot: {e}")

except Exception as e:
    print(f"  ERROR in comparison: {e}")

# ============================================================================
# 9. ENRICHRESULT METHODS
# ============================================================================
print("\n9. Testing EnrichResult methods...")

if go_result and len(go_result.result) > 0:
    # Filter by p-value
    filtered = go_result.filter(pvalue=0.01)
    print(f"  - Filtered results (p<0.01): {len(filtered.result)} terms")

    # Get top results
    top_results = go_result.top(n=3)
    print(f"  - Top 3 results: {len(top_results.result)} terms")

    # Export to CSV
    try:
        go_result.to_csv('go_enrichment_results.csv')
        print("  - Results exported to CSV successfully")
    except Exception as e:
        print(f"  ERROR exporting to CSV: {e}")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("Test Summary")
print("=" * 80)
print(f"GO enrichment: {'✓' if go_result else '✗'}")
print(f"KEGG enrichment: {'✓' if kegg_result else '✗'}")
print(f"GSEA: {'✓' if gsea_result else '✗'}")
print(f"Clustering: {'✓' if clustered_result else '✗'}")
print("\nAll major functions tested!")
print("\nTo save plots, uncomment the fig.savefig() lines in the script.")
print("=" * 80)
