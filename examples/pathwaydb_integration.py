#!/usr/bin/env python3
"""
pathwaydb Integration Example for richAnn

This example demonstrates how to use richAnn with pathwaydb for
gene set enrichment analysis. pathwaydb provides easy access to
GO and KEGG annotations for 12 model species.

Requirements:
    pip install pathwaydb

Usage:
    python pathwaydb_integration.py
"""

import richAnn as ra

# Check if pathwaydb is installed
try:
    from pathwaydb import GO, KEGG, get_supported_species
    PATHWAYDB_AVAILABLE = True
except ImportError:
    PATHWAYDB_AVAILABLE = False
    print("=" * 60)
    print("pathwaydb is not installed!")
    print("Install it with: pip install pathwaydb")
    print("Or from source: git clone https://github.com/guokai8/pathwaydb.git")
    print("=" * 60)


def main():
    if not PATHWAYDB_AVAILABLE:
        print("\nRunning demo with synthetic data instead...\n")
        run_demo_without_pathwaydb()
        return

    print("=" * 60)
    print("richAnn + pathwaydb Integration Example")
    print("=" * 60)

    # Show supported species
    print("\nSupported species in pathwaydb:")
    species_list = get_supported_species()
    print(f"  {', '.join(species_list)}")

    # =========================================================================
    # Step 1: Download/Load GO Annotations
    # =========================================================================
    print("\n" + "=" * 60)
    print("Step 1: Setting up GO annotations")
    print("=" * 60)

    go = GO(storage_path='go_human.db')

    # Check if already downloaded
    try:
        stats = go.storage.stats() if go.storage else None
        if stats and stats.get('total_annotations', 0) > 0:
            print(f"Using existing GO database: {stats['total_annotations']:,} annotations")
        else:
            raise Exception("Empty database")
    except:
        print("Downloading GO annotations for human (this may take a few minutes)...")
        go = GO(storage_path='go_human.db')
        go.download_annotations(species='human')
        print("Download complete!")

    # =========================================================================
    # Step 2: Download/Load KEGG Annotations
    # =========================================================================
    print("\n" + "=" * 60)
    print("Step 2: Setting up KEGG annotations")
    print("=" * 60)

    kegg = KEGG(species='hsa', storage_path='kegg_human.db')

    try:
        stats = kegg.storage.stats() if kegg.storage else None
        if stats and stats.get('total_annotations', 0) > 0:
            print(f"Using existing KEGG database: {stats['total_annotations']:,} annotations")
        else:
            raise Exception("Empty database")
    except:
        print("Downloading KEGG annotations for human...")
        kegg = KEGG(species='hsa', storage_path='kegg_human.db')
        kegg.download_annotations()
        kegg.convert_ids_to_symbols()
        print("Download complete!")

    # =========================================================================
    # Step 3: Convert to richAnn Format
    # =========================================================================
    print("\n" + "=" * 60)
    print("Step 3: Converting to richAnn format")
    print("=" * 60)

    # Convert GO annotations (Biological Process only)
    go_data = ra.from_pathwaydb_go(go.storage, ontology="BP")
    print(f"GO (BP): {go_data['GOterm'].nunique():,} terms, {go_data['GeneID'].nunique():,} genes")

    # Convert KEGG annotations
    kegg_data = ra.from_pathwaydb_kegg(kegg.storage)
    print(f"KEGG: {kegg_data['Pathway'].nunique():,} pathways, {kegg_data['GeneID'].nunique():,} genes")

    # =========================================================================
    # Step 4: Define Gene List
    # =========================================================================
    print("\n" + "=" * 60)
    print("Step 4: Gene list for analysis")
    print("=" * 60)

    # Example: DNA damage response genes
    genes = [
        # Core DDR genes
        "TP53", "BRCA1", "BRCA2", "ATM", "ATR",
        # Checkpoint kinases
        "CHEK1", "CHEK2",
        # DNA repair
        "RAD51", "PALB2", "NBN", "MRE11", "RAD50",
        "XRCC1", "PARP1", "FANCD2", "FANCA",
        # Cell cycle
        "CDKN1A", "MDM2", "GADD45A"
    ]
    print(f"Analyzing {len(genes)} genes: {', '.join(genes[:5])}...")

    # =========================================================================
    # Step 5: GO Enrichment Analysis
    # =========================================================================
    print("\n" + "=" * 60)
    print("Step 5: GO Enrichment Analysis (Biological Process)")
    print("=" * 60)

    go_result = ra.richGO(
        genes=genes,
        godata=go_data,
        ontology="BP",
        pvalue=0.05,
        padj=0.1,
        minSize=5,
        maxSize=500
    )

    print("\nGO Enrichment Summary:")
    go_result.summary()

    if len(go_result.result) > 0:
        print("\nTop 10 Enriched GO Terms:")
        print("-" * 80)
        for i, row in go_result.top(10).result.iterrows():
            print(f"  {row['Annot']}: {row['Term'][:50]}")
            print(f"    Padj={row['Padj']:.2e}, Count={row['Count']}, RichFactor={row['RichFactor']:.2f}")

    # =========================================================================
    # Step 6: KEGG Enrichment Analysis
    # =========================================================================
    print("\n" + "=" * 60)
    print("Step 6: KEGG Pathway Enrichment Analysis")
    print("=" * 60)

    kegg_result = ra.richKEGG(
        genes=genes,
        kodata=kegg_data,
        pvalue=0.05,
        padj=0.1,
        minSize=5,
        maxSize=500
    )

    print("\nKEGG Enrichment Summary:")
    kegg_result.summary()

    if len(kegg_result.result) > 0:
        print("\nTop 10 Enriched KEGG Pathways:")
        print("-" * 80)
        for i, row in kegg_result.top(10).result.iterrows():
            print(f"  {row['Annot']}: {row['Term'][:50]}")
            print(f"    Padj={row['Padj']:.2e}, Count={row['Count']}, RichFactor={row['RichFactor']:.2f}")

    # =========================================================================
    # Step 7: Visualizations
    # =========================================================================
    print("\n" + "=" * 60)
    print("Step 7: Creating Visualizations")
    print("=" * 60)

    if len(go_result.result) > 0:
        # GO dot plot
        print("Creating GO dot plot...")
        fig_go = ra.ggdot(go_result, top=15, usePadj=True)
        fig_go.savefig("pathwaydb_go_dotplot.png", dpi=300, bbox_inches='tight')
        print("  Saved: pathwaydb_go_dotplot.png")

        # GO bar plot
        print("Creating GO bar plot...")
        fig_bar = ra.ggbar(go_result, top=10, colorBy='Padj')
        fig_bar.savefig("pathwaydb_go_barplot.png", dpi=300, bbox_inches='tight')
        print("  Saved: pathwaydb_go_barplot.png")

    if len(kegg_result.result) > 0:
        # KEGG dot plot
        print("Creating KEGG dot plot...")
        fig_kegg = ra.ggdot(kegg_result, top=10, usePadj=True)
        fig_kegg.savefig("pathwaydb_kegg_dotplot.png", dpi=300, bbox_inches='tight')
        print("  Saved: pathwaydb_kegg_dotplot.png")

    # Combined network (if both have results)
    if len(go_result.result) > 0 and len(kegg_result.result) > 0:
        print("Creating combined network plot...")
        try:
            fig_net = ra.ggnetmap(
                {'GO': go_result, 'KEGG': kegg_result},
                top=8
            )
            fig_net.savefig("pathwaydb_combined_network.png", dpi=300, bbox_inches='tight')
            print("  Saved: pathwaydb_combined_network.png")
        except Exception as e:
            print(f"  Network plot skipped: {e}")

    # =========================================================================
    # Step 8: Export Results
    # =========================================================================
    print("\n" + "=" * 60)
    print("Step 8: Exporting Results")
    print("=" * 60)

    if len(go_result.result) > 0:
        go_result.to_csv("pathwaydb_go_results.csv")
        print("Saved: pathwaydb_go_results.csv")

    if len(kegg_result.result) > 0:
        kegg_result.to_csv("pathwaydb_kegg_results.csv")
        print("Saved: pathwaydb_kegg_results.csv")

    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)


def run_demo_without_pathwaydb():
    """Demo with synthetic data when pathwaydb is not installed."""
    import pandas as pd

    print("Creating synthetic annotation data...")

    # Create synthetic GO data
    go_data = pd.DataFrame({
        'GeneID': ['TP53', 'BRCA1', 'BRCA2', 'ATM', 'TP53', 'BRCA1', 'CHEK2', 'RAD51'] * 10,
        'GOterm': ['GO:0006281', 'GO:0006281', 'GO:0006281', 'GO:0006281',
                   'GO:0006915', 'GO:0006915', 'GO:0006915', 'GO:0006915'] * 10,
        'GOname': ['DNA repair', 'DNA repair', 'DNA repair', 'DNA repair',
                   'apoptotic process', 'apoptotic process', 'apoptotic process', 'apoptotic process'] * 10,
        'Ontology': ['BP'] * 80
    }).drop_duplicates()

    genes = ["TP53", "BRCA1", "BRCA2", "ATM", "CHEK2", "RAD51"]

    print(f"Running GO enrichment with {len(genes)} genes...")
    result = ra.richGO(genes, go_data, ontology="BP", pvalue=0.5, padj=0.5)

    print("\nResults:")
    result.summary()

    print("\nTo use real data, install pathwaydb:")
    print("  pip install pathwaydb")


if __name__ == "__main__":
    main()
