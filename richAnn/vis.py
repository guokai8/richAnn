"""
Visualization functions for richAnn package
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from typing import Dict, Optional
import logging

from .core import EnrichResult

logger = logging.getLogger(__name__)


def ggbar(enrich_result: EnrichResult,
          top: int = 10,
          pvalue: float = 0.05,
          padj: Optional[float] = None,
          usePadj: bool = True,
          orderBy: str = 'Padj',
          colorBy: str = 'Padj',
          palette: str = 'RdYlBu_r',
          horiz: bool = True,
          useTerm: bool = True,
          fontsize: int = 10,
          figsize: tuple = (10, 6),
          filename: Optional[str] = None,
          dpi: int = 300,
          **kwargs):
    """
    Create barplot for enrichment results
    
    Parameters:
    -----------
    enrich_result : EnrichResult
        Enrichment result object
    top : int
        Number of top terms to display
    pvalue : float
        P-value cutoff for filtering
    padj : float
        Adjusted p-value cutoff
    usePadj : bool
        Use adjusted p-value for coloring
    orderBy : str
        Column to order by
    colorBy : str
        Column to color by
    palette : str
        Color palette
    horiz : bool
        Horizontal bars
    useTerm : bool
        Use term names vs IDs
    fontsize : int
        Font size
    figsize : tuple
        Figure size
    filename : str
        Save to file
    dpi : int
        Resolution
    
    Returns:
    --------
    matplotlib Figure object
    """
    df = enrich_result.result.copy()
    
    if padj is not None:
        df = df[df['Padj'] <= padj]
    else:
        df = df[df['Pvalue'] <= pvalue]
    
    if len(df) == 0:
        raise ValueError("No terms remain after filtering")
    
    ascending = orderBy in ['Pvalue', 'Padj']
    df = df.sort_values(orderBy, ascending=ascending).head(top)
    df = df.iloc[::-1]
    
    labels = df['Term'].values if useTerm else df['Annot'].values
    
    if colorBy == 'Padj' or (colorBy == 'Pvalue' and usePadj):
        colors = -np.log10(df['Padj'].values)
        cbar_label = '-log10(Adjusted P-value)'
    elif colorBy == 'Pvalue':
        colors = -np.log10(df['Pvalue'].values)
        cbar_label = '-log10(P-value)'
    elif colorBy == 'RichFactor':
        colors = df['RichFactor'].values
        cbar_label = 'Rich Factor'
    else:
        colors = df[colorBy].values
        cbar_label = colorBy
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if horiz:
        bars = ax.barh(range(len(df)), df['Count'].values, height=0.7)
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(labels, fontsize=fontsize)
        ax.set_xlabel('Gene Count', fontsize=fontsize+2, fontweight='bold')
        ax.invert_yaxis()
    else:
        bars = ax.bar(range(len(df)), df['Count'].values, width=0.7)
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=fontsize)
        ax.set_ylabel('Gene Count', fontsize=fontsize+2, fontweight='bold')
    
    # Handle edge case where all colors are identical
    if colors.min() == colors.max():
        # Use a single color for all bars
        single_color = plt.cm.get_cmap(palette)(0.5)
        for bar in bars:
            bar.set_color(single_color)
            bar.set_edgecolor('black')
            bar.set_linewidth(0.5)
        # Add a simple colorbar with single value
        norm = plt.Normalize(vmin=colors.min(), vmax=colors.min() + 1)
        sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap(palette), norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label(f'{cbar_label} (all {colors.min():.2f})', fontsize=fontsize)
    else:
        norm = plt.Normalize(vmin=colors.min(), vmax=colors.max())
        cmap = plt.cm.get_cmap(palette)

        for bar, color_val in zip(bars, colors):
            bar.set_color(cmap(norm(color_val)))
            bar.set_edgecolor('black')
            bar.set_linewidth(0.5)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label(cbar_label, fontsize=fontsize)
    
    ax.set_title(f'Top {len(df)} Enriched Terms', fontsize=fontsize+4, fontweight='bold')
    ax.grid(axis='x' if horiz else 'y', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {filename}")
    
    return fig


def ggdot(enrich_result: EnrichResult,
          top: int = 10,
          pvalue: float = 0.05,
          padj: Optional[float] = None,
          usePadj: bool = True,
          orderBy: str = 'Padj',
          palette: str = 'RdYlBu_r',
          useTerm: bool = True,
          fontsize: int = 10,
          size_range: tuple = (50, 500),
          figsize: tuple = (10, 6),
          filename: Optional[str] = None,
          dpi: int = 300,
          **kwargs):
    """
    Create dot plot for enrichment results
    
    X-axis: RichFactor
    Dot size: Gene count
    Dot color: P-value significance
    """
    df = enrich_result.result.copy()
    
    if padj is not None:
        df = df[df['Padj'] <= padj]
    else:
        df = df[df['Pvalue'] <= pvalue]
    
    if len(df) == 0:
        raise ValueError("No terms remain after filtering")
    
    ascending = orderBy in ['Pvalue', 'Padj']
    df = df.sort_values(orderBy, ascending=ascending).head(top)
    df = df.iloc[::-1]
    
    labels = df['Term'].values if useTerm else df['Annot'].values
    rich_factors = df['RichFactor'].values
    
    if usePadj:
        colors = -np.log10(df['Padj'].values)
        cbar_label = '-log10(Adjusted P-value)'
    else:
        colors = -np.log10(df['Pvalue'].values)
        cbar_label = '-log10(P-value)'
    
    counts = df['Count'].values
    sizes = np.interp(counts, (counts.min(), counts.max()), size_range)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    scatter = ax.scatter(rich_factors, range(len(df)), 
                        s=sizes,
                        c=colors, 
                        cmap=palette,
                        alpha=0.8,
                        edgecolors='black',
                        linewidth=1)
    
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(labels, fontsize=fontsize)
    ax.set_xlabel('Rich Factor', fontsize=fontsize+2, fontweight='bold')
    ax.axvline(x=1.0, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='RichFactor=1')
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(cbar_label, fontsize=fontsize)
    
    size_values = [counts.min(), np.median(counts), counts.max()]
    size_labels = [f'{int(v)}' for v in size_values]
    
    legend_elements = []
    for val, label in zip(size_values, size_labels):
        size = np.interp(val, (counts.min(), counts.max()), size_range)
        legend_elements.append(
            plt.scatter([], [], s=size, c='gray', alpha=0.6,
                       edgecolors='black', linewidth=0.5, label=label)
        )
    
    ax.legend(handles=legend_elements, scatterpoints=1, frameon=True,
             labelspacing=2, title='Gene Count', loc='lower right',
             fontsize=fontsize-2)
    
    ax.set_title(f'Top {len(df)} Enriched Terms', fontsize=fontsize+4, fontweight='bold')
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {filename}")
    
    return fig
# Continue in richAnn/visualization.py


def ggnetplot(enrich_result: EnrichResult,
              top: int = 20,
              pvalue: float = 0.05,
              padj: Optional[float] = None,
              usePadj: bool = True,
              layout: str = 'spring',
              node_size_scale: float = 300,
              label_size: int = 8,
              edge_alpha: float = 0.2,
              figsize: tuple = (14, 10),
              filename: Optional[str] = None,
              dpi: int = 300,
              seed: int = 42,
              **kwargs):
    """
    Create term-gene network plot
    
    Shows bipartite network: terms connected to their genes
    """
    df = enrich_result.result.copy()
    
    if padj is not None:
        df = df[df['Padj'] <= padj]
    else:
        df = df[df['Pvalue'] <= pvalue]
    
    df = df.sort_values('Padj').head(top)
    
    G = nx.Graph()
    
    term_nodes = []
    gene_nodes = set()
    
    for idx, row in df.iterrows():
        term_id = row['Annot']
        term_nodes.append(term_id)
        G.add_node(term_id, node_type='term',
                  pvalue=row['Padj'] if usePadj else row['Pvalue'],
                  count=row['Count'],
                  term=row['Term'])
        
        genes = str(row['GeneID']).replace(',', ';').split(';')
        for gene in genes:
            gene = gene.strip()
            if gene:
                gene_nodes.add(gene)
                G.add_node(gene, node_type='gene')
                G.add_edge(term_id, gene)
    
    np.random.seed(seed)
    if layout == 'spring':
        pos = nx.spring_layout(G, k=1, iterations=50, seed=seed)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)
    elif layout == 'circular':
        pos = nx.circular_layout(G)
    elif layout == 'bipartite':
        pos = nx.bipartite_layout(G, term_nodes)
    else:
        pos = nx.random_layout(G, seed=seed)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    gene_pos = {k: v for k, v in pos.items() if k in gene_nodes}
    nx.draw_networkx_nodes(G, gene_pos,
                          nodelist=list(gene_nodes),
                          node_color='lightgray',
                          node_size=50,
                          alpha=0.6,
                          ax=ax)
    
    term_pos = {k: v for k, v in pos.items() if k in term_nodes}
    term_colors = [-np.log10(G.nodes[t]['pvalue']) for t in term_nodes]
    term_sizes = [G.nodes[t]['count'] * node_size_scale for t in term_nodes]
    
    nodes = nx.draw_networkx_nodes(G, term_pos,
                                   nodelist=term_nodes,
                                   node_color=term_colors,
                                   node_size=term_sizes,
                                   cmap='RdYlBu_r',
                                   alpha=0.8,
                                   edgecolors='black',
                                   linewidths=2,
                                   ax=ax)
    
    nx.draw_networkx_edges(G, pos, alpha=edge_alpha, width=0.5, ax=ax)
    
    term_labels = {t: G.nodes[t]['term'][:30] + '...' if len(G.nodes[t]['term']) > 30
                   else G.nodes[t]['term'] for t in term_nodes}
    nx.draw_networkx_labels(G, term_pos, term_labels,
                           font_size=label_size,
                           font_weight='bold',
                           ax=ax)

    # Handle edge case where all p-values are identical
    if min(term_colors) == max(term_colors):
        vmin, vmax = min(term_colors), min(term_colors) + 1
    else:
        vmin, vmax = min(term_colors), max(term_colors)

    sm = plt.cm.ScalarMappable(cmap='RdYlBu_r',
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('-log10(Adjusted P-value)' if usePadj else '-log10(P-value)',
                  fontsize=10)
    
    ax.set_title('Term-Gene Network', fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {filename}")
    
    return fig


def ggnetwork(enrich_result: EnrichResult,
              top: int = 30,
              pvalue: float = 0.05,
              padj: Optional[float] = None,
              usePadj: bool = True,
              weightcut: float = 0.2,
              layout: str = 'spring',
              node_size_scale: float = 500,
              label_size: int = 8,
              edge_width_scale: float = 3,
              figsize: tuple = (14, 10),
              filename: Optional[str] = None,
              dpi: int = 300,
              seed: int = 42,
              **kwargs):
    """
    Create term-term network based on shared genes
    
    Edges represent gene overlap (Jaccard similarity)
    """
    df = enrich_result.result.copy()
    
    if padj is not None:
        df = df[df['Padj'] <= padj]
    else:
        df = df[df['Pvalue'] <= pvalue]
    
    df = df.sort_values('Padj').head(top)
    
    term_genesets = {}
    for idx, row in df.iterrows():
        genes = set(str(row['GeneID']).replace(',', ';').split(';'))
        genes = {g.strip() for g in genes if g.strip()}
        term_genesets[row['Annot']] = {
            'genes': genes,
            'pvalue': row['Padj'] if usePadj else row['Pvalue'],
            'count': row['Count'],
            'term': row['Term']
        }
    
    G = nx.Graph()
    
    for term_id, data in term_genesets.items():
        G.add_node(term_id, **data)
    
    term_ids = list(term_genesets.keys())
    for i in range(len(term_ids)):
        for j in range(i+1, len(term_ids)):
            term_i = term_ids[i]
            term_j = term_ids[j]
            
            genes_i = term_genesets[term_i]['genes']
            genes_j = term_genesets[term_j]['genes']
            
            intersection = len(genes_i.intersection(genes_j))
            union = len(genes_i.union(genes_j))
            
            if union > 0:
                jaccard = intersection / union
                
                if jaccard >= weightcut:
                    G.add_edge(term_i, term_j, weight=jaccard, overlap=intersection)
    
    np.random.seed(seed)
    if layout == 'spring':
        pos = nx.spring_layout(G, k=2, iterations=50, seed=seed)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)
    elif layout == 'circular':
        pos = nx.circular_layout(G)
    else:
        pos = nx.random_layout(G, seed=seed)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    
    nx.draw_networkx_edges(G, pos,
                          width=[w * edge_width_scale for w in weights],
                          alpha=0.4,
                          edge_color='gray',
                          ax=ax)
    
    node_colors = [-np.log10(G.nodes[n]['pvalue']) for n in G.nodes()]
    node_sizes = [G.nodes[n]['count'] * node_size_scale for n in G.nodes()]
    
    nodes = nx.draw_networkx_nodes(G, pos,
                                   node_color=node_colors,
                                   node_size=node_sizes,
                                   cmap='RdYlBu_r',
                                   alpha=0.8,
                                   edgecolors='black',
                                   linewidths=2,
                                   ax=ax)
    
    labels = {n: G.nodes[n]['term'][:25] + '...' if len(G.nodes[n]['term']) > 25
              else G.nodes[n]['term'] for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels,
                           font_size=label_size,
                           font_weight='bold',
                           ax=ax)

    # Handle edge case where all p-values are identical
    if min(node_colors) == max(node_colors):
        vmin, vmax = min(node_colors), min(node_colors) + 1
    else:
        vmin, vmax = min(node_colors), max(node_colors)

    sm = plt.cm.ScalarMappable(cmap='RdYlBu_r',
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('-log10(Adjusted P-value)' if usePadj else '-log10(P-value)',
                  fontsize=10)
    
    ax.set_title('Term-Term Network (Shared Genes)', fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {filename}")
    
    return fig


def ggnetmap(result_dict: Dict[str, EnrichResult],
             top: int = 15,
             weightcut: float = 0.2,
             layout: str = 'spring',
             node_size_scale: float = 400,
             label_size: int = 8,
             figsize: tuple = (16, 12),
             filename: Optional[str] = None,
             dpi: int = 300,
             seed: int = 42,
             **kwargs):
    """
    Create combined network from multiple enrichment results
    
    Parameters:
    -----------
    result_dict : dict
        Dictionary mapping database names to EnrichResult objects
        e.g., {'GO': go_result, 'KEGG': kegg_result}
    """
    all_terms = {}
    
    for db_name, result in result_dict.items():
        df = result.result.sort_values('Padj').head(top)
        
        for idx, row in df.iterrows():
            term_id = f"{db_name}:{row['Annot']}"
            genes = set(str(row['GeneID']).replace(',', ';').split(';'))
            genes = {g.strip() for g in genes if g.strip()}
            
            all_terms[term_id] = {
                'genes': genes,
                'pvalue': row['Padj'],
                'count': row['Count'],
                'term': row['Term'],
                'database': db_name
            }
    
    G = nx.Graph()
    
    for term_id, data in all_terms.items():
        G.add_node(term_id, **data)
    
    term_ids = list(all_terms.keys())
    for i in range(len(term_ids)):
        for j in range(i+1, len(term_ids)):
            term_i = term_ids[i]
            term_j = term_ids[j]
            
            genes_i = all_terms[term_i]['genes']
            genes_j = all_terms[term_j]['genes']
            
            intersection = len(genes_i.intersection(genes_j))
            union = len(genes_i.union(genes_j))
            
            if union > 0:
                jaccard = intersection / union
                
                if jaccard >= weightcut:
                    G.add_edge(term_i, term_j, weight=jaccard, overlap=intersection)
    
    np.random.seed(seed)
    if layout == 'spring':
        pos = nx.spring_layout(G, k=2, iterations=50, seed=seed)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)
    else:
        pos = nx.circular_layout(G)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    nx.draw_networkx_edges(G, pos,
                          width=[w * 3 for w in weights],
                          alpha=0.3,
                          edge_color='gray',
                          ax=ax)
    
    databases = list(result_dict.keys())
    colors_map = plt.cm.get_cmap('Set3', len(databases))
    db_colors = {db: colors_map(i) for i, db in enumerate(databases)}
    
    for db in databases:
        db_nodes = [n for n in G.nodes() if G.nodes[n]['database'] == db]
        if not db_nodes:
            continue
            
        node_sizes = [G.nodes[n]['count'] * node_size_scale for n in db_nodes]
        
        nx.draw_networkx_nodes(G, pos,
                              nodelist=db_nodes,
                              node_color=[db_colors[db]] * len(db_nodes),
                              node_size=node_sizes,
                              alpha=0.7,
                              edgecolors='black',
                              linewidths=2,
                              label=db,
                              ax=ax)
    
    labels = {n: G.nodes[n]['term'][:20] + '...' if len(G.nodes[n]['term']) > 20
              else G.nodes[n]['term'] for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels,
                           font_size=label_size,
                           ax=ax)
    
    ax.set_title('Combined Network Map', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12, loc='upper right')
    ax.axis('off')
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {filename}")
    
    return fig
# Continue in richAnn/visualization.py


def ggheatmap(compare_df: pd.DataFrame,
              pvalue: float = 0.05,
              top: int = 30,
              figsize: tuple = (12, 10),
              cmap: str = 'RdYlBu_r',
              fontsize: int = 10,
              cluster_rows: bool = True,
              cluster_cols: bool = False,
              filename: Optional[str] = None,
              dpi: int = 300,
              **kwargs):
    """
    Create heatmap for comparing enrichment across samples
    
    Parameters:
    -----------
    compare_df : pd.DataFrame
        DataFrame from compareResult() with columns: Term, Sample, Padj
    cluster_rows : bool
        Cluster rows (terms)
    cluster_cols : bool
        Cluster columns (samples)
    """
    df = compare_df[compare_df['Padj'] <= pvalue].copy()
    
    if len(df) == 0:
        raise ValueError("No terms remain after filtering")
    
    df['NegLogPadj'] = -np.log10(df['Padj'])
    
    pivot = df.pivot_table(index='Term', columns='Sample',
                          values='NegLogPadj', fill_value=0)
    
    pivot['total'] = pivot.sum(axis=1)
    pivot = pivot.sort_values('total', ascending=False).head(top)
    pivot = pivot.drop('total', axis=1)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if cluster_rows or cluster_cols:
        import scipy.cluster.hierarchy as sch
        
        if cluster_rows:
            row_linkage = sch.linkage(pivot.values, method='average', metric='euclidean')
            row_order = sch.leaves_list(row_linkage)
            pivot = pivot.iloc[row_order]
        
        if cluster_cols:
            col_linkage = sch.linkage(pivot.T.values, method='average', metric='euclidean')
            col_order = sch.leaves_list(col_linkage)
            pivot = pivot.iloc[:, col_order]
    
    sns.heatmap(pivot, cmap=cmap,
               cbar_kws={'label': '-log10(Adjusted P-value)'},
               linewidths=0.5, linecolor='gray',
               ax=ax, **kwargs)
    
    ax.set_xlabel('Sample', fontsize=fontsize+2, fontweight='bold')
    ax.set_ylabel('Term', fontsize=fontsize+2, fontweight='bold')
    ax.set_title('Enrichment Heatmap', fontsize=fontsize+4, fontweight='bold')
    
    plt.xticks(fontsize=fontsize, rotation=45, ha='right')
    plt.yticks(fontsize=fontsize-2)
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {filename}")
    
    return fig


def ggGSEA(enrich_result: EnrichResult,
           gene_scores: Dict[str, float],
           term_id: str,
           figsize: tuple = (10, 6),
           filename: Optional[str] = None,
           dpi: int = 300,
           **kwargs):
    """
    Create GSEA enrichment plot for a specific term
    
    Parameters:
    -----------
    enrich_result : EnrichResult
        GSEA result object
    gene_scores : dict
        Dictionary mapping gene IDs to scores (same as used in richGSEA)
    term_id : str
        Term/pathway ID to plot
    """
    df = enrich_result.result
    term_row = df[df['Annot'] == term_id]
    
    if len(term_row) == 0:
        raise ValueError(f"Term {term_id} not found in results")
    
    term_row = term_row.iloc[0]
    
    leading_edge = set(str(term_row['LeadingEdge']).split(';'))
    
    ranked_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)
    ranked_gene_list = [g[0] for g in ranked_genes]
    ranked_scores = np.array([g[1] for g in ranked_genes])
    
    N = len(ranked_gene_list)
    hit_indices = [i for i, g in enumerate(ranked_gene_list) if g in leading_edge]
    
    running_sum = np.zeros(N)
    abs_scores = np.abs(ranked_scores)
    N_R = np.sum(abs_scores[hit_indices])
    N_miss = N - len(hit_indices)
    
    for i in range(N):
        if i in hit_indices:
            if N_R > 0:
                running_sum[i] = abs_scores[i] / N_R
            else:
                running_sum[i] = 1.0 / len(hit_indices)
        else:
            if N_miss > 0:
                running_sum[i] = -1.0 / N_miss
    
    running_sum = np.cumsum(running_sum)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize,
                                   gridspec_kw={'height_ratios': [3, 1]})
    
    ax1.plot(range(N), running_sum, color='green', linewidth=2)
    ax1.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    
    max_idx = np.argmax(np.abs(running_sum))
    ax1.plot(max_idx, running_sum[max_idx], 'ro', markersize=8)
    
    ax1.set_ylabel('Enrichment Score', fontsize=12, fontweight='bold')
    ax1.set_title(f'{term_row["Term"]}\nES={term_row["ES"]:.3f}, NES={term_row["NES"]:.3f}, P={term_row["Pvalue"]:.4f}',
                 fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, N)
    
    for idx in hit_indices:
        ax2.axvline(x=idx, color='black', linewidth=0.5, alpha=0.5)
    
    ax2.fill_between(range(N), 0, 1, color='lightgray', alpha=0.3)
    ax2.set_xlim(0, N)
    ax2.set_ylim(0, 1)
    ax2.set_xlabel('Rank in Ordered Dataset', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Hits', fontsize=12, fontweight='bold')
    ax2.set_yticks([])
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {filename}")
    
    return fig


def ggcluster(data,
              top: int = 20,
              pvalue: float = 0.05,
              padj: Optional[float] = None,
              method: str = 'auto',
              group_col: Optional[str] = None,
              color_low: str = 'pink',
              color_high: str = 'red',
              color_mid: str = 'white',
              size_range: tuple = (20, 200),
              curve_color: str = 'grey',
              curve_linewidth: float = 0.5,
              vertical_line_color: str = 'darkcyan',
              vertical_line_width: float = 1.5,
              dot_line_color: str = 'grey',
              dot_line_width: float = 0.3,
              dot_line_style: str = 'dotted',
              vline_color: str = 'grey',
              vline_style: str = 'dashed',
              label_fontsize: int = 8,
              label_fontweight: str = 'bold',
              pathway_fontsize: int = 7,
              pathway_fontstyle: str = 'italic',
              x_label_fontsize: int = 8,
              figsize: tuple = (12, 10),
              filename: Optional[str] = None,
              dpi: int = 300,
              **kwargs):
    """
    Create KEGG cluster visualization plot

    This function creates a visualization showing pathways grouped by Level2 categories
    with cluster comparisons. Supports both enrichment and GSEA modes.

    Can accept either:
    - EnrichResult from richKEGG() (requires Level2 column in KEGG annotation data)
    - DataFrame with required columns

    Parameters:
    -----------
    data : EnrichResult or pd.DataFrame
        Enrichment results. For EnrichResult from richKEGG, requires Level2 column.
        For DataFrame:
        - Enrichment mode: 'Term', 'Level2', 'Padj', and either 'RichFactor' or
          ('Significant'/'Count' and 'Annotated')
        - GSEA mode: 'pathway', 'Level2', 'padj', 'NES'
    top : int
        Number of top terms to display (default: 20)
    pvalue : float
        P-value cutoff for filtering (default: 0.05)
    padj : float, optional
        Adjusted p-value cutoff (overrides pvalue if provided)
    method : str
        One of 'auto', 'enrich', or 'gsea'. If 'auto', guesses from columns.
    group_col : str, optional
        Column name for grouping (e.g., 'Cluster', 'Sample'). If None, uses 'group'
        column if present, otherwise creates a single default group.
    color_low : str
        Color for low values in gradient (default: 'pink')
    color_high : str
        Color for high values in gradient (default: 'red')
    color_mid : str
        Color for middle values in diverging gradient for GSEA (default: 'white')
    size_range : tuple
        Range for dot sizes (default: (20, 200))
    curve_color : str
        Color for connecting curves (default: 'grey')
    curve_linewidth : float
        Line width for curves (default: 0.5)
    vertical_line_color : str
        Color for vertical lines between Level2 and terms (default: 'darkcyan')
    vertical_line_width : float
        Width of vertical lines (default: 1.5)
    dot_line_color : str
        Color for dotted lines connecting clusters to pathways (default: 'grey')
    dot_line_width : float
        Width of dotted lines (default: 0.3)
    dot_line_style : str
        Style of dotted lines (default: 'dotted')
    vline_color : str
        Color for vertical lines separating clusters (default: 'grey')
    vline_style : str
        Style for cluster separator lines (default: 'dashed')
    label_fontsize : int
        Font size for Level2 labels (default: 8)
    label_fontweight : str
        Font weight for Level2 labels (default: 'bold')
    pathway_fontsize : int
        Font size for pathway labels (default: 7)
    pathway_fontstyle : str
        Font style for pathway labels (default: 'italic')
    x_label_fontsize : int
        Font size for x-axis labels (default: 8)
    figsize : tuple
        Figure size (default: (12, 10))
    filename : str, optional
        Path to save the figure
    dpi : int
        Resolution for saved figure (default: 300)

    Returns:
    --------
    matplotlib Figure object

    Examples:
    ---------
    >>> # Using with richKEGG results (requires Level2 in KEGG annotation)
    >>> kegg_result = ra.richKEGG(genes, kegg_data)
    >>> fig = ra.ggcluster(kegg_result, top=15)

    >>> # Using with DataFrame
    >>> import pandas as pd
    >>> data = pd.DataFrame({
    ...     'Term': ['Peroxisome', 'PPAR signaling pathway', 'Fatty acid elongation'],
    ...     'Level2': ['Transport and catabolism', 'Endocrine system', 'Lipid metabolism'],
    ...     'group': ['Cluster1', 'Cluster1', 'Cluster2'],
    ...     'Padj': [0.015, 0.04, 0.03],
    ...     'Count': [7, 5, 3],
    ...     'RichFactor': [2.1, 1.5, 1.8]
    ... })
    >>> fig = ggcluster(data)
    """
    from matplotlib.patches import FancyBboxPatch
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    # Handle EnrichResult input
    if hasattr(data, 'result'):
        df = data.result.copy()
    else:
        df = data.copy()

    # Filter by p-value
    if padj is not None:
        df = df[df['Padj'] <= padj]
    else:
        df = df[df['Pvalue'] <= pvalue] if 'Pvalue' in df.columns else df[df['Padj'] <= pvalue]

    if len(df) == 0:
        raise ValueError("No terms remain after filtering")

    # Select top terms
    df = df.sort_values('Padj').head(top)

    # Handle group column
    if group_col is not None and group_col in df.columns:
        # If using Cluster column with numeric values, add prefix
        if group_col == 'Cluster' and df[group_col].dtype in ['int64', 'int32', 'float64']:
            df['group'] = df[group_col].apply(lambda x: f'Cluster{int(x)}')
        else:
            df['group'] = df[group_col].astype(str)
    elif 'group' not in df.columns:
        if 'Cluster' in df.columns:
            df['group'] = df['Cluster'].apply(lambda x: f'Cluster{int(x)}' if pd.notna(x) else 'Unknown')
        elif 'Sample' in df.columns:
            df['group'] = df['Sample'].astype(str)
        else:
            df['group'] = 'All'

    # Check for Level2 column
    if 'Level2' not in df.columns:
        raise ValueError(
            "Level2 column not found. For richKEGG results, ensure your KEGG "
            "annotation data includes Level2 hierarchy information. "
            "You can add it manually or use pathwaydb which provides KEGG hierarchy."
        )

    # Auto-detect method based on GSEA-specific columns (NES, ES)
    if method == 'auto':
        if 'NES' in df.columns or 'ES' in df.columns:
            method = 'gsea'
        elif 'Term' in df.columns or 'pathway' in df.columns:
            method = 'enrich'
        else:
            raise ValueError("Cannot auto-detect method. Please specify method='enrich' or 'gsea'.")

    if method == 'enrich':
        return _ggcluster_enrich(
            df, color_low=color_low, color_high=color_high,
            size_range=size_range, curve_color=curve_color,
            curve_linewidth=curve_linewidth,
            vertical_line_color=vertical_line_color,
            vertical_line_width=vertical_line_width,
            dot_line_color=dot_line_color, dot_line_width=dot_line_width,
            dot_line_style=dot_line_style, vline_color=vline_color,
            vline_style=vline_style, label_fontsize=label_fontsize,
            label_fontweight=label_fontweight,
            pathway_fontsize=pathway_fontsize,
            pathway_fontstyle=pathway_fontstyle,
            x_label_fontsize=x_label_fontsize,
            figsize=figsize, filename=filename, dpi=dpi, **kwargs
        )
    else:
        return _ggcluster_gsea(
            df, color_low=color_low, color_high=color_high,
            color_mid=color_mid, size_range=size_range,
            curve_color=curve_color, curve_linewidth=curve_linewidth,
            vertical_line_color=vertical_line_color,
            vertical_line_width=vertical_line_width,
            dot_line_color=dot_line_color, dot_line_width=dot_line_width,
            dot_line_style=dot_line_style, vline_color=vline_color,
            vline_style=vline_style, label_fontsize=label_fontsize,
            label_fontweight=label_fontweight,
            pathway_fontsize=pathway_fontsize,
            pathway_fontstyle=pathway_fontstyle,
            figsize=figsize, filename=filename, dpi=dpi, **kwargs
        )


def _ggcluster_enrich(data: pd.DataFrame,
                      color_low: str = 'pink',
                      color_high: str = 'red',
                      size_range: tuple = (20, 200),
                      curve_color: str = 'grey',
                      curve_linewidth: float = 0.5,
                      vertical_line_color: str = 'darkcyan',
                      vertical_line_width: float = 1.5,
                      dot_line_color: str = 'grey',
                      dot_line_width: float = 0.3,
                      dot_line_style: str = 'dotted',
                      vline_color: str = 'grey',
                      vline_style: str = 'dashed',
                      label_fontsize: int = 8,
                      label_fontweight: str = 'bold',
                      pathway_fontsize: int = 7,
                      pathway_fontstyle: str = 'italic',
                      x_label_fontsize: int = 8,
                      figsize: tuple = (12, 10),
                      filename: Optional[str] = None,
                      dpi: int = 300,
                      **kwargs):
    """
    Internal function for enrichment-mode cluster visualization
    """
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D
    import matplotlib.patches as mpatches

    # Validate required columns - flexible to handle richKEGG output
    required_cols = ['Term', 'Level2', 'group', 'Padj']
    missing = [c for c in required_cols if c not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = data.copy()
    df['Padj'] = pd.to_numeric(df['Padj'], errors='coerce')
    df['neg_log10_Padj'] = -np.log10(df['Padj'].clip(lower=1e-300))

    # Handle RichFactor - use pre-computed if available, otherwise calculate
    if 'RichFactor' in df.columns:
        df['RichFactor'] = pd.to_numeric(df['RichFactor'], errors='coerce')
    elif 'Significant' in df.columns and 'Annotated' in df.columns:
        df['RichFactor'] = df['Significant'] / df['Annotated']
    elif 'Count' in df.columns and 'BgRatio' in df.columns:
        # Extract denominator from BgRatio (e.g., "50/5000" -> 5000)
        def parse_bg_total(ratio):
            try:
                parts = str(ratio).split('/')
                return int(parts[0]) if len(parts) >= 1 else 1
            except:
                return 1
        bg_counts = df['BgRatio'].apply(parse_bg_total)
        df['RichFactor'] = df['Count'] / bg_counts.clip(lower=1)
    else:
        # Default to using Count if nothing else available
        if 'Count' in df.columns:
            df['RichFactor'] = df['Count'] / df['Count'].max()
        else:
            df['RichFactor'] = 1.0
            logger.warning("Could not compute RichFactor, using default value of 1.0")

    # Create term ordering by Level2
    terms = df[['Term', 'Level2']].drop_duplicates()
    terms = terms.sort_values(['Level2', 'Term'])
    terms['y'] = range(len(terms))
    term_to_y = dict(zip(terms['Term'], terms['y']))

    # Get Level2 label positions (mean y of terms in that Level2)
    level2_labels = terms.groupby('Level2')['y'].agg(['mean', 'min', 'max']).reset_index()
    level2_labels.columns = ['Level2', 'y_level2', 'y_min', 'y_max']

    # Merge positions into data
    df['y'] = df['Term'].map(term_to_y)
    df = df.merge(level2_labels[['Level2', 'y_level2']], on='Level2', how='left')

    # Get groups and their x positions
    groups = sorted(df['group'].unique())
    n_groups = len(groups)
    group_x_positions = {g: 2 + i for i, g in enumerate(groups)}
    df['x_group'] = df['group'].map(group_x_positions)

    x_pathway = 2 + n_groups + 0.3

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    y_min_global = terms['y'].min()
    y_max_global = terms['y'].max()

    # Draw curves connecting Level2 labels to their terms
    from matplotlib.patches import FancyArrowPatch

    for _, level2_row in level2_labels.iterrows():
        level2 = level2_row['Level2']
        y_level2 = level2_row['y_level2']
        level2_terms = terms[terms['Level2'] == level2].copy()
        level2_terms = level2_terms.sort_values('y')
        n_terms = len(level2_terms)

        for idx, (_, term_row) in enumerate(level2_terms.iterrows()):
            y_term = term_row['y']
            # Draw curve from (1, y_level2) to (1.5, y_term)
            if n_terms == 1:
                # Single term: straight line
                ax.plot([1, 1.5], [y_level2, y_term],
                       color=curve_color, linewidth=curve_linewidth)
            else:
                # Multiple terms: use FancyArrowPatch for curved connection
                # Determine curvature direction based on position relative to center
                if idx < n_terms / 2:
                    # Upper half: curve bends upward (rad > 0)
                    connectionstyle = f"arc3,rad=0.15"
                else:
                    # Lower half: curve bends downward (rad < 0)
                    connectionstyle = f"arc3,rad=-0.15"

                arrow = FancyArrowPatch(
                    (1, y_level2), (1.5, y_term),
                    connectionstyle=connectionstyle,
                    arrowstyle='-',
                    color=curve_color,
                    linewidth=curve_linewidth,
                    mutation_scale=1
                )
                ax.add_patch(arrow)

    # Draw vertical lines for each Level2 group
    for _, level2_row in level2_labels.iterrows():
        y_min = level2_row['y_min']
        y_max = level2_row['y_max']
        delta = 0.1 if y_min == y_max else 0
        ax.plot([1.5, 1.5], [y_min - delta, y_max + delta],
               color=vertical_line_color, linewidth=vertical_line_width)

    # Draw horizontal dotted lines from cluster dots to pathway labels
    for _, row in df.drop_duplicates(['Term', 'y']).iterrows():
        y = row['y']
        x_start = max(group_x_positions.values())
        ax.plot([x_start, x_pathway], [y, y],
               color=dot_line_color, linewidth=dot_line_width,
               linestyle=dot_line_style)

    # Draw vertical dashed lines separating clusters
    for x in group_x_positions.values():
        ax.axvline(x=x, ymin=0, ymax=1, color=vline_color,
                  linestyle=vline_style, linewidth=0.5, alpha=0.5)

    # Add Level2 labels
    for _, level2_row in level2_labels.iterrows():
        ax.text(0.9, level2_row['y_level2'], level2_row['Level2'],
               ha='right', va='center', fontsize=label_fontsize,
               fontweight=label_fontweight)

    # Plot dots with color by -log10(Padj) and size by RichFactor
    color_values = df['neg_log10_Padj'].values
    size_values = df['RichFactor'].values

    # Handle edge case where all values are identical
    if size_values.min() == size_values.max():
        sizes = np.full_like(size_values, np.mean(size_range))
    else:
        sizes = np.interp(size_values, (size_values.min(), size_values.max()), size_range)

    cmap = LinearSegmentedColormap.from_list('custom', [color_low, color_high])

    if color_values.min() == color_values.max():
        norm = plt.Normalize(vmin=color_values.min(), vmax=color_values.min() + 1)
    else:
        norm = plt.Normalize(vmin=color_values.min(), vmax=color_values.max())

    scatter = ax.scatter(df['x_group'], df['y'],
                        c=color_values, s=sizes,
                        cmap=cmap, norm=norm,
                        edgecolors='black', linewidth=0.5, alpha=0.8)

    # Add pathway labels on the right
    for _, row in terms.iterrows():
        ax.text(x_pathway + 0.1, row['y'], row['Term'],
               ha='left', va='center', fontsize=pathway_fontsize,
               fontstyle=pathway_fontstyle)

    # Set axis properties
    ax.set_xlim(0, x_pathway + 2)
    ax.set_ylim(y_min_global - 1, y_max_global + 1)

    # X-axis: cluster labels
    ax.set_xticks(list(group_x_positions.values()))
    ax.set_xticklabels(groups, rotation=90, ha='center', fontsize=x_label_fontsize)

    # Remove y-axis
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.02, pad=0.15)
    cbar.set_label('-log10(Padj)', fontsize=10)

    # Add size legend
    size_legend_values = [size_values.min(), np.median(size_values), size_values.max()]
    size_legend_sizes = np.interp(size_legend_values,
                                  (size_values.min(), size_values.max()), size_range)

    legend_elements = []
    for val, size in zip(size_legend_values, size_legend_sizes):
        legend_elements.append(
            plt.scatter([], [], s=size, c='gray', alpha=0.6,
                       edgecolors='black', linewidth=0.5,
                       label=f'{val:.2f}')
        )

    ax.legend(handles=legend_elements, scatterpoints=1, frameon=True,
             labelspacing=1.5, title='RichFactor', loc='upper right',
             fontsize=8, bbox_to_anchor=(1.15, 1))

    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {filename}")

    return fig


def _ggcluster_gsea(data: pd.DataFrame,
                    color_low: str = 'blue',
                    color_high: str = 'red',
                    color_mid: str = 'white',
                    size_range: tuple = (20, 200),
                    curve_color: str = 'grey',
                    curve_linewidth: float = 0.5,
                    vertical_line_color: str = 'darkcyan',
                    vertical_line_width: float = 1.5,
                    dot_line_color: str = 'grey',
                    dot_line_width: float = 0.3,
                    dot_line_style: str = 'dotted',
                    vline_color: str = 'grey',
                    vline_style: str = 'dashed',
                    label_fontsize: int = 8,
                    label_fontweight: str = 'bold',
                    pathway_fontsize: int = 7,
                    pathway_fontstyle: str = 'italic',
                    figsize: tuple = (12, 10),
                    filename: Optional[str] = None,
                    dpi: int = 300,
                    **kwargs):
    """
    Internal function for GSEA-mode cluster visualization

    Uses NES (normalized enrichment score) for coloring with diverging colormap
    """
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    from matplotlib.lines import Line2D

    df = data.copy()

    # Normalize column names - handle both richGSEA (Term, Padj) and custom (pathway, padj)
    if 'Term' in df.columns and 'pathway' not in df.columns:
        df['pathway'] = df['Term']
    if 'Padj' in df.columns and 'padj' not in df.columns:
        df['padj'] = df['Padj']

    # Validate required columns
    required_cols = ['pathway', 'Level2', 'group', 'padj', 'NES']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
    df['NES'] = pd.to_numeric(df['NES'], errors='coerce')
    df['neg_log10_padj'] = -np.log10(df['padj'].clip(lower=1e-300))

    # Create term ordering by Level2
    terms = df[['pathway', 'Level2']].drop_duplicates()
    terms = terms.sort_values(['Level2', 'pathway'])
    terms['y'] = range(len(terms))
    term_to_y = dict(zip(terms['pathway'], terms['y']))

    # Get Level2 label positions
    level2_labels = terms.groupby('Level2')['y'].agg(['mean', 'min', 'max']).reset_index()
    level2_labels.columns = ['Level2', 'y_level2', 'y_min', 'y_max']

    # Merge positions into data
    df['y'] = df['pathway'].map(term_to_y)
    df = df.merge(level2_labels[['Level2', 'y_level2']], on='Level2', how='left')

    # Get groups and their x positions
    groups = sorted(df['group'].unique())
    n_groups = len(groups)
    group_x_positions = {g: 2 + i for i, g in enumerate(groups)}
    df['x_group'] = df['group'].map(group_x_positions)

    x_pathway = 2 + n_groups + 0.3

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    y_min_global = terms['y'].min()
    y_max_global = terms['y'].max()

    # Draw curves connecting Level2 labels to their terms
    from matplotlib.patches import FancyArrowPatch

    for _, level2_row in level2_labels.iterrows():
        level2 = level2_row['Level2']
        y_level2 = level2_row['y_level2']
        level2_terms = terms[terms['Level2'] == level2].copy()
        level2_terms = level2_terms.sort_values('y')
        n_terms = len(level2_terms)

        for idx, (_, term_row) in enumerate(level2_terms.iterrows()):
            y_term = term_row['y']
            if n_terms == 1:
                ax.plot([1, 1.5], [y_level2, y_term],
                       color=curve_color, linewidth=curve_linewidth)
            else:
                # Determine curvature direction based on position relative to center
                if idx < n_terms / 2:
                    connectionstyle = f"arc3,rad=0.15"
                else:
                    connectionstyle = f"arc3,rad=-0.15"

                arrow = FancyArrowPatch(
                    (1, y_level2), (1.5, y_term),
                    connectionstyle=connectionstyle,
                    arrowstyle='-',
                    color=curve_color,
                    linewidth=curve_linewidth,
                    mutation_scale=1
                )
                ax.add_patch(arrow)

    # Draw vertical lines for each Level2 group
    for _, level2_row in level2_labels.iterrows():
        y_min = level2_row['y_min']
        y_max = level2_row['y_max']
        delta = 0.1 if y_min == y_max else 0
        ax.plot([1.5, 1.5], [y_min - delta, y_max + delta],
               color=vertical_line_color, linewidth=vertical_line_width)

    # Draw horizontal dotted lines
    for _, row in df.drop_duplicates(['pathway', 'y']).iterrows():
        y = row['y']
        x_start = max(group_x_positions.values())
        ax.plot([x_start, x_pathway], [y, y],
               color=dot_line_color, linewidth=dot_line_width,
               linestyle=dot_line_style)

    # Draw vertical dashed lines separating clusters
    for x in group_x_positions.values():
        ax.axvline(x=x, ymin=0, ymax=1, color=vline_color,
                  linestyle=vline_style, linewidth=0.5, alpha=0.5)

    # Add Level2 labels
    for _, level2_row in level2_labels.iterrows():
        ax.text(0.9, level2_row['y_level2'], level2_row['Level2'],
               ha='right', va='center', fontsize=label_fontsize,
               fontweight=label_fontweight)

    # Plot dots with diverging color by NES and size by -log10(padj)
    color_values = df['NES'].values
    size_values = df['neg_log10_padj'].values

    if size_values.min() == size_values.max():
        sizes = np.full_like(size_values, np.mean(size_range))
    else:
        sizes = np.interp(size_values, (size_values.min(), size_values.max()), size_range)

    # Diverging colormap for NES (negative = blue, positive = red)
    cmap = LinearSegmentedColormap.from_list('custom_diverging',
                                              [color_low, color_mid, color_high])

    # Center the colormap at 0
    vmin, vmax = color_values.min(), color_values.max()
    if vmin >= 0:
        norm = plt.Normalize(vmin=0, vmax=vmax)
    elif vmax <= 0:
        norm = plt.Normalize(vmin=vmin, vmax=0)
    else:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    scatter = ax.scatter(df['x_group'], df['y'],
                        c=color_values, s=sizes,
                        cmap=cmap, norm=norm,
                        edgecolors='black', linewidth=0.5, alpha=0.8)

    # Add pathway labels on the right
    for _, row in terms.iterrows():
        ax.text(x_pathway + 0.1, row['y'], row['pathway'],
               ha='left', va='center', fontsize=pathway_fontsize,
               fontstyle=pathway_fontstyle)

    # Set axis properties
    ax.set_xlim(0, x_pathway + 2)
    ax.set_ylim(y_min_global - 1, y_max_global + 1)

    ax.set_xticks(list(group_x_positions.values()))
    ax.set_xticklabels(groups, rotation=90, ha='center', fontsize=8)

    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Add colorbar for NES
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.02, pad=0.15)
    cbar.set_label('NES', fontsize=10)

    # Add size legend for -log10(padj)
    size_legend_values = [size_values.min(), np.median(size_values), size_values.max()]
    size_legend_sizes = np.interp(size_legend_values,
                                  (size_values.min(), size_values.max()), size_range)

    legend_elements = []
    for val, size in zip(size_legend_values, size_legend_sizes):
        legend_elements.append(
            plt.scatter([], [], s=size, c='gray', alpha=0.6,
                       edgecolors='black', linewidth=0.5,
                       label=f'{val:.2f}')
        )

    ax.legend(handles=legend_elements, scatterpoints=1, frameon=True,
             labelspacing=1.5, title='-log10(padj)', loc='upper right',
             fontsize=8, bbox_to_anchor=(1.15, 1))

    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {filename}")

    return fig


def comparedot(compare_df: pd.DataFrame,
               pvalue: float = 0.05,
               top: int = 10,
               by_sample: bool = True,
               figsize: tuple = (14, 8),
               fontsize: int = 10,
               filename: Optional[str] = None,
               dpi: int = 300,
               **kwargs):
    """
    Create comparative dot plot across samples using RichFactor
    
    Parameters:
    -----------
    compare_df : pd.DataFrame
        DataFrame from compareResult()
    by_sample : bool
        If True, create faceted plot by sample
        If False, overlay all samples
    """
    df = compare_df[compare_df['Padj'] <= pvalue].copy()
    
    if len(df) == 0:
        raise ValueError("No terms remain after filtering")
    
    term_scores = df.groupby('Term')['RichFactor'].mean()
    top_terms = term_scores.nlargest(top).index.tolist()
    
    df = df[df['Term'].isin(top_terms)]
    df['NegLogPadj'] = -np.log10(df['Padj'])
    
    if by_sample:
        samples = df['Sample'].unique()
        n_samples = len(samples)
        
        fig, axes = plt.subplots(1, n_samples, figsize=figsize, sharey=True)
        if n_samples == 1:
            axes = [axes]
        
        for idx, (sample, ax) in enumerate(zip(samples, axes)):
            sample_df = df[df['Sample'] == sample].copy()
            sample_df = sample_df.sort_values('RichFactor', ascending=True)
            
            scatter = ax.scatter(sample_df['RichFactor'],
                               range(len(sample_df)),
                               s=sample_df['Count'] * 20,
                               c=sample_df['NegLogPadj'],
                               cmap='RdYlBu_r',
                               alpha=0.8,
                               edgecolors='black',
                               linewidth=0.5)
            
            if idx == 0:
                ax.set_yticks(range(len(sample_df)))
                ax.set_yticklabels(sample_df['Term'], fontsize=fontsize)
            
            ax.set_xlabel('Rich Factor', fontsize=fontsize+1, fontweight='bold')
            ax.set_title(sample, fontsize=fontsize+2, fontweight='bold')
            ax.axvline(x=1.0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
            ax.grid(True, alpha=0.3, axis='x')
        
        cbar = fig.colorbar(scatter, ax=axes, fraction=0.02, pad=0.04)
        cbar.set_label('-log10(Adjusted P-value)', fontsize=fontsize)
        
    else:
        fig, ax = plt.subplots(figsize=figsize)
        
        samples = df['Sample'].unique()
        colors_map = plt.cm.get_cmap('Set2', len(samples))
        
        term_order = df.groupby('Term')['RichFactor'].mean().sort_values(ascending=True)
        term_positions = {term: i for i, term in enumerate(term_order.index)}
        
        for idx, sample in enumerate(samples):
            sample_df = df[df['Sample'] == sample].copy()
            y_positions = [term_positions[term] for term in sample_df['Term']]
            
            ax.scatter(sample_df['RichFactor'],
                      y_positions,
                      s=sample_df['Count'] * 20,
                      c=[colors_map(idx)] * len(sample_df),
                      alpha=0.7,
                      edgecolors='black',
                      linewidth=0.5,
                      label=sample)
        
        ax.set_yticks(range(len(term_order)))
        ax.set_yticklabels(term_order.index, fontsize=fontsize)
        ax.set_xlabel('Rich Factor', fontsize=fontsize+2, fontweight='bold')
        ax.axvline(x=1.0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left',
                 fontsize=fontsize)
        ax.grid(True, alpha=0.3, axis='x')
    
    plt.suptitle('Comparative Enrichment Analysis',
                fontsize=fontsize+4, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Plot saved to {filename}")
    
    return fig

