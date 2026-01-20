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

