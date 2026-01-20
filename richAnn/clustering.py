"""
Clustering functions for richAnn package
"""

import numpy as np
import pandas as pd
from typing import Set
import logging

from .core import EnrichResult

logger = logging.getLogger(__name__)


def richCluster(enrich_result: EnrichResult,
                cutoff: float = 0.5,
                overlap: float = 0.5,
                minSize: int = 5,
                escore: float = 3) -> EnrichResult:
    """
    Cluster enrichment results using kappa statistics
    
    Parameters:
    -----------
    enrich_result : EnrichResult
        Enrichment result object
    cutoff : float
        Kappa score threshold (0-1) for term similarity
    overlap : float
        Overlap threshold for merging clusters (0-1)
    minSize : int
        Minimum cluster size
    escore : float
        Enrichment score cutoff (-log10(Padj))
        
    Returns:
    --------
    EnrichResult with cluster assignments in 'Cluster' column
    """
    df = enrich_result.result.copy()
    
    # Filter by enrichment score
    df['Escore'] = -np.log10(df['Padj'])
    df = df[df['Escore'] >= escore]
    
    if len(df) < 2:
        df['Cluster'] = 1
        logger.warning("Fewer than 2 terms passed escore filter, all assigned to cluster 1")
        return EnrichResult(df, enrich_result.enrichment_type, enrich_result.parameters)
    
    logger.info(f"Clustering {len(df)} terms with kappa >= {cutoff}")
    
    # Build gene sets
    term_genesets = {}
    for idx, row in df.iterrows():
        genes = set(str(row['GeneID']).replace(',', ';').split(';'))
        genes = {g.strip() for g in genes if g.strip()}
        term_genesets[idx] = genes
    
    # Calculate kappa matrix
    indices = list(term_genesets.keys())
    n_terms = len(indices)
    kappa_matrix = np.zeros((n_terms, n_terms))
    
    for i in range(n_terms):
        for j in range(i+1, n_terms):
            idx_i, idx_j = indices[i], indices[j]
            kappa = _calculate_kappa(term_genesets[idx_i], term_genesets[idx_j])
            kappa_matrix[i, j] = kappa
            kappa_matrix[j, i] = kappa
    
    # Initial clustering
    clusters = []
    assigned = set()
    
    for i in range(n_terms):
        if i in assigned:
            continue
        
        cluster = {i}
        assigned.add(i)
        
        for j in range(n_terms):
            if j != i and j not in assigned and kappa_matrix[i, j] >= cutoff:
                cluster.add(j)
                assigned.add(j)
        
        clusters.append(cluster)
    
    # Merge overlapping clusters
    merged = True
    iteration = 0
    while merged and iteration < 100:
        merged = False
        new_clusters = []
        used = set()
        
        for i, cluster_i in enumerate(clusters):
            if i in used:
                continue
            
            for j, cluster_j in enumerate(clusters[i+1:], i+1):
                if j in used:
                    continue
                
                overlap_size = len(cluster_i.intersection(cluster_j))
                min_size = min(len(cluster_i), len(cluster_j))
                
                if min_size > 0 and overlap_size / min_size >= overlap:
                    cluster_i = cluster_i.union(cluster_j)
                    used.add(j)
                    merged = True
            
            new_clusters.append(cluster_i)
            used.add(i)
        
        clusters = new_clusters
        iteration += 1
    
    # Filter by minimum size
    clusters = [c for c in clusters if len(c) >= minSize]
    
    # Assign cluster IDs
    cluster_assignment = {}
    for cluster_id, cluster in enumerate(clusters, 1):
        for idx_pos in cluster:
            idx = indices[idx_pos]
            cluster_assignment[idx] = cluster_id
    
    df['Cluster'] = df.index.map(cluster_assignment)
    df = df[df['Cluster'].notna()]
    
    logger.info(f"Clustering complete: {len(clusters)} clusters, {len(df)} terms retained")
    
    return EnrichResult(df, enrich_result.enrichment_type, enrich_result.parameters)


def _calculate_kappa(genes_a: Set[str], genes_b: Set[str]) -> float:
    """Calculate Cohen's kappa statistic between two gene sets"""
    all_genes = genes_a.union(genes_b)
    N = len(all_genes)
    
    if N == 0:
        return 0.0
    
    overlap = len(genes_a.intersection(genes_b))
    
    # Observed agreement
    p_o = (overlap + (N - len(genes_a) - len(genes_b) + overlap)) / N
    
    # Expected agreement
    p_a = (len(genes_a) / N) * (len(genes_b) / N)
    p_b = ((N - len(genes_a)) / N) * ((N - len(genes_b)) / N)
    p_e = p_a + p_b
    
    if p_e >= 1:
        return 0.0
    
    kappa = (p_o - p_e) / (1 - p_e)
    
    return max(0.0, min(1.0, kappa))

