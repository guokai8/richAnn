"""
Enrichment analysis functions for richAnn package
"""

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from typing import List, Dict, Optional, Union, Set, Tuple
import logging

from .core import EnrichResult
from .utils import (_normalize_genes, _validate_annotation, _calculate_effect_sizes,
                    _validate_pvalue, _validate_size_params, _validate_positive_int)

logger = logging.getLogger(__name__)


def richGO(genes: Union[List[str], Set[str], np.ndarray], 
           godata: pd.DataFrame,
           ontology: str = "BP",
           pvalue: float = 0.05,
           padj: Optional[float] = None,
           padj_method: str = "BH",
           minSize: int = 2,
           maxSize: int = 500,
           keepRich: bool = False,
           universe: Optional[Union[List[str], Set[str]]] = None,
           case_sensitive: bool = False,
           sep: str = ";") -> EnrichResult:
    """
    Perform GO enrichment analysis using hypergeometric test
    
    Parameters:
    -----------
    genes : list, set, or array
        Gene identifiers to test for enrichment
    godata : pd.DataFrame
        GO annotation with columns: GeneID, GOterm, GOname, Ontology
    ontology : str
        "BP", "MF", or "CC"
    pvalue : float
        Raw p-value cutoff
    padj : float
        Adjusted p-value cutoff
    padj_method : str
        Multiple testing correction method
    minSize : int
        Minimum gene set size
    maxSize : int
        Maximum gene set size
    keepRich : bool
        Keep terms with RichFactor = 1
    universe : list/set
        Background gene universe
    case_sensitive : bool
        Case-sensitive gene matching
    sep : str
        Separator for gene IDs in results
        
    Returns:
    --------
    EnrichResult object
    """
    # Validate inputs
    _validate_annotation(godata, ['GeneID', 'GOterm', 'GOname', 'Ontology'])
    _validate_pvalue(pvalue, "pvalue")
    if padj is not None:
        _validate_pvalue(padj, "padj")
    _validate_size_params(minSize, maxSize)

    if ontology not in ['BP', 'MF', 'CC']:
        raise ValueError("Ontology must be 'BP', 'MF', or 'CC'")

    try:
        query_genes = _normalize_genes(genes, case_sensitive)
    except (TypeError, AttributeError) as e:
        raise ValueError(f"Invalid gene list format: {e}. Expected list, set, or array of gene IDs.")

    if len(query_genes) == 0:
        raise ValueError("No valid genes provided. Check that your gene list is not empty and contains valid identifiers.")
    
    godata = godata.copy()
    if not case_sensitive:
        godata['GeneID'] = godata['GeneID'].str.upper()
    
    godata = godata[godata['Ontology'] == ontology].copy()
    
    if len(godata) == 0:
        raise ValueError(f"No GO terms found for ontology '{ontology}'")
    
    if universe is not None:
        background_genes = _normalize_genes(universe, case_sensitive)
    else:
        background_genes = set(godata['GeneID'].unique())
    
    N = len(background_genes)
    query_genes = query_genes.intersection(background_genes)
    n = len(query_genes)

    if n == 0:
        original_count = len(_normalize_genes(genes, case_sensitive))
        raise ValueError(
            f"No query genes found in background universe. "
            f"Started with {original_count} genes, but none overlap with the {N} genes in the annotation. "
            f"Check: (1) Gene ID format matches annotation, (2) Case sensitivity setting (case_sensitive={case_sensitive}), "
            f"(3) Annotation covers your organism/genes."
        )

    logger.info(f"GO enrichment: {n} genes, {N} background, ontology={ontology}")
    
    go_groups = godata.groupby(['GOterm', 'GOname'])
    
    results = []
    for (go_id, go_name), group in go_groups:
        term_genes = set(group['GeneID'].unique())
        term_genes = term_genes.intersection(background_genes)
        M = len(term_genes)
        
        if M < minSize or M > maxSize:
            continue
        
        overlap_genes = query_genes.intersection(term_genes)
        k = len(overlap_genes)
        
        if k == 0:
            continue
        
        pval = hypergeom.sf(k - 1, N, M, n)
        effect_sizes = _calculate_effect_sizes(k, n, M, N)
        
        if not keepRich and abs(effect_sizes['RichFactor'] - 1.0) < 1e-10:
            continue
        
        results.append({
            'Annot': go_id,
            'Term': go_name,
            'Pvalue': pval,
            'GeneID': sep.join(sorted(overlap_genes)),
            'Count': k,
            'GeneRatio': f"{k}/{n}",
            'BgRatio': f"{M}/{N}",
            **effect_sizes,
            'Annotated': M,
            'Significant': k
        })
    
    if len(results) == 0:
        raise ValueError(
            f"No enriched terms found for {ontology} ontology. "
            f"Tested {n} query genes against {N} background genes. "
            f"Try: (1) Relaxing pvalue cutoff (current: {pvalue}), "
            f"(2) Adjusting gene set size limits (current: {minSize}-{maxSize}), "
            f"(3) Using a different ontology (BP/MF/CC)."
        )

    results_df = pd.DataFrame(results)
    results_df['Padj'] = multipletests(results_df['Pvalue'], method=padj_method)[1]
    
    if padj is not None:
        results_df = results_df[results_df['Padj'] <= padj]
    else:
        results_df = results_df[results_df['Pvalue'] <= pvalue]
    
    results_df = results_df.sort_values('Padj').reset_index(drop=True)
    
    parameters = {
        'ontology': ontology,
        'n_query_genes': n,
        'n_background_genes': N,
        'pvalue_cutoff': pvalue,
        'padj_cutoff': padj,
        'padj_method': padj_method,
        'minSize': minSize,
        'maxSize': maxSize
    }
    
    logger.info(f"Enrichment complete: {len(results_df)} significant terms")
    
    return EnrichResult(results_df, enrichment_type="GO", parameters=parameters)


def richKEGG(genes: Union[List[str], Set[str], np.ndarray], 
             kodata: pd.DataFrame,
             pvalue: float = 0.05,
             padj: Optional[float] = None,
             padj_method: str = "BH",
             minSize: int = 2,
             maxSize: int = 500,
             keepRich: bool = False,
             universe: Optional[Union[List[str], Set[str]]] = None,
             case_sensitive: bool = False,
             sep: str = ";") -> EnrichResult:
    """
    Perform KEGG pathway enrichment analysis
    
    Parameters same as richGO, except:
    kodata : pd.DataFrame
        KEGG annotation with columns: GeneID, Pathway, PathwayName
    """
    # Validate inputs
    _validate_annotation(kodata, ['GeneID', 'Pathway', 'PathwayName'])
    _validate_pvalue(pvalue, "pvalue")
    if padj is not None:
        _validate_pvalue(padj, "padj")
    _validate_size_params(minSize, maxSize)
    
    query_genes = _normalize_genes(genes, case_sensitive)
    
    if len(query_genes) == 0:
        raise ValueError("No valid genes provided")
    
    kodata = kodata.copy()
    if not case_sensitive:
        kodata['GeneID'] = kodata['GeneID'].str.upper()
    
    if universe is not None:
        background_genes = _normalize_genes(universe, case_sensitive)
    else:
        background_genes = set(kodata['GeneID'].unique())
    
    N = len(background_genes)
    query_genes = query_genes.intersection(background_genes)
    n = len(query_genes)
    
    if n == 0:
        raise ValueError("No query genes found in background universe")
    
    logger.info(f"KEGG enrichment: {n} genes, {N} background")
    
    pathway_groups = kodata.groupby(['Pathway', 'PathwayName'])

    # Check if hierarchy columns exist
    has_hierarchy = all(col in kodata.columns for col in ['Level3', 'Level2', 'Level1'])

    results = []
    for (pathway_id, pathway_name), group in pathway_groups:
        term_genes = set(group['GeneID'].unique())
        term_genes = term_genes.intersection(background_genes)
        M = len(term_genes)

        if M < minSize or M > maxSize:
            continue

        overlap_genes = query_genes.intersection(term_genes)
        k = len(overlap_genes)

        if k == 0:
            continue

        pval = hypergeom.sf(k - 1, N, M, n)
        effect_sizes = _calculate_effect_sizes(k, n, M, N)

        if not keepRich and abs(effect_sizes['RichFactor'] - 1.0) < 1e-10:
            continue

        result_dict = {
            'Annot': pathway_id,
            'Term': pathway_name,
            'Pvalue': pval,
            'GeneID': sep.join(sorted(overlap_genes)),
            'Count': k,
            'GeneRatio': f"{k}/{n}",
            'BgRatio': f"{M}/{N}",
            **effect_sizes,
            'Annotated': M,
            'Significant': k,
            'ko': pathway_id  # Add 'ko' as copy of 'Annot' for KEGG compatibility
        }

        # Add hierarchy columns if they exist in the annotation
        if has_hierarchy:
            result_dict['Level3'] = group['Level3'].iloc[0]
            result_dict['Level2'] = group['Level2'].iloc[0]
            result_dict['Level1'] = group['Level1'].iloc[0]

        results.append(result_dict)
    
    if len(results) == 0:
        raise ValueError("No enriched pathways found")
    
    results_df = pd.DataFrame(results)
    results_df['Padj'] = multipletests(results_df['Pvalue'], method=padj_method)[1]
    
    if padj is not None:
        results_df = results_df[results_df['Padj'] <= padj]
    else:
        results_df = results_df[results_df['Pvalue'] <= pvalue]
    
    results_df = results_df.sort_values('Padj').reset_index(drop=True)
    
    parameters = {
        'n_query_genes': n,
        'n_background_genes': N,
        'pvalue_cutoff': pvalue,
        'padj_cutoff': padj,
        'padj_method': padj_method,
        'minSize': minSize,
        'maxSize': maxSize
    }
    
    logger.info(f"Enrichment complete: {len(results_df)} significant pathways")
    
    return EnrichResult(results_df, enrichment_type="KEGG", parameters=parameters)


def richGSEA(gene_scores: Dict[str, float],
             geneset_db: pd.DataFrame,
             nperm: int = 1000,
             min_size: int = 15,
             max_size: int = 500,
             weighted_score: bool = True,
             random_state: int = 42,
             case_sensitive: bool = False) -> EnrichResult:
    """
    Perform Gene Set Enrichment Analysis (GSEA)
    
    Parameters:
    -----------
    gene_scores : dict
        Gene IDs mapped to scores (higher = more important)
    geneset_db : pd.DataFrame
        Gene sets with columns: GeneSet, GeneSetName, GeneID
    nperm : int
        Number of permutations
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size
    weighted_score : bool
        Use weighted enrichment score
    random_state : int
        Random seed
    case_sensitive : bool
        Case-sensitive matching
    """
    # Validate inputs
    _validate_annotation(geneset_db, ['GeneSet', 'GeneSetName', 'GeneID'])
    _validate_positive_int(nperm, "nperm")
    _validate_size_params(min_size, max_size)

    if not gene_scores:
        raise ValueError("gene_scores dictionary cannot be empty")
    if not all(isinstance(v, (int, float)) for v in gene_scores.values()):
        raise TypeError("All gene_scores values must be numbers")
    
    np.random.seed(random_state)
    
    if not case_sensitive:
        gene_scores = {k.upper(): v for k, v in gene_scores.items()}
        geneset_db = geneset_db.copy()
        geneset_db['GeneID'] = geneset_db['GeneID'].str.upper()
    
    ranked_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)
    ranked_gene_list = [g[0] for g in ranked_genes]
    ranked_scores = np.array([g[1] for g in ranked_genes])
    
    logger.info(f"GSEA: {len(ranked_gene_list)} ranked genes")
    
    geneset_groups = geneset_db.groupby(['GeneSet', 'GeneSetName'])
    
    results = []
    for (geneset_id, geneset_name), group in geneset_groups:
        geneset_genes = set(group['GeneID'].unique())
        
        if len(geneset_genes) < min_size or len(geneset_genes) > max_size:
            continue
        
        es, nes, pval, leading_edge = _calculate_gsea_score(
            ranked_gene_list, ranked_scores, geneset_genes, nperm, weighted_score
        )
        
        results.append({
            'Annot': geneset_id,
            'Term': geneset_name,
            'ES': es,
            'NES': nes,
            'Pvalue': pval,
            'LeadingEdge': ';'.join(leading_edge),
            'Count': len(geneset_genes),
            'Significant': len(leading_edge)
        })
    
    if len(results) == 0:
        raise ValueError("No gene sets passed size filters")
    
    results_df = pd.DataFrame(results)
    results_df['Padj'] = multipletests(results_df['Pvalue'], method='BH')[1]
    results_df = results_df.sort_values('Pvalue').reset_index(drop=True)
    
    parameters = {
        'n_genes': len(ranked_gene_list),
        'nperm': nperm,
        'min_size': min_size,
        'max_size': max_size,
        'weighted_score': weighted_score,
        'random_state': random_state
    }
    
    logger.info(f"GSEA complete: {len(results_df)} gene sets tested")
    
    return EnrichResult(results_df, enrichment_type="GSEA", parameters=parameters)


def _calculate_gsea_score(ranked_genes: List[str], 
                          ranked_scores: np.ndarray,
                          geneset: Set[str], 
                          nperm: int,
                          weighted: bool = True) -> Tuple[float, float, float, List[str]]:
    """
    Calculate GSEA enrichment score with permutation testing
    
    Returns:
    --------
    tuple: (ES, NES, pvalue, leading_edge_genes)
    """
    N = len(ranked_genes)

    hit_indices = np.array([i for i, g in enumerate(ranked_genes) if g in geneset])

    if len(hit_indices) == 0:
        return 0, 0, 1.0, []

    # Convert to set for O(1) membership checking instead of O(n) with arrays
    hit_set = set(hit_indices)

    running_sum = np.zeros(N)

    if weighted:
        abs_scores = np.abs(ranked_scores)
        N_R = np.sum(abs_scores[hit_indices])
        N_miss = N - len(hit_indices)

        for i in range(N):
            if i in hit_set:
                if N_R > 0:
                    running_sum[i] = abs_scores[i] / N_R
                else:
                    running_sum[i] = 1.0 / len(hit_indices)
            else:
                if N_miss > 0:
                    running_sum[i] = -1.0 / N_miss
    else:
        N_hit = len(hit_indices)
        N_miss = N - N_hit

        for i in range(N):
            if i in hit_set:
                running_sum[i] = 1.0 / N_hit
            else:
                running_sum[i] = -1.0 / N_miss
    
    running_sum = np.cumsum(running_sum)
    
    max_es = np.max(running_sum)
    min_es = np.min(running_sum)
    es = max_es if abs(max_es) > abs(min_es) else min_es
    
    peak_idx = np.argmax(running_sum) if es > 0 else np.argmin(running_sum)
    leading_edge = [ranked_genes[i] for i in hit_indices if i <= peak_idx]
    
    # Permutation test
    perm_scores = []
    for _ in range(nperm):
        perm_hit_indices = np.random.choice(N, size=len(hit_indices), replace=False)
        perm_hit_set = set(perm_hit_indices)  # Convert to set for fast lookup

        perm_running = np.zeros(N)

        if weighted:
            abs_scores = np.abs(ranked_scores)
            N_R_perm = np.sum(abs_scores[perm_hit_indices])
            N_miss = N - len(hit_indices)

            for i in range(N):
                if i in perm_hit_set:
                    if N_R_perm > 0:
                        perm_running[i] = abs_scores[i] / N_R_perm
                    else:
                        perm_running[i] = 1.0 / len(perm_hit_indices)
                else:
                    perm_running[i] = -1.0 / N_miss
        else:
            for i in range(N):
                if i in perm_hit_set:
                    perm_running[i] = 1.0 / len(perm_hit_indices)
                else:
                    perm_running[i] = -1.0 / N_miss
        
        perm_running = np.cumsum(perm_running)
        perm_max = np.max(perm_running)
        perm_min = np.min(perm_running)
        perm_es = perm_max if abs(perm_max) > abs(perm_min) else perm_min
        perm_scores.append(perm_es)
    
    perm_scores = np.array(perm_scores)
    
    # Calculate NES
    if es >= 0:
        pos_perms = perm_scores[perm_scores >= 0]
        nes = es / np.mean(pos_perms) if len(pos_perms) > 0 and np.mean(pos_perms) > 0 else es
        pval = np.sum(perm_scores >= es) / nperm
    else:
        neg_perms = perm_scores[perm_scores < 0]
        nes = -es / np.mean(-neg_perms) if len(neg_perms) > 0 and np.mean(-neg_perms) > 0 else es
        pval = np.sum(perm_scores <= es) / nperm
    
    pval = max(pval, 1.0 / nperm)
    
    return es, nes, pval, leading_edge

