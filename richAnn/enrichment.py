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


def _convert_padj_method(method: str) -> str:
    """
    Convert R-style method names to statsmodels format

    Parameters:
    -----------
    method : str
        Method name in R format (e.g., 'BH', 'BY', 'bonferroni')
        or statsmodels format (e.g., 'fdr_bh', 'fdr_by', 'bonferroni')

    Returns:
    --------
    str: Method name in statsmodels format
    """
    # Mapping from R-style to statsmodels-style
    method_map = {
        'BH': 'fdr_bh',
        'BY': 'fdr_by',
        'bonferroni': 'bonferroni',
        'holm': 'holm',
        'hommel': 'hommel',
        'hochberg': 'simes-hochberg',
        'fdr': 'fdr_bh',  # Default FDR is Benjamini-Hochberg
    }

    # If already in statsmodels format, return as-is
    if method in ['fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky', 'bonferroni',
                  'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel']:
        return method

    # Convert from R-style
    if method.upper() in method_map:
        return method_map[method.upper()]

    # If not recognized, return as-is and let multipletests raise the error
    return method


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

    # Convert method name to statsmodels format
    method_converted = _convert_padj_method(padj_method)
    results_df['Padj'] = multipletests(results_df['Pvalue'], method=method_converted)[1]
    
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

    # Convert method name to statsmodels format
    method_converted = _convert_padj_method(padj_method)
    results_df['Padj'] = multipletests(results_df['Pvalue'], method=method_converted)[1]
    
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
             permutation_num: int = 1000,
             min_size: int = 15,
             max_size: int = 500,
             weight: float = 1.0,
             seed: int = 123,
             case_sensitive: bool = False,
             padj_method: str = "BH") -> EnrichResult:
    """
    Perform Gene Set Enrichment Analysis (GSEA) - matches GSEApy implementation

    Parameters:
    -----------
    gene_scores : dict
        Gene IDs mapped to scores (higher = more important). Genes will be ranked
        by these scores in descending order.
    geneset_db : pd.DataFrame
        Gene sets with columns: GeneSet, GeneSetName, GeneID
    permutation_num : int
        Number of permutations for statistical testing (default: 1000)
    min_size : int
        Minimum gene set size (default: 15)
    max_size : int
        Maximum gene set size (default: 500)
    weight : float
        Weighting exponent used in the enrichment score calculation (default: 1.0)
        - weight = 0: standard Kolmogorov-Smirnov statistic
        - weight = 1: weighted by correlation (GSEApy default)
        - weight = 2: over-weighted by correlation
    seed : int
        Random seed for permutation testing (default: 123)
    case_sensitive : bool
        Case-sensitive gene matching (default: False)
    padj_method : str
        Multiple testing correction method (default: "BH" for Benjamini-Hochberg)

    Returns:
    --------
    EnrichResult object with columns:
        - Annot: Gene set ID
        - Term: Gene set name
        - ES: Enrichment Score
        - NES: Normalized Enrichment Score
        - Pvalue: Nominal p-value
        - Padj: Adjusted p-value (FDR)
        - LeadingEdge: Genes in the leading edge
        - Count: Total genes in gene set
        - Significant: Number of leading edge genes

    Notes:
    ------
    This implementation matches GSEApy's prerank algorithm:
    1. Genes are ranked by input scores
    2. Running enrichment score is calculated with weighted scoring
    3. Statistical significance determined by permutation testing
    4. NES computed by normalizing against permutation distribution
    """
    # Validate inputs
    _validate_annotation(geneset_db, ['GeneSet', 'GeneSetName', 'GeneID'])
    _validate_positive_int(permutation_num, "permutation_num")
    _validate_size_params(min_size, max_size)

    if not gene_scores:
        raise ValueError("gene_scores dictionary cannot be empty")
    if not all(isinstance(v, (int, float)) for v in gene_scores.values()):
        raise TypeError("All gene_scores values must be numbers")
    if not isinstance(weight, (int, float)) or weight < 0:
        raise ValueError(f"weight must be a non-negative number, got {weight}")

    # Set random seed for reproducibility
    np.random.seed(seed)

    # Handle case sensitivity
    if not case_sensitive:
        gene_scores = {str(k).upper(): v for k, v in gene_scores.items()}
        geneset_db = geneset_db.copy()
        geneset_db['GeneID'] = geneset_db['GeneID'].astype(str).str.upper()

    # Rank genes by score (descending order - higher scores first)
    ranked_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)
    ranked_gene_list = [g[0] for g in ranked_genes]
    ranked_scores = np.array([g[1] for g in ranked_genes])

    logger.info(f"GSEA: {len(ranked_gene_list)} ranked genes, weight={weight}, permutations={permutation_num}")

    # Group gene sets
    geneset_groups = geneset_db.groupby(['GeneSet', 'GeneSetName'])

    results = []
    for (geneset_id, geneset_name), group in geneset_groups:
        geneset_genes = set(group['GeneID'].astype(str).unique())
        geneset_size = len(geneset_genes)

        # Apply size filters
        if geneset_size < min_size or geneset_size > max_size:
            continue

        # Calculate GSEA score
        es, nes, pval, leading_edge = _calculate_gsea_score(
            ranked_gene_list, ranked_scores, geneset_genes, permutation_num, weight
        )

        results.append({
            'Annot': geneset_id,
            'Term': geneset_name,
            'ES': es,
            'NES': nes,
            'Pvalue': pval,
            'LeadingEdge': ';'.join(leading_edge),
            'Count': geneset_size,
            'Significant': len(leading_edge)
        })

    if len(results) == 0:
        raise ValueError(
            f"No gene sets passed size filters (min_size={min_size}, max_size={max_size}). "
            f"Check that your gene sets contain genes present in the ranked list."
        )

    results_df = pd.DataFrame(results)

    # Apply multiple testing correction
    method_converted = _convert_padj_method(padj_method)
    results_df['Padj'] = multipletests(results_df['Pvalue'], method=method_converted)[1]

    # Sort by nominal p-value
    results_df = results_df.sort_values('Pvalue').reset_index(drop=True)

    parameters = {
        'n_genes': len(ranked_gene_list),
        'permutation_num': permutation_num,
        'min_size': min_size,
        'max_size': max_size,
        'weight': weight,
        'seed': seed,
        'padj_method': padj_method
    }

    logger.info(f"GSEA complete: {len(results_df)} gene sets tested")

    return EnrichResult(results_df, enrichment_type="GSEA", parameters=parameters)


def _calculate_gsea_score(ranked_genes: List[str],
                          ranked_scores: np.ndarray,
                          geneset: Set[str],
                          permutation_num: int,
                          weight: float = 1.0) -> Tuple[float, float, float, List[str]]:
    """
    Calculate GSEA enrichment score with permutation testing - matches GSEApy algorithm

    This implements the weighted Kolmogorov-Smirnov-like running sum statistic
    used in GSEA (Subramanian et al., PNAS 2005).

    Parameters:
    -----------
    ranked_genes : list
        Gene identifiers in ranked order (highest to lowest score)
    ranked_scores : np.ndarray
        Corresponding scores for ranked genes
    geneset : set
        Set of gene identifiers in the gene set to test
    permutation_num : int
        Number of permutations for p-value calculation
    weight : float
        Weighting exponent (default: 1.0)
        - weight = 0: standard Kolmogorov-Smirnov
        - weight = 1: weighted by correlation (GSEApy default)
        - weight > 1: over-weighted

    Returns:
    --------
    tuple: (ES, NES, pvalue, leading_edge_genes)
        ES: Enrichment Score (maximum deviation from zero)
        NES: Normalized Enrichment Score (ES normalized by permutation distribution)
        pvalue: Nominal p-value from permutation test
        leading_edge_genes: Genes contributing to the enrichment signal

    Algorithm:
    ----------
    1. Calculate running enrichment score:
       - For hits (genes in set): add weighted score
       - For misses: subtract uniform penalty
    2. ES is the maximum deviation from zero
    3. Leading edge: all hit genes up to the peak
    4. Permutation: randomly shuffle hit positions
    5. NES: normalize by mean of permutation distribution
    6. P-value: fraction of permutations with ES >= observed (or <= for negative ES)
    """
    N = len(ranked_genes)

    # Find positions of genes in the gene set
    hit_indices = np.array([i for i, g in enumerate(ranked_genes) if g in geneset])

    if len(hit_indices) == 0:
        return 0.0, 0.0, 1.0, []

    N_hit = len(hit_indices)
    N_miss = N - N_hit

    # Calculate observed running enrichment score
    running_sum = _calculate_running_sum(ranked_scores, hit_indices, N, N_hit, N_miss, weight)

    # Find enrichment score (maximum deviation from zero)
    max_es = np.max(running_sum)
    min_es = np.min(running_sum)
    es = max_es if abs(max_es) > abs(min_es) else min_es

    # Find leading edge genes (all hits up to the peak)
    peak_idx = np.argmax(running_sum) if es > 0 else np.argmin(running_sum)
    leading_edge = [ranked_genes[i] for i in hit_indices if i <= peak_idx]

    # Permutation test to calculate p-value and NES
    perm_scores = np.zeros(permutation_num)

    for perm_idx in range(permutation_num):
        # Randomly shuffle hit positions
        perm_hit_indices = np.random.choice(N, size=N_hit, replace=False)

        # Calculate running sum for permutation
        perm_running = _calculate_running_sum(
            ranked_scores, perm_hit_indices, N, N_hit, N_miss, weight
        )

        # Get permutation ES
        perm_max = np.max(perm_running)
        perm_min = np.min(perm_running)
        perm_scores[perm_idx] = perm_max if abs(perm_max) > abs(perm_min) else perm_min

    # Calculate Normalized Enrichment Score (NES)
    # Normalize by the mean of positive or negative permutation scores
    if es >= 0:
        pos_perms = perm_scores[perm_scores >= 0]
        mean_pos = np.mean(pos_perms) if len(pos_perms) > 0 else 1.0
        nes = es / mean_pos if mean_pos > 0 else es
    else:
        neg_perms = perm_scores[perm_scores < 0]
        mean_neg = np.mean(np.abs(neg_perms)) if len(neg_perms) > 0 else 1.0
        # Preserve the negative sign of ES in NES
        nes = es / mean_neg if mean_neg > 0 else es

    # Calculate nominal p-value
    if es >= 0:
        pval = np.sum(perm_scores >= es) / permutation_num
    else:
        pval = np.sum(perm_scores <= es) / permutation_num

    # Bound p-value to avoid zeros (minimum is 1/permutation_num)
    pval = max(pval, 1.0 / permutation_num)

    return es, nes, pval, leading_edge


def _calculate_running_sum(ranked_scores: np.ndarray,
                           hit_indices: np.ndarray,
                           N: int,
                           N_hit: int,
                           N_miss: int,
                           weight: float) -> np.ndarray:
    """
    Calculate the running sum statistic for GSEA

    Parameters:
    -----------
    ranked_scores : np.ndarray
        Scores for all ranked genes
    hit_indices : np.ndarray
        Indices of genes in the gene set
    N : int
        Total number of genes
    N_hit : int
        Number of genes in the gene set
    N_miss : int
        Number of genes not in the gene set
    weight : float
        Weighting exponent

    Returns:
    --------
    np.ndarray: Cumulative running sum
    """
    running_sum = np.zeros(N)
    hit_set = set(hit_indices)

    if weight == 0:
        # Unweighted (classic Kolmogorov-Smirnov)
        for i in range(N):
            if i in hit_set:
                running_sum[i] = 1.0 / N_hit
            else:
                running_sum[i] = -1.0 / N_miss if N_miss > 0 else 0.0
    else:
        # Weighted by correlation scores
        abs_scores = np.abs(ranked_scores)

        # Apply weight exponent and calculate normalization factor
        weighted_scores = np.power(abs_scores[hit_indices], weight)
        N_R = np.sum(weighted_scores)

        for i in range(N):
            if i in hit_set:
                if N_R > 0:
                    running_sum[i] = np.power(abs_scores[i], weight) / N_R
                else:
                    # Fallback if all scores are zero
                    running_sum[i] = 1.0 / N_hit
            else:
                running_sum[i] = -1.0 / N_miss if N_miss > 0 else 0.0

    # Calculate cumulative sum
    return np.cumsum(running_sum)

