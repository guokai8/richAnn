"""Utility functions for richAnn package"""

import numpy as np
import pandas as pd
from typing import List, Set, Union
import logging
from .core import EnrichResult

logger = logging.getLogger(__name__)


def _normalize_genes(genes: Union[List[str], Set[str], np.ndarray], 
                     case_sensitive: bool = False) -> Set[str]:
    """Normalize gene list: convert to set, handle case sensitivity"""
    if isinstance(genes, np.ndarray):
        genes = genes.tolist()
    elif not isinstance(genes, (list, set)):
        raise TypeError("Genes must be a list, set, or numpy array")
    
    gene_set = {str(g).strip() for g in genes if str(g).strip()}
    
    if not case_sensitive:
        gene_set = {g.upper() for g in gene_set}
    
    return gene_set


def _validate_annotation(annot_df: pd.DataFrame, 
                        required_cols: List[str]) -> None:
    """Validate annotation dataframe structure"""
    if not isinstance(annot_df, pd.DataFrame):
        raise TypeError("Annotation must be a pandas DataFrame")
    
    missing_cols = set(required_cols) - set(annot_df.columns)
    if missing_cols:
        raise ValueError(f"Annotation missing required columns: {missing_cols}")
    
    if len(annot_df) == 0:
        raise ValueError("Annotation dataframe is empty")


def _calculate_effect_sizes(k: int, n: int, M: int, N: int) -> dict:
    """
    Calculate multiple effect size metrics for enrichment

    Parameters:
    -----------
    k : int - Number of significant genes (overlap)
    n : int - Total query genes
    M : int - Total genes in term/pathway (annotated)
    N : int - Total background genes

    Returns:
    --------
    dict with effect size metrics
    """
    epsilon = 1e-10

    # RichFactor = k/M (proportion of pathway genes that are significant)
    rich_factor = k / max(M, 1)

    # FoldEnrichment = (k/n) / (M/N) (enrichment ratio)
    gene_ratio = k / max(n, 1)
    bg_ratio = M / max(N, 1)
    fold_enrichment = gene_ratio / max(bg_ratio, epsilon)

    # Odds ratio
    numerator = k * (N - M)
    denominator = max((n - k) * M, epsilon)
    odds_ratio = numerator / denominator

    log_odds_ratio = np.log2(max(odds_ratio, epsilon))

    # Calculate z-score for enrichment
    # Expected count under null hypothesis
    expected = n * M / max(N, 1)
    # Variance under hypergeometric distribution
    variance = n * (M / max(N, 1)) * ((N - M) / max(N, 1)) * ((N - n) / max(N - 1, 1))
    # Z-score
    zscore = (k - expected) / max(np.sqrt(variance), epsilon)

    return {
        'RichFactor': rich_factor,
        'FoldEnrichment': fold_enrichment,
        'FoldEnrich': fold_enrichment,  # Keep for backward compatibility
        'OddsRatio': odds_ratio,
        'LogOddsRatio': log_odds_ratio,
        'zscore': zscore
    }


def compareResult(result_dict: dict) -> pd.DataFrame:
    """Compare enrichment results across multiple samples/conditions"""
    if not result_dict:
        raise ValueError("result_dict cannot be empty")

    combined = []
    for sample_name, result in result_dict.items():
        if not isinstance(result, EnrichResult):
            raise TypeError(f"All values must be EnrichResult objects")

        df = result.result.copy()
        df['Sample'] = sample_name
        combined.append(df)

    return pd.concat(combined, ignore_index=True)


def _validate_pvalue(pvalue: float, param_name: str = "pvalue") -> None:
    """Validate p-value parameter is in valid range"""
    if not isinstance(pvalue, (int, float)):
        raise TypeError(f"{param_name} must be a number, got {type(pvalue)}")
    if not 0 < pvalue <= 1:
        raise ValueError(f"{param_name} must be between 0 and 1, got {pvalue}")


def _validate_size_params(min_size: int, max_size: int) -> None:
    """Validate min/max size parameters"""
    if not isinstance(min_size, int) or not isinstance(max_size, int):
        raise TypeError("minSize and maxSize must be integers")
    if min_size < 0:
        raise ValueError(f"minSize must be non-negative, got {min_size}")
    if max_size < min_size:
        raise ValueError(f"maxSize ({max_size}) must be >= minSize ({min_size})")


def _validate_positive_int(value: int, param_name: str) -> None:
    """Validate that a parameter is a positive integer"""
    if not isinstance(value, int):
        raise TypeError(f"{param_name} must be an integer, got {type(value)}")
    if value <= 0:
        raise ValueError(f"{param_name} must be positive, got {value}")

