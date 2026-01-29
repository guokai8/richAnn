"""
richAnn: Functional Enrichment Analysis and Visualization for Python

A comprehensive Python package for gene set enrichment analysis with 
publication-quality visualizations, inspired by the richR package.
"""

from .core import EnrichResult
from .enrichment import richGO, richKEGG, richGSEA
from .clustering import richCluster
from .utils import compareResult
from .vis import (
    ggbar, ggdot, ggnetplot, ggnetwork, ggnetmap,
    ggheatmap, ggGSEA, comparedot, ggcluster
)
from .data import (
    load_gmt, load_go_gaf, load_kegg_mapping,
    create_custom_annotation, validate_annotation_format,
    from_pathwaydb_go, from_pathwaydb_kegg
)

__version__ = "1.0.0"
__author__ = "richAnn Development Team"
__email__ = "support@richann.org"

__all__ = [
    # Core classes
    'EnrichResult',
    # Enrichment analysis
    'richGO',
    'richKEGG',
    'richGSEA',
    'richCluster',
    # Utilities
    'compareResult',
    # Visualizations
    'ggbar',
    'ggdot',
    'ggnetplot',
    'ggnetwork',
    'ggnetmap',
    'ggheatmap',
    'ggGSEA',
    'comparedot',
    'ggcluster',
    # Data loading
    'load_gmt',
    'load_go_gaf',
    'load_kegg_mapping',
    'create_custom_annotation',
    'validate_annotation_format',
    # pathwaydb integration
    'from_pathwaydb_go',
    'from_pathwaydb_kegg',
]
