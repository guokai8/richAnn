"""
Data loading and annotation utilities for richAnn package
"""

import pandas as pd
import logging
from typing import Optional, List, Dict
from pathlib import Path

logger = logging.getLogger(__name__)


def load_gmt(filepath: str,
             annot_col: str = "GeneSet",
             name_col: str = "GeneSetName",
             gene_col: str = "GeneID") -> pd.DataFrame:
    """
    Load gene sets from GMT (Gene Matrix Transposed) format file

    GMT format: Each line contains:
    - Gene set name
    - Description (optional)
    - Tab-separated list of genes

    Parameters:
    -----------
    filepath : str
        Path to GMT file
    annot_col : str
        Name for annotation ID column (default: "GeneSet")
    name_col : str
        Name for annotation name/description column (default: "GeneSetName")
    gene_col : str
        Name for gene ID column (default: "GeneID")

    Returns:
    --------
    pd.DataFrame with columns: [annot_col, name_col, gene_col]

    Examples:
    ---------
    >>> # Load MSigDB GMT file
    >>> geneset_db = load_gmt("c2.cp.kegg.v7.4.symbols.gmt")
    >>> result = richGSEA(gene_scores, geneset_db)
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"GMT file not found: {filepath}")

    records = []

    with open(filepath, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')

            if len(parts) < 3:
                logger.warning(f"Line {line_num}: Expected at least 3 fields, got {len(parts)}. Skipping.")
                continue

            geneset_id = parts[0]
            geneset_name = parts[1] if parts[1] else geneset_id
            genes = parts[2:]

            for gene in genes:
                gene = gene.strip()
                if gene:
                    records.append({
                        annot_col: geneset_id,
                        name_col: geneset_name,
                        gene_col: gene
                    })

    if not records:
        raise ValueError(f"No valid gene sets found in {filepath}")

    df = pd.DataFrame(records)
    logger.info(f"Loaded {df[annot_col].nunique()} gene sets with {df[gene_col].nunique()} unique genes from {filepath}")

    return df


def load_go_gaf(filepath: str,
                ontology: Optional[str] = None,
                evidence_codes: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Load GO annotations from GAF (Gene Association File) format

    GAF format is the standard GO annotation format. This function parses
    the tab-delimited file and extracts gene-GO term associations.

    Parameters:
    -----------
    filepath : str
        Path to GAF file (can be .gaf or .gaf.gz)
    ontology : str, optional
        Filter by ontology: "BP", "MF", or "CC". If None, loads all.
    evidence_codes : list of str, optional
        Filter by evidence codes (e.g., ["IEA", "IDA"]). If None, keeps all.

    Returns:
    --------
    pd.DataFrame with columns: GeneID, GOterm, GOname, Ontology, Evidence

    Examples:
    ---------
    >>> # Load all GO annotations
    >>> go_data = load_go_gaf("goa_human.gaf")
    >>>
    >>> # Load only Biological Process with experimental evidence
    >>> go_bp = load_go_gaf("goa_human.gaf", ontology="BP",
    ...                     evidence_codes=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"])
    """
    import gzip

    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"GAF file not found: {filepath}")

    # Determine if file is gzipped
    open_func = gzip.open if filepath.suffix == '.gz' else open
    mode = 'rt' if filepath.suffix == '.gz' else 'r'

    records = []
    ontology_map = {'P': 'BP', 'F': 'MF', 'C': 'CC'}

    with open_func(filepath, mode, encoding='utf-8') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('!'):
                continue

            parts = line.strip().split('\t')

            # GAF 2.x has at least 15 columns
            if len(parts) < 15:
                continue

            gene_id = parts[2]  # Gene symbol/ID
            go_id = parts[4]    # GO ID
            evidence = parts[6]  # Evidence code
            aspect = parts[8]    # Ontology aspect (P/F/C)
            go_name = parts[9] if len(parts) > 9 else go_id  # GO term name

            # Apply filters
            onto = ontology_map.get(aspect, aspect)
            if ontology and onto != ontology:
                continue
            if evidence_codes and evidence not in evidence_codes:
                continue

            records.append({
                'GeneID': gene_id,
                'GOterm': go_id,
                'GOname': go_name if go_name else go_id,
                'Ontology': onto,
                'Evidence': evidence
            })

    if not records:
        raise ValueError(f"No valid GO annotations found in {filepath}")

    df = pd.DataFrame(records)
    logger.info(f"Loaded {len(df)} GO annotations for {df['GeneID'].nunique()} genes from {filepath}")

    return df


def load_kegg_mapping(filepath: str,
                      organism: str = "hsa") -> pd.DataFrame:
    """
    Load KEGG pathway annotations from a mapping file

    Expected format: Tab-delimited file with columns:
    - Gene ID
    - Pathway ID (e.g., "hsa00010")
    - Pathway Name

    Parameters:
    -----------
    filepath : str
        Path to KEGG mapping file
    organism : str
        Organism code (e.g., "hsa" for human, "mmu" for mouse)

    Returns:
    --------
    pd.DataFrame with columns: GeneID, Pathway, PathwayName

    Examples:
    ---------
    >>> # Load human KEGG pathways
    >>> kegg_data = load_kegg_mapping("kegg_human.txt", organism="hsa")
    >>> result = richKEGG(genes, kegg_data)
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"KEGG mapping file not found: {filepath}")

    try:
        df = pd.read_csv(filepath, sep='\t', header=0)

        # Standardize column names
        expected_cols = ['GeneID', 'Pathway', 'PathwayName']
        if len(df.columns) >= 3:
            df.columns = expected_cols[:len(df.columns)]
        else:
            raise ValueError(f"Expected at least 3 columns, got {len(df.columns)}")

        # Filter by organism if specified
        if organism:
            df = df[df['Pathway'].str.startswith(organism)]

        if len(df) == 0:
            raise ValueError(f"No pathways found for organism '{organism}'")

        logger.info(f"Loaded {df['Pathway'].nunique()} KEGG pathways for {df['GeneID'].nunique()} genes")

        return df[['GeneID', 'Pathway', 'PathwayName']]

    except Exception as e:
        raise ValueError(f"Error parsing KEGG mapping file: {e}")


def create_custom_annotation(gene_sets: Dict[str, List[str]],
                             annot_col: str = "Annot",
                             name_col: str = "Term",
                             gene_col: str = "GeneID") -> pd.DataFrame:
    """
    Create annotation DataFrame from dictionary of gene sets

    Parameters:
    -----------
    gene_sets : dict
        Dictionary mapping annotation IDs to lists of gene IDs
        Keys: Annotation IDs
        Values: Lists of gene IDs
    annot_col : str
        Column name for annotation IDs
    name_col : str
        Column name for annotation names (will use IDs if not provided)
    gene_col : str
        Column name for gene IDs

    Returns:
    --------
    pd.DataFrame suitable for enrichment analysis

    Examples:
    ---------
    >>> # Create custom gene sets
    >>> custom_sets = {
    ...     "DNA_REPAIR": ["BRCA1", "BRCA2", "TP53", "ATM"],
    ...     "CELL_CYCLE": ["CDK1", "CDK2", "CCNA1", "CCNB1"]
    ... }
    >>> annot = create_custom_annotation(custom_sets)
    """
    if not gene_sets:
        raise ValueError("gene_sets dictionary cannot be empty")

    records = []
    for annot_id, genes in gene_sets.items():
        if not genes:
            logger.warning(f"Gene set '{annot_id}' is empty, skipping")
            continue

        for gene in genes:
            gene = str(gene).strip()
            if gene:
                records.append({
                    annot_col: annot_id,
                    name_col: annot_id,  # Use ID as name
                    gene_col: gene
                })

    if not records:
        raise ValueError("No valid genes found in any gene set")

    df = pd.DataFrame(records)
    logger.info(f"Created annotation with {df[annot_col].nunique()} gene sets and {df[gene_col].nunique()} unique genes")

    return df


def from_pathwaydb_go(go_db, ontology: Optional[str] = None) -> pd.DataFrame:
    """
    Convert pathwaydb GOAnnotationDB output to richAnn GO format.

    Parameters:
    -----------
    go_db : pathwaydb.GOAnnotationDB or pd.DataFrame
        Either a GOAnnotationDB instance or a DataFrame from go_db.to_dataframe()
    ontology : str, optional
        Filter by ontology: "BP", "MF", or "CC". If None, loads all.

    Returns:
    --------
    pd.DataFrame with columns: GeneID, GOterm, GOname, Ontology

    Examples:
    ---------
    >>> from pathwaydb import GOAnnotationDB
    >>> go_db = GOAnnotationDB('go_human.db')
    >>> go_data = from_pathwaydb_go(go_db, ontology="BP")
    >>> result = richGO(genes, go_data, ontology="BP")
    """
    # Handle both GOAnnotationDB instance and DataFrame
    if hasattr(go_db, 'to_dataframe'):
        # It's a GOAnnotationDB instance
        records = go_db.to_dataframe()
        df = pd.DataFrame(records)
    elif isinstance(go_db, pd.DataFrame):
        df = go_db.copy()
    elif isinstance(go_db, list):
        # List of dicts from to_dataframe()
        df = pd.DataFrame(go_db)
    else:
        raise TypeError(
            f"Expected GOAnnotationDB, DataFrame, or list of dicts, got {type(go_db)}"
        )

    # Validate required columns from pathwaydb format
    required = ['GeneID', 'TERM', 'Aspect']
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Expected pathwaydb GO format.")

    # Map Aspect (P/F/C) to Ontology (BP/MF/CC)
    aspect_map = {'P': 'BP', 'F': 'MF', 'C': 'CC'}
    df['Ontology'] = df['Aspect'].map(aspect_map).fillna(df['Aspect'])

    # Filter by ontology if specified
    if ontology:
        if ontology not in ['BP', 'MF', 'CC']:
            raise ValueError("ontology must be 'BP', 'MF', or 'CC'")
        df = df[df['Ontology'] == ontology]

    # Get term names if available (pathwaydb may have term_name column)
    if 'term_name' in df.columns:
        df['GOname'] = df['term_name'].fillna(df['TERM'])
    else:
        # Use GO ID as name placeholder - user should populate term names
        df['GOname'] = df['TERM']

    # Rename columns to richAnn format
    result = df[['GeneID', 'TERM', 'GOname', 'Ontology']].copy()
    result.columns = ['GeneID', 'GOterm', 'GOname', 'Ontology']

    # Remove duplicates
    result = result.drop_duplicates()

    logger.info(f"Converted {result['GOterm'].nunique()} GO terms for {result['GeneID'].nunique()} genes")

    return result


def from_pathwaydb_kegg(kegg_db) -> pd.DataFrame:
    """
    Convert pathwaydb KEGGAnnotationDB output to richAnn KEGG format.

    Parameters:
    -----------
    kegg_db : pathwaydb.KEGGAnnotationDB or pd.DataFrame
        Either a KEGGAnnotationDB instance or a DataFrame from kegg_db.to_dataframe()

    Returns:
    --------
    pd.DataFrame with columns: GeneID, Pathway, PathwayName, Level1, Level2, Level3

    The Level columns contain KEGG pathway hierarchy information:
    - Level1: Top-level category (e.g., "Metabolism", "Human Diseases")
    - Level2: Sub-category (e.g., "Carbohydrate metabolism", "Cancer")
    - Level3: Pathway name (same as PathwayName)

    Examples:
    ---------
    >>> from pathwaydb import KEGGAnnotationDB
    >>> kegg_db = KEGGAnnotationDB('kegg_human.db')
    >>> kegg_data = from_pathwaydb_kegg(kegg_db)
    >>> result = richKEGG(genes, kegg_data)
    >>> # Results will include Level1, Level2, Level3 columns
    """
    # Handle both KEGGAnnotationDB instance and DataFrame
    if hasattr(kegg_db, 'to_dataframe'):
        # It's a KEGGAnnotationDB instance
        records = kegg_db.to_dataframe()
        df = pd.DataFrame(records)
    elif isinstance(kegg_db, pd.DataFrame):
        df = kegg_db.copy()
    elif isinstance(kegg_db, list):
        # List of dicts from to_dataframe()
        df = pd.DataFrame(kegg_db)
    else:
        raise TypeError(
            f"Expected KEGGAnnotationDB, DataFrame, or list of dicts, got {type(kegg_db)}"
        )

    # Validate required columns from pathwaydb format
    required = ['GeneID', 'PATH', 'Annot']
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Expected pathwaydb KEGG format.")

    # Check if hierarchy columns exist
    has_hierarchy = all(col in df.columns for col in ['Level1', 'Level2', 'Level3'])

    # Select columns to include
    cols_to_include = ['GeneID', 'PATH', 'Annot']
    if has_hierarchy:
        cols_to_include.extend(['Level1', 'Level2', 'Level3'])

    result = df[cols_to_include].copy()

    # Rename columns to richAnn format
    rename_map = {'PATH': 'Pathway', 'Annot': 'PathwayName'}
    result = result.rename(columns=rename_map)

    # Remove duplicates
    result = result.drop_duplicates()

    # Log info about hierarchy
    if has_hierarchy:
        n_with_hierarchy = result['Level1'].notna().sum()
        logger.info(f"Converted {result['Pathway'].nunique()} KEGG pathways for {result['GeneID'].nunique()} genes ({n_with_hierarchy} with hierarchy)")
    else:
        logger.info(f"Converted {result['Pathway'].nunique()} KEGG pathways for {result['GeneID'].nunique()} genes (no hierarchy)")

    return result


def validate_annotation_format(df: pd.DataFrame,
                               required_cols: List[str],
                               check_duplicates: bool = True,
                               check_empty: bool = True) -> Dict[str, any]:
    """
    Validate annotation DataFrame and return quality metrics

    Parameters:
    -----------
    df : pd.DataFrame
        Annotation DataFrame to validate
    required_cols : list of str
        Required column names
    check_duplicates : bool
        Check for duplicate rows
    check_empty : bool
        Check for empty values

    Returns:
    --------
    dict with validation results and statistics

    Examples:
    ---------
    >>> stats = validate_annotation_format(go_data, ['GeneID', 'GOterm', 'GOname', 'Ontology'])
    >>> print(f"Valid: {stats['valid']}")
    >>> print(f"Unique genes: {stats['n_genes']}")
    """
    results = {'valid': True, 'errors': [], 'warnings': []}

    # Check required columns
    missing_cols = set(required_cols) - set(df.columns)
    if missing_cols:
        results['valid'] = False
        results['errors'].append(f"Missing required columns: {missing_cols}")
        return results

    # Check for empty DataFrame
    if len(df) == 0:
        results['valid'] = False
        results['errors'].append("DataFrame is empty")
        return results

    # Check for empty values
    if check_empty:
        for col in required_cols:
            null_count = df[col].isna().sum()
            if null_count > 0:
                results['warnings'].append(f"Column '{col}' has {null_count} null values")

            empty_count = (df[col].astype(str).str.strip() == '').sum()
            if empty_count > 0:
                results['warnings'].append(f"Column '{col}' has {empty_count} empty strings")

    # Check for duplicates
    if check_duplicates:
        dup_count = df.duplicated().sum()
        if dup_count > 0:
            results['warnings'].append(f"Found {dup_count} duplicate rows")

    # Collect statistics
    annot_col = required_cols[0]  # Assume first col is annotation ID
    gene_col = required_cols[-1] if len(required_cols) > 1 else None  # Assume last col is gene ID

    results['n_rows'] = len(df)
    results['n_annotations'] = df[annot_col].nunique()
    if gene_col:
        results['n_genes'] = df[gene_col].nunique()
        results['avg_genes_per_annot'] = len(df) / results['n_annotations']

    return results
