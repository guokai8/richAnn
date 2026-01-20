"""
Pytest fixtures for richAnn tests
"""

import pytest
import pandas as pd
import numpy as np


@pytest.fixture
def sample_genes():
    """Sample gene list for testing"""
    return ["BRCA1", "BRCA2", "TP53", "ATM", "CHEK2", "RAD51", "PALB2", "BARD1"]


@pytest.fixture
def sample_gene_scores():
    """Sample gene scores for GSEA testing"""
    genes = ["BRCA1", "BRCA2", "TP53", "ATM", "CHEK2", "RAD51", "PALB2", "BARD1",
             "CDK1", "CDK2", "CCNA1", "CCNB1", "MYC", "JUN", "FOS", "EGFR"]
    scores = [2.5, 2.3, 2.1, 1.9, 1.7, 1.5, 1.3, 1.1,
              -0.5, -0.7, -0.9, -1.1, -1.3, -1.5, -1.7, -1.9]
    return dict(zip(genes, scores))


@pytest.fixture
def sample_go_data():
    """Sample GO annotation data for testing"""
    data = {
        'GeneID': ['BRCA1', 'BRCA2', 'TP53', 'ATM', 'CHEK2', 'RAD51',
                   'BRCA1', 'BRCA2', 'TP53',
                   'CDK1', 'CDK2', 'CCNA1', 'CCNB1',
                   'MYC', 'JUN', 'FOS'],
        'GOterm': ['GO:0006281', 'GO:0006281', 'GO:0006281', 'GO:0006281', 'GO:0006281', 'GO:0006281',
                   'GO:0006974', 'GO:0006974', 'GO:0006974',
                   'GO:0007049', 'GO:0007049', 'GO:0007049', 'GO:0007049',
                   'GO:0006355', 'GO:0006355', 'GO:0006355'],
        'GOname': ['DNA repair', 'DNA repair', 'DNA repair', 'DNA repair', 'DNA repair', 'DNA repair',
                   'DNA damage response', 'DNA damage response', 'DNA damage response',
                   'cell cycle', 'cell cycle', 'cell cycle', 'cell cycle',
                   'transcription regulation', 'transcription regulation', 'transcription regulation'],
        'Ontology': ['BP', 'BP', 'BP', 'BP', 'BP', 'BP',
                     'BP', 'BP', 'BP',
                     'BP', 'BP', 'BP', 'BP',
                     'BP', 'BP', 'BP']
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_kegg_data():
    """Sample KEGG pathway annotation data for testing"""
    data = {
        'GeneID': ['BRCA1', 'BRCA2', 'TP53', 'ATM',
                   'CDK1', 'CDK2', 'CCNA1', 'CCNB1',
                   'MYC', 'JUN', 'FOS', 'EGFR'],
        'Pathway': ['hsa03440', 'hsa03440', 'hsa03440', 'hsa03440',
                    'hsa04110', 'hsa04110', 'hsa04110', 'hsa04110',
                    'hsa05200', 'hsa05200', 'hsa05200', 'hsa05200'],
        'PathwayName': ['Homologous recombination', 'Homologous recombination', 'Homologous recombination', 'Homologous recombination',
                        'Cell cycle', 'Cell cycle', 'Cell cycle', 'Cell cycle',
                        'Pathways in cancer', 'Pathways in cancer', 'Pathways in cancer', 'Pathways in cancer']
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_geneset_db():
    """Sample gene set database for GSEA testing"""
    data = {
        'GeneSet': ['DNA_REPAIR', 'DNA_REPAIR', 'DNA_REPAIR', 'DNA_REPAIR', 'DNA_REPAIR', 'DNA_REPAIR',
                    'CELL_CYCLE', 'CELL_CYCLE', 'CELL_CYCLE', 'CELL_CYCLE',
                    'TRANSCRIPTION', 'TRANSCRIPTION', 'TRANSCRIPTION', 'TRANSCRIPTION'],
        'GeneSetName': ['DNA repair pathway', 'DNA repair pathway', 'DNA repair pathway', 'DNA repair pathway', 'DNA repair pathway', 'DNA repair pathway',
                        'Cell cycle regulation', 'Cell cycle regulation', 'Cell cycle regulation', 'Cell cycle regulation',
                        'Transcription factors', 'Transcription factors', 'Transcription factors', 'Transcription factors'],
        'GeneID': ['BRCA1', 'BRCA2', 'TP53', 'ATM', 'CHEK2', 'RAD51',
                   'CDK1', 'CDK2', 'CCNA1', 'CCNB1',
                   'MYC', 'JUN', 'FOS', 'EGFR']
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_enrich_result(sample_go_data, sample_genes):
    """Sample EnrichResult object for testing"""
    from richAnn import richGO

    try:
        result = richGO(sample_genes, sample_go_data, ontology="BP", pvalue=1.0)
        return result
    except ValueError:
        # If no significant results, return minimal result
        from richAnn.core import EnrichResult
        df = pd.DataFrame({
            'Annot': ['GO:0006281'],
            'Term': ['DNA repair'],
            'Pvalue': [0.01],
            'Padj': [0.05],
            'Count': [6],
            'GeneID': ['BRCA1;BRCA2;TP53;ATM;CHEK2;RAD51'],
            'RichFactor': [2.0],
            'GeneRatio': ['6/8'],
            'BgRatio': ['6/16']
        })
        return EnrichResult(df, enrichment_type="GO")
